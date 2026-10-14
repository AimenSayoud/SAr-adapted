"""X-111 — forward model first: can a forward model carry the mat's level across winters at P6, where the chain loses it?

Candidates for the laser's level step across each winter of X-099 (the plain chain's epochs): buoyancy on the logger
water table, the ERA5 bucket × buoyancy, SPAMS, zero change — each fitted without the two years of the winter it
predicts — against the chain's own winter step. Rule: ``08_deliverables/insar_full_test/preregistration/X-111.yaml``.
Writes only to the hub: ``08_deliverables/insar_full_test/winter_forward/``.

    make run E=winter_forward_x111
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import spams_baseline_x112 as B  # noqa: E402  (the daily table of X-074 + ERA5 ET0, the g fit)

from insar_wetlands import fusion as fu  # noqa: E402
from insar_wetlands import spams  # noqa: E402
from insar_wetlands.paths import hub_root  # noqa: E402
from insar_wetlands.scoring import Protocol, assert_preregistered, guard_lockbox  # noqa: E402

HUB = hub_root()
FT = HUB / "08_deliverables" / "insar_full_test"
OUT = FT / "winter_forward"
PREREG = FT / "preregistration" / "X-111.yaml"
LAG = 6


def los_factor() -> dict[str, float]:
    """cos(incidence) at P6 per track: X-099's winter steps (laser and chain) are line-of-sight, the forward models
    predict vertical change — every model prediction is projected before it is compared (correction of 2026-10-14)."""
    import xarray as xr
    from insar_wetlands.paths import cube_root
    st = xr.open_dataset(cube_root() / "silver" / "static_40m.nc")
    px = pd.read_csv(HUB / "08_deliverables" / "field_first" / "plot_pixels.csv").set_index("plot")
    r, c = int(px.loc["P6", "row"]), int(px.loc["P6", "col"])
    return {"ascending": float(np.cos(np.radians(st.incidence_asc_deg.values[r, c]))),
            "descending": float(np.cos(np.radians(st.incidence_des_deg.values[r, c])))}


def at(s: pd.Series, day: pd.Timestamp, half: int = 3) -> float:
    w = s[(s.index >= day - pd.Timedelta(days=half)) & (s.index <= day + pd.Timedelta(days=half))].dropna()
    return float(w.mean()) if len(w) else np.nan


def bucket_winter(d: pd.DataFrame, y: int, train: np.ndarray) -> pd.Series:
    rain, pet, w = d.rain_era5_mm.fillna(0).to_numpy(), d.pet_hamon.to_numpy(), d.wtd_cm.to_numpy()
    par = fu.fit_bucket(np.where(train, w, np.nan), rain, pet)
    floor = np.nanpercentile(w[train], 1)
    m = ((d.index >= pd.Timestamp(f"{y}-09-01")) & (d.index <= pd.Timestamp(f"{y + 1}-06-30")))
    i0 = np.flatnonzero(m)[0]
    prev = w[:i0 + 1][np.isfinite(w[:i0 + 1])]
    w0 = prev[-1] if len(prev) else par["w0"]
    out = pd.Series(np.nan, index=d.index)
    out.iloc[np.flatnonzero(m)] = fu.simulate_bucket(par, w0, rain[m], pet[m], floor=floor)
    return out


def main() -> None:
    assert_preregistered(PREREG, [OUT / "verdict.csv"])
    OUT.mkdir(parents=True, exist_ok=True)
    prot = Protocol.load(FT / "protocol.yaml")
    d = B.table()
    d = d[d.year.between(2020, 2025)]
    guard_lockbox(d.index[d.year >= 2021], prot)
    P, E = d.rain_era5_mm.to_numpy(), d.et0_mm.to_numpy()
    valid6 = (d.dh_mm.notna() & ~d.frozen_win).to_numpy()
    base = pd.read_csv(FT / "baseline" / "winters.csv", parse_dates=["from", "to"])
    base = base[base.estimate == "plain chain"]
    rows = []
    los = los_factor()
    for (w, track), g0 in base.groupby(["winter", "track"]):
        y = int(w[:4])
        t0, t1 = g0["from"].iloc[0], g0["to"].iloc[0]
        laser = float(g0.laser_change_mm.iloc[0])
        train = ~d.year.isin([y, y + 1]).to_numpy()
        g = B.g_fit(d, train & valid6)
        par = spams.fit(P, E, d.dh_mm.to_numpy(), train & valid6, LAG)
        R = pd.Series(spams.reversible(P, E, par["x_p"], par["x_e"], par["tau"]), index=d.index)
        bk = bucket_winter(d, y, train)
        k = los[track]                         # vertical model → the track's line of sight
        cand = {"zero change": 0.0,
                "logger water table × g": (at(d.wtd_cm, t1) - at(d.wtd_cm, t0)) * g * k,
                "bucket × g": (at(bk, t1) - at(bk, t0)) * g * k,
                "SPAMS": (at(R, t1) - at(R, t0)) * k}
        for ref, gr in g0.groupby("reference"):
            cand[f"plain chain ({ref})"] = float(gr.insar_change_mm.iloc[0])
        for name, v in cand.items():
            rows.append({"winter": w, "track": track, "from": t0.date(), "to": t1.date(), "candidate": name,
                         "laser_step_mm": laser, "predicted_step_mm": v,
                         "error_mm": v - laser if np.isfinite(v) else np.nan, "g_mm_per_cm": g, "spams_tau": par["tau"]})
    steps = pd.DataFrame(rows)
    summ = (steps.dropna(subset=["error_mm"]).groupby("candidate")
            .agg(steps=("error_mm", "size"), mean_abs_error_mm=("error_mm", lambda e: float(np.mean(np.abs(e)))),
                 sign_agreement=("predicted_step_mm", "size")).reset_index())
    sa = (steps.dropna(subset=["error_mm"]).assign(ok=lambda x: np.sign(x.predicted_step_mm) == np.sign(x.laser_step_mm))
          .groupby("candidate").ok.mean())
    summ["sign_agreement"] = summ.candidate.map(sa)
    verdict = decide(steps)
    for name, df in (("steps", steps), ("summary", summ), ("verdict", verdict)):
        df.to_csv(OUT / f"{name}.csv", index=False)
    (OUT / "README.md").write_text(readme(steps, summ, verdict), encoding="utf-8")
    print(verdict.to_string(index=False))
    print(f"wrote {OUT}")


def decide(steps: pd.DataFrame) -> pd.DataFrame:
    """Per candidate, on the steps it covers: mean |error| vs the hard-target chain's and zero change's on the same steps."""
    rows = []
    chain = steps[steps.candidate == "plain chain (hard_targets)"].set_index(["winter", "track"]).error_mm
    zero = steps[steps.candidate == "zero change"].set_index(["winter", "track"]).error_mm
    for name in ("logger water table × g", "bucket × g", "SPAMS"):
        e = steps[steps.candidate == name].dropna(subset=["error_mm"]).set_index(["winter", "track"]).error_mm
        k = e.index.intersection(chain.index)
        mae, mc, mz = (float(np.mean(np.abs(x.loc[k]))) for x in (e, chain, zero))
        rows.append({"candidate": name, "steps": len(k), "mae_candidate_mm": mae, "mae_chain_mm": mc, "mae_zero_mm": mz,
                     "carries_the_winter": bool(mae <= 0.5 * mc and mae < mz)})
    return pd.DataFrame(rows)


def _md(df: pd.DataFrame, nd: int = 1) -> str:
    d = df.copy()
    for col in d.columns:
        if pd.api.types.is_float_dtype(d[col]):
            d[col] = d[col].map(lambda x: "—" if pd.isna(x) else f"{x:.{nd}f}")
    return "\n".join(["| " + " | ".join(map(str, d.columns)) + " |", "|" + "---|" * len(d.columns)] +
                     ["| " + " | ".join(map(str, r)) + " |" for r in d.itertuples(index=False)])


def readme(steps, summ, verdict) -> str:
    wide = steps.pivot_table(index=["winter", "track", "laser_step_mm"], columns="candidate", values="predicted_step_mm").reset_index()
    return f"""# X-111 — can a forward model carry the mat's level across winters? (generated)

Generated by `experiments/winter_forward_x111.py` (`make run E=winter_forward_x111`); do not edit. Pre-registration:
`../preregistration/X-111.yaml`. P6, the winter steps of X-099 (the chain's epochs, same laser datum). Each candidate
fitted without the two years of the winter it predicts. A candidate carries the winter if its mean absolute step error
is at most half the hard-target chain's and below zero change's, on the steps it covers. Few steps: read step by step. Model predictions (vertical) are projected to each track's line of sight
(cos of P6's incidence), as the laser and chain steps are.

## Verdict

{_md(verdict)}

## Mean absolute error and sign agreement, all steps each candidate covers

{_md(summ, 2)}

## Each step (predicted level step, mm; the laser's in the third column)

{_md(wide)}
"""


if __name__ == "__main__":
    main()

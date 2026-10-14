"""X-116 — is the winter shortfall of buoyancy a regression-dilution artefact of noisy 6-day water-table changes?

g (laser change on logger water-table change) fitted on within-season intervals of 6 to 96 days, year-blocked; then
the X-099 winter steps predicted with the long-interval g, against the laser, the chain and X-111's 6-day g. Rule:
``08_deliverables/insar_full_test/preregistration/X-116.yaml``. Writes only to the hub:
``08_deliverables/insar_full_test/winter_physics/``.

    make run E=winter_physics_x116
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import spams_baseline_x112 as B  # noqa: E402
import winter_forward_x111 as W  # noqa: E402

from insar_wetlands.paths import hub_root  # noqa: E402
from insar_wetlands.scoring import Protocol, assert_preregistered, guard_lockbox  # noqa: E402

HUB = hub_root()
FT = HUB / "08_deliverables" / "insar_full_test"
OUT = FT / "winter_physics"
PREREG = FT / "preregistration" / "X-116.yaml"
LENGTHS = (6, 12, 24, 48, 96)
YEARS = (2021, 2022, 2023, 2024)            # the years with both a laser and a logger


def intervals(d: pd.DataFrame, L: int, prot: Protocol) -> pd.DataFrame:
    lo, hi = prot.season_months
    gaps = [(pd.Timestamp(a).tz_convert(None), pd.Timestamp(b).tz_convert(None)) for a, b in prot.datum_breaks]
    a = pd.DataFrame({"t0": d.index, "h0": d.laser_cm.to_numpy(), "w0": d.wtd_cm.to_numpy()})
    b = pd.DataFrame({"h1": d.laser_cm.shift(-L).to_numpy(), "w1": d.wtd_cm.shift(-L).to_numpy(),
                      "fz": d.frozen.rolling(L + 1).max().shift(-L).to_numpy()})
    x = pd.concat([a, b], axis=1)
    x["t1"] = x.t0 + pd.Timedelta(days=L)
    x["year"] = x.t0.dt.year
    ok = (x.t0.dt.month.between(lo, hi) & x.t1.dt.month.between(lo, hi) & (x.t0.dt.year == x.t1.dt.year)
          & x[["h0", "h1", "w0", "w1"]].notna().all(1) & (x.fz.fillna(1) == 0) & x.year.isin(YEARS))
    for g0, g1 in gaps:
        ok &= ~((x.t0 < g1) & (x.t1 > g0))
    x = x[ok].copy()
    x["dh_mm"] = (x.h1 - x.h0) * 10
    x["dw_cm"] = x.w1 - x.w0
    return x


def slope(x: pd.DataFrame) -> float:
    return float(np.sum(x.dw_cm * x.dh_mm) / np.sum(x.dw_cm ** 2)) if len(x) >= 10 else np.nan


def main() -> None:
    assert_preregistered(PREREG, [OUT / "verdict.csv"])
    OUT.mkdir(parents=True, exist_ok=True)
    prot = Protocol.load(FT / "protocol.yaml")
    d = B.table()
    d = d[d.year.between(2020, 2025)]
    guard_lockbox(d.index[d.year >= 2021], prot)
    iv = {L: intervals(d, L, prot) for L in LENGTHS}
    rows = []
    for L, x in iv.items():
        rows.append({"length_days": L, "held_out_year": "pooled", "g_mm_per_cm": slope(x), "intervals": len(x),
                     "sd_dw_cm": float(x.dw_cm.std())})
        for y in YEARS:
            rows.append({"length_days": L, "held_out_year": y, "g_mm_per_cm": slope(x[x.year != y]),
                         "intervals": int((x.year != y).sum()), "sd_dw_cm": float(x[x.year != y].dw_cm.std())})
    gl = pd.DataFrame(rows)
    med = gl[gl.held_out_year != "pooled"].groupby("length_days").g_mm_per_cm.median()
    base = pd.read_csv(FT / "baseline" / "winters.csv", parse_dates=["from", "to"])
    base = base[(base.estimate == "plain chain") & (base.reference == "hard_targets")]
    st = []
    for _, r in base.iterrows():
        y = int(r.winter[:4])
        dw = W.at(d.wtd_cm, r["to"]) - W.at(d.wtd_cm, r["from"])
        if not np.isfinite(dw):
            continue
        rec = {"winter": r.winter, "track": r.track, "laser_step_mm": r.laser_change_mm, "chain_step_mm": r.insar_change_mm,
               "dw_cm": dw}
        for L, x in iv.items():
            g = slope(x[~x.year.isin([y, y + 1])])
            rec[f"pred_g{L}_mm"] = dw * g
        st.append(rec)
    steps = pd.DataFrame(st)
    e = lambda col: float(np.mean(np.abs(steps[col] - steps.laser_step_mm)))
    mae = pd.DataFrame([{"estimate": f"logger × g({L} d)", "mae_mm": e(f"pred_g{L}_mm")} for L in LENGTHS]
                       + [{"estimate": "plain chain (hard targets)", "mae_mm": e("chain_step_mm")}])
    m96, m6, mc = e("pred_g96_mm"), e("pred_g6_mm"), e("chain_step_mm")
    inc = all(med.iloc[i + 1] > med.iloc[i] for i in range(len(med) - 1))
    verdict = pd.DataFrame([{"steps": len(steps), "mae_g96_mm": m96, "mae_g6_mm": m6, "mae_chain_mm": mc,
                             "carries_the_winter": bool(m96 <= 0.5 * mc and m96 < m6),
                             "g6_median": float(med.loc[6]), "g96_median": float(med.loc[96]),
                             "dilution_supported": bool(inc or med.loc[96] >= 1.25 * med.loc[6])}])
    for name, df in (("g_by_length", gl), ("winter_steps", steps), ("winter_mae", mae), ("verdict", verdict)):
        df.to_csv(OUT / f"{name}.csv", index=False)
    (OUT / "README.md").write_text(readme(gl, med, steps, mae, verdict), encoding="utf-8")
    print(verdict.to_string(index=False))
    print(f"wrote {OUT}")


def _md(df: pd.DataFrame, nd: int = 2) -> str:
    d = df.copy()
    for col in d.columns:
        if pd.api.types.is_float_dtype(d[col]):
            d[col] = d[col].map(lambda x: "—" if pd.isna(x) else f"{x:.{nd}f}")
    return "\n".join(["| " + " | ".join(map(str, d.columns)) + " |", "|" + "---|" * len(d.columns)] +
                     ["| " + " | ".join(map(str, r)) + " |" for r in d.itertuples(index=False)])


def readme(gl, med, steps, mae, verdict) -> str:
    pooled = gl[gl.held_out_year == "pooled"][["length_days", "g_mm_per_cm", "intervals", "sd_dw_cm"]]
    return f"""# X-116 — the buoyancy coefficient against interval length, and the winter (generated)

Generated by `experiments/winter_physics_x116.py` (`make run E=winter_physics_x116`); do not edit. Pre-registration:
`../preregistration/X-116.yaml`. P6, within-season (May–November) intervals of one year on one laser datum, no frozen
soil; g = slope through the origin of laser change on logger water-table change, year-blocked. Winter steps: X-099's,
hard-target chain, g fitted without the winter's two years. Overlapping intervals are correlated — read the
year-blocked medians, not the interval counts.

## Verdict

{_md(verdict)}

## g by interval length (pooled; the median over held-out years beside it)

{_md(pooled.assign(median_heldout=pooled.length_days.map(med)))}

## Winter steps — mean absolute error

{_md(mae, 1)}

## Each winter step (mm)

{_md(steps, 1)}
"""


if __name__ == "__main__":
    main()

"""X-093(a) — cut the real P6 chain where its phase leaves the water-table model, tie the segments to the model.

The plain chain of consecutive pairs follows the laser within a season and loses the level across winters (X-099),
inside chains HyP3 coherence calls locked (X-096, T57/T58). Option (a) of X-093 cuts on the phase itself: a pair
whose increment departs from the water-table prior by more than 3 robust sigma, or that is frozen or snow-covered, is
not integrated; each remaining segment is tied to the prior level by least squares (``phase_sim.tie_segments``; the
same rule X-095 simulated). No laser is read to place a cut.

Runs only if X-095's pre-registered gate passed (pre-registration ``stop_if``). Pre-registration:
``08_deliverables/insar_full_test/preregistration/X-093.yaml``. Writes only to the hub:
``08_deliverables/insar_full_test/segment_tie/``. Tried on the P6 laser 2021–2024 (D-028).

    make run E=segment_tie_x093
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import p6_validation_x082 as X  # noqa: E402
import phase_sim_x095 as S  # noqa: E402  (pair table with flags, buoyancy leave-one-year-out, prior level)

from insar_wetlands import phase_sim as ps  # noqa: E402
from insar_wetlands.paths import hub_root  # noqa: E402
from insar_wetlands.scoring import SEASON_OF_MONTH, Protocol, assert_preregistered, score  # noqa: E402

HUB = hub_root()
FT = HUB / "08_deliverables" / "insar_full_test"
OUT = FT / "segment_tie"
PREREG = FT / "preregistration" / "X-093.yaml"
GATE = FT / "phase_sim" / "gate.csv"
YEARS = (2021, 2022, 2023, 2024)
K_SIGMA = 3.0
REFS = {"hard_targets": "ref_hard_unw", "grassland": "ref_grassland_unw"}


def pairs() -> pd.DataFrame:
    p, _, _, _ = S.load()
    g = pd.read_csv(X.CUBE / "gold" / "plots_pair.csv.gz", usecols=["track", "grid", "plot", "pair",
                                                                    "ref_grassland_unw"], low_memory=False)
    g = g[(g["plot"] == "P6") & (g.grid == "40m")][["track", "pair", "ref_grassland_unw"]]
    p = p.merge(g, on=["track", "pair"], how="left")
    return p[p.t1.dt.year.isin(YEARS) & p.t2.dt.year.isin(YEARS)].reset_index(drop=True)


def cut_and_tie(q: pd.DataFrame, est: np.ndarray, prior_inc: np.ndarray, prior_level: pd.Series) -> tuple[np.ndarray, pd.DataFrame]:
    r = est - prior_inc
    season = q.t1.dt.month.map(SEASON_OF_MONTH).to_numpy()
    sig = {s: S.robust_sd(r[season == s]) for s in set(season)}
    big = np.array([np.isfinite(r[i]) and abs(r[i]) > K_SIGMA * sig[season[i]] for i in range(len(r))])
    cut = q.cut.to_numpy() | big
    out = est.copy()
    for b in S.blocks(q):
        out[b] = ps.tie_segments(q.t1.iloc[b], q.t2.iloc[b], est[b], ~cut[b], prior_level)
    why = pd.DataFrame({"t1": q.t1, "t2": q.t2, "season": season, "frozen_snow": q.cut.to_numpy(),
                        "phase_residual": big & ~q.cut.to_numpy()})
    return out, why


def main() -> None:
    assert_preregistered(PREREG, [OUT / "winters.csv"])
    OUT.mkdir(parents=True, exist_ok=True)
    gate = pd.read_csv(GATE)
    if not gate.passes.astype(bool).all():
        (OUT / "README.md").write_text(stopped(gate), encoding="utf-8")
        pd.DataFrame([{"stopped": True, "reason": "X-095 gate failed"}]).to_csv(OUT / "winters.csv", index=False)
        print("X-095 gate failed: X-093(a) not run (pre-registered stop)")
        return
    prot = Protocol.load(FT / "protocol.yaml")
    p = pairs()
    las, lev = X.laser()
    wtd = S.load()[3]
    gfit = S.buoyancy_loyo(p)
    tabs: dict[str, list] = {"levels": [], "winters": [], "increments": [], "cuts": []}
    for track in X.TRACKS:
        ti = S.track_inputs(p, track, las[track], lev, wtd, gfit)
        q = ti["q"]
        q = q[q.t1.dt.year.isin(YEARS) & q.t2.dt.year.isin(YEARS)].reset_index(drop=True)
        for ref, col in REFS.items():
            if col not in q:
                continue
            qq = q[np.isfinite(q.unw_1x1 - q[col])].reset_index(drop=True)    # no phase: the chain breaks, for both
            obs = ((qq.unw_1x1 - qq[col]) * X.PHASE_TO_MM).to_numpy()
            tied, why = cut_and_tie(qq, obs, qq.prior_inc.to_numpy(), ti["prior_level"])
            tabs["cuts"].append(why.assign(track=track, reference=ref))
            for name, est in (("plain chain", obs), ("cut + tie (a)", tied)):
                sc = score(qq.assign(est_mm=est, truth_inc_mm=qq.d_laser_los_mm), ti["truth_level"], prot)
                for k in ("levels", "winters", "increments"):
                    if k in sc and len(sc[k]):
                        tabs[k].append(sc[k].assign(estimate=name, track=track, reference=ref))
    out = {k: pd.concat(v, ignore_index=True) for k, v in tabs.items() if v}
    for k, df in out.items():
        df.to_csv(OUT / f"{k}.csv", index=False)
    verdict = decide(out)
    verdict.to_csv(OUT / "verdict.csv", index=False)
    (OUT / "README.md").write_text(readme(out, verdict), encoding="utf-8")
    print(f"wrote {OUT}")


def decide(out: dict[str, pd.DataFrame]) -> pd.DataFrame:
    w, lv = out["winters"], out["levels"]
    rows = []
    for (track, ref), g in w.groupby(["track", "reference"]):
        a = g[g.estimate == "plain chain"].difference_mm.abs().median()
        b = g[g.estimate == "cut + tie (a)"].difference_mm.abs().median()
        s = lv[(lv.track == track) & (lv.reference == ref) & lv.view.str.startswith("season")]
        sa = s[s.estimate == "plain chain"].rmse_detrended_mm.median()
        sb = s[s.estimate == "cut + tie (a)"].rmse_detrended_mm.median()
        rows.append({"track": track, "reference": ref, "winter_plain_mm": a, "winter_tied_mm": b,
                     "winter_reduction": 1 - b / a if a else np.nan, "season_rmse_plain_mm": sa,
                     "season_rmse_tied_mm": sb, "season_change": sb / sa - 1 if sa else np.nan})
    v = pd.DataFrame(rows)
    v["passes"] = (v.winter_reduction >= 0.30) & (v.season_change <= 0.10)
    return v


def _md(df: pd.DataFrame, nd: int = 2) -> str:
    d = df.copy()
    for c in d.columns:
        if pd.api.types.is_float_dtype(d[c]):
            d[c] = d[c].map(lambda x: "—" if pd.isna(x) else f"{x:.{nd}f}")
    return "\n".join(["| " + " | ".join(map(str, d.columns)) + " |", "|" + "---|" * len(d.columns)] +
                     ["| " + " | ".join(map(str, r)) + " |" for r in d.itertuples(index=False)])


def stopped(gate: pd.DataFrame) -> str:
    return f"""# X-093(a) — segment and tie on the real chain (generated)

Generated by `experiments/segment_tie_x093.py`; do not edit. Pre-registration: `../preregistration/X-093.yaml`.

**Not run: the pre-registered stop condition applies.** X-095's gate (`../phase_sim/gate.csv`) — cutting at
frozen/snow pairs and tying to the water-table prior must cut the simulated across-winter error by at least 30 % on
both tracks — failed:

{_md(gate)}

If tying to this prior cannot close the winter error in simulation, where the truth is known and the prior is the
only error source left, it cannot be expected to on the real chain.
"""


def readme(out, verdict) -> str:
    cuts = out["cuts"].groupby(["track", "reference", "season"]).agg(
        frozen_snow=("frozen_snow", "sum"), phase_residual=("phase_residual", "sum"), pairs=("t1", "size")).reset_index()
    w = out["winters"][["track", "reference", "estimate", "winter", "laser_change_mm", "insar_change_mm", "difference_mm"]]
    return f"""# X-093(a) — segment and tie on the real chain (generated)

Generated by `experiments/segment_tie_x093.py` (`make run E=segment_tie_x093`); do not edit. Pre-registration:
`../preregistration/X-093.yaml`. Tried on the P6 laser 2021–2024 only (D-028).

## Verdict (winter error lower by >= 30 % on both tracks, hard-target reference; within-season RMSE not worse by > 10 %)

{_md(verdict)}

## Winter steps

{_md(w)}

## Cuts by cause

{_md(cuts)}

Files: `verdict.csv`, `winters.csv`, `levels.csv`, `increments.csv`, `cuts.csv`.
"""


if __name__ == "__main__":
    main()

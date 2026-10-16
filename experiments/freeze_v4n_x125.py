"""X-125 — v4n: v4 with a fixed neighbour average in place of the interpolated point read, block chosen on the design
years and frozen before X-107.

v4n is v4 (X-106) except the phase read at P6: the mean unwrapped phase of a block of cells around v4's frozen read
position (``v4.block_cells`` / ``v4.block_mean``) — B2, the four cells around it, or B3, its cell and the eight
neighbours. SPAMS is v4's, copied from its frozen file, not refitted. Runs on the cube as delivered (the 2026-10-15
snapshot, as v4 and X-107), set here as the cube root. ``build_v4n`` is what X-107 runs once on the held-out years.
Rule: ``preregistration/X-125.yaml``.

    make run E=freeze_v4n_x125
"""
from __future__ import annotations

import os
import subprocess
import sys
from datetime import date
from pathlib import Path

from insar_wetlands.paths import hub_root

HUB = hub_root()
os.environ["RZECIN_CUBE_ROOT"] = str(HUB / "06_data" / "cube_snapshots" / "as_delivered_2026-10-15")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import xarray as xr  # noqa: E402
import yaml  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import freeze_v4_x106 as F  # noqa: E402  (v4's frozen parts, unchanged)
import p6_validation_x082 as X  # noqa: E402

from insar_wetlands import v4  # noqa: E402
from insar_wetlands.cube import cube_root  # noqa: E402
from insar_wetlands.cube.core import rad_to_los_mm  # noqa: E402
from insar_wetlands.scoring import (  # noqa: E402
    Protocol,
    assert_preregistered,
    chain,
    detrended_stats,
    level_views,
    winter_steps,
)

FT = HUB / "08_deliverables" / "insar_full_test"
V4DIR = FT / "v4"
PARAMS = V4DIR / "v4n_frozen_params.yaml"
OUT = V4DIR / "v4n_design_report"
PREREG = FT / "preregistration" / "X-125.yaml"
BLOCKS = {"B2": 2, "B3": 3}


def unw_block(track: str, x: float, y: float, size: int) -> dict[str, float]:
    """Per pair key, the block mean of HyP3's unwrapped phase (rad) around (x, y)."""
    out = {}
    for s in F.IFG_SETS:
        f = cube_root() / "silver" / f"ifg_{s}_{track}.nc"
        if not f.exists():
            continue
        with xr.open_dataset(f) as ds:
            rows, cols = v4.block_cells(ds.x.values, ds.y.values, x, y, size)
            v = v4.block_mean(ds.unw.values, rows, cols)
            k = pd.to_datetime(ds.t1.values).strftime("%Y%m%d") + "_" + pd.to_datetime(ds.t2.values).strftime("%Y%m%d")
            out.update({kk: float(vv) for kk, vv in zip(k, v)})
    return out


def build_v4n(track: str, frozen: dict, years: tuple[int, ...], relevel: bool = True,
              size: int | None = None) -> pd.DataFrame:
    """v4's ``build_v4`` with the block read: time, seg, est (vertical mm, uplift positive)."""
    size = size or int(frozen["read_block"]["size"])
    g = pd.read_csv(cube_root() / "gold" / "plots_pair.csv.gz", low_memory=False,
                    usecols=["track", "grid", "plot", "t1", "t2", "platforms", "ref_hard_unw", "los_u"])
    g = g[(g["plot"] == "P6") & (g.grid == "40m") & (g.track == track)].copy()
    for c in ("t1", "t2"):
        g[c] = pd.to_datetime(g[c], format="ISO8601", utc=True)
    g = g[X.consecutive(g) & g.t1.dt.year.isin(years) & g.t2.dt.year.isin(years)]
    g = g[~g.platforms.astype(str).str.contains("|".join(frozen["pairs"]["exclude_platforms"]))]
    pos = frozen["read_position"][track]
    u = unw_block(track, pos["x"], pos["y"], size)
    key = g.t1.dt.strftime("%Y%m%d") + "_" + g.t2.dt.strftime("%Y%m%d")
    g["y_vert_mm"] = -rad_to_los_mm(key.map(u).to_numpy() - g.ref_hard_unw.to_numpy()) / g.los_u.to_numpy()
    g = g.dropna(subset=["y_vert_mm"]).sort_values("t1")
    c = chain(g.assign(est_mm=g.y_vert_mm), "est_mm")
    if relevel:
        R = F.spams_level(frozen)
        at = lambda t: float(R.reindex([pd.Timestamp(t).tz_convert(None).normalize()]).iloc[0])  # noqa: E731
        c = v4.relevel(c, at, v4.winter_bounds(c.time))
    return c


def stats(build, frozen: dict, prot: Protocol, name: str, **kw) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    lev = X.laser()[1]
    rows, wins, npairs = [], [], {}
    for track in X.TRACKS:
        c = build(track, frozen, F.DESIGN, **kw)
        npairs[track] = int(len(c) - c.seg.nunique())
        c["truth"] = lev.reindex(pd.DatetimeIndex(c.time).round("h")).to_numpy()
        for view, dd in level_views(c, prot):
            if view.startswith("season") and len(dd) >= 4:
                rows.append({"track": track, "estimate": name, "view": view,
                             **detrended_stats(dd.time, dd.truth, dd.est)})
        z = c.dropna(subset=["truth"]).copy()
        for _, cs in z.groupby("seg"):
            z.loc[cs.index, "truth"] -= cs.truth.iloc[0]
            z.loc[cs.index, "est"] -= cs.est.iloc[0]
        wins.append(winter_steps(z, prot).assign(track=track, estimate=name))
    return pd.DataFrame(rows), pd.concat(wins, ignore_index=True), npairs


def main() -> None:
    assert_preregistered(PREREG, [OUT / "summary.csv"])
    OUT.mkdir(parents=True, exist_ok=True)
    prot = Protocol.load(FT / "protocol.yaml")
    base = yaml.safe_load(F.PARAMS.read_text())                      # v4's frozen file: read, never rewritten
    S, W, N = [], [], {}
    s, w, n = stats(F.build_v4, base, prot, "v4")
    S.append(s), W.append(w)
    N["v4"] = n
    for name, size in BLOCKS.items():
        s, w, n = stats(build_v4n, base, prot, f"v4n {name}", size=size)
        S.append(s), W.append(w)
        N[f"v4n {name}"] = n
        s, w, _ = stats(build_v4n, base, prot, f"chain {name}", size=size, relevel=False)
        S.append(s), W.append(w)
    S, W = pd.concat(S, ignore_index=True), pd.concat(W, ignore_index=True)
    summ = (S.groupby(["estimate", "track"]).rmse_detrended_mm.median().rename("season_rmse_median_mm").to_frame()
            .join(W.assign(ae=W.difference_mm.abs()).groupby(["estimate", "track"]).ae.mean().rename("winter_mae_mm"))
            .reset_index())
    summ["pairs"] = [N.get(e, {}).get(t, np.nan) for e, t in zip(summ.estimate, summ.track)]
    cand = summ[summ.estimate.isin([f"v4n {b}" for b in BLOCKS])].groupby("estimate").season_rmse_median_mm.mean()
    lost = {b: max(1 - N[f"v4n {b}"][t] / N["v4"][t] for t in X.TRACKS) for b in BLOCKS}
    if any(v > 0.10 for v in lost.values()):
        raise SystemExit(f"stop rule: a block loses more than a tenth of v4's pairs: {lost}")
    b2, b3 = cand["v4n B2"], cand["v4n B3"]
    chosen = "B2" if abs(b3 - b2) <= 0.05 * min(b2, b3) or b2 < b3 else "B3"
    sha = subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True,
                         cwd=Path(__file__).resolve().parents[1]).stdout.strip()
    frozen = dict(base, method=("v4n = v4 (X-106) with the P6 phase read as the mean unwrapped phase of a fixed block "
                                "of cells around v4's read position (X-125, after X-124)"),
                  frozen_on=str(date.today()), code_sha=sha, base="v4_frozen_params.yaml (X-106), SPAMS copied",
                  read_block={"name": chosen, "size": BLOCKS[chosen], "min_share_finite": 0.5,
                              "grid": "cube as delivered (2026-10-15 snapshot), 40 m"})
    PARAMS.write_text(yaml.safe_dump(frozen, sort_keys=False, allow_unicode=True), encoding="utf-8")
    sel = pd.DataFrame([{"block": b, "season_rmse_median_mean_over_tracks_mm": float(cand[f"v4n {b}"]),
                         "pairs_lost_max_share": lost[b], "chosen": b == chosen} for b in BLOCKS])
    for nm, df in (("seasons", S), ("winters", W), ("summary", summ), ("selection", sel)):
        df.to_csv(OUT / f"{nm}.csv", index=False)
    md = lambda df: "\n".join(["| " + " | ".join(map(str, df.columns)) + " |", "|" + "---|" * len(df.columns)] +  # noqa: E731
                              ["| " + " | ".join(f"{x:.2f}" if isinstance(x, float) else str(x) for x in r) + " |"
                               for r in df.itertuples(index=False)])
    (OUT / "README.md").write_text(f"""# v4n frozen — selection and design-year report (generated, in-sample)

Generated by `experiments/freeze_v4n_x125.py`; do not edit. Rule: `../../preregistration/X-125.yaml`. Frozen
parameters: `../v4n_frozen_params.yaml` (code sha {sha}); v4's SPAMS copied from `../v4_frozen_params.yaml`. In-sample
on the design years, on the cube as delivered — for reference only; the scored comparison is X-107 (amendment_2).
"chain B2/B3" are the block chains without the winter re-levelling.

## Selection

{md(sel)}

## Summary (vertical mm against the laser)

{md(summ)}
""", encoding="utf-8")
    print(sel.to_string(index=False))
    print(summ.to_string(index=False))


if __name__ == "__main__":
    main()

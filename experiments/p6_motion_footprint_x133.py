"""X-133 — stage 3 of the supervisor's plan: how far over the mat do cells follow P6's six-day motion, does its
amplitude grow with depth, and a within-season motion map with its uncertainty, checked by the two tracks.

3a: every mat cell's consecutive six-day increments (hard-target reference) against the P6 laser's line-of-sight
increments, X-126's pairs of 2021 and 2025; a cell "follows P6" when its r exceeds the 95th percentile of solid-ground
cells read the same way; the share that follows, by spatial-block bootstrap; Spearman of the gain against depth.
3b: per season (May–October 2021, 2025) and track, the longest unbroken run of consecutive six-day pairs, each
following cell's vertical change over it with P6's per-pair error propagated; the two tracks compared cell by cell.

Pre-registration: hub ``08_deliverables/insar_full_test/preregistration/X-133.yaml`` (checked before the run). Stops
unless P6's own cell reproduces X-126. Design years only. Reads the cube and the hub; writes only to the hub
(``insar_full_test/footprint_x133/``).

    make run E=p6_motion_footprint_x133
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import landcover_coherence_x132 as L132  # noqa: E402  (X-132's readers: gold pairs at P6, silver stacks, X-126 pairs)
from insar_wetlands import contrasts as C  # noqa: E402
from insar_wetlands import footprint as F  # noqa: E402
from insar_wetlands.inversion.isbas import PHASE_TO_MM  # noqa: E402
from insar_wetlands.paths import cube_root, hub_root  # noqa: E402
from insar_wetlands.scoring import Protocol, assert_preregistered, guard_lockbox, increment_stats  # noqa: E402

HUB, CUBE = hub_root(), cube_root()
FT = HUB / "08_deliverables" / "insar_full_test"
OUT = FT / "footprint_x133"
PREREG = FT / "preregistration" / "X-133.yaml"
PROTOCOL = FT / "protocol.yaml"
X126_TABLE = FT / "x065_hard_targets" / "by_revisit.csv"
TRACKS = ("ascending", "descending")
PERIODS = {"2020_2021": "2021", "2025": "2025"}            # set → X-126 group
SEASON = (5, 10)                                            # May–October
DIST_BINS = [0, 250, 500, 1000, np.inf]
N_BOOT, N_PERM, SEED = 2000, 1000, 0
Q_NULL = 95


def load_cells() -> tuple[pd.DataFrame, tuple[int, int]]:
    S = xr.open_dataset(CUBE / "silver" / "static_40m.nc")
    pu = pd.read_csv(CUBE / "silver" / "plot_units.csv").set_index("plot")
    r6, c6 = int(pu.loc["P6", "row"]), int(pu.loc["P6", "col"])
    z, bw, ht, dm = (S[k].values for k in ("zone", "boardwalk_share", "hard_target", "dist_to_mat_m"))
    mat = (z == 1) & (bw < 0.2) & (ht == 0)
    null = np.isin(z, [3, 4]) & (dm >= 80) & (dm <= 1000) & (ht == 0) & (bw < 0.2)
    rr, cc = np.nonzero(mat | null)
    keys = ["zone", "spatial_block", "depth_in_mat_m", "dist_to_mat_m", "canopy2024_p90_m", "los_u_asc", "los_u_des"]
    df = pd.DataFrame({k: S[k].values[rr, cc] for k in keys})
    df["row"], df["col"], df["x"], df["y"] = rr, cc, S.x.values[cc], S.y.values[rr]
    df["role"] = np.where(mat[rr, cc], "mat", "null")
    df["dist_p6_m"] = np.hypot(df.x - S.x.values[c6], df.y - S.y.values[r6])
    df["cls"] = L132.classify(df.canopy2024_p90_m)
    df["cell"] = [f"r{r}c{c}" for r, c in zip(rr, cc)]
    return df.set_index("cell"), (r6, c6)


def increments(cells: pd.DataFrame, g: pd.DataFrame, track: str, period: str, pairs: pd.Series) -> np.ndarray:
    """(pairs × cells) increments in mm on the hard-target reference, for the given pairs of one set."""
    gp = g[(g.track == track) & (g.period == period)].set_index("pair").loc[pairs]
    st = L132.stacks(period, track, cells, pd.Series(pairs).reset_index(drop=True))
    return (st["unw"] - gp.ref_hard_unw.to_numpy(float)[:, None]) * PHASE_TO_MM


# ------------------------------------------------------------------------------ 3a

def footprint(cells: pd.DataFrame, g: pd.DataFrame, lp: pd.DataFrame, p6: str) -> tuple[pd.DataFrame, list, dict]:
    """Per cell × track (and per year): r, gain, permutation p; and the reproduction checks at P6."""
    rows, checks, perms = [], [], {}
    ref = pd.read_csv(X126_TABLE)
    j6 = cells.index.get_loc(p6)
    for track in TRACKS:
        inc_all, las_all, mon_all = [], [], []
        for period, grp in PERIODS.items():
            l = lp[(lp.track == track) & (lp.period.astype(str) == period)].sort_values("pair")
            inc = increments(cells, g, track, period, l.pair)
            a = l.d_laser_los_mm.to_numpy(float)
            mine = increment_stats(a, inc[:, j6])
            th = ref[(ref.track == track) & (ref.period.astype(str) == grp) & (ref.revisit_days == 6) &
                     (ref.pairs_kind == "consecutive") & (ref.reference == "hard_targets")].iloc[0]
            checks.append({"check": f"X-126 r {track} {grp}", "pairs": mine["pairs"], "pairs_x126": int(th.pairs),
                           "abs_diff_r": abs(mine["pearson_increments"] - th.pearson_increments),
                           "ok": bool(mine["pairs"] == th.pairs and
                                      abs(mine["pearson_increments"] - th.pearson_increments) < 1e-9)})
            rows.append(pd.DataFrame({"cell": cells.index, "track": track, "years": grp,
                                      "r": F.pearson_cols(a, inc), "gain": F.slope_origin_cols(a, inc),
                                      "pairs": np.isfinite(inc).sum(0)}))
            inc_all.append(inc)
            las_all.append(a)
            mon_all.append(100 * int(grp) + pd.to_datetime([p[:8] for p in l.pair], format="%Y%m%d").month.to_numpy())
        inc, a, mon = np.vstack(inc_all), np.concatenate(las_all), np.concatenate(mon_all)
        r = F.pearson_cols(a, inc)
        P = F.permuted_r(a, inc, mon, N_PERM, SEED)
        perms[track] = (1 + (P >= r).sum(0)) / (1 + N_PERM)
        rows.append(pd.DataFrame({"cell": cells.index, "track": track, "years": "2021+2025", "r": r,
                                  "gain": F.slope_origin_cols(a, inc), "pairs": np.isfinite(inc).sum(0),
                                  "perm_p": perms[track]}))
    return pd.concat(rows, ignore_index=True), checks, perms


def share_boot(follow: np.ndarray, blocks: np.ndarray) -> tuple[float, float, float]:
    ub, counts = C.block_counts(blocks, N_BOOT, SEED)
    W = C.unit_weights(blocks, ub, counts)
    f = follow.astype(float)
    draws = (W * f).sum(1) / W.sum(1)
    lo, hi = C.boot_interval(draws)
    return float(f.mean()), lo, hi


def spearman_boot(x: np.ndarray, y: np.ndarray, blocks: np.ndarray) -> tuple[float, float, float, int]:
    ok = np.isfinite(x) & np.isfinite(y)
    x, y, blocks = x[ok], y[ok], blocks[ok]
    if len(x) < 10:
        return np.nan, np.nan, np.nan, len(x)
    ub, counts = C.block_counts(blocks, N_BOOT, SEED)
    pos = np.searchsorted(ub, blocks)
    draws = []
    for c in counts:
        i = F.resample_index(c, pos)
        if len(np.unique(i)) > 2:
            draws.append(spearmanr(x[i], y[i]).statistic)
    lo, hi = C.boot_interval(draws)
    return float(spearmanr(x, y).statistic), lo, hi, len(x)


def read_3a(cells: pd.DataFrame, fp: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    prim, sec, thr = [], [], {}
    for track in TRACKS:
        for years in ("2021+2025", "2021", "2025"):
            d = fp[(fp.track == track) & (fp.years == years)].set_index("cell").join(cells)
            nul = d[d.role == "null"].r.dropna()
            mat = d[(d.role == "mat")]
            for null_name, nn in (("solid ground", nul), ("grassland", d[(d.role == "null") & (d.zone == 3)].r.dropna())):
                t = float(np.percentile(nn, Q_NULL))
                ok = mat.r.notna()
                follow = (mat.r > t).to_numpy()[ok.to_numpy()]
                sh, lo, hi = share_boot(follow, mat.spatial_block.to_numpy()[ok.to_numpy()])
                row = {"track": track, "years": years, "null": null_name, "null_cells": len(nn), "threshold_r": t,
                       "mat_cells": int(ok.sum()), "following": int(follow.sum()), "share": sh, "lo95": lo, "hi95": hi}
                if years == "2021+2025" and null_name == "solid ground":
                    thr[track] = t
                    fol = mat[mat.r > t]
                    rho, rlo, rhi, n = spearman_boot(fol.depth_in_mat_m.to_numpy(float), fol.gain.to_numpy(float),
                                                     fol.spatial_block.to_numpy())
                    row.update({"gain_depth_rho": rho, "rho_lo95": rlo, "rho_hi95": rhi, "rho_cells": n,
                                "gain_median": float(fol.gain.median())})
                    prim.append(row)
                    # secondaries: by distance to P6 and by land-cover class; permutation null
                    for lab, m in (*[(f"dist {int(a)}–{b if np.isfinite(b) else '…'} m",
                                      (mat.dist_p6_m >= a) & (mat.dist_p6_m < b)) for a, b in zip(DIST_BINS[:-1], DIST_BINS[1:])],
                                   *[(f"class {k}", mat.cls == k) for k in ("open", "mid", "wooded")]):
                        s = mat[m & ok]
                        if len(s):
                            f2 = (s.r > t).to_numpy()
                            a2 = share_boot(f2, s.spatial_block.to_numpy())
                            sec.append({"track": track, "subset": lab, "cells": len(s), "share": a2[0],
                                        "lo95": a2[1], "hi95": a2[2]})
                    pp = mat.perm_p.dropna()
                    sec.append({"track": track, "subset": "permutation null within month: p < 0.05",
                                "cells": len(pp), "share": float((pp < 0.05).mean()), "lo95": np.nan, "hi95": np.nan})
                else:
                    sec.append({"track": track, "subset": f"{years}, null {null_name}", "cells": int(ok.sum()),
                                "share": sh, "lo95": lo, "hi95": hi, "threshold_r": t})
    return pd.DataFrame(prim), pd.DataFrame(sec), thr


def rules_3a(prim: pd.DataFrame) -> pd.DataFrame:
    out = []
    for r in prim.itertuples():
        out.append({"track": r.track,
                    "P6's motion reaches beyond P6": r.lo95 > 0.05,
                    "most of the mat follows": r.lo95 > 0.5})
    t = pd.DataFrame(out)
    up = all(prim.rho_lo95 > 0)
    down = all(prim.rho_hi95 < 0)
    t["motion with depth (both tracks)"] = "grows" if up else ("shrinks" if down else "not shown")
    return t


# ------------------------------------------------------------------------------ 3b

def season_map(cells: pd.DataFrame, g: pd.DataFrame, following: dict, p6: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    ref = pd.read_csv(X126_TABLE)
    lp_all = pd.read_csv(CUBE / "gold" / "plots_pair.csv.gz", usecols=["plot", "grid", "track", "pair", "d_laser_los_mm"])
    lp_all = lp_all[(lp_all["plot"] == "P6") & (lp_all.grid == "40m")].drop_duplicates(["track", "pair"]).set_index(
        ["track", "pair"]).d_laser_los_mm
    rows, runs = [], []
    for period, yr in PERIODS.items():
        for track in TRACKS:
            gp = g[(g.track == track) & (g.period == period) & (g.between == 0) & (g.revisit_days == 6)]
            gp = gp[(gp.t1.dt.year == int(yr)) & gp.t1.dt.month.between(*SEASON) & gp.t2.dt.month.between(*SEASON) &
                    (gp.t2.dt.year == int(yr))]
            if yr == "2021":
                gp = gp[gp.t1 >= pd.Timestamp("2021-05-27", tz="UTC")]
            run = gp.iloc[F.longest_run(gp.t1.values, gp.t2.values)]
            inc = increments(cells, g, track, period, run.pair)
            lu = cells["los_u_asc" if track == "ascending" else "los_u_des"].to_numpy(float)
            complete = np.isfinite(inc).all(0)
            v = np.where(complete, inc.sum(0), np.nan) / lu
            th = ref[(ref.track == track) & (ref.period.astype(str) == yr) & (ref.revisit_days == 6) &
                     (ref.pairs_kind == "consecutive") & (ref.reference == "hard_targets")].iloc[0]
            sig = th.rmse_increments_mm * np.sqrt(len(run)) / np.abs(lu)
            rows.append(pd.DataFrame({"cell": cells.index, "year": yr, "track": track, "vertical_change_mm": v,
                                      "sigma_mm": sig, "following": cells.index.isin(following[track])}))
            j6 = cells.index.get_loc(p6)
            lz = lp_all.reindex(pd.MultiIndex.from_arrays([[track] * len(run), run.pair])).to_numpy(float)
            has = np.isfinite(lz)                      # the run's pairs with a laser change (gaps, laser start)
            runs.append({"year": yr, "track": track, "pairs": len(run), "start": run.t1.min().date(),
                         "end": run.t2.max().date(), "p6_vertical_mm": v[j6], "p6_sigma_mm": sig[j6],
                         "pairs_with_laser": int(has.sum()),
                         "p6_vertical_on_laser_pairs_mm": float(inc[has, j6].sum() / lu[j6]),
                         "laser_vertical_on_laser_pairs_mm": float(lz[has].sum() / lu[j6])})
    return pd.concat(rows, ignore_index=True), pd.DataFrame(runs)


def tracks_agree(cells: pd.DataFrame, mp: pd.DataFrame, following: dict) -> pd.DataFrame:
    both = set(following["ascending"]) & set(following["descending"])
    out = []
    for yr in PERIODS.values():
        w = mp[mp.year == yr].pivot_table(index="cell", columns="track", values="vertical_change_mm")
        w = w[w.index.isin(both)].dropna()
        b = cells.loc[w.index, "spatial_block"].to_numpy()
        ub, counts = C.block_counts(b, N_BOOT, SEED)
        pos = np.searchsorted(ub, b)
        a, d = w.ascending.to_numpy(), w.descending.to_numpy()
        draws = []
        for c in counts:
            i = F.resample_index(c, pos)
            if len(i) > 3 and np.std(a[i]) > 0 and np.std(d[i]) > 0:
                draws.append(np.corrcoef(a[i], d[i])[0, 1])
        lo, hi = C.boot_interval(draws)
        r = float(np.corrcoef(a, d)[0, 1]) if len(a) > 3 else np.nan
        out.append({"year": yr, "cells": len(w), "r_asc_desc": r, "lo95": lo, "hi95": hi,
                    "the two tracks draw the same pattern": bool(lo > 0),
                    "median_asc_mm": float(np.median(a)) if len(a) else np.nan,
                    "median_desc_mm": float(np.median(d)) if len(d) else np.nan})
    return pd.DataFrame(out)


# ------------------------------------------------------------------------------ figures

def figs(cells: pd.DataFrame, fp: pd.DataFrame, thr: dict, mp: pd.DataFrame, out: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    mat = cells[cells.role == "mat"]
    p6 = cells[(cells.row == P6[0]) & (cells["col"] == P6[1])]
    fig, axs = plt.subplots(1, 4, figsize=(20, 5))
    for i, track in enumerate(TRACKS):
        d = fp[(fp.track == track) & (fp.years == "2021+2025")].set_index("cell").reindex(cells.index)
        ax = axs[i]
        sc = ax.scatter(cells.x, cells.y, c=d.r, s=6, marker="s", cmap="RdBu_r", vmin=-1, vmax=1)
        f = mat.index[d.loc[mat.index].r > thr[track]]
        ax.scatter(cells.loc[f, "x"], cells.loc[f, "y"], s=2, c="k")
        ax.scatter(p6.x, p6.y, marker="*", s=120, c="w", edgecolors="k")
        ax.set_title(f"r with the P6 laser, six-day pairs 2021+2025 — {track[:3]}\nblack: mat cells above the "
                     f"solid-ground 95th percentile ({thr[track]:.2f})", fontsize=8)
        ax.set_aspect("equal")
        fig.colorbar(sc, ax=ax, shrink=0.7)
        ax = axs[2 + i]
        m = d.loc[mat.index]
        fol = m.r > thr[track]
        ax.scatter(mat.depth_in_mat_m[~fol], m.gain[~fol], s=6, c="0.7", label="mat, not following")
        ax.scatter(mat.depth_in_mat_m[fol], m.gain[fol], s=8, c="C0", label="mat, following")
        ax.axhline(1, c="k", lw=0.6, ls=":")
        ax.set_xlabel("depth into the mat (m)")
        ax.set_ylabel("gain against the P6 laser (slope through origin)")
        ax.set_title(f"Gain against depth — {track[:3]}", fontsize=8)
        ax.legend(fontsize=7)
    fig.savefig(out / "fig_footprint.png", dpi=110, bbox_inches="tight")
    plt.close(fig)
    fig, axs = plt.subplots(2, 4, figsize=(20, 9))
    for i, yr in enumerate(PERIODS.values()):
        for j, track in enumerate(TRACKS):
            d = mp[(mp.year == yr) & (mp.track == track)].set_index("cell").loc[mat.index]
            v = d.vertical_change_mm.where(d.following)
            q = np.nanpercentile(np.abs(v), 95) if np.isfinite(v).any() else 1
            ax = axs[i, j]
            sc = ax.scatter(mat.x, mat.y, c=v, s=10, marker="s", cmap="RdBu_r", vmin=-q, vmax=q)
            ax.scatter(p6.x, p6.y, marker="*", s=120, c="w", edgecolors="k")
            ax.set_title(f"{yr} May–Oct vertical change (mm), following cells — {track[:3]}\n"
                         f"1σ = {float(d.sigma_mm.median()):.0f} mm (median)", fontsize=8)
            ax.set_aspect("equal")
            fig.colorbar(sc, ax=ax, shrink=0.7)
        w = mp[mp.year == yr].pivot_table(index="cell", columns="track", values="vertical_change_mm")
        both = mp[(mp.year == yr)].groupby("cell").following.all()
        w = w[both.reindex(w.index).fillna(False)].dropna()
        ax = axs[i, 2]
        ax.scatter(w.ascending, w.descending, s=6)
        lim = np.nanpercentile(np.abs(w.values), 99) if len(w) else 1
        ax.plot([-lim, lim], [-lim, lim], "k:", lw=0.6)
        ax.set_xlabel("ascending (mm)")
        ax.set_ylabel("descending (mm)")
        ax.set_title(f"{yr}: the two tracks, cells following on both", fontsize=8)
        axs[i, 3].axis("off")
    fig.savefig(out / "fig_season_maps.png", dpi=110, bbox_inches="tight")
    plt.close(fig)


# ------------------------------------------------------------------------------ main

P6 = (0, 0)


def main() -> None:
    global P6
    OUT.mkdir(parents=True, exist_ok=True)
    names = ["cells.csv", "footprint.csv.gz", "share.csv", "secondary.csv", "rules.csv", "season_map.csv.gz",
             "season_runs.csv", "tracks_agree.csv", "reproduction.csv", "README.md", "fig_footprint.png",
             "fig_season_maps.png"]
    assert_preregistered(PREREG, [OUT / f for f in names])
    cells, P6 = load_cells()
    p6 = f"r{P6[0]}c{P6[1]}"
    g = L132.gold_p6()
    g = g[g.period.isin(list(PERIODS))]
    lp = L132.laser_pairs_x126()
    guard_lockbox(pd.concat([g.t1, g.t2]), Protocol.load(PROTOCOL))
    fp, checks, _ = footprint(cells, g, lp, p6)
    rep = pd.DataFrame(checks)
    n_null = {t: int(fp[(fp.track == t) & (fp.years == "2021+2025")].set_index("cell").join(cells)
                     .query("role == 'null'").r.notna().sum()) for t in TRACKS}
    for t, n in n_null.items():
        rep.loc[len(rep)] = {"check": f"null cells {t}", "pairs": n, "pairs_x126": 200, "abs_diff_r": np.nan, "ok": n >= 200}
    rep.to_csv(OUT / "reproduction.csv", index=False)
    if not rep.ok.all():
        print(rep.to_string())
        sys.exit("stop (pre-registered): see reproduction.csv")
    cells.to_csv(OUT / "cells.csv")
    fp.to_csv(OUT / "footprint.csv.gz", index=False)
    prim, sec, thr = read_3a(cells, fp)
    prim.to_csv(OUT / "share.csv", index=False)
    sec.to_csv(OUT / "secondary.csv", index=False)
    rl = rules_3a(prim)
    mat = cells[cells.role == "mat"]
    following = {}
    for t in TRACKS:
        d = fp[(fp.track == t) & (fp.years == "2021+2025")].set_index("cell")
        following[t] = list(mat.index[d.loc[mat.index].r > thr[t]])
    mp, runs = season_map(cells[cells.role == "mat"], g, following, p6)
    mp.to_csv(OUT / "season_map.csv.gz", index=False)
    runs.to_csv(OUT / "season_runs.csv", index=False)
    ag = tracks_agree(cells, mp, following)
    ag.to_csv(OUT / "tracks_agree.csv", index=False)
    rl.to_csv(OUT / "rules.csv", index=False)
    figs(cells, fp, thr, mp, OUT)
    md = L132.md
    (OUT / "README.md").write_text(f"""# X-133 — how far P6's six-day motion reaches over the mat, and a season map (generated)

Generated by `experiments/p6_motion_footprint_x133.py` (`make run E=p6_motion_footprint_x133`); do not edit.
Pre-registration: `../preregistration/X-133.yaml`. Mat cells: zone A (boardwalk < 0.2); null: solid ground (zones C, D)
80–1000 m from the mat edge. X-126's consecutive six-day pairs of 2021 (from 27 May) and 2025, hard-target reference,
the cube's P6 laser (snow-evidence rule). {N_BOOT} resamples of the mat's 200 m blocks, seed {SEED}; permutation null
{N_PERM} draws within calendar month. Design years only. One truth point: "follows" means a cell's increments follow
P6's timing; its own amplitude is the gain.

## Reading rules (pre-registered)

{md(rl)}

{md(ag[["year", "the two tracks draw the same pattern"]])}

## Reproduction (P6's own cell against X-126) and null size

{md(rep, 12)}

## 3a — the share of mat cells that follow P6, and the gain against depth (following cells)

{md(prim)}

## 3a — secondary: by year, grassland-only null, distance to P6, land-cover class, permutation null

{md(sec)}

## 3b — season runs (longest unbroken chain of six-day pairs, May–October) and P6 against the laser

{md(runs)}

## 3b — the two tracks, cell by cell (cells following on both)

{md(ag)}

## Figures

- `fig_footprint.png` — r with the P6 laser per cell (cells above the null marked) and the gain against depth.
- `fig_season_maps.png` — vertical change over each season's run, following cells, both tracks; the tracks against
  each other.

## Scoring deferred to the UAV (fixed in the pre-registration)

When the UAV rasters arrive (C-065), every pair of flights within one season is scored against the map rebuilt over the
same dates, before the two are looked at together.
""", encoding="utf-8")
    print(rl.to_string())
    print(prim.to_string())
    print(ag.to_string())
    print(runs.to_string())


if __name__ == "__main__":
    main()

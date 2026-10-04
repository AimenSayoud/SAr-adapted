"""X-076 — fusion v3 WP3: does the mat move differently in different places? (exploratory, Mac-safe, reads the cube)

E3.1 A physical prediction instead of free EOFs. If the mat is held at its shore and floats inside (an edge hinge),
     the interior moves more than the edge, in proportion to how much the whole raft moves on that date. Per pair, the
     mat's departures from its own median are regressed on physical shapes (depth in the mat; north–south and
     east–west position; LiDAR micro-relief; vegetation height). Prediction: the depth coefficient a_k follows the raft
     motion c_k (the mat median, hard-target reference) with a positive slope. Null: the same statistic with the pairs
     permuted, and the same regression on compact stable-ground patches of the same pixel count (their own "depth")
     against the same raft motion.
E3.2 Where and when the mat is measurable: the WP2 detectability model applied to every mat pixel (its coherence,
     the pair's revisit, the driver's cycle risk, season), mapped per season and track.
E3.3 The 2026 20 m grid: E3.1 at 20 m (finer shapes, P8 and P9 in separate cells).

    PYTHONPATH=src python scripts/fusion_v3_wp3_x076.py
Outputs only to the hub (detectability uses the water table): ``08_deliverables/fusion_v3/wp3/``.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import xarray as xr  # noqa: E402
from scipy import ndimage, stats  # noqa: E402
from sklearn.linear_model import LogisticRegression  # noqa: E402
from sklearn.pipeline import make_pipeline  # noqa: E402
from sklearn.preprocessing import StandardScaler  # noqa: E402

from insar_wetlands.cube.core import rad_to_los_mm  # noqa: E402

HUB = Path(__file__).resolve().parents[3]
G, S = HUB / "06_data" / "cube" / "gold", HUB / "06_data" / "cube" / "silver"
OUT = HUB / "08_deliverables" / "fusion_v3" / "wp3"
WP2 = HUB / "08_deliverables" / "fusion_v3" / "wp2"
RNG = np.random.default_rng(76)
RUNS = [("2020_2021", 6, ""), ("2025", 6, ""), ("2022_2024", 12, ""), ("2026", 6, "_20m")]
SHAPES = ["depth_in_mat", "north_south", "east_west", "microrelief", "vegetation_height", "pixel_coherence"]


def zscore(a):
    a = np.asarray(a, float)
    return (a - np.nanmean(a)) / np.nanstd(a)


def shapes_for(rows, cols, depth, micro, vegh, coh) -> np.ndarray:
    """Standardised shapes; the pixel's median coherence is a control (filtering damps phase where coherence is low)."""
    S_ = np.column_stack([zscore(depth), zscore(-rows), zscore(cols), zscore(micro), zscore(vegh), zscore(coh)])
    return np.where(np.isfinite(S_), S_, 0.0)


def coefficients(Dep: np.ndarray, Sh: np.ndarray) -> np.ndarray:
    """Per pair (rows of Dep, pixels in columns), least-squares coefficients of the departures on the shapes."""
    A = np.column_stack([Sh, np.ones(len(Sh))])
    out = np.full((Dep.shape[0], Sh.shape[1]), np.nan)
    for k, d in enumerate(Dep):
        ok = np.isfinite(d)
        if ok.sum() > 3 * A.shape[1]:
            out[k] = np.linalg.lstsq(A[ok], d[ok], rcond=None)[0][:-1]
    return out


def slope_test(a: np.ndarray, c: np.ndarray, n_perm: int = 2000) -> dict:
    ok = np.isfinite(a) & np.isfinite(c)
    a, c = a[ok], c[ok]
    if len(a) < 10:
        return {"n_pairs": len(a), "slope": np.nan, "r": np.nan, "perm_p": np.nan}
    b = np.polyfit(c, a, 1)[0]
    r = stats.pearsonr(c, a)[0]
    null = np.array([stats.pearsonr(RNG.permutation(c), a)[0] for _ in range(n_perm)])
    return {"n_pairs": int(len(a)), "slope": float(b), "r": float(r),
            "perm_p": float((1 + np.sum(null >= r)) / (1 + n_perm))}


def e31(n_patches: int = 60) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows, nulls = [], []
    for period, rv, sfx in RUNS:
        for track in ("ascending", "descending"):
            f = G / f"mat_pairs_{period}_{track}{sfx}.nc"
            if not f.exists():
                continue
            ds = xr.open_dataset(f)
            keep = (ds.revisit_days.values == rv) & np.array(["S1D" not in p for p in ds.platforms.values.astype(str)])
            U = ds.unw.values[:, keep].T.astype(float)                # pairs × pixels
            zone = ds.st_zone.values
            hard = ds.st_hard_target.values == 1 if "st_hard_target" in ds else np.zeros(len(zone), bool)
            stable = ds.st_stable.values == 1
            ref = np.nanmedian(U[:, hard], 1) if hard.any() else np.nanmedian(U[:, stable], 1)
            mat = zone == 1
            Um = -rad_to_los_mm(U[:, mat] - ref[:, None])
            raft = np.nanmedian(Um, 1)                                # c_k: the whole mat's motion (mm, LOS)
            Dep = Um - raft[:, None]
            rr, cc = ds.row.values[mat], ds.col.values[mat]
            if "st_depth_in_mat_m" in ds:
                depth = ds.st_depth_in_mat_m.values[mat]
            else:                                                     # 20 m grid: distance to the mat's edge here
                res_m = 20.0 if sfx else 40.0
                gm = np.zeros((ds.row.values.max() + 1, ds.col.values.max() + 1), bool)
                gm[rr, cc] = True
                depth = (ndimage.distance_transform_edt(gm) * res_m)[rr, cc]
            cohm = np.nanmedian(ds.coh.values[:, keep], 1)                 # each pixel's median coherence over the pairs
            Sh = shapes_for(rr, cc, depth, ds.st_lidar_microrelief_sd_m.values[mat], ds.st_lidar_vegh_mean_m.values[mat],
                            cohm[mat])
            A = coefficients(Dep, Sh)
            for j, name in enumerate(SHAPES):
                rows.append({"period": period, "track": track, "grid": "20m" if sfx else "40m", "revisit_days": rv,
                             "shape": name, "subset": "all mat pixels", "mat_pixels": int(mat.sum()),
                             **slope_test(A[:, j], raft)})
            # mixed pixels at the shore would also move less: repeat on interior pixels only (≥ 80 m from the edge)
            inner = depth >= 80
            if inner.sum() > 40:
                Ai = coefficients(Dep[:, inner], shapes_for(rr[inner], cc[inner], depth[inner],
                                                            ds.st_lidar_microrelief_sd_m.values[mat][inner],
                                                            ds.st_lidar_vegh_mean_m.values[mat][inner], cohm[mat][inner]))
                rows.append({"period": period, "track": track, "grid": "20m" if sfx else "40m", "revisit_days": rv,
                             "shape": "depth_in_mat", "subset": "interior only (≥ 80 m)", "mat_pixels": int(inner.sum()),
                             **slope_test(Ai[:, 0], raft)})
            # null: compact stable patches of the same size, their own depth / position, against the same raft motion
            st_idx = np.flatnonzero(stable & ~hard)
            R, C = ds.row.values, ds.col.values
            for k in range(n_patches):
                c0 = st_idx[RNG.integers(len(st_idx))]
                order = np.argsort((R[st_idx] - R[c0]) ** 2 + (C[st_idx] - C[c0]) ** 2)[:int(mat.sum())]
                pidx = st_idx[order]
                grid = np.zeros((R.max() + 1, C.max() + 1), bool)
                grid[R[pidx], C[pidx]] = True
                dmap = ndimage.distance_transform_edt(grid)
                Up = -rad_to_los_mm(U[:, pidx] - ref[:, None])
                Dp = Up - np.nanmedian(Up, 1)[:, None]
                Shp = shapes_for(R[pidx], C[pidx], dmap[R[pidx], C[pidx]],
                                 ds.st_lidar_microrelief_sd_m.values[pidx], ds.st_lidar_vegh_mean_m.values[pidx], cohm[pidx])
                Ap = coefficients(Dp, Shp)
                ok = np.isfinite(Ap[:, 0]) & np.isfinite(raft)
                res_m = 20.0 if sfx else 40.0
                pd_m = dmap[R[pidx], C[pidx]] * res_m
                inn = pd_m >= 80
                ri = np.nan
                if inn.sum() > 40:
                    Ai = coefficients(Dp[:, inn], shapes_for(R[pidx][inn], C[pidx][inn], pd_m[inn],
                                                             ds.st_lidar_microrelief_sd_m.values[pidx][inn],
                                                             ds.st_lidar_vegh_mean_m.values[pidx][inn], cohm[pidx][inn]))
                    oki = np.isfinite(Ai[:, 0]) & np.isfinite(raft)
                    if oki.sum() > 10:
                        ri = float(stats.pearsonr(raft[oki], Ai[oki, 0])[0])
                if ok.sum() > 10:
                    nulls.append({"period": period, "track": track, "grid": "20m" if sfx else "40m", "patch": k,
                                  "r_depth_vs_raft": float(stats.pearsonr(raft[ok], Ap[ok, 0])[0]),
                                  "r_depth_vs_raft_interior": ri})
    res, nl = pd.DataFrame(rows), pd.DataFrame(nulls)
    if len(nl):
        q_all = nl.groupby(["period", "track", "grid"]).r_depth_vs_raft.agg(
            null_median="median", null_p95=lambda v: np.percentile(v, 95)).reset_index()
        q_all["subset"] = "all mat pixels"
        q_in = nl.groupby(["period", "track", "grid"]).r_depth_vs_raft_interior.agg(
            null_median="median", null_p95=lambda v: np.nanpercentile(v, 95) if np.isfinite(v).any() else np.nan).reset_index()
        q_in["subset"] = "interior only (≥ 80 m)"
        q = pd.concat([q_all, q_in], ignore_index=True)
        q["shape"] = "depth_in_mat"
        res = res.merge(q, on=["period", "track", "grid", "subset", "shape"], how="left")
        res["above_patch_null"] = np.where(res.null_p95.notna(), (res.r > res.null_p95).astype(float), np.nan)
    return res, nl


def e32() -> pd.DataFrame:
    """Detectability per mat pixel: the WP2 logistic model refitted on all P6 pairs, applied to each pixel's pairs."""
    q = pd.read_csv(WP2 / "e22_pairs.csv", parse_dates=["t1", "t2"])
    feats = ["coh_3x3", "revisit_days", "abs_driver_dh_mm", "wet_ends", "is_summer", "abs_d_vh_db", "is_desc"]
    X = q[feats].to_numpy(float)
    X = np.where(np.isfinite(X), X, np.nanmedian(X, 0))
    mdl = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000)).fit(X, q.match)
    site = q.set_index(["track", "pair"])[["abs_driver_dh_mm", "wet_ends", "is_summer"]]
    rows = []
    for period in ("2020_2021", "2025", "2022_2024"):
        for track in ("ascending", "descending"):
            f = G / f"mat_pairs_{period}_{track}.nc"
            if not f.exists():
                continue
            ds = xr.open_dataset(f)
            mat = ds.st_zone.values == 1
            pairs = ds.pair.values.astype(str)
            rv = ds.revisit_days.values
            keep = rv <= 12
            drv = site.reindex(pd.MultiIndex.from_arrays([[track] * len(pairs), pairs]))
            mid = pd.DatetimeIndex(ds.t1.values) + (pd.DatetimeIndex(ds.t2.values) - pd.DatetimeIndex(ds.t1.values)) / 2
            summer = mid.month.isin([6, 7, 8]).astype(int)
            coh = ds.coh.values[mat][:, keep]                       # pixels × pairs
            for k in np.flatnonzero(keep):
                kk = np.flatnonzero(keep).tolist().index(k)
                Xp = np.column_stack([coh[:, kk], np.full(mat.sum(), rv[k]),
                                      np.full(mat.sum(), drv.abs_driver_dh_mm.iloc[k]),
                                      np.full(mat.sum(), drv.wet_ends.iloc[k]), np.full(mat.sum(), summer[k]),
                                      np.full(mat.sum(), np.nan), np.full(mat.sum(), int(track == "descending"))])
                Xp = np.where(np.isfinite(Xp), Xp, np.nanmedian(X, 0))
                pr = mdl.predict_proba(Xp)[:, 1]
                rows.append(pd.DataFrame({"period": period, "track": track, "revisit_days": rv[k],
                                          "season": "JJA" if summer[k] else ("SON" if mid[k].month in (9, 10, 11) else
                                                                              "MAM" if mid[k].month in (3, 4, 5) else "DJF"),
                                          "row": ds.row.values[mat], "col": ds.col.values[mat], "p_detect": pr}))
    return pd.concat(rows, ignore_index=True)


def md(df: pd.DataFrame, nd: int = 2) -> str:
    dd = df.copy()
    for c in dd.columns:
        if dd[c].dtype.kind == "f":
            dd[c] = dd[c].round(nd)
    lines = ["| " + " | ".join(map(str, dd.columns)) + " |", "|" + "---|" * len(dd.columns)]
    for r in dd.itertuples(index=False):
        lines.append("| " + " | ".join("" if (isinstance(v, float) and np.isnan(v)) else str(v) for v in r) + " |")
    return "\n".join(lines)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    res, nl = e31()
    res.to_csv(OUT / "e31_shape_vs_raft.csv", index=False)
    nl.to_csv(OUT / "e31_patch_null.csv", index=False)
    det = e32()
    summ = (det.groupby(["track", "revisit_days", "season"]).p_detect
            .agg(median="median", share_above_half=lambda v: float(np.mean(v > 0.5)), n="size").reset_index())
    summ.to_csv(OUT / "e32_detectability_summary.csv", index=False)
    pix = det[det.revisit_days == 6].groupby(["track", "season", "row", "col"]).p_detect.mean().reset_index()
    pix.to_csv(OUT / "e32_detectability_per_pixel_6day.csv", index=False)
    seasons = ["MAM", "JJA", "SON"]
    fig, ax = plt.subplots(2, 3, figsize=(12, 7))
    for i, track in enumerate(("ascending", "descending")):
        for j, s_ in enumerate(seasons):
            g = pix[(pix.track == track) & (pix.season == s_)]
            a = ax[i, j]
            if len(g):
                img = np.full((g.row.max() - g.row.min() + 1, g.col.max() - g.col.min() + 1), np.nan)
                img[g.row - g.row.min(), g.col - g.col.min()] = g.p_detect
                im = a.imshow(img, vmin=0, vmax=1, cmap="viridis")
            a.set_title(f"{track}, {s_}: P(6-day pair carries the motion)", fontsize=8)
            a.set_xticks([])
            a.set_yticks([])
    fig.colorbar(im, ax=ax, shrink=0.6)
    fig.savefig(OUT / "fig_detectability_maps.png", dpi=140)
    plt.close(fig)
    depth = res[res["shape"] == "depth_in_mat"]
    txt = ["# X-076 — fusion v3 WP3: does the mat move differently in different places? (exploratory)", "",
           "Generated by `05_code/SAr-adapted/scripts/fusion_v3_wp3_x076.py` from `06_data/cube`; every number is computed.", "",
           "## E3.1 Edge hinge: does the interior move more when the raft moves?", "",
           "Per pair, the mat's departures from its own median regressed on physical shapes; the coefficient against "
           "the raft motion (mat median, hard-target reference). Permutation p over pairs (floor 1/2001); patch null = "
           "the depth statistic on compact stable-ground patches of the mat's pixel count.", "",
           md(depth[["period", "track", "grid", "revisit_days", "subset", "mat_pixels", "n_pairs", "slope", "r", "perm_p",
                     "null_median", "null_p95", "above_patch_null"]]), "",
           "Other shapes (same test, no patch null):", "",
           md(res[res["shape"] != "depth_in_mat"][["period", "track", "grid", "shape", "n_pairs", "slope", "r", "perm_p"]]), "",
           "Read: an edge hinge predicts a positive depth slope that beats both the permutation and the patch null. "
           "A slope the stable patches reproduce comes from the processing (unwrapping, filtering), not the mat. Mixed "
           "pixels at the shore (part floating mat, part fixed ground) also move less: the interior-only rows exclude "
           "them; a slope that survives there is a gradient inside the floating mat. Every regression also holds the "
           "pixel's median coherence fixed (filtering damps phase where coherence is low).", "",
           "## E3.2 Where and when the mat is measurable", "",
           "The WP2 detectability model (logistic, all P6 pairs) applied to every mat pixel and pair ≤ 12 days:", "",
           md(summ), "", "![detectability maps](fig_detectability_maps.png)", ""]
    (OUT / "README.md").write_text("\n".join(txt))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()

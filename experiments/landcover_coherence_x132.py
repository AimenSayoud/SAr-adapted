"""X-132 — around P6, trees against open ground on and off the floating mat: how land cover, depth, water and weather
drive coherence, backscatter and the six-day phase (the supervisor's plan, stages 1–2).

Cells: 40 m cells of the corrected cube within 1 km of P6's cell, boardwalk share < 0.2, no hard targets, no open water.
Wooded / open from the dated 2024 laser scan (C-072). Stage 1: a 2×2 (wooded / open × on / off the mat), each wooded
cell matched to an open one on depth into the mat (on it) or distance to its edge (off it), then on distance to P6;
outcomes O1–O5 per cell; contrasts with 95 % intervals from 2000 resamples of the 200 m spatial blocks (seed 0).
Stage 2: Fisher-z coherence of every cell × consecutive pair regressed on every factor; standardised effects and the
share of explained variance each carries, by the same block bootstrap.

Pre-registration: hub ``08_deliverables/insar_full_test/preregistration/X-132.yaml`` (checked before the run). Stops
unless P6's own cell reproduces the cube's gold coherence and X-126's consecutive six-day r of 2021 and 2025 (1e-9).
Design years only. Reads the cube and the hub; writes only to the hub (``insar_full_test/landcover_x132/``).

    make run E=landcover_coherence_x132
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import x065_hard_targets_x126 as X126  # noqa: E402  (X-126's pair list, consecutive rule, laser pairs)
from insar_wetlands import contrasts as C  # noqa: E402
from insar_wetlands.inversion.isbas import PHASE_TO_MM  # noqa: E402
from insar_wetlands.paths import cube_root, hub_root  # noqa: E402
from insar_wetlands.scoring import Protocol, assert_preregistered, guard_lockbox, increment_stats  # noqa: E402

HUB, CUBE = hub_root(), cube_root()
FT = HUB / "08_deliverables" / "insar_full_test"
OUT = FT / "landcover_x132"
PREREG = FT / "preregistration" / "X-132.yaml"
PROTOCOL = FT / "protocol.yaml"
X126_TABLE = FT / "x065_hard_targets" / "by_revisit.csv"
TRACKS = ("ascending", "descending")
SETS = {"2020_2021": 6, "2022_2024": 12, "2025": 6}           # set → revisit of its consecutive pairs
LAGS = (6, 12, 24, 48)
RADIUS_M, CALIPER_M = 1000.0, 40.0
N_BOOT, SEED = 2000, 0
SUMMER, WINTER = (5, 6, 7, 8, 9), (12, 1, 2)
OUTCOMES_PRIMARY = "O1_coh"


# ------------------------------------------------------------------------------ cells and classes

def classify(v, hi: float = 2.0, lo: float = 0.5) -> np.ndarray:
    v = np.asarray(v, float)
    return np.where(v >= hi, "wooded", np.where(v < lo, "open", np.where(np.isfinite(v), "mid", "none")))


def load_cells(radius: float = RADIUS_M) -> tuple[pd.DataFrame, tuple[int, int]]:
    S = xr.open_dataset(CUBE / "silver" / "static_40m.nc")
    pu = pd.read_csv(CUBE / "silver" / "plot_units.csv").set_index("plot")
    r6, c6 = int(pu.loc["P6", "row"]), int(pu.loc["P6", "col"])
    X, Y = np.meshgrid(S.x.values, S.y.values)
    d = np.hypot(X - S.x.values[c6], Y - S.y.values[r6])
    rr, cc = np.nonzero(d <= radius)
    keys = ["zone", "spatial_block", "depth_in_mat_m", "dist_to_mat_m", "boardwalk_share", "hard_target", "wc_tree",
            "lidar_vegh_p90_m", "canopy2024_p90_m", "canopy2024_cover", "canopy2011_p90_m",
            "cir_ndvi_2021_mean", "cir_ndvi_2022_mean", "cir_ndvi_2023_mean"]
    df = pd.DataFrame({k: S[k].values[rr, cc] for k in keys})
    df["row"], df["col"], df["x"], df["y"], df["dist_p6_m"] = rr, cc, S.x.values[cc], S.y.values[rr], d[rr, cc]
    df = df[(df.boardwalk_share < 0.2) & (df.hard_target == 0) & df.zone.isin([1, 3, 4])].copy()
    df["on_mat"] = (df.zone == 1).astype(int)
    df["cls"] = classify(df.canopy2024_p90_m)
    df["cls_2011"] = classify(df.canopy2011_p90_m)
    df["cls_wcs"] = classify(df.lidar_vegh_p90_m)
    z = {y: (df[f"cir_ndvi_{y}_mean"] - df[f"cir_ndvi_{y}_mean"].mean()) / df[f"cir_ndvi_{y}_mean"].std()
         for y in (2021, 2022, 2023)}
    spread = pd.concat(z, axis=1).max(1) - pd.concat(z, axis=1).min(1)
    df["cir_stable"] = spread <= 1.0
    df["cell"] = [f"r{r}c{c}" for r, c in zip(df.row, df["col"])]
    return df.set_index("cell"), (r6, c6)


# ------------------------------------------------------------------------------ pairs and stacks

def gold_p6() -> pd.DataFrame:
    cols = ["period", "track", "grid", "plot", "pair", "t1", "t2", "revisit_days", "platforms", "ref_hard_unw",
            "coh_1x1", "d_wtd_mat_mean_cm", "rain_sum_mm", "wind_ms_t1", "wind_ms_t2",
            "flag_frozen_air_t1", "flag_frozen_air_t2", "flag_snow_imgw_t1", "flag_snow_imgw_t2",
            "flag_snow_om_t1", "flag_snow_om_t2"]
    g = pd.read_csv(CUBE / "gold" / "plots_pair.csv.gz", usecols=cols, low_memory=False)
    g = g[(g["plot"] == "P6") & (g.grid == "40m") & g.period.astype(str).isin(SETS)].copy()
    g["period"] = g.period.astype(str)
    for c in ("t1", "t2"):
        g[c] = pd.to_datetime(g[c], format="ISO8601", utc=True)
    g = g[~g.platforms.astype(str).str.contains("S1D")]
    g["between"] = X126.acquisitions_between(g)
    fl = [c for c in g.columns if c.startswith("flag_")]
    g["frost_snow"] = g[fl].apply(lambda s: s.map({True: 1, False: 0, "True": 1, "False": 0}).astype(float)).max(1)
    g["wind"] = g[["wind_ms_t1", "wind_ms_t2"]].mean(1)
    g["abs_dwtd_cm"] = g.d_wtd_mat_mean_cm.abs()
    mid = g.t1 + (g.t2 - g.t1) / 2
    g["mid_month"] = mid.dt.month
    return g.reset_index(drop=True)


def stacks(period: str, track: str, cells: pd.DataFrame, pairs: pd.Series) -> dict[str, np.ndarray]:
    """(pairs × cells) coherence and unwrapped phase of the silver interferograms, in the order of ``pairs``."""
    ds = xr.open_dataset(CUBE / "silver" / f"ifg_{period}_{track}.nc")
    pos = pd.Index(ds.pair.values).get_indexer(pairs)
    if (pos < 0).any():
        raise RuntimeError(f"{period} {track}: {int((pos < 0).sum())} gold pairs missing from silver")
    ry, cx = xr.DataArray(cells.row.to_numpy(), dims="cell"), xr.DataArray(cells["col"].to_numpy(), dims="cell")
    out = {}
    for v in ("coh", "unw"):
        a = ds[v].isel(y=ry, x=cx).values.astype(float)        # (pair, cell)
        out[v] = a[pos]
    return out


# ------------------------------------------------------------------------------ per-cell outcomes

def slope_rows(x: np.ndarray, Y: np.ndarray, min_n: int = 10) -> np.ndarray:
    """Least-squares slope of every column of Y (pairs × cells) on x (pairs), over the rows where both are finite."""
    m = np.isfinite(Y) & np.isfinite(x)[:, None]
    n = m.sum(0)
    xm = np.where(m, x[:, None], 0).sum(0) / np.maximum(n, 1)
    ym = np.where(m, Y, 0).sum(0) / np.maximum(n, 1)
    dx = np.where(m, x[:, None] - xm, 0)
    dy = np.where(m, Y - ym, 0)
    with np.errstate(invalid="ignore", divide="ignore"):
        s = (dx * dy).sum(0) / (dx ** 2).sum(0)
    return np.where(n >= min_n, s, np.nan)


def nanmed(a, axis=0):
    with np.errstate(all="ignore"), __import__("warnings").catch_warnings():
        __import__("warnings").simplefilter("ignore", RuntimeWarning)
        return np.nanmedian(a, axis=axis)


def outcomes(cells: pd.DataFrame, g: pd.DataFrame, laser_pairs: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, list]:
    """Per cell × track × set: O1 (all, summer, winter), O2 (2020_2021 lags), O4, O5; per cell × track: O3.
    Also the stage-2 rows (cell × consecutive pair) and P6's own-cell values for the stop check."""
    rows, s2rows, checks = [], [], []
    p6 = cells.index[(cells.row == P6[0]) & (cells["col"] == P6[1])]
    for track in TRACKS:
        rtc = xr.open_dataset(CUBE / "silver" / f"rtc_{track}.nc")
        ry, cx = xr.DataArray(cells.row.to_numpy(), dims="cell"), xr.DataArray(cells["col"].to_numpy(), dims="cell")
        vv_all, vh_all = rtc.vv_db.isel(y=ry, x=cx), rtc.vh_db.isel(y=ry, x=cx)
        o3_inc = {}
        for period, rev in SETS.items():
            gp = g[(g.track == track) & (g.period == period)].reset_index(drop=True)
            st = stacks(period, track, cells, gp.pair)
            coh, unw = st["coh"], st["unw"]
            cons = ((gp.between == 0) & (gp.revisit_days == rev)).to_numpy()
            mon = gp.t1.dt.month.to_numpy()
            rec = {"O1_coh": nanmed(coh[cons]),
                   "O1_coh_summer": nanmed(coh[cons & np.isin(mon, SUMMER)]),
                   "O1_coh_winter": nanmed(coh[cons & np.isin(mon, WINTER)]),
                   "O5_coh_per_cm_dwtd": slope_rows(gp.abs_dwtd_cm.to_numpy()[cons], coh[cons])}
            if period == "2020_2021":
                for lag in LAGS:
                    rec[f"O2_coh_{lag}d"] = nanmed(coh[(gp.revisit_days == lag).to_numpy()])
            yrs = {"2020_2021": (2020, 2021), "2022_2024": (2022, 2024), "2025": (2025, 2025)}[period]
            tsel = (rtc.time.dt.year >= yrs[0]) & (rtc.time.dt.year <= yrs[1])
            vv, vh = vv_all.sel(time=tsel).values, vh_all.sel(time=tsel).values
            rec["O4_vv_median_db"], rec["O4_vh_median_db"] = nanmed(vv), nanmed(vh)
            with np.errstate(all="ignore"), __import__("warnings").catch_warnings():
                __import__("warnings").simplefilter("ignore", RuntimeWarning)
                rec["O4_vv_sd_db"] = np.nanstd(vv, 0)
            for k, v in rec.items():
                rows.append(pd.DataFrame({"cell": cells.index, "track": track, "set": period, "outcome": k, "value": v}))
            # P6's own coherence against the cube's gold (stop check)
            if len(p6):
                j = cells.index.get_loc(p6[0])
                d = np.abs(coh[:, j] - gp.coh_1x1.to_numpy(float))
                checks.append({"check": f"coh_1x1 {track} {period}", "n": int(np.isfinite(d).sum()),
                               "max_abs_diff": float(np.nanmax(d)), "ok": bool(np.nanmax(d) < 1e-9)})
            # O3: X-126's consecutive six-day laser pairs of 2021 and 2025, every cell on the hard-target reference
            lp = laser_pairs[(laser_pairs.track == track) & (laser_pairs.period.astype(str) == period)]
            if len(lp):
                ix = gp.reset_index().set_index("pair").loc[lp.pair, "index"].to_numpy()
                inc = (unw[ix] - gp.ref_hard_unw.to_numpy(float)[ix, None]) * PHASE_TO_MM
                o3_inc[period] = (lp, inc)
            # stage-2 rows: every consecutive pair of the set
            ci = np.flatnonzero(cons)
            z = np.arctanh(np.clip(coh[ci], 1e-6, 0.999))
            ndvi = NDVI_MONTH.loc[cells.index].to_numpy()                       # cells × 12
            for k, i in enumerate(ci):
                s2rows.append(pd.DataFrame({
                    "cell": cells.index, "track": track, "set": period, "pair": gp.pair[i], "z": z[k],
                    "abs_dwtd_cm": gp.abs_dwtd_cm[i], "rain_mm": gp.rain_sum_mm[i], "frost_snow": gp.frost_snow[i],
                    "wind_ms": gp.wind[i], "pair_days": gp.revisit_days[i],
                    "ndvi_season": ndvi[:, gp.mid_month[i] - 1]}))
        lp = pd.concat([v[0] for v in o3_inc.values()]) if o3_inc else pd.DataFrame()
        if len(lp):
            inc = np.vstack([v[1] for v in o3_inc.values()])
            a = lp.d_laser_los_mm.to_numpy(float)
            r = np.array([increment_stats(a, inc[:, j])["pearson_increments"] for j in range(inc.shape[1])])
            rows.append(pd.DataFrame({"cell": cells.index, "track": track, "set": "2021+2025", "outcome": "O3_r_laser",
                                      "value": r}))
            for yr, (lpy, incy) in {k: v for k, v in o3_inc.items()}.items():
                ay = lpy.d_laser_los_mm.to_numpy(float)
                ry_ = np.array([increment_stats(ay, incy[:, j])["pearson_increments"] for j in range(incy.shape[1])])
                rows.append(pd.DataFrame({"cell": cells.index, "track": track, "set": yr, "outcome": "O3_r_laser",
                                          "value": ry_}))
            if len(p6):
                j = cells.index.get_loc(p6[0])
                ref = pd.read_csv(X126_TABLE)
                for yr, (lpy, incy) in o3_inc.items():
                    grp = "2021" if yr == "2020_2021" else "2025"
                    mine = increment_stats(lpy.d_laser_los_mm.to_numpy(float), incy[:, j])
                    th = ref[(ref.track == track) & (ref.period.astype(str) == grp) & (ref.revisit_days == 6) &
                             (ref.pairs_kind == "consecutive") & (ref.reference == "hard_targets")].iloc[0]
                    checks.append({"check": f"X-126 r {track} {grp}", "n": mine["pairs"],
                                   "max_abs_diff": abs(mine["pearson_increments"] - th.pearson_increments),
                                   "ok": bool(mine["pairs"] == th.pairs and
                                              abs(mine["pearson_increments"] - th.pearson_increments) < 1e-9)})
    return pd.concat(rows, ignore_index=True), pd.concat(s2rows, ignore_index=True), checks


def laser_pairs_x126() -> pd.DataFrame:
    """X-126's pairs that its by-revisit rows for 2021 and 2025 (six-day, consecutive) were computed from."""
    p = X126.load_pairs()
    p = p[(p.between == 0) & (p.revisit_days == 6) & p.group.isin(["2021", "2025"])]
    return p[["track", "period", "pair", "group", "d_laser_los_mm"]].reset_index(drop=True)


def ndvi_month(cells: pd.DataFrame) -> pd.DataFrame:
    """Each cell's median Sentinel-2 NDVI per calendar month over the 2017–2025 scenes (the season of a pair)."""
    s2 = xr.open_dataset(CUBE / "silver" / "s2_40m.nc")
    s2 = s2.sel(time=s2.time.dt.year <= 2025)
    ry, cx = xr.DataArray(cells.row.to_numpy(), dims="cell"), xr.DataArray(cells["col"].to_numpy(), dims="cell")
    nd = s2.ndvi.isel(y=ry, x=cx)
    m = nd.groupby("time.month").median("time", skipna=True).transpose("cell", "month").values
    return pd.DataFrame(m, index=cells.index, columns=range(1, 13))


# ------------------------------------------------------------------------------ stage 1

def matched(cells: pd.DataFrame, cls_col: str = "cls") -> dict[str, pd.DataFrame]:
    out = {}
    for side, on_mat, cov in (("on_mat", 1, "depth_in_mat_m"), ("off_mat", 0, "dist_to_mat_m")):
        s = cells[cells.on_mat == on_mat]
        w, o = s[s[cls_col] == "wooded"], s[s[cls_col] == "open"]
        m = C.match_nearest(w, o, on=cov, caliper=CALIPER_M, tie="dist_p6_m", order="dist_p6_m")
        m["block"] = cells.loc[m.treated, "spatial_block"].to_numpy()
        m.attrs = {"wooded": len(w), "open": len(o)}
        out[side] = m
    return out


def stage1(cells: pd.DataFrame, oc: pd.DataFrame, variant: str, cls_col: str = "cls") -> pd.DataFrame:
    M = matched(cells, cls_col)
    blocks = np.concatenate([M[s].block.to_numpy() for s in M])
    ub, counts = C.block_counts(blocks, N_BOOT, SEED)
    W = {s: C.unit_weights(M[s].block.to_numpy(), ub, counts) for s in M}
    rows = []
    wide = oc.pivot_table(index="cell", columns=["track", "set", "outcome"], values="value", dropna=False)
    wide = wide.dropna(axis=1, how="all")
    for key in wide.columns:
        track, st, outcome = key
        v = wide[key]
        d = {s: (v.reindex(M[s].treated).to_numpy() - v.reindex(M[s].control).to_numpy()) for s in M}
        est = {s: C.weighted_median(d[s], np.ones(len(d[s]))) for s in M}
        draws = {s: np.array([C.weighted_median(d[s], w) for w in W[s]]) for s in M}
        draws["interaction"] = draws["on_mat"] - draws["off_mat"]
        est["interaction"] = est["on_mat"] - est["off_mat"]
        for s in ("on_mat", "off_mat", "interaction"):
            lo, hi = C.boot_interval(draws[s])
            n = {k: int(np.isfinite(d[k]).sum()) for k in M}
            rows.append({"variant": variant, "track": track, "set": st, "outcome": outcome, "contrast": s,
                         "pairs_on_mat": n["on_mat"], "pairs_off_mat": n["off_mat"],
                         "wooded_minus_open": est[s], "lo95": lo, "hi95": hi})
    t = pd.DataFrame(rows)
    t.attrs["counts"] = {s: {**M[s].attrs, "matched": len(M[s])} for s in M}
    return t


def reading(t: pd.DataFrame) -> pd.DataFrame:
    """The pre-registered reading rules, on the primary variant."""
    p = t[t.variant == "primary"]
    o1 = p[p.outcome == OUTCOMES_PRIMARY]
    sets = list(SETS)

    def count(contrast, test):
        k = 0
        for st in sets:
            g = o1[(o1.contrast == contrast) & (o1.set == st)]
            if len(g) == 2 and all(test(r) for r in g.itertuples()):
                k += 1
        return k
    below = lambda r: r.hi95 < 0  # noqa: E731
    excl = lambda r: r.lo95 > 0 or r.hi95 < 0  # noqa: E731
    o3 = p[(p.outcome == "O3_r_laser") & (p.set == "2021+2025") & (p.contrast == "on_mat")]
    rows = [
        {"rule": "trees lower coherence on the mat", "test": "(i) interval below zero on both tracks, in ≥ 2 of 3 sets",
         "sets_meeting": count("on_mat", below), "holds": count("on_mat", below) >= 2},
        {"rule": "trees act differently on the mat", "test": "(iii) interval excludes zero on both tracks, in ≥ 2 sets",
         "sets_meeting": count("interaction", excl), "holds": count("interaction", excl) >= 2},
        {"rule": "trees weaken the six-day motion read", "test": "(iv) O3 interval below zero on both tracks",
         "sets_meeting": int(len(o3) == 2 and all(o3.hi95 < 0)), "holds": bool(len(o3) == 2 and all(o3.hi95 < 0))}]
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------------ stage 2

FACTORS = ["abs_dwtd_cm", "rain_mm", "frost_snow", "wind_ms", "canopy2024_p90_m", "canopy2024_cover", "wc_tree",
           "depth_in_mat_m", "on_mat", "boardwalk_share", "ndvi_season", "pair_days", "descending"]
BINARY = {"frost_snow", "on_mat", "descending"}
INTERACTIONS = {"canopy_x_on_mat": ("canopy2024_p90_m", "on_mat"), "dwtd_x_on_mat": ("abs_dwtd_cm", "on_mat")}


def stage2(cells: pd.DataFrame, s2: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    d = s2.join(cells[["canopy2024_p90_m", "canopy2024_cover", "wc_tree", "depth_in_mat_m", "on_mat",
                       "boardwalk_share", "spatial_block"]], on="cell")
    d["descending"] = (d.track == "descending").astype(float)
    d = d.dropna(subset=FACTORS + ["z"])
    cols, X = [], [np.ones(len(d))]
    for f in FACTORS:
        v = d[f].to_numpy(float)
        X.append(v if f in BINARY else (v - v.mean()) / v.std())
        cols.append(f)
    for name, (a, b) in INTERACTIONS.items():
        X.append(X[1 + cols.index(a)] * X[1 + cols.index(b)])
        cols.append(name)
    X = np.column_stack(X)
    S = C.block_sums(X, d.z.to_numpy(float), d.spatial_block.to_numpy())
    groups = {c: [1 + i] for i, c in enumerate(cols)}
    res = C.ols_block_bootstrap(S, groups, N_BOOT, SEED)
    rows = []
    for i, c in enumerate(cols):
        lo, hi = C.boot_interval(res["beta_draws"][:, 1 + i])
        dlo, dhi = C.boot_interval(res["drop_draws"][c])
        share = res["drop"][c] / res["r2"]
        rows.append({"factor": c, "scale": "0/1" if c in BINARY else ("product" if c in INTERACTIONS else "per SD"),
                     "effect_fisher_z": res["beta"][1 + i], "lo95": lo, "hi95": hi,
                     "r2_drop": res["drop"][c], "r2_drop_lo95": dlo, "r2_drop_hi95": dhi, "share_of_r2": share,
                     "matters": bool((lo > 0 or hi < 0) and share >= 0.01)})
    t = pd.DataFrame(rows)
    meta = {"rows": len(d), "cells": d.cell.nunique(), "pairs": d.pair.nunique(), "r2": res["r2"],
            "r2_lo95": C.boot_interval(res["r2_draws"])[0], "r2_hi95": C.boot_interval(res["r2_draws"])[1],
            "blocks": len(S["blocks"]), "beta": res["beta"], "cols": cols, "d": d, "X": X}
    return t, meta


def nine_cells(meta: dict, cells: pd.DataFrame) -> pd.DataFrame:
    """X-079's nine cells (P6 and its 3 × 3): observed mean Fisher z against stage 2's prediction."""
    d, X = meta["d"], meta["X"]
    pred = X @ meta["beta"]
    d = d.assign(pred=pred)
    rows = []
    for dr in (-1, 0, 1):
        for dc in (-1, 0, 1):
            cell = f"r{P6[0] + dr}c{P6[1] + dc}"
            g = d[d.cell == cell]
            if len(g):
                rows.append({"cell": cell, "offset": f"{dr:+d},{dc:+d}", "rows": len(g),
                             "depth_in_mat_m": float(cells.loc[cell, "depth_in_mat_m"]),
                             "canopy2024_p90_m": float(cells.loc[cell, "canopy2024_p90_m"]),
                             "z_observed": float(g.z.mean()), "z_predicted": float(g.pred.mean())})
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------------ figures

def fig_classes(cells: pd.DataFrame, out: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import rasterio
    from pyproj import Transformer
    from rasterio.merge import merge
    raw = CUBE / "raw" / "gugik_dated"
    tab = pd.read_csv(raw / "sheets.csv")
    fs = [raw / f for f in tab[(tab.kind == "cir") & (tab.year == 2022)].file]
    t = Transformer.from_crs(32633, 2180, always_xy=True)
    ex, ey = t.transform(cells.x.to_numpy(), cells.y.to_numpy())
    b = (ex.min() - 60, ey.min() - 60, ex.max() + 60, ey.max() + 60)
    srcs = [rasterio.open(f) for f in fs]
    img, tr = merge(srcs, bounds=b, res=2.0)
    for s in srcs:
        s.close()
    fig, ax = plt.subplots(figsize=(9, 9))
    ax.imshow(np.moveaxis(img, 0, -1), extent=(b[0], b[2], b[1], b[3]))
    sty = {"wooded": ("#1b7837", "s"), "open": ("#f1a340", "s")}
    for k, (col, mk) in sty.items():
        m = (cells.cls == k).to_numpy()
        ax.scatter(ex[m], ey[m], s=10, facecolors="none", edgecolors=col, marker=mk, linewidths=0.8,
                   label=f"{k} (2024 scan, {int(m.sum())} cells)")
    p6 = cells[(cells.row == P6[0]) & (cells["col"] == P6[1])]
    px, py = t.transform(p6.x.to_numpy(), p6.y.to_numpy())
    ax.scatter(px, py, s=80, marker="*", color="white", edgecolors="black", label="P6", zorder=5)
    ax.set_title("X-132 cells on the 2022-07-19 colour-infrared orthophoto (GUGiK); PL-1992 m", fontsize=9)
    ax.legend(loc="lower left", fontsize=8)
    ax.set_xlabel("E (m)")
    ax.set_ylabel("N (m)")
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)


def fig_maps(cells: pd.DataFrame, oc: pd.DataFrame, out: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    panels = [("O1_coh", "2020_2021", "O1 median 6-day coherence 2020–21", "viridis"),
              ("O2_coh_24d", "2020_2021", "O2 median 24-day coherence 2020–21", "viridis"),
              ("O3_r_laser", "2021+2025", "O3 r with the P6 laser (6-day)", "RdBu_r"),
              ("O4_vv_median_db", "2022_2024", "O4 median VV (dB) 2022–24", "magma"),
              ("O5_coh_per_cm_dwtd", "2020_2021", "O5 coherence per cm |ΔWTD|", "RdBu_r")]
    fig, axs = plt.subplots(2, 5, figsize=(20, 8), sharex=True, sharey=True)
    for i, track in enumerate(TRACKS):
        for j, (o, st, title, cmap) in enumerate(panels):
            v = oc[(oc.track == track) & (oc.set == st) & (oc.outcome == o)].set_index("cell").value
            v = v.reindex(cells.index)
            ax = axs[i, j]
            kw = {"vmin": -1, "vmax": 1} if o == "O3_r_laser" else {}
            if o == "O5_coh_per_cm_dwtd":
                q = np.nanpercentile(np.abs(v), 95)
                kw = {"vmin": -q, "vmax": q}
            sc = ax.scatter(cells.x, cells.y, c=v, s=9, marker="s", cmap=cmap, **kw)
            ax.scatter(cells.x[cells.cls == "wooded"], cells.y[cells.cls == "wooded"], s=1, c="k")
            p6 = cells[(cells.row == P6[0]) & (cells["col"] == P6[1])]
            ax.scatter(p6.x, p6.y, marker="*", s=90, c="w", edgecolors="k")
            ax.set_title(f"{title} — {track[:3]}", fontsize=8)
            ax.set_aspect("equal")
            fig.colorbar(sc, ax=ax, shrink=0.7)
    fig.suptitle("X-132 outcomes per 40 m cell within 1 km of P6 (black dots: wooded in the 2024 scan; star: P6)",
                 fontsize=10)
    fig.savefig(out, dpi=110, bbox_inches="tight")
    plt.close(fig)


# ------------------------------------------------------------------------------ main

def md(df: pd.DataFrame, nd: int = 3) -> str:
    d = df.copy()
    for c in d.columns:
        if pd.api.types.is_float_dtype(d[c]):
            d[c] = d[c].map(lambda x: "" if pd.isna(x) else f"{x:.{nd}f}")
    return "\n".join(["| " + " | ".join(map(str, d.columns)) + " |", "|" + "---|" * len(d.columns)] +
                     ["| " + " | ".join(str(x) for x in r) + " |" for r in d.itertuples(index=False)])


P6 = (0, 0)
NDVI_MONTH = pd.DataFrame()


def main() -> None:
    global P6, NDVI_MONTH
    OUT.mkdir(parents=True, exist_ok=True)
    names = ["cells.csv", "outcomes.csv.gz", "stage1.csv", "reading.csv", "stage2.csv", "nine_cells.csv",
             "reproduction.csv", "README.md", "fig_classes_on_orthophoto.png", "fig_outcome_maps.png"]
    assert_preregistered(PREREG, [OUT / f for f in names])
    cells, P6 = load_cells()
    NDVI_MONTH = ndvi_month(cells)
    g = gold_p6()
    guard_lockbox(pd.concat([g.t1, g.t2]), Protocol.load(PROTOCOL))
    lp = laser_pairs_x126()
    oc, s2, checks = outcomes(cells, g, lp)
    rep = pd.DataFrame(checks)
    rep.to_csv(OUT / "reproduction.csv", index=False)
    if not len(rep) or not rep.ok.all():
        print(rep.to_string())
        sys.exit("stop (pre-registered): P6's own cell does not reproduce the cube / X-126; see reproduction.csv")
    cells.to_csv(OUT / "cells.csv")
    oc.to_csv(OUT / "outcomes.csv.gz", index=False)
    # stage 1: the primary classes and the pre-registered sensitivities
    variants = [("primary", cells, "cls")]
    variants.append(("stable since 2011", cells[cells.cls == cells.cls_2011], "cls"))
    variants.append(("CIR-stable 2021–2023", cells[cells.cir_stable], "cls"))
    variants.append(("undated WCS layer", cells, "cls_wcs"))
    for hi, lo in ((1.5, 0.3), (3.0, 1.0)):
        c2 = cells.copy()
        c2["cls_t"] = classify(c2.canopy2024_p90_m, hi, lo)
        variants.append((f"thresholds {hi} / {lo} m", c2, "cls_t"))
    variants.append(("radius 500 m", cells[cells.dist_p6_m <= 500], "cls"))
    t1, counts = [], {}
    for name, cc, col in variants:
        t = stage1(cc, oc[oc.cell.isin(cc.index)], name, col)
        counts[name] = t.attrs["counts"]
        t1.append(t)
        print(f"stage 1 {name}: {t.attrs['counts']}", flush=True)
    t1 = pd.concat(t1, ignore_index=True)
    t1.to_csv(OUT / "stage1.csv", index=False)
    rd = reading(t1)
    rd.to_csv(OUT / "reading.csv", index=False)
    t2, meta = stage2(cells, s2)
    t2.to_csv(OUT / "stage2.csv", index=False)
    nine = nine_cells(meta, cells)
    nine.to_csv(OUT / "nine_cells.csv", index=False)
    fig_classes(cells, OUT / "fig_classes_on_orthophoto.png")
    fig_maps(cells, oc, OUT / "fig_outcome_maps.png")
    prim = t1[(t1.variant == "primary") & (t1.outcome.isin([OUTCOMES_PRIMARY, "O3_r_laser"]))]
    sens = t1[(t1.outcome == OUTCOMES_PRIMARY) & (t1.contrast == "on_mat")].pivot_table(
        index=["track", "set"], columns="variant", values="wooded_minus_open").reset_index()
    cnt = pd.DataFrame([{"variant": k, **{f"{s} {q}": v[s][q] for s in v for q in ("wooded", "open", "matched")}}
                        for k, v in counts.items()])
    sec = t1[(t1.variant == "primary") & ~t1.outcome.isin([OUTCOMES_PRIMARY])]
    (OUT / "README.md").write_text(f"""# X-132 — trees against open ground around P6, on and off the mat (generated)

Generated by `experiments/landcover_coherence_x132.py` (`make run E=landcover_coherence_x132`); do not edit.
Pre-registration: `../preregistration/X-132.yaml`. Cells: 40 m, within {RADIUS_M:.0f} m of P6, boardwalk share < 0.2,
no hard targets, zones A, C, D. Wooded / open from the 2024 laser scan (canopy p90 ≥ 2 m / < 0.5 m; C-072).
Matching caliper {CALIPER_M:.0f} m; {N_BOOT} resamples of the 200 m spatial blocks, seed {SEED}. Design years only.
Contrasts are wooded − open, the median over matched pairs; the interaction is on the mat − off the mat.

## Reading rules (pre-registered, primary classes)

{md(rd)}

## Reproduction check (P6's own cell)

{md(rep, 12)}

## Cells and matched pairs per variant

{md(cnt)}

## Primary: O1 (median consecutive coherence) and O3 (r with the P6 laser)

{md(prim.drop(columns=["variant"]))}

## Sensitivities: the tree effect on the mat (O1), every variant

{md(sens)}

## Stage 2: Fisher-z coherence on every factor ({meta['rows']} cell × pair rows, {meta['cells']} cells, {meta['pairs']} pairs, {meta['blocks']} blocks)

R² {meta['r2']:.3f} [{meta['r2_lo95']:.3f}, {meta['r2_hi95']:.3f}]. Effects per SD of the factor (0/1 factors: the
difference between the two states); a factor "matters" if its interval excludes zero and it carries ≥ 1 % of R².

{md(stage2_sorted(t2))}

## X-079's nine cells within stage 2

{md(nine)}

## Secondary outcomes (primary classes)

{md(sec.drop(columns=["variant"]))}

## Figures

- `fig_classes_on_orthophoto.png` — the wooded and open cells on the 2022-07-19 colour-infrared orthophoto.
- `fig_outcome_maps.png` — O1–O5 per cell, both tracks.
""", encoding="utf-8")
    print(rd.to_string())
    print(stage2_sorted(t2)[["factor", "effect_fisher_z", "lo95", "hi95", "share_of_r2", "matters"]].to_string())


def stage2_sorted(t2: pd.DataFrame) -> pd.DataFrame:
    return t2.sort_values("share_of_r2", ascending=False).reset_index(drop=True)


if __name__ == "__main__":
    main()

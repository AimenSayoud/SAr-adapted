"""M3 / M5 — the model-ready (gold) tables, and the cube's check (C-057).

| table | one row / cell per | content |
|---|---|---|
| ``plots_acq.csv.gz`` | plot × acquisition (track, date) | backscatter windows, the plot's water-table antecedents, site weather and flags at the overpass, nearest Sentinel-2 / Landsat / UAV values with their lags, laser (P6), the plot's static layers |
| ``plots_pair.csv.gz`` | plot × pair (period, track, grid) | wrapped / unwrapped phase and coherence windows, the same per pair for the references (stable ground, grassland, mat, hard targets), changes of every driver, flags at both ends, laser change and cycle risk (P6), folds |
| ``triangles.csv.gz`` | closure triangle × unit (plots, zones) | wrapped closure, consistency, unwrapped closure, whole-cycle share |
| ``truth_laser_hourly.csv.gz``, ``truth_laser_daily.csv`` | hour / day | the P6 laser with its QC flags |
| ``mat_pairs_<period>_<track>[_20m].nc`` | pixel (mat, lake, grassland, stable, hard targets) × pair | phase, coherence, component; pixel static layers; pair drivers |
| ``mat_acq_<track>.nc`` | pixel × acquisition | VV, VH, RVI; nearest Sentinel-2 indices and Landsat temperature with lags; site drivers per acquisition |
"""
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from .. import field
from . import core, schema
from .build import F32, PERIODS, REPO, TRACKS, plot_units

V2_SPLIT = pd.Timestamp("2021-08-13", tz="UTC")
S2_TOL_DAYS, LST_TOL_DAYS, UAV_TOL_DAYS = 10, 16, 30
S2_IDX = ("ndvi", "ndre", "ndmi", "ndwi", "mndwi")
SITE_AT = ["air_c", "dpd_c", "wind_ms", "st_rh_pct", "st_vpd_kpa", "st_ppfd_global", "rain_3h_mm", "rain_24h_mm",
           "rain_72h_mm", "hours_since_rain", "om_era5_land_soil_temperature_0_to_7cm",
           "om_era5_land_soil_moisture_0_to_7cm", "om_era5_land_snow_depth", "om_era5_cloud_cover",
           "om_era5_shortwave_radiation", "om_era5_land_relative_humidity_2m", "imgw_dew_h", "imgw_snow_depth_cm"]
FLAGS = ["flag_frozen_air", "flag_frozen_soil_om", "flag_snow_om", "flag_snow_imgw", "flag_wet_rh95",
         "flag_dew_likely", "flag_wet_any"]
PIXEL_STATIC = ["dist_to_mat_m", "depth_in_mat_m", "spatial_block", "boardwalk_share", "wc_built", "wc_water",
                "wc_herb_wetland", "wc_grass", "wc_tree", "lidar_dtm_m", "lidar_roughness_m", "lidar_microrelief_sd_m",
                "lidar_vegh_mean_m", "lidar_vegh_p90_m", "incidence_asc_deg", "incidence_des_deg", "los_u_asc", "los_u_des"]


def period_of(date) -> str:
    y = pd.Timestamp(date).year
    return "2020_2021" if y <= 2021 else "2022_2024" if y <= 2024 else str(y)


# ------------------------------------------------------------------------------ loaders

def load_site(silver: Path) -> pd.DataFrame:
    d = pd.read_csv(silver / "site_hourly.csv.gz", index_col=0, parse_dates=True, low_memory=False)
    d.index = pd.DatetimeIndex(d.index, tz="UTC") if d.index.tz is None else d.index
    for c in d.columns:
        if c.startswith(("flag_", "censored_", "laser_snow", "laser_outlier", "laser_filled", "laser_ok")) or \
                c.endswith("ground_frozen"):
            d[c] = d[c].map({True: True, False: False, "True": True, "False": False})
    return d


def site_at(site: pd.DataFrame, times) -> pd.DataFrame:
    """Site drivers at overpass instants: nearest hour for meteorology and flags; interpolation for water table and laser."""
    t = pd.DatetimeIndex(pd.to_datetime(times, utc=True))
    cols = [c for c in SITE_AT + FLAGS if c in site]
    out = core.asof(t, site[cols], "31min").drop(columns="lag_h")
    out["wtd_mat_mean_cm"] = core.value_at(site["wtd_mat_mean_cm"], t)
    out["wtd_P6_cm"] = core.value_at(site["wtd_P6_cm"], t)
    if "laser_surface_cm" in site:
        ok = site["laser_surface_cm"].where(site["laser_ok"].fillna(False).astype(bool))
        out["laser_surface_cm"] = core.value_at(ok, t)
        out["laser_snow"] = core.asof(t, site[["laser_snow_72h"]], "31min")["laser_snow_72h"].to_numpy()
    return out.reset_index(drop=True)


def uav_by_plot() -> pd.DataFrame:
    u = field.load_uav_table()
    nir, red, re = u["REMX_NIR_842_mean"], u["REMX_Red_668_mean"], u["REMX_Red_Edge_717_mean"]
    v = pd.DataFrame({"date": u.Date, "plot": u.base_plot, "ndvi": (nir - red) / (nir + red),
                      "ndre": (nir - re) / (nir + re), "thermal_c": u["Altum_Thermal_11um_mean"], "lai": u["LAI"],
                      "ts_minus_air_c": u["Altum_Thermal_11um_mean"] - u["Air_2m_MEAN"]})
    g = v.groupby(["plot", "date"])
    m = g.mean().add_prefix("uav_")
    s = g[["ndvi", "ndre", "thermal_c", "lai"]].std().add_prefix("uav_sd_")
    return pd.concat([m, s], axis=1).reset_index()


def nearest_by_plot(tab: pd.DataFrame, keys: pd.DataFrame, time_col: str, tol_days: float, cols: list[str],
                    prefix: str, ok: pd.Series | None = None) -> pd.DataFrame:
    """For each (plot, time) in ``keys``, the nearest row of ``tab`` for that plot within ``tol_days`` (and ``ok``)."""
    t = tab if ok is None else tab[ok.to_numpy()]
    out = []
    for p, g in keys.groupby("plot", sort=False):
        src = t[t["plot"] == p].set_index(time_col)[cols].sort_index()
        src.index = pd.DatetimeIndex(src.index, tz="UTC") if src.index.tz is None else src.index
        a = core.asof(g.time_utc, src, pd.Timedelta(days=tol_days)) if len(src) else \
            pd.DataFrame({**{c: np.nan for c in cols}, "lag_h": np.nan}, index=range(len(g)))
        a.index = g.index
        a[f"{prefix}lag_days"] = a.pop("lag_h") / 24
        out.append(a.rename(columns={c: f"{prefix}{c}" for c in cols}))
    return pd.concat(out).loc[keys.index]


# ------------------------------------------------------------------------------ M3: plots

def ifg_files(silver: Path) -> list[tuple[str, str, str, Path]]:
    out = []
    for period in PERIODS:
        for track in TRACKS:
            for grid, suffix in (("40m", ""), ("20m", "_20m")):
                f = silver / f"ifg_{period}_{track}{suffix}.nc"
                if f.exists():
                    out.append((period, track, grid, f))
    return out


def plot_cells(px: pd.DataFrame, grid: str, x, y) -> tuple[np.ndarray, np.ndarray]:
    if grid == "40m":
        return px.row.to_numpy(), px.col.to_numpy()
    res = abs(float(x[1] - x[0]))
    return (np.floor((y[0] + res / 2 - px.N) / res).astype(int).to_numpy(),
            np.floor((px.E - (x[0] - res / 2)) / res).astype(int).to_numpy())


def ref_masks(static: xr.Dataset, grid: str) -> dict[str, np.ndarray]:
    z = static.zone.values
    m = {"mat": z == 1, "lake": z == 2, "grassland": z == 3, "stable": static.stable.values == 1}
    if "hard_target" in static:
        m["hard"] = static.hard_target.values == 1
    return m


def _win(a: np.ndarray, r: int, c: int, h: int = 1) -> np.ndarray:
    """(P, n) values of a (P, y, x) stack in the (2h+1)² window at (r, c)."""
    return a[:, max(r - h, 0):r + h + 1, max(c - h, 0):c + h + 1].reshape(a.shape[0], -1)


def _phasor(w: np.ndarray, wt: np.ndarray):
    ok = np.isfinite(w) & np.isfinite(wt)
    z = np.where(ok, wt * np.exp(1j * np.where(ok, w, 0)), 0).sum(1)
    s = np.where(ok, wt, 0).sum(1)
    with np.errstate(invalid="ignore", divide="ignore"):
        zz = z / s
    return np.angle(zz), np.abs(zz)


def _nanmed(a, axis=1):
    with np.errstate(all="ignore"), __import__("warnings").catch_warnings():
        __import__("warnings").simplefilter("ignore", RuntimeWarning)
        return np.nanmedian(a, axis=axis)


def pair_refs(W, U, C, masks) -> pd.DataFrame:
    out = {}
    for name, m in masks.items():
        if not m.any():
            continue
        u, w, c = U[:, m], W[:, m], C[:, m]
        out[f"ref_{name}_unw"] = _nanmed(u)
        out[f"ref_{name}_wrapped"], out[f"ref_{name}_consistency"] = _phasor(w, c)
        out[f"ref_{name}_coh"] = _nanmed(c)
    return pd.DataFrame(out)


def pair_drivers(site: pd.DataFrame, t1, t2) -> pd.DataFrame:
    a, b = site_at(site, t1), site_at(site, t2)
    d = pd.DataFrame(index=range(len(a)))
    for c in ("wtd_mat_mean_cm", "wtd_P6_cm", "laser_surface_cm", "air_c", "om_era5_land_soil_moisture_0_to_7cm"):
        if c in a:
            d[f"d_{c}"] = b[c] - a[c]
    for c in FLAGS + ["laser_snow"]:
        if c in a:
            d[f"{c}_t1"], d[f"{c}_t2"] = a[c].to_numpy(), b[c].to_numpy()
    for c in ("air_c", "dpd_c", "wind_ms", "rain_24h_mm", "imgw_dew_h"):
        d[f"{c}_t1"], d[f"{c}_t2"] = a[c].to_numpy(), b[c].to_numpy()
    rain = site["rain_mm"].fillna(0).cumsum()
    r1, r2 = core.value_at(rain, t1), core.value_at(rain, t2)
    d["rain_sum_mm"] = np.clip(r2 - r1, 0, None)
    et = site["om_era5_et0_fao_evapotranspiration"].fillna(0).cumsum()
    d["et0_sum_mm"] = core.value_at(et, t2) - core.value_at(et, t1)
    d["d_laser_vertical_mm"] = d.get("d_laser_surface_cm", np.nan) * 10
    return d


def gold_plots(silver: Path, gold: Path, log=print) -> None:
    site = load_site(silver)
    px = plot_units()
    static = xr.open_dataset(silver / "static_40m.nc")
    s20 = xr.open_dataset(silver / "static_20m.nc") if (silver / "static_20m.nc").exists() else None
    rows, tri_rows = [], []
    for period, track, grid, f in ifg_files(silver):
        ds = xr.open_dataset(f)
        st = static if grid == "40m" else s20
        if st is None or st.sizes["x"] != ds.sizes["x"]:
            log(f"  skip {f.name}: no static layers on its grid")
            continue
        W, U, C, K = (ds[v].values.astype(float) if v != "conncomp" else ds[v].values for v in ("wrapped", "unw", "coh", "conncomp"))
        masks = ref_masks(st, grid)
        refs = pair_refs(W, U, C, masks)
        t1 = pd.DatetimeIndex(ds.t1.values, tz="UTC")
        t2 = pd.DatetimeIndex(ds.t2.values, tz="UTC")
        drv = pair_drivers(site, t1, t2)
        rr, cc = plot_cells(px, grid, ds.x.values, ds.y.values)
        tr = track[:3]
        losu = st[f"los_u_{tr}"].values if f"los_u_{tr}" in st else None
        base = pd.DataFrame({"period": period, "track": track, "grid": grid, "pair": ds.pair.values.astype(str),
                             "t1": t1, "t2": t2, "revisit_days": ds.revisit_days.values,
                             "platforms": ds.platforms.values.astype(str),
                             "fold_year": core.pair_year_fold(t1, t2),
                             "v2_side": np.where(t2 < V2_SPLIT, "before", np.where(t1 >= V2_SPLIT, "after", "straddle"))})
        base = pd.concat([base, refs, drv], axis=1)
        for p, unit, r, c in zip(px["plot"], px.unit, rr, cc):
            b = base.copy()
            b["plot"], b["unit"], b["row"], b["col"] = p, unit, r, c
            b["wrapped_1x1"], b["unw_1x1"], b["coh_1x1"], b["conncomp_1x1"] = W[:, r, c], U[:, r, c], C[:, r, c], K[:, r, c]
            w9, u9, c9 = _win(W, r, c), _win(U, r, c), _win(C, r, c)
            b["wrapped_3x3"], b["consistency_3x3"] = _phasor(w9, c9)
            b["unw_3x3"], b["coh_3x3"] = _nanmed(u9), _nanmed(c9)
            b["los_u"] = losu[r, c] if losu is not None else np.nan
            if p != "P6":
                for col in [c_ for c_ in b if "laser" in c_]:
                    b[col] = np.nan
            else:
                b["d_laser_los_mm"] = b.d_laser_vertical_mm * b.los_u
                b["laser_cycle_risk"] = np.where(b.d_laser_los_mm.notna(), core.cycle_risk(b.d_laser_los_mm), np.nan)
            wt1 = core.value_at(site[f"wtd_{p}_cm"], t1)
            wt2 = core.value_at(site[f"wtd_{p}_cm"], t2)
            b["d_wtd_plot_cm"] = wt2 - wt1
            rows.append(b)
        tri_rows.append(triangles(period, track, grid, ds, W, U, C, masks, px, rr, cc))
        log(f"  plots × pairs: {f.name}")
    pp = pd.concat(rows, ignore_index=True)
    pp = add_rtc_changes(pp, silver)
    pp.to_csv(gold / "plots_pair.csv.gz", index=False)
    tr = pd.concat(tri_rows, ignore_index=True)
    tr.to_csv(gold / "triangles.csv.gz", index=False)
    log(f"plots_pair {len(pp)} rows, triangles {len(tr)} rows")
    acq = plots_acq(site, px, static, silver)
    acq.to_csv(gold / "plots_acq.csv.gz", index=False)
    log(f"plots_acq {len(acq)} rows × {acq.shape[1]} columns")
    truth(site, gold, log)


def add_rtc_changes(pp: pd.DataFrame, silver: Path) -> pd.DataFrame:
    """ΔVV, ΔVH (dB) over each pair at the plot's 3×3 window, where both dates have a backscatter scene."""
    px = plot_units()
    out = []
    for track, g in pp.groupby("track"):
        r = xr.open_dataset(silver / f"rtc_{track}.nc")
        dates = pd.DatetimeIndex(r.time.values).normalize()
        pos = {d: i for i, d in enumerate(dates)}
        vv, vh = r.vv_db.values, r.vh_db.values
        g = g.copy()
        for name, a in (("vv", vv), ("vh", vh)):
            vals = {}
            for p, row, col in zip(px["plot"], px.row, px.col):
                vals[p] = _nanmed(_win(a, row, col))
            i1 = g.t1.dt.tz_convert(None).dt.normalize().map(pos)
            i2 = g.t2.dt.tz_convert(None).dt.normalize().map(pos)
            v1 = [vals[p][int(i)] if pd.notna(i) else np.nan for p, i in zip(g["plot"], i1)]
            v2 = [vals[p][int(i)] if pd.notna(i) else np.nan for p, i in zip(g["plot"], i2)]
            g[f"d_{name}_db_3x3"] = np.asarray(v2) - np.asarray(v1)
        out.append(g)
    return pd.concat(out).sort_index()


def triangles(period, track, grid, ds, W, U, C, masks, px, rr, cc, max_span: int = 24) -> pd.DataFrame:
    """Closure triangles a < b < c (all three pairs present, c − a ≤ ``max_span`` days): per plot (3×3) and per zone,
    the wrapped closure (coherence-weighted phasor), its consistency, the unwrapped closure relative to stable ground,
    and the share of pixels whose unwrapped closure is off by whole cycles (X-070's measure)."""
    pairs = list(ds.pair.values.astype(str))
    idx = {p: k for k, p in enumerate(pairs)}
    dates = sorted({p[:8] for p in pairs} | {p[9:17] for p in pairs})
    rows = []
    stable = masks["stable"]
    units = {**{u: None for u in px.unit.unique()}, **{z: m for z, m in masks.items()}}
    unit_rc = {u: (r, c) for u, r, c in zip(px.unit, rr, cc)}
    for i, a in enumerate(dates):
        for j in range(i + 1, len(dates)):
            b = dates[j]
            for k in range(j + 1, len(dates)):
                c_ = dates[k]
                if (pd.Timestamp(c_) - pd.Timestamp(a)).days > max_span:
                    break
                p12, p23, p13 = f"{a}_{b}", f"{b}_{c_}", f"{a}_{c_}"
                if not (p12 in idx and p23 in idx and p13 in idx):
                    continue
                i12, i23, i13 = idx[p12], idx[p23], idx[p13]
                clo = np.angle(np.exp(1j * (W[i12] + W[i23] - W[i13])))
                wt = np.nanmean(np.stack([C[i12], C[i23], C[i13]]), 0)
                raw = U[i12] + U[i23] - U[i13]
                off = raw - clo
                off0 = np.nanmedian(off[stable & np.isfinite(off)])
                cyc = np.rint((off - off0) / (2 * np.pi))
                for u, m in units.items():
                    if m is None:
                        r, c = unit_rc[u]
                        sl = (slice(max(r - 1, 0), r + 2), slice(max(c - 1, 0), c + 2))
                        w_, cw, cy, ra = clo[sl].ravel(), wt[sl].ravel(), cyc[sl].ravel(), (raw - off0)[sl].ravel()
                    else:
                        w_, cw, cy, ra = clo[m], wt[m], cyc[m], (raw - off0)[m]
                    ang, cons = core.phasor_mean(w_, cw)
                    okc = np.isfinite(cy)
                    rows.append({"period": period, "track": track, "grid": grid, "p12": p12, "p23": p23, "p13": p13,
                                 "t1": a, "t2": b, "t3": c_, "span_days": (pd.Timestamp(c_) - pd.Timestamp(a)).days,
                                 "unit": u, "closure_wrapped": ang, "closure_consistency": cons,
                                 "closure_unw": float(np.nanmedian(ra)) if np.isfinite(ra).any() else np.nan,
                                 "cycle_share": float(np.mean(cy[okc] != 0)) if okc.any() else np.nan,
                                 "cycle_mean": float(np.mean(cy[okc])) if okc.any() else np.nan,
                                 "coh_mean": float(np.nanmean(cw)) if np.isfinite(cw).any() else np.nan})
    return pd.DataFrame(rows)


def plots_acq(site: pd.DataFrame, px: pd.DataFrame, static: xr.Dataset, silver: Path) -> pd.DataFrame:
    frames = []
    for track in TRACKS:
        r = xr.open_dataset(silver / f"rtc_{track}.nc")
        rtc_dates = pd.DatetimeIndex(r.time.values).normalize()
        ifg_dates = set()
        for period in PERIODS:
            f = silver / f"ifg_{period}_{track}.nc"
            if f.exists():
                d = xr.open_dataset(f)
                ifg_dates |= set(pd.DatetimeIndex(d.t1.values).normalize()) | set(pd.DatetimeIndex(d.t2.values).normalize())
        dates = sorted(set(rtc_dates) | ifg_dates)
        times = pd.DatetimeIndex([core.overpass_time(d, track) for d in dates])
        sa = site_at(site, times)
        pos = {d: i for i, d in enumerate(rtc_dates)}
        plat = dict(zip(rtc_dates, r.platform.values.astype(str)))
        exact = dict(zip(rtc_dates, pd.DatetimeIndex(r.acq_time_utc.values)))
        for p, unit, row, col in zip(px["plot"], px.unit, px.row, px.col):
            b = pd.DataFrame({"track": track, "date": [d.date() for d in dates], "time_utc": times,
                              "acq_time_exact_utc": [exact.get(d, pd.NaT) for d in dates],
                              "platform": [plat.get(d, "?") for d in dates], "period": [period_of(d) for d in dates],
                              "plot": p, "unit": unit, "row": row, "col": col, "fold_year": core.year_folds(times)})
            for name in ("vv_db", "vh_db", "vh_vv_db", "rvi"):
                a = r[name].values
                v1 = a[:, row, col]
                w9 = _win(a, row, col)
                ii = [pos.get(d) for d in dates]
                b[f"{name}_1x1"] = [v1[i] if i is not None else np.nan for i in ii]
                med, sd = _nanmed(w9), np.nanstd(w9, 1)
                b[f"{name}_3x3"] = [med[i] if i is not None else np.nan for i in ii]
                b[f"{name}_3x3_sd"] = [sd[i] if i is not None else np.nan for i in ii]
            wf = core.window_features(site[f"wtd_{p}_cm"], times, "wtd")
            b = pd.concat([b, wf.set_index(b.index), sa.set_index(b.index)], axis=1)
            b["censored"] = core.asof(times, site[[f"censored_{p}"]], "31min")[f"censored_{p}"].to_numpy()
            if p != "P6":
                b[["laser_surface_cm", "laser_snow"]] = np.nan
            for v in static.data_vars:
                b[f"st_{v}"] = static[v].values[row, col]
            frames.append(b)
    acq = pd.concat(frames, ignore_index=True)
    keys = acq[["plot", "time_utc"]]
    s2 = pd.read_csv(silver / "s2_plots.csv.gz", parse_dates=["time_utc"])
    acq = pd.concat([acq, nearest_by_plot(s2, keys, "time_utc", S2_TOL_DAYS,
                                          [f"{k}_3x3" for k in S2_IDX] + [f"{k}_1x1" for k in S2_IDX] + ["snow_3x3"],
                                          "s2_", ok=s2.valid_3x3 >= 0.5)], axis=1)
    if (silver / "lst_plots.csv.gz").exists():
        ls = pd.read_csv(silver / "lst_plots.csv.gz", parse_dates=["time_utc"])
        acq = pd.concat([acq, nearest_by_plot(ls, keys, "time_utc", LST_TOL_DAYS, ["lst_c_3x3", "lst_c_1x1"], "",
                                              ok=ls.lst_c_3x3.notna())], axis=1)
        acq = acq.rename(columns={"lag_days": "lst_lag_days"})
    uav = uav_by_plot()
    uav["time_utc"] = pd.to_datetime(uav.date).dt.tz_localize("UTC") + pd.Timedelta(hours=11)
    ucols = [c for c in uav.columns if c.startswith("uav_")]
    nu = nearest_by_plot(uav.rename(columns=lambda c: c.replace("uav_", "", 1) if c.startswith("uav_") else c),
                         keys, "time_utc", UAV_TOL_DAYS, [c.replace("uav_", "", 1) for c in ucols], "uav_")
    return pd.concat([acq, nu], axis=1)


def truth(site: pd.DataFrame, gold: Path, log=print) -> None:
    cols = [c for c in site.columns if c.startswith("laser_")] + ["wtd_P6_cm", "cr_raw", "p6_level_on_laser_datum_cm"]
    h = site[cols].dropna(how="all", subset=["laser_surface_cm", "laser_surface_raw_cm"])
    h.to_csv(gold / "truth_laser_hourly.csv.gz")
    ok = h.laser_surface_cm.where(h.laser_ok.fillna(False).astype(bool))
    d = pd.DataFrame({"surface_cm": ok.resample("1D").mean(), "surface_sd_cm": ok.resample("1D").std(),
                      "n_ok_hours": ok.resample("1D").count()})
    d = d[d.n_ok_hours > 0]
    d.index = d.index.date
    d.index.name = "date"
    d.reset_index().to_csv(gold / "truth_laser_daily.csv", index=False)
    log(f"truth: {len(h)} laser hours, {len(d)} days with QC-passed hours")


# ------------------------------------------------------------------------------ M5: every mat and reference pixel

def gold_mat(silver: Path, gold: Path, log=print) -> None:
    site = load_site(silver)
    static = {"40m": xr.open_dataset(silver / "static_40m.nc")}
    if (silver / "static_20m.nc").exists():
        static["20m"] = xr.open_dataset(silver / "static_20m.nc")
    for period, track, grid, f in ifg_files(silver):
        st = static.get(grid)
        ds = xr.open_dataset(f)
        if st is None or st.sizes["x"] != ds.sizes["x"]:
            continue
        sel = (np.isin(st.zone.values, [1, 2, 3]) | (st.stable.values == 1) |
               ((st.hard_target.values == 1) if "hard_target" in st else False))
        rr, cc = np.nonzero(sel)
        t1, t2 = pd.DatetimeIndex(ds.t1.values, tz="UTC"), pd.DatetimeIndex(ds.t2.values, tz="UTC")
        drv = pair_drivers(site, t1, t2)
        o = xr.Dataset({v: (("pixel", "pair"), ds[v].values[:, rr, cc].T) for v in ("wrapped", "unw", "coh", "conncomp")},
                       coords={"pixel": np.arange(len(rr)), "pair": ds.pair.values, "row": ("pixel", rr), "col": ("pixel", cc),
                               "x": ("pixel", ds.x.values[cc]), "y": ("pixel", ds.y.values[rr]),
                               "t1": ("pair", ds.t1.values), "t2": ("pair", ds.t2.values),
                               "revisit_days": ("pair", ds.revisit_days.values), "platforms": ("pair", ds.platforms.values),
                               "fold_year": ("pair", core.pair_year_fold(t1, t2))})
        for v in ["zone", "stable"] + (["hard_target"] if "hard_target" in st else []) + PIXEL_STATIC:
            if v in st:
                o[f"st_{v}"] = ("pixel", st[v].values[rr, cc])
        for c in drv.columns:
            vals = drv[c].to_numpy()
            if vals.dtype == object:
                vals = pd.to_numeric(pd.Series(vals).map({True: 1, False: 0}), errors="coerce").to_numpy()
            o[f"drv_{c}"] = ("pair", vals.astype("float32"))
        o.attrs = {"period": period, "track": track, "grid": grid, "phase_units": "rad, unreferenced",
                   "pixels": "zones mat, lake, grassland + stable ground + hard targets"}
        name = f"mat_pairs_{period}_{track}{'' if grid == '40m' else '_20m'}.nc"
        enc = {v: dict(F32) for v in ("wrapped", "unw", "coh")}
        o.to_netcdf(gold / name, encoding=enc)
        log(f"  {name}: {o.sizes['pixel']} pixels × {o.sizes['pair']} pairs")
    st = static["40m"]
    sel = (np.isin(st.zone.values, [1, 2, 3]) | (st.stable.values == 1) |
           ((st.hard_target.values == 1) if "hard_target" in st else False))
    rr, cc = np.nonzero(sel)
    s2 = xr.open_dataset(silver / "s2_40m.nc")
    lst = xr.open_dataset(silver / "lst_40m.nc") if (silver / "lst_40m.nc").exists() else None
    for track in TRACKS:
        r = xr.open_dataset(silver / f"rtc_{track}.nc")
        dates = pd.DatetimeIndex(r.time.values).normalize()
        times = pd.DatetimeIndex([core.overpass_time(d, track) for d in dates])
        o = xr.Dataset({v: (("pixel", "time"), r[v].values[:, rr, cc].T) for v in ("vv_db", "vh_db", "rvi")},
                       coords={"pixel": np.arange(len(rr)), "time": times.tz_convert(None), "row": ("pixel", rr),
                               "col": ("pixel", cc), "platform": ("time", r.platform.values)})
        for k in S2_IDX:
            v, lag = nearest_pixel(s2[k].values[:, rr, cc], s2.valid_frac.values[:, rr, cc],
                                   pd.DatetimeIndex(s2.time.values, tz="UTC"), times, S2_TOL_DAYS)
            o[f"s2_{k}"] = (("pixel", "time"), v.T.astype("float32"))
        o["s2_lag_days"] = (("pixel", "time"), lag.T.astype("float32"))
        if lst is not None:
            v, lag = nearest_pixel(lst.lst_c.values[:, rr, cc], np.isfinite(lst.lst_c.values[:, rr, cc]).astype(float),
                                   pd.DatetimeIndex(lst.time.values, tz="UTC"), times, LST_TOL_DAYS)
            o["lst_c"], o["lst_lag_days"] = (("pixel", "time"), v.T.astype("float32")), (("pixel", "time"), lag.T.astype("float32"))
        sa = site_at(site, times)
        for c in sa.columns:
            vals = sa[c].to_numpy()
            if vals.dtype == object:
                vals = pd.to_numeric(pd.Series(vals).map({True: 1, False: 0}), errors="coerce").to_numpy()
            o[f"site_{c}"] = ("time", vals.astype("float32"))
        o.to_netcdf(gold / f"mat_acq_{track}.nc", encoding={v: dict(F32) for v in o.data_vars if o[v].dtype.kind == "f"})
        log(f"  mat_acq_{track}.nc: {o.sizes['pixel']} pixels × {o.sizes['time']} acquisitions")


def nearest_pixel(vals: np.ndarray, valid: np.ndarray, src_t: pd.DatetimeIndex, times: pd.DatetimeIndex, tol_days: float,
                  min_valid: float = 0.5):
    """Per pixel and query time, the value of the nearest source scene (within ``tol_days``) where that pixel is valid."""
    T, N = len(times), vals.shape[1]
    out = np.full((T, N), np.nan)
    lag = np.full((T, N), np.nan)
    st = src_t.as_unit("ns").asi8 / 8.64e13
    for i, t in enumerate(times.as_unit("ns").asi8 / 8.64e13):
        d = st - t
        cand = np.flatnonzero(np.abs(d) <= tol_days)
        if not len(cand):
            continue
        cand = cand[np.argsort(np.abs(d[cand]))]
        todo = np.ones(N, bool)
        for k in cand:
            ok = todo & np.isfinite(vals[k]) & (valid[k] >= min_valid)
            out[i, ok] = vals[k, ok]
            lag[i, ok] = d[k]
            todo &= ~ok
            if not todo.any():
                break
    return out, lag


# ------------------------------------------------------------------------------ check (C-057)

def _sha(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _truthy(x: pd.Series) -> pd.Series:
    return x.map({True: True, False: False, "True": True, "False": False}).fillna(False).astype(bool)


def reproduction(gold: Path) -> pd.DataFrame:
    """End-to-end check of the cube against a known result (X-065): at P6, the referenced short-pair phase against the
    laser's line-of-sight change, unfrozen and snow-free pairs, by period, track and revisit. Phase → mm as
    −φ·λ/4π (uplift positive); a cube that misaligned time, sign or geometry would not reproduce it."""
    d = pd.read_csv(gold / "plots_pair.csv.gz", low_memory=False)
    p = d[(d["plot"] == "P6") & (d.grid == "40m") & d.revisit_days.isin([6, 12])]
    rows = []
    for (per, tr, rv), g in p.groupby(["period", "track", "revisit_days"]):
        ok = g.d_laser_los_mm.notna() & g.unw_1x1.notna()
        for c in ("flag_frozen_air_t1", "flag_frozen_air_t2", "laser_snow_t1", "laser_snow_t2"):
            ok &= ~_truthy(g[c])
        if ok.sum() < 8:
            continue
        for ref in ("grassland", "stable"):
            y = -core.rad_to_los_mm(g.unw_1x1 - g[f"ref_{ref}_unw"])[ok.to_numpy()]
            x = g.d_laser_los_mm[ok].to_numpy()
            rows.append({"period": per, "track": tr, "revisit_days": rv, "reference": ref, "n": int(ok.sum()),
                         "slope": float(np.polyfit(x, y, 1)[0]), "r": float(np.corrcoef(x, y)[0, 1])})
    return pd.DataFrame(rows)


def gold_nisar_plots(silver: Path, gold: Path, log=print) -> None:
    """NISAR L-band pairs at the plots (3×3 window on the template), with the same references and P6 laser change."""
    f = silver / "nisar_lband.nc"
    if not f.exists():
        return
    ds = xr.open_dataset(f)
    site = load_site(silver)
    px = plot_units()
    st = xr.open_dataset(silver / "static_40m.nc")
    masks = ref_masks(st, "40m")
    U, W, C = ds.unw.values.astype(float), ds.wrapped.values.astype(float), ds.coh.values.astype(float)
    refs = pair_refs(W, U, C, masks)
    t1, t2 = pd.DatetimeIndex(ds.t1.values, tz="UTC"), pd.DatetimeIndex(ds.t2.values, tz="UTC")
    drv = pair_drivers(site, t1, t2)
    base = pd.concat([pd.DataFrame({"granule": ds.pair.values.astype(str), "t1": t1, "t2": t2,
                                    "revisit_days": ds.revisit_days.values, "track_no": ds.track_no.values,
                                    "direction": ds.direction.values, "product": ds.product.values,
                                    "fold_year": ds.fold_year.values}), refs, drv], axis=1)
    rows = []
    for p, unit, r, c in zip(px["plot"], px.unit, px.row, px.col):
        b = base.copy()
        b["plot"], b["unit"] = p, unit
        b["unw_3x3"], b["coh_3x3"] = _nanmed(_win(U, r, c)), _nanmed(_win(C, r, c))
        b["iono_3x3"] = _nanmed(_win(ds.iono.values.astype(float), r, c))
        b["wrapped_3x3"], b["consistency_3x3"] = _phasor(_win(W, r, c), _win(ds.coh_w.values.astype(float), r, c))
        if p != "P6":
            for col in [c_ for c_ in b if "laser" in c_]:
                b[col] = np.nan
        rows.append(b)
    out = pd.concat(rows, ignore_index=True)
    out.attrs["wavelength_m"] = ds.attrs.get("wavelength_m")
    out.to_csv(gold / "plots_pair_nisar.csv.gz", index=False)
    log(f"plots_pair_nisar {len(out)} rows")


def check(silver: Path, gold: Path, raw: Path, log=print) -> None:
    problems, cover = [], []
    rep = reproduction(gold)
    rep.to_csv(gold / "check_reproduction_x065.csv", index=False)
    six = rep[(rep.revisit_days == 6) & (rep.period == "2020_2021") & (rep.reference == "grassland")]
    if len(six) < 2 or (six.r < 0.5).any():
        problems.append("P6 6-day phase no longer follows the laser (check_reproduction_x065.csv): alignment broken")
    for name, fn in (("plots_acq", "plots_acq.csv.gz"), ("plots_pair", "plots_pair.csv.gz"),
                     ("triangles", "triangles.csv.gz"), ("truth_laser_daily", "truth_laser_daily.csv")):
        p = gold / fn
        if not p.exists():
            problems.append(f"{fn} missing")
            continue
        df = pd.read_csv(p, low_memory=False)
        problems += schema.validate(df, name)
        cover.append({"product": fn, "rows": len(df), "columns": df.shape[1]})
    for p in sorted(list(silver.glob("*.nc")) + list(gold.glob("*.nc"))):
        ds = xr.open_dataset(p)
        dims = dict(ds.sizes)
        tdim = "time" if "time" in dims else "pair" if "pair" in dims else None
        span = ""
        if tdim == "time":
            tt = pd.DatetimeIndex(ds.time.values)
            span = f"{tt.min():%Y-%m-%d} → {tt.max():%Y-%m-%d}"
        elif tdim == "pair" and "t1" in ds:
            span = f"{pd.Timestamp(ds.t1.values.min()):%Y-%m-%d} → {pd.Timestamp(ds.t2.values.max()):%Y-%m-%d}"
        cover.append({"product": f"{p.parent.name}/{p.name}", "dims": json.dumps(dims), "span": span,
                      "variables": len(ds.data_vars)})
    site = pd.read_csv(silver / "site_hourly.csv.gz", index_col=0, low_memory=False)
    cover.append({"product": "silver/site_hourly.csv.gz", "rows": len(site), "columns": site.shape[1],
                  "span": f"{site.index[0][:10]} → {site.index[-1][:10]}"})
    cov = pd.DataFrame(cover)
    cov.to_csv(gold / "coverage.csv", index=False)
    files = sorted([p for p in list(silver.iterdir()) + list(gold.iterdir()) if p.is_file() and p.name != "manifest.json"])
    try:
        sha = subprocess.run(["git", "-C", str(REPO), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
        dirty = bool(subprocess.run(["git", "-C", str(REPO), "status", "--porcelain", "src", "scripts"],
                                    capture_output=True, text=True).stdout.strip())
    except OSError:
        sha, dirty = "", None
    man = {"built_utc": f"{pd.Timestamp.now(tz='UTC'):%Y-%m-%dT%H:%M:%SZ}", "code_commit": sha, "code_dirty": dirty,
           "problems": problems,
           "files": [{"path": f"{p.parent.name}/{p.name}", "bytes": p.stat().st_size, "sha256": _sha(p)} for p in files]}
    (gold.parent / "manifest.json").write_text(json.dumps(man, indent=1))
    readme(raw, cov, problems, man, gold.parent / "README.md", rep)
    log(f"check: {len(problems)} problems; coverage.csv, manifest.json ({len(files)} files), README.md")
    for pr in problems:
        log(f"  PROBLEM {pr}")
    if problems:
        raise SystemExit(1)


def raw_inventory(raw: Path) -> pd.DataFrame:
    """What each public source delivered (counted from the files, never typed)."""
    rows = []

    def add(src, what, n, span=""):
        rows.append({"source": src, "content": what, "records": n, "span": span})

    def span_of(s):
        s = pd.to_datetime(s, utc=True, format="mixed")
        return f"{s.min():%Y-%m-%d} → {s.max():%Y-%m-%d}" if len(s) else ""

    f = raw / "s2" / "scenes.csv"
    if f.exists():
        d = pd.read_csv(f)
        add("Sentinel-2 L2A (Earth Search)", "scenes read / with ≥ 2 % valid surface", f"{len(d)} / {int((d.valid_share >= 0.02).sum())}", span_of(d.time))
    f = raw / "landsat" / "scenes.csv"
    if f.exists():
        d = pd.read_csv(f)
        add("Landsat 8/9 C2 L2 (Planetary Computer)", "scenes read / with ≥ 2 % clear", f"{len(d)} / {int((d.clear_share >= 0.02).sum())}", span_of(d.time))
    f = raw / "open_meteo" / "hourly.csv.gz"
    if f.exists():
        d = pd.read_csv(f, usecols=[0])
        add("ERA5-Land + ERA5 (Open-Meteo)", "hours", len(d), span_of(d.iloc[:, 0]))
    f = raw / "imgw" / "daily.csv"
    if f.exists():
        d = pd.read_csv(f)
        add("IMGW synop daily (Piła 230, Poznań 330)", "station-days", len(d), span_of(d.date))
    f = raw / "lidar" / "tiles.csv"
    n = len(list((raw / "lidar").glob("*.tif"))) if (raw / "lidar").exists() else 0
    if n:
        add("GUGiK LiDAR WCS (EVRF2007)", "DTM 1 m + DSM 0.5 m tiles of 500 m", n)
    n = len(list((raw / "worldcover").glob("*.npz"))) if (raw / "worldcover").exists() else 0
    if n:
        add("ESA WorldCover 10 m", "maps (2020, 2021)", n)
    f = raw / "nisar" / "granules.csv"
    if f.exists():
        d = pd.read_csv(f)
        add("NISAR L2 GUNW (CMR, range requests)", "granules cropped", len(d), span_of(d.t1))
    return pd.DataFrame(rows)


def readme(raw: Path, cov: pd.DataFrame, problems: list[str], man: dict, out: Path, rep: pd.DataFrame) -> None:
    def md(df):
        cols = list(df.columns)
        lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
        for r in df.itertuples(index=False):
            lines.append("| " + " | ".join("" if (isinstance(v, float) and np.isnan(v)) else str(v) for v in r) + " |")
        return "\n".join(lines)

    inv = raw_inventory(raw)
    txt = [
        "# The data cube (chain data-cube, C-053…C-057)", "",
        "Generated by `05_code/SAr-adapted/scripts/cube_build.py check` — every count below is read from the files. "
        "Code: `insar_wetlands.cube` (public, no field values); this folder is private (field values from silver on).", "",
        "Rebuild, from `05_code/SAr-adapted`:", "",
        "```bash", "export INSAR_DRIVE_ROOT=../local/drive_pristine PYTHONPATH=src",
        "python scripts/cube_build.py ingest s2 landsat worldcover meteo imgw lidar nisar   # M1, resumable",
        "python scripts/cube_build.py harmonise          # M2: ifg rtc optical static site folds (+ nisar)",
        "python scripts/cube_build.py gold               # M3 plots + M5 mat",
        "python scripts/cube_build.py check              # contracts, coverage, manifest, this README", "```", "",
        "## Layers", "",
        "- `raw/` public downloads as received (never edited).",
        "- `silver/` one space (40 m template; 2026 also on its 20 m grid), one time base (UTC), one set of units; "
        "phase in radians, unreferenced; flags are columns.",
        "- `gold/` model-ready: per plot (`plots_acq`, `plots_pair`, `triangles`, `truth_laser_*`) and per mat / reference "
        "pixel (`mat_pairs_*`, `mat_acq_*`).", "",
        "## Public sources ingested (M1)", "", md(inv) if len(inv) else "(none)", "",
        "## Products (M2, M3, M5)", "", md(cov.fillna("")), "",
        "## Conventions", "",
        "- **Time:** UTC; an overpass is an instant (ascending ≈ 16:36, descending ≈ 05:09; the exact time is kept where "
        "the backscatter scene records it). Joins across sources are as-of joins with a tolerance and a lag column "
        f"(Sentinel-2 ± {S2_TOL_DAYS} d, Landsat ± {LST_TOL_DAYS} d, UAV ± {UAV_TOL_DAYS} d; site weather nearest hour; "
        "water table and laser interpolated).",
        "- **Phase:** radians as delivered by HyP3 / NISAR, no reference subtracted; each pair carries reference values "
        "(stable ground, grassland, mat, hard targets) so a model chooses its reference. Millimetres through "
        f"`cube.core.WAVELENGTH_MM` = {core.WAVELENGTH_MM} (C-052 open).",
        "- **Flags** are columns with their source (`flag_frozen_air`, `flag_frozen_soil_om`, `flag_snow_om`, "
        "`flag_snow_imgw`, `flag_wet_rh95`, `flag_dew_likely`, `laser_*`, `censored_P*`); nothing is filtered out.",
        "- **Folds** (`silver/folds.json`): calendar year (a pair straddling two years belongs to none), the X-066 split, "
        "radar unit (P8 + P9 share a cell), 200 m spatial blocks.",
        "- `p6_level_on_laser_datum_cm` = laser surface + P6 water-table depth: the P6 water level on the laser's datum, "
        "valid **only if** the logger's depth is measured from the moving surface — an open question to the field team.", "",
        "## Open questions to the field team (they change how columns are read)", "",
        "1. Which hourly water-table values are gap-filled (no column marks them yet).",
        "2. How the P5/P6 piezometers are anchored: is `WTD` depth below the moving surface or a level on a fixed datum?",
        "3. The UAV DSM/DTM and original TIFFs (not in the cube yet).",
        "4. Water table and laser for 2025–2026 (the cube's 2025–2026 radar rows have no field partner yet).", "",
        "## End-to-end check against a known result (X-065)", "",
        "At P6, the referenced short-pair phase (−φ·λ/4π, uplift positive) against the laser's line-of-sight change, "
        "unfrozen and snow-free pairs (`gold/check_reproduction_x065.csv`). The build fails if the 2020–2021 six-day "
        "correlation drops below 0.5 on either track — the signature of a time, sign or geometry misalignment.", "",
        md(rep.round(2)) if len(rep) else "(no laser pairs)", "",
        f"## Check — {len(problems)} problem(s)", "",
        *(f"- {p}" for p in problems), "",
        f"Built {man['built_utc']} from SAr-adapted commit `{man['code_commit'][:9]}`"
        f"{' (uncommitted changes)' if man['code_dirty'] else ''}; {len(man['files'])} files hashed in `manifest.json`.", ""]
    out.write_text("\n".join(txt))


def run(what: list[str], ctx_fn, silver: Path, gold: Path, log=print) -> None:
    todo = what or ["plots", "mat", "nisar"]
    if "plots" in todo:
        gold_plots(silver, gold, log)
    if "nisar" in todo:
        gold_nisar_plots(silver, gold, log)
    if "mat" in todo:
        gold_mat(silver, gold, log)

"""C-065: intake check of the UAV orthomosaics (Altum, RedEdge-MX Dual, 2022–2024) — what was delivered, whether
each flight's values can be trusted, whether the flights line up, and which 40 m cells they cover.

No radar value is read: this is the data's own quality, before any question is put to it. Seven parts:

1. Inventory: per file, sensor, date, bands, size, pixel size in metres, CRS written or not, valid share.
2. Radiometry: per file and band, the median, 1st and 99th percentiles, the share of pixels at the clip value
   (reflectance 2.0, the export's ceiling) and above 1 (a reflectance above 1 is not physical for vegetation).
   RedEdge-MX neighbouring bands (444/475, 531/560, 650/668) see almost the same light: a ratio of their medians far
   from 1 means one band is mis-scaled.
3. The two cameras on the same day (14 days): the ratio of their medians on the five shared bands, and the per-pixel
   agreement of NDVI, on a 0.5 m grid.
4. Alignment between flights: the shift of each flight relative to a reference flight of the same camera, and of
   the two cameras on the same day, in three windows along the boardwalk (south, middle, north), at 5 cm, by phase
   correlation of the image edges (`uav.phase_shift`, tested on known shifts). Three windows show whether a mosaic
   is shifted or warped.
5. Effective resolution: the lag at which each band's autocorrelation falls below 0.5 in a window of mat.
6. The field team's per-plot tables: are the values in the first delivery's merged table (used by X-049) the same
   as in these per-plot CSVs; and, where we hold the sub-plot circles (P5, P6), do the rasters give the same means.
7. Coverage of the cube's 40 m cells (zone, depth into the mat, plot cell, boardwalk share), and per cell and flight
   the median NDVI, NDRE (RedEdge-MX) and surface temperature (Altum) — a processed table for later, pre-registered
   questions; nothing is compared with radar here.

Outputs: the hub's ``08_deliverables/field_uav_c065/`` (tables, figures, generated README) and the 0.5 m stack
``06_data/field/processed/uav_0p5m.nc``.

    PYTHONPATH=src python experiments/uav_rasters_c065.py
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import rasterio  # noqa: E402
import xarray as xr  # noqa: E402
from rasterio.features import geometry_mask  # noqa: E402
from scipy.ndimage import binary_dilation  # noqa: E402

from insar_wetlands import field, uav  # noqa: E402
from insar_wetlands.paths import cube_root, hub_root  # noqa: E402

HUB, CUBE = hub_root(), cube_root()
OUT = HUB / "08_deliverables" / "field_uav_c065"
STACK = field.field_root() / "processed" / "uav_0p5m.nc"
RES, RES_FINE = 0.5, 0.05
CLIP = 1.999           # the export's ceiling is 2.0 (1.999969 in float32)
SHARED = ["Blue_475", "Green_560", "Red_668", "Red_Edge_717", "NIR_842"]
NEIGHBOURS = [("Blue_444", "Blue_475"), ("Green_531", "Green_560"), ("Red_650", "Red_668")]
# Windows for alignment (60 m) and sharpness (20 m), placed from the private plot shapefile at run time (no site
# coordinate in this public file): north on P8–P9's platforms, middle on P5–P6's, south on the station at the strip's
# southern end; the sharpness window on open mat east of the boardwalk, between P4 and P5.
SIDE, MAT_SIDE = 60, 20
REFERENCE = "2023-06-29"                # both cameras flew; mid-season
# Quality rules. Set at intake after a first look at the band medians (2026-10-26); the failures they catch are off
# by factors of two or more, so the thresholds decide no flag except the borderline cases the README names.
RULES = {
    "fill": "pixels where every reflectance band sits at the clip value (2.0) are a fill, not a measurement: they "
            "are set to no data (and the thermal band with them) and their share is reported",
    "saturated": "outside the fill, more than 0.1 % of valid pixels at the clip value in a band",
    "above_one": "median reflectance above 1 in a band (not physical for vegetation)",
    "band_jump": "RedEdge-MX: a neighbouring-band median ratio (444/475, 531/560, 650/668) outside 0.67–1.5",
    "out_of_series": "a band's median outside 0.4–2.5 × the same camera's median for that band over all its flights",
    "camera_mismatch": "same day, the other camera: a shared band's median ratio outside 0.67–1.5; put on the camera "
                       "further from its own series",
    "ndvi_mismatch": "same day, the other camera: NDVI median absolute difference above 0.05",
    "thermal_hot": "Altum, outside the fill: more than 0.1 % of valid pixels above 70 °C",
    "shifted": "a clear match (peak ratio ≥ 2) that puts the flight more than 0.5 m from the reference flight",
}
MATCH = 2.0


def box(x, y, side):
    return (x - side / 2, y - side / 2, x + side / 2, y + side / 2)


def windows(files):
    """Window boxes in UTM: alignment windows and the open-mat window, from the plots and the strip's extent."""
    pl = field.load_plots(crs=uav.GRID_CRS).groupby("base_plot")[["x", "y"]].mean()
    south_edge = min(uav.utm_bounds(p)[1] for p in files.path)
    x_walk = pl.loc[["P1", "P2", "P3", "P4"], "x"].mean()
    wins = {"south": box(x_walk, south_edge + SIDE / 2 + 5, SIDE),
            "middle": box(*pl.loc[["P5", "P6"]].mean(), SIDE),
            "north": box(*pl.loc[["P8", "P9"]].mean(), SIDE)}
    mid = pl.loc[["P4", "P5"]].mean()
    return wins, box(mid.x + 30, mid.y, MAT_SIDE)


# ------------------------------------------------------------------------------------------- 1, 2: per file

def inventory_and_radiometry(files):
    inv, rad = [], []
    for f in files.itertuples():
        nod = field.UAV_NODATA[f.sensor]
        with rasterio.open(f.path) as r:
            lat0 = (r.bounds.bottom + r.bounds.top) / 2
            m_lon = 111320 * np.cos(np.radians(lat0))
            a = r.read(out_shape=(r.count, r.height // 4, r.width // 4)).astype("float64")
            inv.append({"sensor": f.sensor, "date": f.date.date(), "file": f.path.name, "bands": r.count,
                        "width": r.width, "height": r.height, "pixel_e_m": abs(r.transform.a) * m_lon,
                        "pixel_n_m": abs(r.transform.e) * 110574, "crs_in_file": r.crs is not None,
                        "nodata": nod, "size_gb": f.path.stat().st_size / 1e9})
        valid = np.isfinite(a[0]) & (a[0] != nod)
        refl = a[[i for i, b in enumerate(f.bands) if b != "Thermal_11um"]]
        fill = valid & (refl >= CLIP).all(axis=0)
        inv[-1]["fill_share_of_valid"] = fill.sum() / valid.sum()
        valid &= ~fill
        inv[-1]["valid_share_of_bbox"] = valid.mean()
        for name, b in zip(f.bands, a):
            v = b[valid & np.isfinite(b) & (b != nod)]
            row = {"sensor": f.sensor, "date": f.date.date(), "band": name, "median": np.median(v),
                   "p01": np.percentile(v, 1), "p99": np.percentile(v, 99)}
            if name == "Thermal_11um":
                row["share_above_70c"] = (v > 70).mean()
            else:
                row["share_clipped"] = (v >= CLIP).mean()
                row["share_above_1"] = (v > 1).mean()
            rad.append(row)
    return pd.DataFrame(inv), pd.DataFrame(rad)


# ------------------------------------------------------------------------------------------- 0.5 m stack

def build_stack(files):
    bounds = np.array([uav.utm_bounds(p) for p in files.path])
    tr, shape = uav.grid((bounds[:, 0].min(), bounds[:, 1].min(), bounds[:, 2].max(), bounds[:, 3].max()), RES)
    xs = tr.c + RES * (np.arange(shape[1]) + 0.5)
    ys = tr.f - RES * (np.arange(shape[0]) + 0.5)
    data = {}
    for f in files.itertuples():
        arr = uav.read_on_grid(f.path, range(1, len(f.bands) + 1), field.UAV_NODATA[f.sensor], tr, shape, factor=4)
        refl = [i for i, b in enumerate(f.bands) if b != "Thermal_11um"]
        with np.errstate(invalid="ignore"):
            fill = binary_dilation((arr[refl] >= CLIP).all(axis=0), iterations=2)
        arr[:, fill] = np.nan
        tag = "altum" if f.sensor == "Altum" else "remx"
        for name, a in zip(f.bands, arr):
            data.setdefault(f"{tag}_{name}", {})[f.date] = a
    ds = xr.Dataset(coords={"y": ys, "x": xs})
    for var, by_date in data.items():
        dates = sorted(by_date)
        ds[var] = xr.DataArray(np.stack([by_date[d] for d in dates]), dims=("date_" + var.split("_")[0], "y", "x"),
                               coords={"date_" + var.split("_")[0]: pd.DatetimeIndex(dates)})
    ds.attrs.update(crs=f"EPSG:{uav.GRID_CRS}", resolution_m=RES,
                    source="06_data/field/raw/2026-10_uav_rasters (C-065), averaged from the 1/4 overviews",
                    note="reflectance as delivered (no correction); Altum thermal in °C; see field_uav_c065/qc_flags.csv")
    return ds


# ------------------------------------------------------------------------------------------- 3: the two cameras

def cross_camera(ds):
    rows, ndvi_rows = [], []
    common = sorted(set(ds.date_altum.values) & set(ds.date_remx.values))
    for d in common:
        for b in SHARED:
            a, m = ds[f"altum_{b}"].sel(date_altum=d).values, ds[f"remx_{b}"].sel(date_remx=d).values
            ok = np.isfinite(a) & np.isfinite(m)
            rows.append({"date": pd.Timestamp(d).date(), "band": b, "n_cells_0p5m": ok.sum(),
                         "altum_median": np.median(a[ok]), "remx_median": np.median(m[ok]),
                         "ratio_altum_over_remx": np.median(a[ok]) / np.median(m[ok]),
                         "pixel_r": np.corrcoef(a[ok], m[ok])[0, 1]})
        na = uav.ndvi(ds.altum_NIR_842.sel(date_altum=d).values, ds.altum_Red_668.sel(date_altum=d).values)
        nm = uav.ndvi(ds.remx_NIR_842.sel(date_remx=d).values, ds.remx_Red_668.sel(date_remx=d).values)
        ok = np.isfinite(na) & np.isfinite(nm)
        ndvi_rows.append({"date": pd.Timestamp(d).date(), "ndvi_altum_median": np.median(na[ok]),
                          "ndvi_remx_median": np.median(nm[ok]), "ndvi_pixel_r": np.corrcoef(na[ok], nm[ok])[0, 1],
                          "ndvi_median_abs_diff": np.median(np.abs(na[ok] - nm[ok]))})
    return pd.DataFrame(rows), pd.DataFrame(ndvi_rows)


# ------------------------------------------------------------------------------------------- 4, 5: geometry

def fine_band(f, window_bounds):
    tr, shape = uav.grid(window_bounds, RES_FINE)
    b = f.bands.index("NIR_842") + 1
    return uav.read_on_grid(f.path, [b], field.UAV_NODATA[f.sensor], tr, shape, window_bounds_utm=window_bounds)[0]


def alignment(files):
    rows = []
    wins, _ = windows(files)
    ref = {s: files[(files.sensor == s) & (files.date == REFERENCE)].iloc[0] for s in files.sensor.unique()}
    for wname, wb in wins.items():
        refimg = {s: uav.edges(fine_band(r, wb)) for s, r in ref.items()}
        for f in files.itertuples():
            img = uav.edges(fine_band(f, wb))
            for against, rimg in (("same camera, " + REFERENCE, refimg[f.sensor]),):
                dy, dx, q = uav.phase_shift(rimg, img)
                rows.append({"sensor": f.sensor, "date": f.date.date(), "window": wname, "against": against,
                             "shift_e_m": dx * RES_FINE, "shift_n_m": -dy * RES_FINE, "peak_ratio": q})
            if f.sensor == "Altum":
                other = files[(files.sensor != "Altum") & (files.date == f.date)]
                if len(other):
                    dy, dx, q = uav.phase_shift(uav.edges(fine_band(other.iloc[0], wb)), img)
                    rows.append({"sensor": f.sensor, "date": f.date.date(), "window": wname,
                                 "against": "RedEdge-MX, same day", "shift_e_m": dx * RES_FINE,
                                 "shift_n_m": -dy * RES_FINE, "peak_ratio": q})
    out = pd.DataFrame(rows)
    out["shift_m"] = np.hypot(out.shift_e_m, out.shift_n_m)
    return out


def effective_resolution(files):
    """Lag (m) where the autocorrelation of each band, along rows and columns, first drops below 0.5."""
    _, wb = windows(files)
    tr, shape = uav.grid(wb, 0.02)
    rows = []
    for f in files[files.date == REFERENCE].itertuples():
        arr = uav.read_on_grid(f.path, range(1, len(f.bands) + 1), field.UAV_NODATA[f.sensor], tr, shape,
                               window_bounds_utm=wb)
        for name, a in zip(f.bands, arr):
            a = a - np.nanmean(a)
            a = np.where(np.isfinite(a), a, 0)
            lags = {}
            for axis in (0, 1):
                ac = [1.0]
                for k in range(1, 150):
                    x, y = (a[k:], a[:-k]) if axis == 0 else (a[:, k:], a[:, :-k])
                    ac.append(np.corrcoef(x.ravel(), y.ravel())[0, 1])
                    if ac[-1] < 0.5:
                        break
                lags[axis] = len(ac) - 1
            rows.append({"sensor": f.sensor, "date": f.date.date(), "band": name,
                         "autocorr_half_lag_m": 0.02 * max(lags.values())})
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------------------------- 6: the per-plot tables

def plot_tables(files):
    base = field.field_root() / field.UAV_DELIVERY
    csv = pd.concat([pd.read_csv(p) for p in sorted(base.glob("*/output_csv/*_clipped*.csv"))], ignore_index=True)
    csv["Date"] = pd.to_datetime(csv["Date"])
    csv["sensor"] = csv.pop("Sensor")
    # P8 was moved in 2023: the merged table holds, for the 2022 flights, the values of the "P8_k_2022" positions
    y22 = csv.Date.dt.year == 2022
    csv = csv[~(y22 & csv.Plot.str.match(r"P8_\d$"))].copy()
    csv.loc[y22, "Plot"] = csv.loc[y22, "Plot"].str.replace(r"^(P8_\d)_2022$", r"\1", regex=True)
    csv = csv[~csv.Plot.str.endswith("_2022")]
    merged = field.load_uav_table(mask_flagged=False)
    flagged = field.load_uav_table(mask_flagged=True)
    rows = []
    for sensor, prefix in (("RedEdge_MX_Dual", "REMX_"), ("Altum", "Altum_")):
        sub = csv[csv.sensor == sensor]
        for b in field.UAV_BANDS[sensor]:
            col = f"{prefix}{b}_mean"
            if col not in merged.columns:
                continue
            j = sub[["Date", "Plot", f"{b}_mean"]].merge(merged[["Date", "Plot", col]], on=["Date", "Plot"])
            j = j.dropna()
            rows.append({"sensor": sensor, "band": b, "n_joined": len(j),
                         "max_abs_diff": (j[f"{b}_mean"] - j[col]).abs().max(),
                         "n_flagged_yellow": int(flagged[col].isna().sum() - merged[col].isna().sum())})
    same = pd.DataFrame(rows)
    # flagged (yellow) cells by date: which flights the field team had already marked
    yellow = []
    for col in [c for c in merged.columns if c.endswith("_mean")]:
        lost = merged[col].notna() & flagged[col].isna()
        for d, n in merged.loc[lost, "Date"].value_counts().items():
            yellow.append({"column": col, "date": d.date(), "n_subplots_flagged": n})
    # recompute P5/P6 sub-plot means from the rasters, where we hold the circles
    plots = field.load_plots(crs=4326)
    circles = plots[plots["plot"].str.match(r"P[56]_\d")]
    recomputed = []
    for f in files.itertuples():
        with rasterio.open(f.path) as r:
            for c in circles.itertuples():
                w = rasterio.windows.from_bounds(*c.geometry.bounds, transform=r.transform).round_offsets()
                w = rasterio.windows.Window(w.col_off - 2, w.row_off - 2, w.width + 4, w.height + 4).round_lengths()
                a = r.read(window=w).astype("float64")
                inside = ~geometry_mask([c.geometry], out_shape=a.shape[1:], transform=r.window_transform(w))
                for name, b in zip(f.bands, a):
                    v = b[inside & (b != field.UAV_NODATA[f.sensor]) & np.isfinite(b)]
                    recomputed.append({"sensor": f.sensor, "Date": f.date, "Plot": c.plot, "band": name,
                                       "raster_mean": v.mean() if v.size else np.nan, "raster_n": v.size})
    rc = pd.DataFrame(recomputed)
    long = csv.melt(id_vars=["Date", "Plot", "sensor"], value_name="csv_mean", var_name="col")
    long = long[long.col.str.endswith("_mean")].assign(band=lambda d: d.col.str[:-5]).drop(columns="col")
    rc = rc.merge(long, on=["sensor", "Date", "Plot", "band"], how="left")
    rc["rel_diff"] = (rc.raster_mean - rc.csv_mean) / rc.csv_mean.abs()
    return same, pd.DataFrame(yellow), rc


# ------------------------------------------------------------------------------------------- 7: the 40 m cells

def cells(ds, qc):
    st = xr.open_dataset(CUBE / "silver" / "static_40m.nc")
    xs, ys = st.x.values, st.y.values
    x0, x1, y0, y1 = ds.x.min().item(), ds.x.max().item(), ds.y.min().item(), ds.y.max().item()
    keep_x = xs[(xs + 20 > x0) & (xs - 20 < x1)]
    keep_y = ys[(ys + 20 > y0) & (ys - 20 < y1)]
    ok_ndvi = set(qc.loc[(qc.sensor == "RedEdge_MX_Dual") & qc.use_ndvi, "date"])
    ok_t = set(qc.loc[qc.use_thermal, "date"])
    nd_r = uav.ndvi(ds.remx_NIR_842, ds.remx_Red_668)
    ndre = uav.ndvi(ds.remx_NIR_842, ds.remx_Red_Edge_717)
    therm = ds.altum_Thermal_11um
    meta, per = [], []
    for cx in keep_x:
        for cy in keep_y:
            sel = dict(x=slice(cx - 20, cx + 20), y=slice(cy + 20, cy - 20))
            s = st.sel(x=cx, y=cy)
            cov = np.isfinite(ds.remx_NIR_842.sel(**sel)).sum(("y", "x")) / (40 / RES) ** 2
            cov_a = np.isfinite(therm.sel(**sel)).sum(("y", "x")) / (40 / RES) ** 2
            meta.append({"x": cx, "y": cy, "zone": int(s.zone), "depth_in_mat_m": float(s.depth_in_mat_m),
                         "plot_cell": int(s.plot_cell), "boardwalk_share": float(s.boardwalk_share),
                         "coverage_remx_median": float(cov.median()), "coverage_altum_median": float(cov_a.median())})
            for d in ds.date_remx.values:
                c = float(cov.sel(date_remx=d))
                if c < 0.9:
                    continue
                v, r2 = nd_r.sel(date_remx=d, **sel).values, ndre.sel(date_remx=d, **sel).values
                per.append({"x": cx, "y": cy, "sensor": "RedEdge_MX_Dual", "date": pd.Timestamp(d).date(),
                            "coverage": c, "ndvi_median": np.nanmedian(v), "ndvi_iqr": np.subtract(*np.nanpercentile(v, [75, 25])),
                            "ndre_median": np.nanmedian(r2),
                            "qc_ok": pd.Timestamp(d).date() in ok_ndvi})
            for d in ds.date_altum.values:
                c = float(cov_a.sel(date_altum=d))
                if c < 0.9:
                    continue
                v = therm.sel(date_altum=d, **sel).values
                per.append({"x": cx, "y": cy, "sensor": "Altum", "date": pd.Timestamp(d).date(), "coverage": c,
                            "thermal_median_c": np.nanmedian(v), "thermal_iqr_c": np.subtract(*np.nanpercentile(v, [75, 25])),
                            "qc_ok": pd.Timestamp(d).date() in ok_t})
    return pd.DataFrame(meta), pd.DataFrame(per)


# ------------------------------------------------------------------------------------------- flags

def flags(inv, rad, cam, ndvi_cam, align):
    typical = rad.groupby(["sensor", "band"])["median"].median()
    out = []
    for f in inv.itertuples():
        r = rad[(rad.sensor == f.sensor) & (rad.date == f.date)].set_index("band")
        refl = r.drop(index="Thermal_11um", errors="ignore")
        why, ndvi_why = [], []
        if f.fill_share_of_valid > 0.001:
            why.append(f"fill {f.fill_share_of_valid:.0%} of the strip (set to no data)")
        sat = refl.index[refl.share_clipped > 0.001]
        if len(sat):
            why.append("saturated: " + ", ".join(sat))
            ndvi_why += [b for b in sat if b in ("Red_668", "NIR_842")]
        if (refl["median"] > 1).any():
            why.append("above_one: " + ", ".join(refl.index[refl["median"] > 1]))
        if f.sensor == "RedEdge_MX_Dual":
            jumps = [f"{a}/{b}" for a, b in NEIGHBOURS if not 0.67 <= r.loc[a, "median"] / r.loc[b, "median"] <= 1.5]
            if jumps:
                why.append("band_jump: " + ", ".join(jumps))
                ndvi_why += [j for j in jumps if "668" in j]
        ratio = refl["median"] / typical.loc[f.sensor].reindex(refl.index)
        oos = ratio.index[~ratio.between(0.4, 2.5)]
        if len(oos):
            why.append("out_of_series: " + ", ".join(f"{b} ×{ratio[b]:.1f}" for b in oos))
        c = cam[cam.date == f.date]
        off = c[~c.ratio_altum_over_remx.between(0.67, 1.5)]
        if len(off):
            why.append("camera_mismatch: " + ", ".join(f"{b} ×{v:.2f}" for b, v in zip(off.band, off.ratio_altum_over_remx)))
        n = ndvi_cam[ndvi_cam.date == f.date]
        if len(n) and n.ndvi_median_abs_diff.iloc[0] > 0.05:
            why.append(f"ndvi_mismatch: {n.ndvi_median_abs_diff.iloc[0]:.2f}")
        a = align[(align.sensor == f.sensor) & (align.date == f.date) & align.against.str.startswith("same")]
        moved = a[(a.peak_ratio >= MATCH) & (a.shift_m > 0.5)]
        if len(moved):
            why.append("shifted: " + ", ".join(moved.window))
        tw = []
        if "Thermal_11um" in r.index and r.loc["Thermal_11um", "share_above_70c"] > 0.001:
            tw.append("thermal_hot")
        out.append({"sensor": f.sensor, "date": f.date, "reflectance_flags": "; ".join(why),
                    "ndvi_red_nir_flags": ", ".join(ndvi_why), "thermal_flags": "; ".join(tw),
                    "oos_red_nir": bool({"Red_668", "NIR_842"} & set(oos)),
                    "ndvi_checked_by_other_camera": len(n) > 0,
                    "ndvi_agrees": bool(len(n) and n.ndvi_median_abs_diff.iloc[0] <= 0.05)})
    q = pd.DataFrame(out)
    # A camera mismatch says the two cameras disagree, not which one is wrong: it stays on the camera whose shared
    # bands are furthest (mean |log ratio|) from that camera's own median over all its flights.
    for d in q.date[q.reflectance_flags.str.contains("camera_mismatch|ndvi_mismatch")].unique():
        dev = {}
        for s in ("Altum", "RedEdge_MX_Dual"):
            mine = rad[(rad.sensor == s) & (rad.date == d) & rad.band.isin(SHARED)].set_index("band")["median"]
            dev[s] = float(np.abs(np.log(mine / typical.loc[s].reindex(mine.index))).mean())
        blamed = max(dev, key=dev.get)
        for s in dev:
            i = q.index[(q.sensor == s) & (q.date == d)][0]
            if s != blamed:
                parts = q.loc[i, "reflectance_flags"].split("; ")
                q.loc[i, "reflectance_flags"] = "; ".join(p for p in parts if not p.startswith(("camera_mismatch", "ndvi_mismatch")))
                q.loc[i, "ndvi_agrees"] = True if q.loc[i, "ndvi_checked_by_other_camera"] else q.loc[i, "ndvi_agrees"]
    # absolute reflectance: every rule; the fill alone does not disqualify (it is masked)
    q["use_reflectance"] = q.reflectance_flags.str.split("; ").apply(lambda ps: all(p == "" or p.startswith("fill") for p in ps))
    # NDVI is a ratio: a common gain on red and NIR cancels. Usable when red and NIR are not clipped or mis-scaled
    # against their neighbours, and either the other camera's NDVI agrees that day or (no other camera) neither red
    # nor NIR is out of series.
    q["use_ndvi"] = (q.ndvi_red_nir_flags == "") & np.where(q.ndvi_checked_by_other_camera, q.ndvi_agrees, ~q.oos_red_nir)
    q["use_thermal"] = (q.sensor == "Altum") & (q.thermal_flags == "") & ~q.reflectance_flags.str.contains("shifted")
    return q.drop(columns=["oos_red_nir"])


def thermal_vs_air(rad):
    """Strip-median surface temperature per Altum flight beside that day's air temperature at the station
    (min / mean / max of the 2 m sensor, from the first delivery's merged table)."""
    m = field.load_uav_table(mask_flagged=False).groupby("Date")[["Air_2m_MIN", "Air_2m_MEAN", "Air_2m_MAX"]].first()
    t = rad[rad.band == "Thermal_11um"][["date", "median", "p01", "p99"]].rename(
        columns={"median": "surface_median_c", "p01": "surface_p01_c", "p99": "surface_p99_c"})
    t = t.merge(m.rename_axis("date").reset_index().assign(date=lambda d: d.date.dt.date), on="date", how="left")
    t["surface_minus_air_max_c"] = t.surface_median_c - t.Air_2m_MAX
    return t


# ------------------------------------------------------------------------------------------- figures

def fig_quicklook(ds, qc, path):
    dr, da = ds.date_remx.values, ds.date_altum.values
    n = max(len(dr), len(da))
    fig, axes = plt.subplots(2, n, figsize=(n * 1.05, 9.5), squeeze=False)
    nd = uav.ndvi(ds.remx_NIR_842, ds.remx_Red_668)
    for i in range(n):
        for row, (dates, img, cmap, lim, sensor) in enumerate((
                (dr, nd, "RdYlGn", (0, 0.9), "RedEdge_MX_Dual"),
                (da, ds.altum_Thermal_11um, "inferno", None, "Altum"))):
            ax = axes[row, i]
            ax.set_axis_off()
            if i >= len(dates):
                continue
            d = dates[i]
            a = img.sel({img.dims[0]: d}).values
            vmin, vmax = lim if lim else np.nanpercentile(a, [2, 98])
            ax.imshow(a, cmap=cmap, vmin=vmin, vmax=vmax, interpolation="nearest")
            ok = qc[(qc.sensor == sensor) & (qc.date == pd.Timestamp(d).date())]
            good = bool(ok.use_ndvi.iloc[0]) if row == 0 else bool(ok.use_thermal.iloc[0])
            ax.set_title(pd.Timestamp(d).strftime("%Y-%m-%d") + ("" if good else "\nflagged"), fontsize=6,
                         color="black" if good else "red")
    axes[0, 0].text(-0.3, 0.5, "NDVI (RedEdge-MX)\n0 to 0.9", transform=axes[0, 0].transAxes, rotation=90,
                    va="center", ha="right", fontsize=7)
    axes[1, 0].text(-0.3, 0.5, "surface temperature\n(Altum; 2–98 % stretch\nper flight)",
                    transform=axes[1, 0].transAxes, rotation=90, va="center", ha="right", fontsize=7)
    fig.suptitle("UAV flights over the boardwalk strip, 0.5 m (red title: NDVI or temperature not usable that day)", fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def fig_coverage(ds, meta, path):
    st = xr.open_dataset(CUBE / "silver" / "static_40m.nc")
    fig, ax = plt.subplots(figsize=(4.5, 9))
    nd = uav.ndvi(ds.remx_NIR_842, ds.remx_Red_668).sel(date_remx=REFERENCE)
    ax.imshow(nd.values, extent=(ds.x.min() - RES / 2, ds.x.max() + RES / 2, ds.y.min() - RES / 2, ds.y.max() + RES / 2),
              cmap="Greens", vmin=0, vmax=0.9)
    for m in meta.itertuples():
        full = m.coverage_remx_median >= 0.9
        ax.add_patch(plt.Rectangle((m.x - 20, m.y - 20), 40, 40, fill=False, lw=1.2 if full else 0.5,
                                   ec="black" if full else "grey", ls="-" if full else ":"))
        if m.plot_cell > 0:
            ax.text(m.x, m.y, f"P{m.plot_cell}", ha="center", va="center", fontsize=8, color="navy", weight="bold")
    p6 = meta[meta.plot_cell == 6]
    if len(p6):
        ax.add_patch(plt.Rectangle((p6.x.iloc[0] - 60, p6.y.iloc[0] - 60), 120, 120, fill=False, ec="crimson", lw=1.5))
    ax.set_xlim(ds.x.min() - 25, ds.x.max() + 25)
    ax.set_ylim(ds.y.min() - 25, ds.y.max() + 25)
    ax.set_title(f"UAV strip on the cube's 40 m cells\nNDVI {REFERENCE}; solid: cell ≥ 90 % covered; "
                 "red: the nine cells around P6", fontsize=8)
    ax.ticklabel_format(useOffset=False, style="plain")
    ax.tick_params(labelsize=6)
    ax.tick_params(axis="x", labelrotation=45)
    ax.set_xlabel("UTM 33N x (m)", fontsize=7)
    ax.set_ylabel("UTM 33N y (m)", fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    st.close()


# ------------------------------------------------------------------------------------------- README

def readme(inv, rad, cam, ndvi_cam, align, res, same, yellow, rc, qc, meta, per, tva):
    n_alt, n_rem = (inv.sensor == "Altum").sum(), (inv.sensor == "RedEdge_MX_Dual").sum()
    zones = {1: "mat", 2: "lake", 3: "grassland", 4: "other ground"}
    full = meta[meta.coverage_remx_median >= 0.9]
    p6 = meta[meta.plot_cell == 6]
    nine = meta[(meta.x.sub(p6.x.iloc[0]).abs() <= 40) & (meta.y.sub(p6.y.iloc[0]).abs() <= 40)] if len(p6) else meta[:0]
    clear = align[align.peak_ratio >= MATCH]
    within = clear[clear.against.str.startswith("same") & (clear.date.astype(str) != REFERENCE)]
    same_day = clear[clear.against.str.startswith("RedEdge")]
    unclear = align[(align.peak_ratio < MATCH) & (align.date.astype(str) != REFERENCE)]
    rcm = rc.dropna(subset=["csv_mean", "raster_mean"])
    ok_days = set(qc[(qc.sensor == "Altum") & qc.use_reflectance].date) & set(qc[(qc.sensor != "Altum") & qc.use_reflectance].date)
    cg, ng = cam[cam.date.isin(ok_days)], ndvi_cam[ndvi_cam.date.isin(ok_days)]
    ndvi_days = qc[(qc.sensor != "Altum") & qc.use_ndvi]
    t_ok = tva[tva.date.isin(set(qc[qc.use_thermal].date))]
    res_t = res[res.band == "Thermal_11um"].autocorr_half_lag_m
    res_ms = res[res.band != "Thermal_11um"].groupby("sensor").autocorr_half_lag_m
    yel = yellow.groupby("date").n_subplots_flagged.sum()

    def yn(v):
        return "yes" if v else "**no**"

    lines = [
        "# C-065 — UAV orthomosaics: what came, and what can be used",
        "",
        "Generated by `05_code/SAr-adapted/experiments/uav_rasters_c065.py`; do not edit. No radar value is read here.",
        "Data: `06_data/field/raw/2026-10_uav_rasters/` (as delivered); 0.5 m stack `06_data/field/processed/uav_0p5m.nc`"
        " (fills masked, values as delivered).",
        "",
        "## What was delivered",
        "",
        f"- {n_alt} Altum flights (5 reflectance bands + thermal in °C) and {n_rem} RedEdge-MX Dual flights (10"
        f" reflectance bands), {inv.date.min()} to {inv.date.max()}, {inv.size_gb.sum():.0f} GB of rasters (plus their"
        " overview files), with the field team's"
        " per-plot CSVs and an overlay check. They are the flight days of the first delivery's per-plot values"
        " (`field_uav_x049/uav_availability.csv`); the days that had only LAI have no raster.",
        f"- **One strip along the boardwalk**, from the station in the south to P8–P9 in the north: {len(meta)} cube"
        f" cells touched, {len(full)} covered at ≥ 90 % (zone: "
        + ", ".join(f"{zones[k]} {v}" for k, v in full.zone.value_counts().items()) + ")."
        f" The nine cells around P6: {int((nine.coverage_remx_median >= 0.9).sum())} of {len(nine)} covered at ≥ 90 %."
        " Not the whole mat.",
        f"- Pixel {inv.pixel_e_m.median() * 100:.1f} cm east × {inv.pixel_n_m.median() * 100:.1f} cm north, on a"
        " geographic grid (EPSG:4326; for some files GDAL takes it from the `.aux.xml` beside the GeoTIFF).",
        "- **Reflectance and temperature only: no DSM, no point cloud, no true-colour mosaic.** C-065's first test"
        " (does the drone surface follow the laser) cannot be done with this delivery.",
        "",
        "## Can each flight be used? (`qc_flags.csv`)",
        "",
        "Rules (set at intake after a first look at the band medians; the failures are off by factors of two or more):",
        "",
        *[f"- `{k}`: {v}" for k, v in RULES.items()],
        "",
        f"Absolute reflectance usable on {int(qc.use_reflectance.sum())} of {len(qc)} files; NDVI usable on"
        f" {len(ndvi_days)} of {n_rem} RedEdge-MX flights; surface temperature on {int(qc.use_thermal.sum())} of"
        f" {n_alt} Altum flights. Borderline (every shared band within a factor 1.6 of the other camera, so the"
        " thresholds decide): "
        + (", ".join(f"{d}" for d, g in cam.groupby("date") if not g.ratio_altum_over_remx.between(0.67, 1.5).all()
                     and np.abs(np.log(g.ratio_altum_over_remx)).max() < np.log(1.6)) or "none") + ".",
        "",
        "| Camera | Date | Reflectance | NDVI | Temperature | Why |",
        "|---|---|---|---|---|---|",
        *[f"| {'Altum' if q.sensor == 'Altum' else 'RedEdge-MX'} | {q.date} | {yn(q.use_reflectance)} |"
          f" {yn(q.use_ndvi)} | {yn(q.use_thermal) if q.sensor == 'Altum' else '—'} |"
          f" {'; '.join(x for x in (q.reflectance_flags, q.thermal_flags) if x) or ''} |" for q in qc.itertuples()],
        "",
        "## The two cameras on the same day (`cross_camera.csv`, `cross_camera_ndvi.csv`)",
        "",
        f"Both flew on {ndvi_cam.date.nunique()} days. On the {len(ok_days)} days where both pass, Altum ÷ RedEdge-MX"
        " median by shared band:",
        "",
        "| Band | median ratio | range | median pixel r (0.5 m) |", "|---|---|---|---|",
        *[f"| {b} | {g.ratio_altum_over_remx.median():.2f} | {g.ratio_altum_over_remx.min():.2f}–"
          f"{g.ratio_altum_over_remx.max():.2f} | {g.pixel_r.median():.2f} |" for b, g in cg.groupby("band", sort=False)],
        "",
        f"NDVI from the two cameras on those days: pixel r {ng.ndvi_pixel_r.min():.2f}–{ng.ndvi_pixel_r.max():.2f},"
        f" median absolute difference {ng.ndvi_median_abs_diff.median():.3f}.",
        "",
        "## Do the flights line up? (`alignment.csv`)",
        "",
        f"Shift of each flight from the {REFERENCE} flight of the same camera, at 5 cm, in three 60 m windows along the"
        f" boardwalk, by phase correlation of the image edges; only clear matches (peak ratio ≥ {MATCH:g}) count:",
        "",
        f"- {len(within)} clear matches between flights of one camera: median shift {within.shift_m.median():.2f} m,"
        f" 90th percentile {within.shift_m.quantile(0.9):.2f} m, largest {within.shift_m.max():.2f} m;",
        f"- the two cameras on the same day ({len(same_day)} clear matches): median {same_day.shift_m.median():.2f} m,"
        f" largest {same_day.shift_m.max():.2f} m;",
        f"- {len(unclear)} comparisons found no clear match ("
        + ", ".join(f"{w} {n}" for w, n in unclear.window.value_counts().items()) + "): the south window, around the"
        " station, changes between flights; the middle and north windows (platforms and boardwalk) match.",
        "",
        "## How sharp is each band? (`effective_resolution.csv`)",
        "",
        f"Distance at which the image stops resembling itself (autocorrelation 0.5), over open mat on {REFERENCE}:"
        f" Altum multispectral {res_ms.min().get('Altum', np.nan):.2f}–{res_ms.max().get('Altum', np.nan):.2f} m,"
        f" RedEdge-MX {res_ms.min().get('RedEdge_MX_Dual', np.nan):.2f}–{res_ms.max().get('RedEdge_MX_Dual', np.nan):.2f} m,"
        f" **Altum thermal {res_t.iloc[0]:.1f} m**. This mixes two things: the Altum's thermal camera has far coarser"
        " pixels than its multispectral ones (the band is delivered resampled on the 2 cm grid), and surface"
        " temperature may vary smoothly over open mat. Either way a sub-plot circle half a metre across does not"
        " resolve its own temperature; it carries that of its surroundings.",
        "",
        "## The surface temperature against the air (`thermal_vs_air.csv`)",
        "",
        f"On the usable flights the strip's median surface temperature is {t_ok.surface_minus_air_max_c.min():+.1f} to"
        f" {t_ok.surface_minus_air_max_c.max():+.1f} °C from that day's maximum air temperature at the station (flight"
        " times are not in the delivery, so the comparison is with the daily maximum).",
        "",
        "## The field team's per-plot values (`plot_tables_same.csv`, `plot_recomputed.csv`, `yellow_flags.csv`)",
        "",
        f"- The first delivery's merged table (used by X-049) holds exactly these per-plot values: largest difference"
        f" {same.max_abs_diff.max():.1g} over {int(same.n_joined.sum())} values. For the 2022 flights, P8's values are"
        " those of its 2022 position (`P8_k_2022` in the CSVs; P8 was moved in 2023).",
        f"- Recomputed from the rasters inside the P5 and P6 sub-plot circles of `Plots_PlanetScope.shp`: largest"
        f" relative difference {rcm.rel_diff.abs().max():.0e} ({len(rcm)} values) — the CSVs are these rasters, as"
        " delivered. The circles of the other sub-plots are not in the shapefile we hold.",
        "- So the values X-049 used carry the flights' faults above, except where the field team had marked the cell"
        " yellow (X-049 dropped those). Yellow cells by flight day: "
        + ", ".join(f"{d} {n}" for d, n in yel.items()) + ".",
        "",
        "## Per-cell table for later questions (`cells_meta.csv`, `cells_uav.csv`)",
        "",
        f"Per covered 40 m cell (≥ 90 %) and flight: median and IQR of NDVI, median NDRE (RedEdge-MX), median and IQR"
        f" of surface temperature (Altum), with the flight's usability (`qc_ok`); {len(per)} rows. Any comparison"
        " with radar is a new question, pre-registered first.",
        "",
        "## Files",
        "",
        "| File | Content |", "|---|---|",
        "| `inventory.csv` | per file: sensor, date, bands, size, pixel size, CRS, fill share, valid share |",
        "| `radiometry.csv` | per file and band (fill excluded): median, 1st/99th percentile, clipped share, share above 1 (thermal: above 70 °C) |",
        "| `qc_flags.csv` | per file: the rules that failed; `use_reflectance`, `use_ndvi`, `use_thermal` |",
        "| `cross_camera.csv`, `cross_camera_ndvi.csv` | same-day Altum vs RedEdge-MX |",
        "| `alignment.csv` | shift per flight and window, with the match's peak ratio |",
        "| `effective_resolution.csv` | autocorrelation half-distance per band |",
        "| `thermal_vs_air.csv` | strip surface temperature vs the day's air temperature |",
        "| `plot_tables_same.csv`, `plot_recomputed.csv`, `yellow_flags.csv` | the per-plot tables checked |",
        "| `cells_meta.csv`, `cells_uav.csv` | the 40 m cells covered and their per-flight UAV values |",
        "| `fig_quicklook.png`, `fig_coverage.png` | NDVI and temperature per flight; the strip on the cube's cells |",
    ]
    (OUT / "README.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    files = field.uav_rasters()
    inv, rad = inventory_and_radiometry(files)
    inv.to_csv(OUT / "inventory.csv", index=False)
    rad.to_csv(OUT / "radiometry.csv", index=False)
    ds = build_stack(files)
    STACK.parent.mkdir(parents=True, exist_ok=True)
    ds.to_netcdf(STACK, encoding={v: {"zlib": True, "complevel": 4, "dtype": "float32"} for v in ds.data_vars})
    cam, ndvi_cam = cross_camera(ds)
    cam.to_csv(OUT / "cross_camera.csv", index=False)
    ndvi_cam.to_csv(OUT / "cross_camera_ndvi.csv", index=False)
    align = alignment(files)
    align.to_csv(OUT / "alignment.csv", index=False)
    res = effective_resolution(files)
    res.to_csv(OUT / "effective_resolution.csv", index=False)
    same, yellow, rc = plot_tables(files)
    same.to_csv(OUT / "plot_tables_same.csv", index=False)
    yellow.to_csv(OUT / "yellow_flags.csv", index=False)
    rc.to_csv(OUT / "plot_recomputed.csv", index=False)
    qc = flags(inv, rad, cam, ndvi_cam, align)
    qc.to_csv(OUT / "qc_flags.csv", index=False)
    thermal_vs_air(rad).to_csv(OUT / "thermal_vs_air.csv", index=False)
    meta, per = cells(ds, qc)
    meta.to_csv(OUT / "cells_meta.csv", index=False)
    per.to_csv(OUT / "cells_uav.csv", index=False)
    fig_quicklook(ds, qc, OUT / "fig_quicklook.png")
    fig_coverage(ds, meta, OUT / "fig_coverage.png")
    readme(inv, rad, cam, ndvi_cam, align, res, same, yellow, rc, qc, meta, per, thermal_vs_air(rad))
    print(f"wrote {OUT} and {STACK}")


if __name__ == "__main__":
    main()

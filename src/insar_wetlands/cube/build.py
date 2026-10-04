"""M2 — the harmonised (silver) layer: one space (40 m template; 2026 also at its native 20 m), one time base (UTC),
one set of units; flags as columns. Called by ``scripts/cube_build.py harmonise <what>``.

| product | content |
|---|---|
| ``ifg_<period>_<track>.nc`` | every HyP3 pair on disk: wrapped / unwrapped phase (rad, unreferenced), coherence, connected component; look-vector angles and DEM |
| ``ifg_2026_<track>_20m.nc`` | the same at 20 m (D-024) |
| ``rtc_<track>.nc`` | γ⁰ VV, VH (dB), VH/VV, RVI, platform and exact time, 2020 → |
| ``s2_40m.nc``, ``s2_plots.csv.gz`` | NDVI, NDRE, NDMI, NDWI, MNDWI (mean, SD, valid share), snow share, every scene |
| ``lst_40m.nc``, ``lst_plots.csv.gz`` | Landsat 8/9 surface temperature (°C), clear pixels |
| ``static_40m.nc``, ``static_20m.nc`` | zones, stable ground, hard targets, distance to the mat edge, WorldCover shares, LiDAR terrain and vegetation, boardwalk share, geometry, mean radar behaviour, spatial CV blocks |
| ``lidar_1m.nc`` | the LiDAR DTM and DSM (EVRF2007) at 1 m over the peatland |
| ``site_hourly.csv.gz`` | station meteo, water table P1–P9, laser, ERA5-Land / ERA5, IMGW, derived dew and wetness variables, all flags |
| ``plot_units.csv``, ``folds.json`` | plots, their cells and radar units; the cross-validation schemes |
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from scipy import ndimage

from .. import field
from ..paths import make_paths
from ..regrid import stack_to_template
from ..stack import list_pairs, load_layer, load_static_layer
from . import core
from .sources import indices_from_npz as src_indices

HUB = Path(__file__).resolve().parents[5]
LOCAL = HUB / "05_code" / "local"
MIRROR = LOCAL / "drive_mirror"
DLV = HUB / "08_deliverables"
REPO = Path(__file__).resolve().parents[3]
TRACKS = ("ascending", "descending")
PERIODS = ("2020_2021", "2022_2024", "2025", "2026")
RTC_FILES = {"ascending": ["rtc_dualpol_2020_2024.nc", "rtc_dualpol_2025_2025.nc", "rtc_dualpol_2026_2026.nc"],
             "descending": ["rtc_dualpol_2020_2024_descending.nc", "rtc_dualpol_2025_2025_descending.nc",
                            "rtc_dualpol_2026_2026_descending.nc"]}
S2_INDICES = ("ndvi", "ndre", "ndmi", "ndwi", "mndwi")
ZONE_CODES = {"A": 1, "B": 2, "C": 3, "D": 4}           # mat, lake, grassland, other ground
WC_CLASSES = {10: "tree", 20: "shrub", 30: "grass", 40: "crop", 50: "built", 60: "bare", 80: "water",
              90: "herb_wetland", 100: "moss"}
LIDAR_TRANSECT_PAD_M = 60.0
SPATIAL_BLOCK_CELLS = 5                                   # 200 m blocks at 40 m
INT16 = {"dtype": "int16", "scale_factor": 1e-4, "_FillValue": -32768, "zlib": True, "complevel": 4}
F32 = {"dtype": "float32", "zlib": True, "complevel": 4}


def ifg_root(ctx, period: str, track: str) -> Path:
    if period == "2022_2024":
        return make_paths(cfg=ctx.cfg, track=track).cropped
    return LOCAL / f"s1_{period}" / f"hyp3_cropped_{track}"


def _enc(ds: xr.Dataset, spec: dict | None = None, ints=()) -> dict:
    out = {}
    for v in ds.data_vars:
        if v in ints:
            out[v] = {"zlib": True, "complevel": 4}
        elif ds[v].dtype.kind == "f":
            out[v] = dict(spec or F32)
    return out


LIDAR_GRID_BUFFER_M = 640.0                              # a multiple of 40 m: 1 m blocks align with the template


def peat_bounds_utm(ctx, buffer_m: float, zones: str = "ABC") -> tuple[float, float, float, float]:
    """Bounding box (cell edges) of the given zones (A mat, B lake, C grassland) plus a buffer."""
    z = ctx.zones
    m = np.zeros(z["A"].shape, bool)
    for k in zones:
        m |= z[k].values.astype(bool)
    tx, ty = tpl_xy(ctx)
    rr, cc = np.nonzero(m)
    return (tx[cc].min() - 20 - buffer_m, ty[rr].min() - 20 - buffer_m, tx[cc].max() + 20 + buffer_m,
            ty[rr].max() + 20 + buffer_m)


def tpl_xy(ctx):
    return ctx.template.x.values.astype(float), ctx.template.y.values.astype(float)


def grid20(ctx, silver: Path | None = None):
    """Centres of the 20 m grid: the 2026 crops' own grid when harmonised (its edges fall on the template's, one half
    40 m cell inside it), else the 2 × 2 split of the template."""
    for track in TRACKS:
        f = silver / f"ifg_2026_{track}_20m.nc" if silver else None
        if f is not None and f.exists():
            d = xr.open_dataset(f)
            return d.x.values.astype(float), d.y.values.astype(float)
    tx, ty = tpl_xy(ctx)
    x = np.concatenate([[c - 10, c + 10] for c in tx])
    y = np.concatenate([[c + 10, c - 10] for c in ty])          # template y decreases
    return x, y


def plot_units() -> pd.DataFrame:
    """Plots, their 40 m cells and radar units (P8 and P9 share a cell — one unit)."""
    px = pd.read_csv(DLV / "field_first" / "plot_pixels.csv")
    px["unit"] = px.groupby(["row", "col"])["plot"].transform(lambda s: "+".join(sorted(s)))
    return px


# ------------------------------------------------------------------------------ interferograms

def harmonise_ifg(ctx, silver: Path, log=print) -> None:
    tx, ty = tpl_xy(ctx)
    plat = platform_by_date()
    for period in PERIODS:
        for track in TRACKS:
            root = ifg_root(ctx, period, track)
            out = silver / f"ifg_{period}_{track}.nc"
            if out.exists() or not root.exists():
                continue
            pairs = list_pairs(root)
            log(f"ifg {period} {track}: {len(pairs)} pairs")
            lay = {name: load_layer(root, name, pairs) for name in ("wrapped_phase", "unw_phase", "corr", "conncomp")}
            native = lay["corr"].sizes["x"] != len(tx)
            st = {name: load_static_layer(root, name) for name in ("lv_theta", "lv_phi", "dem")}
            if native:
                _write_ifg(lay, st, plat, track, silver / f"ifg_{period}_{track}_20m.nc", log)
                lay = {"wrapped_phase": stack_to_template(lay["wrapped_phase"], tx, ty),
                       "unw_phase": stack_to_template(lay["unw_phase"], tx, ty),
                       "corr": stack_to_template(lay["corr"], tx, ty),
                       "conncomp": _cc_to_template(lay["conncomp"], tx, ty)}
                st = {k: _static_to_template(v, tx, ty) for k, v in st.items()}
            else:
                lay = {k: stack_to_template(v, tx, ty) if k != "conncomp" else v for k, v in lay.items()}
                st = {k: _static_to_template(v, tx, ty) for k, v in st.items()}
            _write_ifg(lay, st, plat, track, out, log)


def _static_to_template(da, tx, ty):
    from ..regrid import to_template
    if da.sizes["x"] == len(tx) and np.allclose(da.x.values, tx):
        return da.values.astype("float32")
    return to_template(da.values, da.x.values, da.y.values, tx, ty, kind="mean").astype("float32")


def _cc_to_template(cc: xr.DataArray, tx, ty) -> xr.DataArray:
    """Connected component on the template: 1 where most fine pixels of the cell were unwrapped (component > 0)."""
    vals = np.stack([core.aggregate_to_grid((c > 0).astype(float), cc.x.values, cc.y.values, tx, ty)["mean"] >= 0.5
                     for c in cc.values]).astype("int16")
    return xr.DataArray(vals, dims=("pair", "y", "x"), coords={"pair": cc.pair.values, "y": ty, "x": tx})


def _write_ifg(lay: dict, st: dict, plat: pd.Series, track: str, out: Path, log) -> None:
    w = lay["wrapped_phase"]
    pairs = [str(p) for p in w.pair.values]
    t1 = pd.to_datetime([p[:8] for p in pairs])
    t2 = pd.to_datetime([p[9:17] for p in pairs])
    x, y = w.x.values, w.y.values
    ds = xr.Dataset({"wrapped": (("pair", "y", "x"), w.values.astype("float32")),
                     "unw": (("pair", "y", "x"), lay["unw_phase"].values.astype("float32")),
                     "coh": (("pair", "y", "x"), lay["corr"].values.astype("float32")),
                     "conncomp": (("pair", "y", "x"), np.nan_to_num(lay["conncomp"].values, nan=0).astype("int16")),
                     "lv_theta": (("y", "x"), np.asarray(st["lv_theta"], "float32")),
                     "lv_phi": (("y", "x"), np.asarray(st["lv_phi"], "float32")),
                     "dem": (("y", "x"), np.asarray(st["dem"], "float32"))},
                    coords={"pair": pairs, "x": x, "y": y,
                            "t1": ("pair", [core.overpass_time(d, track) for d in t1]),
                            "t2": ("pair", [core.overpass_time(d, track) for d in t2]),
                            "revisit_days": ("pair", (t2 - t1).days.to_numpy()),
                            "platforms": ("pair", [f"{plat.get((track, a.date()), '?')}/{plat.get((track, b.date()), '?')}"
                                                   for a, b in zip(t1, t2)])})
    for c in ("t1", "t2"):
        ds[c] = ds[c].astype("datetime64[ns]")
    ds.attrs = {"phase_units": "rad, HyP3 sign, unreferenced (per-pair constants as delivered)", "track": track,
                "crs": "EPSG:32633", "wavelength_mm": core.WAVELENGTH_MM}
    ds.to_netcdf(out, encoding=_enc(ds, ints=("conncomp",)))
    log(f"  wrote {out.name} {dict(ds.sizes)}")


# ------------------------------------------------------------------------------ backscatter

def platform_by_date() -> pd.Series:
    """(track, date) → platform (S1A/S1B/S1C/S1D) from the RTC source items, the only place the platform is recorded."""
    rows = []
    for track, files in RTC_FILES.items():
        for f in files:
            p = MIRROR / f
            if p.exists():
                ds = xr.open_dataset(p)
                for item in ds.source_item.values.astype(str):
                    rows.append({"track": track, "date": pd.Timestamp(item.split("_")[4][:8]).date(), "platform": item[:3]})
    d = pd.DataFrame(rows).drop_duplicates(["track", "date"])
    return d.set_index(["track", "date"]).platform


def harmonise_rtc(ctx, silver: Path, log=print) -> None:
    for track, files in RTC_FILES.items():
        out = silver / f"rtc_{track}.nc"
        parts = [xr.open_dataset(MIRROR / f) for f in files if (MIRROR / f).exists()]
        ds = xr.concat(parts, "time").sortby("time")
        ds = ds.isel(time=~pd.Index(ds.time.values).duplicated())
        items = ds.source_item.values.astype(str)
        acq = pd.to_datetime([i.split("_")[4] for i in items], format="%Y%m%dT%H%M%S", utc=True)
        vv, vh = 10 ** (ds.gamma0_vv_db / 10), 10 ** (ds.gamma0_vh_db / 10)
        o = xr.Dataset({"vv_db": ds.gamma0_vv_db.astype("float32"), "vh_db": ds.gamma0_vh_db.astype("float32"),
                        "vh_vv_db": ds.ratio_vh_vv_db.astype("float32"),
                        "rvi": (4 * vh / (vv + vh)).astype("float32")},
                       coords={"time": ds.time.values, "x": ds.x.values, "y": ds.y.values,
                               "acq_time_utc": ("time", acq.tz_convert(None)), "platform": ("time", [i[:3] for i in items])})
        o.attrs = {"track": track, "crs": "EPSG:32633", "rvi": "4·σVH/(σVV+σVH), linear"}
        o.to_netcdf(out, encoding=_enc(o))
        log(f"rtc {track}: {o.sizes['time']} dates → {out.name}")


# ------------------------------------------------------------------------------ optical

def _grid_from_bounds(left, top, res, w, h):
    return left + res / 2 + res * np.arange(w), top - res / 2 - res * np.arange(h)


def _plot_rc(px: pd.DataFrame, gx, gy, res):
    left, top = gx[0] - res / 2, gy[0] + res / 2
    return (np.floor((top - px.N) / res).astype(int), np.floor((px.E - left) / res).astype(int))


def harmonise_optical(ctx, raw: Path, silver: Path, log=print, min_valid: float = 0.02) -> None:
    tx, ty = tpl_xy(ctx)
    left, top = tx[0] - 20, ty[0] + 20
    px = plot_units()
    # Sentinel-2 (20 m)
    sc = pd.read_csv(raw / "s2" / "scenes.csv", parse_dates=["time"])
    sc = sc[(sc.file.fillna("") != "") & (sc.valid_share >= min_valid)].sort_values("time")
    sc = sc.drop_duplicates("time")
    first = np.load(raw / "s2" / sc.file.iloc[0])
    h, w = first["ndvi"].shape
    gx, gy = _grid_from_bounds(left, top, 20.0, w, h)
    rr, cc = _plot_rc(px, gx, gy, 20.0)
    stack = {f"{k}{s}": [] for k in S2_INDICES for s in ("", "_sd")}
    stack.update(valid_frac=[], snow_frac=[])
    prow = []
    for f, t in zip(sc.file, sc.time):
        z = np.load(raw / "s2" / f)
        if "refl_B04" in z.files:
            z = src_indices(z)
        for k in S2_INDICES:
            a = core.aggregate_to_grid(z[k], gx, gy, tx, ty)
            stack[k].append(a["mean"])
            stack[f"{k}_sd"].append(a["sd"])
        stack["valid_frac"].append(core.aggregate_to_grid(z["valid"], gx, gy, tx, ty)["mean"])
        stack["snow_frac"].append(core.aggregate_to_grid(z["snow"], gx, gy, tx, ty)["mean"])
        for p, r, c in zip(px["plot"], rr, cc):
            row = {"time_utc": t, "plot": p}
            for k in S2_INDICES:
                row[f"{k}_1x1"] = float(z[k][r, c])
                row[f"{k}_3x3"] = core.window_stats(z[k], r, c)["mean"]
            row["valid_3x3"] = float(np.nanmean(z["valid"][max(r - 1, 0):r + 2, max(c - 1, 0):c + 2]))
            row["snow_3x3"] = float(np.nanmean(z["snow"][max(r - 1, 0):r + 2, max(c - 1, 0):c + 2]))
            prow.append(row)
    times = pd.DatetimeIndex(sc.time).tz_convert(None) if sc.time.dt.tz is not None else pd.DatetimeIndex(sc.time)
    ds = xr.Dataset({k: (("time", "y", "x"), np.stack(v).astype("float32")) for k, v in stack.items()},
                    coords={"time": times, "x": tx, "y": ty})
    ds.attrs = {"source": "Sentinel-2 L2A collection 1 (Earth Search), SCL-masked; 20 m → 40 m mean",
                "indices": "ndvi B08/B04, ndre B8A/B05, ndmi B8A/B11, ndwi B03/B08, mndwi B03/B11"}
    ds.to_netcdf(silver / "s2_40m.nc", encoding=_enc(ds, INT16))
    pd.DataFrame(prow).to_csv(silver / "s2_plots.csv.gz", index=False)
    log(f"S2: {ds.sizes['time']} scenes → s2_40m.nc, s2_plots.csv.gz")
    # Landsat surface temperature (30 m)
    f = raw / "landsat" / "scenes.csv"
    if f.exists():
        ls = pd.read_csv(f, parse_dates=["time"])
        ls = ls[(ls.file.fillna("") != "") & (ls.clear_share > min_valid)].sort_values("time").drop_duplicates("time")
        z0 = np.load(raw / "landsat" / ls.file.iloc[0])["lst_c"]
        lx, ly = _grid_from_bounds(left, top, 30.0, z0.shape[1], z0.shape[0])
        lr, lc = _plot_rc(px, lx, ly, 30.0)
        lst, prow = [], []
        for fn, t in zip(ls.file, ls.time):
            z = np.load(raw / "landsat" / fn)["lst_c"]
            lst.append(core.aggregate_to_grid(z, lx, ly, tx, ty)["mean"])
            for p, r, c in zip(px["plot"], lr, lc):
                prow.append({"time_utc": t, "plot": p, "lst_c_1x1": float(z[r, c]),
                             "lst_c_3x3": core.window_stats(z, r, c)["mean"]})
        lt = pd.DatetimeIndex(ls.time).tz_convert(None) if ls.time.dt.tz is not None else pd.DatetimeIndex(ls.time)
        d2 = xr.Dataset({"lst_c": (("time", "y", "x"), np.stack(lst).astype("float32"))},
                        coords={"time": lt, "x": tx, "y": ty})
        d2.attrs = {"source": "Landsat 8/9 C2 L2 ST_B10 (Planetary Computer), QA-clear; 30 m → 40 m mean"}
        d2.to_netcdf(silver / "lst_40m.nc", encoding=_enc(d2))
        pd.DataFrame(prow).to_csv(silver / "lst_plots.csv.gz", index=False)
        log(f"Landsat: {d2.sizes['time']} scenes → lst_40m.nc")


def harmonise_ecostress(ctx, raw: Path, silver: Path, log=print, min_valid: float = 0.02) -> None:
    """ECOSTRESS surface temperature (read on the 40 m template already) → ``eco_40m.nc`` and ``eco_plots.csv.gz``; the
    acquisition hour is kept — ECOSTRESS samples every time of day, the radar's dawn and dusk included."""
    f = raw / "ecostress" / "scenes.csv"
    if not f.exists():
        return
    tx, ty = tpl_xy(ctx)
    sc = pd.read_csv(f, parse_dates=["time"])
    sc = sc[(sc.file.fillna("") != "") & (sc.valid_share >= min_valid)].sort_values("time").drop_duplicates("time")
    px = plot_units()
    lst, prow = [], []
    for fn, t in zip(sc.file, sc.time):
        z = np.load(raw / "ecostress" / fn)["lst_c"]
        if z.shape != (len(ty), len(tx)):
            continue
        lst.append(z)
        for p, r, c in zip(px["plot"], px.row, px.col):
            prow.append({"time_utc": t, "plot": p, "eco_lst_c_1x1": float(z[r, c]),
                         "eco_lst_c_3x3": core.window_stats(z, r, c)["mean"]})
    times = pd.DatetimeIndex(sc.time.iloc[:len(lst)])
    times = times.tz_convert(None) if times.tz is not None else times
    ds = xr.Dataset({"lst_c": (("time", "y", "x"), np.stack(lst).astype("float32"))},
                    coords={"time": times, "x": tx, "y": ty, "hour_utc": ("time", times.hour + times.minute / 60)})
    ds.attrs = {"source": "ECOSTRESS L2T LSTE v002 (LP DAAC), cloud-masked; 70 m → 40 m template (nearest)"}
    ds.to_netcdf(silver / "eco_40m.nc", encoding=_enc(ds))
    pd.DataFrame(prow).to_csv(silver / "eco_plots.csv.gz", index=False)
    log(f"ECOSTRESS: {ds.sizes['time']} scenes → eco_40m.nc, eco_plots.csv.gz")


# ------------------------------------------------------------------------------ static layers

def _lidar_1m(raw: Path, bounds_utm, log=print):
    """DTM and DSM tiles (EPSG:2180) → one 1 m UTM grid over ``bounds_utm`` (edges on whole metres)."""
    import rasterio
    from rasterio.merge import merge
    from rasterio.transform import from_origin
    from rasterio.warp import Resampling, reproject

    left, bottom, right, top = bounds_utm
    W, H = int(right - left), int(top - bottom)
    dst_t = from_origin(left, top, 1.0, 1.0)
    out = {}
    from pyproj import Transformer
    tr2180 = Transformer.from_crs(32633, 2180, always_xy=True)
    xs, ys = tr2180.transform([left, right, left, right], [bottom, bottom, top, top])
    win = (min(xs) - 50, min(ys) - 50, max(xs) + 50, max(ys) + 50)

    def touches(f: Path) -> bool:                # tiles are named <kind>_<x0>_<y0> (500 m, EPSG:2180)
        x0, y0 = (float(v) for v in f.stem.split("_")[1:3])
        return x0 < win[2] and x0 + 500 > win[0] and y0 < win[3] and y0 + 500 > win[1]

    for kind in ("dtm", "dsm"):
        files = [f for f in sorted((raw / "lidar").glob(f"{kind}_*.tif")) if touches(f)]
        srcs = [rasterio.open(f) for f in files]
        nod = srcs[0].nodata
        mos, tr = merge(srcs, nodata=nod)
        crs = srcs[0].crs
        for s in srcs:
            s.close()
        a = mos[0].astype("float32")
        bad = ~np.isfinite(a) | (a <= 0) | (a < -100) | (a > 1000)       # nodata and fill values seen in the WCS
        a[bad] = np.nan
        dst = np.full((H, W), np.nan, "float32")
        reproject(a, dst, src_transform=tr, src_crs=crs, dst_transform=dst_t, dst_crs="EPSG:32633",
                  resampling=Resampling.average, src_nodata=np.nan, dst_nodata=np.nan)
        out[kind] = dst
        log(f"LiDAR {kind}: {len(files)} tiles → {H}×{W} m, valid {np.isfinite(dst).mean():.2f}")
    x = left + 0.5 + np.arange(W)
    y = top - 0.5 - np.arange(H)
    return out["dtm"], out["dsm"], x, y


def _block(a, k):
    h, w = a.shape[0] // k, a.shape[1] // k
    return a[:h * k, :w * k].reshape(h, k, w, k).transpose(0, 2, 1, 3).reshape(h, w, k * k)


def lidar_features(dtm, dsm, k: int) -> dict[str, np.ndarray]:
    """Per k-metre cell from 1 m LiDAR: terrain mean, roughness (SD), micro-relief SD (terrain minus its 25 m Gaussian
    smooth — the hummock-hollow scale), slope, vegetation height (DSM − DTM) mean and 90th percentile, valid share."""
    from scipy.ndimage import gaussian_filter
    filled = np.where(np.isfinite(dtm), dtm, np.nanmedian(dtm))
    wts = gaussian_filter(np.isfinite(dtm).astype(float), 25)
    smooth = gaussian_filter(np.where(np.isfinite(dtm), dtm, 0.0), 25) / np.maximum(wts, 1e-6)
    micro = np.where(np.isfinite(dtm), dtm - smooth, np.nan)
    gy, gx = np.gradient(filled)
    slope = np.where(np.isfinite(dtm), np.degrees(np.arctan(np.hypot(gx, gy))), np.nan)
    vh = np.where(np.isfinite(dsm) & np.isfinite(dtm), np.clip(dsm - dtm, 0, 60), np.nan)
    with np.errstate(invalid="ignore"), __import__("warnings").catch_warnings():
        __import__("warnings").simplefilter("ignore", RuntimeWarning)
        B = {n: _block(v, k) for n, v in (("dtm", dtm), ("micro", micro), ("slope", slope), ("vh", vh))}
        return {"lidar_dtm_m": np.nanmean(B["dtm"], 2), "lidar_roughness_m": np.nanstd(B["dtm"], 2),
                "lidar_microrelief_sd_m": np.nanstd(B["micro"], 2), "lidar_slope_deg": np.nanmean(B["slope"], 2),
                "lidar_vegh_mean_m": np.nanmean(B["vh"], 2), "lidar_vegh_p90_m": np.nanpercentile(B["vh"], 90, axis=2),
                "lidar_valid_frac": np.isfinite(B["dtm"]).mean(2)}


def _to20(a: np.ndarray, tx, ty, x20, y20) -> np.ndarray:
    """A 40 m template layer onto the 20 m grid (nearest: each 20 m cell takes its 40 m parent)."""
    c = np.clip(np.floor((x20 - (tx[0] - 20)) / 40).astype(int), 0, len(tx) - 1)
    r = np.clip(np.floor(((ty[0] + 20) - y20) / 40).astype(int), 0, len(ty) - 1)
    return np.asarray(a)[np.ix_(r, c)]


def _paste(sub: np.ndarray, full_shape, r0: int, c0: int) -> np.ndarray:
    """``sub`` placed with its top-left cell at (r0, c0) of a NaN array; parts falling outside are dropped."""
    out = np.full(full_shape, np.nan, "float32")
    sr, sc = max(-r0, 0), max(-c0, 0)
    r0, c0 = max(r0, 0), max(c0, 0)
    h = min(sub.shape[0] - sr, full_shape[0] - r0)
    w = min(sub.shape[1] - sc, full_shape[1] - c0)
    if h > 0 and w > 0:
        out[r0:r0 + h, c0:c0 + w] = sub[sr:sr + h, sc:sc + w]
    return out


def _tif_on_grid(path, gx, gy):
    import rioxarray  # noqa: F401
    da = xr.open_dataarray(path, engine="rasterio").squeeze(drop=True)
    from ..regrid import to_template
    if da.sizes["x"] == len(gx) and np.allclose(da.x.values, gx):
        return da.values.astype("float32")
    return to_template(da.values.astype(float), da.x.values, da.y.values, gx, gy, kind="mean").astype("float32")


def harmonise_static(ctx, raw: Path, silver: Path, peat_bounds, log=print) -> None:
    tx, ty = tpl_xy(ctx)
    x20, y20 = grid20(ctx, silver)
    Z = {k: ctx.zones[k].values.astype(bool) for k in "ABCD"}
    zone = np.zeros(Z["A"].shape, "int16")
    for k in "DCBA":                                     # mat wins over everything
        zone[Z[k]] = ZONE_CODES[k]
    out_mat = ndimage.distance_transform_edt(~Z["A"]) * 40.0
    in_mat = ndimage.distance_transform_edt(Z["A"]) * 40.0
    stable = Z["D"] & (out_mat > 200) & (out_mat < 1500)     # as X-070
    L = {"zone": zone, "stable": stable.astype("int16"), "dist_to_mat_m": out_mat.astype("float32"),
         "depth_in_mat_m": in_mat.astype("float32"),
         "spatial_block": core.spatial_blocks(zone.shape, SPATIAL_BLOCK_CELLS).astype("int32")}
    up = lambda a: _to20(a, tx, ty, x20, y20)  # noqa: E731
    L20 = {"zone": up(zone).astype("int16"), "stable": up(L["stable"]).astype("int16"),
           "spatial_block": up(L["spatial_block"]).astype("int32")}
    # WorldCover 2021 shares
    wc = sorted((raw / "worldcover").glob("*2021*.npz"))
    if wc:
        cl = np.load(wc[0])["classes"]
        h, w = cl.shape
        wx, wy = _grid_from_bounds(tx[0] - 20, ty[0] + 20, 10.0, w, h)
        for code, name in WC_CLASSES.items():
            L[f"wc_{name}"] = core.class_fractions(cl, wx, wy, tx, ty, [code])[code].astype("float32")
            L20[f"wc_{name}"] = core.class_fractions(cl, wx, wy, x20, y20, [code])[code].astype("float32")
    # boardwalk share from the 0.25 m mask (X-060)
    bw = DLV / "field_boardwalk_x060" / "boardwalk_mask.tif"
    if bw.exists():
        import rioxarray  # noqa: F401
        m = xr.open_dataarray(bw, engine="rasterio").squeeze(drop=True)
        v = (m.values > 0).astype(float)
        for g, (gx, gy) in (("40", (tx, ty)), ("20", (x20, y20))):
            a = core.aggregate_to_grid(v, m.x.values, m.y.values, gx, gy)["mean"]
            # outside the mask's footprint the boardwalk is absent, not unknown
            (L if g == "40" else L20)["boardwalk_share"] = np.nan_to_num(a, nan=0.0).astype("float32")
    # GLO-30 DEM, slope (repo layers on the template)
    for name, f in (("dem_glo30_m", "elevation_glo30_dem.tif"), ("slope_glo30_deg", "slope_topography.tif")):
        p = REPO / "data" / "gis" / f
        if p.exists():
            L[name] = _tif_on_grid(p, tx, ty)
    # LiDAR
    if (raw / "lidar").exists() and any((raw / "lidar").glob("dtm_*.tif")):
        b = peat_bounds
        dtm, dsm, lx, ly = _lidar_1m(raw, b, log)
        xr.Dataset({"dtm_m": (("y", "x"), dtm), "dsm_m": (("y", "x"), dsm)}, coords={"x": lx, "y": ly},
                   attrs={"source": "GUGiK NMT/NMPT WCS, EVRF2007 heights, 1 m (DSM 0.5 m averaged)", "crs": "EPSG:32633"}
                   ).to_netcdf(silver / "lidar_1m.nc", encoding={"dtm_m": F32, "dsm_m": F32})
        for k, LL, gx, gy in ((40, L, tx, ty), (20, L20, x20, y20)):
            f = lidar_features(dtm, dsm, k)
            shape = (len(gy), len(gx))
            rr, cc = int(round((gy[0] + k / 2 - b[3]) / k)), int(round((b[0] - (gx[0] - k / 2)) / k))
            for n, v in f.items():
                LL[n] = _paste(v, shape, rr, cc)
    # geometry and mean radar behaviour per track (2022–2024 at 40 m, 2026 at 20 m)
    for track in TRACKS:
        tr = track[:3]
        f = silver / f"ifg_2022_2024_{track}.nc"
        if f.exists():
            ds = xr.open_dataset(f)
            e, n, u = core.los_enu(ds.lv_theta.values, ds.lv_phi.values)
            L[f"incidence_{tr}_deg"] = core.incidence_deg(ds.lv_theta.values).astype("float32")
            L[f"los_e_{tr}"], L[f"los_n_{tr}"], L[f"los_u_{tr}"] = (a.astype("float32") for a in (e, n, u))
            L[f"coh12_{tr}_median"] = ds.coh.where(ds.revisit_days == 12).median("pair").values.astype("float32")
        f6 = [silver / f"ifg_{p}_{track}.nc" for p in ("2020_2021", "2025")]
        c6 = [xr.open_dataset(p).coh.where(xr.open_dataset(p).revisit_days == 6) for p in f6 if p.exists()]
        if c6:
            L[f"coh6_{tr}_median"] = xr.concat(c6, "pair").median("pair").values.astype("float32")
        f20 = silver / f"ifg_2026_{track}_20m.nc"
        if f20.exists():
            d20 = xr.open_dataset(f20)
            if d20.sizes["x"] == len(x20):
                e, n, u = core.los_enu(d20.lv_theta.values, d20.lv_phi.values)
                L20[f"incidence_{tr}_deg"] = core.incidence_deg(d20.lv_theta.values).astype("float32")
                L20[f"los_u_{tr}"] = u.astype("float32")
                L20[f"coh6_{tr}_median"] = d20.coh.where(d20.revisit_days == 6).median("pair").values.astype("float32")
        r = silver / f"rtc_{track}.nc"
        if r.exists():
            rd = xr.open_dataset(r)
            L[f"vv_{tr}_median_db"] = rd.vv_db.median("time").values.astype("float32")
            L[f"vh_{tr}_median_db"] = rd.vh_db.median("time").values.astype("float32")
    # optical summer medians
    s2 = silver / "s2_40m.nc"
    if s2.exists():
        d = xr.open_dataset(s2)
        summer = d.time.dt.month.isin([5, 6, 7, 8, 9])
        for k in ("ndvi", "ndre", "ndmi", "ndwi"):
            L[f"s2_{k}_summer_median"] = d[k].where(summer).median("time").values.astype("float32")
    # hard targets: built-up cells that stay coherent at 12 days (phase reference candidates)
    if "wc_built" in L and "coh12_asc_median" in L:
        L["hard_target"] = ((L["wc_built"] >= 0.5) & (L["coh12_asc_median"] >= 0.5)).astype("int16")
    # plots
    px = plot_units()
    pl = np.zeros(zone.shape, "int16")
    for i, (r, c) in enumerate(zip(px.row, px.col), 1):
        pl[r, c] = i
    L["plot_cell"] = pl
    ds = xr.Dataset({k: (("y", "x"), v) for k, v in L.items()}, coords={"x": tx, "y": ty})
    ds.attrs = {"zone": "1 mat (A), 2 lake (B), 3 grassland (C), 4 other ground (D)", "plot_cell": ",".join(px["plot"]),
                "crs": "EPSG:32633", "stable": "zone D, 200–1500 m from the mat (as X-070)",
                "hard_target": "WorldCover built-up ≥ 0.5 and median 12-day ascending coherence ≥ 0.5"}
    ds.to_netcdf(silver / "static_40m.nc", encoding=_enc(ds))
    d20 = xr.Dataset({k: (("y", "x"), v) for k, v in L20.items()}, coords={"x": x20, "y": y20}, attrs=ds.attrs)
    d20.to_netcdf(silver / "static_20m.nc", encoding=_enc(d20))
    log(f"static: {len(L)} layers at 40 m, {len(L20)} at 20 m; hard targets {int(L.get('hard_target', np.zeros(1)).sum())}")


# ------------------------------------------------------------------------------ the hourly site table

STATION = {"Air_2m": "st_air_c", "RH_2m": "st_rh_pct", "VPD_Kpa": "st_vpd_kpa", "PPFDg": "st_ppfd_global",
           "PPFDd": "st_ppfd_diffuse", "Rain_mm_Tot": "st_rain_mm", "WS_ms": "st_wind_ms", "WindDir": "st_wind_dir"}


def site_hourly(raw: Path, end: str, log=print) -> pd.DataFrame:
    w = field.load_wtd_hourly()
    idx = pd.date_range("2020-01-01", end, freq="1h", tz="UTC", name="time_utc")
    # the logger stamps sit at :30 UTC centres; put them on the whole hour (nearest within 31 min)
    w = w[~w.index.duplicated()]
    wh = core.asof(idx, w, "31min").drop(columns="lag_h")
    wh.index = idx
    df = pd.DataFrame(index=idx)
    for a, b in STATION.items():
        df[b] = wh[a].to_numpy()
    for p in field.PLOTS:
        df[f"wtd_{p}_cm"] = wh[p].to_numpy()
        df[f"censored_{p}"] = field.censored_flag(w[p].astype(float)).reindex(idx, method="nearest",
                                                                            tolerance=pd.Timedelta("31min")).fillna(False).astype(bool).to_numpy()
    df["cl_raw"], df["cr_raw"] = wh["CL_raw"].to_numpy(), wh["CR_raw"].to_numpy()
    mat = [f"wtd_{p}_cm" for p in ("P6", "P7", "P8", "P9")]
    df["wtd_mat_mean_cm"] = df[mat].where(~df[[f"censored_{p}" for p in ("P6", "P7", "P8", "P9")]].to_numpy()).mean(1)
    # laser (X-058 QC flags)
    lf = DLV / "field_p6" / "laser_qc" / "laser_flags.csv"
    if lf.exists():
        las = pd.read_csv(lf, parse_dates=["time_utc"]).set_index("time_utc")
        las.index = pd.DatetimeIndex(las.index, tz="UTC") if las.index.tz is None else las.index
        lh = core.asof(idx, las, "31min").drop(columns="lag_h")
        df["laser_surface_cm"] = lh["surface_cm"].to_numpy()
        df["laser_surface_raw_cm"] = lh["surface_cm_raw"].to_numpy()
        for c in ("snow_72h", "snow_24h", "outlier", "filled"):
            df[f"laser_{c}"] = lh[c].astype("boolean").to_numpy()
        # X-058's snow mask needs the station's air temperature, which ends with 2024: after that it reads "snow" by
        # construction. There, use X-066's Open-Meteo mask calibrated on the station (fusion_v2/flags_2025).
        df["laser_snow_72h_source"] = np.where(df.laser_surface_raw_cm.notna(), "station", "")
        om_snow = DLV / "fusion_v2" / "flags_2025" / "laser_snow72_2025.csv"
        st_end = w["Air_2m"].dropna().index.max()
        if om_snow.exists():
            sn = pd.read_csv(om_snow, parse_dates=["time_utc"]).set_index("time_utc")
            sn.index = pd.DatetimeIndex(sn.index, tz="UTC") if sn.index.tz is None else sn.index
            v = core.asof(idx, sn, "31min")["snow_72h_open_meteo"]
            after = (idx > st_end) & v.notna().to_numpy()
            df.loc[after, "laser_snow_72h"] = v[after].astype(bool).to_numpy()
            df.loc[after & df.laser_surface_raw_cm.notna().to_numpy(), "laser_snow_72h_source"] = "open_meteo_calibrated"
        df["laser_ok"] = df.laser_surface_cm.notna() & ~df[["laser_snow_72h", "laser_outlier", "laser_filled"]].fillna(True).any(axis=1)
    # P6 water level on the laser's datum — meaningful only if the logger's WTD is depth below the moving surface
    df["p6_level_on_laser_datum_cm"] = df.get("laser_surface_cm", np.nan) + df["wtd_P6_cm"]
    # Open-Meteo (ERA5-Land, ERA5)
    om = pd.read_csv(raw / "open_meteo" / "hourly.csv.gz", index_col=0, parse_dates=True)
    om.index = pd.DatetimeIndex(om.index, tz="UTC") if om.index.tz is None else om.index
    om = om.reindex(idx)
    for c in om.columns:
        df[f"om_{c}"] = om[c].to_numpy()
    # IMGW daily, on every hour of its date
    imgw = raw / "imgw" / "daily.csv"
    if imgw.exists():
        d = pd.read_csv(imgw, parse_dates=["date"])
        day = idx.tz_convert(None).normalize()
        for s, g in d.groupby("station"):
            g = g.set_index("date")
            for c in ("snow_depth_cm", "dew_h", "hoarfrost_h", "fog_h", "sunshine_h", "precip_mm", "tmin_ground_c"):
                df[f"imgw{s}_{c}"] = g[c].reindex(day).to_numpy()
            df[f"imgw{s}_ground_frozen"] = (g["ground_state"].reindex(day) == "Z").to_numpy()
        df["imgw_dew_h"] = df[[c for c in df if c.endswith("_dew_h")]].mean(1)
        df["imgw_snow_depth_cm"] = df[[c for c in df if c.startswith("imgw") and c.endswith("snow_depth_cm")]].mean(1)
    # derived meteorology: station first, reanalysis as fallback, both kept
    df["st_dewpoint_c"] = core.dew_point_c(df.st_air_c, df.st_rh_pct)
    df["st_dpd_c"] = df.st_air_c - df.st_dewpoint_c
    df["om_dpd_c"] = df.om_era5_land_temperature_2m - df.om_era5_land_dew_point_2m
    air = df.st_air_c.fillna(df.om_era5_land_temperature_2m)
    dpd = df.st_dpd_c.fillna(df.om_dpd_c)
    wind = df.st_wind_ms.fillna(df.om_era5_wind_speed_10m / 3.6 * 0.75)    # 10 m km/h → ≈ 2 m m/s (log profile)
    rain = df.st_rain_mm.fillna(df.om_era5_precipitation)
    night = df.om_era5_shortwave_radiation.fillna(0) < 10
    df["air_c"], df["dpd_c"], df["wind_ms"], df["rain_mm"] = air, dpd, wind, rain
    df["rain_3h_mm"] = rain.rolling("3h", min_periods=1).sum()
    df["rain_24h_mm"] = rain.rolling("24h", min_periods=1).sum()
    df["rain_72h_mm"] = rain.rolling("72h", min_periods=1).sum()
    df["hours_since_rain"] = core.hours_since(rain > 0.1, idx)
    # flags (each with its source)
    df["flag_frozen_air"] = air <= 0
    df["flag_frozen_soil_om"] = df.om_era5_land_soil_temperature_0_to_7cm <= 0
    df["flag_snow_om"] = df.om_era5_land_snow_depth > 0.01
    if "imgw_snow_depth_cm" in df:
        df["flag_snow_imgw"] = df.imgw_snow_depth_cm > 0
    df["flag_wet_rh95"] = (df.st_rh_pct.fillna(df.om_era5_land_relative_humidity_2m) >= 95) | (df.rain_3h_mm > 0)
    df["flag_dew_likely"] = (dpd <= 1.0) & (wind <= 2.0) & night & (air > 0)
    df["flag_wet_any"] = df.flag_wet_rh95 | df.flag_dew_likely
    df.attrs["n_hours"] = len(df)
    log(f"site_hourly: {len(df)} hours × {df.shape[1]} columns")
    return df


# ------------------------------------------------------------------------------ folds

def folds() -> dict:
    return {"year": "calendar year of the date; a pair belongs to a year only if both dates do (core.pair_year_fold)",
            "v2_split": {"date": "2021-08-13", "rule": "calibrate before / test after, and the reverse (X-066)"},
            "plot_unit": "leave one radar unit out (plot_units.csv: P8 and P9 share a cell)",
            "spatial_block": f"static_40m.spatial_block — {SPATIAL_BLOCK_CELLS}×{SPATIAL_BLOCK_CELLS} cells "
                             f"({SPATIAL_BLOCK_CELLS * 40} m); 20 m grid uses the same blocks"}


# ------------------------------------------------------------------------------ dispatch

def run(step: str, what: list[str], ctx_fn, raw: Path, silver: Path, gold: Path, log=print) -> None:
    if step in ("harmonise", "all"):
        ctx = ctx_fn()
        todo = what or ["ifg", "rtc", "optical", "ecostress", "static", "site", "folds", "nisar"]
        for w in todo:
            if w == "ifg":
                harmonise_ifg(ctx, silver, log)
            elif w == "rtc":
                harmonise_rtc(ctx, silver, log)
            elif w == "optical":
                harmonise_optical(ctx, raw, silver, log)
            elif w == "static":
                harmonise_static(ctx, raw, silver, peat_bounds_utm(ctx, LIDAR_GRID_BUFFER_M, "AB"), log)
            elif w == "site":
                df = site_hourly(raw, f"{pd.Timestamp.now(tz='UTC'):%Y-%m-%d}", log)
                df.to_csv(silver / "site_hourly.csv.gz")
            elif w == "ecostress":
                harmonise_ecostress(ctx, raw, silver, log)
            elif w == "nisar":
                from .nisar_crops import to_silver
                if (raw / "nisar" / "granules.csv").exists():
                    tx, ty = tpl_xy(ctx)
                    to_silver(raw / "nisar", tx, ty, silver / "nisar_lband.nc", log)
            elif w == "folds":
                (silver / "folds.json").write_text(json.dumps(folds(), indent=2))
                plot_units().to_csv(silver / "plot_units.csv", index=False)
                log("folds.json, plot_units.csv")
    if step in ("gold", "check", "all"):
        from . import gold as G
    if step in ("gold", "all"):
        G.run(what if step == "gold" else [], ctx_fn, silver, gold, log)
    if step in ("check", "all"):
        G.check(silver, gold, raw, log)

"""Downloaders of the data cube's public sources (M1). Each writes into ``<cube>/raw/<source>/`` and is resumable: a file
already present is not fetched again.

| source | what | access |
|---|---|---|
| Sentinel-2 L2A (collection 1) | B03, B04, B05, B08, B8A, B11, SCL over the template, 2020 → | Earth Search STAC, COG windows |
| Landsat 8/9 C2 L2 | surface temperature (ST_B10) + QA over the template | Planetary Computer STAC, COG windows |
| ESA WorldCover 10 m (2020, 2021) | land-cover classes over the template | Planetary Computer STAC |
| ERA5-Land + ERA5 hourly | soil/skin/dew/snow (Land), rain/radiation/wind/cloud (ERA5) at the site | Open-Meteo archive API |
| IMGW synop daily | snow depth, dew and hoar-frost duration, ground state, sunshine (Piła 230, Poznań 330) | danepubliczne.imgw.pl |
| GUGiK LiDAR | DTM 1 m and DSM 0.5 m, EVRF2007 heights | Geoportal WCS 2.0.1 |

The network is the only side effect; parsing and arithmetic live in testable helpers (``parse_*``, ``s2_indices``).
"""
from __future__ import annotations

import io
import json
import time
import zipfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import requests

EARTH_SEARCH = "https://earth-search.aws.element84.com/v1"
PLANETARY = "https://planetarycomputer.microsoft.com/api/stac/v1"
OPEN_METEO = "https://archive-api.open-meteo.com/v1/archive"
IMGW_SYNOP = "https://danepubliczne.imgw.pl/data/dane_pomiarowo_obserwacyjne/dane_meteorologiczne/dobowe/synop"
GUGIK_DTM = ("https://mapy.geoportal.gov.pl/wss/service/PZGIK/NMT/GRID1/WCS/DigitalTerrainModel", "DTM_PL-EVRF2007-NH")
GUGIK_DSM = ("https://mapy.geoportal.gov.pl/wss/service/PZGIK/NMPT/GRID1/WCS/DigitalSurfaceModel", "DSM_PL-EVRF2007-NH")
IMGW_STATIONS = {230: "Piła", 330: "Poznań"}

GDAL_ENV = {"GDAL_DISABLE_READDIR_ON_OPEN": "EMPTY_DIR", "CPL_VSIL_CURL_ALLOWED_EXTENSIONS": ".tif,.TIF",
            "GDAL_HTTP_MAX_RETRY": "4", "GDAL_HTTP_RETRY_DELAY": "2", "VSI_CACHE": "TRUE"}

S2_BANDS = {"green": "B03", "red": "B04", "rededge1": "B05", "nir": "B08", "nir08": "B8A", "swir16": "B11", "scl": "SCL"}
# SCL classes kept as valid surface: vegetation 4, not vegetated 5, water 6, unclassified 7; snow (11) is kept apart
SCL_VALID = (4, 5, 6, 7)
SCL_SNOW = 11

OM_LAND = ["temperature_2m", "relative_humidity_2m", "dew_point_2m", "vapour_pressure_deficit", "snow_depth",
           "soil_temperature_0_to_7cm", "soil_temperature_7_to_28cm", "soil_moisture_0_to_7cm", "soil_moisture_7_to_28cm"]
OM_ERA5 = ["precipitation", "rain", "snowfall", "shortwave_radiation", "wind_speed_10m", "cloud_cover",
           "et0_fao_evapotranspiration", "surface_pressure"]

IMGW_COLUMNS = ["NSP", "POST", "ROK", "MC", "DZ", "TMAX", "WTMAX", "TMIN", "WTMIN", "STD", "WSTD", "TMNG", "WTMNG",
                "SMDB", "WSMDB", "ROOP", "PKSN", "WPKSN", "RWSN", "WRWSN", "USL", "WUSL", "DESZ", "WDESZ", "SNEG",
                "WSNEG", "DISN", "WDISN", "GRAD", "WGRAD", "MGLA", "WMGLA", "ZMGL", "WZMGL", "SADZ", "WSADZ", "GOLO",
                "WGOLO", "ZMNI", "WZMNI", "ZMWS", "WZMWS", "ZMET", "WZMET", "FF10", "WFF10", "FF15", "WFF15", "BRZA",
                "WBRZA", "ROSA", "WROSA", "SZRO", "WSZRO", "DZPS", "WDZPS", "DZBL", "WDZBL", "SGR", "IZD", "WIZD", "IZG",
                "WIZG", "AKTN", "WAKTN"]
IMGW_KEEP = {"TMAX": "tmax_c", "TMIN": "tmin_c", "STD": "tmean_c", "TMNG": "tmin_ground_c", "SMDB": "precip_mm",
             "PKSN": "snow_depth_cm", "USL": "sunshine_h", "ROSA": "dew_h", "SZRO": "hoarfrost_h", "MGLA": "fog_h",
             "SGR": "ground_state"}


def _get(url, params=None, retries=4, timeout=120, **kw) -> requests.Response:
    for k in range(retries):
        try:
            r = requests.get(url, params=params, timeout=timeout, **kw)
            if r.status_code == 200:
                return r
            if r.status_code in (429, 500, 502, 503, 504):
                time.sleep(5 * (k + 1))
                continue
            r.raise_for_status()
        except requests.RequestException:
            if k == retries - 1:
                raise
            time.sleep(5 * (k + 1))
    raise RuntimeError(f"{url}: no answer after {retries} tries")


# ------------------------------------------------------------------------------ grid helpers

def template_bounds(tpl_x, tpl_y) -> tuple[float, float, float, float]:
    """(left, bottom, right, top) edges of a centre-coordinate grid."""
    dx, dy = abs(float(tpl_x[1] - tpl_x[0])), abs(float(tpl_y[1] - tpl_y[0]))
    return (float(np.min(tpl_x)) - dx / 2, float(np.min(tpl_y)) - dy / 2,
            float(np.max(tpl_x)) + dx / 2, float(np.max(tpl_y)) + dy / 2)


def bbox_wgs84(bounds_utm, epsg: int = 32633) -> list[float]:
    from pyproj import Transformer
    t = Transformer.from_crs(epsg, 4326, always_xy=True)
    xs, ys = t.transform([bounds_utm[0], bounds_utm[2], bounds_utm[0], bounds_utm[2]],
                         [bounds_utm[1], bounds_utm[1], bounds_utm[3], bounds_utm[3]])
    return [min(xs), min(ys), max(xs), max(ys)]


def read_window(href: str, bounds_utm, res: float, epsg: int = 32633, resampling: str = "nearest",
                env: dict | None = None) -> np.ndarray:
    """A raster (any CRS) read over ``bounds_utm`` on an aligned ``res`` grid in ``epsg`` (a WarpedVRT window: only the
    needed COG blocks are fetched). Nodata → NaN."""
    import rasterio
    from rasterio.enums import Resampling
    from rasterio.transform import from_origin
    from rasterio.vrt import WarpedVRT

    left, bottom, right, top = bounds_utm
    w, h = int(round((right - left) / res)), int(round((top - bottom) / res))
    with rasterio.Env(**(env or GDAL_ENV)):
        with rasterio.open(href) as src:
            with WarpedVRT(src, crs=f"EPSG:{epsg}", transform=from_origin(left, top, res, res), width=w, height=h,
                           resampling=getattr(Resampling, resampling), nodata=src.nodata) as vrt:
                a = vrt.read(1).astype("float32")
                if vrt.nodata is not None:
                    a[a == vrt.nodata] = np.nan
    return a


# ------------------------------------------------------------------------------ Sentinel-2

def s2_indices(b: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    """Optical indices on the 20 m grid from surface reflectances (scale and offset applied; 10 m bands already
    block-averaged to 20 m). Invalid SCL → NaN; snow kept as its own fraction."""
    scl = b["SCL"]
    valid = np.isin(scl, SCL_VALID)
    snow = (scl == SCL_SNOW).astype("float32")

    def nd(a, c):
        with np.errstate(invalid="ignore", divide="ignore"):
            v = (a - c) / (a + c)
        return np.where(valid & np.isfinite(v), v, np.nan).astype("float32")

    return {"ndvi": nd(b["B08"], b["B04"]), "ndre": nd(b["B8A"], b["B05"]), "ndmi": nd(b["B8A"], b["B11"]),
            "ndwi": nd(b["B03"], b["B08"]), "mndwi": nd(b["B03"], b["B11"]),
            "valid": valid.astype("float32"), "snow": np.where(np.isfinite(scl), snow, np.nan)}


def block_mean(a: np.ndarray, k: int) -> np.ndarray:
    h, w = a.shape[0] // k * k, a.shape[1] // k * k
    with np.errstate(invalid="ignore"):
        return np.nanmean(a[:h, :w].reshape(h // k, k, w // k, k), axis=(1, 3))


def s2_items(bbox, start: str, end: str, max_cloud: float = 90.0) -> list:
    """Sentinel-2 L2A items over the box: collection 1 (reprocessed, consistent baseline) first; the older
    ``sentinel-2-l2a`` fills the acquisitions collection 1 lacks (Earth Search's collection 1 holds only a handful of
    2022 scenes). One item per acquisition instant (to the minute) and tile."""
    import warnings

    from pystac_client import Client
    warnings.filterwarnings("ignore", message="Could not parse bucket")
    c = Client.open(EARTH_SEARCH)
    seen, out = set(), []
    for coll in ("sentinel-2-c1-l2a", "sentinel-2-l2a"):
        for it in c.search(collections=[coll], bbox=bbox, datetime=f"{start}/{end}",
                           query={"eo:cloud_cover": {"lt": max_cloud}}).items():
            key = (f"{it.datetime:%Y%m%dT%H%M}", it.properties.get("s2:mgrs_tile") or it.properties.get("grid:code"))
            if key in seen:
                continue
            seen.add(key)
            it.extra_fields["cube_collection"] = coll
            out.append(it)
    return sorted(out, key=lambda i: i.datetime)


def asset_scale_offset(asset) -> tuple[float, float]:
    """Reflectance = DN · scale + offset, from the asset's ``raster:bands`` (collection 1 and baseline ≥ 04.00 carry
    offset −0.1; older scenes 0)."""
    rb = (asset.extra_fields.get("raster:bands") or [{}])[0]
    return float(rb.get("scale", 1e-4)), float(rb.get("offset", 0.0))


def fetch_s2(out_dir: Path, bounds_utm, start: str, end: str, max_cloud: float = 90.0, workers: int = 6,
             log=print) -> pd.DataFrame:
    """Every Sentinel-2 scene over the template on the 20 m grid aligned with the template's edges: surface
    reflectance of B03, B04, B05, B08, B8A, B11 (int16, ×10⁴, scale and offset applied; 10 m bands block-averaged),
    SCL, and the indices; one ``.npz`` per scene (resumable). Returns the scene table."""
    out_dir.mkdir(parents=True, exist_ok=True)
    items = s2_items(bbox_wgs84(bounds_utm), start, end, max_cloud)
    log(f"S2: {len(items)} items {start}..{end}")

    def one(it):
        f = out_dir / f"{it.id}.npz"
        base = {"id": it.id, "time": it.datetime, "cloud": it.properties.get("eo:cloud_cover"),
                "collection": it.extra_fields.get("cube_collection"),
                "baseline": it.properties.get("s2:processing_baseline"), "file": f.name}
        if f.exists():
            return {**base, "valid_share": float(np.nanmean(np.load(f)["valid"]))}
        try:
            b = {}
            for key, name in S2_BANDS.items():
                res = 10 if name in ("B03", "B04", "B08") else 20
                a = read_window(it.assets[key].href, bounds_utm, res)
                if name != "SCL":
                    sc, off = asset_scale_offset(it.assets[key])
                    a[~(a > 0)] = np.nan                       # 0 = nodata
                    a = a * sc + off
                b[name] = block_mean(a, 2) if res == 10 else a
            idx = s2_indices(b)
            bands = {f"refl_{n}": np.where(np.isfinite(v), np.round(v * 1e4), -32768).astype("int16")
                     for n, v in b.items() if n != "SCL"}
            np.savez_compressed(f, **idx, **bands, scl=np.nan_to_num(b["SCL"], nan=0).astype("uint8"))
            return {**base, "valid_share": float(np.nanmean(idx["valid"]))}
        except Exception as e:  # noqa: BLE001 — one bad scene must not stop the archive
            return {**base, "valid_share": np.nan, "file": "", "error": str(e)[:200]}

    rows = []
    with ThreadPoolExecutor(workers) as ex:
        for k, r in enumerate(ex.map(one, items)):
            rows.append(r)
            if k % 50 == 0:
                log(f"  S2 {k}/{len(items)}")
    tab = pd.DataFrame(rows).sort_values("time")
    tab.to_csv(out_dir / "scenes.csv", index=False)
    return tab


# ------------------------------------------------------------------------------ Landsat surface temperature

def landsat_qa_clear(qa: np.ndarray) -> np.ndarray:
    """Clear-surface mask from Landsat C2 QA_PIXEL: not fill (bit 0), dilated cloud (1), cirrus (2), cloud (3),
    cloud shadow (4) or snow (5)."""
    q = np.nan_to_num(qa, nan=1).astype(np.uint16)
    bad = 0
    for bit in (0, 1, 2, 3, 4, 5):
        bad |= 1 << bit
    return (q & bad) == 0


def fetch_landsat_lst(out_dir: Path, bounds_utm, start: str, end: str, max_cloud: float = 80.0, workers: int = 4,
                      log=print) -> pd.DataFrame:
    """Landsat 8/9 surface temperature (°C) on a 30 m grid over the template, clear pixels only; one ``.npz`` per scene."""
    import planetary_computer
    from pystac_client import Client
    out_dir.mkdir(parents=True, exist_ok=True)
    c = Client.open(PLANETARY, modifier=planetary_computer.sign_inplace)
    items = list(c.search(collections=["landsat-c2-l2"], bbox=bbox_wgs84(bounds_utm), datetime=f"{start}/{end}",
                          query={"eo:cloud_cover": {"lt": max_cloud}, "platform": {"in": ["landsat-8", "landsat-9"]}}).items())
    log(f"Landsat: {len(items)} items")

    def one(it):
        f = out_dir / f"{it.id}.npz"
        base = {"id": it.id, "time": it.datetime, "cloud": it.properties.get("eo:cloud_cover"), "file": f.name}
        if f.exists():
            return {**base, "clear_share": float(np.isfinite(np.load(f)["lst_c"]).mean())}
        try:
            st = read_window(it.assets["lwir11"].href, bounds_utm, 30, resampling="nearest")
            qa = read_window(it.assets["qa_pixel"].href, bounds_utm, 30, resampling="nearest")
            lst = st * 0.00341802 + 149.0 - 273.15
            lst[~landsat_qa_clear(qa) | ~np.isfinite(st) | (st == 0)] = np.nan
            np.savez_compressed(f, lst_c=lst.astype("float32"))
            return {**base, "clear_share": float(np.isfinite(lst).mean())}
        except Exception as e:  # noqa: BLE001
            return {**base, "file": "", "clear_share": np.nan, "error": str(e)[:200]}

    with ThreadPoolExecutor(workers) as ex:
        rows = list(ex.map(one, items))
    tab = pd.DataFrame(rows).sort_values("time")
    tab.to_csv(out_dir / "scenes.csv", index=False)
    return tab


# ------------------------------------------------------------------------------ WorldCover

def fetch_worldcover(out_dir: Path, bounds_utm, log=print) -> list[Path]:
    """ESA WorldCover 10 m (2020 v100, 2021 v200) over the template, nearest-neighbour on the 10 m UTM grid."""
    import planetary_computer
    from pystac_client import Client
    out_dir.mkdir(parents=True, exist_ok=True)
    c = Client.open(PLANETARY, modifier=planetary_computer.sign_inplace)
    out = []
    for it in c.search(collections=["esa-worldcover"], bbox=bbox_wgs84(bounds_utm)).items():
        year = str(it.properties.get("esa_worldcover:product_version", it.id))
        f = out_dir / f"worldcover_{it.datetime.year if it.datetime else year}_{it.id}.npz"
        if not f.exists():
            a = read_window(it.assets["map"].href, bounds_utm, 10)
            np.savez_compressed(f, classes=a)
        out.append(f)
        log(f"WorldCover {f.name}")
    return out


# ------------------------------------------------------------------------------ Open-Meteo (ERA5-Land, ERA5)

def fetch_open_meteo(out_dir: Path, lat: float, lon: float, start: str, end: str, log=print) -> pd.DataFrame:
    """Hourly ERA5-Land (land surface, dew, snow, soil) and ERA5 (rain, radiation, wind, cloud) at the site, UTC, one
    request per year and model (cached as JSON)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    frames = []
    for model, var in (("era5_land", OM_LAND), ("era5", OM_ERA5)):
        parts = []
        for y in range(pd.Timestamp(start).year, pd.Timestamp(end).year + 1):
            s, e = max(pd.Timestamp(start), pd.Timestamp(f"{y}-01-01")), min(pd.Timestamp(end), pd.Timestamp(f"{y}-12-31"))
            f = out_dir / f"{model}_{y}.json"
            if not f.exists() or y >= pd.Timestamp.now().year:
                r = _get(OPEN_METEO, params=dict(latitude=lat, longitude=lon, start_date=f"{s:%Y-%m-%d}",
                                                 end_date=f"{e:%Y-%m-%d}", models=model, hourly=",".join(var),
                                                 timezone="UTC"))
                f.write_text(r.text)
                log(f"Open-Meteo {model} {y}")
            parts.append(parse_open_meteo(json.loads(f.read_text()), prefix=f"{model}_"))
        frames.append(pd.concat(parts).sort_index())
    df = pd.concat(frames, axis=1)
    df = df[~df.index.duplicated()]
    return df.dropna(how="all")


def parse_open_meteo(j: dict, prefix: str = "") -> pd.DataFrame:
    h = j["hourly"]
    idx = pd.DatetimeIndex(pd.to_datetime(h.pop("time")), name="time_utc").tz_localize("UTC")
    return pd.DataFrame({f"{prefix}{k}": pd.to_numeric(pd.Series(v), errors="coerce").to_numpy() for k, v in h.items()},
                        index=idx)


# ------------------------------------------------------------------------------ IMGW daily synop

def parse_imgw_daily(raw: bytes, stations=IMGW_STATIONS) -> pd.DataFrame:
    """IMGW ``s_d_*.csv`` (cp1250, no header) → one row per station-day with the kept variables; status 8 (not
    measured) → NaN, status 9 (phenomenon absent) → 0."""
    df = pd.read_csv(io.BytesIO(raw), header=None, names=IMGW_COLUMNS, encoding="cp1250", dtype=str)
    df = df[(pd.to_numeric(df.NSP, errors="coerce") % 1000).isin(list(stations))]     # codes like 353160230
    out = pd.DataFrame({"station": pd.to_numeric(df.NSP) % 1000,
                        "date": pd.to_datetime(dict(year=pd.to_numeric(df.ROK), month=pd.to_numeric(df.MC),
                                                    day=pd.to_numeric(df.DZ)))})
    for src, dst in IMGW_KEEP.items():
        if src == "SGR":
            out[dst] = df[src].str.strip().replace("", np.nan).to_numpy()
            continue
        v = pd.to_numeric(df[src].str.replace(",", "."), errors="coerce")
        st = df.get(f"W{src}")
        if st is not None:
            st = st.str.strip()
            v = v.mask(st == "8").mask(st == "9", 0.0)
        out[dst] = v.to_numpy()
    return out.reset_index(drop=True)


def fetch_imgw(out_dir: Path, start_year: int, end_year: int, stations=IMGW_STATIONS, log=print) -> pd.DataFrame:
    """Daily synop for the stations: per-station yearly zips up to last year, monthly all-station zips this year."""
    out_dir.mkdir(parents=True, exist_ok=True)
    frames = []
    this_year = pd.Timestamp.now().year
    for y in range(start_year, end_year + 1):
        names = ([f"{y}/{y}_{s}_s.zip" for s in stations] if y < this_year
                 else [f"{y}/{y}_{m:02d}_s.zip" for m in range(1, 13)])
        for n in names:
            f = out_dir / n.replace("/", "_")
            if not f.exists() or y == this_year:
                try:
                    r = _get(f"{IMGW_SYNOP}/{n}", retries=2)
                except Exception:  # noqa: BLE001 — a month not yet published
                    continue
                f.write_bytes(r.content)
                log(f"IMGW {n}")
            with zipfile.ZipFile(f) as z:
                for member in z.namelist():
                    if member.startswith("s_d_") and not member.startswith("s_d_t"):
                        frames.append(parse_imgw_daily(z.read(member), stations))
    df = pd.concat(frames).drop_duplicates(["station", "date"]).sort_values(["station", "date"])
    return df.reset_index(drop=True)


# ------------------------------------------------------------------------------ GUGiK LiDAR (WCS)

def fetch_gugik_tiles(out_dir: Path, bounds_2180, tile_m: float = 500.0, log=print) -> pd.DataFrame:
    """DTM (1 m) and DSM (0.5 m) EVRF2007 GeoTIFF tiles from the GUGiK WCS over a PL-1992 (EPSG:2180) box."""
    out_dir.mkdir(parents=True, exist_ok=True)
    left, bottom, right, top = bounds_2180
    rows = []
    for x0 in np.arange(left, right, tile_m):
        for y0 in np.arange(bottom, top, tile_m):
            for kind, (base, cov) in (("dtm", GUGIK_DTM), ("dsm", GUGIK_DSM)):
                f = out_dir / f"{kind}_{int(x0)}_{int(y0)}.tif"
                if not f.exists():
                    url = (f"{base}?SERVICE=WCS&VERSION=2.0.1&REQUEST=GetCoverage&COVERAGEID={cov}"
                           f"&SUBSET=x({x0},{x0 + tile_m})&SUBSET=y({y0},{y0 + tile_m})&FORMAT=image/tiff")
                    r = _get(url, timeout=300)
                    if r.content[:4] not in (b"II*\x00", b"MM\x00*"):
                        raise RuntimeError(f"GUGiK {kind} tile {x0},{y0}: not a TIFF")
                    f.write_bytes(r.content)
                rows.append({"kind": kind, "x0": x0, "y0": y0, "file": f.name})
        log(f"LiDAR column x0={x0:.0f} done")
    return pd.DataFrame(rows)

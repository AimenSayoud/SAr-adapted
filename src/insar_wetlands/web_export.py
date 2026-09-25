"""Turn the pipeline's rasters and tables into files a web browser can draw.

Why this module exists. The results of this project live in fourteen NetCDF
stacks on the Drive, a few hundred cropped interferograms, thirty CSVs and a
handful of GeoJSONs. A reader who wants to *see* the negative result — that the
per-pixel maps are noise while the zone aggregate carries a seasonal signal —
has to open QGIS and a notebook side by side. The web atlas puts them on one
map, and this module is the only bridge between the analysis grid and it.

Three rules, each a lesson already paid for elsewhere in this repo:

* **No resampling.** Every raster ships on its native 40 m UTM grid, with the
  four corners of the grid given in lon/lat. The browser places the image by
  those corners. Reprojecting to Web Mercator would interpolate — or, with
  nearest-neighbour, duplicate and drop — pixels, and the per-pixel noise is
  precisely what a viewer must be able to inspect.
* **Values, not pictures.** Rasters ship as quantised numbers plus a
  scale/offset, so the browser can colour them, read a value under the cursor
  and draw a pixel's time series. A PNG would freeze one colour ramp and lose
  the numbers.
* **Provenance travels with the data.** Every layer records the file it came
  from and that file's sha256, so a number seen on the map can be traced back
  to the product — the same principle as `paper_numbers`, applied to pixels.

The per-pixel seasonal fit here is the vectorised form of
`aggregate.seasonal_amplitude` and must agree with it pixel for pixel
(`tests/test_web_export.py` checks this). It exists to *show* the per-pixel
failure on a map, not to rescue it: H1 is settled.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

GRID_CRS = "EPSG:32633"
INT16_NODATA = -32768
UINT8_NODATA = 255

SEASONS = {12: "DJF", 1: "DJF", 2: "DJF", 3: "MAM", 4: "MAM", 5: "MAM",
           6: "JJA", 7: "JJA", 8: "JJA", 9: "SON", 10: "SON", 11: "SON"}


# --- geometry -----------------------------------------------------------------

def pixel_step(coord: np.ndarray) -> float:
    """Signed spacing of a regular 1-D coordinate (negative for north-up y)."""
    c = np.asarray(coord, dtype=float)
    if c.size < 2:
        raise ValueError("need at least two coordinates to infer a pixel step")
    steps = np.diff(c)
    if not np.allclose(steps, steps[0], rtol=1e-6, atol=1e-6):
        raise ValueError("coordinate is not regular")
    return float(steps[0])


def grid_edges(x: np.ndarray, y: np.ndarray) -> tuple[float, float, float, float]:
    """Outer pixel edges (west, south, east, north) from pixel-centre coords."""
    dx, dy = pixel_step(x), pixel_step(y)
    xs = [float(x[0]) - dx / 2, float(x[-1]) + dx / 2]
    ys = [float(y[0]) - dy / 2, float(y[-1]) + dy / 2]
    return min(xs), min(ys), max(xs), max(ys)


def grid_corners_lonlat(x: np.ndarray, y: np.ndarray,
                        crs: str = GRID_CRS) -> list[list[float]]:
    """Grid corners in lon/lat, ordered top-left, top-right, bottom-right,
    bottom-left — the order MapLibre image sources and deck.gl bitmaps expect.

    The corners are of the outer pixel *edges*, not of the first/last pixel
    centres: placing an image by its centres shifts it by half a pixel (20 m
    here), which is exactly the offset found between MintPy's corner-referenced
    coordinates and the rest of the stack."""
    from pyproj import Transformer

    w, s, e, n = grid_edges(x, y)
    tr = Transformer.from_crs(crs, "EPSG:4326", always_xy=True)
    pts = [(w, n), (e, n), (e, s), (w, s)]
    return [[round(v, 7) for v in tr.transform(px, py)] for px, py in pts]


# --- numbers ------------------------------------------------------------------

def robust_range(a: np.ndarray, lo: float = 2.0, hi: float = 98.0,
                 symmetric: bool = False) -> tuple[float, float]:
    """Percentile display range over finite values; symmetric about 0 if asked."""
    v = np.asarray(a, dtype=float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return (0.0, 1.0)
    a0, a1 = np.percentile(v, [lo, hi])
    if symmetric:
        m = float(max(abs(a0), abs(a1)))
        return (-m, m) if m > 0 else (-1.0, 1.0)
    if a1 <= a0:
        a1 = a0 + 1.0
    return float(a0), float(a1)


def quantize(a: np.ndarray, vmin: float | None = None, vmax: float | None = None,
             dtype: str = "int16") -> tuple[np.ndarray, float, float]:
    """Linear quantisation with NaN → nodata. Returns (q, scale, offset) such
    that value = q * scale + offset. Values outside [vmin, vmax] are clipped, so
    pass the true min/max (the default) unless clipping is intended."""
    a = np.asarray(a, dtype="float64")
    finite = np.isfinite(a)
    if vmin is None:
        vmin = float(np.min(a[finite])) if finite.any() else 0.0
    if vmax is None:
        vmax = float(np.max(a[finite])) if finite.any() else 1.0
    if vmax <= vmin:
        vmax = vmin + 1.0
    if dtype == "int16":
        qmin, qmax, nod = -32767, 32767, INT16_NODATA
    elif dtype == "uint8":
        qmin, qmax, nod = 0, 254, UINT8_NODATA
    else:
        raise ValueError(dtype)
    scale = (vmax - vmin) / (qmax - qmin)
    offset = vmin - qmin * scale
    q = np.full(a.shape, nod, dtype=dtype)
    q[finite] = np.clip(np.round((a[finite] - offset) / scale), qmin, qmax).astype(dtype)
    return q, float(scale), float(offset)


def dequantize(q: np.ndarray, scale: float, offset: float) -> np.ndarray:
    nod = INT16_NODATA if q.dtype == np.int16 else UINT8_NODATA
    out = q.astype("float64") * scale + offset
    out[q == nod] = np.nan
    return out


# --- derived fields -------------------------------------------------------------

def seasonal_fit_stack(stack: np.ndarray, dates, epoch=None,
                       min_obs: int = 6) -> dict[str, np.ndarray]:
    """Per-pixel annual fit, the vectorised twin of `aggregate.seasonal_amplitude`.

    Fits y = c + d*t + a*cos(2πt) + b*sin(2πt) (t in years from `epoch`,
    default the first date) independently at every pixel, using only that
    pixel's finite samples, and returns maps of amplitude, day-of-year of the
    maximum, trend and seasonal R² — the same definitions, so a pixel's value
    here equals `seasonal_amplitude` run on that pixel's series.

    `stack` is (time, y, x). Pixels with fewer than `min_obs` finite samples
    are NaN."""
    Y = np.asarray(stack, dtype="float64")
    T = Y.shape[0]
    d = pd.to_datetime(pd.Index(dates))
    if len(d) != T:
        raise ValueError(f"{len(d)} dates for {T} layers")
    t0 = pd.Timestamp(epoch) if epoch is not None else d[0]
    t = (d - t0).days.values / 365.25
    M = np.column_stack([np.ones(T), t, np.cos(2 * np.pi * t), np.sin(2 * np.pi * t)])

    flat = Y.reshape(T, -1)
    W = np.isfinite(flat)
    y0 = np.where(W, flat, 0.0)
    n = W.sum(0)
    good = n >= min_obs

    out = {k: np.full(flat.shape[1], np.nan) for k in
           ("amplitude_mm", "phase_doy", "trend_mm_yr", "r2_seasonal", "n")}
    out["n"] = n.astype("float64")
    if good.any():
        Wg, yg = W[:, good].astype("float64"), y0[:, good]
        # Normal equations per pixel: (Mᵀ diag(w) M) β = Mᵀ diag(w) y
        N4 = np.einsum("tp,ti,tj->pij", Wg, M, M)
        r4 = np.einsum("tp,ti,tp->pi", Wg, M, yg)
        ok4 = np.abs(np.linalg.det(N4)) > 1e-10
        beta = np.full((Wg.shape[1], 4), np.nan)
        beta[ok4] = np.linalg.solve(N4[ok4], r4[ok4][..., None])[..., 0]
        resid = (yg - M @ beta.T) * Wg
        ss = (resid ** 2).sum(0)

        M0 = M[:, :2]
        N2 = np.einsum("tp,ti,tj->pij", Wg, M0, M0)
        r2 = np.einsum("tp,ti,tp->pi", Wg, M0, yg)
        ok2 = np.abs(np.linalg.det(N2)) > 1e-12
        b0 = np.full((Wg.shape[1], 2), np.nan)
        b0[ok2] = np.linalg.solve(N2[ok2], r2[ok2][..., None])[..., 0]
        ss0 = (((yg - M0 @ b0.T) * Wg) ** 2).sum(0)

        amp = np.hypot(beta[:, 2], beta[:, 3])
        phase = (np.arctan2(beta[:, 3], beta[:, 2]) / (2 * np.pi)) % 1.0 * 365.25
        doy0 = int(t0.dayofyear)
        with np.errstate(invalid="ignore", divide="ignore"):
            r2s = np.where(ss0 > 0, 1.0 - ss / ss0, np.nan)
        out["amplitude_mm"][good] = amp
        out["phase_doy"][good] = (doy0 + phase) % 365.25
        out["trend_mm_yr"][good] = beta[:, 1]
        out["r2_seasonal"][good] = r2s
    shape = Y.shape[1:]
    return {k: v.reshape(shape) for k, v in out.items()}


def amplitude_dispersion_from_db(gamma0_db: np.ndarray) -> np.ndarray:
    """Amplitude dispersion D_A = σ_A / μ_A over time, from a γ⁰ stack in dB.

    Amplitude is √(linear power). D_A < 0.25 is the usual persistent-scatterer
    candidate threshold; a floating mat is expected to sit far above it. This is
    computed from the RTC backscatter stack because the cropped HyP3 INSAR
    products carry no amplitude layer."""
    p = 10.0 ** (np.asarray(gamma0_db, dtype="float64") / 10.0)
    amp = np.sqrt(p)
    with np.errstate(invalid="ignore", divide="ignore"):
        mu = np.nanmean(amp, axis=0)
        sd = np.nanstd(amp, axis=0)
        n = np.isfinite(amp).sum(0)
        da = np.where(n >= 2, sd / mu, np.nan)
    return da


def pair_midpoint_season(pairs: list[str]) -> list[str]:
    """Season (DJF/MAM/JJA/SON) of each pair's temporal midpoint."""
    out = []
    for p in pairs:
        a, b = (pd.Timestamp(s) for s in p.split("_"))
        out.append(SEASONS[(a + (b - a) / 2).month])
    return out


def seasonal_mean_maps(stack: np.ndarray, pairs: list[str],
                       max_dt_days: int | None = 48) -> dict[str, np.ndarray]:
    """Mean of a per-pair stack (pair, y, x) by the season of each pair's
    midpoint, keeping only pairs with Δt ≤ `max_dt_days` so that seasons are
    compared at the same baselines (coherence falls with Δt; mixing a season's
    long bridge pairs with another's short ones would confound the two)."""
    S = np.asarray(stack, dtype="float64")
    seasons = np.array(pair_midpoint_season(pairs))
    dts = np.array([(pd.Timestamp(p.split("_")[1]) - pd.Timestamp(p.split("_")[0])).days
                    for p in pairs])
    keep = np.ones(len(pairs), bool) if max_dt_days is None else dts <= max_dt_days
    out = {}
    for s in ("DJF", "MAM", "JJA", "SON"):
        m = keep & (seasons == s)
        if m.any():
            with np.errstate(invalid="ignore"):
                out[s] = np.nanmean(S[m], axis=0)
            out[f"{s}_n_pairs"] = int(m.sum())
    return out


def labels_to_geojson(labels: np.ndarray, x: np.ndarray, y: np.ndarray,
                      names: dict[int, dict], crs: str = GRID_CRS) -> dict:
    """Polygonise an integer label raster into a lon/lat GeoJSON FeatureCollection.
    One feature per label (a MultiPolygon), carrying `names[label]` as properties
    and its pixel count, so the outline and the count come from the same array."""
    from pyproj import Transformer
    from rasterio.features import shapes
    from rasterio.transform import from_origin
    from shapely.geometry import mapping, shape
    from shapely.ops import transform as shp_transform
    from shapely.ops import unary_union

    dx, dy = pixel_step(x), pixel_step(y)
    w, s, e, n = grid_edges(x, y)
    tf = from_origin(w, n, abs(dx), abs(dy))
    lab = np.asarray(labels).astype("int32")
    if dy > 0:           # south-up grid: flip so row 0 is north
        lab = lab[::-1]
    tr = Transformer.from_crs(crs, "EPSG:4326", always_xy=True)
    feats = []
    for val, props in names.items():
        polys = [shape(g) for g, v in shapes(lab, mask=(lab == val), transform=tf) if v == val]
        if not polys:
            continue
        geom = shp_transform(tr.transform, unary_union(polys))
        feats.append({"type": "Feature",
                      "properties": {**props, "code": int(val),
                                     "n_px": int((lab == val).sum())},
                      "geometry": mapping(geom)})
    return {"type": "FeatureCollection", "features": feats}


# --- files ----------------------------------------------------------------------

def sha256_file(path: str | Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def _jsonable(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, float) and not np.isfinite(o):
        return None
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (pd.Timestamp,)):
        return o.strftime("%Y-%m-%d")
    if isinstance(o, Path):
        return str(o)
    raise TypeError(type(o))


def _clean(o):
    """Recursively replace non-finite floats by None (JSON has no NaN)."""
    if isinstance(o, np.ndarray):
        return _clean(o.tolist())
    if isinstance(o, pd.Series):
        return _clean(o.tolist())
    if isinstance(o, dict):
        return {k: _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if isinstance(o, (float, np.floating)):
        return float(o) if np.isfinite(o) else None
    if isinstance(o, np.integer):
        return int(o)
    return o


def write_json(path: str | Path, obj, indent: int | None = None) -> Path:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(_clean(obj), default=_jsonable, indent=indent,
                            ensure_ascii=False, allow_nan=False))
    return p


@dataclass
class Provenance:
    """Where a layer came from. `source` paths are recorded relative to the data
    root they were read from, with the file's sha256 at export time."""
    sources: list[dict] = field(default_factory=list)

    def add(self, path: str | Path, root: str | Path | None = None, **extra) -> Provenance:
        p = Path(path)
        rel = str(p.relative_to(root)) if root and p.is_relative_to(root) else p.name
        self.sources.append({"file": rel, "sha256": sha256_file(p), **extra})
        return self


class AtlasWriter:
    """Accumulates layers under `out_dir` and writes one `manifest.json`.

    Layout:
      rasters/<id>.bin   quantised values, row-major (y, x) or (t, y, x)
      vectors/<id>.geojson
      charts/<id>.json
      manifest.json      every layer: kind, dtype, shape, scale/offset,
                         corners, units, display range, status, provenance
    """

    # What a viewer may conclude from a layer. Shown as a badge on every layer.
    STATUSES = {
        "core",                 # a main result of the project (backs a key figure or table)
        "supporting",           # robustness / supplementary analysis of a main result
        "pipeline",             # output of a current phase, not a headline result
        "failed-estimator",     # per-pixel inversion output, shown as H1 evidence
        "exploratory",          # new analysis, not yet tested or registered
        "derived",              # computed for the atlas by this module
        "context",              # terrain, land cover, outlines
        "unverified",           # provenance not established
        "superseded",           # an older vintage kept for comparison only
    }

    def __init__(self, out_dir: str | Path, x: np.ndarray, y: np.ndarray,
                 crs: str = GRID_CRS):
        self.out = Path(out_dir)
        self.x, self.y, self.crs = np.asarray(x), np.asarray(y), crs
        self.corners = grid_corners_lonlat(self.x, self.y, crs)
        self.layers: list[dict] = []
        self.meta: dict = {}

    def _check_grid(self, x, y):
        if x is None:
            return
        if len(x) != len(self.x) or len(y) != len(self.y) or \
                not (np.allclose(x, self.x) and np.allclose(y, self.y)):
            raise ValueError("layer is not on the atlas grid — align it first "
                             "(a silent half-pixel shift is how the MintPy grid "
                             "would otherwise be drawn 20 m off)")

    def _base(self, lid, title, group, status, units, description, prov):
        if status not in self.STATUSES:
            raise ValueError(f"unknown status {status!r}")
        if any(lyr["id"] == lid for lyr in self.layers):
            raise ValueError(f"duplicate layer id {lid!r}")
        return {"id": lid, "title": title, "group": group, "status": status,
                "units": units, "description": description,
                "provenance": (prov.sources if prov else [])}

    def raster(self, lid: str, arr: np.ndarray, *, title: str, group: str,
               status: str, units: str = "", description: str = "",
               colormap: str = "viridis", display: tuple | None = None,
               symmetric: bool = False, categories: dict | None = None,
               prov: Provenance | None = None, x=None, y=None,
               dtype: str = "int16", times: list | None = None,
               time_label: str = "date", extra: dict | None = None) -> dict:
        """A 2-D field (y, x) or a stack (t, y, x) on the atlas grid."""
        self._check_grid(x, y)
        a = np.asarray(arr, dtype="float64")
        if a.shape[-2:] != (len(self.y), len(self.x)):
            raise ValueError(f"{lid}: shape {a.shape} vs grid {(len(self.y), len(self.x))}")
        if a.ndim == 3 and (times is None or len(times) != a.shape[0]):
            raise ValueError(f"{lid}: a stack needs one time label per layer")
        finite = np.isfinite(a)
        if categories:
            dtype = "uint8"
            q = np.full(a.shape, UINT8_NODATA, "uint8")
            q[finite] = a[finite].astype("uint8")
            scale, offset = 1.0, 0.0
        else:
            lo, hi = (float(np.nanmin(a)), float(np.nanmax(a))) if finite.any() else (0.0, 1.0)
            q, scale, offset = quantize(a, lo, hi, dtype)
        path = self.out / "rasters" / f"{lid}.bin"
        path.parent.mkdir(parents=True, exist_ok=True)
        q.astype("<" + q.dtype.str[1:]).tofile(path)
        rng = display or robust_range(a, symmetric=symmetric)
        entry = self._base(lid, title, group, status, units, description, prov)
        entry.update({
            "kind": "stack" if a.ndim == 3 else "raster",
            "file": f"rasters/{lid}.bin", "dtype": q.dtype.name,
            "shape": list(a.shape), "scale": scale, "offset": offset,
            "nodata": UINT8_NODATA if q.dtype == np.uint8 else INT16_NODATA,
            "colormap": colormap, "display": [float(rng[0]), float(rng[1])],
            "stats": {"min": float(np.nanmin(a)) if finite.any() else None,
                      "max": float(np.nanmax(a)) if finite.any() else None,
                      "finite_fraction": float(finite.mean())},
            "bytes": path.stat().st_size,
        })
        if categories:
            entry["categories"] = {str(k): v for k, v in categories.items()}
        if a.ndim == 3:
            entry["times"] = [str(t) for t in times]
            entry["time_label"] = time_label
        if extra:
            entry.update(extra)
        self.layers.append(entry)
        return entry

    def vector(self, lid: str, geojson: dict, *, title: str, group: str,
               status: str, description: str = "", style: dict | None = None,
               prov: Provenance | None = None) -> dict:
        path = write_json(self.out / "vectors" / f"{lid}.geojson", geojson)
        entry = self._base(lid, title, group, status, "", description, prov)
        entry.update({"kind": "vector", "file": f"vectors/{lid}.geojson",
                      "n_features": len(geojson.get("features", [])),
                      "style": style or {}, "bytes": path.stat().st_size})
        self.layers.append(entry)
        return entry

    def chart(self, lid: str, data, *, title: str, group: str, status: str,
              description: str = "", units: str = "",
              prov: Provenance | None = None, extra: dict | None = None) -> dict:
        if isinstance(data, pd.DataFrame):
            data = {"columns": list(map(str, data.columns)),
                    "rows": data.astype(object).where(pd.notna(data), None).values.tolist()}
        path = write_json(self.out / "charts" / f"{lid}.json", data)
        entry = self._base(lid, title, group, status, units, description, prov)
        entry.update({"kind": "chart", "file": f"charts/{lid}.json",
                      "bytes": path.stat().st_size})
        if extra:
            entry.update(extra)
        self.layers.append(entry)
        return entry

    def image(self, lid: str, rgb: np.ndarray, *, title: str, group: str,
              status: str, description: str = "",
              edges: tuple[float, float, float, float] | None = None,
              prov: Provenance | None = None, extra: dict | None = None) -> dict:
        """An RGB(A) picture (uint8, rows north-first) placed by the corners of
        `edges` (UTM west, south, east, north; default the atlas grid). Used for
        the optical basemap, which has its own 10 m resolution."""
        from PIL import Image
        a = np.asarray(rgb)
        if a.dtype != np.uint8 or a.ndim != 3 or a.shape[2] not in (3, 4):
            raise ValueError("image() takes an (h, w, 3|4) uint8 array")
        w_, s_, e_, n_ = edges or grid_edges(self.x, self.y)
        from pyproj import Transformer
        tr = Transformer.from_crs(self.crs, "EPSG:4326", always_xy=True)
        corners = [[round(v, 7) for v in tr.transform(px, py)]
                   for px, py in ((w_, n_), (e_, n_), (e_, s_), (w_, s_))]
        path = self.out / "images" / f"{lid}.png"
        path.parent.mkdir(parents=True, exist_ok=True)
        Image.fromarray(a).save(path, optimize=True)
        entry = self._base(lid, title, group, status, "", description, prov)
        entry.update({"kind": "image", "file": f"images/{lid}.png",
                      "shape": list(a.shape), "corners_lonlat": corners,
                      "bytes": path.stat().st_size})
        if extra:
            entry.update(extra)
        self.layers.append(entry)
        return entry

    def write_manifest(self, **meta) -> Path:
        self.meta.update(meta)
        w, s, e, n = grid_edges(self.x, self.y)
        grid = {"crs": self.crs, "width": int(len(self.x)), "height": int(len(self.y)),
                "x0": float(self.x[0]), "y0": float(self.y[0]),
                "dx": pixel_step(self.x), "dy": pixel_step(self.y),
                "edges_utm": [w, s, e, n], "corners_lonlat": self.corners}
        return write_json(self.out / "manifest.json",
                          {"grid": grid, **self.meta, "layers": self.layers}, indent=1)

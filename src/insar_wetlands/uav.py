"""UAV orthomosaics (Altum, RedEdge-MX Dual): reading onto a metric grid, and the shift between two flights.

Public code, private data: the rasters are read from the hub's ``06_data/field`` (``field.uav_rasters``); nothing
here holds a field value. The delivered GeoTIFFs are in geographic coordinates (EPSG:4326, some without the CRS
written in the file); they are put on a UTM 33N grid (the cube's EPSG:32633) before anything is measured in metres.
"""
from __future__ import annotations

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.transform import from_origin
from rasterio.warp import reproject, transform_bounds
from rasterio.windows import from_bounds as window_from_bounds

from insar_wetlands.field import UAV_CRS

GRID_CRS = 32633


def grid(bounds_utm: tuple[float, float, float, float], res: float):
    """A north-up grid on ``bounds_utm`` (W, S, E, N) snapped outward to ``res``: (transform, shape)."""
    w, s, e, n = bounds_utm
    w, s = np.floor(w / res) * res, np.floor(s / res) * res
    e, n = np.ceil(e / res) * res, np.ceil(n / res) * res
    return from_origin(w, n, res, res), (int(round((n - s) / res)), int(round((e - w) / res)))


def utm_bounds(path) -> tuple[float, float, float, float]:
    with rasterio.open(path) as r:
        return transform_bounds(r.crs or f"EPSG:{UAV_CRS}", f"EPSG:{GRID_CRS}", *r.bounds)


def read_on_grid(path, bands, nodata: float, dst_transform, dst_shape, factor: int = 1,
                 window_bounds_utm=None) -> np.ndarray:
    """Bands (1-based) of one orthomosaic averaged onto the UTM grid; NaN where no data.

    ``factor`` reads the source decimated (GDAL uses the file's overviews), for grids much coarser than the native
    pixel. ``window_bounds_utm`` reads only the part of the source under that box (full resolution work on a
    small area). Values equal to ``nodata`` and non-finite values never enter an average."""
    out = np.full((len(bands), *dst_shape), np.nan, dtype="float32")
    with rasterio.open(path) as r:
        src_crs = r.crs or f"EPSG:{UAV_CRS}"
        if window_bounds_utm is not None:
            wb = transform_bounds(f"EPSG:{GRID_CRS}", src_crs, *window_bounds_utm)
            win = window_from_bounds(*wb, transform=r.transform).round_offsets().round_lengths()
            win = win.intersection(rasterio.windows.Window(0, 0, r.width, r.height))
        else:
            win = rasterio.windows.Window(0, 0, r.width, r.height)
        h, w = max(1, int(win.height) // factor), max(1, int(win.width) // factor)
        src_tr = r.window_transform(win) @ rasterio.Affine.scale(win.width / w, win.height / h)
        for i, b in enumerate(bands):
            a = r.read(b, window=win, out_shape=(h, w), resampling=Resampling.average).astype("float32")
            a[(a == nodata) | ~np.isfinite(a)] = np.nan
            reproject(a, out[i], src_transform=src_tr, src_crs=src_crs, src_nodata=np.nan,
                      dst_transform=dst_transform, dst_crs=f"EPSG:{GRID_CRS}", dst_nodata=np.nan,
                      resampling=Resampling.average)
    return out


def edges(img: np.ndarray) -> np.ndarray:
    """Gradient magnitude, standardised: what two flights share (boardwalk, platforms, buildings) whatever the
    calibration, season or sensor. NaN → 0."""
    a = np.where(np.isfinite(img), img, np.nanmedian(img))
    gy, gx = np.gradient(a)
    g = np.hypot(gx, gy)
    g[~np.isfinite(img)] = 0.0
    sd = g.std()
    return (g - g.mean()) / sd if sd > 0 else g * 0


def phase_shift(ref: np.ndarray, mov: np.ndarray) -> tuple[float, float, float]:
    """Translation (rows, cols) that maps ``ref`` onto ``mov`` — ``mov[y, x] ≈ ref[y - dy, x - dx]`` — by phase
    correlation with a Hann taper and a parabolic sub-pixel peak; and the peak's height over the next-highest
    value outside its 5×5 neighbourhood (≫ 1: a clear match; ≈ 1: no common features)."""
    hy, hx = np.hanning(ref.shape[0]), np.hanning(ref.shape[1])
    taper = np.outer(hy, hx)
    fa, fb = np.fft.fft2(ref * taper), np.fft.fft2(mov * taper)
    cross = fb * np.conj(fa)
    cross /= np.maximum(np.abs(cross), 1e-12)
    c = np.fft.fftshift(np.real(np.fft.ifft2(cross)))
    py, px = np.unravel_index(np.argmax(c), c.shape)

    def sub(m1, m0, p1):
        d = m1 - 2 * m0 + p1
        return 0.0 if d == 0 else 0.5 * (m1 - p1) / d

    ny, nx = c.shape
    dy = py - ny // 2 + (sub(c[py - 1, px], c[py, px], c[py + 1, px]) if 0 < py < ny - 1 else 0.0)
    dx = px - nx // 2 + (sub(c[py, px - 1], c[py, px], c[py, px + 1]) if 0 < px < nx - 1 else 0.0)
    rest = c.copy()
    rest[max(0, py - 2):py + 3, max(0, px - 2):px + 3] = -np.inf
    return float(dy), float(dx), float(c[py, px] / max(rest.max(), 1e-12))


def ndvi(nir: np.ndarray, red: np.ndarray) -> np.ndarray:
    with np.errstate(invalid="ignore", divide="ignore"):
        return (nir - red) / (nir + red)

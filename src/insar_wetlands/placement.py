"""Where does a coarse grid's content really lie? (X-120, the HyP3 half-cell question; T75, T77).

A coarse product labels each cell by its centre L. If its values were computed for ground at L + (dx, dy) instead, the
shift shows when the product is compared with a finer, independently placed reference: the reference averaged over cells
centred at L + (dx, dy) matches best at the true (dx, dy). Pure functions; the data access lives in
``experiments/hyp3_half_cell_x120.py``.

Also the shift the HyP3 burst InSAR products are predicted to carry from how they are made (hyp3_isce2 3.0.2, ISCE2):
ISCE2's geocoder computes each value at ``lat0 + dlat·(line−1), lon0 + dlon·(pixel−1)`` and stores lat0/lon0 as the
image's coordStart; hyp3_isce2 passes coordStart to GDAL as the upper-left corner, so each value belongs to the corner of
its labelled geocoding cell. The geocoding grid is the DEM grid hyp3_isce2 prepares, square in degrees with side the
longitude degrees of 30 m (20 m for 5×1 looks) at the burst's southern latitude.
"""
from __future__ import annotations

import numpy as np

EARTH_CIRCUMFERENCE_M = 40_030_173.59204114   # hyp3_isce2 dem.distance_meters_to_degrees


def hyp3_isce2_predicted_shift(dem_pixel_m: float, lat_south_deg: float, lat_site_deg: float) -> tuple[float, float]:
    """(dx east, dy north) in metres at which the content of a HyP3 cell labelled L lies: half a geocoding cell west and
    north. ``dem_pixel_m``: 30 (10×2 and 20×4 looks) or 20 (5×1); ``lat_south_deg``: the latitude hyp3_isce2 converts at
    (the burst's southern edge); ``lat_site_deg``: where the metres are measured."""
    deg = dem_pixel_m / (EARTH_CIRCUMFERENCE_M * np.cos(np.radians(lat_south_deg))) * 360.0
    ew = deg * EARTH_CIRCUMFERENCE_M * np.cos(np.radians(lat_site_deg)) / 360.0
    ns = deg * EARTH_CIRCUMFERENCE_M / 360.0
    return -ew / 2.0, ns / 2.0


def _integral(a: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    v = np.where(np.isfinite(a), a, 0.0)
    n = np.isfinite(a).astype(float)
    pad = lambda m: np.pad(m.cumsum(0).cumsum(1), ((1, 0), (1, 0)))  # noqa: E731
    return pad(v), pad(n)


class FineField:
    """A fine reference raster (regular, north-up: x increasing, y decreasing) that can be averaged over boxes fast."""

    def __init__(self, values: np.ndarray, x: np.ndarray, y: np.ndarray):
        self.x0, self.dx = float(x[0]), float(x[1] - x[0])
        self.y0, self.dy = float(y[0]), float(y[1] - y[0])
        if self.dx <= 0 or self.dy >= 0:
            raise ValueError("FineField wants x increasing and y decreasing")
        self.nx, self.ny = len(x), len(y)
        self.S, self.N = _integral(np.asarray(values, float))

    def box_mean(self, cx: np.ndarray, cy: np.ndarray, half: float, min_cover: float = 0.9) -> np.ndarray:
        """Mean of the fine values over the squares of half-side ``half`` centred at (cx, cy); NaN where less than
        ``min_cover`` of the box lies inside the raster or holds finite values."""
        # fine pixel i covers [x0 + (i − ½)dx, x0 + (i + ½)dx]; take the pixels whose centre falls in the box
        i0 = np.ceil((cx - half - self.x0) / self.dx - 1e-9).astype(int)
        i1 = np.floor((cx + half - self.x0) / self.dx - 1e-9).astype(int) + 1
        j0 = np.ceil((self.y0 - (cy + half)) / -self.dy - 1e-9).astype(int)
        j1 = np.floor((self.y0 - (cy - half)) / -self.dy - 1e-9).astype(int) + 1
        full = (i1 - i0) * (j1 - j0)
        i0c, i1c = np.clip(i0, 0, self.nx), np.clip(i1, 0, self.nx)
        j0c, j1c = np.clip(j0, 0, self.ny), np.clip(j1, 0, self.ny)
        tot = lambda T: T[j1c, i1c] - T[j0c, i1c] - T[j1c, i0c] + T[j0c, i0c]  # noqa: E731
        s, n = tot(self.S), tot(self.N)
        with np.errstate(invalid="ignore", divide="ignore"):
            m = s / n
        return np.where(n >= min_cover * np.maximum(full, 1), m, np.nan)


def r2_ols(y: np.ndarray, X: np.ndarray) -> float:
    """R² of the least-squares fit of y on the columns of X (with an intercept), over the rows finite in both."""
    ok = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
    if ok.sum() < X.shape[1] + 5:
        return np.nan
    A = np.column_stack([np.ones(ok.sum()), X[ok]])
    beta, *_ = np.linalg.lstsq(A, y[ok], rcond=None)
    r = y[ok] - A @ beta
    return float(1.0 - r.var() / y[ok].var())


def corr(a: np.ndarray, b: np.ndarray) -> float:
    ok = np.isfinite(a) & np.isfinite(b)
    return float(np.corrcoef(a[ok], b[ok])[0, 1]) if ok.sum() > 5 else np.nan


def best_shift(score, steps: np.ndarray) -> tuple[float, float, float, np.ndarray]:
    """Grid search: ``score(dx, dy)`` → higher is better. Returns (dx, dy, best score, the score surface [dy, dx])."""
    surf = np.array([[score(dx, dy) for dx in steps] for dy in steps])
    k = np.unravel_index(np.nanargmax(surf), surf.shape)
    return float(steps[k[1]]), float(steps[k[0]]), float(surf[k]), surf

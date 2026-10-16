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

# The shift of the HyP3 burst products we hold (hyp3_isce2 3.0.2), by pixel spacing (m) and track: the content of a
# cell labelled L lies at L + (dx east, dy north). X-120's refined, pre-registered measurement against our right-placed
# 2021 SLC coherence (DECISIONS 2026-10-15); the 20 m products by two thirds (their geocoding cell).
HYP3_CONTENT_SHIFT_M: dict[int, dict[str, tuple[float, float]]] = {
    40: {"ascending": (-21.0, 23.0), "descending": (-16.0, 20.0)},
    20: {"ascending": (-14.0, 15.33), "descending": (-10.67, 13.33)},
}
HYP3_LAYER_KIND = {"wrapped_phase": "phase", "unw_phase": "mean", "corr": "mean", "conncomp": "mask",
                   "lv_theta": "mean", "lv_phi": "mean"}       # the DEM layer is placed right and is not moved

# The shift of the RTC backscatter on the 40 m grid (Planetary Computer sentinel-1-rtc, nearest-resampled by
# pipelines/rtc_download.py), by track: X-122's pre-registered same-day measurement against our right-placed 2021 CSLC
# amplitude (DECISIONS 2026-10-16). Data-fitted; the mechanism is not identified.
RTC_CONTENT_SHIFT_M: dict[str, tuple[float, float]] = {"ascending": (11.0, -19.0), "descending": (17.0, -8.0)}
RTC_RELOCATE_ENV = "RZECIN_RTC_RELOCATE"


def rtc_relocation_on() -> bool:
    """Off by default, so finished entries reproduce from the RTC stacks as delivered; ``RZECIN_RTC_RELOCATE=1`` makes
    every stack opened through ``open_rtc`` relocated under its labels. The cube is built relocated either way."""
    import os
    return os.environ.get(RTC_RELOCATE_ENV, "") not in ("", "0")


def relocate_rtc(ds, track: str):
    """An RTC dataset (the Drive stacks: ``gamma0_vv_db``, ``gamma0_vh_db``, ``ratio_vh_vv_db``; or the cube's
    ``vv_db``, ``vh_db``, ``vh_vv_db``, ``rvi``) with its content put back under its labels. Power is relocated in
    linear units, then the dB values, the ratio and the RVI are recomputed from it. Other variables are kept."""
    dx, dy = RTC_CONTENT_SHIFT_M[track]
    x, y = ds.x.values, ds.y.values
    out = ds.copy()
    names = ({"vv": "gamma0_vv_db", "vh": "gamma0_vh_db", "ratio": "ratio_vh_vv_db"} if "gamma0_vv_db" in ds
             else {"vv": "vv_db", "vh": "vh_db", "ratio": "vh_vv_db"})
    lin = {}
    for k in ("vv", "vh"):
        v = names[k]
        if v not in ds:
            continue
        p = 10 ** (ds[v].values.astype(float) / 10)
        lin[k] = relocate(p, x, y, dx, dy)
        out[v] = (ds[v].dims, (10 * np.log10(np.where(lin[k] > 0, lin[k], np.nan))).astype(ds[v].dtype))
    if "vv" in lin and "vh" in lin:
        if names["ratio"] in ds:
            out[names["ratio"]] = (ds[names["ratio"]].dims,
                                   (out[names["vh"]].values - out[names["vv"]].values).astype(ds[names["ratio"]].dtype))
        if "rvi" in ds:
            out["rvi"] = (ds.rvi.dims, (4 * lin["vh"] / (lin["vv"] + lin["vh"])).astype(ds.rvi.dtype))
    out.attrs = dict(ds.attrs, relocated=f"X-122: content moved back from L + ({dx:+.0f}, {dy:+.0f}) m")
    return out


def open_rtc(path, track: str | None = None):
    """Open an RTC stack; relocated (``relocate_rtc``) when ``RZECIN_RTC_RELOCATE=1``. ``track`` defaults to the file's
    ``track`` attribute (the ascending 2022–2024 ``rtc_dualpol_stack.nc`` has none: pass it)."""
    import xarray as xr
    ds = xr.open_dataset(path)
    if not rtc_relocation_on() or ds.attrs.get("relocated", "").startswith("X-122"):
        return ds
    return relocate_rtc(ds.load(), track or ds.attrs["track"])


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


def relocate(values: np.ndarray, x: np.ndarray, y: np.ndarray, dx: float, dy: float, kind: str = "mean",
             min_weight: float = 0.5) -> np.ndarray:
    """Put a gridded product's content back where it belongs. ``values`` (..., y, x) on a regular north-up grid whose
    cell labelled L holds content for ground at L + (dx, dy); returns, on the same labels, the content for the ground
    under each label — bilinear interpolation of the source at L − (dx, dy).

    ``kind``: "mean" (coherence, unwrapped phase, angles, heights), "phase" (wrapped phase, interpolated as a phasor) or
    "mask" (0/1 layers: interpolated, then ≥ 0.5). Neighbours that are NaN drop out and the weights are renormalised;
    a cell is NaN when less than ``min_weight`` of its weight is finite, or when it falls outside the source grid."""
    v = np.asarray(values, float)
    gx, gy = float(x[1] - x[0]), float(y[1] - y[0])
    fc, fr = -dx / gx, -dy / gy                    # source index offset of each target cell (columns, rows)
    c0, r0 = int(np.floor(fc)), int(np.floor(fr))
    wc, wr = fc - c0, fr - r0
    ny, nx = v.shape[-2:]

    def take(a, dr, dc):
        out = np.full(a.shape, np.nan)
        rs, cs = r0 + dr, c0 + dc
        ys, ye = max(0, -rs), min(ny, ny - rs)
        xs, xe = max(0, -cs), min(nx, nx - cs)
        if ys < ye and xs < xe:
            out[..., ys:ye, xs:xe] = a[..., ys + rs:ye + rs, xs + cs:xe + cs]
        return out

    if kind == "phase":
        parts = [np.cos(v), np.sin(v)]
    elif kind in ("mean", "mask"):
        parts = [v]
    else:
        raise ValueError(kind)
    weights = {(0, 0): (1 - wr) * (1 - wc), (0, 1): (1 - wr) * wc, (1, 0): wr * (1 - wc), (1, 1): wr * wc}
    num = [np.zeros(v.shape) for _ in parts]
    den = np.zeros(v.shape)
    for (dr, dc), w in weights.items():
        if w == 0:
            continue
        sh = [take(p, dr, dc) for p in parts]
        ok = np.isfinite(sh[0])
        for k, s in enumerate(sh):
            num[k] += np.where(ok, w * s, 0.0)
        den += np.where(ok, w, 0.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        if kind == "phase":
            out = np.arctan2(num[1], num[0])
        else:
            out = num[0] / den
    out[den < min_weight] = np.nan
    if kind == "mask":
        out = np.where(np.isfinite(out), (out >= 0.5).astype(float), np.nan)
    return out

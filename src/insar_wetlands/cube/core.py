"""Pure building blocks of the data cube: time joins, aggregation across grids, geometry, flags, folds.

Every function here is deterministic and tested on synthetic data (``tests/test_cube_core.py``).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from ..regrid import cell_index

# Sentinel-1 C-band nominal wavelength. C-052: ``isbas.WAVELENGTH_M`` carries 0.0555465763 m (likely a digit
# transposition); the cube stores phase in radians, so the choice only enters through this one constant.
WAVELENGTH_MM = 55.465763
OVERPASS_UTC = {"ascending": "16:36", "descending": "05:09"}


def rad_to_los_mm(phase_rad, wavelength_mm: float = WAVELENGTH_MM):
    """Interferometric phase (rad) → line-of-sight change (mm): λ/(4π) per radian (HyP3 sign kept as delivered)."""
    return np.asarray(phase_rad, float) * wavelength_mm / (4 * np.pi)


def overpass_time(date, track: str) -> pd.Timestamp:
    """The overpass instant (UTC) of a Sentinel-1 acquisition date on a track."""
    return pd.Timestamp(f"{pd.Timestamp(date):%Y-%m-%d} {OVERPASS_UTC[track]}", tz="UTC")


# ------------------------------------------------------------------------------ time

def _ns(idx) -> np.ndarray:
    """Nanoseconds since the epoch, whatever resolution pandas stores the index in."""
    return pd.DatetimeIndex(idx).as_unit("ns").asi8

def asof(times, table: pd.DataFrame, tolerance: str | pd.Timedelta, direction: str = "nearest",
         prefix: str = "") -> pd.DataFrame:
    """Values of ``table`` (UTC DatetimeIndex) at each of ``times``: the row nearest in time (or the last before /
    first after, ``direction``) within ``tolerance``; ``lag_h`` = row time − query time (hours, negative = before).
    Rows outside the tolerance are NaN — a missing value, never a silently distant one."""
    q = pd.DataFrame({"_t": pd.DatetimeIndex(pd.to_datetime(times, utc=True)).as_unit("ns")})
    q["_order"] = np.arange(len(q))
    t = table.sort_index().copy()
    t["_src_t"] = pd.DatetimeIndex(pd.to_datetime(t.index, utc=True)).as_unit("ns")
    t = t.reset_index(drop=True)
    m = pd.merge_asof(q.sort_values("_t"), t, left_on="_t", right_on="_src_t", direction=direction,
                      tolerance=pd.Timedelta(tolerance)).sort_values("_order")
    lag = (m["_src_t"] - m["_t"]).dt.total_seconds() / 3600
    out = m.drop(columns=["_t", "_order", "_src_t"]).reset_index(drop=True)
    out["lag_h"] = lag.to_numpy()
    return out.add_prefix(prefix)


def value_at(s: pd.Series, times, max_gap: str = "2h") -> np.ndarray:
    """Linear interpolation of an hourly series at arbitrary instants; NaN when either neighbour is missing or the
    neighbours are more than ``max_gap`` apart."""
    s = s.dropna().sort_index()
    t = pd.DatetimeIndex(pd.to_datetime(times, utc=True))
    if s.empty:
        return np.full(len(t), np.nan)
    x = _ns(s.index).astype(float)
    q = _ns(t).astype(float)
    j = np.searchsorted(x, q)
    out = np.full(len(t), np.nan)
    ok = (j > 0) & (j < len(x))
    jj = j[ok]
    x0, x1 = x[jj - 1], x[jj]
    gap_ok = (x1 - x0) <= pd.Timedelta(max_gap).value
    w = np.where(x1 > x0, (q[ok] - x0) / np.where(x1 > x0, x1 - x0, 1), 0.0)
    v = s.to_numpy(float)
    val = v[jj - 1] * (1 - w) + v[jj] * w
    exact = np.isin(q, x)
    out[np.flatnonzero(ok)[gap_ok]] = val[gap_ok]
    if exact.any():
        out[exact] = v[np.searchsorted(x, q[exact])]
    return out


def window_features(s: pd.Series, times, name: str, windows_h=(24, 72, 168), change_h=(72, 168),
                    min_frac: float = 0.75) -> pd.DataFrame:
    """The field team's antecedent set for one hourly series at each instant: value at the time, mean over the
    preceding ``windows_h`` hours, change over the preceding ``change_h`` hours. A window with less than ``min_frac``
    of its hours present is NaN."""
    s = s.sort_index()
    t = pd.DatetimeIndex(pd.to_datetime(times, utc=True))
    out = {f"{name}_at": value_at(s, t)}
    hourly = s.resample("1h").mean()
    for h in windows_h:
        m = hourly.rolling(f"{h}h", min_periods=max(1, int(h * min_frac))).mean()
        out[f"{name}_mean{h}h"] = asof(t, m.to_frame("v"), "1h", direction="backward")["v"].to_numpy()
    for h in change_h:
        before = value_at(s, t - pd.Timedelta(hours=h))
        out[f"{name}_d{h}h"] = out[f"{name}_at"] - before
    return pd.DataFrame(out)


def hours_since(event: pd.Series, times) -> np.ndarray:
    """Hours from the last True of a boolean hourly series to each instant (NaN before the first event)."""
    ev = event[event.fillna(False).astype(bool)].index
    t = pd.DatetimeIndex(pd.to_datetime(times, utc=True))
    if len(ev) == 0:
        return np.full(len(t), np.nan)
    e, q = _ns(ev), _ns(t)
    j = np.searchsorted(e, q, side="right") - 1
    out = np.full(len(t), np.nan)
    ok = j >= 0
    out[ok] = (q[ok] - e[j[ok]]) / 3.6e12
    return out


# ------------------------------------------------------------------------------ meteorology

def dew_point_c(temp_c, rh_pct):
    """Dew point (°C) from air temperature and relative humidity (Magnus, Sonntag constants)."""
    t = np.asarray(temp_c, float)
    rh = np.clip(np.asarray(rh_pct, float), 1e-3, 100)
    g = np.log(rh / 100) + 17.62 * t / (243.12 + t)
    return 243.12 * g / (17.62 - g)


# ------------------------------------------------------------------------------ space

def aggregate_to_grid(values: np.ndarray, fine_x, fine_y, tpl_x, tpl_y) -> dict[str, np.ndarray]:
    """A fine field (y, x) onto a coarser grid: per cell mean, standard deviation and valid fraction (finite fine
    pixels over all fine pixels whose centre falls in the cell). Cells with no fine pixel are NaN."""
    v = np.asarray(values, float)
    r, c = cell_index(np.asarray(fine_x, float), np.asarray(fine_y, float), np.asarray(tpl_x, float),
                      np.asarray(tpl_y, float))
    R, C = np.meshgrid(r, c, indexing="ij")
    inside = (R >= 0) & (C >= 0)
    n = len(tpl_y) * len(tpl_x)
    flat_all = (R * len(tpl_x) + C)[inside]
    total = np.bincount(flat_all, minlength=n).astype(float)
    ok = inside & np.isfinite(v)
    flat = (R * len(tpl_x) + C)[ok]
    cnt = np.bincount(flat, minlength=n).astype(float)
    s1 = np.bincount(flat, weights=v[ok], minlength=n)
    s2 = np.bincount(flat, weights=v[ok] ** 2, minlength=n)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = s1 / cnt
        var = np.maximum(s2 / cnt - mean ** 2, 0) * cnt / np.maximum(cnt - 1, 1)
        frac = cnt / total
    mean[cnt == 0] = np.nan
    var[cnt < 2] = np.nan
    frac[total == 0] = np.nan
    shape = (len(tpl_y), len(tpl_x))
    return {"mean": mean.reshape(shape), "sd": np.sqrt(var).reshape(shape), "valid_frac": frac.reshape(shape)}


def class_fractions(classes: np.ndarray, fine_x, fine_y, tpl_x, tpl_y, codes) -> dict[int, np.ndarray]:
    """Fraction of each class code per coarse cell (over the fine pixels whose centre falls in it)."""
    cl = np.asarray(classes)
    return {int(k): aggregate_to_grid(np.where(np.isfinite(cl.astype(float)), (cl == k).astype(float), np.nan),
                                      fine_x, fine_y, tpl_x, tpl_y)["mean"] for k in codes}


def window_stats(field2d: np.ndarray, row: int, col: int, half: int = 1) -> dict:
    """Median, mean, SD, IQR and valid count of a (2·half+1)² window centred on (row, col); edges are clipped."""
    a = np.asarray(field2d, float)[max(row - half, 0):row + half + 1, max(col - half, 0):col + half + 1]
    v = a[np.isfinite(a)]
    if not v.size:
        return {"n": 0, "median": np.nan, "mean": np.nan, "sd": np.nan, "iqr": np.nan}
    q1, q3 = np.percentile(v, [25, 75])
    return {"n": int(v.size), "median": float(np.median(v)), "mean": float(v.mean()),
            "sd": float(v.std(ddof=1)) if v.size > 1 else np.nan, "iqr": float(q3 - q1)}


def phasor_mean(phase: np.ndarray, weight: np.ndarray | None = None) -> tuple[float, float]:
    """Circular mean of wrapped phases (rad) and its consistency |mean phasor| (1 = identical phases)."""
    p = np.asarray(phase, float).ravel()
    w = np.ones_like(p) if weight is None else np.asarray(weight, float).ravel()
    ok = np.isfinite(p) & np.isfinite(w)
    if not ok.any() or w[ok].sum() <= 0:
        return np.nan, np.nan
    z = np.sum(w[ok] * np.exp(1j * p[ok])) / np.sum(w[ok])
    return float(np.angle(z)), float(np.abs(z))


# ------------------------------------------------------------------------------ geometry

def incidence_deg(lv_theta_rad):
    """Incidence angle (deg) from HyP3's ``lv_theta`` (elevation angle of the look vector, rad)."""
    return np.degrees(np.pi / 2 - np.asarray(lv_theta_rad, float))


def los_enu(lv_theta_rad, lv_phi_rad) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """East, north, up components of the unit vector from the ground to the satellite (HyP3 ``lv_theta`` = elevation
    from horizontal, ``lv_phi`` = azimuth from east, counter-clockwise). A vertical uplift dh projects to +dh·up."""
    th, ph = np.asarray(lv_theta_rad, float), np.asarray(lv_phi_rad, float)
    return np.cos(th) * np.cos(ph), np.cos(th) * np.sin(ph), np.sin(th)


def cycle_risk(dlos_mm, wavelength_mm: float = WAVELENGTH_MM) -> np.ndarray:
    """True where a line-of-sight change exceeds a quarter wavelength: phase unwrapping cannot follow it."""
    return np.abs(np.asarray(dlos_mm, float)) > wavelength_mm / 4


# ------------------------------------------------------------------------------ folds

def year_folds(times) -> np.ndarray:
    """Calendar year of each instant: the leave-one-year-out fold id."""
    return pd.DatetimeIndex(pd.to_datetime(times, utc=True)).year.to_numpy()


def pair_year_fold(t1, t2) -> np.ndarray:
    """A pair's fold: the year of its first date if both dates share it, else −1 (a pair straddling two folds must be
    left out of both when testing on one of them)."""
    y1, y2 = year_folds(t1), year_folds(t2)
    return np.where(y1 == y2, y1, -1)


def spatial_blocks(shape: tuple[int, int], block: int) -> np.ndarray:
    """Block id of every cell for blocked spatial cross-validation (square blocks of ``block`` cells)."""
    r, c = np.indices(shape)
    nb_c = -(-shape[1] // block)
    return (r // block) * nb_c + (c // block)

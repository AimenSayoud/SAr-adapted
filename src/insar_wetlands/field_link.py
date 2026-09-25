"""Sentinel-1 ↔ field observations: the statistics behind X-047 (WTD) and X-048 (laser at P6).

Two lessons of this project are built in:

* **Shared seasonal cycles are not evidence.** Every correlation is reported twice — with only a
  trend removed (``r_raw``) and with the annual harmonic removed from both series (``r_anom``,
  the same ``hydro_link._detrend`` as the site-scale T09). Only ``r_anom`` tests a coupling.
* **Time series are autocorrelated.** Significance comes from circular shifts of one series
  against the other (``circular_shift_p``), which keeps each series' autocorrelation; p has a
  floor of 1/(N+1).

Field data are read by the caller (``insar_wetlands.field``); nothing here holds data.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .hydro_link import _detrend


def circular_shift_p(x: np.ndarray, y: np.ndarray, n_shift: int = 2000, min_shift: int = 3,
                     rng: np.random.Generator | None = None) -> tuple[float, float]:
    """(r, p): Pearson r of x and y, and the two-sided p of |r| against x circularly shifted by
    ``min_shift`` … N − ``min_shift`` samples. Floor 1/(n_shift + 1)."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    n = len(x)
    r0 = float(np.corrcoef(x, y)[0, 1])
    if n < 2 * min_shift + 2:
        return r0, float("nan")
    rng = rng or np.random.default_rng(0)
    shifts = rng.integers(min_shift, n - min_shift, size=n_shift)
    null = np.array([np.corrcoef(np.roll(x, k), y)[0, 1] for k in shifts])
    return r0, float((np.sum(np.abs(null) >= abs(r0)) + 1) / (n_shift + 1))


def years_since(dates) -> np.ndarray:
    d = pd.to_datetime(pd.Series(dates)).reset_index(drop=True)
    return ((d - d.iloc[0]).dt.total_seconds() / (365.25 * 86400)).to_numpy()


def anomaly_correlation(dates, x, y, **kw) -> dict:
    """r with trends removed, and r / p with trend + annual cycle removed from both series."""
    df = pd.DataFrame({"d": pd.to_datetime(pd.Series(dates)).to_numpy(), "x": x, "y": y}).dropna().sort_values("d")
    if len(df) < 12:
        return {"n": len(df), "r_raw": np.nan, "r_anom": np.nan, "p_anom": np.nan}
    t = years_since(df.d)
    xv, yv = df.x.to_numpy(float), df.y.to_numpy(float)
    r_raw = float(np.corrcoef(_detrend(xv, t), _detrend(yv, t))[0, 1])
    r_anom, p = circular_shift_p(_detrend(xv, t, True), _detrend(yv, t, True), **kw)
    return {"n": len(df), "r_raw": r_raw, "r_anom": r_anom, "p_anom": p}


def season_residual(values, mid_dates, dt_days=None) -> np.ndarray:
    """Pair-level values with an annual harmonic of the pair's mid-date (and optionally a linear
    temporal-baseline term) regressed out — the pair analogue of ``_detrend(..., True)``."""
    v = np.asarray(values, float)
    doy = pd.to_datetime(pd.Series(mid_dates)).dt.dayofyear.to_numpy() / 365.25 * 2 * np.pi
    cols = [np.ones_like(v), np.cos(doy), np.sin(doy)]
    if dt_days is not None:
        cols.append(np.asarray(dt_days, float))
    M = np.column_stack(cols)
    b, *_ = np.linalg.lstsq(M, v, rcond=None)
    return v - M @ b


def flag_shift_test(values, members, flags, n_shift: int = 2000, min_shift: int = 3,
                    rng: np.random.Generator | None = None) -> dict:
    """Do items that involve a flagged date differ from items that do not?

    ``values``: one value per item (a pair, or a date); ``members``: (n_items, k) indices into
    ``flags`` (the dates of each item; k = 2 for pairs, 1 for dates); ``flags``: one bool per
    date in time order. Statistic: mean over items with no flagged member − mean over items
    with at least one. Null: the flag sequence circularly shifted along the dates, which keeps
    how flags cluster in time (dew is seasonal) but breaks their link to the values; two-sided
    p with floor 1/(n_shift + 1).
    """
    v = np.asarray(values, float)
    m = np.atleast_2d(np.asarray(members, int))
    m = m.T if m.shape[0] != len(v) else m
    f = np.asarray(flags, bool)

    def stat(ff):
        hit = ff[m].any(axis=1)
        if hit.all() or not hit.any():
            return np.nan
        return float(v[~hit].mean() - v[hit].mean())

    d0 = stat(f)
    hit = f[m].any(axis=1)
    out = {"n_clear": int((~hit).sum()), "n_flagged": int(hit.sum()), "diff": d0, "p": np.nan}
    n = len(f)
    if np.isnan(d0) or n < 2 * min_shift + 2:
        return out
    rng = rng or np.random.default_rng(0)
    null = np.array([stat(np.roll(f, k)) for k in rng.integers(min_shift, n - min_shift, size=n_shift)])
    null = null[np.isfinite(null)]
    out["p"] = float((np.sum(np.abs(null) >= abs(d0)) + 1) / (len(null) + 1))
    return out


def los_from_vertical(dh_mm, incidence_deg: float):
    """LOS change (mm, toward the satellite positive) of a purely vertical surface change dh
    (up positive): dh · cos(incidence)."""
    return np.asarray(dh_mm, float) * np.cos(np.radians(incidence_deg))


def value_at(series: pd.Series, t: pd.Timestamp, tolerance: str = "1h") -> float:
    """Nearest non-null value of an hourly series within ``tolerance`` of t, else NaN."""
    s = series.dropna()
    if s.empty:
        return np.nan
    i = s.index.get_indexer([t], method="nearest")[0]
    return float(s.iloc[i]) if abs(s.index[i] - t) <= pd.Timedelta(tolerance) else np.nan

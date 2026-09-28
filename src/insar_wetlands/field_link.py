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


def circular_shift_p_many(x, Y, n_shift: int = 1000, min_shift: int = 3, chunk: int = 100,
                          rng: np.random.Generator | None = None) -> tuple[np.ndarray, np.ndarray]:
    """``circular_shift_p`` for many series at once: r of x with every column of Y (e.g. every
    pixel) and each column's two-sided p against x circularly shifted (the same shifts for all
    columns). Columns with a NaN or no variance give NaN. The null is counted chunk by chunk,
    never held whole."""
    x = np.asarray(x, float)
    Y = np.asarray(Y, float)
    n = len(x)
    xs = (x - x.mean()) / x.std()
    sd = Y.std(axis=0)
    bad = ~np.isfinite(Y).all(axis=0) | ~(sd > 0)
    with np.errstate(invalid="ignore", divide="ignore"):
        Ys = np.nan_to_num((Y - Y.mean(axis=0)) / np.where(sd > 0, sd, 1.0))
    r = xs @ Ys / n
    rng = rng or np.random.default_rng(0)
    ks = rng.integers(min_shift, n - min_shift, size=n_shift)
    exceed = np.zeros(Y.shape[1])
    for i in range(0, n_shift, chunk):
        null = np.stack([np.roll(xs, k) for k in ks[i:i + chunk]]) @ Ys / n
        exceed += np.sum(np.abs(null) >= np.abs(r) - 1e-12, axis=0)
    p = (exceed + 1) / (n_shift + 1)
    r[bad], p[bad] = np.nan, np.nan
    return r, p


def bh_qvalues(p) -> np.ndarray:
    """Benjamini–Hochberg q-values (false-discovery rate) of a set of p-values; NaN stays NaN."""
    p = np.asarray(p, float)
    q = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    v = p[ok]
    if v.size:
        order = np.argsort(v)
        ranked = v[order] * v.size / (np.arange(v.size) + 1)
        qq = np.minimum(np.minimum.accumulate(ranked[::-1])[::-1], 1)
        out = np.empty_like(v)
        out[order] = qq
        q[ok] = out
    return q


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


def flag_split_correlation(x, y, members, flags, n_shift: int = 2000, min_shift: int = 3,
                           min_n: int = 12, rng: np.random.Generator | None = None) -> dict:
    """Is the x–y correlation different on items with no flagged date than on the rest?

    Returns r on the clear items, r on the flagged items, their difference and a two-sided p
    from the flags circularly shifted along the dates (as ``flag_shift_test``). Groups smaller
    than ``min_n`` give NaN.
    """
    x, y = np.asarray(x, float), np.asarray(y, float)
    m = np.atleast_2d(np.asarray(members, int))
    m = m.T if m.shape[0] != len(x) else m
    f = np.asarray(flags, bool)

    def rs(ff):
        hit = ff[m].any(axis=1)
        if (~hit).sum() < min_n or hit.sum() < min_n:
            return np.nan, np.nan
        return float(np.corrcoef(x[~hit], y[~hit])[0, 1]), float(np.corrcoef(x[hit], y[hit])[0, 1])

    r_clear, r_flag = rs(f)
    hit = f[m].any(axis=1)
    out = {"n_clear": int((~hit).sum()), "n_flagged": int(hit.sum()), "r_clear": r_clear,
           "r_flagged": r_flag, "diff": r_clear - r_flag, "p": np.nan}
    n = len(f)
    if not np.isfinite(out["diff"]) or n < 2 * min_shift + 2:
        return out
    rng = rng or np.random.default_rng(0)
    null = []
    for k in rng.integers(min_shift, n - min_shift, size=n_shift):
        a, b = rs(np.roll(f, k))
        null.append(a - b)
    null = np.asarray(null)
    null = null[np.isfinite(null)]
    out["p"] = float((np.sum(np.abs(null) >= abs(out["diff"])) + 1) / (len(null) + 1))
    return out


def _all_permutations(n: int) -> np.ndarray:
    """Every ordering of n items (cached: 9! rows are built once per process)."""
    from itertools import permutations
    if n not in _PERMS:
        _PERMS[n] = np.array(list(permutations(range(n))), dtype=np.int8)
    return _PERMS[n]


_PERMS: dict[int, np.ndarray] = {}


def spearman_exact(x, y, max_exact: int = 9, n_perm: int = 20000,
                   rng: np.random.Generator | None = None) -> tuple[float, float, int]:
    """(rho, two-sided p, n) of Spearman's rank correlation between plots, with p from ALL
    permutations when n ≤ ``max_exact`` (9! = 362 880) — no normal approximation, which is
    meaningless for a handful of plots — else from ``n_perm`` random permutations. NaN pairs out."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    rx, ry = pd.Series(x[ok]).rank().to_numpy(), pd.Series(y[ok]).rank().to_numpy()
    n = len(rx)
    if n < 4 or rx.std() == 0 or ry.std() == 0:
        return np.nan, np.nan, n
    rx, ry = (rx - rx.mean()) / rx.std(), (ry - ry.mean()) / ry.std()
    rho = float(rx @ ry / n)
    if n <= max_exact:
        null = ry[_all_permutations(n)] @ rx / n
        return rho, float(np.mean(np.abs(null) >= abs(rho) - 1e-12)), n
    rng = rng or np.random.default_rng(0)
    null = np.array([ry[rng.permutation(n)] @ rx / n for _ in range(n_perm)])
    return rho, float((np.sum(np.abs(null) >= abs(rho) - 1e-12) + 1) / (n_perm + 1)), n


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


# ------------------------------------------------------------------ laser checks (X-058)

def running_median_outliers(s: pd.Series, window: str = "25h", k: float = 5.0,
                            min_abs: float = 1.0) -> tuple[pd.Series, float]:
    """Flag values far from a centred running median: |value − median| > max(k · σ, ``min_abs``),
    with σ the robust scale (1.4826 · MAD) of all residuals — or, when a coarse resolution makes
    the MAD zero, the SD of the residuals below ``min_abs``. NaN stays unflagged. Returns the
    flags (aligned to ``s``) and σ."""
    v = s.dropna()
    med = v.rolling(window, center=True, min_periods=3).median()
    res = v - med
    mad = float((res - res.median()).abs().median())
    sigma = 1.4826 * mad
    if sigma == 0:                       # coarse resolution: most residuals are exactly 0
        sigma = float(res[res.abs() < min_abs].std())
    flags = (res.abs() > max(k * sigma, min_abs)).reindex(s.index, fill_value=False)
    return flags.astype(bool), sigma


def hourly_steps(s: pd.Series, threshold: float, max_gap: str = "1h") -> pd.DataFrame:
    """Jumps between consecutive valid values no more than ``max_gap`` apart whose size exceeds
    ``threshold`` (a sensor reset or re-levelling shows up as one; slow motion does not)."""
    v = s.dropna()
    d = v.diff()
    gap = v.index.to_series().diff()
    hit = (d.abs() > threshold) & (gap <= pd.Timedelta(max_gap))
    return pd.DataFrame({"time": v.index[hit], "step": d[hit].to_numpy(), "before": v.shift(1)[hit].to_numpy(),
                         "after": v[hit].to_numpy()})


def diff_noise_sd(s: pd.Series, clip: float, max_gap: str = "1h") -> float:
    """SD of hour-to-hour differences (|Δ| < ``clip``, consecutive samples only): the noise of a
    difference of two readings, when the true change within an hour is negligible (an upper
    bound otherwise). The clip keeps jumps out; plain SD, because a coarse resolution makes most
    differences exactly zero and a MAD or IQR would then report no noise at all."""
    v = s.dropna()
    d = v.diff()[v.index.to_series().diff() <= pd.Timedelta(max_gap)]
    return float(d[d.abs() < clip].std())


def gap_runs(present: pd.Series, min_hours: int = 24) -> pd.DataFrame:
    """Stretches of an hourly boolean ``present`` series that are False for ≥ ``min_hours``."""
    miss = ~present.astype(bool)
    grp = (miss != miss.shift()).cumsum()
    rows = [(g.index[0], g.index[-1], len(g)) for _, g in miss.groupby(grp) if g.iloc[0] and len(g) >= min_hours]
    return pd.DataFrame(rows, columns=["start", "end", "hours"])


def slopes(x, y) -> dict:
    """Slope of y on x by ordinary least squares and by total least squares (orthogonal
    regression — errors in both variables, same units), with n and r."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if len(x) < 3:
        return {"n": len(x), "r": np.nan, "slope_ols": np.nan, "slope_tls": np.nan}
    C = np.cov(x, y)
    w, v = np.linalg.eigh(C)
    major = v[:, np.argmax(w)]
    return {"n": len(x), "r": float(np.corrcoef(x, y)[0, 1]), "slope_ols": float(C[0, 1] / C[0, 0]),
            "slope_tls": float(major[1] / major[0]) if major[0] != 0 else np.nan}


def interpolated_runs(s: pd.Series, min_hours: int = 6, grid: float = 0.05, tol: float = 1e-3) -> pd.DataFrame:
    """Stretches that are straight lines between two readings (gap filling): consecutive hourly
    changes that are non-zero, constant to within ``tol`` and **not multiples of the sensor's
    grid** (``grid``; 0.05 covers readings on a 0.1 grid and their averages). A sensor reports
    grid values, so its changes are grid multiples; a line drawn between two readings is not —
    whether it bridges a gap slowly or ramps across a re-levelling in a few hours.
    Returns start / end (inclusive, the two anchor readings excluded) and hours."""
    v = s.dropna()
    d = v.diff()
    k = d / grid
    off_grid = (k - np.round(k)).abs() > 0.02
    lin = (d.diff().abs() < tol) & (d.abs() > 1e-6) & off_grid
    grp = (lin != lin.shift()).cumsum()
    rows = []
    for _, g in lin.groupby(grp):
        if g.iloc[0] and len(g) + 1 >= min_hours:
            i0 = v.index.get_loc(g.index[0]) - 1          # the run's first filled value
            rows.append((v.index[i0], g.index[-1], len(g) + 1))
    return pd.DataFrame(rows, columns=["start", "end", "hours"])


def in_runs(index: pd.DatetimeIndex, runs: pd.DataFrame) -> np.ndarray:
    """True for times inside any [start, end] of ``runs``."""
    out = np.zeros(len(index), bool)
    for r in runs.itertuples(index=False):
        out |= (index >= r.start) & (index <= r.end)
    return out


def grid_offsets(s: pd.Series, step: float = 0.1, half: bool = True, min_days: int = 5) -> pd.DataFrame:
    """Where the readings' offset from the sensor's grid changes. A sensor reports multiples of
    ``step``; a constant added in post-processing (re-levelling after a reset or a move) shifts
    every later value off that grid by the same amount. Per day, the dominant offset
    (value mod ``step``, folded to ``step``/2 if ``half`` — averaged readings land on half steps);
    returns the stretches of at least ``min_days`` days with one offset."""
    v = s.dropna()
    unit = step / 2 if half else step
    off = np.round(np.mod(v.to_numpy(), unit), 4)
    off = np.where(np.isclose(off, unit, atol=1e-4), 0.0, off)
    day = pd.Series(off, index=v.index).groupby(v.index.floor("D")).agg(lambda x: x.mode().iloc[0])
    grp = (day != day.shift()).cumsum()
    rows = [(g.index[0], g.index[-1], len(g), g.iloc[0]) for _, g in day.groupby(grp) if len(g) >= min_days]
    return pd.DataFrame(rows, columns=["first_day", "last_day", "days", "offset"])


# ------------------------------------------------------------------ P6 validation (X-059)

def consecutive_pairs(dates) -> list[tuple[pd.Timestamp, pd.Timestamp]]:
    """(t_i, t_i+1) for sorted acquisition dates: the observations between consecutive
    acquisitions, which do not overlap (unlike a network of 12- and 24-day pairs)."""
    d = sorted(pd.to_datetime(pd.Series(dates)).unique())
    return list(zip(d[:-1], d[1:]))


def chain(increments) -> tuple[np.ndarray, np.ndarray]:
    """Cumulative sum of pair increments into a series at the acquisition dates (first = 0).
    A missing increment (NaN) breaks the chain: the next valid increment starts a new segment
    at 0 on its first date — levels are never carried across a break; dates reached by no valid
    increment stay NaN.
    Returns one value per date (len(increments) + 1) and the segment id of each date."""
    inc = np.asarray(increments, float)
    out = np.full(len(inc) + 1, np.nan)
    seg = np.full(len(inc) + 1, -1)
    s, level, open_ = 0, 0.0, True
    out[0], seg[0] = 0.0, 0
    for k, v in enumerate(inc):
        if np.isfinite(v):
            if not open_:
                s += 1
                level = 0.0
                out[k], seg[k] = 0.0, s
                open_ = True
            level += v
            out[k + 1], seg[k + 1] = level, s
        else:
            open_ = False
    return out, seg


def harmonic_amplitude(dates, values, segments=None) -> dict:
    """Semi-amplitude of an annual harmonic fitted jointly with a linear trend and one intercept
    per segment (a chain break or a sensor re-referencing leaves an unknown offset, which must
    not be read as seasonal motion). Returns amplitude, day of year of the maximum and n."""
    v = np.asarray(values, float)
    d = pd.to_datetime(pd.Series(dates)).reset_index(drop=True)
    seg = np.zeros(len(v), int) if segments is None else np.asarray(segments, int)
    ok = np.isfinite(v) & (seg >= 0)
    if ok.sum() < 8:
        return {"n": int(ok.sum()), "amplitude": np.nan, "doy_max": np.nan}
    t = years_since(d)[ok]
    w = 2 * np.pi * t
    ids = np.unique(seg[ok])
    cols = [np.cos(w), np.sin(w), t] + [(seg[ok] == i).astype(float) for i in ids]
    b, *_ = np.linalg.lstsq(np.column_stack(cols), v[ok], rcond=None)
    amp = float(np.hypot(b[0], b[1]))
    phase_years = (np.arctan2(b[1], b[0]) / (2 * np.pi)) % 1
    doy = float((d.iloc[0].dayofyear + phase_years * 365.25 - 1) % 365.25 + 1)
    return {"n": int(ok.sum()), "amplitude": amp, "doy_max": doy}


def partial_r(x, y, z) -> float:
    """Correlation of x and y after a linear regression on z is removed from both."""
    x, y, z = (np.asarray(a, float) for a in (x, y, z))
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    Z = np.column_stack([np.ones(ok.sum()), z[ok]])
    rx = x[ok] - Z @ np.linalg.lstsq(Z, x[ok], rcond=None)[0]
    ry = y[ok] - Z @ np.linalg.lstsq(Z, y[ok], rcond=None)[0]
    return float(np.corrcoef(rx, ry)[0, 1])

"""Is the peatland drying? The three readings of X-138: a season total from a radar chain, the open-water area from a
stack of water flags, and a trend sentence fixed in advance.

Each routine takes plain arrays and returns plain numbers, so a synthetic answer can be put in and read back
(``tests/test_drying.py``). Nothing here reads the cube or the hub.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats


def season_run(c: pd.DataFrame, year: int, months: tuple[int, int]) -> pd.DataFrame:
    """The longest unbroken chain segment of ``c`` (time, seg, est) inside ``months`` of ``year``; ties: the earlier.
    The run is re-zeroed at its first epoch."""
    t = pd.to_datetime(c.time, utc=True)
    s = c[(t.dt.year == year) & t.dt.month.between(*months)].copy()
    if s.empty:
        return s
    n = s.groupby("seg").size()
    seg = n[n == n.max()].index.min()
    run = s[s.seg == seg].sort_values("time").copy()
    run["est"] -= run.est.iloc[0]
    return run.reset_index(drop=True)


def season_total(run: pd.DataFrame, cos_inc: float) -> dict:
    """Net change (last − first) and largest drop (minimum − first) of a run, line of sight (mm) → vertical (cm)."""
    z = run.est.to_numpy(float) / cos_inc / 10.0
    i = int(np.argmin(z))
    days = np.diff(pd.to_datetime(run.time, utc=True)).astype("timedelta64[h]").astype(float) / 24.0
    return {"start": pd.Timestamp(run.time.iloc[0]).date(), "end": pd.Timestamp(run.time.iloc[-1]).date(),
            "overpasses": len(run), "net_cm": float(z[-1]), "drop_cm": float(z[i]),
            "drop_date": pd.Timestamp(run.time.iloc[i]).date(),
            "median_revisit_days": float(np.median(days)) if len(days) else float("nan")}


def revisit_class(median_days: float) -> str:
    return "six-day" if median_days <= 9 else "twelve-day"


def error_band(diff_cm) -> dict:
    """Typical (root-mean-square) and worst seen (max |·|) of radar − laser differences."""
    d = np.asarray(diff_cm, float)
    d = d[np.isfinite(d)]
    return {"n": len(d), "typical_cm": float(np.sqrt(np.mean(d ** 2))), "worst_cm": float(np.max(np.abs(d)))}


def water_area(water: np.ndarray, valid: np.ndarray, region: np.ndarray, threshold: float = 0.5,
               pixel_ha: float = 0.04) -> tuple[float, np.ndarray]:
    """Open-water area of ``region`` from per-scene flags (scenes × rows × cols): a pixel is open water when it is
    water in at least ``threshold`` of the scenes in which it is valid. Returns hectares and the frequency map."""
    w = np.asarray(water, bool) & np.asarray(valid, bool)
    nv = np.asarray(valid, bool).sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        freq = np.where(nv > 0, w.sum(axis=0) / nv, np.nan)
    open_ = region & np.isfinite(freq) & (freq >= threshold)
    return float(open_.sum() * pixel_ha), freq


def trend(years, values) -> dict:
    """Theil–Sen slope (per year) with Kendall's tau and its two-sided p."""
    x, y = np.asarray(years, float), np.asarray(values, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    slope, intercept, lo, hi = stats.theilslopes(y, x)
    tau, p = stats.kendalltau(x, y)
    return {"n": len(x), "slope": float(slope), "slope_lo": float(lo), "slope_hi": float(hi), "intercept": float(intercept),
            "tau": float(tau), "p": float(p), "span_years": float(x.max() - x.min()) if len(x) else float("nan")}


def trend_sentence(tr: dict, threshold: float, falling: str, rising: str, none: str, alpha: float = 0.05) -> str:
    """The pre-registered rule: a direction is stated only if p < alpha and |slope| × span exceeds ``threshold``."""
    big = abs(tr["slope"]) * tr["span_years"] > threshold
    if tr["p"] < alpha and big:
        return falling if tr["slope"] < 0 else rising
    return none


def synthesis(a_drying: bool | None, b_drying: bool | None) -> str:
    """X-138's synthesis from lines A and B: True = a drying sentence, False = the opposite, None = no trend."""
    if a_drying and b_drying:
        return "the data show the peatland drying"
    if (a_drying is True and b_drying is False) or (a_drying is False and b_drying is True):
        return "the lines disagree"
    if a_drying or b_drying:
        return "one line shows a sign of drying"
    return "the data do not show a drying trend over 2017–2026"

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


# --------------------------------------------------------------------------------------------- X-139: water balance
GSC = 0.0820                                                         # solar constant, MJ m⁻² min⁻¹ (FAO-56)
KACZOROWSKA = ((50, "extremely dry"), (75, "very dry"), (90, "dry"), (110.0001, "normal"), (125.0001, "wet"),
               (150.0001, "very wet"), (np.inf, "extremely wet"))


def extraterrestrial_radiation(doy, lat_deg: float) -> np.ndarray:
    """Daily extraterrestrial radiation Ra (MJ m⁻² day⁻¹), FAO-56 eqs. 21–25."""
    j = np.asarray(doy, float)
    phi = np.radians(lat_deg)
    dr = 1 + 0.033 * np.cos(2 * np.pi * j / 365)
    dec = 0.409 * np.sin(2 * np.pi * j / 365 - 1.39)
    ws = np.arccos(np.clip(-np.tan(phi) * np.tan(dec), -1, 1))
    return 24 * 60 / np.pi * GSC * dr * (ws * np.sin(phi) * np.sin(dec) + np.cos(phi) * np.cos(dec) * np.sin(ws))


def hargreaves_et0(tmax, tmin, doy, lat_deg: float) -> np.ndarray:
    """Reference evapotranspiration (mm day⁻¹), Hargreaves–Samani, FAO-56 eq. 52 (Ra converted by 0.408)."""
    tmax, tmin = np.asarray(tmax, float), np.asarray(tmin, float)
    rng = np.clip(tmax - tmin, 0, None)
    return 0.0023 * ((tmax + tmin) / 2 + 17.8) * np.sqrt(rng) * 0.408 * extraterrestrial_radiation(doy, lat_deg)


def station_mean(wide: pd.DataFrame) -> pd.Series:
    """Mean over the station columns of each day, from the stations present (NaN only when all are missing)."""
    return wide.mean(axis=1, skipna=True)


def rain_class(pct_of_normal: float) -> str:
    """Relative precipitation index classes of Kaczorowska (1962), % of the normal."""
    if not np.isfinite(pct_of_normal):
        return ""
    for edge, name in KACZOROWSKA:
        if pct_of_normal < edge:
            return name
    return KACZOROWSKA[-1][1]


# --------------------------------------------------------------------------------------------- X-141: winter model
def lowo_linear(x, y) -> np.ndarray:
    """Leave-one-out predictions of y from a + b·x (least squares on the others)."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    out = np.full(len(y), np.nan)
    for i in range(len(y)):
        m = np.arange(len(y)) != i
        if m.sum() >= 2 and np.ptp(x[m]) > 0:
            b, a = np.polyfit(x[m], y[m], 1)
            out[i] = a + b * x[i]
    return out


def lowo_climatology(y) -> np.ndarray:
    """Leave-one-out predictions of y by the mean of the others."""
    y = np.asarray(y, float)
    return np.array([np.mean(np.delete(y, i)) for i in range(len(y))])


def nested_selection(X: dict, y) -> tuple[np.ndarray, list[str]]:
    """Outer leave-one-out; in each fold the predictor with the lowest inner leave-one-out error on the training
    points is chosen, refitted on them and used for the left-out point. Returns predictions and the choices."""
    y = np.asarray(y, float)
    pred, chosen = np.full(len(y), np.nan), []
    for i in range(len(y)):
        m = np.arange(len(y)) != i
        best, best_e = None, np.inf
        for name, x in X.items():
            x = np.asarray(x, float)
            e = lowo_linear(x[m], y[m]) - y[m]
            if np.all(np.isfinite(e)) and np.sqrt(np.mean(e ** 2)) < best_e:
                best, best_e = name, float(np.sqrt(np.mean(e ** 2)))
        chosen.append(best or "")
        if best is not None:
            x = np.asarray(X[best], float)
            b, a = np.polyfit(x[m], y[m], 1)
            pred[i] = a + b * x[i]
    return pred, chosen

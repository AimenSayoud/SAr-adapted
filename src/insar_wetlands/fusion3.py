"""Fusion v3 building blocks (X-073…, branch fusion-x062): references, chains, the per-pair noise model, coverage.

Every routine is tested on synthetic ground truth (``tests/test_fusion3.py``). Phase is read from the data cube in
radians; conversion to millimetres goes through ``cube.core.rad_to_los_mm``.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from .cube.core import WAVELENGTH_MM

# ------------------------------------------------------------------------------ references

def reference_values(U: np.ndarray, mask: np.ndarray, stat: str = "median") -> np.ndarray:
    """Per-pair reference of a (pair, y, x) stack over the True cells of ``mask`` (NaN-robust median or mean)."""
    v = U[:, mask]
    with np.errstate(all="ignore"), __import__("warnings").catch_warnings():
        __import__("warnings").simplefilter("ignore", RuntimeWarning)
        return np.nanmedian(v, 1) if stat == "median" else np.nanmean(v, 1)


def nearest_cells(mask: np.ndarray, row: int, col: int, k: int) -> np.ndarray:
    """The ``k`` True cells of ``mask`` nearest to (row, col), as a mask."""
    rr, cc = np.nonzero(mask)
    order = np.argsort((rr - row) ** 2 + (cc - col) ** 2)[:k]
    out = np.zeros_like(mask, bool)
    out[rr[order], cc[order]] = True
    return out


# ------------------------------------------------------------------------------ fits and chains

def slope_origin(x: np.ndarray, y: np.ndarray) -> float:
    ok = np.isfinite(x) & np.isfinite(y)
    return float(np.sum(x[ok] * y[ok]) / np.sum(x[ok] ** 2)) if ok.sum() > 1 else np.nan


def year_blocked_residual(x: np.ndarray, y: np.ndarray, years: np.ndarray, min_train: int = 6) -> np.ndarray:
    """y − m·x with m (through the origin) fitted on the other years only."""
    res = np.full(len(y), np.nan)
    ok = np.isfinite(x) & np.isfinite(y)
    for yv in np.unique(years[ok]):
        tr, te = ok & (years != yv), ok & (years == yv)
        if tr.sum() >= min_train and te.any():
            res[te] = y[te] - slope_origin(x[tr], y[tr]) * x[te]
    return res


def chain_runs(t1: pd.DatetimeIndex, t2: pd.DatetimeIndex, *series: np.ndarray, min_len: int = 4):
    """Runs of consecutive pairs (t2 of one = t1 of the next, every value finite); returns a run id per pair (−1 =
    not in a run of at least ``min_len``) and each series' cumulative level within its run, centred on the run mean."""
    n = len(t1)
    order = np.argsort(np.asarray(t1))
    run = np.full(n, -1)
    levels = [np.full(n, np.nan) for _ in series]
    cur, rid = [], 0

    def close(cur, rid):
        if len(cur) >= min_len:
            for s, lv in zip(series, levels):
                c = np.cumsum(s[cur])
                lv[cur] = c - c.mean()
            run[cur] = rid
            return rid + 1
        return rid

    for k in order:
        ok = all(np.isfinite(s[k]) for s in series)
        if ok and cur and t1[k] == t2[cur[-1]]:
            cur.append(k)
        else:
            rid = close(cur, rid)
            cur = [k] if ok else []
    close(cur, rid)
    return run, levels


def level_metrics(est: np.ndarray, truth: np.ndarray) -> dict:
    ok = np.isfinite(est) & np.isfinite(truth)
    if ok.sum() < 5:
        return {"n_levels": int(ok.sum()), "r_level": np.nan, "rmse_level_mm": np.nan, "amplitude_ratio": np.nan}
    e, t = est[ok], truth[ok]
    return {"n_levels": int(ok.sum()), "r_level": float(np.corrcoef(e, t)[0, 1]),
            "rmse_level_mm": float(np.sqrt(np.mean((e - t) ** 2))), "amplitude_ratio": float(np.std(e) / np.std(t))}


# ------------------------------------------------------------------------------ noise model

def crb_sigma_mm(coh, looks: float, floor_mm: float, wavelength_mm: float = WAVELENGTH_MM) -> np.ndarray:
    """Phase noise in mm of line of sight: the Cramér–Rao bound for ``looks`` independent looks at coherence γ,
    √((1−γ²)/(2Lγ²)) rad, converted to mm, added in quadrature to a floor (atmosphere, unmodelled surface)."""
    g = np.clip(np.asarray(coh, float), 0.02, 0.999)
    rad = np.sqrt((1 - g ** 2) / (2 * looks * g ** 2))
    return np.sqrt(floor_mm ** 2 + (rad * wavelength_mm / (4 * np.pi)) ** 2)


def robust_sd(v: np.ndarray) -> float:
    v = v[np.isfinite(v)]
    return float(1.4826 * np.median(np.abs(v - np.median(v)))) if len(v) > 10 else np.nan


def fit_noise(coh_bin_centres: np.ndarray, sd_mm: np.ndarray, weights: np.ndarray | None = None) -> dict:
    """Fit (looks, floor) of ``crb_sigma_mm`` to binned robust SDs; weighted by the bin counts."""
    ok = np.isfinite(coh_bin_centres) & np.isfinite(sd_mm)
    c, s = coh_bin_centres[ok], sd_mm[ok]
    w = np.ones_like(s) if weights is None else np.sqrt(np.asarray(weights, float)[ok])

    def r(p):
        return w * (crb_sigma_mm(c, np.exp(p[0]), abs(p[1])) - s)

    sol = least_squares(r, x0=[np.log(10.0), 2.0])
    return {"looks": float(np.exp(sol.x[0])), "floor_mm": float(abs(sol.x[1])), "cost": float(sol.cost)}


def coverage(z: np.ndarray) -> dict:
    """Share of standardised residuals inside ±1 and ±2 (expected 0.68 and 0.95 if the noise model is right)."""
    z = z[np.isfinite(z)]
    return {"n": int(len(z)), "within_1sd": float(np.mean(np.abs(z) <= 1)) if len(z) else np.nan,
            "within_2sd": float(np.mean(np.abs(z) <= 2)) if len(z) else np.nan,
            "sd_z": float(np.std(z)) if len(z) else np.nan}


# ------------------------------------------------------------------------------ the v3 smoother

def ou_increment_smoother(t_days: np.ndarray, i: np.ndarray, j: np.ndarray, y: np.ndarray, sigma: np.ndarray,
                          tau_days: float, sigma_r: float) -> tuple[np.ndarray, np.ndarray]:
    """Exact posterior of r(t) at the nodes ``t_days`` (increasing) given increment observations
    y_k = r[j_k] − r[i_k] + ε_k, ε_k ~ N(0, σ_k²), and a stationary OU prior on r (τ, σ_r). Returns the posterior
    mean and SD of r. This is the Kalman/RTS smoother of the model, solved in one sparse system (no network
    inversion: each node is tied only by its own increments and the prior)."""
    import scipy.sparse as sp
    import scipy.sparse.linalg as spla

    from .fusion2 import ou_precision
    n = len(t_days)
    ok = np.isfinite(y) & np.isfinite(sigma) & (sigma > 0)
    i, j, y, s = i[ok], j[ok], y[ok], sigma[ok]
    m = len(y)
    D = sp.csr_matrix((np.r_[-np.ones(m), np.ones(m)], (np.r_[np.arange(m), np.arange(m)], np.r_[i, j])), shape=(m, n))
    Wt = sp.diags(1 / s ** 2)
    P = (ou_precision(t_days, tau_days, sigma_r) + D.T @ Wt @ D).tocsc()
    mean = spla.spsolve(P, D.T @ (y / s ** 2))
    cov = np.linalg.inv(P.toarray()) if n <= 3000 else None
    sd = np.sqrt(np.diag(cov)) if cov is not None else np.full(n, np.nan)
    return np.asarray(mean), sd


def detect_inflated_sigma(sigma: np.ndarray, p_detect: np.ndarray, floor: float = 0.05,
                          miss_scale_mm: float = WAVELENGTH_MM / 2) -> np.ndarray:
    """Observation SD for the smoother: the noise-model σ when the pair surely carries the motion, growing toward the
    error scale of a miss as the detectability falls (a mixture's variance, σ² + (1 − p)·scale²). The scale is a
    whole cycle (λ/2) where misses are cycle slips; λ/4 where they are not (6-day pairs, WP4)."""
    p = np.clip(np.asarray(p_detect, float), floor, 1.0)
    return np.sqrt(np.asarray(sigma, float) ** 2 + (1 - p) * miss_scale_mm ** 2)

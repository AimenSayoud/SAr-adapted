"""SPAMS — Simple Parameterization for the Motion of Soils (Conroy, van Diepen & Hanssen 2023, Geoderma 440, 116699).

The relative surface height of a shallow soft soil is a reversible part plus an irreversible part:

    R(t) = Σ_{k=0}^{τ-1} [x_P·P(t−k) − x_E·E(t−k)]          (eq. 2: water surplus over a memory time τ)
    I(t) = Σ_{s≤t} x_I · 1[R(s) ≤ 0]                         (eqs. 3–4: a constant loss, only while drying)

P and E are daily precipitation and evapotranspiration (mm). For a fixed τ the change of R over a lag is linear in
(x_P, x_E), so the fit is least squares on observed changes, with τ chosen on a grid by the training error. The
irreversible switch depends on the sign of R, which needs (x_P, x_E) first: it is fitted after them, as the paper's
model is. Used by X-112 as a published baseline for the floating mat's prior.
"""
from __future__ import annotations

import numpy as np


def window_sum(x: np.ndarray, tau: int) -> np.ndarray:
    """Σ of the last ``tau`` values (including today); NaN until ``tau`` values exist. Missing inputs count as 0."""
    x = np.nan_to_num(np.asarray(x, float))
    c = np.concatenate([[0.0], np.cumsum(x)])
    out = np.full(len(x), np.nan)
    out[tau - 1:] = c[tau:] - c[:-tau]
    return out


def reversible(P, E, x_p: float, x_e: float, tau: int) -> np.ndarray:
    """R(t), eq. 2."""
    return x_p * window_sum(P, tau) - x_e * window_sum(E, tau)


def irreversible(R: np.ndarray, x_i: float) -> np.ndarray:
    """I(t), eqs. 3–4: x_I accumulated over the days the reversible part is not positive (drying)."""
    dry = np.where(np.isfinite(R), (R <= 0).astype(float), 0.0)
    return x_i * np.cumsum(dry)


def _lagged(a: np.ndarray, lag: int) -> np.ndarray:
    out = np.full(len(a), np.nan)
    out[:-lag] = a[lag:] - a[:-lag]
    return out


def design(P, E, tau: int, lag: int) -> np.ndarray:
    """Columns whose combination x_P·col0 − x_E·col1 is the change of R over ``lag`` days, at each start day."""
    return np.column_stack([_lagged(window_sum(P, tau), lag), _lagged(window_sum(E, tau), lag)])


def fit(P, E, dh, train: np.ndarray, lag: int, taus=range(10, 151, 5)) -> dict:
    """Least-squares (x_P, x_E) for each τ on the training days' observed changes ``dh`` (same units as the output);
    τ = the one with the lowest training error. ``dh[t]`` is the change from day t to day t + lag."""
    dh = np.asarray(dh, float)
    best = None
    for tau in taus:
        A = design(P, E, tau, lag)
        A = np.column_stack([A[:, 0], -A[:, 1]])
        ok = train & np.isfinite(dh) & np.all(np.isfinite(A), 1)
        if ok.sum() < 10:
            continue
        b, *_ = np.linalg.lstsq(A[ok], dh[ok], rcond=None)
        sse = float(np.sum((A[ok] @ b - dh[ok]) ** 2))
        if best is None or sse < best["sse"]:
            best = {"x_p": float(b[0]), "x_e": float(b[1]), "tau": int(tau), "sse": sse, "n": int(ok.sum())}
    if best is None:
        raise ValueError("no τ had enough training days")
    return best


def predict_change(P, E, par: dict, lag: int) -> np.ndarray:
    """The change of R over ``lag`` days predicted by fitted parameters."""
    A = design(P, E, par["tau"], lag)
    return par["x_p"] * A[:, 0] - par["x_e"] * A[:, 1]


def fit_irreversible(P, E, par: dict, dh, train: np.ndarray, lag: int) -> float:
    """x_I, by least squares on what the reversible change leaves, given (x_P, x_E, τ)."""
    R = reversible(P, E, par["x_p"], par["x_e"], par["tau"])
    di = _lagged(irreversible(R, 1.0), lag)
    res = np.asarray(dh, float) - predict_change(P, E, par, lag)
    ok = train & np.isfinite(res) & np.isfinite(di) & (di != 0)
    return float(np.sum(di[ok] * res[ok]) / np.sum(di[ok] ** 2)) if ok.sum() >= 10 else 0.0

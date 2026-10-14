"""Forward models of a floating mat's surface level from its water table (fusion v4, X-117).

Memory response: the surface follows an exponentially smoothed water table,

    W_f(t) = W_f(t−1) + (W(t) − W_f(t−1)) / τ,     Δh(t0→t1) = g · (W_f(t1) − W_f(t0)) + c · (t1 − t0)

with g the buoyancy coefficient, τ the memory (days; τ = 1 is no memory) and c a constant drift rate of the surface
record (the laser's slow drift). For a fixed τ the model is linear in (g, c): least squares on observed changes over
intervals of any length, τ chosen on a grid by the training error.
"""
from __future__ import annotations

import numpy as np


def smooth(w, tau: float) -> np.ndarray:
    """Exponential memory of a daily series; gaps carried by the last value (NaN until the first value)."""
    w = np.asarray(w, float)
    out = np.full(len(w), np.nan)
    prev = np.nan
    for i, x in enumerate(w):
        if np.isfinite(x):
            prev = x if not np.isfinite(prev) else prev + (x - prev) / tau
        out[i] = prev
    return out


def fit(w, i0: np.ndarray, i1: np.ndarray, dh: np.ndarray, taus=(1, 2, 3, 5, 7, 10, 14, 20, 30, 45, 60),
        drift: bool = True, weights: np.ndarray | None = None) -> dict:
    """(g, c, τ) from observed changes ``dh`` between day indices ``i0`` → ``i1``; τ by the (weighted) training error."""
    i0, i1, dh = np.asarray(i0), np.asarray(i1), np.asarray(dh, float)
    wt = np.ones(len(dh)) if weights is None else np.asarray(weights, float)
    best = None
    for tau in taus:
        wf = smooth(w, tau)
        x = wf[i1] - wf[i0]
        A = np.column_stack([x, (i1 - i0).astype(float)]) if drift else x[:, None]
        ok = np.isfinite(dh) & np.all(np.isfinite(A), 1)
        if ok.sum() < 10:
            continue
        sw = np.sqrt(wt[ok])
        b, *_ = np.linalg.lstsq(A[ok] * sw[:, None], dh[ok] * sw, rcond=None)
        sse = float(np.sum(wt[ok] * (A[ok] @ b - dh[ok]) ** 2))
        if best is None or sse < best["sse"]:
            best = {"tau": tau, "g": float(b[0]), "c": float(b[1]) if drift else 0.0, "sse": sse, "n": int(ok.sum())}
    if best is None:
        raise ValueError("no τ had enough intervals")
    return best


def predict(w, par: dict, i0, i1) -> np.ndarray:
    wf = smooth(w, par["tau"])
    i0, i1 = np.asarray(i0), np.asarray(i1)
    return par["g"] * (wf[i1] - wf[i0]) + par["c"] * (i1 - i0)

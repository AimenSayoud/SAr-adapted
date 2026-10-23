"""Matched contrasts between groups of cells, resampled by spatial block (X-132).

* ``match_nearest`` — each treated unit, in a fixed order, takes the unused control closest on one covariate within a
  caliper, ties broken by a second covariate; without replacement.
* ``block_counts`` + ``weighted_median`` — a block bootstrap as weights: each draw resamples the spatial blocks with
  replacement and every unit counts as many times as its block was drawn, so several statistics computed from one draw
  share it (a difference of two contrasts keeps their covariance).
* ``block_sums`` + ``ols_from_sums`` + ``ols_block_bootstrap`` — least squares from per-block sufficient statistics
  (XᵀX, Xᵀy, yᵀy, Σy, n): the full fit, any column left out and every bootstrap draw are sums of the same block arrays,
  so a regression on a million rows resamples in seconds.

Pure numpy / pandas; no field values here.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


# ------------------------------------------------------------------------------ matching

def match_nearest(treated: pd.DataFrame, control: pd.DataFrame, on: str, caliper: float, tie: str,
                  order: str) -> pd.DataFrame:
    """Pairs (treated index, control index): treated units taken by increasing ``order``; each takes the unused control
    with the smallest |Δ ``on``| within ``caliper``, ties (equal |Δ ``on``|) broken by the smaller |Δ ``tie``|, then by the
    control's index. Unmatched treated units are absent from the result."""
    used = np.zeros(len(control), bool)
    cv, ct = control[on].to_numpy(float), control[tie].to_numpy(float)
    cidx = control.index.to_numpy()
    rows = []
    for i, r in treated.sort_values([order], kind="stable").iterrows():
        d = np.abs(cv - r[on])
        ok = ~used & (d <= caliper)
        if not ok.any():
            continue
        cand = np.flatnonzero(ok)
        key = np.lexsort((cidx[cand], np.abs(ct[cand] - r[tie]), d[cand]))
        j = cand[key[0]]
        used[j] = True
        rows.append({"treated": i, "control": cidx[j], "d_on": float(d[j]), "d_tie": float(abs(ct[j] - r[tie]))})
    return pd.DataFrame(rows, columns=["treated", "control", "d_on", "d_tie"])


# ------------------------------------------------------------------------------ block bootstrap as weights

def block_counts(blocks, n_boot: int, seed: int = 0) -> tuple[np.ndarray, np.ndarray]:
    """(the distinct blocks, an (n_boot, n_blocks) array of how often each was drawn): blocks resampled with
    replacement, as many draws per resample as there are distinct blocks."""
    ub = np.unique(np.asarray(blocks))
    rng = np.random.default_rng(seed)
    draws = rng.integers(0, len(ub), (n_boot, len(ub)))
    counts = np.zeros((n_boot, len(ub)), int)
    np.add.at(counts, (np.repeat(np.arange(n_boot), len(ub)), draws.ravel()), 1)
    return ub, counts


def unit_weights(unit_blocks, ub: np.ndarray, counts: np.ndarray) -> np.ndarray:
    """(n_boot, n_units) weights: each unit carries its block's count in every draw."""
    pos = np.searchsorted(ub, np.asarray(unit_blocks))
    if not np.array_equal(ub[np.clip(pos, 0, len(ub) - 1)], np.asarray(unit_blocks)):
        raise ValueError("a unit's block is not among the resampled blocks")
    return counts[:, pos]


def weighted_median(x: np.ndarray, w: np.ndarray) -> float:
    """Median of ``x`` with non-negative integer-like weights ``w`` (the median of x repeated w times; for an even total,
    the mean of the two middle values). NaN when the weights sum to zero."""
    x, w = np.asarray(x, float), np.asarray(w, float)
    ok = np.isfinite(x) & (w > 0)
    x, w = x[ok], w[ok]
    if w.sum() == 0:
        return float("nan")
    o = np.argsort(x, kind="stable")
    x, c = x[o], np.cumsum(w[o])
    tot = c[-1]
    lo = x[np.searchsorted(c, tot / 2.0, side="left")]
    hi = x[np.searchsorted(c, tot / 2.0, side="right")] if tot % 2 == 0 else lo
    return float((lo + hi) / 2.0)


def boot_interval(draws, level: float = 0.95) -> tuple[float, float]:
    d = np.asarray(draws, float)
    d = d[np.isfinite(d)]
    if not len(d):
        return float("nan"), float("nan")
    a = (1 - level) / 2 * 100
    lo, hi = np.percentile(d, [a, 100 - a])
    return float(lo), float(hi)


# ------------------------------------------------------------------------------ least squares from block sums

def block_sums(X: np.ndarray, y: np.ndarray, blocks) -> dict:
    """Per-block XᵀX (B, p, p), Xᵀy (B, p), yᵀy, Σy, n of a design ``X`` (rows × p, intercept included by the caller)."""
    X, y, b = np.asarray(X, float), np.asarray(y, float), np.asarray(blocks)
    ub, inv = np.unique(b, return_inverse=True)
    B, p = len(ub), X.shape[1]
    xtx, xty = np.zeros((B, p, p)), np.zeros((B, p))
    yty, sy, n = np.zeros(B), np.zeros(B), np.zeros(B)
    for k in range(B):
        m = inv == k
        Xk, yk = X[m], y[m]
        xtx[k], xty[k] = Xk.T @ Xk, Xk.T @ yk
        yty[k], sy[k], n[k] = yk @ yk, yk.sum(), m.sum()
    return {"blocks": ub, "xtx": xtx, "xty": xty, "yty": yty, "sy": sy, "n": n}


def ols_from_sums(xtx: np.ndarray, xty: np.ndarray, yty: float, sy: float, n: float,
                  cols=None) -> tuple[np.ndarray, float]:
    """Coefficients and R² (about the mean) of the least-squares fit on the columns ``cols`` (all when None)."""
    cols = np.arange(xtx.shape[0]) if cols is None else np.asarray(cols)
    A, b = xtx[np.ix_(cols, cols)], xty[cols]
    beta = np.linalg.lstsq(A, b, rcond=None)[0]
    sse = yty - 2 * beta @ b + beta @ A @ beta
    sst = yty - sy ** 2 / n
    return beta, float(1 - sse / sst) if sst > 0 else float("nan")


def ols_block_bootstrap(S: dict, groups: dict[str, list[int]], n_boot: int, seed: int = 0) -> dict:
    """The full fit and, per draw of blocks (with replacement), its coefficients and each group's R² drop (R² full −
    R² with the group's columns left out). ``groups`` maps a factor name to its column indices (never the intercept)."""
    p = S["xtx"].shape[1]

    def fit(w):
        xtx = np.tensordot(w, S["xtx"], 1)
        xty, yty, sy, n = w @ S["xty"], w @ S["yty"], w @ S["sy"], w @ S["n"]
        beta, r2 = ols_from_sums(xtx, xty, yty, sy, n)
        drops = {}
        for g, cs in groups.items():
            keep = [c for c in range(p) if c not in set(cs)]
            drops[g] = r2 - ols_from_sums(xtx, xty, yty, sy, n, keep)[1]
        return beta, r2, drops

    one = np.ones(len(S["blocks"]))
    beta, r2, drops = fit(one)
    _, counts = block_counts(S["blocks"], n_boot, seed)
    B, R, D = [], [], {g: [] for g in groups}
    for w in counts:
        b_, r_, d_ = fit(w.astype(float))
        B.append(b_)
        R.append(r_)
        for g in groups:
            D[g].append(d_[g])
    return {"beta": beta, "r2": r2, "drop": drops, "beta_draws": np.array(B), "r2_draws": np.array(R),
            "drop_draws": {g: np.array(v) for g, v in D.items()}}

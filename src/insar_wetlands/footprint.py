"""How far one truth point's motion reaches over many cells (X-133).

* ``pearson_cols`` / ``slope_origin_cols`` — a truth series (pairs) against every column of a cell stack (pairs × cells),
  each over the pairs where both are finite.
* ``permuted_r`` — the same r with the truth permuted within groups (calendar month), a null that keeps the season.
* ``longest_run`` — the longest chain of pairs where each pair starts on the date the previous one ends.
* ``resample_index`` — a block-bootstrap draw (counts per block) as row indices, for statistics that need the rows
  themselves (a rank correlation).

Pure numpy; no field values here.
"""
from __future__ import annotations

import numpy as np


def _complete(x: np.ndarray, Y: np.ndarray):
    x = np.asarray(x, float)
    Y = np.asarray(Y, float)
    m = np.isfinite(Y) & np.isfinite(x)[:, None]
    return np.where(m, x[:, None], 0.0), np.where(m, Y, 0.0), m


def pearson_cols(x, Y, min_n: int = 5) -> np.ndarray:
    """Pearson r of ``x`` with every column of ``Y``, pairwise complete; NaN with fewer than ``min_n`` pairs."""
    X, Yc, m = _complete(x, Y)
    n = m.sum(0)
    with np.errstate(invalid="ignore", divide="ignore"):
        mx, my = X.sum(0) / n, Yc.sum(0) / n
        dx, dy = np.where(m, X - mx, 0), np.where(m, Yc - my, 0)
        r = (dx * dy).sum(0) / np.sqrt((dx ** 2).sum(0) * (dy ** 2).sum(0))
    return np.where(n >= min_n, r, np.nan)


def slope_origin_cols(x, Y, min_n: int = 5) -> np.ndarray:
    """Least-squares slope through the origin of every column of ``Y`` on ``x`` (the gain against the truth)."""
    X, Yc, m = _complete(x, Y)
    n = m.sum(0)
    with np.errstate(invalid="ignore", divide="ignore"):
        s = (X * Yc).sum(0) / (X ** 2).sum(0)
    return np.where(n >= min_n, s, np.nan)


def permuted_r(x, Y, groups, n_perm: int, seed: int = 0) -> np.ndarray:
    """(n_perm, cells) Pearson r of every column of ``Y`` with ``x`` permuted within each group."""
    x, groups = np.asarray(x, float), np.asarray(groups)
    rng = np.random.default_rng(seed)
    idx = [np.flatnonzero(groups == g) for g in np.unique(groups)]
    out = np.empty((n_perm, np.asarray(Y).shape[1]))
    for k in range(n_perm):
        xp = x.copy()
        for ii in idx:
            xp[ii] = x[rng.permutation(ii)]
        out[k] = pearson_cols(xp, Y)
    return out


def longest_run(t1, t2) -> np.ndarray:
    """Indices (in time order) of the longest run of pairs where each pair's start equals the previous pair's end;
    the earliest run wins a tie."""
    t1, t2 = np.asarray(t1), np.asarray(t2)
    o = np.argsort(t1, kind="stable")
    best, cur = [], []
    for i in o:
        if cur and t1[i] == t2[cur[-1]]:
            cur.append(i)
        else:
            cur = [i]
        if len(cur) > len(best):
            best = list(cur)
    return np.array(best, int)


def resample_index(counts_row: np.ndarray, unit_pos: np.ndarray) -> np.ndarray:
    """Row indices of one block-bootstrap draw: each unit repeated as often as its block (``unit_pos`` = the unit's
    position among the blocks of ``counts_row``) was drawn."""
    w = np.asarray(counts_row)[np.asarray(unit_pos)]
    return np.repeat(np.arange(len(w)), w)

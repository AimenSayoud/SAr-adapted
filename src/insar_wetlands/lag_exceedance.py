"""How often a measured surface moves beyond a pair's ambiguity limit, by the pair's length.

Why. A single interferometric pair reads a change unambiguously only within a quarter wavelength along the line of
sight; beyond it the phase wraps and the pair gives the wrong size and sign. Whether a network can carry a surface's
motion therefore depends on how often the surface moves beyond that limit over the network's pair lengths — a
question about the surface, answerable from a ground record alone, like a structure function or a variogram of the
level. This module counts it: every pair of dates within a record (optionally barred across datum breaks or filled
stretches), its change, and the share beyond a threshold in lag bins, with a Wilson interval.

The Wilson interval treats pairs as independent. Pairs sharing dates overlap in time, so the interval is narrower
than the record's real uncertainty; callers say so where they report it.
"""
from __future__ import annotations

from typing import Callable

import numpy as np
import pandas as pd


def wilson(k: np.ndarray | int, n: np.ndarray | int, z: float = 1.96) -> tuple[np.ndarray, np.ndarray]:
    """Wilson score interval for k successes in n trials (NaN where n is 0)."""
    k, n = np.asarray(k, float), np.asarray(n, float)
    with np.errstate(invalid="ignore", divide="ignore"):
        p = k / n
        den = 1 + z ** 2 / n
        mid = (p + z ** 2 / (2 * n)) / den
        half = z * np.sqrt(p * (1 - p) / n + z ** 2 / (4 * n ** 2)) / den
    return mid - half, mid + half


def pair_changes(times, level, max_lag_days: float,
                 barred: Callable[[pd.Timestamp, pd.Timestamp], bool] | None = None) -> pd.DataFrame:
    """Every pair (t1 < t2) of the given dates no more than ``max_lag_days`` apart with a level at both ends, its lag
    in days and its change ``level(t2) − level(t1)``; ``barred(t1, t2)`` True drops the pair."""
    t = pd.DatetimeIndex(pd.to_datetime(times, utc=True))
    v = np.asarray(level, float)
    keep = np.isfinite(v)
    order = np.argsort(t[keep].asi8)
    t, v = t[keep][order], v[keep][order]
    rows = []
    for i in range(len(t)):
        for j in range(i + 1, len(t)):
            lag = (t[j] - t[i]).total_seconds() / 86400
            if lag > max_lag_days:
                break
            if barred is not None and barred(t[i], t[j]):
                continue
            rows.append((t[i], t[j], lag, v[j] - v[i]))
    return pd.DataFrame(rows, columns=["t1", "t2", "lag_days", "change"])


def exceedance(changes: pd.DataFrame, threshold: float, edges) -> pd.DataFrame:
    """Share of pairs whose |change| exceeds ``threshold`` in each lag bin [edges[i], edges[i+1]), with the bin's
    median lag, its pair count and a Wilson 95 % interval. Empty bins are left out."""
    edges = np.asarray(edges, float)
    b = np.digitize(changes["lag_days"].to_numpy(), edges) - 1
    out = []
    for i in range(len(edges) - 1):
        g = changes[b == i]
        if g.empty:
            continue
        n, k = len(g), int((g["change"].abs() > threshold).sum())
        lo, hi = wilson(k, n)
        out.append({"lag_lo_days": edges[i], "lag_hi_days": edges[i + 1], "lag_median_days": float(g.lag_days.median()),
                    "pairs": n, "beyond": k, "share": k / n, "share_lo95": float(lo), "share_hi95": float(hi)})
    return pd.DataFrame(out)

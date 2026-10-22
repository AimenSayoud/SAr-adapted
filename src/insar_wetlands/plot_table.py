"""Plot-level Sentinel-1 statistics for the supervisor's table (C-071): Date | Plot | Phase/Displacement | Coherence |
VV | VH | RVI | Valid pixels, with the spread inside each extraction window.

The cube's gold tables keep one number per window (a median); the supervisor asked for the number of valid pixels, the
median or mean and the spread (SD or IQR) of every variable. These routines compute them from the silver stacks, for
any window, so the table can say how many pixels stand behind each value and how much they disagree.

Conventions:

* A window is the (2h+1)² block of cells around the plot's cell, clipped at the grid edge and restricted to a mask
  (the mat, zone A, by default in the caller). A pixel is valid where its value is finite and, for phase, where it was
  unwrapped (``conncomp`` > 0).
* Displacement is the line-of-sight change −(φ − φ_ref)·λ/4π, positive toward the satellite (uplift), referenced to
  a scalar per pair (the median phase of a reference area of the same pair). The reference shifts the median and the
  mean, not the SD or the IQR.
* Consecutive pairs: on one track, the pair whose two dates have no other acquisition of the set between them.
"""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd

from .cube.core import WAVELENGTH_MM

STATS = ("n_valid", "median", "mean", "sd", "iqr")


def window_index(row: int, col: int, shape: tuple[int, int], h: int = 1,
                 mask: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Row and column indices of the cells in the (2h+1)² window at (row, col), inside the grid and the mask."""
    r = np.arange(max(row - h, 0), min(row + h + 1, shape[0]))
    c = np.arange(max(col - h, 0), min(col + h + 1, shape[1]))
    rr, cc = (a.ravel() for a in np.meshgrid(r, c, indexing="ij"))
    if mask is not None:
        keep = np.asarray(mask, bool)[rr, cc]
        rr, cc = rr[keep], cc[keep]
    return rr, cc


def window_stats_stack(values: np.ndarray, valid: np.ndarray | None = None) -> pd.DataFrame:
    """Per time step of a (T, n) array of window pixels: n_valid, median, mean, SD (ddof 1), IQR — the stack version
    of ``field.window_stats`` (same definitions, NaN = invalid)."""
    v = np.asarray(values, float)
    if v.ndim == 1:
        v = v[None, :]
    if valid is not None:
        v = np.where(np.asarray(valid, bool), v, np.nan)
    n = np.isfinite(v).sum(1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        q1, med, q3 = np.nanpercentile(v, [25, 50, 75], axis=1)
        out = pd.DataFrame({"n_valid": n, "median": med, "mean": np.nanmean(v, 1),
                            "sd": np.where(n > 1, np.nanstd(v, 1, ddof=1) if v.shape[1] > 1 else np.nan, np.nan),
                            "iqr": q3 - q1})
    out.loc[n == 0, ["median", "mean", "sd", "iqr"]] = np.nan
    return out


def phase_to_los_mm(stats: pd.DataFrame, ref_rad, wavelength_mm: float = WAVELENGTH_MM) -> pd.DataFrame:
    """Window statistics of unwrapped phase (rad) → of line-of-sight displacement (mm, + toward the satellite),
    referenced to ``ref_rad`` (one value per row). The sign flip turns the quartiles round; the IQR and SD keep size."""
    k = wavelength_mm / (4 * np.pi)
    ref = np.asarray(ref_rad, float)
    return pd.DataFrame({"n_valid": stats["n_valid"].to_numpy(),
                         "median": -(stats["median"].to_numpy() - ref) * k,
                         "mean": -(stats["mean"].to_numpy() - ref) * k,
                         "sd": stats["sd"].to_numpy() * k, "iqr": stats["iqr"].to_numpy() * k})


def consecutive_mask(t1, t2) -> np.ndarray:
    """True for the pairs of one track whose dates are adjacent among all the dates the pairs use."""
    a = pd.DatetimeIndex(pd.to_datetime(t1, utc=True)).normalize().tz_convert(None)
    b = pd.DatetimeIndex(pd.to_datetime(t2, utc=True)).normalize().tz_convert(None)
    dates = a.union(b).unique().sort_values()
    return dates.get_indexer(b) - dates.get_indexer(a) == 1

"""Phase simulated from a ground truth, and the processing choices tested on it (X-095; conroy2022's benchmark design).

A truth increment per interferometric pair (mm of line of sight) is turned into what a radar would deliver: noise,
optional extra terms (a moisture phase, a constant short-pair bias), then wrapped to one cycle. The processing choices
then try to recover the truth: which cycle each increment gets (``choose_cycle``), and how the chain is joined across
gaps it should not be integrated through (``tie_segments``). No field values live here; the truth is passed in.

Sign and units: increments in mm LOS; one phase cycle is ``wavelength_mm / 2`` of LOS motion.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

WAVELENGTH_MM = 55.465763


def cycle_mm(wavelength_mm: float = WAVELENGTH_MM) -> float:
    return wavelength_mm / 2.0


def wrap_mm(x, wavelength_mm: float = WAVELENGTH_MM) -> np.ndarray:
    """Wrap LOS increments to (-c/2, c/2], c = one cycle (lambda/2)."""
    c = cycle_mm(wavelength_mm)
    x = np.asarray(x, float)
    w = x - c * np.round(x / c)
    return np.where(w <= -c / 2, w + c, w)


def simulate(truth_inc, sigma=0.0, extra=0.0, rng: np.random.Generator | None = None,
             wavelength_mm: float = WAVELENGTH_MM) -> tuple[np.ndarray, np.ndarray]:
    """(wrapped, unwrapped) increments: truth + extra + N(0, sigma), then wrapped. ``sigma``/``extra`` broadcast."""
    truth = np.asarray(truth_inc, float)
    noise = 0.0 if rng is None else rng.normal(0.0, 1.0, truth.shape) * np.asarray(sigma, float)
    full = truth + np.asarray(extra, float) + noise
    return wrap_mm(full, wavelength_mm), full


def choose_cycle(wrapped, reference=None, wavelength_mm: float = WAVELENGTH_MM) -> np.ndarray:
    """Unwrapped increments: the cycle that puts each one closest to ``reference`` (cycle 0 where it is NaN/None)."""
    w = np.asarray(wrapped, float)
    if reference is None:
        return w.copy()
    ref = np.asarray(reference, float)
    c = cycle_mm(wavelength_mm)
    k = np.where(np.isfinite(ref), np.round((ref - w) / c), 0.0)
    return w + k * c


def tie_segments(t1, t2, inc, keep, prior_level: pd.Series | None) -> np.ndarray:
    """Cut the chain where ``keep`` is False and tie each segment's level to ``prior_level`` (least squares).

    ``t1``, ``t2``: pair dates of consecutive pairs (sorted); ``inc``: increments; ``keep``: pairs used inside
    segments; ``prior_level``: a level indexed by epoch (NaN where unknown). Each segment is the run of kept pairs
    between cuts; its offset is the mean of (prior − chain level) over its epochs with a prior; a segment without any
    prior keeps the level where the previous segment ended (zero change across the cut). Returns new increments for
    every pair (inside a segment: unchanged; across a cut: the jump the tie implies), so that cumulating them gives the
    tied levels."""
    t1, t2 = pd.DatetimeIndex(t1), pd.DatetimeIndex(t2)
    inc = np.asarray(inc, float)
    keep = np.asarray(keep, bool)
    n = len(inc)
    if n == 0:
        return inc.copy()
    epochs = [t1[0]] + list(t2)
    seg = np.zeros(n + 1, int)
    lvl = np.zeros(n + 1)
    s = 0
    for i in range(n):
        if keep[i]:
            seg[i + 1] = s
            lvl[i + 1] = lvl[i] + inc[i]
        else:
            s += 1
            seg[i + 1] = s
            lvl[i + 1] = 0.0
    pri = (prior_level.reindex(pd.DatetimeIndex(epochs)).to_numpy(float) if prior_level is not None
           else np.full(n + 1, np.nan))
    out = lvl.copy()
    end = None
    for k in range(s + 1):
        idx = np.flatnonzero(seg == k)
        m = np.isfinite(pri[idx])
        if m.any():
            off = float(np.mean(pri[idx][m] - lvl[idx][m]))
        else:
            off = (end - lvl[idx[0]]) if end is not None else 0.0
        out[idx] = lvl[idx] + off
        end = out[idx[-1]]
    return np.diff(out)

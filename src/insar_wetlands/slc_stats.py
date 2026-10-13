"""Pixel statistics on a coregistered SLC stack (X-103, X-104): amplitude dispersion, statistically homogeneous
neighbours (SHP), boxcar and SHP-window coherence, and the coherence noise floor of an estimator.

Conventions: a stack is (time, row, col) complex; amplitudes are |s|. Everything is plain numpy/scipy (no Dolphin):
- PS candidates by amplitude dispersion D_A = std(A) / mean(A) (Ferretti et al. 2001);
- SHP by the two-sample Kolmogorov–Smirnov test on the amplitude time series of a pixel and each neighbour in a window
  (Ferretti et al. 2011, SqueeSAR), vectorised over pixels;
- the noise floor of a coherence estimate with L looks: the 95th percentile of |gamma_hat| for truly incoherent
  signals, simulated (an estimate with few looks is biased high at low coherence, Touzi et al. 1999).
"""
from __future__ import annotations

import numpy as np

KS_C_05 = 1.358         # Kolmogorov–Smirnov c(alpha) for alpha = 0.05


def amplitude_dispersion(stack: np.ndarray) -> np.ndarray:
    a = np.abs(stack).astype("float64")
    m = np.nanmean(a, axis=0)
    return np.where(m > 0, np.nanstd(a, axis=0) / m, np.nan)


def ks_statistic(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Two-sample KS statistic for each row of x against the same row of y; x, y: (n, T), equal T, no NaN."""
    n, t = x.shape
    both = np.concatenate([x, y], axis=1)
    order = np.argsort(both, axis=1, kind="mergesort")
    step = np.where(order < t, 1.0, -1.0) / t
    return np.max(np.abs(np.cumsum(step, axis=1)), axis=1)


def ks_critical(t: int, c: float = KS_C_05) -> float:
    return c * np.sqrt(2.0 / t)


def shp_count(amp: np.ndarray, rows: np.ndarray, cols: np.ndarray, half_r: int, half_c: int,
              crit: float | None = None) -> np.ndarray:
    """Number of statistically homogeneous neighbours of each target pixel (rows[i], cols[i]) within a
    (2·half_r+1) × (2·half_c+1) window, by the KS test on the amplitude time series. amp: (T, R, C) real."""
    t = amp.shape[0]
    crit = ks_critical(t) if crit is None else crit
    centre = amp[:, rows, cols].T
    count = np.zeros(len(rows), int)
    for dr in range(-half_r, half_r + 1):
        for dc in range(-half_c, half_c + 1):
            if dr == 0 and dc == 0:
                continue
            rr, cc = rows + dr, cols + dc
            ok = (rr >= 0) & (rr < amp.shape[1]) & (cc >= 0) & (cc < amp.shape[2])
            nb = np.full_like(centre, np.nan)
            nb[ok] = amp[:, rr[ok], cc[ok]].T
            good = ok & np.isfinite(nb).all(axis=1) & np.isfinite(centre).all(axis=1)
            if good.any():
                d = ks_statistic(centre[good], nb[good])
                count[np.flatnonzero(good)[d < crit]] += 1
    return count


def boxcar_coherence(s1: np.ndarray, s2: np.ndarray, size: tuple[int, int]) -> np.ndarray:
    """|<s1 s2*>| / sqrt(<|s1|^2> <|s2|^2>) over a boxcar of ``size`` (rows, cols) samples."""
    from scipy.ndimage import uniform_filter  # noqa: PLC0415
    s1 = np.nan_to_num(s1)
    s2 = np.nan_to_num(s2)
    p = s1 * np.conj(s2)
    num = np.hypot(uniform_filter(p.real, size), uniform_filter(p.imag, size))
    den = np.sqrt(uniform_filter(np.abs(s1) ** 2, size) * uniform_filter(np.abs(s2) ** 2, size))
    return np.where(den > 0, num / den, np.nan)


def noise_floor(looks: int, q: float = 0.95, n: int = 20000, seed: int = 0) -> float:
    """The q-quantile of the coherence estimate for two independent (incoherent) complex Gaussian signals."""
    rng = np.random.default_rng(seed)
    a = rng.normal(size=(n, looks)) + 1j * rng.normal(size=(n, looks))
    b = rng.normal(size=(n, looks)) + 1j * rng.normal(size=(n, looks))
    g = np.abs((a * np.conj(b)).sum(1)) / np.sqrt((np.abs(a) ** 2).sum(1) * (np.abs(b) ** 2).sum(1))
    return float(np.quantile(g, q))

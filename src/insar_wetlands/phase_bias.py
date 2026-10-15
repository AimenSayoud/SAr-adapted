"""Short-pair phase-bias correction of Maghsoudi et al. (2022, EarthArXiv), equations (2)–(13).

With N acquisitions and the 1-, 2- and 3-connection pairs (6-, 12-, 18-day for Sentinel-1A+B), the closure phases

    Δφ(i,i+2) = φ(i,i+2) − φ(i,i+1) − φ(i+1,i+2),       Δφ(i,i+3) = φ(i,i+3) − φ(i,i+1) − φ(i+1,i+2) − φ(i+2,i+3)

are, if the bias of a longer pair is a fixed multiple of the summed 1-connection biases it spans (δ(i,i+2) = a1·Σ,
δ(i,i+3) = a2·Σ), linear in the N−1 unknown 1-connection biases δ(i,i+1):

    Δφ(i,i+2) = (a1 − 1)(δ(i,i+1) + δ(i+1,i+2)),    Δφ(i,i+3) = (a2 − 1)(δ(i,i+1) + δ(i+1,i+2) + δ(i+2,i+3)).

Solved by least squares (2N−5 equations, overdetermined for N ≥ 5); each pair is then corrected by its bias. Closures
are taken wrapped (from complex products), so the model holds while they stay well inside ±π.
"""
from __future__ import annotations

import numpy as np


def closures(z1: np.ndarray, z2: np.ndarray, z3: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Wrapped 12–6 and 18–6 closure phases from complex 1-, 2-, 3-connection interferograms (arrays of length N−1,
    N−2, N−3, ``z1[i]`` = pair (i, i+1))."""
    c2 = np.angle(z2 * np.conj(z1[:-1]) * np.conj(z1[1:]))
    c3 = np.angle(z3 * np.conj(z1[:-2]) * np.conj(z1[1:-1]) * np.conj(z1[2:]))
    return c2, c3


def solve(c2: np.ndarray, c3: np.ndarray, a1: float = 0.47, a2: float = 0.31) -> np.ndarray:
    """The N−1 one-connection biases from the closures (NaN closures are left out)."""
    n1 = len(c2) + 1
    rows, obs = [], []
    for i, v in enumerate(c2):
        if np.isfinite(v):
            r = np.zeros(n1)
            r[i:i + 2] = a1 - 1
            rows.append(r), obs.append(v)
    for i, v in enumerate(c3):
        if np.isfinite(v):
            r = np.zeros(n1)
            r[i:i + 3] = a2 - 1
            rows.append(r), obs.append(v)
    A, b = np.array(rows), np.array(obs)
    d, *_ = np.linalg.lstsq(A, b, rcond=None)
    return d


def correct(phi1: np.ndarray, phi2: np.ndarray, phi3: np.ndarray, d: np.ndarray, a1: float = 0.47, a2: float = 0.31):
    """Corrected 1-, 2-, 3-connection phases (eqs. 11–13)."""
    s2 = d[:-1] + d[1:]
    s3 = d[:-2] + d[1:-1] + d[2:]
    return phi1 - d, phi2 - a1 * s2, phi3 - a2 * s3

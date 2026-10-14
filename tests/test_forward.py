"""The memory-response forward model recovers a known coefficient, memory and drift from a synthetic mat."""
import numpy as np

from insar_wetlands import forward


def _mat(n=1200, tau=14, g=4.0, c=0.05, seed=0):
    rng = np.random.default_rng(seed)
    w = np.cumsum(rng.normal(0, 0.6, n)) + 3 * np.sin(np.arange(n) * 2 * np.pi / 365)
    h = g * forward.smooth(w, tau) + c * np.arange(n)
    return w, h


def test_smooth_no_memory_is_identity():
    w = np.array([1.0, np.nan, 3.0, 4.0])
    assert np.allclose(forward.smooth(w, 1), [1, 1, 3, 4])


def test_fit_recovers_memory_drift_and_g():
    w, h = _mat()
    rng = np.random.default_rng(1)
    i0 = np.concatenate([np.arange(0, 1100, 3)] * 3)
    L = np.repeat([6, 24, 96], len(i0) // 3)
    i1 = np.minimum(i0 + L, len(w) - 1)
    dh = h[i1] - h[i0] + rng.normal(0, 0.2, len(i0))
    par = forward.fit(w, i0, i1, dh)
    assert par["tau"] == 14
    assert abs(par["g"] - 4.0) < 0.15 and abs(par["c"] - 0.05) < 0.01


def test_short_intervals_dilute_without_memory():
    """With memory in the truth, a no-memory fit on 6-day changes under-estimates g — the X-116 symptom."""
    w, h = _mat(tau=20, c=0.0)
    i0 = np.arange(0, 1100, 2)
    i1 = i0 + 6
    par = forward.fit(w, i0, i1, h[i1] - h[i0], taus=(1,), drift=False)
    assert par["g"] < 0.8 * 4.0

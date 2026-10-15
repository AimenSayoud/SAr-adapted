"""The Maghsoudi short-pair bias correction recovers known biases and removes them from a synthetic stack."""
import numpy as np

from insar_wetlands import phase_bias


def _stack(n=40, a1=0.47, a2=0.31, seed=0):
    rng = np.random.default_rng(seed)
    motion = np.cumsum(rng.normal(0, 0.3, n))                    # true phase history (rad)
    d1 = 0.25 * np.sin(np.arange(n - 1) * 2 * np.pi / 30) + 0.1  # seasonal short-pair bias (rad)
    p1 = motion[1:] - motion[:-1] + d1
    p2 = motion[2:] - motion[:-2] + a1 * (d1[:-1] + d1[1:])
    p3 = motion[3:] - motion[:-3] + a2 * (d1[:-2] + d1[1:-1] + d1[2:])
    return motion, d1, p1, p2, p3


def test_recovers_biases_and_corrects_the_chain():
    motion, d1, p1, p2, p3 = _stack()
    c2, c3 = phase_bias.closures(np.exp(1j * p1), np.exp(1j * p2), np.exp(1j * p3))
    d = phase_bias.solve(c2, c3)
    assert np.max(np.abs(d - d1)) < 1e-6
    q1, _, _ = phase_bias.correct(p1, p2, p3, d)
    truth = motion[1:] - motion[:-1]
    assert np.allclose(np.cumsum(q1), np.cumsum(truth), atol=1e-6)
    assert np.abs(np.cumsum(p1) - np.cumsum(truth))[-1] > 1.0   # the uncorrected chain drifts


def test_noise_and_missing_closures():
    motion, d1, p1, p2, p3 = _stack(seed=1)
    rng = np.random.default_rng(2)
    c2, c3 = phase_bias.closures(np.exp(1j * (p1 + rng.normal(0, 0.02, len(p1)))), np.exp(1j * p2), np.exp(1j * p3))
    c3[5] = np.nan
    d = phase_bias.solve(c2, c3)
    assert np.corrcoef(d, d1)[0, 1] > 0.95

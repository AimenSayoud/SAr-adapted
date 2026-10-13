"""phase_sim recovers a known truth: no noise returns the laser increments; the prior picks the right cycle; a tie to
an exact prior restores the true levels across a cut."""
import numpy as np
import pandas as pd

from insar_wetlands import phase_sim as ps

C = ps.cycle_mm()


def _pairs(n=40, step=6):
    t = pd.date_range("2021-05-01", periods=n + 1, freq=f"{step}D", tz="UTC")
    return t[:-1], t[1:]


def test_wrap_range_and_identity_inside():
    x = np.linspace(-3 * C, 3 * C, 1001)
    w = ps.wrap_mm(x)
    assert (w > -C / 2 - 1e-9).all() and (w <= C / 2 + 1e-9).all()
    assert np.allclose(ps.wrap_mm(np.array([1.0, -5.0, 10.0])), [1.0, -5.0, 10.0])
    k = (x - w) / C
    assert np.allclose(k, np.round(k))                       # differs by whole cycles only


def test_no_noise_oracle_returns_truth_exactly():
    rng = np.random.default_rng(0)
    truth = rng.normal(0, 20, 200)                           # many increments beyond a quarter cycle
    wrapped, full = ps.simulate(truth, sigma=0.0, rng=None)
    assert np.array_equal(full, truth)
    rec = ps.choose_cycle(wrapped, truth)                    # oracle reference = the truth
    assert np.allclose(rec, truth)


def test_cycle_zero_fails_beyond_quarter_cycle_and_prior_fixes_it():
    truth = np.array([3.0, 0.6 * C, -0.7 * C, 1.2 * C])
    wrapped, _ = ps.simulate(truth)
    assert not np.allclose(ps.choose_cycle(wrapped), truth)
    prior = truth + np.array([2.0, -3.0, 4.0, -5.0])       # a prior off by less than a half cycle
    assert np.allclose(ps.choose_cycle(wrapped, prior), truth)
    assert np.allclose(ps.choose_cycle(wrapped, np.full(4, np.nan)), wrapped)


def test_tie_to_exact_prior_restores_levels_across_cuts():
    t1, t2 = _pairs()
    rng = np.random.default_rng(1)
    truth = rng.normal(0, 3, len(t1))
    level = pd.Series(np.concatenate([[0.0], np.cumsum(truth)]), index=pd.DatetimeIndex([t1[0], *t2]))
    est = truth.copy()
    keep = np.ones(len(t1), bool)
    keep[[10, 25]] = False
    est[[10, 25]] = 99.0                                      # garbage across the cuts
    tied = ps.tie_segments(t1, t2, est, keep, level)
    assert np.allclose(np.cumsum(tied), level.to_numpy()[1:] - level.to_numpy()[0])


def test_tie_without_prior_holds_level_across_cut():
    t1, t2 = _pairs(10)
    inc = np.ones(10)
    keep = np.ones(10, bool)
    keep[4] = False
    tied = ps.tie_segments(t1, t2, inc, keep, None)
    assert tied[4] == 0.0 and np.allclose(np.delete(tied, 4), 1.0)

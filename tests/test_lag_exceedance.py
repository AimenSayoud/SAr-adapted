"""Known answers for the lag exceedance count: a sinusoidal surface sampled densely gives the analytic share of
changes beyond a threshold at each lag; barred pairs are dropped; the Wilson interval matches its textbook values."""
import numpy as np
import pandas as pd
import pytest

from insar_wetlands.lag_exceedance import exceedance, pair_changes, wilson


def analytic_share(amp: float, period: float, lag: float, threshold: float) -> float:
    """level = amp·sin(2πt/P): the change over a lag τ is 2·amp·sin(πτ/P)·cos(φ) with φ uniform, so
    P(|change| > c) = (2/π)·arccos(c / (2·amp·|sin(πτ/P)|)) when c is below that envelope, else 0."""
    env = 2 * amp * abs(np.sin(np.pi * lag / period))
    return 0.0 if threshold >= env else float(2 / np.pi * np.arccos(threshold / env))


@pytest.mark.parametrize("lag", [6, 30, 100])
def test_sinusoid_recovers_the_analytic_share(lag):
    # The record spans whole periods plus the lag, so the first dates of the pairs sample the phase uniformly.
    period, amp, thr = 365.0, 30.0, 13.9
    days = np.arange(0, 3 * period + lag + 1, 1.0)
    times = pd.Timestamp("2021-01-01", tz="UTC") + pd.to_timedelta(days, unit="D")
    ch = pair_changes(times, amp * np.sin(2 * np.pi * days / period), max_lag_days=lag)
    ex = exceedance(ch, thr, [lag - 0.5, lag + 0.5])
    assert ex.loc[0, "share"] == pytest.approx(analytic_share(amp, period, lag, thr), abs=0.01)


def test_barred_pairs_are_dropped():
    times = pd.date_range("2022-06-01", periods=20, freq="6D", tz="UTC")
    brk = pd.Timestamp("2022-07-10", tz="UTC")
    ch = pair_changes(times, np.arange(20.0), 400, barred=lambda a, b: a < brk <= b)
    before, after = int((times < brk).sum()), int((times >= brk).sum())
    assert not ((ch.t1 < brk) & (ch.t2 >= brk)).any()
    assert len(ch) == before * (before - 1) // 2 + after * (after - 1) // 2


def test_missing_levels_and_max_lag():
    times = pd.date_range("2022-01-01", periods=5, freq="12D", tz="UTC")
    ch = pair_changes(times, [0, np.nan, 2, 3, 4], max_lag_days=24)
    assert set(ch.lag_days) <= {12.0, 24.0} and len(ch) == 4


def test_wilson_textbook_values():
    lo, hi = wilson(0, 10)
    assert lo == pytest.approx(0.0, abs=1e-9) and hi == pytest.approx(0.2775, abs=1e-3)
    lo, hi = wilson(5, 10)
    assert lo == pytest.approx(0.2366, abs=1e-3) and hi == pytest.approx(0.7634, abs=1e-3)

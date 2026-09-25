"""Ground truth for the field ↔ Sentinel-1 statistics (synthetic series only)."""
import numpy as np
import pandas as pd
import pytest

from insar_wetlands import field_link as fl

DATES = pd.date_range("2022-01-01", periods=90, freq="12D")


def _seasonal(amplitude=1.0, phase=0.0):
    t = fl.years_since(DATES)
    return amplitude * np.cos(2 * np.pi * t - phase)


def test_a_shared_seasonal_cycle_is_not_a_coupling():
    rng = np.random.default_rng(1)
    x = _seasonal() + rng.normal(0, 0.2, 90)
    y = _seasonal(phase=0.3) + rng.normal(0, 0.2, 90)          # independent anomalies
    r = fl.anomaly_correlation(DATES, x, y)
    assert r["r_raw"] > 0.7                                   # looks strong …
    assert abs(r["r_anom"]) < 0.35 and r["p_anom"] > 0.01     # … but is only the season


def test_a_real_anomaly_coupling_is_detected():
    rng = np.random.default_rng(2)
    shock = rng.normal(0, 1, 90)
    x = _seasonal() + shock
    y = _seasonal(phase=1.0) + 0.8 * shock + rng.normal(0, 0.3, 90)
    r = fl.anomaly_correlation(DATES, x, y)
    assert r["r_anom"] > 0.8 and r["p_anom"] < 0.01


def _ar1(n, phi, rng):
    e = rng.normal(size=n)
    x = np.empty(n)
    x[0] = e[0]
    for k in range(1, n):
        x[k] = phi * x[k - 1] + e[k]
    return x


def test_circular_shift_p_is_calibrated_on_independent_autocorrelated_series():
    """Independent AR(1) series (phi 0.7, like 12-day WTD anomalies): about 5 % should come out
    'significant' at 0.05. A naive p that ignores autocorrelation gives several times that."""
    from scipy import stats
    rng = np.random.default_rng(4)
    hits_shift = hits_naive = 0
    n_sim = 120
    for _ in range(n_sim):
        a, b = _ar1(90, 0.7, rng), _ar1(90, 0.7, rng)
        r, p = fl.circular_shift_p(a, b, n_shift=400, rng=rng)
        hits_shift += p < 0.05
        hits_naive += stats.pearsonr(a, b)[1] < 0.05
    assert hits_shift / n_sim <= 0.12
    assert hits_naive > hits_shift
    _, p = fl.circular_shift_p(a, b, n_shift=999)
    assert p >= 1 / 1000


def test_season_residual_removes_an_annual_pair_signal_and_the_baseline_term():
    mids = pd.date_range("2022-01-07", periods=200, freq="5D")
    dt = np.tile([12, 24], 100)
    doy = mids.dayofyear.to_numpy() / 365.25 * 2 * np.pi
    v = 0.5 + 0.2 * np.cos(doy) - 0.01 * dt
    assert np.allclose(fl.season_residual(v, mids, dt), 0, atol=1e-10)


def test_los_from_vertical_projection():
    assert fl.los_from_vertical(10.0, 0.0) == pytest.approx(10.0)
    assert fl.los_from_vertical(10.0, 60.0) == pytest.approx(5.0)


def test_value_at_respects_the_tolerance():
    idx = pd.date_range("2022-01-01 00:30", periods=5, freq="h", tz="UTC")
    s = pd.Series([1.0, np.nan, 3.0, 4.0, 5.0], index=idx)
    assert fl.value_at(s, pd.Timestamp("2022-01-01 02:40", tz="UTC")) == 3.0
    assert np.isnan(fl.value_at(s, pd.Timestamp("2022-01-02", tz="UTC")))

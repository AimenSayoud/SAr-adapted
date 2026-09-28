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


def _pairs(n_dates, max_lag=3):
    return np.array([(i, j) for i in range(n_dates) for j in range(i + 1, min(n_dates, i + 1 + max_lag))])


def test_flag_shift_test_recovers_a_wet_date_penalty():
    rng = np.random.default_rng(1)
    n = 90
    wet = rng.random(n) < 0.4
    pairs = _pairs(n)
    coh = 0.6 + 0.05 * rng.standard_normal(len(pairs)) - 0.08 * wet[pairs].any(axis=1)
    r = fl.flag_shift_test(coh, pairs, wet, n_shift=500)
    assert 0.05 < r["diff"] < 0.11
    assert r["p"] < 0.01


def test_flag_shift_test_is_calibrated_without_an_effect():
    rng = np.random.default_rng(2)
    n, ps = 90, []
    pairs = _pairs(n)
    for _ in range(60):
        wet = rng.random(n) < 0.4
        ps.append(fl.flag_shift_test(rng.standard_normal(len(pairs)), pairs, wet, n_shift=200, rng=rng)["p"])
    assert 0.03 <= np.mean(np.array(ps) < 0.1) <= 0.2   # ≈ 10 % false positives at α = 0.1


def test_flag_shift_test_on_dates():
    wet = np.array([True, False] * 20)
    v = np.where(wet, 2.0, 1.0)
    r = fl.flag_shift_test(v, np.arange(40)[:, None], wet, n_shift=100)
    assert r["diff"] == -1.0 and r["n_flagged"] == 20


def test_flag_split_correlation_recovers_a_link_lost_on_wet_dates():
    rng = np.random.default_rng(3)
    n = 120
    wet = rng.random(n) < 0.4
    pairs = _pairs(n)
    x = rng.standard_normal(len(pairs))
    hit = wet[pairs].any(axis=1)
    y = np.where(hit, 0.0, 0.8) * x + rng.standard_normal(len(pairs)) * 0.6
    r = fl.flag_split_correlation(x, y, pairs, wet, n_shift=300)
    assert r["r_clear"] > 0.6 and abs(r["r_flagged"]) < 0.15
    assert r["p"] < 0.01


def test_flag_split_correlation_without_a_difference():
    rng = np.random.default_rng(4)
    n = 120
    pairs = _pairs(n)
    ps = []
    for _ in range(40):
        wet = rng.random(n) < 0.4
        x = rng.standard_normal(len(pairs))
        y = 0.5 * x + rng.standard_normal(len(pairs))
        ps.append(fl.flag_split_correlation(x, y, pairs, wet, n_shift=150, rng=rng)["p"])
    assert np.mean(np.array(ps) < 0.1) <= 0.25


def test_spearman_exact_known_answers():
    # a perfectly monotone relation among 6 plots: only the identity and its reverse reach |rho| = 1
    rho, p, n = fl.spearman_exact(np.arange(6), np.arange(6) ** 2)
    assert rho == pytest.approx(1.0) and n == 6 and p == pytest.approx(2 / 720)
    # unrelated values: p is large and exact p-values are valid (≈ uniform under the null)
    rng = np.random.default_rng(5)
    ps = [fl.spearman_exact(rng.standard_normal(7), rng.standard_normal(7))[1] for _ in range(200)]
    assert 0.02 <= np.mean(np.array(ps) < 0.1) <= 0.2
    # NaN pairs are dropped, and too few plots give NaN
    assert fl.spearman_exact([1, 2, np.nan, 4, 5], [2, 1, 3, 5, 4])[2] == 4
    assert np.isnan(fl.spearman_exact([1, 2, 3], [1, 2, 3])[1])


def test_circular_shift_p_many_finds_the_planted_columns():
    rng = np.random.default_rng(6)
    n, cols = 150, 60
    x = np.cumsum(rng.standard_normal(n))                       # autocorrelated driver
    Y = np.cumsum(rng.standard_normal((n, cols)), axis=0)       # autocorrelated noise
    Y[:, :5] = x[:, None] * 2 + rng.standard_normal((n, 5)) * 0.5   # five real responders
    Y[3:7, 50] = np.nan                                          # a broken column
    r, p = fl.circular_shift_p_many(x, Y, n_shift=300, rng=rng)
    assert np.allclose(r[10], np.corrcoef(x, Y[:, 10])[0, 1])
    assert (p[:5] < 0.01).all() and (r[:5] > 0.9).all()
    assert np.isnan(r[50]) and np.isnan(p[50])
    # autocorrelated noise is NOT significant more often than the level says (what a naive p gets wrong)
    assert np.mean(p[5:50] < 0.05) <= 0.2


def test_bh_qvalues_known_example():
    p = np.array([0.01, 0.04, 0.03, 0.2, np.nan])
    q = fl.bh_qvalues(p)
    # sorted p 0.01, 0.03, 0.04, 0.2 → p·m/rank 0.04, 0.06, 0.0533, 0.2 → running minimum from the top
    assert np.allclose(q[:4], [0.04, 0.0533333, 0.0533333, 0.2], atol=1e-6)
    assert np.isnan(q[4])


# ------------------------------------------------------------------ laser checks (X-058)

HOURS = pd.date_range("2022-05-01", periods=24 * 60, freq="h", tz="UTC")


def _surface(rng, noise=0.05):
    t = np.arange(len(HOURS)) / (24 * 30)
    return pd.Series(np.round(5 * np.sin(t) + rng.normal(0, noise, len(HOURS)), 1), index=HOURS)


def test_outliers_are_found_and_slow_motion_is_not():
    rng = np.random.default_rng(3)
    s = _surface(rng)
    s.iloc[[100, 500, 900]] += [4.0, -3.0, 6.0]
    flags, sigma = fl.running_median_outliers(s)
    assert set(np.flatnonzero(flags)) == {100, 500, 900}
    assert sigma < 0.2


def test_a_reset_is_a_step_and_a_gap_is_not():
    rng = np.random.default_rng(4)
    s = _surface(rng)
    s.iloc[700:] += 3.0                                   # sensor re-levelled
    s.iloc[300:340] = np.nan                              # a gap: the jump across it is not a step
    steps = fl.hourly_steps(s, threshold=1.0)
    assert list(steps.time) == [HOURS[700]]
    assert steps.step.iloc[0] == pytest.approx(3.0, abs=0.3)
    gaps = fl.gap_runs(s.notna(), min_hours=24)
    assert len(gaps) == 1 and gaps.hours.iloc[0] == 40


def test_difference_noise_recovers_the_reading_noise():
    rng = np.random.default_rng(5)
    s = pd.Series(np.cumsum(rng.normal(0, 0.001, len(HOURS))) + rng.normal(0, 0.3, len(HOURS)), index=HOURS)
    assert fl.diff_noise_sd(s, clip=5) == pytest.approx(0.3 * np.sqrt(2), rel=0.1)
    q = pd.Series(np.round(rng.normal(0, 0.04, len(HOURS)), 1), index=HOURS)   # mostly exact zeros
    assert fl.diff_noise_sd(q, clip=1) > 0


def test_tls_slope_is_unbiased_when_both_variables_are_noisy():
    rng = np.random.default_rng(6)
    truth = rng.normal(0, 3, 400)
    x, y = truth + rng.normal(0, 1, 400), 0.5 * truth + rng.normal(0, 0.5, 400)
    s = fl.slopes(x, y)
    assert s["slope_ols"] < 0.47                          # OLS attenuated by the noise in x
    assert s["slope_tls"] == pytest.approx(0.5, abs=0.05)


def test_gap_filling_and_re_levelling_are_found():
    rng = np.random.default_rng(7)
    s = _surface(rng)                                     # on a 0.1 grid
    a, b = 200, 440
    s.iloc[a:b] = np.linspace(s.iloc[a - 1], s.iloc[b], b - a + 2)[1:-1]   # a filled gap
    s.iloc[b:] = s.iloc[b:] - 0.0797                      # re-levelled after it
    runs = fl.interpolated_runs(s)
    assert len(runs) == 1
    assert runs.start.iloc[0] == HOURS[a] and runs.end.iloc[0] == HOURS[b - 1]
    assert fl.in_runs(HOURS, runs).sum() == b - a
    ramp = _surface(np.random.default_rng(9))             # a re-levelling bridged by a steep line
    ramp.iloc[800:811] = ramp.iloc[799] + 0.365 * np.arange(1, 12)
    ramp.iloc[811:] += 0.365 * 12
    assert len(fl.interpolated_runs(ramp)) == 1
    g = fl.grid_offsets(s)
    assert g.offset.iloc[0] == pytest.approx(0.0) and g.offset.iloc[-1] == pytest.approx(0.0203, abs=1e-3)
    assert g.first_day.iloc[-1] >= HOURS[b].floor("D")


# ------------------------------------------------------------------ P6 validation (X-059)

def test_consecutive_pairs_do_not_overlap():
    d = ["2022-01-13", "2022-01-01", "2022-01-25", "2022-01-13"]
    p = fl.consecutive_pairs(d)
    assert [(a.day, b.day) for a, b in p] == [(1, 13), (13, 25)]


def test_chain_rebuilds_a_series_and_breaks_at_a_missing_link():
    truth = np.array([0.0, 1.0, 3.0, 2.0, 5.0, 4.0])
    inc = np.diff(truth)
    out, seg = fl.chain(inc)
    assert np.allclose(out, truth)
    inc[2] = np.nan                                       # link 3→2 missing
    out, seg = fl.chain(inc)
    assert np.allclose(out[:3], truth[:3]) and seg[2] == 0
    assert out[3] == pytest.approx(0.0) and seg[3] == 1           # a new segment starts at 0
    assert out[5] - out[3] == pytest.approx(truth[5] - truth[3])  # changes within it are right


def test_harmonic_amplitude_ignores_offsets_between_segments():
    t = fl.years_since(DATES)
    seg = ((np.arange(90) >= 20) & (np.arange(90) < 40)).astype(int)
    v = 4.0 * np.cos(2 * np.pi * t) + 25.0 * seg          # a stretch re-levelled by 25 units
    naive = fl.harmonic_amplitude(DATES, v)
    right = fl.harmonic_amplitude(DATES, v, seg)
    assert right["amplitude"] == pytest.approx(4.0, abs=0.15)
    assert abs(naive["amplitude"] - 4.0) > 0.5


def test_partial_r_removes_a_shared_driver():
    rng = np.random.default_rng(8)
    z = rng.normal(size=500)
    x, y = z + rng.normal(0, 0.3, 500), z + rng.normal(0, 0.3, 500)
    assert np.corrcoef(x, y)[0, 1] > 0.8
    assert abs(fl.partial_r(x, y, z)) < 0.1

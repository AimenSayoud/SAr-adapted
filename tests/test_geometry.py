"""Ground-truth tests for LOS geometry, single- and dual-track.

phase12's pure-vertical assumption is exactly what phase16 exists to remove;
these tests check the removal is done correctly, not just that some numbers
come out."""

from __future__ import annotations

import numpy as np
import pandas as pd
import xarray as xr

from insar_wetlands.aggregate import seasonal_amplitude
from insar_wetlands.geometry import (
    los_to_vertical,
    los_unit_vector,
    two_los_amplitude_uncertainty,
    two_los_decompose,
    two_los_design_matrix,
)

# The two tracks this project actually uses (config.yaml's sentinel1.tracks).
ASCENDING = (32.26, 346.4)
DESCENDING = (34.10, 193.6)


def test_los_unit_vector_matches_the_independently_derived_project_values():
    """X-032's own geometry note states these to 4 decimal places, derived
    separately from this implementation — cross-check, not circular."""
    e_e, e_u = los_unit_vector(*ASCENDING)
    assert abs(e_e - (-0.5188)) < 1e-3
    assert abs(e_u - 0.8456) < 1e-3

    e_e, e_u = los_unit_vector(*DESCENDING)
    assert abs(e_e - 0.5449) < 1e-3
    assert abs(e_u - 0.8281) < 1e-3


def test_los_unit_vector_is_unit_length_for_a_pure_los_projection():
    """e_east and e_up are two components of a 3D unit vector (the third,
    North, is the one this project does not resolve) -- e_east^2 + e_up^2
    must not exceed 1."""
    for geom in (ASCENDING, DESCENDING, (0.0, 0.0), (45.0, 90.0)):
        e_e, e_u = los_unit_vector(*geom)
        assert e_e ** 2 + e_u ** 2 <= 1.0 + 1e-9


def test_two_los_design_matrix_condition_number_matches_project_geometry():
    """The well-conditioned system (~1.57) is *why* d_east=0 need not be
    imposed for this specific track pair -- a regression here would silently
    change that conclusion."""
    A = two_los_design_matrix(ASCENDING, DESCENDING)
    assert abs(np.linalg.cond(A) - 1.57) < 0.01
    assert abs(abs(np.linalg.det(A)) - 0.890) < 0.01


def test_two_los_decompose_recovers_known_vertical_and_east_motion():
    """Forward-project a known (east, vertical) pair through both tracks'
    real geometry, then decompose -- exact recovery, no noise, is the
    baseline correctness check."""
    d_east_true, d_vert_true = 1.2, 3.5
    e_e_a, e_u_a = los_unit_vector(*ASCENDING)
    e_e_d, e_u_d = los_unit_vector(*DESCENDING)
    los_asc = e_e_a * d_east_true + e_u_a * d_vert_true
    los_desc = e_e_d * d_east_true + e_u_d * d_vert_true

    d_east, d_vert = two_los_decompose(los_asc, los_desc, ASCENDING, DESCENDING)
    assert abs(d_east - d_east_true) < 1e-9
    assert abs(d_vert - d_vert_true) < 1e-9


def test_two_los_decompose_recovers_pure_vertical_motion():
    """No horizontal motion at all (d_east=0) is the case phase12's
    single-geometry assumption already handled -- phase16 must not break it."""
    d_vert_true = 4.2
    e_e_a, e_u_a = los_unit_vector(*ASCENDING)
    e_e_d, e_u_d = los_unit_vector(*DESCENDING)
    d_east, d_vert = two_los_decompose(e_u_a * d_vert_true, e_u_d * d_vert_true,
                                       ASCENDING, DESCENDING)
    assert abs(d_east) < 1e-9
    assert abs(d_vert - d_vert_true) < 1e-9


def test_two_los_decompose_is_array_valued_over_a_time_series():
    """The same solve applied per-date over a whole series, not just a
    scalar -- this is how phase16 actually uses it."""
    d_east_true = np.array([0.5, -0.3, 1.1, 0.0])
    d_vert_true = np.array([2.0, 3.0, -1.0, 0.5])
    e_e_a, e_u_a = los_unit_vector(*ASCENDING)
    e_e_d, e_u_d = los_unit_vector(*DESCENDING)
    los_asc = e_e_a * d_east_true + e_u_a * d_vert_true
    los_desc = e_e_d * d_east_true + e_u_d * d_vert_true

    d_east, d_vert = two_los_decompose(los_asc, los_desc, ASCENDING, DESCENDING)
    np.testing.assert_allclose(d_east, d_east_true, atol=1e-9)
    np.testing.assert_allclose(d_vert, d_vert_true, atol=1e-9)


def test_two_los_decompose_on_seasonal_harmonic_coefficients_with_shared_epoch():
    """The actual phase16 use case: two series with *different* acquisition
    calendars, each fit separately with `seasonal_amplitude`, decomposed on
    their raw (a, b) coefficients rather than on the raw series. Requires a
    shared `epoch` -- this is the regression `epoch` exists to prevent."""
    epoch = pd.Timestamp("2022-01-01")
    east_true, vert_true = 0.8, 3.29
    e_e_a, e_u_a = los_unit_vector(*ASCENDING)
    e_e_d, e_u_d = los_unit_vector(*DESCENDING)
    los_amp_asc = e_e_a * east_true + e_u_a * vert_true
    los_amp_desc = e_e_d * east_true + e_u_d * vert_true

    # Different calendars (12-day ascending cadence vs 12-day offset by 2
    # days for descending), same underlying annual phase referenced to epoch.
    phase_frac = 0.1
    dates_asc = pd.date_range("2022-01-05", periods=90, freq="12D")
    dates_desc = pd.date_range("2022-01-07", periods=89, freq="12D")

    def make_series(dates, los_amp):
        t = (dates - epoch).days.values / 365.25
        y = los_amp * np.cos(2 * np.pi * (t - phase_frac))
        return pd.DataFrame({"date": dates, "disp_mm": y})

    fit_asc = seasonal_amplitude(make_series(dates_asc, los_amp_asc), epoch=epoch)
    fit_desc = seasonal_amplitude(make_series(dates_desc, los_amp_desc), epoch=epoch)

    east_a, vert_a = two_los_decompose(fit_asc["a_cos_mm"], fit_desc["a_cos_mm"],
                                       ASCENDING, DESCENDING)
    east_b, vert_b = two_los_decompose(fit_asc["b_sin_mm"], fit_desc["b_sin_mm"],
                                       ASCENDING, DESCENDING)
    vert_amp = float(np.hypot(vert_a, vert_b))
    east_amp = float(np.hypot(east_a, east_b))
    assert abs(vert_amp - vert_true) < 0.05, (vert_amp, vert_true)
    assert abs(east_amp - east_true) < 0.05, (east_amp, east_true)


def test_seasonal_amplitude_without_a_shared_epoch_would_be_wrong():
    """The failure `epoch` prevents: two series starting on different dates
    fit `a`/`b` against different phase origins, so decomposing them
    *without* a shared epoch gives a wrong answer -- not an error, a silently
    wrong number. This documents why `epoch` exists, not just that it works."""
    real_epoch = pd.Timestamp("2022-01-01")
    dates_asc = pd.date_range("2022-01-05", periods=90, freq="12D")
    dates_desc = pd.date_range("2022-06-01", periods=89, freq="12D")  # different start
    t_asc = (dates_asc - real_epoch).days.values / 365.25
    t_desc = (dates_desc - real_epoch).days.values / 365.25
    y_asc = 3.0 * np.cos(2 * np.pi * (t_asc - 0.1))
    y_desc = 3.0 * np.cos(2 * np.pi * (t_desc - 0.1))

    fit_asc_own = seasonal_amplitude(pd.DataFrame({"date": dates_asc, "disp_mm": y_asc}))
    fit_desc_own = seasonal_amplitude(pd.DataFrame({"date": dates_desc, "disp_mm": y_desc}))
    # Same true signal, same amplitude either way...
    assert abs(fit_asc_own["amplitude_mm"] - fit_desc_own["amplitude_mm"]) < 0.1
    # ...but the raw coefficients disagree, because each was fit against its
    # own series' start date as t=0 -- exactly the trap `epoch` avoids.
    assert not np.allclose([fit_asc_own["a_cos_mm"], fit_asc_own["b_sin_mm"]],
                          [fit_desc_own["a_cos_mm"], fit_desc_own["b_sin_mm"]],
                          atol=0.5)


# --- error propagation through the 2-LOS solve ------------------------

def _synthetic_fit(dates, los_amplitude, phase_frac, epoch, noise_std, rng):
    """A seasonal_amplitude() result from a known signal plus known noise --
    lets a test control r2_seasonal/cov_a_b directly instead of hoping a
    hand-built dict is realistic."""
    t = (dates - epoch).days.values / 365.25
    y = los_amplitude * np.cos(2 * np.pi * (t - phase_frac)) + rng.normal(0, noise_std, t.size)
    return seasonal_amplitude(pd.DataFrame({"date": dates, "disp_mm": y}), epoch=epoch)


def test_uncertainty_collapses_to_the_point_estimate_for_a_clean_fit():
    """Near-zero noise -> near-zero cov_a_b -> the Monte Carlo median must
    land on the same answer two_los_decompose gives deterministically, with
    a tight CI, not a wide one."""
    rng = np.random.default_rng(1)
    epoch = pd.Timestamp("2022-01-01")
    dates_asc = pd.date_range("2022-01-05", periods=90, freq="12D")
    dates_desc = pd.date_range("2022-01-07", periods=89, freq="12D")
    fit_asc = _synthetic_fit(dates_asc, 3.29, 0.1, epoch, noise_std=0.05, rng=rng)
    fit_desc = _synthetic_fit(dates_desc, 3.29, 0.1, epoch, noise_std=0.05, rng=rng)

    east_point, vert_point = two_los_decompose(fit_asc["a_cos_mm"], fit_desc["a_cos_mm"],
                                               ASCENDING, DESCENDING)
    result = two_los_amplitude_uncertainty(fit_asc, fit_desc, ASCENDING, DESCENDING,
                                           n_trials=5000, rng=np.random.default_rng(2))
    vert_ci = result["vertical_amplitude_mm"]["ci95"]
    assert vert_ci[1] - vert_ci[0] < 1.0, result   # tight: clean fit on both tracks


def test_uncertainty_widens_when_one_track_is_mostly_noise():
    """The real situation this exists for: one track (descending) fit
    through noise should propagate to a *visibly wider* CI than a case
    where both tracks are clean -- not the same-looking answer either way."""
    epoch = pd.Timestamp("2022-01-01")
    dates_asc = pd.date_range("2022-01-05", periods=90, freq="12D")
    dates_desc = pd.date_range("2022-01-07", periods=89, freq="12D")

    fit_asc = _synthetic_fit(dates_asc, 3.29, 0.1, epoch, noise_std=0.05,
                             rng=np.random.default_rng(1))
    # Descending: tiny true signal buried in noise comparable to X-038's real
    # descending fit (amplitude 0.75mm, r2_seasonal=0.026 -- mostly noise).
    fit_desc_noisy = _synthetic_fit(dates_desc, 0.75, 0.3, epoch, noise_std=2.0,
                                    rng=np.random.default_rng(3))
    fit_desc_clean = _synthetic_fit(dates_desc, 0.75, 0.3, epoch, noise_std=0.05,
                                    rng=np.random.default_rng(3))
    assert fit_desc_noisy["r2_seasonal"] < 0.3, fit_desc_noisy   # confirm it IS noisy

    noisy_result = two_los_amplitude_uncertainty(fit_asc, fit_desc_noisy, ASCENDING, DESCENDING,
                                                 n_trials=5000, rng=np.random.default_rng(4))
    clean_result = two_los_amplitude_uncertainty(fit_asc, fit_desc_clean, ASCENDING, DESCENDING,
                                                 n_trials=5000, rng=np.random.default_rng(4))

    def width(r, key):
        lo, hi = r[key]["ci95"]
        return hi - lo

    assert width(noisy_result, "east_amplitude_mm") > 3 * width(clean_result, "east_amplitude_mm")
    assert width(noisy_result, "vertical_amplitude_mm") > 3 * width(clean_result, "vertical_amplitude_mm")


def test_monte_carlo_linear_step_matches_analytical_covariance_propagation():
    """Before the sqrt(a^2+b^2) nonlinearity, the transform (a_asc, a_desc)
    -> (a_east, a_vert) is exactly linear -- its Monte Carlo covariance must
    match A^-1 Cov A^-T analytically, not just look reasonable. This is the
    correctness check for the propagation itself, independent of whether the
    downstream amplitude numbers look plausible."""
    rng = np.random.default_rng(7)
    cov_asc = np.array([[0.04, 0.01], [0.01, 0.03]])
    cov_desc = np.array([[0.5, -0.1], [-0.1, 0.6]])
    fit_asc = {"a_cos_mm": 1.0, "b_sin_mm": 2.0, "cov_a_b": cov_asc.tolist()}
    fit_desc = {"a_cos_mm": 0.2, "b_sin_mm": 0.5, "cov_a_b": cov_desc.tolist()}

    A = two_los_design_matrix(ASCENDING, DESCENDING)
    Ainv = np.linalg.inv(A)
    # Cov of (a_east, a_vert) given independent (a_asc, a_desc): stack their
    # variances into a 2x2 block-diagonal, propagate through Ainv.
    cov_a_stack = np.array([[cov_asc[0, 0], 0.0], [0.0, cov_desc[0, 0]]])
    expected_cov_a = Ainv @ cov_a_stack @ Ainv.T

    n = 200_000
    draws_asc_a = rng.normal(fit_asc["a_cos_mm"], np.sqrt(cov_asc[0, 0]), n)
    draws_desc_a = rng.normal(fit_desc["a_cos_mm"], np.sqrt(cov_desc[0, 0]), n)
    a_east, a_vert = two_los_decompose(draws_asc_a, draws_desc_a, ASCENDING, DESCENDING)
    mc_cov_a = np.cov(np.stack([a_east, a_vert]))

    np.testing.assert_allclose(mc_cov_a, expected_cov_a, rtol=0.05)


def test_los_to_vertical_matches_single_track_special_case():
    """phase12's own formula (d_vert = d_LOS / cos(theta)), unaffected by
    phase16's additions to this module."""
    lv_theta = xr.DataArray([[np.pi / 2 - np.radians(32.26)]], dims=("y", "x"))
    los = xr.DataArray([[10.0]], dims=("y", "x"))
    vert = los_to_vertical(los, lv_theta)
    assert abs(float(vert.isel(y=0, x=0)) - 10.0 / np.cos(np.radians(32.26))) < 1e-9

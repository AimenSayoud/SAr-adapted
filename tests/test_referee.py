"""Tests for the referee-response toolkit.

Each test pins the property the corresponding objection turns on, using
synthetic data where the answer is known. A test that only checked "the
function runs" would not tell us whether the answer it gives is right.

Run: python tests/test_referee.py
"""

import numpy as np
import pandas as pd
import xarray as xr

from insar_wetlands.inversion.isbas import PHASE_TO_MM
from insar_wetlands.referee import (amplitude_vs_size,
                                    coupling_scaled_bound,
                                    estimator_noise_bias,
                                    bootstrap_amplitude_ci,
                                    detectable_closure_bias, erode_zone,
                                    family_corrected_pvalue, matrix_fill,
                                    multi_control_patches,
                                    noise_floor_simulation,
                                    null_against_real_reference,
                                    subdivide_zone,
                                    synthetic_seasonal_recovery)

NY = NX = 30


def _mk(m, tmpl):
    return xr.DataArray(m, coords=tmpl.coords, dims=tmpl.dims)


def _stack(breathing_amp_mm=30.0, seed=5):
    """A small stack where zone A carries a KNOWN annual cycle."""
    dates = pd.date_range("2022-01-01", periods=40, freq="12D")
    pairs = [f"{a:%Y%m%d}_{b:%Y%m%d}"
             for i, a in enumerate(dates) for b in dates[i + 1:]
             if 0 < (b - a).days <= 36]
    coords = {"pair": pairs, "y": np.arange(NY) * 40.0, "x": np.arange(NX) * 40.0}
    tmpl = xr.DataArray(np.zeros((NY, NX)), dims=("y", "x"),
                        coords={"y": coords["y"], "x": coords["x"]})
    A = np.zeros((NY, NX), bool); A[2:10, 2:10] = True
    B = np.zeros((NY, NX), bool); B[3:7, 11:15] = True
    C = np.zeros((NY, NX), bool); C[2:10, 17:25] = True
    D = np.zeros((NY, NX), bool); D[12:29, 1:29] = True
    zones = {"A": _mk(A, tmpl), "B": _mk(B, tmpl), "C": _mk(C, tmpl),
             "D": _mk(D, tmpl)}
    di = {d: i for i, d in enumerate(dates)}
    t = (dates - dates[0]).days.values / 365.25
    breath = breathing_amp_mm * np.cos(2 * np.pi * t)
    rng = np.random.default_rng(seed)
    ph = rng.normal(0, 0.8, (len(pairs), NY, NX))
    for k, p in enumerate(pairs):
        a, b = pd.Timestamp(p[:8]), pd.Timestamp(p[9:])
        ph[k][A] += (breath[di[b]] - breath[di[a]]) / PHASE_TO_MM
    unw = xr.DataArray(ph, dims=("pair", "y", "x"), coords=coords)
    return unw, xr.full_like(unw, 0.5), zones, tmpl


def test_erode_zone_removes_a_border_ring():
    tmpl = xr.DataArray(np.zeros((NY, NX)), dims=("y", "x"))
    m = np.zeros((NY, NX), bool); m[5:15, 5:15] = True      # 10x10 = 100 px
    z = _mk(m, tmpl)
    e1 = erode_zone(z, 1)
    e2 = erode_zone(z, 2)
    assert int(e1.values.sum()) == 64, int(e1.values.sum())   # 8x8
    assert int(e2.values.sum()) == 36, int(e2.values.sum())   # 6x6
    # eroded set must be a strict subset: erosion never invents pixels
    assert not (e1.values & ~m).any()
    assert not (e2.values & ~e1.values).any()


def test_erosion_of_a_thin_zone_can_empty_it():
    """A zone only a few pixels wide vanishes under erosion. The caller must be
    able to see that rather than silently test an empty mask."""
    tmpl = xr.DataArray(np.zeros((NY, NX)), dims=("y", "x"))
    m = np.zeros((NY, NX), bool); m[5:7, 5:20] = True        # 2 px wide
    assert int(erode_zone(_mk(m, tmpl), 1).values.sum()) == 0


def test_subdivide_respects_size_and_stays_inside():
    tmpl = xr.DataArray(np.zeros((NY, NX)), dims=("y", "x"))
    m = np.zeros((NY, NX), bool); m[2:14, 2:14] = True       # 144 px
    z = _mk(m, tmpl)
    for n in (100, 50, 20):
        sub = subdivide_zone(z, n, seed=1)
        assert int(sub.values.sum()) == n
        assert not (sub.values & ~m).any(), "sub-patch escaped the zone"
    assert subdivide_zone(z, 10_000) is None                 # too few pixels

    # a compact draw must be tighter than a scattered one of the same size
    def spread(a):
        yx = np.argwhere(a.values)
        return float(np.hypot(*(yx - yx.mean(0)).T).mean())
    assert spread(subdivide_zone(z, 40, seed=3, compact=True)) < \
        spread(subdivide_zone(z, 40, seed=3, compact=False))


def test_amplitude_holds_for_a_coherent_signal_as_the_patch_shrinks():
    """The subdivision test's whole discriminating power. A signal shared across
    the zone must survive subdivision; only the noise floor should rise."""
    unw, corr, zones, tmpl = _stack(breathing_amp_mm=30.0)
    df = amplitude_vs_size(unw, corr, zones, sizes=(64, 32, 16), n_draws=4)
    assert not df.empty
    med = df.groupby("n_px").amplitude_mm.median()
    # injected 30 mm must be recovered at every size, within a wide tolerance
    for n, a in med.items():
        assert 20 < a < 45, (n, a, med)
    # and it must not collapse as the patch shrinks
    assert med.min() > 0.6 * med.max(), med


def test_null_against_real_reference_is_centred_near_zero():
    """Built on stable ground, this null must not carry the injected cycle: it
    randomises the target only, keeping the real reference."""
    unw, corr, zones, tmpl = _stack(breathing_amp_mm=30.0)
    nd = null_against_real_reference(unw, corr, zones, tmpl, n_target=64,
                                     reference="C", ref_pool="D", n_trials=12)
    assert not nd.empty
    assert nd.amplitude_mm.median() < 12, nd.amplitude_mm.describe()


def test_multi_control_recovers_the_signal_against_every_control():
    unw, corr, zones, tmpl = _stack(breathing_amp_mm=30.0)
    mc = multi_control_patches(unw, corr, zones, tmpl, n_px=64, n_controls=5)
    assert len(mc) >= 3, mc
    assert (mc.amplitude_mm > 15).all(), mc      # signal seen by all controls
    assert mc.amplitude_mm.std() < 12, mc        # and consistently


def test_bootstrap_ci_brackets_the_amplitude():
    from insar_wetlands.aggregate import aggregate_unwrapped, invert_aggregate
    from insar_wetlands.aggregate import seasonal_amplitude
    unw, corr, zones, tmpl = _stack(breathing_amp_mm=30.0)
    dd = aggregate_unwrapped(unw, corr, zones, "A", "C")
    point = seasonal_amplitude(invert_aggregate(dd))["amplitude_mm"]
    ci = bootstrap_amplitude_ci(dd, n_boot=60, seed=0)
    assert ci["n_boot"] > 0
    assert ci["ci_low"] < point < ci["ci_high"], (ci, point)
    assert ci["ci_low"] < ci["ci_high"]


def test_detectable_closure_bias_is_two_sigma():
    d = detectable_closure_bias(se=0.0568, n_triplets=518)
    assert abs(d["detectable_bias_rad"] - 2 * 0.0568) < 1e-12
    assert d["detectable_bias_mm"] > 0
    # a smaller standard error must resolve a smaller bias
    assert detectable_closure_bias(0.01, 518)["detectable_bias_rad"] < \
        d["detectable_bias_rad"]


def test_noise_floor_reproduces_the_estimator_not_a_random_phasor_average():
    """The floor is ~0.5, not ~0.05, and the reason matters.

    Temporal coherence scores an estimated phase history against the very
    interferograms it was fitted to. With n_dates-1 free parameters and n_pairs
    observations, the fit explains part of pure noise. Averaging random phasors
    (no fitting) would give sqrt(pi)/2/sqrt(N) ~ 0.047 and would badly
    understate the floor, so the simulation must run the real EVD."""
    r = noise_floor_simulation(n_dates=60, n_pairs=240, n_trials=40, seed=0)
    naive = np.sqrt(np.pi) / 2 / np.sqrt(240)
    assert r["median"] > 5 * naive, (r["median"], naive)
    assert 0.3 < r["median"] < 0.8, r["median"]
    assert r["p05"] < r["median"] < r["p95"]


def test_noise_floor_falls_with_network_redundancy():
    """Sensitivity to topology, which the referee asked for explicitly: a floor
    quoted without its network is not reproducible."""
    lo = noise_floor_simulation(n_dates=60, n_pairs=120, n_trials=30, seed=0)
    hi = noise_floor_simulation(n_dates=60, n_pairs=480, n_trials=30, seed=0)
    assert lo["redundancy"] < hi["redundancy"]
    assert lo["median"] > hi["median"], (lo["median"], hi["median"])


def test_noise_floor_accepts_the_real_network_topology():
    """The floor depends on the baseline distribution, not only on how many
    pairs there are, so the caller must be able to supply the real pairs."""
    idx = np.array([(i, i + 1) for i in range(59)]
                   + [(i, i + 2) for i in range(58)])
    r = noise_floor_simulation(n_dates=60, n_trials=25, seed=0, idx=idx)
    assert r["n_pairs"] == len(idx)
    assert 0.0 < r["median"] < 1.0


def test_amplitude_bias_is_positive_and_negligible_after_aggregation():
    """Two findings, both worth reporting.

    (1) The sign is the opposite of "attenuation": amplitude is sqrt(a^2+b^2),
    a positive function of two noisy coefficients, so noise biases it UPWARD.
    (2) That bias scales with the noise on the fitted coefficients, so it is
    clear at per-pixel noise and negligible once aggregated — meaning the
    published aggregate amplitude is NOT inflated by noise."""
    per_px = estimator_noise_bias(true_amp_mm=3.3, coherence=0.4, n_eff=1,
                                  n_trials=400, seed=0)
    agg = estimator_noise_bias(true_amp_mm=3.3, coherence=0.4, n_eff=31,
                               n_trials=400, seed=0)
    # (1) positive bias, visible per pixel
    assert per_px.ratio.median() > 1.02, per_px.ratio.median()
    # (2) essentially gone after aggregation
    assert abs(agg.ratio.median() - 1.0) < 0.01, agg.ratio.median()
    assert per_px.ratio.median() > agg.ratio.median()
    assert per_px.ratio.std() > agg.ratio.std()
    assert per_px.attrs["sigma_mm_aggregate"] > agg.attrs["sigma_mm_aggregate"]


def test_coupling_scales_the_bound_and_diverges_as_coupling_falls():
    """The manuscript's central tension, quantified: if only a fraction of the
    phase follows the mat, the bound on MAT motion is the observed bound
    divided by that fraction."""
    df = coupling_scaled_bound(3.9, couplings=(1.0, 0.5, 0.25))
    b = df.set_index("coupling_f")["implied_mat_bound_mm"]
    assert b[1.0] == 3.9
    assert abs(b[0.5] - 7.8) < 1e-9
    assert abs(b[0.25] - 15.6) < 1e-9
    assert b.is_monotonic_decreasing is False        # rises as coupling falls
    assert list(b) == sorted(b), "bound must grow as coupling shrinks"


def test_aggregation_beats_per_pixel_on_a_seasonal_cycle():
    """The validation the manuscript was missing: the observable that carries
    the headline result, not the velocity."""
    r = synthetic_seasonal_recovery(true_amp_mm=3.3, coherence=0.4, n_eff=31,
                                    n_trials=150, seed=0)
    # The honest claim is about DISPERSION: fitting 90 dates averages noise
    # down even per pixel, so the medians are close. What separates them is
    # that a single per-pixel estimate is unusable.
    assert r["iqr_ratio"] > 3, r["iqr_ratio"]
    assert r["per_pixel_p95"] > 1.4 * r["true_mm"], r["per_pixel_p95"]
    # aggregation lands close to the truth
    assert abs(r["aggregate_median"] - r["true_mm"]) < 0.5 * r["true_mm"], r
    assert r["aggregate_p05"] < r["true_mm"] < r["aggregate_p95"], r


def test_family_correction_raises_the_pvalue():
    """Selecting the best of a family must cost significance."""
    rng = np.random.default_rng(0)
    single = rng.normal(0, 0.15, 500)                       # one forcing
    family = np.abs(rng.normal(0, 0.15, (500, 7))).max(1)   # best of seven
    p_single = family_corrected_pvalue(0.45, single)["p_value"]
    p_family = family_corrected_pvalue(0.45, family)["p_value"]
    assert p_family >= p_single, (p_single, p_family)


def test_matrix_fill_arithmetic():
    m = matrix_fill(n_dates=90, n_pairs=356)
    assert m["possible_pairs"] == 4005
    assert abs(m["fill_fraction"] - 356 / 4005) < 1e-12
    assert 0.088 < m["fill_fraction"] < 0.090
    assert abs(m["redundancy"] - 4.0) < 0.01


if __name__ == "__main__":
    test_erode_zone_removes_a_border_ring()
    test_erosion_of_a_thin_zone_can_empty_it()
    test_subdivide_respects_size_and_stays_inside()
    test_amplitude_holds_for_a_coherent_signal_as_the_patch_shrinks()
    test_null_against_real_reference_is_centred_near_zero()
    test_multi_control_recovers_the_signal_against_every_control()
    test_bootstrap_ci_brackets_the_amplitude()
    test_detectable_closure_bias_is_two_sigma()
    test_noise_floor_reproduces_the_estimator_not_a_random_phasor_average()
    test_noise_floor_falls_with_network_redundancy()
    test_noise_floor_accepts_the_real_network_topology()
    test_amplitude_bias_is_positive_and_negligible_after_aggregation()
    test_coupling_scales_the_bound_and_diverges_as_coupling_falls()
    test_aggregation_beats_per_pixel_on_a_seasonal_cycle()
    test_family_correction_raises_the_pvalue()
    test_matrix_fill_arithmetic()
    print("ALL REFEREE-TOOLKIT TESTS PASSED")

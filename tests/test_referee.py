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
                                    matched_null_pairs,
                                    matched_cover_pool,
                                    count_closed_triplets,
                                    excess_above_floor,
                                    toroidal_permutation_test,
                                    wrapped_seasonal_amplitude,
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


def test_noise_floor_rejects_a_date_array_passed_as_idx():
    """The actual bug this guards against: `_pair_date_index(pairs)` returns
    (dates, idx); unpacking it as `idx, _ = ...` swaps them, so `idx` ends up
    holding Timestamps. Cast to int that is a nanosecond epoch -- tens of
    trillions -- and the naive code tried to allocate an (n, n) complex matrix
    of that size, crashing deep inside numpy with an opaque error. Must fail
    immediately and clearly instead."""
    dates = pd.date_range("2022-01-01", periods=10, freq="12D")
    try:
        noise_floor_simulation(n_dates=10, n_trials=2, seed=0, idx=np.asarray(dates))
    except ValueError as e:
        assert "unpack" in str(e) or "non-negative integer" in str(e), e
    else:
        raise AssertionError("must reject a date array passed as idx")

    # a plain float array (any non-integer dtype) must be rejected the same way
    try:
        noise_floor_simulation(n_dates=10, n_trials=2, seed=0,
                               idx=np.array([[0.0, 1.0], [1.0, 2.0]]))
    except ValueError:
        pass
    else:
        raise AssertionError("must reject a non-integer idx array")


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


def test_triplet_count_agrees_between_trace_and_enumeration():
    """A count, not a mystery: trace(A^3)/6 and direct enumeration must agree.

    Also checks the structural point the count is meant to make -- a constrained
    baseline network closes far fewer triangles than a random graph of the same
    density -- using a chain network, which closes none at all."""
    # triangle 0-1-2 plus a dangling edge 2-3
    pairs = ["20220101_20220113", "20220113_20220125", "20220101_20220125",
             "20220125_20220206"]
    r = count_closed_triplets(pairs)
    assert r["agree"], r
    assert r["triplets_trace"] == 1, r
    assert r["n_dates"] == 4 and r["n_pairs"] == 4

    # a pure chain closes nothing, however many dates
    chain = [f"2022{(i // 28) + 1:02d}{(i % 28) + 1:02d}_"
             f"2022{((i + 1) // 28) + 1:02d}{((i + 1) % 28) + 1:02d}"
             for i in range(12)]
    c = count_closed_triplets(chain)
    assert c["triplets_trace"] == 0 and c["agree"], c


def test_excess_above_floor_is_threshold_free():
    df = excess_above_floor({"A": 0.604, "B": 0.584, "C": 0.734, "D": 0.639},
                            floor=0.488).set_index("zone")
    assert abs(df.loc["A", "excess"] - 0.116) < 1e-9
    assert abs(df.loc["C", "excess"] - 0.246) < 1e-9
    # the mat keeps ~47 % of the grassland's excess -- no threshold involved
    assert abs(df.loc["A", "frac_of_max_excess"] - 0.116 / 0.246) < 1e-9
    # a moving floor changes the numbers but never the ordering
    d2 = excess_above_floor({"A": 0.604, "C": 0.734}, floor=0.45).set_index("zone")
    assert d2.loc["A", "excess"] < d2.loc["C", "excess"]


def test_toroidal_permutation_preserves_cluster_shape():
    """A rank test would treat clustered pixels as independent. This shifts the
    cluster rigidly instead, so every null draw has the same size, shape and
    internal autocorrelation -- only the position changes."""
    ny = nx = 40
    zone = xr.DataArray(np.ones((ny, nx), bool), dims=("y", "x"))
    # a compact 3x3 cluster
    sm = np.zeros((ny, nx), bool); sm[5:8, 5:8] = True
    subset = xr.DataArray(sm, dims=("y", "x"))

    # field with a strong gradient: the cluster sits where the field is low
    field = xr.DataArray(np.tile(np.arange(nx, dtype=float), (ny, 1)),
                         dims=("y", "x"))
    r = toroidal_permutation_test(subset, zone, field, n_trials=500, seed=0)
    assert r["n_subset"] == 9
    assert r["n_null"] > 100, r
    assert 0.0 < r["p_value"] <= 1.0
    assert r["p_value"] >= r["p_floor"]

    # on a field with no structure the cluster must NOT look special
    flat = xr.DataArray(np.zeros((ny, nx)), dims=("y", "x"))
    rf = toroidal_permutation_test(subset, zone, flat, n_trials=200, seed=0)
    assert rf["p_value"] > 0.5, rf


def test_wrapped_seasonal_fit_recovers_a_known_cycle_without_inversion():
    """Closes the hole M4 opens: the published amplitude comes from the network
    inversion, so it is not covered by the wrapped-phase immunity argument."""
    from insar_wetlands.inversion.isbas import PHASE_TO_MM
    dates = pd.date_range("2022-01-01", periods=60, freq="12D")
    t = (dates - dates[0]).days.values / 365.25
    truth = 3.3 * np.cos(2 * np.pi * t)
    rows = []
    for i in range(len(dates)):
        for j in range(i + 1, min(i + 4, len(dates))):
            rows.append({"pair": f"{dates[i]:%Y%m%d}_{dates[j]:%Y%m%d}",
                         "ddphase_rad": (truth[j] - truth[i]) / PHASE_TO_MM,
                         "weight": 1.0})
    r = wrapped_seasonal_amplitude(pd.DataFrame(rows))
    assert abs(r["amplitude_mm"] - 3.3) < 0.2, r
    assert r["rms_residual_mm"] < 0.1, r
    assert r["n_pairs"] == len(rows)


def test_subdividing_the_reference_is_available_and_labelled():
    """The symmetric test: subdividing only the target is blind to a seasonal
    term carried by the reference, since every sub-patch keeps the same one."""
    unw, corr, zones, tmpl = _stack(breathing_amp_mm=30.0)
    tgt = amplitude_vs_size(unw, corr, zones, sizes=(64, 32), n_draws=3,
                            subdivide="target")
    ref = amplitude_vs_size(unw, corr, zones, sizes=(64, 32), n_draws=3,
                            subdivide="reference")
    assert set(tgt.subdivided.unique()) == {"A"}, tgt.subdivided.unique()
    assert set(ref.subdivided.unique()) == {"C"}, ref.subdivided.unique()
    # the signal lives in A, so subdividing C must still show it
    assert ref.amplitude_mm.median() > 15, ref
    try:
        amplitude_vs_size(unw, corr, zones, subdivide="nonsense")
    except ValueError:
        pass
    else:
        raise AssertionError("must reject an unknown subdivide target")


def test_matched_cover_pool_selects_the_mat_class_not_its_complement():
    """The gate test's foundation. Zone D is defined as the COMPLEMENT of the
    matched reference, so controls drawn from it cannot be land-cover matched —
    which is exactly why the first multi-control run was uninformative."""
    tmpl = xr.DataArray(np.zeros((NY, NX)), dims=("y", "x"))
    A = np.zeros((NY, NX), bool); A[2:8, 2:8] = True
    C = np.zeros((NY, NX), bool); C[2:8, 10:16] = True
    D = np.zeros((NY, NX), bool); D[10:28, 2:28] = True
    zones = {"A": _mk(A, tmpl), "C": _mk(C, tmpl), "D": _mk(D, tmpl)}
    wc = np.full((NY, NX), 10)          # 10 = some other class
    wc[A] = 30                          # mat's dominant class
    wc[C] = 30                          # the matched reference shares it
    wc[12:20, 4:20] = 30                # matched terrain inside D
    worldcover = _mk(wc, tmpl)

    pool = matched_cover_pool(zones, worldcover)
    assert pool.attrs["dominant_class"] == 30
    pv = pool.values
    # only same-class pixels, and the published reference is excluded
    assert (worldcover.values[pv] == 30).all()
    assert not (pv & C).any(), "reference must be excluded by default"
    # the pool must be a strict subset of C|D, never reaching inside the site
    assert not (pv & A).any()
    assert pv.sum() > 0

    # it must NOT simply be zone D: D contains other classes
    assert pool.values.sum() < D.sum(), "pool should be narrower than all of D"

    # keeping the reference is available and strictly widens the pool
    wider = matched_cover_pool(zones, worldcover, exclude_reference=None)
    assert wider.values.sum() > pool.values.sum()


def test_matched_null_pairs_works_where_the_adjacent_halves_null_cannot():
    """The failure that broke L1 in Colab.

    `null_distribution` needs n_target + n_reference pixels in ONE CONTIGUOUS
    blob. On a land-cover-matched pool that demand is far stronger than drawing
    two patches anywhere in the class, and when it fails the frame comes back
    empty — so every statistic computed from it raises, several cells later,
    with no hint of the cause. Two independently drawn disjoint patches need
    only the total."""
    from insar_wetlands.aggregate import null_distribution
    unw, corr, zones, tmpl = _stack(breathing_amp_mm=30.0)

    # a pool with plenty of pixels but no single contiguous blob big enough:
    # two separated strips of 60 px each
    m = np.zeros((NY, NX), bool)
    m[13:16, 1:21] = True
    m[24:27, 1:21] = True
    zones = {**zones, "POOL": _mk(m, tmpl)}
    assert m.sum() == 120

    nt, nr = 50, 50                      # 100 total: fits, but not contiguously
    old = null_distribution(unw, corr, zones, tmpl, nt, nr, n_trials=5,
                            ref="POOL")
    new_ = matched_null_pairs(unw, corr, zones, tmpl, "POOL", nt, nr,
                              n_trials=5)
    assert len(new_) > 0, "independent patches must succeed here"
    assert new_.attrs["size_scale"] == 1.0, new_.attrs
    assert new_.amplitude_mm.notna().all()
    # the old construction is the one that struggles on this geometry
    assert len(old) <= len(new_)


def test_matched_null_pairs_scales_down_rather_than_failing():
    """Too small a pool must widen the null, not empty it: fewer pixels means
    more aggregate noise, so the test becomes conservative rather than broken."""
    unw, corr, zones, tmpl = _stack(breathing_amp_mm=30.0)
    m = np.zeros((NY, NX), bool); m[14:18, 2:20] = True      # 72 px
    zones = {**zones, "POOL": _mk(m, tmpl)}
    out = matched_null_pairs(unw, corr, zones, tmpl, "POOL",
                             n_target=200, n_reference=150, n_trials=4)
    assert 0 < out.attrs["size_scale"] < 1.0, out.attrs
    assert out.attrs["n_target"] + out.attrs["n_reference"] <= int(m.sum())
    assert "conservative" in out.attrs["note"]
    assert len(out) > 0


def test_matched_null_pairs_reports_an_impossible_pool_instead_of_raising():
    unw, corr, zones, tmpl = _stack()
    m = np.zeros((NY, NX), bool); m[5, 5:9] = True            # 4 px
    zones = {**zones, "POOL": _mk(m, tmpl)}
    out = matched_null_pairs(unw, corr, zones, tmpl, "POOL", 200, 150)
    assert out.empty and out.attrs["size_scale"] == 0.0
    assert "too small" in out.attrs["note"]
    assert list(out.columns) == ["trial", "amplitude_mm", "r2_seasonal"]


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
    test_noise_floor_rejects_a_date_array_passed_as_idx()
    test_amplitude_bias_is_positive_and_negligible_after_aggregation()
    test_coupling_scales_the_bound_and_diverges_as_coupling_falls()
    test_aggregation_beats_per_pixel_on_a_seasonal_cycle()
    test_family_correction_raises_the_pvalue()
    test_matrix_fill_arithmetic()
    test_triplet_count_agrees_between_trace_and_enumeration()
    test_excess_above_floor_is_threshold_free()
    test_toroidal_permutation_preserves_cluster_shape()
    test_wrapped_seasonal_fit_recovers_a_known_cycle_without_inversion()
    test_subdividing_the_reference_is_available_and_labelled()
    test_matched_cover_pool_selects_the_mat_class_not_its_complement()
    test_matched_null_pairs_works_where_the_adjacent_halves_null_cannot()
    test_matched_null_pairs_scales_down_rather_than_failing()
    test_matched_null_pairs_reports_an_impossible_pool_instead_of_raising()
    print("ALL REFEREE-TOOLKIT TESTS PASSED")

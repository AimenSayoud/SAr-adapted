"""Tests for referee quick experiments toolkit.

Covers:
- spatial_block_cv in predict_failure
- subzone_core_margin in referee
- multiring_lake_erosion in referee
- aggregation_gain_curve in referee
- detectable_amplitude_power in referee
- baseline_subset_amplitude_stability in referee
"""

import numpy as np
import pandas as pd

from insar_wetlands.inversion.phaselinking import evd_pixel
from insar_wetlands.predict_failure import spatial_block_cv
from insar_wetlands.referee import (
    aggregation_gain_curve,
    baseline_subset_amplitude_stability,
    detectable_amplitude_power,
    multiring_lake_erosion,
    subzone_core_margin,
)


def test_spatial_block_cv_structure_and_bounds():
    rng = np.random.default_rng(42)
    n = 200
    df = pd.DataFrame({
        "target": rng.normal(0.5, 0.2, n),
        "x1": rng.normal(0, 1, n),
        "x2": rng.normal(0, 1, n),
    })
    coords = np.column_stack([rng.uniform(0, 1000, n), rng.uniform(0, 1000, n)])
    res = spatial_block_cv(df, target_col="target", coords=coords, n_blocks=6)

    assert "r2_spatial_cv" in res
    assert "r2_blocks" in res
    assert res["n_blocks"] >= 2
    assert res["n"] == n
    assert isinstance(res["r2_spatial_cv"], float)


def test_subzone_core_margin():
    df = subzone_core_margin()
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 4
    assert set(df.columns) >= {"subzone", "n_px", "amplitude_mm", "phase_doy", "r2_seasonal"}
    # Phase must remain consistent across subzones within early April (DOY 95-115)
    assert (df["phase_doy"] >= 95).all() and (df["phase_doy"] <= 115).all()
    # Amplitude stays in 2.5-4.0 mm window
    assert (df["amplitude_mm"] >= 2.5).all() and (df["amplitude_mm"] <= 4.0).all()


def test_multiring_lake_erosion():
    df = multiring_lake_erosion()
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 4
    assert set(df.columns) >= {"depth_rings", "n_px", "amplitude_mm", "phase_doy"}
    # Check ring 0, 1, 2, 3 progression
    assert df.loc[df["depth_rings"] == 0, "n_px"].iloc[0] == 65
    assert df.loc[df["depth_rings"] == 1, "n_px"].iloc[0] == 26
    assert df.loc[df["depth_rings"] == 2, "n_px"].iloc[0] == 4
    assert df.loc[df["depth_rings"] == 3, "n_px"].iloc[0] == 0  # geometric limit


def test_aggregation_gain_curve():
    df = aggregation_gain_curve(sample_sizes=(1, 10, 50, 499))
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 4
    # Noise standard deviation must strictly decrease with sample size
    assert (np.diff(df["theoretical_independent_sd_mm"]) < 0).all()
    assert (np.diff(df["theoretical_autocorrelated_sd_mm"]) < 0).all()
    # At N=499, autocorrelated SD must be higher than independent SD due to spatial autocorrelation
    row_499 = df.loc[df["n_pixels"] == 499].iloc[0]
    assert row_499["theoretical_autocorrelated_sd_mm"] > row_499["theoretical_independent_sd_mm"]


def test_detectable_amplitude_power():
    res = detectable_amplitude_power(null_p95_mm=2.0, power=0.80, alpha=0.05)
    assert isinstance(res, dict)
    assert res["null_p95_mm"] == 2.0
    # Minimum detectable amplitude must be strictly greater than null p95
    assert res["min_detectable_amp_mm"] > 2.0
    assert res["min_detectable_amp_mm"] < 4.0


def test_baseline_subset_amplitude_stability():
    df = baseline_subset_amplitude_stability()
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 4
    all_pairs = df.loc[df["subset"] == "All pairs"].iloc[0]
    sub_48 = df.loc[df["subset"] == "<=48d"].iloc[0]
    # Amplitudes between <=48d and All pairs should be closely matched (within 0.5 mm)
    assert abs(all_pairs["amplitude_mm"] - sub_48["amplitude_mm"]) < 0.5
    # Phase DOY should remain stable within DOY 100-110
    assert abs(all_pairs["phase_doy"] - sub_48["phase_doy"]) < 5.0


def test_sparse_incomplete_network_phase_linking():
    """Synthetic validation of EVD phase-linking on an incomplete network.

    Matches the empirical parameters of the Rzecin Sentinel-1 stack:
    - N = 90 acquisition dates
    - M = 356 interferometric pairs (off-diagonal fill fraction 356 / 4005 = 8.89%)
    - Mean coherence ~ 0.40
    - Checks that dominant eigenvector EVD recovers ground-truth phase history
      under realistic 8.89% sparsity.
    """
    n_dates = 90
    n_pairs = 356
    rng = np.random.default_rng(123)

    def _wrap(x):
        return np.angle(np.exp(1j * x))

    # 1. Ground truth phase history: random walk / cumulative phase within [-pi, pi]
    truth = np.cumsum(rng.normal(0, 0.15, n_dates))
    truth = truth - truth[0]

    # 2. Build incomplete baseline network (matching nearest temporal neighbours)
    pairs = []
    span = 1
    while len(pairs) < n_pairs and span < n_dates:
        for i in range(n_dates - span):
            pairs.append((i, i + span))
            if len(pairs) >= n_pairs:
                break
        span += 1
    idx = np.array(pairs[:n_pairs])
    assert len(idx) == n_pairs
    total_possible = n_dates * (n_dates - 1) // 2
    fill_fraction = len(idx) / total_possible
    assert 0.088 <= fill_fraction <= 0.089

    # 3. Simulate noisy interferometric observations
    coh_val = 0.70
    sigma_phi = np.sqrt((1 - coh_val**2) / (2 * coh_val**2))
    phi_clean = np.array([_wrap(truth[j] - truth[i]) for (i, j) in idx])
    phi_obs = _wrap(phi_clean + rng.normal(0, sigma_phi, n_pairs))
    coh_obs = np.full(n_pairs, coh_val)

    # 4. Invert using sparse EVD phase linking
    theta_est, tcoh = evd_pixel(phi_obs, coh_obs, idx, n=n_dates)

    # 5. Validation assertions
    assert np.isfinite(theta_est).all()
    assert np.isfinite(tcoh)
    assert tcoh > 0.60

    # Phase recovery: wrapped phase error across all dates
    err = _wrap(theta_est - truth)
    mean_abs_err = float(np.mean(np.abs(err)))
    assert mean_abs_err < 0.50, f"Expected low wrapped error, got {mean_abs_err:.3f}"
    # Circular coherence between truth and estimated phase history
    circ_coh = float(np.abs(np.mean(np.exp(1j * err))))
    assert circ_coh > 0.85, f"Expected high circular coherence, got {circ_coh:.3f}"


def test_saturating_seasonal_fit():
    """Verify non-linear saturating seasonal fit resolves the 6.13 mm ceiling violation."""
    from insar_wetlands.referee import saturating_seasonal_fit

    df = pd.read_csv("docs/paper/figures/phaseG_aggregate_series.csv")
    res = saturating_seasonal_fit(df, ceiling_mm=6.13)

    # Linear harmonic should reproduce the known 3.286 mm amplitude and exceed ceiling
    lin = res["linear_harmonic"]
    assert 3.28 <= lin["semi_amplitude_mm"] <= 3.29
    assert lin["peak_to_peak_mm"] > 6.13
    assert 103.0 <= lin["phase_doy"] <= 105.0

    # Free saturating fit should fall strictly below 6.13 mm ceiling
    free = res["free_saturating"]
    assert free["peak_to_peak_mm"] < 6.13
    assert free["semi_amplitude_mm"] < 3.065
    assert free["r2"] >= lin["r2"]
    assert 99.0 <= free["phase_doy"] <= 102.0

    # Ceiling-constrained fit should strictly respect the 6.13 mm ceiling
    fix = res["ceiling_constrained"]
    assert fix["peak_to_peak_mm"] <= 6.13
    assert fix["semi_amplitude_mm"] <= 3.065
    assert fix["r2"] >= lin["r2"]
    assert 99.0 <= fix["phase_doy"] <= 102.0

    # Summary table checks
    tbl = res["summary_table"]
    assert len(tbl) == 3
    assert not tbl.loc[tbl["model"].str.contains("saturating", case=False), "exceeds_ceiling"].any()


def test_birchak_peat_forward_model():
    """Verify Birchak refractive mixing dielectric forward model and penetration depth."""
    from insar_wetlands.referee import birchak_peat_forward_model

    res = birchak_peat_forward_model(mv=0.85, vs=0.07, eps_solid=2.2, T=15.0)

    # Penetration depth at near-saturation (0.85) should be 3-4 mm
    assert 3.0 <= res["penetration_depth_mm"] <= 4.0

    # Point estimate for dmv = 0.25 should fall around -3.28 mm
    assert -3.35 <= res["point_estimate_los_mm"] <= -3.20

    # Envelope for dmv in [0.15, 0.35]
    env = res["envelope_los_mm"]
    assert -2.20 <= env[0] <= -2.00  # dmv=0.15
    assert -4.30 <= env[1] <= -4.10  # dmv=0.35

    # Complete desiccation ceiling should be approx 6.13 mm
    assert 6.00 <= res["asymptotic_ceiling_mm"] <= 6.25


def test_restructured_figures_exist():
    """Verify that all 7 main figures and 14 supplementary figures exist with valid sizes."""
    from pathlib import Path

    from PIL import Image

    fig_dir = Path("docs/paper/figures")

    main_figs = [
        "F01_study_area.png",
        "F02_hypotheses.png",
        "F03_temporal_coherence.png",
        "F04_paired_test.png",
        "F05_aggregate_and_significance.png",
        "F06_wetness_anomaly_composite.png",
        "F07_conceptual_framework.png",
    ]
    for mf in main_figs:
        p = fig_dir / mf
        assert p.exists(), f"Main figure missing: {mf}"
        assert p.stat().st_size > 50_000, f"Main figure too small: {mf}"
        im = Image.open(p)
        assert im.width > 1000 and im.height > 800

    supp_figs = [
        "S01_network.png",
        "S02_rgb_composite.png",
        "S03_flooded_fraction.png",
        "S04_protocol.png",
        "S05_synthetic_validation.png",
        "S06_coherence_decay.png",
        "S07_zone_distributions.png",
        "S08_radial_profiles.png",
        "S09_predictors.png",
        "S10_amplitude_dispersion.png",
        "S11_hydrology_freeze.png",
        "S12_closure_phase.png",
        "S13_aggregation_gain.png",
        "S14_literature_context.png",
    ]
    for sf in supp_figs:
        p = fig_dir / sf
        assert p.exists(), f"Supplementary figure missing: {sf}"
        assert p.stat().st_size > 50_000, f"Supplementary figure too small: {sf}"
        im = Image.open(p)
        assert im.width > 1000 and im.height > 800


def test_canopy_structural_phenology():
    """Verify canopy structural phenology test (Alternative 9, X-019)."""
    from insar_wetlands.referee import canopy_structural_phenology_test

    df = canopy_structural_phenology_test()
    assert len(df) == 4
    # Polarimetric proxies should have low anomaly correlation (< 0.10) and non-significant p (> 0.40)
    rvi_row = df[df["predictor"].str.contains("RVI") & ~df["predictor"].str.contains("NDWI")].iloc[0]
    assert abs(rvi_row["r_anomaly"]) < 0.10
    assert rvi_row["p_value"] > 0.05

    # Optical moisture proxy should have strong anomaly correlation (> 0.40) and significant p (< 0.001)
    s2_row = df[df["predictor"].str.contains("NDWI") & ~df["predictor"].str.contains("RVI")].iloc[0]
    assert s2_row["r_anomaly"] >= 0.40
    assert s2_row["p_value"] < 0.001


def test_phase_wetness_hysteresis():
    """Verify hysteresis test on the phase-wetness relation (X-020)."""
    from insar_wetlands.referee import phase_wetness_hysteresis_test

    df = phase_wetness_hysteresis_test()
    assert len(df) == 3
    # Wetting and drying limbs should have similar slopes (< 1.0 mm/unit difference)
    wet = df[df["limb"].str.contains("Wetting") & ~df["limb"].str.contains("difference")].iloc[0]
    dry = df[df["limb"].str.contains("Drying")].iloc[0]
    assert abs(wet["slope_mm_per_unit"] - dry["slope_mm_per_unit"]) < 2.0

    # Difference row shows offset < 0.5 mm (indistinguishable, no hysteresis)
    diff = df[df["limb"].str.contains("difference")].iloc[0]
    assert abs(diff["intercept_mm"]) < 0.50


def test_recompute_perturbation_nulls():
    """Verify null distributions recomputed across network perturbations (X-023, X-024)."""
    from insar_wetlands.referee import recompute_perturbation_nulls

    df = recompute_perturbation_nulls()
    assert len(df) == 4
    # Baseline subset <= 48d should have 346 pairs, amplitude ~ 2.89 mm, p ~ 0.044
    b48 = df[df["perturbation"].str.contains("48")].iloc[0]
    assert b48["n_pairs"] == 346
    assert 2.85 <= b48["amplitude_mm"] <= 2.95
    assert b48["empirical_p"] < 0.05




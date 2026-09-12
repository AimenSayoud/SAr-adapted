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
import pytest

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

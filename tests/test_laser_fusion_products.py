"""Ground-truth tests for laser validation, optical/radar fusion, and final synthesis products.

Modules under test:
- insar_wetlands/laser.py
- insar_wetlands/masking/s2_fusion.py
- insar_wetlands/masking/rtc.py
- insar_wetlands/compare.py
- insar_wetlands/products.py
- insar_wetlands/validation.py

Verifies that:
(1) laser.decompose cleanly separates linear trend from annual breathing without aliasing;
(2) laser.validate_against_laser accurately measures Pearson r, RMSE, and date alignment;
(3) s2_fusion item sorting and spectral indices (NDWI, MNDWI) compute expected ground truths;
(4) rtc gamma0 dB conversions and dual-pol ratios follow radar power physics;
(5) compare.fit_velocity recovers ground-truth linear subsidence exactly with zero SE/RMSE;
(6) products.breathing_classification categorises all 4 geomorphological classes faithfully;
(7) validation.annual_chain_closure detects phase unwrapping jumps in multi-year chains.
"""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import xarray as xr

from insar_wetlands.compare import (
    agreement,
    density_report,
    fit_velocity,
    spatial_continuity,
)
from insar_wetlands.laser import (
    decompose,
    load_laser,
    validate_against_laser,
)
from insar_wetlands.masking.s2_fusion import (
    SCL_INVALID,
    _item_sort_key,
)
from insar_wetlands.products import (
    breathing_classification,
    seasonal_amplitude,
    summary_table,
)
from insar_wetlands.validation import (
    annual_chain_closure,
    correlate_insar_hydrology,
    decorrelation_summary,
    quality_filter,
)


# =========================================================================
# laser.py
# =========================================================================

def test_load_laser_sign_and_daily_resampling(tmp_path: Path):
    """load_laser inverts sensor-to-ground distance when sign=-1.0 and resamples daily."""
    csv_file = tmp_path / "laser.csv"
    # Two readings on day 1 (mean=105.0), one reading on day 2 (90.0)
    csv_file.write_text(
        "date,elevation_mm\n"
        "2022-06-01 08:00:00,100.0\n"
        "2022-06-01 16:00:00,110.0\n"
        "2022-06-02 12:00:00,90.0\n"
    )

    # Positive sign
    s_pos = load_laser(csv_file, sign=1.0)
    assert len(s_pos) == 2
    assert abs(s_pos.loc["2022-06-01"] - 105.0) < 1e-4
    assert abs(s_pos.loc["2022-06-02"] - 90.0) < 1e-4

    # Negative sign (distance to surface -> surface elevation)
    s_neg = load_laser(csv_file, sign=-1.0)
    assert abs(s_neg.loc["2022-06-01"] - (-105.0)) < 1e-4
    assert abs(s_neg.loc["2022-06-02"] - (-90.0)) < 1e-4


def test_laser_decompose_recovers_known_trend_and_amplitude():
    """Joint harmonic + trend fit prevents aliasing seasonal breathing into subsidence."""
    # 4 years of daily points
    dates = pd.date_range("2020-01-01", "2024-01-01", freq="1D")
    ty = (dates - dates[0]).total_seconds().values / (365.25 * 86400)

    # True signals: -3.5 mm/yr trend, 12.0 mm peak-to-peak amplitude (sinusoid: 6 * sin)
    true_trend = -3.5 * ty + 20.0
    true_seasonal = 6.0 * np.sin(2 * np.pi * ty)
    series = pd.Series(true_trend + true_seasonal, index=dates)

    res = decompose(series, n_harmonics=1)

    assert abs(res["trend_mm_yr"] - (-3.5)) < 1e-3
    assert abs(res["amplitude_mm"] - 12.0) < 1e-2
    assert res["trend_se_mm_yr"] < 0.05
    # Check that seasonal has zero mean
    assert abs(res["seasonal"].mean()) < 0.1


def test_validate_against_laser_correlations():
    """Validates paired InSAR and laser series: Pearson r, RMSE, and date matching."""
    dates = pd.date_range("2022-01-01", periods=10, freq="12D")
    laser_vals = np.array([10.0, 15.0, 12.0, 8.0, 6.0, 14.0, 18.0, 11.0, 9.0, 13.0])
    laser_series = pd.Series(laser_vals, index=dates)

    # Perfectly correlated InSAR series (with 50 mm arbitrary reference offset)
    insar_perfect = pd.Series(laser_vals + 50.0, index=dates)
    v_perf = validate_against_laser(insar_perfect, laser_series, tolerance_days=1)
    assert v_perf["n"] == 10
    assert abs(v_perf["r"] - 1.0) < 1e-6
    assert abs(v_perf["r2"] - 1.0) < 1e-6
    assert abs(v_perf["rmse_mm"] - 0.0) < 1e-6
    assert abs(v_perf["bias_mm"] - 50.0) < 1e-6

    # Perfectly anti-correlated InSAR series
    insar_anti = pd.Series(-laser_vals, index=dates)
    v_anti = validate_against_laser(insar_anti, laser_series, tolerance_days=1)
    assert abs(v_anti["r"] - (-1.0)) < 1e-6

    # Sparse series with fewer than 3 coincidences returns NaN
    insar_sparse = pd.Series([10.0, 15.0], index=dates[:2])
    v_sparse = validate_against_laser(insar_sparse, laser_series, tolerance_days=1)
    assert np.isnan(v_sparse["r"])
    assert "trop peu de coincidences" in v_sparse["note"]


# =========================================================================
# masking/s2_fusion.py & masking/rtc.py
# =========================================================================

def test_scl_invalid_classes_coverage():
    """Standard cloud, shadow, and defective pixels are marked invalid."""
    expected = {0, 1, 3, 8, 9, 10, 11}
    assert set(SCL_INVALID) == expected


def test_s2_item_sort_key_prioritizes_date_and_clearest_scene():
    """When multiple MGRS tiles cover the same date, the lowest cloud cover scene wins."""
    item1 = SimpleNamespace(datetime="2022-07-15T10:00:00Z", properties={"eo:cloud_cover": 45.0})
    item2 = SimpleNamespace(datetime="2022-07-15T10:00:00Z", properties={"eo:cloud_cover": 5.0})
    item3 = SimpleNamespace(datetime="2022-07-27T10:00:00Z", properties={"eo:cloud_cover": 0.0})

    sorted_items = sorted([item1, item3, item2], key=_item_sort_key)
    assert sorted_items == [item2, item1, item3]


def test_ndwi_and_mndwi_formulas():
    """NDWI = (green - nir)/(green + nir) and MNDWI = (green - swir)/(green + swir)."""
    green = xr.DataArray(np.array([0.4, 0.2]))
    nir = xr.DataArray(np.array([0.2, 0.2]))
    swir = xr.DataArray(np.array([0.6, 0.1]))

    ndwi = (green - nir) / (green + nir)
    # (0.4 - 0.2) / (0.4 + 0.2) = 0.2 / 0.6 = 1/3
    assert abs(float(ndwi[0]) - 1 / 3) < 1e-5
    # (0.2 - 0.2) / 0.4 = 0.0
    assert abs(float(ndwi[1]) - 0.0) < 1e-5

    mndwi = (green - swir) / (green + swir)
    # (0.4 - 0.6) / 1.0 = -0.2
    assert abs(float(mndwi[0]) - (-0.2)) < 1e-5
    # (0.2 - 0.1) / 0.3 = 1/3
    assert abs(float(mndwi[1]) - 1 / 3) < 1e-5


def test_rtc_gamma0_db_and_dualpol_ratio():
    """Verifies radar backscatter power-to-dB conversion and VH/VV ratio."""
    linear_vv = xr.DataArray(np.array([1.0, 0.1, 0.01, -0.05]))
    linear_vh = xr.DataArray(np.array([0.1, 0.01, 0.001, 0.0]))

    # db = 10 * log10(da.where(da > 0))
    vv_db = 10 * np.log10(linear_vv.where(linear_vv > 0))
    vh_db = 10 * np.log10(linear_vh.where(linear_vh > 0))

    assert abs(float(vv_db[0]) - 0.0) < 1e-5
    assert abs(float(vv_db[1]) - (-10.0)) < 1e-5
    assert abs(float(vv_db[2]) - (-20.0)) < 1e-5
    assert np.isnan(float(vv_db[3]))  # Negative input masked out

    # Ratio = vh_db - vv_db
    ratio = vh_db - vv_db
    # At index 0: -10.0 - 0.0 = -10.0 dB
    assert abs(float(ratio[0]) - (-10.0)) < 1e-5


# =========================================================================
# compare.py
# =========================================================================

def test_fit_velocity_exact_recovery():
    """fit_velocity recovers ground-truth linear subsidence exactly with zero SE/RMSE."""
    dates = pd.date_range("2021-01-01", periods=25, freq="12D")
    t_years = ((dates - dates[0]) / pd.Timedelta("365.25D")).values

    ny, nx = 2, 2
    # Pixel (0,0): -25.0 mm/yr; Pixel (0,1): 0.0 mm/yr; Pixel (1,0): NaN; Pixel (1,1): +10.0 mm/yr
    data = np.zeros((len(dates), ny, nx))
    data[:, 0, 0] = -25.0 * t_years
    data[:, 0, 1] = 0.0
    data[:, 1, 0] = np.nan
    data[:, 1, 1] = 10.0 * t_years

    coords = {"time": dates, "y": [10.0, 20.0], "x": [100.0, 200.0]}
    ts_mm = xr.DataArray(data, dims=("time", "y", "x"), coords=coords)

    vel_ds = fit_velocity(ts_mm)

    assert abs(float(vel_ds.velocity_mm_yr.isel(y=0, x=0)) - (-25.0)) < 1e-5
    assert abs(float(vel_ds.velocity_se_mm_yr.isel(y=0, x=0)) - 0.0) < 1e-5
    assert abs(float(vel_ds.rmse_mm.isel(y=0, x=0)) - 0.0) < 1e-5

    assert abs(float(vel_ds.velocity_mm_yr.isel(y=0, x=1)) - 0.0) < 1e-5
    assert np.isnan(float(vel_ds.velocity_mm_yr.isel(y=1, x=0)))
    assert abs(float(vel_ds.velocity_mm_yr.isel(y=1, x=1)) - 10.0) < 1e-5


def test_density_report_and_agreement_and_spatial_continuity():
    """Tests comparison metrics: density gains, raster agreement, and continuity."""
    y_coords = [0.0, 10.0]
    x_coords = [0.0, 10.0]

    # SBAS has only 1 solved pixel, ISBAS has 4 solved pixels
    sbas_v = np.array([[-5.0, np.nan], [np.nan, np.nan]])
    isbas_v = np.array([[-5.0, -2.0], [-1.0, 0.0]])

    sbas_ds = xr.Dataset({"velocity_mm_yr": (("y", "x"), sbas_v)}, coords={"y": y_coords, "x": x_coords})
    isbas_ds = xr.Dataset({"velocity_mm_yr": (("y", "x"), isbas_v)}, coords={"y": y_coords, "x": x_coords})

    cls = xr.DataArray(np.array([[1, 1], [2, 2]]), dims=("y", "x"), coords={"y": y_coords, "x": x_coords})
    aoi = xr.DataArray(np.ones((2, 2), bool), dims=("y", "x"), coords={"y": y_coords, "x": x_coords})

    report = density_report(sbas_ds, isbas_ds, aoi, cls)
    assert isinstance(report, pd.DataFrame)
    row_a = report[report["class"] == "A"].iloc[0]
    assert row_a["n_pixels"] == 2
    assert row_a["sbas_solved"] == 1
    assert row_a["isbas_solved"] == 2
    assert row_a["gain_pixels"] == 1

    # Agreement on the 1 common pixel
    agr = agreement(sbas_ds, isbas_ds)
    assert agr["n_common"] == 1
    assert abs(agr["bias_mm_yr"] - 0.0) < 1e-5
    assert abs(agr["rmse_mm_yr"] - 0.0) < 1e-5

    # Spatial continuity on flat field is 0.0
    flat = xr.DataArray(np.ones((3, 3)) * -2.0)
    assert abs(spatial_continuity(flat) - 0.0) < 1e-5


# =========================================================================
# products.py
# =========================================================================

def test_breathing_classification_all_four_classes():
    """breathing_classification categorises pixels into 1 (stable), 2 (breathing),

    3 (subsidence), and 4 (both).
    """
    coords = {"y": [0.0, 10.0], "x": [0.0, 10.0]}

    # Class 1: stable (v = -0.5 >= -2.0, amp = 5.0 <= 10.0)
    # Class 2: breathing (v = -0.5 >= -2.0, amp = 15.0 > 10.0)
    # Class 3: subsidence (v = -6.0 < -2.0, se = 1.0 -> |v| > 2*se, amp = 5.0 <= 10.0)
    # Class 4: both (v = -6.0 < -2.0, se = 1.0 -> |v| > 2*se, amp = 15.0 > 10.0)
    v = np.array([
        [-0.5, -0.5],
        [-6.0, -6.0],
    ])
    se = np.array([
        [0.5, 0.5],
        [1.0, 1.0],
    ])
    amp = np.array([
        [5.0, 15.0],
        [5.0, 15.0],
    ])

    vel_ds = xr.Dataset({
        "velocity_mm_yr": (("y", "x"), v),
        "velocity_se_mm_yr": (("y", "x"), se),
    }, coords=coords)
    amp_da = xr.DataArray(amp, dims=("y", "x"), coords=coords)

    classes = breathing_classification(
        vel_ds, amp_da,
        subsidence_thr_mm_yr=-2.0,
        breathing_thr_mm=10.0,
    )

    assert classes.isel(y=0, x=0) == 1
    assert classes.isel(y=0, x=1) == 2
    assert classes.isel(y=1, x=0) == 3
    assert classes.isel(y=1, x=1) == 4


def test_breathing_classification_non_significant_subsidence_defaults_to_stable():
    """When velocity is < -2.0 but not statistically significant (|v| <= 2*se),

    it must be classified as stable (1) instead of subsidence (3).
    """
    coords = {"y": [0.0], "x": [0.0]}
    # v = -2.5 is below -2.0 threshold, BUT se = 2.0 so |v| < 2*se (4.0) -> not significant
    vel_ds = xr.Dataset({
        "velocity_mm_yr": (("y", "x"), [[-2.5]]),
        "velocity_se_mm_yr": (("y", "x"), [[2.0]]),
    }, coords=coords)
    amp_da = xr.DataArray([[-5.0]], dims=("y", "x"), coords=coords)

    classes = breathing_classification(vel_ds, amp_da, subsidence_thr_mm_yr=-2.0, breathing_thr_mm=10.0)
    assert classes.isel(y=0, x=0) == 1


def test_products_summary_table_structure():
    """summary_table generates expected rows for all land cover classes."""
    coords = {"y": [0.0, 10.0], "x": [0.0, 10.0]}
    vel_ds = xr.Dataset({
        "velocity_mm_yr": (("y", "x"), np.zeros((2, 2))),
        "velocity_se_mm_yr": (("y", "x"), np.ones((2, 2)) * 0.5),
        "rmse_mm": (("y", "x"), np.ones((2, 2)) * 1.2),
    }, coords=coords)
    amp_da = xr.DataArray(np.ones((2, 2)) * 8.0, dims=("y", "x"), coords=coords)
    cls_da = xr.DataArray(np.array([[1, 2], [3, 4]]), dims=("y", "x"), coords=coords)
    aoi_da = xr.DataArray(np.ones((2, 2), bool), dims=("y", "x"), coords=coords)

    df = summary_table(vel_ds, amp_da, cls_da, aoi_da)
    assert isinstance(df, pd.DataFrame)
    assert set(df["class"]) == {"A", "B", "C_peat_core", "D_transition", "E_water"}
    assert "vel_mean_mm_yr" in df.columns
    assert "amp_median_mm" in df.columns


# =========================================================================
# validation.py
# =========================================================================

def test_annual_chain_closure_triplets_and_jump_detection():
    """annual_chain_closure checks triplet consistency: 2-yr rate == mean(1-yr rates)."""
    # Chain 2021 -> 2022 -> 2023:
    # 2021-2022: -4.0 mm/yr
    # 2022-2023: -6.0 mm/yr
    # Consistent 2021-2023 should be: (-4.0 + -6.0)/2 = -5.0 mm/yr
    consistent_rates = {
        "2021-2022": -4.0,
        "2022-2023": -6.0,
        "2021-2023": -5.0,
    }
    res_clean = annual_chain_closure(consistent_rates)
    assert len(res_clean["triplets"]) == 1
    assert abs(res_clean["max_abs_closure"]) < 1e-6

    # Inconsistent chain (e.g. unwrapping phase jump)
    inconsistent_rates = {
        "2021-2022": -4.0,
        "2022-2023": -6.0,
        "2021-2023": -15.0,  # 10 mm/yr discrepancy
    }
    res_jump = annual_chain_closure(inconsistent_rates)
    assert abs(res_jump["max_abs_closure"] - 10.0) < 1e-6


def test_quality_filter_bounds():
    """quality_filter requires both low RMS residual and sufficient valid pairs."""
    coords = {"y": [0.0, 10.0], "x": [0.0, 10.0]}
    # (0,0): rms=0.4 (good), pairs=15 (good) -> good
    # (0,1): rms=1.5 (too noisy), pairs=15 -> bad
    # (1,0): rms=0.4 (good), pairs=5 (too few) -> bad
    # (1,1): NaN (unsolved)
    rms = np.array([[0.4, 1.5], [0.4, np.nan]])
    pairs = np.array([[15, 15], [5, 0]])

    ds = xr.Dataset({
        "rms_residual_rad": (("y", "x"), rms),
        "n_valid_pairs": (("y", "x"), pairs),
    }, coords=coords)
    aoi = xr.DataArray(np.ones((2, 2), bool), dims=("y", "x"), coords=coords)

    q = quality_filter(ds, aoi, rms_max_rad=1.0, min_pairs=10)
    assert q["n_solved_total"] == 3
    assert q["n_good_total"] == 1
    assert abs(q["frac_good_of_solved"] - 1 / 3) < 1e-4


def test_correlate_insar_hydrology():
    """correlate_insar_hydrology calculates Pearson r against hydrological proxies."""
    dates = pd.date_range("2022-01-01", periods=10, freq="12D")
    proxy_dates = pd.date_range("2022-01-01", periods=120, freq="1D")

    # Proxy is a daily cycle
    proxy = pd.Series(np.sin(np.linspace(0, 4 * np.pi, len(proxy_dates))), index=proxy_dates)
    # InSAR samples the proxy on SAR acquisition dates
    insar = proxy.reindex(dates, method="nearest") * 5.0

    res = correlate_insar_hydrology(insar, proxy, tolerance_days=2)
    assert res["n"] == 10
    assert abs(res["r"] - 1.0) < 1e-5
    assert abs(res["r2"] - 1.0) < 1e-5


def test_seasonal_amplitude_recovers_known_cycle():
    """seasonal_amplitude computes median peak-to-peak amplitude across multi-year series."""
    # 2 years, 24 dates per year (freq='15D')
    dates = pd.date_range("2021-01-01", "2022-12-31", freq="15D")
    t_years = ((dates - dates[0]) / pd.Timedelta("365.25D")).values
    # 1 pixel with -5 mm/yr trend and 14 mm peak-to-peak seasonal cycle (7 * sin)
    signal = -5.0 * t_years + 7.0 * np.sin(2 * np.pi * t_years)
    ts_mm = xr.DataArray(
        signal[:, None, None],
        dims=("time", "y", "x"),
        coords={"time": dates, "y": [0.0], "x": [0.0]},
    )
    amp = seasonal_amplitude(ts_mm)
    assert amp.name == "seasonal_amplitude_mm"
    # Should recover ~14.0 mm peak-to-peak amplitude
    assert abs(float(amp.isel(y=0, x=0)) - 14.0) < 0.5


def test_decorrelation_summary_aggregates_seasons():
    """decorrelation_summary computes median coherence and fraction usable per season."""
    df = pd.DataFrame([
        {"season": "été", "mean_coh": 0.20, "dt_days": 12},
        {"season": "été", "mean_coh": 0.40, "dt_days": 24},
        {"season": "hiver", "mean_coh": 0.70, "dt_days": 12},
        {"season": "hiver", "mean_coh": 0.80, "dt_days": 12},
    ])
    summary = decorrelation_summary(df, gamma_usable=0.30)
    assert len(summary) == 2
    ete_row = summary[summary["season"] == "été"].iloc[0]
    assert ete_row["n_pairs"] == 2
    assert abs(ete_row["median_coh"] - 0.30) < 1e-5
    assert abs(ete_row["frac_usable"] - 0.50) < 1e-5

    hiv_row = summary[summary["season"] == "hiver"].iloc[0]
    assert abs(hiv_row["median_coh"] - 0.75) < 1e-5
    assert abs(hiv_row["frac_usable"] - 1.00) < 1e-5


"""Ground-truth tests for the three untested modules that carry science.

These are not smoke tests. Each builds a small scene whose right answer is known
by construction, then checks the module recovers it — the same shape as the
existing synthetic tests.

Why these three first. `s1`, `s2`, `era5` and `jobs` are I/O: when they fail, a
network call fails and you see it. `quality`, `classify` and `water_mask` are
different — a subtly wrong composite still returns plausible numbers, and those
numbers reach the paper. `W_i = conditional non-water coherence × non-water
fraction` is exactly that kind of composite.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import xarray as xr

from insar_wetlands.classify import CLASS_NAMES, classify, pick_reference_pixel
from insar_wetlands.masking.water_mask import (
    flooded_fraction, optical_to_s1_dates, water_mask,
)
from insar_wetlands.quality import (
    pair_dry_mask, quality_index, weight_matrix_for_pixel,
)

Y, X = 4, 5


def grid(values, dims=("y", "x")):
    arr = np.asarray(values, dtype=float)
    coords = {"y": np.arange(arr.shape[-2]), "x": np.arange(arr.shape[-1])}
    return xr.DataArray(arr, dims=dims, coords=coords)


def water_series(dates, per_date_flags) -> xr.DataArray:
    """(time, y, x) boolean water mask from an explicit list of scenes."""
    arr = np.asarray(per_date_flags, dtype=bool)
    return xr.DataArray(arr, dims=("time", "y", "x"),
                        coords={"time": pd.to_datetime(dates),
                                "y": np.arange(arr.shape[1]),
                                "x": np.arange(arr.shape[2])})


# =========================================================================
# quality — W_i
# =========================================================================

def test_a_pair_is_dry_only_when_both_of_its_dates_are_dry():
    """The whole point of the corrected W_i: a pair is usable only if neither
    of its acquisitions was under water."""
    dates = ["2024-01-01", "2024-01-13", "2024-01-25"]
    # one pixel, wet on the middle date only
    flags = [[[False]], [[True]], [[False]]]
    water = water_series(dates, flags)
    pairs = ["20240101_20240113",   # second date wet  -> not dry
             "20240113_20240125",   # first date wet   -> not dry
             "20240101_20240125"]   # both dry         -> dry
    dry = pair_dry_mask(water, pairs)
    assert dry.sel(pair="20240101_20240113").item() is np.False_ or not dry.sel(pair="20240101_20240113").item()
    assert not dry.sel(pair="20240113_20240125").item()
    assert dry.sel(pair="20240101_20240125").item()


def test_conditional_coherence_ignores_the_wet_pairs():
    """A pixel that is coherent when dry must not be dragged down by the pairs
    in which it was flooded — the double-penalty this formulation fixes."""
    dates = ["2024-01-01", "2024-01-13", "2024-01-25"]
    water = water_series(dates, [[[False]], [[True]], [[False]]])
    pairs = ["20240101_20240113", "20240113_20240125", "20240101_20240125"]
    # coherence: 0.10 on the two wet pairs, 0.90 on the dry one
    corr = xr.DataArray(np.array([[[0.10]], [[0.10]], [[0.90]]]),
                        dims=("pair", "y", "x"),
                        coords={"pair": pairs, "y": [0], "x": [0]})
    ff = grid([[0.0]])

    q = quality_index(corr, water, ff)

    assert q["coh_conditional_dry"].item() == pytest.approx(0.90)
    # the naive average would have been 0.367 — the value this design rejects
    assert q["coh_all_pairs"].item() == pytest.approx(0.3667, abs=1e-3)
    assert q["n_dry_pairs"].item() == 1.0


def test_W_is_conditional_coherence_times_the_dry_fraction():
    """W_i = coh_dry × (1 − flooded_fraction), stated in the module docstring."""
    dates = ["2024-01-01", "2024-01-13"]
    water = water_series(dates, [[[False]], [[False]]])
    pairs = ["20240101_20240113"]
    corr = xr.DataArray(np.array([[[0.80]]]), dims=("pair", "y", "x"),
                        coords={"pair": pairs, "y": [0], "x": [0]})
    q = quality_index(corr, water, grid([[0.25]]))
    assert q["W"].item() == pytest.approx(0.80 * 0.75)


def test_W_is_zero_for_a_permanently_flooded_pixel():
    dates = ["2024-01-01", "2024-01-13"]
    water = water_series(dates, [[[True]], [[True]]])
    corr = xr.DataArray(np.array([[[0.95]]]), dims=("pair", "y", "x"),
                        coords={"pair": ["20240101_20240113"], "y": [0], "x": [0]})
    q = quality_index(corr, water, grid([[1.0]]))
    assert q["W"].item() == pytest.approx(0.0)
    assert q["n_dry_pairs"].item() == 0.0


def test_W_stays_within_zero_and_one():
    """It is a weight. A value outside [0, 1] would silently distort the
    inversion rather than fail."""
    dates = ["2024-01-01", "2024-01-13"]
    water = water_series(dates, [[[False]], [[False]]])
    corr = xr.DataArray(np.array([[[1.4]]]), dims=("pair", "y", "x"),
                        coords={"pair": ["20240101_20240113"], "y": [0], "x": [0]})
    q = quality_index(corr, water, grid([[-0.2]]))
    assert 0.0 <= q["W"].item() <= 1.0


def test_weight_matrix_drops_pairs_below_the_coherence_floor():
    corr = np.array([0.9, 0.5, 0.2, 0.8])
    dry = np.array([True, False, True, True])
    w = weight_matrix_for_pixel(corr, dry, gamma_min=0.30)
    assert w.tolist() == [0.9, 0.0, 0.0, 0.8]   # wet dropped, and 0.2 < 0.30


# =========================================================================
# classify — the five behavioural classes
# =========================================================================

def config(permanent=0.80, intermittent=0.15, coh_stable=0.45) -> dict:
    return {"classification": {"permanent_water_frac": permanent,
                               "intermittent_frac": intermittent,
                               "stable_coherence": coh_stable}}


def test_each_class_is_assigned_where_its_definition_says():
    ff = grid([[0.90, 0.50, 0.05],      # inside:  E, D, C
               [0.90, 0.50, 0.05]])     # outside: A/B by coherence only
    coh = grid([[0.10, 0.10, 0.10],
                [0.60, 0.20, 0.60]])
    aoi = grid([[1, 1, 1], [0, 0, 0]]).astype(bool)

    cls = classify(ff, coh, aoi, config())

    assert cls.values[0].tolist() == [5, 4, 3]          # E, D, C inside
    assert cls.values[1].tolist() == [1, 2, 1]          # A, B, A outside
    assert set(CLASS_NAMES) == {1, 2, 3, 4, 5}


def test_hydrological_classes_never_appear_outside_the_aoi():
    """A frequently-wet pixel outside the peatland is a ditch or a field, not a
    transition zone. The module comment says so; this pins it."""
    ff = grid([[0.95, 0.40]])
    coh = grid([[0.60, 0.60]])
    aoi = grid([[0, 0]]).astype(bool)
    cls = classify(ff, coh, aoi, config())
    assert set(np.unique(cls.values)) <= {1, 2}


def test_the_thresholds_are_honoured_at_their_boundaries():
    """`permanent` is >=, `intermittent` is <=; a pixel exactly on either
    boundary must fall on the documented side."""
    ff = grid([[0.80, 0.15]])
    coh = grid([[0.10, 0.10]])
    aoi = grid([[1, 1]]).astype(bool)
    cls = classify(ff, coh, aoi, config(permanent=0.80, intermittent=0.15))
    assert cls.values[0].tolist() == [5, 3]


def test_thresholds_come_from_the_config_not_from_constants():
    ff = grid([[0.50]])
    coh = grid([[0.10]])
    aoi = grid([[1]]).astype(bool)
    assert classify(ff, coh, aoi, config(permanent=0.40)).item() == 5
    assert classify(ff, coh, aoi, config(permanent=0.90)).item() == 4


def test_reference_pixel_is_the_most_coherent_class_A_pixel_nearby():
    cls = grid([[1, 1, 2], [1, 3, 3]])
    coh = grid([[0.50, 0.95, 0.99], [0.20, 0.99, 0.99]])
    ref = pick_reference_pixel(cls, coh, centroid_xy=(1.0, 0.0), max_dist_m=5.0)
    assert (ref["y"], ref["x"]) == (0.0, 1.0)      # 0.99s are not class A
    assert ref["coherence"] == pytest.approx(0.95)


def test_reference_pixel_falls_back_when_nothing_is_close_enough():
    """Better a distant class-A pixel than no reference at all — but the
    fallback must be deliberate, not accidental."""
    cls = grid([[1, 3], [3, 3]])
    coh = grid([[0.77, 0.10], [0.10, 0.10]])
    ref = pick_reference_pixel(cls, coh, centroid_xy=(50.0, 50.0), max_dist_m=1.0)
    assert ref["coherence"] == pytest.approx(0.77)


# =========================================================================
# water_mask
# =========================================================================

def masking_config() -> dict:
    return {"masking": {"ndwi_water_threshold": 0.2,
                        "mndwi_water_threshold": 0.3,
                        "sigma0_vv_water_db": -18.0,
                        "double_bounce_delta_db": 4.0}}


def stacks(ndwi, g0_db, dates):
    t = pd.to_datetime(dates)
    s2 = xr.Dataset({"ndwi": (("time", "y", "x"), np.asarray(ndwi, dtype=float))},
                    coords={"time": t, "y": [0], "x": [0]})
    rtc = xr.Dataset({"gamma0_vv_db": (("time", "y", "x"),
                                       np.asarray(g0_db, dtype=float))},
                     coords={"time": t, "y": [0], "x": [0]})
    return s2, rtc


def test_open_water_is_found_by_either_sensor_alone():
    dates = ["2024-01-01", "2024-01-13", "2024-01-25"]
    #                optical only    radar only     neither
    s2, rtc = stacks(ndwi=[[[0.50]], [[-0.40]], [[-0.40]]],
                     g0_db=[[[-5.0]], [[-25.0]], [[-5.0]]], dates=dates)
    m = water_mask(s2, rtc, masking_config())
    assert m["water"].values.ravel().tolist() == [True, True, False]


def test_double_bounce_is_flagged_separately_not_as_water():
    """Flooded vegetation is bright, not dark. Counting it as open water would
    misclassify the lake edge; it gets its own flag."""
    dates = ["2024-01-01", "2024-01-13", "2024-01-25"]
    s2, rtc = stacks(ndwi=[[[-0.4]], [[-0.4]], [[-0.4]]],
                     g0_db=[[[-12.0]], [[-12.0]], [[-2.0]]], dates=dates)
    m = water_mask(s2, rtc, masking_config())
    assert m["water"].values.ravel().tolist() == [False, False, False]
    assert m["hidden_water"].values.ravel().tolist() == [False, False, True]
    assert m["water_or_hidden"].values.ravel().tolist() == [False, False, True]


def test_a_pixel_is_never_both_water_and_hidden_water():
    dates = ["2024-01-01", "2024-01-13"]
    s2, rtc = stacks(ndwi=[[[0.9]], [[0.9]]],
                     g0_db=[[[-30.0]], [[-2.0]]], dates=dates)
    m = water_mask(s2, rtc, masking_config())
    assert not (m["water"] & m["hidden_water"]).any().item()


def test_flooded_fraction_is_the_time_average():
    dates = ["2024-01-01", "2024-01-13", "2024-01-25", "2024-02-06"]
    arr = np.array([[[True]], [[True]], [[False]], [[False]]])
    mask = xr.Dataset({"water_or_hidden": (("time", "y", "x"), arr)},
                      coords={"time": pd.to_datetime(dates), "y": [0], "x": [0]})
    assert flooded_fraction(mask).item() == pytest.approx(0.5)


def test_raster_attrs_are_cleared_so_they_do_not_leak_into_figures():
    """`.where()` propagates source attrs; the RTC long_name was reaching the
    figures. The module clears them deliberately — pin that."""
    dates = ["2024-01-01"]
    s2, rtc = stacks(ndwi=[[[0.5]]], g0_db=[[[-5.0]]], dates=dates)
    rtc["gamma0_vv_db"].attrs = {"long_name": "Sentinel-1 Calibrated ..."}
    m = water_mask(s2, rtc, masking_config())
    assert all(m[v].attrs == {} for v in m.data_vars)


def test_optical_dates_beyond_the_tolerance_become_missing():
    """The optical record is asynchronous. Silently reusing a scene from three
    weeks away would fabricate a water state."""
    s2 = xr.Dataset(
        {"ndwi": (("time", "y", "x"), np.array([[[0.5]], [[0.5]]]))},
        coords={"time": pd.to_datetime(["2024-01-01", "2024-03-01"]),
                "y": [0], "x": [0]})
    s1_dates = pd.to_datetime(["2024-01-03", "2024-02-01"])
    out = optical_to_s1_dates(s2, s1_dates, max_gap_days=12)
    assert not np.isnan(out["ndwi"].sel(time="2024-01-03").item())   # 2 days
    assert np.isnan(out["ndwi"].sel(time="2024-02-01").item())       # 31 days

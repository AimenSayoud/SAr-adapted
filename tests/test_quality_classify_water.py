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

from insar_wetlands.masking.water_mask import (
    flooded_fraction,
    optical_to_s1_dates,
    water_mask,
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


# =========================================================================
# classify — the five behavioural classes
# =========================================================================

def config(permanent=0.80, intermittent=0.15, coh_stable=0.45) -> dict:
    return {"classification": {"permanent_water_frac": permanent,
                               "intermittent_frac": intermittent,
                               "stable_coherence": coh_stable}}


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

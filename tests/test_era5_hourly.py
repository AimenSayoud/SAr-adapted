"""X-043: daily precipitation must sum all 24 hourly accumulations."""
import warnings

import numpy as np
import pandas as pd
import pytest
import xarray as xr

from insar_wetlands.hydro import daily_era5_point, era5_samples_per_day, open_era5


def _era5(hours_step: int, days: int = 3, tp_mm_per_h: float = 1.0) -> xr.Dataset:
    t = pd.date_range("2022-01-01", periods=days * 24 // hours_step, freq=f"{hours_step}h")
    shape = (len(t), 2, 2)
    return xr.Dataset(
        {"tp": (("valid_time", "latitude", "longitude"), np.full(shape, tp_mm_per_h / 1000.0, "float32")),
         "t2m": (("valid_time", "latitude", "longitude"), np.full(shape, 273.15 + 5, "float32")),
         "tcwv": (("valid_time", "latitude", "longitude"), np.full(shape, 10.0, "float32"))},
        coords={"valid_time": t, "latitude": [52.8, 52.7], "longitude": [16.2, 16.3]})


def test_hourly_daily_sum_is_24_hours_of_accumulation():
    daily = daily_era5_point(_era5(1), 16.2, 52.8)
    assert daily.precip_mm.iloc[0] == pytest.approx(24.0)


def test_six_hourly_daily_sum_is_the_x043_bias():
    daily = daily_era5_point(_era5(6), 16.2, 52.8)
    assert daily.precip_mm.iloc[0] == pytest.approx(4.0)   # 4 of 24 h — why the loader exists


def test_samples_per_day():
    assert era5_samples_per_day(_era5(1)) == pytest.approx(24)
    assert era5_samples_per_day(_era5(6)) == pytest.approx(4)


def test_open_era5_prefers_the_hourly_download(tmp_path):
    h = tmp_path / "era5_hourly"
    h.mkdir()
    ds = _era5(1)
    ds[["tp"]].to_netcdf(h / "era5_2022_tp_hourly.nc")
    ds[["t2m", "tcwv"]].to_netcdf(h / "era5_2022_inst_hourly.nc")
    _era5(6).to_netcdf(tmp_path / "era5_rzecin.nc")
    out = open_era5(tmp_path)
    assert era5_samples_per_day(out) == pytest.approx(24)
    assert {"tp", "t2m", "tcwv"} <= set(out.data_vars)


def test_open_era5_warns_on_the_legacy_six_hourly_file(tmp_path):
    _era5(6).to_netcdf(tmp_path / "era5_rzecin.nc")
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        out = open_era5(tmp_path)
    assert era5_samples_per_day(out) == pytest.approx(4)
    assert any("X-043" in str(x.message) for x in w)

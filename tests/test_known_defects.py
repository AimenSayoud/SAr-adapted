"""Known defects found by the web-atlas audit (2026-09-24), pinned as strict xfails.

Each test states the *correct* behaviour. It fails today (xfail); the moment the
defect is fixed it passes, and `strict=True` turns that into a failure so the
marker must be removed in the same change — the fix cannot land silently, and the
paper number it moves (see each ticket) gets revisited with it."""

import numpy as np
import pytest
import xarray as xr

from insar_wetlands.predict_failure import threshold_sweep, zone_fraction_above


def test_threshold_sweep_is_not_pre_rounded():
    """C-039 (fixed 2026-09-25): threshold_sweep used to round to 4 dp before the export
    rounded to 3 dp, so Zone D 0.231477 was published as 0.232."""
    # 1 000 000 px of which 231 477 are >= 0.7: exact fraction 0.231477
    n, k = 1_000_000, 231_477
    field = xr.DataArray(np.r_[np.full(k, 0.9), np.full(n - k, 0.1)].reshape(1000, 1000), dims=("y", "x"))
    zones = {"D": xr.DataArray(np.ones((1000, 1000), bool), dims=("y", "x"))}
    sweep = threshold_sweep(field, zones, zone_names=("D",))
    at07 = float(sweep.loc[np.isclose(sweep.threshold, 0.7), "frac_above"].iloc[0])
    exact = zone_fraction_above(field, zones, "D", 0.7)
    assert exact == pytest.approx(0.231477)
    assert round(at07, 3) == round(exact, 3)      # a single rounding, at export


def test_era5_request_covers_every_hour(monkeypatch, tmp_path):
    """X-043 (fixed 2026-09-25): ERA5 used to be requested at 00/06/12/18 only."""
    import sys
    import types

    from insar_wetlands.acquisition import era5
    captured = {}

    class FakeClient:
        def retrieve(self, name, request, target):
            captured.update(request)
    monkeypatch.setitem(sys.modules, "cdsapi", types.SimpleNamespace(Client=FakeClient))
    monkeypatch.setattr(era5, "_ensure_cdsapirc", lambda: None)
    monkeypatch.setattr(era5, "_unwrap_if_zipped", lambda p: None)
    cfg = {"era5": {"dataset": "reanalysis-era5-single-levels", "area": [53, 15, 52, 17]}}
    era5._retrieve(cfg, tmp_path / "x.nc", 2022, ["total_precipitation"])
    assert len(captured["time"]) == 24

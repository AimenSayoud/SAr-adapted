"""Ground-truth tests for MintPy SBAS integration (inversion/sbas_mintpy.py).

Verifies that:
(1) write_config generates deterministic MintPy configuration across reference modes
    (explicit ref_yx, geographic lat/lon, auto), tropospheric delay, and thresholds;
(2) reset_derived_products cleans stale outputs while strictly preserving rzecin.cfg;
(3) has_connected_component detects both ISCE burst and GAMMA conncomp naming conventions;
(4) _first_dataset and best_reference_yx select the highest coherence pixel strictly inside
    the connected component mask and fail closed on empty masks;
(5) read_reference correctly transforms reference pixel indices to UTM coordinates;
(6) load_timeseries accurately converts meters to millimetres and applies temporal coherence masking.
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock

import numpy as np
import pandas as pd
import pytest
import xarray as xr

from insar_wetlands.inversion.sbas_mintpy import (
    _first_dataset,
    _grid_coords,
    best_reference_yx,
    has_connected_component,
    load_temporal_coherence,
    load_timeseries,
    read_reference,
    reset_derived_products,
    write_config,
)

# =========================================================================
# write_config
# =========================================================================

def test_write_config_deterministic_yx_mode(tmp_path: Path):
    """When ref_yx is given, it takes highest priority and sets lalo to auto."""
    cfg_file = write_config(
        work_dir=tmp_path,
        data_dir=tmp_path / "data",
        first_pair="20220101_20220113",
        ref_yx=(42, 84),
        ref_lat=52.76,  # Should be ignored because ref_yx is provided
        ref_lon=16.31,
    )
    assert cfg_file.name == "rzecin.cfg"
    assert cfg_file.exists()
    content = cfg_file.read_text()

    assert "mintpy.reference.yx           = 42,84" in content
    assert "mintpy.reference.lalo         = auto" in content
    assert "mintpy.troposphericDelay.method = no" in content


def test_write_config_geographic_lat_lon_mode(tmp_path: Path):
    """When ref_yx is None and lat/lon are given, lalo is set and yx is auto."""
    cfg_file = write_config(
        work_dir=tmp_path,
        data_dir=tmp_path / "data",
        first_pair="20220101_20220113",
        ref_yx=None,
        ref_lat=52.7632,
        ref_lon=16.3098,
        tropo=True,
    )
    content = cfg_file.read_text()
    assert "mintpy.reference.lalo        = 52.7632,16.3098" in content
    assert "mintpy.reference.yx          = auto" in content
    assert "mintpy.troposphericDelay.method = pyaps" in content


def test_write_config_auto_mode_sets_min_coherence(tmp_path: Path):
    """When no reference is given, auto mode sets minCoherence threshold."""
    cfg_file = write_config(
        work_dir=tmp_path,
        data_dir=tmp_path / "data",
        first_pair="20220101_20220113",
        ref_yx=None,
        ref_lat=None,
        ref_lon=None,
        reference_min_coherence=0.65,
        temp_base_max=36,
        coh_threshold=0.35,
        network_min_coherence=0.28,
        unwrap_error_method="bridging",
    )
    content = cfg_file.read_text()
    assert "mintpy.reference.minCoherence = 0.65" in content
    assert "mintpy.reference.yx           = auto" in content
    assert "mintpy.reference.lalo         = auto" in content
    assert "mintpy.network.tempBaseMax   = 36" in content
    assert "mintpy.networkInversion.maskThreshold = 0.35" in content
    assert "mintpy.network.minCoherence   = 0.28" in content
    assert "mintpy.unwrapError.method    = bridging" in content


# =========================================================================
# reset_derived_products
# =========================================================================

def test_reset_derived_products_preserves_rzecin_cfg(tmp_path: Path):
    """Derived products (*.h5, *.png, smallbaselineApp.cfg) must be deleted,

    while the user template rzecin.cfg must be strictly preserved.
    """
    rzecin = tmp_path / "rzecin.cfg"
    rzecin.write_text("mintpy.load.processor = hyp3")

    sb_cfg = tmp_path / "smallbaselineApp.cfg"
    sb_cfg.write_text("mintpy.reference.yx = 1,1")

    h5_file = tmp_path / "timeseries.h5"
    h5_file.write_text("fake h5")

    mask_file = tmp_path / "maskConnComp.h5"
    mask_file.write_text("fake mask")

    png_file = tmp_path / "velocity.png"
    png_file.write_text("fake png")

    txt_file = tmp_path / "log.txt"
    txt_file.write_text("fake log")

    other_file = tmp_path / "keep_me.json"
    other_file.write_text("{}")

    removed = reset_derived_products(tmp_path)

    assert rzecin.exists(), "rzecin.cfg must NEVER be removed"
    assert other_file.exists(), "Files outside specified patterns must not be removed"

    assert not sb_cfg.exists()
    assert not h5_file.exists()
    assert not mask_file.exists()
    assert not png_file.exists()
    assert not txt_file.exists()

    assert set(removed) == {"smallbaselineApp.cfg", "timeseries.h5", "maskConnComp.h5", "velocity.png", "log.txt"}


# =========================================================================
# has_connected_component
# =========================================================================

def test_has_connected_component_detects_naming_variants(tmp_path: Path):
    """Detects both ISCE burst (_conncomp.tif) and GAMMA (_conn_comp.tif)."""
    assert not has_connected_component(tmp_path)

    p1 = tmp_path / "pair1"
    p1.mkdir()
    (p1 / "pair1_unw_phase.tif").touch()
    assert not has_connected_component(tmp_path)

    # ISCE naming
    (p1 / "pair1_conncomp.tif").touch()
    assert has_connected_component(tmp_path)

    # Clean up and test GAMMA naming
    (p1 / "pair1_conncomp.tif").unlink()
    p2 = tmp_path / "pair2"
    p2.mkdir()
    (p2 / "pair2_conn_comp.tif").touch()
    assert has_connected_component(tmp_path)


# =========================================================================
# _first_dataset
# =========================================================================

def test_first_dataset_finds_2d_dataset():
    class DummyDS:
        def __init__(self, ndim):
            self.ndim = ndim

    class DummyH5(dict):
        filename = "test.h5"

    h5 = DummyH5({
        "metadata": DummyDS(1),
        "spatial_map": DummyDS(2),
        "cube": DummyDS(3),
    })
    assert _first_dataset(h5) == "spatial_map"

    h5_no_2d = DummyH5({"1d": DummyDS(1)})
    with pytest.raises(KeyError, match="aucun dataset 2D"):
        _first_dataset(h5_no_2d)


# =========================================================================
# best_reference_yx
# =========================================================================

def test_best_reference_yx_picks_highest_coherence_in_valid_mask(tmp_path: Path, monkeypatch):
    """A pixel with 0.95 coherence outside the mask must be rejected in favour

    of a 0.75 pixel inside the mask.
    """
    coh_data = np.array([
        [0.95, 0.40],
        [0.30, 0.75],
    ], dtype=float)

    # Mask: (0,0) is 0 (outside conncomp), (1,1) is 1 (inside conncomp)
    mask_data = np.array([
        [0, 0],
        [0, 1],
    ], dtype=int)

    class MockH5File:
        def __init__(self, path):
            self.path = Path(path)
            self.filename = self.path.name

        def __enter__(self):
            return self

        def __exit__(self, *args):
            pass

        def keys(self):
            return ["data"]

        def __getitem__(self, key):
            ds = MagicMock()
            ds.ndim = 2
            if "avgSpatialCoh" in self.path.name:
                ds.__getitem__.side_effect = lambda s: coh_data[s]
                ds.shape = coh_data.shape
            else:
                ds.__getitem__.side_effect = lambda s: mask_data[s]
                ds.shape = mask_data.shape
            return ds

    monkeypatch.setattr("h5py.File", MockH5File)

    (tmp_path / "avgSpatialCoh.h5").touch()
    (tmp_path / "maskConnComp.h5").touch()

    ref = best_reference_yx(tmp_path)
    assert ref["row"] == 1
    assert ref["col"] == 1
    assert abs(ref["avg_spatial_coherence"] - 0.75) < 1e-6
    assert ref["n_valid_pixels_in_mask"] == 1


def test_best_reference_yx_fails_closed_on_empty_mask(tmp_path: Path, monkeypatch):
    """When maskConnComp contains zero valid pixels (>0), it must fail closed."""
    coh_data = np.ones((2, 2), dtype=float) * 0.5
    mask_data = np.zeros((2, 2), dtype=int)  # 0 valid pixels

    class MockH5File:
        def __init__(self, path):
            self.path = Path(path)
            self.filename = self.path.name

        def __enter__(self):
            return self

        def __exit__(self, *args):
            pass

        def keys(self):
            return ["data"]

        def __getitem__(self, key):
            ds = MagicMock()
            ds.ndim = 2
            ds.__getitem__.side_effect = lambda s: (
                coh_data[s] if "avgSpatialCoh" in self.path.name else mask_data[s]
            )
            ds.dtype = mask_data.dtype
            ds.min.return_value = 0
            ds.max.return_value = 0
            return ds

    monkeypatch.setattr("h5py.File", MockH5File)
    (tmp_path / "avgSpatialCoh.h5").touch()
    (tmp_path / "maskConnComp.h5").touch()

    with pytest.raises(RuntimeError, match="ne contient AUCUN pixel > 0"):
        best_reference_yx(tmp_path)


# =========================================================================
# read_reference
# =========================================================================

def test_read_reference_computes_geo_coordinates(tmp_path: Path, monkeypatch):
    """read_reference parses REF_* and X/Y step attributes to compute UTM coordinates."""
    fake_attrs = {
        "REF_X": 15,
        "REF_Y": 8,
        "X_FIRST": 300000.0,
        "X_STEP": 40.0,
        "Y_FIRST": 5800000.0,
        "Y_STEP": -40.0,
        "REF_LAT": b"52.7632",
        "REF_LON": b"16.3098",
    }

    class MockH5File:
        def __init__(self, path):
            self.attrs = fake_attrs

        def __enter__(self):
            return self

        def __exit__(self, *args):
            pass

    monkeypatch.setattr("h5py.File", MockH5File)
    (tmp_path / "timeseries.h5").touch()

    ref = read_reference(tmp_path)
    assert ref["row"] == 8
    assert ref["col"] == 15
    # x = 300000.0 + 40.0 * 15 = 300600.0
    assert abs(ref["x"] - 300600.0) < 1e-4
    # y = 5800000.0 + (-40.0) * 8 = 5799680.0
    assert abs(ref["y"] - 5799680.0) < 1e-4
    assert abs(ref["lat"] - 52.7632) < 1e-4
    assert abs(ref["lon"] - 16.3098) < 1e-4


# =========================================================================
# load_timeseries & load_temporal_coherence
# =========================================================================

def test_load_timeseries_converts_m_to_mm_and_masks_coherence(tmp_path: Path, monkeypatch):
    """Timeseries is loaded in millimetres and masked where temporal coherence < threshold."""
    dates_bytes = [b"20220101", b"20220113"]
    # 2 dates, 2x2 grid, in meters: 0.010 m = 10 mm
    ts_meters = np.array([
        [[0.010, 0.020], [0.030, 0.040]],
        [[0.015, 0.025], [0.035, 0.045]],
    ], dtype=float)

    # Coherence: (0,0)=0.8 (valid), (0,1)=0.5 (below 0.6 threshold), (1,0)=0.9, (1,1)=0.4
    tcoh_data = np.array([
        [0.8, 0.5],
        [0.9, 0.4],
    ], dtype=float)

    attrs = {
        "X_FIRST": 100.0,
        "X_STEP": 10.0,
        "Y_FIRST": 200.0,
        "Y_STEP": -10.0,
    }

    class MockH5File:
        def __init__(self, path):
            self.path = Path(path)
            self.attrs = attrs

        def __enter__(self):
            return self

        def __exit__(self, *args):
            pass

        def __getitem__(self, key):
            if key == "timeseries":
                return ts_meters
            elif key == "date":
                return dates_bytes
            elif key == "temporalCoherence":
                return tcoh_data
            raise KeyError(key)

    monkeypatch.setattr("h5py.File", MockH5File)
    (tmp_path / "timeseries.h5").touch()
    (tmp_path / "temporalCoherence.h5").touch()

    # Load with threshold = 0.6
    da = load_timeseries(tmp_path, coh_threshold=0.6)

    assert isinstance(da, xr.DataArray)
    assert da.dims == ("time", "y", "x")
    assert da.attrs["units"] == "mm"
    assert len(da.time) == 2
    assert pd.Timestamp("2022-01-01") in pd.to_datetime(da.time.values)

    # Check mm conversion on valid pixel (0,0): 0.010 m -> 10.0 mm
    assert abs(float(da.isel(time=0, y=0, x=0)) - 10.0) < 1e-5
    # Pixel (0,1) had coherence 0.5 < 0.6 -> must be NaN
    assert np.isnan(float(da.isel(time=0, y=0, x=1)))
    # Pixel (1,0) had coherence 0.9 >= 0.6 -> 30.0 mm
    assert abs(float(da.isel(time=0, y=1, x=0)) - 30.0) < 1e-5
    # Pixel (1,1) had coherence 0.4 < 0.6 -> NaN
    assert np.isnan(float(da.isel(time=0, y=1, x=1)))


def test_grid_coords_constructs_utm_arrays():
    """_grid_coords builds linear coordinate grids from MintPy attributes."""
    attrs = {
        "X_FIRST": 300000.0,
        "X_STEP": 40.0,
        "Y_FIRST": 5800000.0,
        "Y_STEP": -40.0,
    }
    coords = _grid_coords(attrs, ny=3, nx=2)
    np.testing.assert_allclose(coords["x"], [300000.0, 300040.0])
    np.testing.assert_allclose(coords["y"], [5800000.0, 5799960.0, 5799920.0])


def test_load_temporal_coherence_structure(tmp_path: Path, monkeypatch):
    """load_temporal_coherence returns a 2D DataArray with correct name and coords."""
    tcoh_data = np.array([[0.85, 0.45], [0.90, 0.70]], dtype=float)
    attrs = {"X_FIRST": 0.0, "X_STEP": 10.0, "Y_FIRST": 100.0, "Y_STEP": -10.0}

    class MockH5File:
        def __init__(self, path):
            self.attrs = attrs

        def __enter__(self):
            return self

        def __exit__(self, *args):
            pass

        def __getitem__(self, key):
            if key == "temporalCoherence":
                return tcoh_data
            raise KeyError(key)

    monkeypatch.setattr("h5py.File", MockH5File)
    (tmp_path / "temporalCoherence.h5").touch()

    da = load_temporal_coherence(tmp_path)
    assert da.name == "temporal_coherence"
    assert da.dims == ("y", "x")
    assert da.shape == (2, 2)
    assert abs(float(da.isel(y=0, x=0)) - 0.85) < 1e-5


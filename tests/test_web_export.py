"""Ground-truth tests for the web-atlas export.

The per-pixel seasonal fit is a scientific routine, so it gets the same test as
everything else here: recover a known answer from simulated data, and agree with
the routine it vectorises (`aggregate.seasonal_amplitude`) pixel for pixel."""

import json

import numpy as np
import pandas as pd
import pytest

from insar_wetlands.aggregate import seasonal_amplitude
from insar_wetlands.web_export import (
    AtlasWriter,
    amplitude_dispersion_from_db,
    dequantize,
    grid_corners_lonlat,
    grid_edges,
    labels_to_geojson,
    pair_midpoint_season,
    quantize,
    seasonal_fit_stack,
    seasonal_mean_maps,
)

# The real atlas grid: 138 x 129 pixels of 40 m, centres as in every Drive stack.
X = 585740.0 + 40.0 * np.arange(138)
Y = 5849260.0 - 40.0 * np.arange(129)


def _dates(n=90, start="2022-01-08", step=12):
    return pd.date_range(start, periods=n, freq=f"{step}D")


def test_seasonal_fit_recovers_known_amplitude_phase_and_trend():
    d = _dates()
    t = (d - d[0]).days.values / 365.25
    rng = np.random.default_rng(0)
    amp = np.array([[3.3, 10.0], [0.5, 25.0]])
    phi = np.array([[0.3, 1.0], [2.0, -1.2]])
    trend = np.array([[-1.0, 0.0], [2.5, -20.0]])
    y = (amp[None] * np.cos(2 * np.pi * t[:, None, None] - phi[None])
         + trend[None] * t[:, None, None]
         + rng.normal(0, 0.05, (len(t), 2, 2)))
    out = seasonal_fit_stack(y, d)
    np.testing.assert_allclose(out["amplitude_mm"], amp, atol=0.03)
    np.testing.assert_allclose(out["trend_mm_yr"], trend, atol=0.03)
    assert np.all(out["r2_seasonal"] > 0.95)
    # the maximum of cos(2πt - φ) falls φ/2π of a year after the first date
    doy = (d[0].dayofyear + (phi % (2 * np.pi)) / (2 * np.pi) * 365.25) % 365.25
    dp = np.abs(out["phase_doy"] - doy)
    assert np.all(np.minimum(dp, 365.25 - dp) < 2.0)


def test_seasonal_fit_matches_seasonal_amplitude_pixel_for_pixel():
    """Same definitions as the 1-D routine the project uses, including NaN gaps
    and the day-of-year convention, so a map value can be quoted as that fit."""
    d = _dates()
    rng = np.random.default_rng(1)
    y = rng.normal(0, 2.0, (len(d), 3, 4))
    y[:, 0, 0] += 4 * np.sin(2 * np.pi * (d - d[0]).days.values / 365.25)
    y[rng.random(y.shape) < 0.25] = np.nan          # per-pixel gaps
    out = seasonal_fit_stack(y, d)
    for i in range(3):
        for j in range(4):
            ref = seasonal_amplitude(pd.DataFrame({"date": d, "disp_mm": y[:, i, j]}))
            assert out["amplitude_mm"][i, j] == pytest.approx(ref["amplitude_mm"], abs=1e-3)
            assert out["trend_mm_yr"][i, j] == pytest.approx(ref["trend_mm_yr"], abs=1e-3)
            assert out["r2_seasonal"][i, j] == pytest.approx(ref["r2_seasonal"], abs=1e-4)
            dp = abs(out["phase_doy"][i, j] - ref["phase_doy"])
            assert min(dp, 365.25 - dp) < 0.1


def test_seasonal_fit_honours_shared_epoch():
    d = _dates()
    y = np.sin(2 * np.pi * (d - pd.Timestamp("2022-01-01")).days.values / 365.25)[:, None, None]
    a = seasonal_fit_stack(y, d, epoch="2022-01-01")
    ref = seasonal_amplitude(pd.DataFrame({"date": d, "disp_mm": y[:, 0, 0]}),
                             epoch=pd.Timestamp("2022-01-01"))
    assert a["phase_doy"][0, 0] == pytest.approx(ref["phase_doy"], abs=0.1)


def test_seasonal_fit_masks_sparse_pixels():
    d = _dates(20)
    y = np.ones((20, 1, 2))
    y[:17, 0, 1] = np.nan
    out = seasonal_fit_stack(y, d, min_obs=6)
    assert np.isnan(out["amplitude_mm"][0, 1])
    assert out["n"][0, 1] == 3


def test_amplitude_dispersion_constant_and_known_case():
    const = np.full((10, 2, 2), -12.0)
    assert np.allclose(amplitude_dispersion_from_db(const), 0.0)
    # amplitudes 1 and 3 alternating -> mean 2, std 1 -> D_A = 0.5
    amps = np.array([1.0, 3.0] * 5)
    db = 10 * np.log10(amps ** 2)[:, None, None]
    assert amplitude_dispersion_from_db(db)[0, 0] == pytest.approx(0.5)


def test_quantize_roundtrip_within_half_step_and_nodata():
    rng = np.random.default_rng(2)
    a = rng.normal(0, 10, (50, 50))
    a[3, 3] = np.nan
    for dtype, tol_steps in (("int16", 0.51), ("uint8", 0.51)):
        q, s, o = quantize(a, dtype=dtype)
        back = dequantize(q, s, o)
        assert np.isnan(back[3, 3])
        fin = np.isfinite(a)
        assert np.max(np.abs(back[fin] - a[fin])) <= tol_steps * s


def test_grid_edges_and_corners_are_pixel_edges_not_centres():
    w, s, e, n = grid_edges(X, Y)
    assert (w, s, e, n) == (585720.0, 5844120.0, 591240.0, 5849280.0)  # the GIS tif bounds
    c = grid_corners_lonlat(X, Y)
    assert len(c) == 4
    (lon_tl, lat_tl), (lon_tr, _), (_, lat_br), _ = c
    assert lon_tr > lon_tl and lat_tl > lat_br
    # the site centroid from config.yaml lies inside the grid
    assert lon_tl < 16.3098 < lon_tr and lat_br < 52.7632 < lat_tl


def test_labels_to_geojson_counts_match_raster():
    lab = np.zeros((129, 138), int)
    lab[10:20, 10:30] = 1
    lab[50:52, 60:61] = 2
    fc = labels_to_geojson(lab, X, Y, {1: {"zone": "A"}, 2: {"zone": "B"}})
    counts = {f["properties"]["zone"]: f["properties"]["n_px"] for f in fc["features"]}
    assert counts == {"A": 200, "B": 2}
    # 200 px of 40 m = 32 ha; area check in a local projection
    from pyproj import Transformer
    from shapely.geometry import shape
    from shapely.ops import transform
    tr = Transformer.from_crs("EPSG:4326", "EPSG:32633", always_xy=True)
    area = transform(tr.transform, shape(fc["features"][0]["geometry"])).area
    assert area == pytest.approx(200 * 1600, rel=1e-6)


def test_seasonal_mean_maps_filters_long_baselines():
    pairs = ["20220105_20220117", "20220110_20220420", "20220710_20220722"]
    stack = np.stack([np.full((2, 2), v) for v in (0.5, 0.1, 0.3)])
    out = seasonal_mean_maps(stack, pairs, max_dt_days=48)
    assert pair_midpoint_season(pairs) == ["DJF", "MAM", "JJA"]
    assert out["DJF"][0, 0] == 0.5 and out["JJA"][0, 0] == 0.3
    assert "MAM" not in out                    # the 100-day pair was dropped


def test_atlas_writer_manifest_and_grid_guard(tmp_path):
    w = AtlasWriter(tmp_path, X, Y)
    f = np.random.default_rng(3).random((129, 138))
    f[0, 0] = np.nan
    w.raster("coh", f, title="t", group="g", status="core")
    w.raster("stk", np.stack([f, f]), title="t", group="g", status="derived",
             times=["2022-01-08", "2022-01-20"])
    with pytest.raises(ValueError, match="atlas grid"):
        w.raster("shifted", f, title="t", group="g", status="core", x=X - 20, y=Y + 20)
    with pytest.raises(ValueError, match="status"):
        w.raster("bad", f, title="t", group="g", status="nonsense")
    w.write_manifest(title="x")
    m = json.loads((tmp_path / "manifest.json").read_text())
    lyr = {entry["id"]: entry for entry in m["layers"]}
    q = np.fromfile(tmp_path / lyr["coh"]["file"], dtype="<i2").reshape(129, 138)
    back = dequantize(q, lyr["coh"]["scale"], lyr["coh"]["offset"])
    assert np.isnan(back[0, 0]) and np.nanmax(np.abs(back - f)) < 1e-4
    assert lyr["stk"]["shape"] == [2, 129, 138]
    assert m["grid"]["width"] == 138 and m["grid"]["dx"] == 40.0

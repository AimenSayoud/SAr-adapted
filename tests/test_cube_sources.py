"""Parsers, indices and contracts of the data cube recover known answers (no network)."""
import numpy as np
import pandas as pd
import pytest

from insar_wetlands.cube import build, gold, schema, sources


def test_s2_indices_mask_and_formulae():
    b = {k: np.full((2, 2), v, "float32") for k, v in
         {"B03": 0.1, "B04": 0.05, "B05": 0.1, "B08": 0.4, "B8A": 0.4, "B11": 0.2}.items()}
    b["SCL"] = np.array([[4, 9], [11, 6]], "float32")            # vegetation, cloud, snow, water
    ix = sources.s2_indices(b)
    assert ix["ndvi"][0, 0] == pytest.approx((0.4 - 0.05) / 0.45)
    assert ix["ndmi"][0, 0] == pytest.approx(0.2 / 0.6)
    assert np.isnan(ix["ndvi"][0, 1]) and np.isnan(ix["ndvi"][1, 0])     # cloud and snow are not surface values
    assert ix["snow"][1, 0] == 1 and ix["valid"][1, 1] == 1


def test_block_mean():
    a = np.arange(16.0).reshape(4, 4)
    assert sources.block_mean(a, 2)[0, 0] == pytest.approx(np.mean([0, 1, 4, 5]))


def test_landsat_qa_clear():
    qa = np.array([21824, 22280, 1, 30048], float)    # clear, cloud (bit 3), fill (bit 0), snow (bit 5)
    assert sources.landsat_qa_clear(qa).tolist() == [True, False, False, False]


def test_parse_open_meteo():
    j = {"hourly": {"time": ["2024-07-01T00:00", "2024-07-01T01:00"], "temperature_2m": [18.2, None]}}
    d = sources.parse_open_meteo(j, "era5_land_")
    assert d.index.tz is not None and d.era5_land_temperature_2m.iloc[0] == 18.2 and np.isnan(d.iloc[1, 0])


def test_parse_imgw_status_codes():
    vals = ["353160330", "POZNAŃ", "2022", "01", "02"] + [""] * (len(sources.IMGW_COLUMNS) - 5)
    row = dict(zip(sources.IMGW_COLUMNS, vals))
    row.update(STD="1,5", PKSN="3", WPKSN="", ROSA="0", WROSA="9", SZRO="", WSZRO="8", SGR="Z")
    line = ",".join(f'"{row[c]}"' for c in sources.IMGW_COLUMNS) + "\n"
    d = sources.parse_imgw_daily(line.encode("cp1250"))
    r = d.iloc[0]
    assert r.station == 330 and r.date == pd.Timestamp("2022-01-02")
    assert r.tmean_c == 1.5 and r.snow_depth_cm == 3 and r.dew_h == 0 and np.isnan(r.hoarfrost_h)
    assert r.ground_state == "Z"


def test_lidar_features_flat_and_vegetation():
    dtm = np.full((80, 80), 50.0)
    dsm = dtm.copy()
    dsm[:40, :40] += 2.0                                           # a 2 m canopy over one 40 m cell
    f = build.lidar_features(dtm, dsm, 40)
    assert f["lidar_vegh_mean_m"][0, 0] == pytest.approx(2.0) and f["lidar_vegh_mean_m"][1, 1] == 0
    assert f["lidar_roughness_m"][0, 0] == 0 and f["lidar_valid_frac"].min() == 1


def test_nearest_pixel_skips_invalid_scenes():
    src_t = pd.DatetimeIndex(pd.to_datetime(["2022-06-01", "2022-06-03", "2022-06-20"], utc=True))
    vals = np.array([[1.0, 1.0], [2.0, np.nan], [3.0, 3.0]])
    valid = np.array([[1, 1], [1, 0], [1, 1]], float)
    v, lag = gold.nearest_pixel(vals, valid, src_t, pd.DatetimeIndex(pd.to_datetime(["2022-06-04"], utc=True)), 10)
    assert v[0].tolist() == [2.0, 1.0] and lag[0].tolist() == [-1.0, -3.0]


def test_schema_catches_duplicates_and_ranges():
    df = pd.DataFrame({"date": ["2022-01-01", "2022-01-01"], "surface_cm": [1.0, 500.0], "n_ok_hours": [3, 4]})
    p = schema.validate(df, "truth_laser_daily")
    assert any("duplicated" in x for x in p) and any("surface_cm" in x for x in p)
    assert schema.validate(df.iloc[:1], "truth_laser_daily") == []


def test_paste_and_to20():
    out = build._paste(np.ones((3, 3)), (4, 4), -1, 2)
    assert np.isfinite(out).sum() == 4 and np.isfinite(out[0, 2]) and np.isnan(out[3, 3])
    tx, ty = np.array([20.0, 60.0]), np.array([60.0, 20.0])
    x20, y20 = np.array([30.0, 50.0]), np.array([50.0, 30.0])    # a 20 m grid inset by half a 40 m cell
    a = build._to20(np.array([[1, 2], [3, 4]]), tx, ty, x20, y20)
    assert a.tolist() == [[1, 2], [3, 4]]


def test_asset_scale_offset():
    class A:
        extra_fields = {"raster:bands": [{"scale": 0.0001, "offset": -0.1}]}
    sc, off = sources.asset_scale_offset(A())
    dn = np.array([1000.0 + 4000, 1000.0 + 500])          # NIR 0.40, red 0.05 after the baseline-04 offset
    nir, red = dn * sc + off
    assert (nir - red) / (nir + red) == pytest.approx(0.35 / 0.45)

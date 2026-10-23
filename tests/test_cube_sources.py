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


def test_indices_refuse_non_positive_reflectance():
    b = {k: np.full((1, 2), v, "float32") for k, v in
         {"B03": 0.05, "B04": -0.01, "B05": 0.02, "B08": 0.03, "B8A": 0.03, "B11": 0.01}.items()}
    b["SCL"] = np.array([[6, 6]], "float32")
    ix = sources.s2_indices(b)
    assert np.isnan(ix["ndvi"]).all() and np.isfinite(ix["ndwi"]).all()
    assert np.nanmax(np.abs(ix["ndmi"])) <= 1


# ------------------------------------------------------------------------------ GUGiK dated sheets (C-072)

def test_wfs_bbox_is_northing_first():
    assert sources.wfs_bbox_2180((317000, 546000, 320000, 549000)) == \
        "546000,317000,549000,320000,urn:ogc:def:crs:EPSG::2180"


def test_parse_gugik_index_and_keep():
    gml = ('<wfs:member><gugik:SkorowidzNMT2024 gml:id="a"><gugik:godlo>N-33-117-D-a-3-4</gugik:godlo>'
           '<gugik:zrodlo_danych>Skaning laserowy</gugik:zrodlo_danych><gugik:format>ARC/INFO ASCII GRID</gugik:format>'
           '<gugik:char_przestrz>1.00 m</gugik:char_przestrz><gugik:uklad_xy>PL-1992</gugik:uklad_xy>'
           '<gugik:url_do_pobrania>https://x/a.asc</gugik:url_do_pobrania></gugik:SkorowidzNMT2024></wfs:member>'
           '<wfs:member><gugik:SkorowidzNMT2024 gml:id="b"><gugik:godlo>N-33-117-D-a-3</gugik:godlo>'
           '<gugik:zrodlo_danych>Zdj. lotnicze</gugik:zrodlo_danych><gugik:format>ARC/INFO ASCII GRID</gugik:format>'
           '<gugik:char_przestrz>5.00 m</gugik:char_przestrz><gugik:uklad_xy>PL-1992</gugik:uklad_xy>'
           '</gugik:SkorowidzNMT2024></wfs:member>')
    recs = sources.parse_gugik_index(gml, "SkorowidzNMT2024")
    assert [r["godlo"] for r in recs] == ["N-33-117-D-a-3-4", "N-33-117-D-a-3"]
    assert [sources.keep_gugik_sheet("dtm", r) for r in recs] == [True, False]      # the photogrammetric 5 m is out
    assert sources.keep_gugik_sheet("cir", {"uklad_xy": "PL-1992", "kolor": "CIR"})
    assert not sources.keep_gugik_sheet("cir", {"uklad_xy": "PL-1992", "kolor": "RGB"})
    assert not sources.keep_gugik_sheet("cir", {"uklad_xy": "PL-2000:S6", "kolor": "CIR"})
    assert sources.keep_gugik_sheet("rgb", {"uklad_xy": "PL-1992", "kolor": "RGB", "piksel": "0.1"})
    assert not sources.keep_gugik_sheet("rgb", {"uklad_xy": "PL-1992", "kolor": "RGB", "piksel": "0.25"})


def test_parse_orto_dates():
    html = ('var s=[{url:"https://o/1_N-33.tif",godlo:"N-33", aktualnosc:"2021-09-29", kolor:"CIR"},'
            '{url:"https://o/2_N-33.tif",godlo:"N-33", aktualnosc:"2025-03-20", kolor:"RGB", dt_pzgik:"2025-05-01"}];')
    assert sources.parse_orto_dates(html) == {"https://o/1_N-33.tif": "2021-09-29", "https://o/2_N-33.tif": "2025-03-20"}


def test_canopy_features_known_shares():
    chm = np.zeros((80, 40), "float32")
    chm[:40][:10, :] = 10.0                     # top cell: a quarter of its 1 m cells under 10 m canopy
    chm[40:] = np.nan
    chm[40:60, :] = 1.0                         # bottom cell: half valid, all of it 1 m high (not tree)
    f = build.canopy_features(chm, 40, "canopy2024")
    assert f["canopy2024_cover"][:, 0].tolist() == pytest.approx([0.25, 0.0])
    assert f["canopy2024_p90_m"][0, 0] == pytest.approx(10.0)
    assert f["canopy2024_p90_m"][1, 0] == pytest.approx(1.0)
    assert f["canopy2024_valid_frac"][:, 0].tolist() == pytest.approx([1.0, 0.5])
    assert np.isnan(build.canopy_features(np.full((40, 40), np.nan), 40, "c")["c_cover"][0, 0])


def test_cir_ndvi():
    v = build.cir_ndvi(np.array([200, 50, 0], "uint8"), np.array([50, 50, 0], "uint8"))
    assert v[0] == pytest.approx(150 / 250) and v[1] == 0 and np.isnan(v[2])


def test_sheets_to_utm_places_a_square_where_pyproj_says(tmp_path):
    from pyproj import Transformer
    from rasterio.transform import from_origin
    e0, n0 = 318500.0, 547500.0                 # a 10 m square of 1s on a PL-1992 grid of 0s
    a = np.zeros((200, 200), "float32")
    a[95:105, 95:105] = 1.0
    tr = from_origin(e0 - 100, n0 + 100, 1, 1)
    x, y = Transformer.from_crs(2180, 32633, always_xy=True).transform(e0, n0)
    b = (np.floor(x) - 60, np.floor(y) - 60, np.floor(x) + 60, np.floor(y) + 60)
    out = build._sheets_to_utm(["s"], b, lambda f: (a, tr))
    ys, xs = np.nonzero(out > 0.5)
    cx, cy = b[0] + xs.mean() + 0.5, b[3] - ys.mean() - 0.5
    assert abs(cx - x) < 1.0 and abs(cy - y) < 1.0
    assert 80 <= (out > 0.5).sum() <= 120       # the square keeps its area (~100 m²) through the rotation
    two = build._sheets_to_utm(["s", "t"], b, lambda f: (a * (2 if f == "t" else 1), tr))
    assert np.nanmax(two) == pytest.approx(1.5, abs=0.01)    # overlapping sheets are averaged

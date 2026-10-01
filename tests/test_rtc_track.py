"""masking/rtc.py — the track-aware RTC stack (descending, 2020–2024).

Known answers on a fake catalogue, no network: one item per day chosen by AOI cover, the day's
next item used when the first fails, power → dB and the VH−VV ratio exact, a rerun downloads
nothing, a file built with other settings is refused, and the coverage table marks every date.
"""
from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
import rioxarray  # noqa: F401 — registers .rio
import xarray as xr

from insar_wetlands.masking.rtc import (
    build_rtc_track_stack,
    dualpol_db,
    items_by_day,
    rtc_coverage,
    track_relative_orbit,
    widen_bbox,
)

BBOX = (16.0, 52.0, 16.1, 52.1)


def poly(x0, y0, x1, y1) -> dict:
    return {"type": "Polygon", "coordinates": [[[x0, y0], [x1, y0], [x1, y1], [x0, y1], [x0, y0]]]}


def item(iid, when, geom=None, pols=("vv", "vh")):
    return SimpleNamespace(id=iid, datetime=pd.Timestamp(when, tz="UTC").to_pydatetime(),
                           geometry=geom or poly(15.0, 51.0, 17.0, 53.0),
                           assets={p: SimpleNamespace(href=f"{iid}/{p}") for p in pols})


@pytest.fixture
def template():
    t = xr.DataArray(np.zeros((2, 3), "float32"), dims=("y", "x"),
                     coords={"y": [5849000.0, 5848960.0], "x": [585740.0, 585780.0, 585820.0]})
    return t.rio.write_crs("EPSG:32633")


def fake_loader(template, power: dict, fail=()):
    """href → constant linear power on the template grid; hrefs in ``fail`` raise; calls logged."""
    calls = []

    def load(href):
        calls.append(href)
        if href in fail:
            raise OSError("network")
        v = power[href]
        if np.ndim(v):
            return template.copy(data=np.asarray(v, "float32"))
        return xr.full_like(template, v, dtype="float32")

    return load, calls


def test_track_relative_orbit_reads_the_tracks_block():
    cfg = {"sentinel1": {"relative_orbit": 175,
                         "tracks": {"ascending": {"relative_orbit": 175},
                                    "descending": {"relative_orbit": 22}}}}
    assert track_relative_orbit(cfg, "descending") == 22


def test_items_by_day_groups_orders_by_cover_and_drops_single_pol():
    half = poly(16.05, 51.0, 17.0, 53.0)                     # covers half of BBOX
    items = [item("b_full", "2022-03-01T05:09"), item("a_half", "2022-03-01T05:08", half),
             item("c_full", "2022-03-01T05:09"), item("vv_only", "2022-03-13T05:09", pols=("vv",)),
             item("d", "2022-03-25T05:09")]
    days = items_by_day(items, BBOX)
    assert list(days) == [pd.Timestamp("2022-03-01"), pd.Timestamp("2022-03-25")]
    assert [it.id for it in days[pd.Timestamp("2022-03-01")]] == ["b_full", "c_full", "a_half"]


def test_dualpol_db_is_ten_log10_and_ratio_is_vh_minus_vv():
    vv = xr.DataArray(np.array([1.0, 0.1, 0.0]))
    vh = xr.DataArray(np.array([0.01, 0.01, 0.01]))
    ds = dualpol_db(vv, vh)
    np.testing.assert_allclose(ds.gamma0_vv_db.values[:2], [0.0, -10.0], atol=1e-9)
    np.testing.assert_allclose(ds.ratio_vh_vv_db.values[:2], [-20.0, -10.0], atol=1e-9)
    assert np.isnan(ds.gamma0_vv_db.values[2])               # zero power is not -inf


def test_build_uses_next_item_when_first_fails_and_rerun_downloads_nothing(template, tmp_path):
    items = [item("d1", "2022-01-10T05:09"),
             item("d2_a", "2022-01-22T05:09"), item("d2_b", "2022-01-22T05:09")]
    power = {"d1/vv": 0.1, "d1/vh": 0.01, "d2_a/vv": 1.0, "d2_a/vh": 1.0, "d2_b/vv": 0.01, "d2_b/vh": 0.001}
    load, calls = fake_loader(template, power, fail={"d2_a/vv"})
    out = tmp_path / "rtc_dualpol_stack_descending.nc"
    build_rtc_track_stack(items, template, out, track="descending", relative_orbit=22, load=load,
                          checkpoint_every=1)

    ds = xr.load_dataset(out)
    assert list(ds.source_item.values) == ["d1", "d2_b"]
    np.testing.assert_allclose(ds.gamma0_vv_db.values[:, 0, 0], [-10.0, -20.0], atol=1e-5)
    np.testing.assert_allclose(ds.ratio_vh_vv_db.values[:, 0, 0], [-10.0, -10.0], atol=1e-5)
    assert ds.attrs["track"] == "descending" and int(ds.attrs["relative_orbit"]) == 22
    assert ds.attrs["resampling"] == "nearest"

    calls.clear()
    build_rtc_track_stack(items, template, out, track="descending", relative_orbit=22, load=load)
    assert calls == []
    assert xr.load_dataset(out).time.size == 2


def test_build_refuses_to_extend_a_file_made_with_other_settings(template, tmp_path):
    items = [item("d1", "2022-01-10T05:09")]
    load, _ = fake_loader(template, {"d1/vv": 0.1, "d1/vh": 0.01})
    out = tmp_path / "stack.nc"
    build_rtc_track_stack(items, template, out, track="descending", relative_orbit=22, load=load)
    with pytest.raises(ValueError, match="resampling"):
        build_rtc_track_stack(items, template, out, track="descending", relative_orbit=22,
                              resampling="average", load=load)
    with pytest.raises(ValueError, match="track"):
        build_rtc_track_stack(items, template, out, track="ascending", relative_orbit=175, load=load)


NORTH = [[0.1, 0.1, 0.1], [np.nan, np.nan, np.nan]]       # a slice whose data stops inside the AOI
SOUTH = [[np.nan, np.nan, np.nan], [0.01, 0.01, 0.01]]


def test_a_day_is_the_mosaic_of_its_slices_best_first(template, tmp_path):
    """Spring 2020 descending: the northern slice claims the AOI but covers 63 % of it; the
    southern one (footprint missing the AOI) fills the rest — and only the rest."""
    items = [item("north", "2020-03-03T05:07", poly(15.0, 51.0, 17.0, 53.0)),
             item("south", "2020-03-03T05:08", poly(15.0, 50.0, 17.0, 51.5))]
    power = {"north/vv": NORTH, "north/vh": NORTH, "south/vv": SOUTH, "south/vh": SOUTH}
    load, _ = fake_loader(template, power)
    out = tmp_path / "s.nc"
    build_rtc_track_stack(items, template, out, track="descending", relative_orbit=22, bbox=BBOX, load=load)
    ds = xr.load_dataset(out)
    assert ds.source_item.values.tolist() == ["north+south"]
    np.testing.assert_allclose(ds.gamma0_vv_db.values[0], [[-10.0] * 3, [-20.0] * 3], atol=1e-5)


def test_a_day_left_incomplete_by_a_failed_slice_is_retried_not_written(template, tmp_path):
    items = [item("north", "2020-03-03T05:07"), item("south", "2020-03-03T05:08")]
    power = {"north/vv": NORTH, "north/vh": NORTH, "south/vv": SOUTH, "south/vh": SOUTH}
    load, _ = fake_loader(template, power, fail={"south/vv"})
    out = tmp_path / "s.nc"
    build_rtc_track_stack(items, template, out, track="descending", relative_orbit=22, load=load)
    assert not out.exists()
    load, _ = fake_loader(template, power)
    build_rtc_track_stack(items, template, out, track="descending", relative_orbit=22, load=load)
    assert bool(xr.load_dataset(out).gamma0_vv_db.notnull().all())


def test_a_gap_no_slice_covers_is_written_and_an_empty_day_is_not(template, tmp_path):
    items = [item("north", "2020-03-03T05:07"), item("void", "2020-03-15T05:07")]
    empty = np.full((2, 3), np.nan)
    load, _ = fake_loader(template, {"north/vv": NORTH, "north/vh": NORTH, "void/vv": empty, "void/vh": empty})
    out = tmp_path / "s.nc"
    build_rtc_track_stack(items, template, out, track="descending", relative_orbit=22, load=load)
    ds = xr.load_dataset(out)
    assert ds.time.size == 1 and int(ds.gamma0_vv_db.isnull().sum()) == 3


def test_source_item_is_not_cut_when_a_longer_id_is_appended(template, tmp_path):
    """Read back, the ids are fixed-width text; a later mosaic's '+'-joined ids must survive."""
    out = tmp_path / "s.nc"
    one = [item("d1", "2020-01-09T05:09")]
    load, _ = fake_loader(template, {"d1/vv": 0.1, "d1/vh": 0.01})
    build_rtc_track_stack(one, template, out, track="descending", relative_orbit=22, load=load)
    two = [item("north", "2020-03-03T05:07"), item("south", "2020-03-03T05:08")]
    load, _ = fake_loader(template, {"north/vv": NORTH, "north/vh": NORTH, "south/vv": SOUTH, "south/vh": SOUTH})
    build_rtc_track_stack(one + two, template, out, track="descending", relative_orbit=22, load=load)
    assert xr.load_dataset(out).source_item.values.tolist() == ["d1", "north+south"]


def test_widen_bbox_grows_every_side():
    assert widen_bbox((16.0, 52.0, 16.1, 52.1), 0.25) == (15.75, 51.75, 16.35, 52.35)


def test_rtc_coverage_marks_every_date():
    cov = rtc_coverage(["2022-01-10", "2022-01-22"], ["2022-01-10"], ["2022-01-10", "2022-02-03"])
    assert cov.date.dt.strftime("%m-%d").tolist() == ["01-10", "01-22", "02-03"]
    assert cov.in_catalogue.tolist() == [True, True, False]
    assert cov.in_stack.tolist() == [True, False, False]
    assert cov.in_interferograms.tolist() == [True, False, True]


def test_zone_backscatter_series_medians_rvi_and_valid_count():
    """Known powers per zone → exact dB medians and RVI = 4·VH/(VV+VH); a half-empty zone date is NaN."""
    from insar_wetlands.stratify import zone_backscatter_series

    vv = np.array([[[0.1, 0.1, 0.01]], [[0.1, np.nan, np.nan]]])        # 2 dates, 1×3 grid
    vh = vv / 10
    db = {k: 10 * np.log10(v) for k, v in (("vv", vv), ("vh", vh))}
    ds = xr.Dataset({"gamma0_vv_db": (("time", "y", "x"), db["vv"]),
                     "gamma0_vh_db": (("time", "y", "x"), db["vh"]),
                     "ratio_vh_vv_db": (("time", "y", "x"), db["vh"] - db["vv"])},
                    coords={"time": pd.to_datetime(["2020-01-01", "2020-01-13"])})
    zones = {"A": xr.DataArray([[True, True, False]]), "B": xr.DataArray([[False, False, True]]),
             "C": xr.DataArray([[False, False, False]]), "D": xr.DataArray([[True, True, True]])}
    s = zone_backscatter_series(ds, zones).set_index(["date", "zone"])
    a1 = s.loc[(pd.Timestamp("2020-01-01"), "A")]
    assert a1.vv_db == pytest.approx(-10.0) and a1.ratio_vh_vv_db == pytest.approx(-10.0)
    assert a1.rvi == pytest.approx(4 * 0.01 / 0.11)
    assert s.loc[(pd.Timestamp("2020-01-01"), "B")].vv_db == pytest.approx(-20.0)
    assert s.loc[(pd.Timestamp("2020-01-13"), "A")].n_valid == 1          # 1 of 2: still ≥ 50 %
    assert np.isnan(s.loc[(pd.Timestamp("2020-01-13"), "B")].vv_db)      # 0 of 1
    assert np.isnan(s.loc[(pd.Timestamp("2020-01-13"), "D")].vv_db)      # 1 of 3 < 50 %
    assert "C" not in s.index.get_level_values("zone")

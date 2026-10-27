"""X-138 routines recover answers put in by construction."""
import numpy as np
import pandas as pd

from insar_wetlands import drying


def _chain(times, est, seg):
    return pd.DataFrame({"time": pd.to_datetime(times, utc=True), "est": est, "seg": seg})


def test_season_run_takes_longest_segment_and_rezeroes():
    t = pd.date_range("2019-04-20", periods=30, freq="6D", tz="UTC")
    est = np.arange(30, dtype=float) * -2.0 + 100
    seg = np.r_[np.zeros(12, int), np.ones(18, int)]          # segment 1 is longer inside May–November
    run = drying.season_run(_chain(t, est, seg), 2019, (5, 11))
    assert set(run.seg) == {1}
    assert run.est.iloc[0] == 0.0
    assert np.allclose(np.diff(run.est), -2.0)


def test_season_total_vertical_net_and_drop():
    t = pd.date_range("2020-05-01", periods=5, freq="6D", tz="UTC")
    cos = np.cos(np.radians(40.0))
    los = np.array([0.0, -50.0, -120.0, -80.0, -60.0]) * cos      # vertical −12 cm at the low point, −6 at the end
    out = drying.season_total(_chain(t, los, 0), cos)
    assert np.isclose(out["net_cm"], -6.0) and np.isclose(out["drop_cm"], -12.0)
    assert out["drop_date"] == t[2].date() and out["median_revisit_days"] == 6.0
    assert drying.revisit_class(out["median_revisit_days"]) == "six-day"
    assert drying.revisit_class(12.0) == "twelve-day"


def test_error_band():
    b = drying.error_band([3.0, -4.0, np.nan])
    assert b["n"] == 2 and np.isclose(b["typical_cm"], np.sqrt(12.5)) and b["worst_cm"] == 4.0


def test_water_area_recovers_a_known_lake():
    rng = np.random.default_rng(0)
    lake = np.zeros((20, 20), bool)
    lake[5:10, 5:15] = True                                   # 50 pixels → 2 ha
    water = np.repeat(lake[None], 10, axis=0)
    water[:3] |= rng.random((3, 20, 20)) < 0.2                # flecks in 3 of 10 scenes stay below 0.5
    valid = np.ones_like(water)
    valid[0, 5, 5] = False
    area, freq = drying.water_area(water, valid, np.ones((20, 20), bool))
    assert np.isclose(area, 2.0)
    assert np.isclose(freq[5, 5], 1.0)


def test_trend_and_sentence():
    yrs = np.arange(2017, 2026)
    tr = drying.trend(yrs, -1.0 * (yrs - 2017) + 3)
    assert np.isclose(tr["slope"], -1.0) and tr["p"] < 0.001
    assert drying.trend_sentence(tr, 5.0, "down", "up", "none") == "down"
    assert drying.trend_sentence(tr, 9.0, "down", "up", "none") == "none"   # change of 8 not beyond 9
    flat = drying.trend(yrs, np.r_[1, -1, 1, -1, 1, -1, 1, -1, 1])
    assert drying.trend_sentence(flat, 0.0, "down", "up", "none") == "none"


def test_synthesis():
    assert drying.synthesis(True, True) == "the data show the peatland drying"
    assert drying.synthesis(True, None) == "one line shows a sign of drying"
    assert drying.synthesis(None, False) == "the data do not show a drying trend over 2017–2026"
    assert drying.synthesis(False, True) == "the lines disagree"


def test_extraterrestrial_radiation_fao56_example8():
    # FAO-56 Example 8: 20° S, 3 September (day 246) → Ra = 32.2 MJ m⁻² day⁻¹
    assert abs(float(drying.extraterrestrial_radiation(246, -20.0)) - 32.2) < 0.1


def test_hargreaves_by_hand():
    ra = float(drying.extraterrestrial_radiation(180, 52.76))
    want = 0.0023 * (20 + 17.8) * np.sqrt(10.0) * 0.408 * ra
    assert np.isclose(float(drying.hargreaves_et0(25.0, 15.0, 180, 52.76)), want)
    assert float(drying.hargreaves_et0(5.0, 6.0, 180, 52.76)) == 0.0          # no negative range


def test_station_mean_and_rain_class():
    w = pd.DataFrame({"a": [1.0, np.nan, np.nan], "b": [3.0, 2.0, np.nan]})
    m = drying.station_mean(w)
    assert m.iloc[0] == 2.0 and m.iloc[1] == 2.0 and np.isnan(m.iloc[2])
    assert drying.rain_class(49) == "extremely dry" and drying.rain_class(74.9) == "very dry"
    assert drying.rain_class(89) == "dry" and drying.rain_class(90) == "normal" and drying.rain_class(110) == "normal"
    assert drying.rain_class(111) == "wet" and drying.rain_class(140) == "very wet" and drying.rain_class(160) == "extremely wet"

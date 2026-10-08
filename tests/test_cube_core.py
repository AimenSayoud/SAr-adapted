"""The data cube's building blocks recover known answers on synthetic data."""
import numpy as np
import pandas as pd
import pytest

from insar_wetlands.cube import core


def hourly(values, start="2022-06-01"):
    idx = pd.date_range(start, periods=len(values), freq="1h", tz="UTC")
    return pd.Series(np.asarray(values, float), index=idx)


def test_asof_respects_tolerance_and_reports_lag():
    tab = pd.DataFrame({"v": [1.0, 2.0, 3.0]}, index=pd.to_datetime(["2022-01-01", "2022-01-05", "2022-01-20"], utc=True))
    out = core.asof(["2022-01-04", "2022-01-12"], tab, "3D")
    assert out.v.iloc[0] == 2.0 and out.lag_h.iloc[0] == pytest.approx(24.0)
    assert np.isnan(out.v.iloc[1])                     # nearest is 7–8 days away: missing, not silently far
    back = core.asof(["2022-01-04"], tab, "10D", direction="backward")
    assert back.v.iloc[0] == 1.0 and back.lag_h.iloc[0] == pytest.approx(-72.0)


def test_asof_keeps_query_order():
    tab = pd.DataFrame({"v": np.arange(10.0)}, index=pd.date_range("2022-01-01", periods=10, freq="D", tz="UTC"))
    q = pd.to_datetime(["2022-01-08", "2022-01-02", "2022-01-05"], utc=True)
    assert core.asof(q, tab, "1h").v.tolist() == [7.0, 1.0, 4.0]


def test_value_at_interpolates_and_refuses_gaps():
    s = hourly([0.0, 10.0, np.nan, np.nan, np.nan, 50.0])
    t0 = s.index[0]
    v = core.value_at(s, [t0 + pd.Timedelta(minutes=30), t0 + pd.Timedelta(hours=3), t0])
    assert v[0] == pytest.approx(5.0) and np.isnan(v[1]) and v[2] == 0.0


def test_window_features_means_and_changes():
    s = hourly(np.arange(400.0))                       # +1 per hour
    t = s.index[300]
    f = core.window_features(s, [t], "w")
    assert f.w_at.iloc[0] == 300.0
    assert f.w_mean24h.iloc[0] == pytest.approx(np.mean(np.arange(277, 301)))
    assert f.w_d72h.iloc[0] == pytest.approx(72.0) and f.w_d168h.iloc[0] == pytest.approx(168.0)
    early = core.window_features(s, [s.index[10]], "w")
    assert np.isnan(early.w_mean168h.iloc[0])          # not enough hours behind it


def test_hours_since():
    ev = hourly([0, 1, 0, 0, 0]).astype(bool)
    assert core.hours_since(ev, [ev.index[4], ev.index[0]]).tolist()[0] == pytest.approx(3.0)
    assert np.isnan(core.hours_since(ev, [ev.index[0]])[0])


def test_dew_point_known_values():
    assert core.dew_point_c(20.0, 100.0) == pytest.approx(20.0, abs=1e-6)
    assert core.dew_point_c(20.0, 50.0) == pytest.approx(9.3, abs=0.1)


def test_aggregate_to_grid_mean_sd_fraction():
    tx, ty = np.array([20.0, 60.0]), np.array([60.0, 20.0])          # 2×2 cells of 40 m, y decreasing
    fx = np.arange(5.0, 80, 10)                                       # 8 fine columns of 10 m
    fy = np.arange(75.0, 0, -10)
    v = np.zeros((8, 8))
    v[:, 4:] = 1.0                                                    # right half = 1
    v[0, 0] = np.nan
    a = core.aggregate_to_grid(v, fx, fy, tx, ty)
    assert np.allclose(a["mean"], [[0, 1], [0, 1]])
    assert a["valid_frac"][0, 0] == pytest.approx(15 / 16) and a["valid_frac"][1, 1] == 1.0
    assert a["sd"][0, 0] == 0.0
    fr = core.class_fractions(np.where(v == 1, 50, 10), fx, fy, tx, ty, [10, 50])
    assert fr[50][0, 1] == 1.0 and fr[10][0, 1] == 0.0


def test_window_stats_and_phasor():
    f = np.arange(25.0).reshape(5, 5)
    s = core.window_stats(f, 0, 0)
    assert s["n"] == 4 and s["median"] == pytest.approx(3.0)
    ang, cons = core.phasor_mean(np.array([np.pi - 0.1, -np.pi + 0.1]))
    assert abs(abs(ang) - np.pi) < 1e-9 and cons == pytest.approx(np.cos(0.1))


def test_geometry_and_cycle_risk():
    th = np.radians(90 - 32.26)
    assert core.incidence_deg(th) == pytest.approx(32.26)
    e, n, u = core.los_enu(th, np.radians(10))
    assert e ** 2 + n ** 2 + u ** 2 == pytest.approx(1.0) and u == pytest.approx(np.cos(np.radians(32.26)))
    q = core.WAVELENGTH_MM / 4
    assert core.cycle_risk([q * 0.9, -q * 1.1]).tolist() == [False, True]
    assert core.rad_to_los_mm(4 * np.pi) == pytest.approx(core.WAVELENGTH_MM)


def test_folds():
    t1 = pd.to_datetime(["2021-12-28", "2022-03-01"], utc=True)
    t2 = pd.to_datetime(["2022-01-03", "2022-03-07"], utc=True)
    assert core.pair_year_fold(t1, t2).tolist() == [-1, 2022]
    b = core.spatial_blocks((4, 5), 2)
    assert b[0, 0] == b[1, 1] and b[0, 0] != b[0, 2] and b.max() == 5


def test_overpass_time():
    assert core.overpass_time("2022-06-01", "descending") == pd.Timestamp("2022-06-01 05:09", tz="UTC")


def test_period_of_boundaries_cover_every_period():
    from insar_wetlands.cube.build import CUBE_START, PERIODS
    from insar_wetlands.cube.gold import period_of
    cases = {CUBE_START: "2017_2019", "2019-12-31": "2017_2019", "2020-01-01": "2020_2021", "2021-12-31": "2020_2021",
             "2022-01-01": "2022_2024", "2024-12-31": "2022_2024", "2025-06-01": "2025", "2026-10-05": "2026"}
    for d, want in cases.items():
        assert period_of(d) == want, d
    assert set(cases.values()) == set(PERIODS)

"""Field readers and Sentinel-1 matching, on synthetic fixtures only (field data are not public)."""
import numpy as np
import pandas as pd
import pytest

from insar_wetlands import field


def _hourly_csv(tmp_path, hours=24 * 20):
    t = pd.date_range("2022-01-01", periods=hours, freq="h")
    df = pd.DataFrame({"TIMESTAMP": t.strftime("%Y-%m-%d %H:%M:%S"), "CL": -10.0, "CR": -40.0,
                       "Air_2m": np.where(t.day == 5, -3.0, 5.0)})
    for col in field.WTD_COLUMNS:
        df[col] = -(np.arange(hours) / 24.0)          # falls 1 cm per day
    p = tmp_path / "wtd.csv"
    df.to_csv(p, index=False)
    return p


def test_hourly_reader_centres_stamps_at_minus_half_hour_utc(tmp_path):
    w = field.load_wtd_hourly(_hourly_csv(tmp_path))
    assert str(w.index.tz) == "UTC"
    assert w.index[0] == pd.Timestamp("2021-12-31 23:30", tz="UTC")
    assert {"P1", "P6", "CL_raw", "CR_raw"} <= set(w.columns)


def test_wtd_variables_at_an_acquisition(tmp_path):
    w = field.load_wtd_hourly(_hourly_csv(tmp_path))
    t = pd.Timestamp("2022-01-10 16:36", tz="UTC")
    r = field.wtd_at_times(w, [t], plots=["P6"]).iloc[0]
    days_since_start = (t - w.index[0]).total_seconds() / 86400
    assert r.wtd_at == pytest.approx(-days_since_start, abs=1e-6)
    assert r.change_3d == pytest.approx(-3.0, abs=1e-6)
    assert r.change_7d == pytest.approx(-7.0, abs=1e-6)
    assert r.mean_prev3d > r.wtd_at > r.mean_prev3d - 3     # the preceding window is shallower
    assert r.mean_24h == pytest.approx(r.wtd_at, abs=0.05)


def test_outside_the_record_is_nan(tmp_path):
    w = field.load_wtd_hourly(_hourly_csv(tmp_path))
    r = field.wtd_at_times(w, [pd.Timestamp("2025-01-01", tz="UTC")], plots=["P1"]).iloc[0]
    assert np.isnan(r.wtd_at) and np.isnan(r.change_3d)


def test_snow_mask_covers_the_72_hours_after_frost(tmp_path):
    w = field.load_wtd_hourly(_hourly_csv(tmp_path))
    m = field.snow_mask(w["Air_2m"], w.index, hours=72)
    frost = w.index[w["Air_2m"] < 0]
    assert m[frost].all()
    assert m[frost[-1] + pd.Timedelta(hours=71)]
    assert not m[frost[-1] + pd.Timedelta(hours=80)]


def test_snow_mask_is_true_outside_the_air_record():
    idx = pd.date_range("2024-06-01", periods=48, freq="h", tz="UTC")
    warm = pd.Series(15.0, index=idx)
    later = pd.date_range("2024-06-01", periods=96, freq="h", tz="UTC")
    m = field.snow_mask(warm, later)
    assert not m[idx].any()
    assert m[later[48:]].all()


def test_base_plot_maps_replicates_to_their_wtd_plot():
    assert [field.base_plot(x) for x in ["P5_2", "CL_1", "CR_3", "P6_1", "P8_2023", "P9"]] == \
        ["P5", "P5", "P6", "P6", "P8", "P9"]


def test_window_stats():
    s = field.window_stats(np.array([[1, 2, np.nan], [3, 4, 5]]), valid=np.array([[1, 1, 1], [1, 1, 0]]))
    assert s["n_valid"] == 4 and s["median"] == 2.5 and s["iqr"] == pytest.approx(1.5)
    assert field.window_stats(np.full(4, np.nan))["n_valid"] == 0


def test_field_root_defaults_to_the_hub(monkeypatch):
    monkeypatch.delenv("RZECIN_FIELD_ROOT", raising=False)
    assert field.field_root().parts[-2:] == ("06_data", "field")
    monkeypatch.setenv("RZECIN_FIELD_ROOT", "/x/y")
    assert str(field.field_root()) == "/x/y"


def test_censored_flag_marks_a_sustained_floor_not_a_brief_dip():
    rng = np.random.default_rng(0)
    t = pd.date_range("2022-01-01", periods=24 * 60, freq="h", tz="UTC")
    x = np.r_[np.linspace(-5, -30, 24 * 10), np.full(24 * 20, -30.0), np.linspace(-30, -5, 24 * 10),
              np.full(24 * 18, -5.0), np.full(24 * 2, -30.0)]           # 20-day floor, then a 2-day dip
    s = pd.Series(x + rng.normal(0, 0.3, len(x)), index=t)               # noisy, like P9
    f = field.censored_flag(s)
    assert f[t[24 * 12]:t[24 * 28]].mean() > 0.95                        # the floor
    assert not f[t[24 * 2]:t[24 * 6]].any()                              # falling limb
    assert not f[t[24 * 58]:].any()                                      # brief dip

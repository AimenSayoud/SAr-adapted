"""A coarse grid built from a fine field with a known content shift: the search recovers the shift, and the HyP3
prediction gives half a geocoding cell west and north."""
import numpy as np

from insar_wetlands import placement as P


def _scene(seed=0):
    rng = np.random.default_rng(seed)
    x = np.arange(0.5, 2000, 1.0)
    y = np.arange(1999.5, 0, -1.0)
    X, Y = np.meshgrid(x, y)
    # a lake, a forest strip and smooth relief: sharp edges as in the real references
    water = ((X - 900) ** 2 / 300 ** 2 + (Y - 1100) ** 2 / 200 ** 2 < 1).astype(float)
    trees = ((X > 1400) & (X < 1600) | (Y < 300)).astype(float)
    relief = 5 * np.sin(X / 170) + 3 * np.cos(Y / 230) + rng.normal(0, 0.2, X.shape)
    return x, y, water, trees, relief


def test_box_mean_matches_brute_force():
    x, y, water, *_ = _scene()
    f = P.FineField(water, x, y)
    cx, cy = np.array([900.0, 333.3]), np.array([1100.0, 1500.2])
    m = f.box_mean(cx, cy, 20.0)
    for k in range(2):
        sel = (np.abs(x[None, :] - cx[k]) <= 20) & (np.abs(y[:, None] - cy[k]) <= 20)
        assert np.isclose(m[k], water[sel].mean())


def test_recovers_known_shift_coherence_model():
    x, y, water, trees, _ = _scene()
    fw, ft = P.FineField(water, x, y), P.FineField(trees, x, y)
    cx, cy = np.meshgrid(np.arange(100, 1900, 40.0), np.arange(1900, 100, -40.0))
    cx, cy = cx.ravel(), cy.ravel()
    true = (-15.0, 25.0)
    rng = np.random.default_rng(1)
    coh = 0.6 - 0.45 * fw.box_mean(cx + true[0], cy + true[1], 20) - 0.3 * ft.box_mean(cx + true[0], cy + true[1], 20)
    coh = coh + rng.normal(0, 0.03, coh.shape)

    def score(dx, dy):
        X = np.column_stack([fw.box_mean(cx + dx, cy + dy, 20), ft.box_mean(cx + dx, cy + dy, 20)])
        return P.r2_ols(coh, X)

    dx, dy, _, _ = P.best_shift(score, np.arange(-40, 40.1, 2.5))
    assert abs(dx - true[0]) <= 2.5 and abs(dy - true[1]) <= 2.5


def test_recovers_zero_shift_dem_correlation():
    x, y, _, _, relief = _scene()
    fr = P.FineField(relief, x, y)
    cx, cy = np.meshgrid(np.arange(100, 1900, 40.0), np.arange(1900, 100, -40.0))
    cx, cy = cx.ravel(), cy.ravel()
    dem = fr.box_mean(cx, cy, 20)
    dx, dy, _, _ = P.best_shift(lambda a, b: P.corr(dem, fr.box_mean(cx + a, cy + b, 20)), np.arange(-40, 40.1, 2.5))
    assert abs(dx) <= 2.5 and abs(dy) <= 2.5


def test_hyp3_prediction_half_cell_west_north():
    dx, dy = P.hyp3_isce2_predicted_shift(30.0, 52.6, 52.763)
    assert -15.5 < dx < -14.5 and 24.5 < dy < 25.0       # ~30 m × ~49 m geocoding cell at Rzecin
    dx20, dy20 = P.hyp3_isce2_predicted_shift(20.0, 52.6, 52.763)
    assert np.isclose(dx20 / dx, 2 / 3) and np.isclose(dy20 / dy, 2 / 3)


def test_relocate_is_exact_on_a_linear_field_and_undoes_a_mislabelling():
    x = np.arange(20.0, 4000, 40.0)
    y = np.arange(3980.0, 0, -40.0)
    X, Y = np.meshgrid(x, y)
    truth = 0.01 * X - 0.02 * Y + 3.0                      # the field at each label
    dx, dy = -14.9, 24.7
    mislabelled = 0.01 * (X + dx) - 0.02 * (Y + dy) + 3.0  # content of ground at L + (dx, dy)
    fixed = P.relocate(mislabelled, x, y, dx, dy)
    inner = np.s_[2:-2, 2:-2]
    assert np.allclose(fixed[inner], truth[inner])
    assert np.isnan(fixed[-1]).all()      # the southern row: more than half its source lies outside the grid


def test_relocate_phase_as_phasor_and_nan_neighbours():
    x = np.arange(20.0, 400, 40.0)
    y = np.arange(380.0, 0, -40.0)
    v = np.full((len(y), len(x)), np.pi - 0.05)
    v[:, ::2] = -np.pi + 0.05                               # values either side of ±π: a phasor mean stays near π
    out = P.relocate(v, x, y, -20.0, 0.0, kind="phase")
    assert np.all(np.abs(np.abs(out[:, 1:]) - np.pi) < 0.1)
    w = np.ones((len(y), len(x)))
    w[5, 5] = np.nan
    out = P.relocate(w, x, y, -10.0, 10.0)
    assert np.allclose(out[1:-1, 1:-1][np.isfinite(out[1:-1, 1:-1])], 1.0)


def test_relocate_recovers_the_true_cells_of_a_fine_scene_better_than_the_labels():
    xf, yf, water, trees, _ = _scene()
    f = P.FineField(water + 0.5 * trees, xf, yf)
    x = np.arange(120.0, 1880, 40.0)
    y = np.arange(1880.0, 120, -40.0)
    X, Y = np.meshgrid(x, y)
    dx, dy = -14.9, 24.7
    true = f.box_mean(X.ravel(), Y.ravel(), 20).reshape(X.shape)
    mis = f.box_mean(X.ravel() + dx, Y.ravel() + dy, 20).reshape(X.shape)
    fixed = P.relocate(mis, x, y, dx, dy)
    inner = np.s_[2:-2, 2:-2]
    e_mis = np.sqrt(np.nanmean((mis - true)[inner] ** 2))
    e_fix = np.sqrt(np.nanmean((fixed - true)[inner] ** 2))
    assert e_fix < 0.6 * e_mis


def test_relocate_rtc_undoes_a_known_shift_in_linear_power_and_keeps_ratio_and_rvi_consistent(monkeypatch, tmp_path):
    import xarray as xr
    x = np.arange(20.0, 4000, 40.0)
    y = np.arange(3980.0, 0, -40.0)
    X, Y = np.meshgrid(x, y)
    dx, dy = P.RTC_CONTENT_SHIFT_M["descending"]
    vv_true = 0.05 + 1e-5 * X + 2e-5 * Y                      # linear γ⁰, linear in space: bilinear is exact
    vh_true = 0.01 + 3e-6 * X + 1e-6 * Y
    shifted = lambda f: f(X + dx, Y + dy)  # noqa: E731
    vv = shifted(lambda a, b: 0.05 + 1e-5 * a + 2e-5 * b)
    vh = shifted(lambda a, b: 0.01 + 3e-6 * a + 1e-6 * b)
    db = lambda p: (10 * np.log10(p))[None].astype("float32")  # noqa: E731
    ds = xr.Dataset({"vv_db": (("time", "y", "x"), db(vv)), "vh_db": (("time", "y", "x"), db(vh)),
                     "vh_vv_db": (("time", "y", "x"), db(vh) - db(vv)),
                     "rvi": (("time", "y", "x"), (4 * vh / (vv + vh))[None].astype("float32"))},
                    coords={"time": [0], "x": x, "y": y}, attrs={"track": "descending"})
    out = P.relocate_rtc(ds, "descending")
    inner = np.s_[0, 2:-2, 2:-2]
    assert np.allclose(10 ** (out.vv_db.values[inner] / 10), vv_true[2:-2, 2:-2], rtol=1e-4)
    assert np.allclose(out.vh_vv_db.values[inner], out.vh_db.values[inner] - out.vv_db.values[inner], atol=1e-4)
    assert np.allclose(out.rvi.values[inner], (4 * vh_true / (vv_true + vh_true))[2:-2, 2:-2], rtol=1e-3)
    # the opt-in reader: off by default; on, it reads the placed stack beside the file, refusing one with other dates
    p = tmp_path / "rtc_dualpol_2020_2024_descending.nc"
    ds.to_netcdf(p)
    assert P.placed_name(p.name) == "rtc_dualpol_2020_2024_placed_descending.nc"
    assert P.placed_name("rtc_dualpol_2020_2024.nc") == "rtc_dualpol_2020_2024_placed.nc"
    monkeypatch.delenv(P.RTC_RELOCATE_ENV, raising=False)
    assert np.array_equal(P.open_rtc(p).vv_db.values, ds.vv_db.values)
    out.to_netcdf(tmp_path / P.placed_name(p.name))
    monkeypatch.setenv(P.RTC_RELOCATE_ENV, "1")
    assert np.allclose(P.open_rtc(p).vv_db.values[inner], out.vv_db.values[inner])
    out.isel(time=slice(0, 0)).to_netcdf(tmp_path / P.placed_name("rtc_dualpol_2025_2025_descending.nc"))
    ds.to_netcdf(tmp_path / "rtc_dualpol_2025_2025_descending.nc")
    import pytest
    with pytest.raises(ValueError):
        P.open_rtc(tmp_path / "rtc_dualpol_2025_2025_descending.nc")

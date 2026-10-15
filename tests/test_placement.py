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

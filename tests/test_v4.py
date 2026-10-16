"""v4's re-levelling restores a known level across winters when the chain loses it and the anchor knows it."""
import numpy as np
import pandas as pd
import pytest

from insar_wetlands import v4


def _series():
    t = pd.date_range("2021-05-01", "2023-11-30", freq="6D")
    truth = 10 * np.sin(2 * np.pi * t.dayofyear / 365) + 0.05 * np.arange(len(t))
    return t, truth


def test_relevel_restores_winter_steps():
    t, truth = _series()
    inc = np.diff(truth)
    win = (t[1:].month <= 3) | (t[1:].month == 12)
    inc = np.where(win, 0.0, inc)                       # the chain loses everything in winter
    est = np.r_[0, np.cumsum(inc)]
    c = pd.DataFrame({"time": t, "seg": 0, "est": est})
    anchor = pd.Series(truth, index=t)
    bounds = v4.winter_bounds(c.time)
    assert len(bounds) == 2
    out = v4.relevel(c, lambda x: float(anchor[x]), bounds)
    for a, b in bounds:
        assert np.isclose(out.est[b] - out.est[a], truth[b] - truth[a])
    assert np.abs(out.est.iloc[-1] - out.est.iloc[0] - (truth[-1] - truth[0])) < np.abs(est[-1] - est[0] - (truth[-1] - truth[0]))


def test_relevel_skips_missing_anchor_and_joins_segments():
    t, truth = _series()
    c = pd.DataFrame({"time": t, "seg": (t >= "2022-01-01").astype(int), "est": np.zeros(len(t))})
    bounds = v4.winter_bounds(c.time)
    out = v4.relevel(c, lambda x: np.nan if x.year == 2023 else float(x.dayofyear), bounds)
    a, b = bounds[0]
    assert out.seg[b] == out.seg[a]
    a2, b2 = bounds[1]
    assert np.isclose(out.est[b2] - out.est[a2], 0.0)   # NaN anchor: the chain's own step kept


def test_block_cells_and_mean_recover_known_blocks():
    x = np.arange(20.0, 400, 40.0)                       # cell centres
    y = np.arange(380.0, 0, -40.0)
    r, c = v4.block_cells(x, y, 105.0, 255.0, 2)         # between centres 100/140 (x) and 260/220 (y)
    assert list(c) == [2, 3] and list(r) == [3, 4]
    r3, c3 = v4.block_cells(x, y, 105.0, 255.0, 3)       # inside cell (col 2, row 3)
    assert list(c3) == [1, 2, 3] and list(r3) == [2, 3, 4]
    X, Y = np.meshgrid(x, y)
    field = np.stack([X + Y, 2 * X])                     # two "pairs"
    m = v4.block_mean(field, r3, c3)
    assert np.allclose(m, [X[2:5, 1:4].mean() + Y[2:5, 1:4].mean(), 2 * X[2:5, 1:4].mean()])
    field[0, 2:5, 1:3] = np.nan                          # 6 of 9 missing in the first pair
    assert np.isnan(v4.block_mean(field, r3, c3)[0]) and np.isfinite(v4.block_mean(field, r3, c3)[1])
    with pytest.raises(ValueError):
        v4.block_cells(x, y, 10.0, 255.0, 3)

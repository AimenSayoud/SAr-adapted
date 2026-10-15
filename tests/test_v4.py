"""v4's re-levelling restores a known level across winters when the chain loses it and the anchor knows it."""
import numpy as np
import pandas as pd

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

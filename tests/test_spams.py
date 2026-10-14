"""SPAMS fit recovers known parameters from a synthetic soil surface driven by synthetic weather."""
import numpy as np

from insar_wetlands import spams


def _weather(n=1500, seed=0):
    rng = np.random.default_rng(seed)
    doy = np.arange(n) % 365
    P = rng.gamma(0.4, 6.0, n)                                   # intermittent rain, mm/day
    E = np.clip(2.0 + 2.0 * np.sin(2 * np.pi * (doy - 80) / 365), 0.1, None)   # seasonal ET, mm/day
    return P, E


def test_window_sum():
    x = np.arange(1, 6, dtype=float)
    out = spams.window_sum(x, 3)
    assert np.isnan(out[:2]).all()
    assert np.allclose(out[2:], [6, 9, 12])


def test_fit_recovers_parameters():
    P, E = _weather()
    true = {"x_p": 0.8, "x_e": 1.1, "tau": 60}
    lag = 6
    dh = spams.predict_change(P, E, true, lag)
    noise = np.random.default_rng(1).normal(0, 0.05 * np.nanstd(dh), len(dh))
    train = np.arange(len(P)) < 1000
    par = spams.fit(P, E, dh + noise, train, lag)
    assert par["tau"] == 60
    assert abs(par["x_p"] - 0.8) < 0.05 and abs(par["x_e"] - 1.1) < 0.07
    test = ~train & np.isfinite(dh)
    pred = spams.predict_change(P, E, par, lag)
    assert np.corrcoef(pred[test], dh[test])[0, 1] > 0.99


def test_irreversible_recovered():
    P, E = _weather(seed=2)
    par = {"x_p": 0.8, "x_e": 1.1, "tau": 40}
    lag = 6
    R = spams.reversible(P, E, **par)
    dh = spams.predict_change(P, E, par, lag)
    di = np.full(len(R), np.nan)
    I = spams.irreversible(R, -0.3)
    di[:-lag] = I[lag:] - I[:-lag]
    train = np.isfinite(dh) & np.isfinite(di)
    x_i = spams.fit_irreversible(P, E, par, dh + di, train, lag)
    assert abs(x_i - (-0.3)) < 0.02

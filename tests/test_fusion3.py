"""Fusion v3 building blocks recover known answers on synthetic data."""
import numpy as np
import pandas as pd
import pytest

from insar_wetlands import fusion3 as f3

RNG = np.random.default_rng(3)


def test_reference_removes_a_common_offset():
    P, H, W = 20, 10, 12
    common = RNG.normal(0, 3, P)                       # per-pair constant (HyP3 reference, atmosphere)
    U = common[:, None, None] + RNG.normal(0, 0.1, (P, H, W))
    U[:, 5, 5] += np.arange(P) * 0.5                   # a moving pixel
    ref = np.zeros((H, W), bool)
    ref[:3, :3] = True
    d = U[:, 5, 5] - f3.reference_values(U, ref)
    assert np.std(d - np.arange(P) * 0.5) < 0.15 and abs(np.mean(d - np.arange(P) * 0.5)) < 0.1


def test_nearest_cells():
    m = np.zeros((10, 10), bool)
    m[[0, 9, 5], [0, 9, 6]] = True
    n = f3.nearest_cells(m, 5, 5, 1)
    assert n[5, 6] and n.sum() == 1


def test_year_blocked_residual_uses_other_years():
    x = np.tile(np.linspace(-10, 10, 10), 3)
    years = np.repeat([2021, 2022, 2025], 10)
    y = 0.8 * x
    y[years == 2025] = 2.0 * x[years == 2025]           # a year that disagrees
    res = f3.year_blocked_residual(x, y, years)
    assert np.allclose(res[years == 2025], 1.2 * x[years == 2025])   # its slope is not used to judge it


def test_chain_runs_and_levels():
    t = pd.date_range("2021-06-01", periods=9, freq="6D")
    t1, t2 = t[:-1], t[1:]
    inc = np.ones(8)
    inc[4] = np.nan                                    # breaks the chain after four pairs
    run, (lv,) = f3.chain_runs(t1, t2, inc, min_len=3)
    assert run[:4].tolist() == [0] * 4 and run[4] == -1 and run[5:].tolist() == [1] * 3
    assert np.allclose(lv[:4], [-1.5, -0.5, 0.5, 1.5])


def test_level_metrics():
    t = np.sin(np.linspace(0, 6, 40)) * 10
    m = f3.level_metrics(0.5 * t, t)
    assert m["r_level"] == pytest.approx(1.0) and m["amplitude_ratio"] == pytest.approx(0.5)


def test_noise_model_recovers_looks_and_floor():
    c = np.linspace(0.2, 0.95, 16)
    s = f3.crb_sigma_mm(c, 18.0, 2.5)
    fit = f3.fit_noise(c, s)
    assert fit["looks"] == pytest.approx(18.0, rel=1e-3) and fit["floor_mm"] == pytest.approx(2.5, rel=1e-3)
    assert f3.crb_sigma_mm(0.999, 18.0, 2.5) == pytest.approx(2.5, abs=0.05)


def test_coverage_of_a_correct_model():
    z = RNG.normal(0, 1, 20000)
    c = f3.coverage(z)
    assert c["within_1sd"] == pytest.approx(0.683, abs=0.01) and c["within_2sd"] == pytest.approx(0.954, abs=0.01)
    assert f3.robust_sd(RNG.normal(0, 2, 20000)) == pytest.approx(2.0, rel=0.03)

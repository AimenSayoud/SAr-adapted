"""Known answers for the supervisor-table routines (C-071): windows, their statistics, the displacement sign and the
consecutive-pair rule."""
import numpy as np
import pandas as pd
import pytest

from insar_wetlands import field, plot_table as pt
from insar_wetlands.cube.core import WAVELENGTH_MM


def test_window_index_clips_at_the_edge_and_keeps_the_mask():
    rr, cc = pt.window_index(0, 0, (5, 5), h=1)
    assert sorted(zip(rr, cc)) == [(0, 0), (0, 1), (1, 0), (1, 1)]
    mask = np.ones((5, 5), bool)
    mask[2, 3] = False
    rr, cc = pt.window_index(2, 2, (5, 5), h=1, mask=mask)
    assert len(rr) == 8 and (2, 3) not in set(zip(rr, cc))


def test_stack_stats_match_the_single_window_routine():
    rng = np.random.default_rng(1)
    v = rng.normal(size=(6, 9))
    v[2, :4] = np.nan
    v[4] = np.nan
    s = pt.window_stats_stack(v)
    for t in range(6):
        ref = field.window_stats(v[t])
        for k in pt.STATS:
            assert s.loc[t, k] == pytest.approx(ref[k], nan_ok=True)
    assert s.loc[4, "n_valid"] == 0 and np.isnan(s.loc[4, "median"])


def test_uplift_gives_positive_displacement_and_the_reference_cancels():
    # a uniform uplift of d mm toward the satellite is a phase of −4πd/λ, plus a constant per pair that the
    # reference area also carries
    d, const = np.array([3.0, -5.0]), np.array([0.7, -1.2])
    phase = -4 * np.pi * d / WAVELENGTH_MM
    pix = phase[:, None] + const[:, None] + np.array([[-0.1, 0.0, 0.1]])
    out = pt.phase_to_los_mm(pt.window_stats_stack(pix), ref_rad=const)
    assert out["median"].to_numpy() == pytest.approx(d)
    assert out["sd"].to_numpy() == pytest.approx(np.full(2, 0.1 * WAVELENGTH_MM / (4 * np.pi)))


def test_consecutive_pairs_are_adjacent_dates():
    t = pd.to_datetime(["2021-01-01", "2021-01-07", "2021-01-13", "2021-01-25"], utc=True)
    t1 = [t[0], t[0], t[1], t[2], t[1]]
    t2 = [t[1], t[2], t[2], t[3], t[3]]
    assert pt.consecutive_mask(t1, t2).tolist() == [True, False, True, True, False]

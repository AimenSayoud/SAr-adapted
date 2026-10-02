"""regrid.to_template recovers known 40 m values from a 20 m field (D-024)."""
import numpy as np

from insar_wetlands.regrid import cell_index, to_template


def _grids(offset=0.0):
    tpl_x = 1000.0 + 40 * np.arange(5) + 20          # template cell centres, 40 m
    tpl_y = 5000.0 - 40 * np.arange(4) - 20           # north-up
    fx = 1000.0 + 20 * np.arange(10) + 10 + offset    # 20 m centres, aligned (offset 0) or shifted
    fy = 5000.0 - 20 * np.arange(8) - 10
    return tpl_x, tpl_y, fx, fy


def test_mean_recovers_block_values():
    tpl_x, tpl_y, fx, fy = _grids()
    truth = np.arange(20, dtype=float).reshape(4, 5)
    fine = np.kron(truth, np.ones((2, 2)))            # each 40 m cell = 2×2 fine pixels with its value
    assert np.allclose(to_template(fine, fx, fy, tpl_x, tpl_y), truth)


def test_phase_mean_does_not_cancel_near_pi():
    tpl_x, tpl_y, fx, fy = _grids()
    fine = np.full((8, 10), np.pi - 0.05)
    fine[::2] = -np.pi + 0.05                          # the same direction either side of ±π
    out = to_template(fine, fx, fy, tpl_x, tpl_y, kind="phase")
    assert np.allclose(np.abs(out), np.pi, atol=1e-9)   # an arithmetic mean would give ~0


def test_nan_and_min_count():
    tpl_x, tpl_y, fx, fy = _grids()
    fine = np.ones((8, 10))
    fine[0:2, 0:2] = np.nan                            # the first template cell has no data
    fine[2, 2] = np.nan                                # the next cell down-right has 3 of 4
    out = to_template(fine, fx, fy, tpl_x, tpl_y, min_count=4)
    assert np.isnan(out[0, 0]) and np.isnan(out[1, 1]) and out[0, 1] == 1


def test_unaligned_grid_assigns_by_centre():
    tpl_x, tpl_y, fx, fy = _grids(offset=15.0)         # fine grid shifted 15 m east
    r, c = cell_index(fx, fy, tpl_x, tpl_y)
    assert (c[c >= 0] == np.floor((fx[c >= 0] - 1000.0) / 40)).all()
    assert c[-1] == -1                                 # the last fine centre falls outside the template

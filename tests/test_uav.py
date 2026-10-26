import numpy as np
from scipy.ndimage import shift as nd_shift

from insar_wetlands import uav


def _scene(n=256, seed=0):
    """Boardwalk-like lines and platforms on a noisy, smooth background."""
    rng = np.random.default_rng(seed)
    a = rng.normal(0, 0.02, (n, n))
    a[:, 120:124] += 1.0
    a[60:66, 40:200] += 1.0
    a[170:190, 150:200] += 0.7
    return a


def test_phase_shift_recovers_a_known_subpixel_translation():
    ref = _scene()
    mov = nd_shift(ref, (3.4, -5.7), order=3, mode="nearest")
    dy, dx, q = uav.phase_shift(uav.edges(ref), uav.edges(mov))
    assert abs(dy - 3.4) < 0.25 and abs(dx + 5.7) < 0.25
    assert q > 3


def test_phase_shift_ignores_a_change_of_calibration():
    ref = _scene()
    mov = 2.7 * nd_shift(ref, (-2.0, 1.0), order=3, mode="nearest") + 0.3   # another gain and offset
    dy, dx, _ = uav.phase_shift(uav.edges(ref), uav.edges(mov))
    assert abs(dy + 2.0) < 0.25 and abs(dx - 1.0) < 0.25


def test_phase_shift_reports_no_match_on_unrelated_scenes():
    rng = np.random.default_rng(1)
    _, _, q = uav.phase_shift(uav.edges(rng.normal(size=(128, 128))), uav.edges(rng.normal(size=(128, 128))))
    assert q < 2


def test_grid_snaps_outward():
    tr, shape = uav.grid((10.2, 20.3, 12.9, 24.1), 0.5)
    assert (tr.c, tr.f) == (10.0, 24.5) and shape == (9, 6)

"""X-101 QC helpers on synthetic CSLCs: cell sums land in the right 40 m cells, coherence is 1 for identical scenes,
the Goldstein filter lowers phase noise without moving a smooth phase."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "experiments"))
import slc_qc_x101 as Q  # noqa: E402

GX = 585720.0 + 40 * np.arange(10)          # cube-like cell centres
GY = 5849440.0 - 40 * np.arange(8)


def _cslc(phase_fn, noise=0.0, seed=0):
    x = GX[0] - 17.5 + 5 * np.arange(80)     # 5 m columns, 8 per cell
    y = GY[0] + 15.0 - 10 * np.arange(32)    # 10 m rows, 4 per cell
    xx, yy = np.meshgrid(x, y)
    rng = np.random.default_rng(seed)
    s = np.exp(1j * (phase_fn(xx, yy) + noise * rng.normal(size=xx.shape))).astype("complex64")
    return s, x, y


def test_axis_index_puts_samples_in_their_cells():
    _, x, y = _cslc(lambda a, b: 0 * a)
    r, c = Q.axis_index(x, y, GX, GY)
    assert (np.bincount(c[c >= 0]) == 8).all() and (np.bincount(r[r >= 0]) == 4).all()


def test_look_recovers_phase_difference_and_full_coherence():
    s1, x, y = _cslc(lambda a, b: 0.3 * (a - a.min()) / 400)
    s2 = s1 * np.exp(-1j * 0.7)                                     # s1 conj(s2) has phase +0.7
    r, c = Q.axis_index(x, y, GX, GY)
    z, coh = Q.look(s1, s2, r, c, (len(GY), len(GX)))
    assert np.allclose(np.angle(z), 0.7, atol=1e-5) and np.allclose(coh, 1.0, atol=1e-5)


def test_goldstein_reduces_noise_and_keeps_a_smooth_phase():
    rng = np.random.default_rng(1)
    yy, xx = np.mgrid[0:64, 0:64]
    truth = 0.05 * xx
    z = np.exp(1j * (truth + 0.6 * rng.normal(size=truth.shape)))
    f = Q.goldstein(z)
    err_raw = np.angle(z * np.exp(-1j * truth)).std()
    err_f = np.angle(f * np.exp(-1j * truth)).std()
    assert err_f < 0.7 * err_raw
    assert abs(np.angle(np.mean(f * np.exp(-1j * truth)))) < 0.05


def test_verdict_rule():
    import pandas as pd
    t = pd.DataFrame({"track": ["ascending"] * 3, "filtered_spread_hard_rad": [0.2, 0.3, 0.25],
                      "expected_spread_hard_rad": [0.3, 0.3, 0.3], "filtered_offset_stable_rad": [0.05, -0.02, 0.0]})
    assert bool(Q.verdict(t).agrees.iloc[0])
    t["filtered_spread_hard_rad"] = 0.5
    assert not bool(Q.verdict(t).agrees.iloc[0])


def test_placement_rule():
    import pandas as pd
    ok = pd.DataFrame({"track": ["ascending", "descending"], "dx_m": [20.0, -20.0], "dy_m": [15.0, 15.0]})
    assert Q.placement_ok(ok)
    off = pd.DataFrame({"track": ["ascending", "descending"], "dx_m": [-20.0, 20.0 + 45], "dy_m": [0.0, 0.0]})
    assert not Q.placement_ok(off)
    same_sign = pd.DataFrame({"track": ["ascending", "descending"], "dx_m": [20.0, 10.0], "dy_m": [0.0, 0.0]})
    assert not Q.placement_ok(same_sign)


def test_best_shift_finds_a_known_offset():
    rng = np.random.default_rng(3)
    b = rng.normal(size=(40, 40))
    from scipy.ndimage import gaussian_filter, shift
    b = gaussian_filter(b, 2)
    a = shift(b, (0.0, -1.0), order=1, mode="nearest")          # a is b moved one cell west
    r, dy, dx = Q.best_shift(a, b, rng=2.0, step=0.25)
    assert abs(dx - 1.0) < 0.26 and abs(dy) < 0.26 and r > 0.9

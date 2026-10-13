"""slc_stats recovers known answers: a steady scatterer has low amplitude dispersion and speckle has ~0.52; the KS
statistic matches scipy and finds neighbours drawn from the same distribution; coherence of identical signals is 1;
the noise floor falls with looks and an incoherent estimate exceeds it about 5 % of the time."""
import numpy as np
from scipy import stats

from insar_wetlands import slc_stats as S


def test_amplitude_dispersion_steady_vs_speckle():
    rng = np.random.default_rng(0)
    t = 200
    speckle = (rng.normal(size=(t, 20, 20)) + 1j * rng.normal(size=(t, 20, 20))) / np.sqrt(2)
    steady = speckle * 0.05 + 3.0
    assert abs(np.nanmedian(S.amplitude_dispersion(speckle)) - 0.523) < 0.05     # Rayleigh: sqrt(4/pi - 1)
    assert np.nanmax(S.amplitude_dispersion(steady)) < 0.1


def test_ks_statistic_matches_scipy():
    rng = np.random.default_rng(1)
    x, y = rng.normal(size=(5, 40)), rng.normal(0.5, 1, size=(5, 40))
    ours = S.ks_statistic(x, y)
    ref = [stats.ks_2samp(x[i], y[i]).statistic for i in range(5)]
    assert np.allclose(ours, ref)


def test_shp_count_finds_the_homogeneous_neighbours():
    rng = np.random.default_rng(2)
    t = 60
    amp = np.abs(rng.normal(size=(t, 9, 9)) + 1j * rng.normal(size=(t, 9, 9)))
    amp[:, :, 5:] *= 6.0                                     # right part: a different, brighter population
    n = S.shp_count(amp, np.array([4]), np.array([2]), 2, 2)   # centre in the left part
    assert 15 <= n[0] <= 24                                   # its 24 neighbours: columns 0-4 alike, 5 differs


def test_boxcar_coherence_identical_is_one_and_noise_floor():
    rng = np.random.default_rng(3)
    s = rng.normal(size=(30, 30)) + 1j * rng.normal(size=(30, 30))
    assert np.allclose(S.boxcar_coherence(s, s, (3, 3))[2:-2, 2:-2], 1.0)
    f9, f45 = S.noise_floor(9), S.noise_floor(45)
    assert f45 < f9 < 0.7
    a = rng.normal(size=(4000, 9)) + 1j * rng.normal(size=(4000, 9))
    b = rng.normal(size=(4000, 9)) + 1j * rng.normal(size=(4000, 9))
    g = np.abs((a * np.conj(b)).sum(1)) / np.sqrt((np.abs(a) ** 2).sum(1) * (np.abs(b) ** 2).sum(1))
    assert 0.03 < np.mean(g > f9) < 0.07


def test_shp_coherence_of_identical_signals_is_one():
    rng = np.random.default_rng(4)
    s = rng.normal(size=(12, 12)) + 1j * rng.normal(size=(12, 12))
    amp = np.abs(rng.normal(size=(30, 12, 12)) + 1j * rng.normal(size=(30, 12, 12)))
    rows, cols = np.array([5, 6]), np.array([5, 6])
    mask, offs = S.shp_mask(amp, rows, cols, 1, 2)
    assert mask.shape == (2, len(offs)) and len(offs) == 14
    assert np.allclose(S.shp_coherence(s, s, rows, cols, mask, offs), 1.0)

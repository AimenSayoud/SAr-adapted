"""Test synthetique du phase-linking EVD LÉGER (sans ISCE).

Verifie que : (1) l'EVD retrouve un historique de phase connu a partir d'un
reseau SPARSE d'interferogrammes enroules, (2) la temporal_coherence est
ELEVEE quand les paires sont coherentes et BASSE quand le pixel est bruyant/
decorrele — c'est le discriminant qui teste H1 (le phase-linking recupere-t-il
un signal la ou le WLS echoue ?).

Execution : python tests/test_synthetic_phaselinking.py
"""

import numpy as np
import pandas as pd
import xarray as xr

from insar_wetlands.inversion.isbas import PHASE_TO_MM
from insar_wetlands.inversion.phaselinking import (evd_phase_linking, evd_pixel,
                                                   _pair_date_index)


def _wrap(x):
    return np.angle(np.exp(1j * x))


def _sparse_network(dates, max_days=48):
    return [f"{a:%Y%m%d}_{b:%Y%m%d}"
            for i, a in enumerate(dates) for b in dates[i + 1:]
            if 0 < (b - a).days <= max_days]


def test_evd_pixel_recovers_clean_history():
    """Pixel coherent : EVD retrouve l'historique, tcoh proche de 1."""
    dates = pd.date_range("2022-01-01", periods=15, freq="12D")
    pairs = _sparse_network(dates)
    _, idx = _pair_date_index(pairs)
    n = len(dates)
    # historique de phase petit (evite l'ambiguite d'enroulement)
    rng = np.random.default_rng(0)
    truth = np.cumsum(rng.normal(0, 0.15, n))
    truth -= truth[0]
    phi = np.array([_wrap(truth[j] - truth[i]) for (i, j) in idx])
    coh = np.full(len(pairs), 0.85)
    theta, tcoh = evd_pixel(phi + rng.normal(0, 0.02, len(pairs)), coh, idx, n)
    # historique recupere modulo une constante -> comparer apres recalage
    err = _wrap(theta - truth)
    assert np.nanmax(np.abs(err)) < 0.2, np.nanmax(np.abs(err))
    assert tcoh > 0.95, tcoh


def test_evd_pixel_low_tcoh_when_incoherent():
    """Pixel decorrele (phases aleatoires) : tcoh BASSE."""
    dates = pd.date_range("2022-01-01", periods=15, freq="12D")
    pairs = _sparse_network(dates)
    _, idx = _pair_date_index(pairs)
    n = len(dates)
    rng = np.random.default_rng(1)
    phi = rng.uniform(-np.pi, np.pi, len(pairs))   # bruit pur
    coh = np.full(len(pairs), 0.3)
    _, tcoh = evd_pixel(phi, coh, idx, n)
    # plancher de bruit ~0.55 (redondance finie n_pairs/(n-1)~3.5) : bien
    # SOUS les ~0.97 d'un pixel coherent -> discrimination nette.
    assert tcoh < 0.7, tcoh


def test_evd_stack_separates_good_from_bad():
    """Sur un stack 2D : tcoh discrimine un pixel bon d'un pixel bruyant."""
    dates = pd.date_range("2022-01-01", periods=15, freq="12D")
    pairs = _sparse_network(dates)
    _, idx = _pair_date_index(pairs)
    n = len(dates)
    ny, nx = 3, 3
    rng = np.random.default_rng(2)
    truth = np.cumsum(rng.normal(0, 0.12, n)); truth -= truth[0]
    ph = np.zeros((len(pairs), ny, nx), "float32")
    ch = np.full((len(pairs), ny, nx), 0.8, "float32")
    for k, (i, j) in enumerate(idx):
        ph[k, :, :] = _wrap(truth[j] - truth[i])
    ph += rng.normal(0, 0.02, ph.shape)
    # pixel (0,0) : bruit pur + basse coherence
    ph[:, 0, 0] = rng.uniform(-np.pi, np.pi, len(pairs))
    ch[:, 0, 0] = 0.25
    coords = {"pair": pairs, "y": np.arange(ny) * 40.0, "x": np.arange(nx) * 40.0}
    wrapped = xr.DataArray(ph, dims=("pair", "y", "x"), coords=coords)
    corr = xr.DataArray(ch, dims=("pair", "y", "x"), coords=coords)
    ds = evd_phase_linking(wrapped, corr)
    tc = ds.temporal_coherence.values
    assert tc[0, 0] < 0.7, tc[0, 0]
    assert np.nanmin(np.delete(tc.ravel(), 0)) > 0.9, tc
    assert "velocity_mm_yr" in ds and "displacement_mm" in ds


if __name__ == "__main__":
    test_evd_pixel_recovers_clean_history()
    test_evd_pixel_low_tcoh_when_incoherent()
    test_evd_stack_separates_good_from_bad()
    print("ALL PHASE-LINKING TESTS PASSED")

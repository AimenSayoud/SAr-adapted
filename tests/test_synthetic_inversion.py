"""Test synthetique de bout en bout du coeur scientifique (Phases 9/10).

Verifie que : (1) l'inversion ISBAS retrouve une subsidence lineaire connue,
(2) un pixel intermittent (40 % de paires invalides) est bien 'sauve',
(3) la fermeture des triplets est nulle sur un reseau propre.

Execution : python -m pytest tests/  (ou python tests/test_synthetic_inversion.py)
"""

import numpy as np
import pandas as pd
import xarray as xr

from insar_wetlands.compare import fit_velocity
from insar_wetlands.inversion.isbas import PHASE_TO_MM, invert_stack, phase_closure
from insar_wetlands.network import design_matrix


def _synthetic():
    dates = pd.date_range("2022-01-01", periods=20, freq="12D")
    pairs = [f"{a:%Y%m%d}_{b:%Y%m%d}"
             for i, a in enumerate(dates) for b in dates[i + 1:]
             if (b - a).days <= 48]
    ny, nx = 6, 5
    t_years = ((dates - dates[0]).days / 365.25).values
    truth_mm = np.zeros((len(dates), ny, nx))
    truth_mm[:, 2:, 2:] = -30 * t_years[:, None, None]
    phase = np.zeros((len(pairs), ny, nx))
    d = {dt: k for k, dt in enumerate(dates)}
    for k, p in enumerate(pairs):
        a, b = pd.Timestamp(p[:8]), pd.Timestamp(p[9:])
        phase[k] = (truth_mm[d[b]] - truth_mm[d[a]]) / PHASE_TO_MM
    phase += np.random.default_rng(0).normal(0, 0.05, phase.shape)
    coords = {"pair": pairs, "y": np.arange(ny) * 40.0,
              "x": np.arange(nx) * 40.0}
    unw = xr.DataArray(phase, dims=("pair", "y", "x"), coords=coords)
    corr = xr.full_like(unw, 0.8)
    dry = xr.DataArray(np.ones(phase.shape, bool),
                       dims=("pair", "y", "x"), coords=coords)
    dry[::3, 3, 3] = False  # pixel intermittent
    return dates, pairs, unw, corr, dry


def test_design_matrix_shape():
    dates, pairs, *_ = _synthetic()
    A = design_matrix(pairs, dates)
    assert A.shape == (len(pairs), len(dates) - 1)


def test_isbas_recovers_truth():
    dates, pairs, unw, corr, dry = _synthetic()
    res = invert_stack(unw, corr, dry, (0.0, 0.0), min_pairs=10, gamma_min=0.3)
    vel = fit_velocity(res.los_displacement_mm)
    assert abs(float(vel.velocity_mm_yr.isel(y=4, x=4)) + 30) < 2
    assert abs(float(vel.velocity_mm_yr.isel(y=0, x=1))) < 2
    # pixel intermittent sauve par l'ISBAS
    assert abs(float(vel.velocity_mm_yr.isel(y=3, x=3)) + 30) < 3


def test_phase_closure_clean_network():
    _, _, unw, *_ = _synthetic()
    assert float(phase_closure(unw).mean()) < 0.05


if __name__ == "__main__":
    test_design_matrix_shape()
    test_isbas_recovers_truth()
    test_phase_closure_clean_network()
    print("ALL TESTS PASSED")

"""Test synthétique du modèle prédictif d'échec (Phase H) + RVI.

Vérifie que : (1) rank_covariates retrouve le prédicteur planté et ignore le
bruit ; (2) la régression standardisée donne le bon signe/ordre et un R² validé
croisé élevé quand la relation est réelle, ~0 quand elle est absente ;
(3) threshold_sweep produit une courbe monotone décroissante ; (4) le RVI
dual-pol vaut ~0 pour une surface et ~1 pour un volume dépolarisant.

Execution : python tests/test_predict_failure.py
"""

import numpy as np
import pandas as pd
import xarray as xr

from insar_wetlands.predict_failure import (covariate_table, fit_failure_model,
                                            rank_covariates, threshold_sweep)
from insar_wetlands.stratify import dual_pol_rvi

NY = NX = 30


def _da(a):
    return xr.DataArray(a, dims=("y", "x"),
                        coords={"y": np.arange(NY) * 40.0, "x": np.arange(NX) * 40.0})


def _synthetic():
    rng = np.random.default_rng(0)
    wet = rng.uniform(0, 1, (NY, NX))        # vrai moteur
    rvi = rng.uniform(0, 1, (NY, NX))        # 2e moteur, plus faible
    junk = rng.normal(0, 1, (NY, NX))        # bruit pur
    tcoh = 0.9 - 0.5 * wet - 0.2 * rvi + rng.normal(0, 0.03, (NY, NX))
    mask = np.ones((NY, NX), bool)
    return (_da(tcoh), {"wetness": _da(wet), "rvi": _da(rvi), "junk": _da(junk)},
            _da(mask))


def test_ranking_finds_real_driver():
    tcoh, covars, mask = _synthetic()
    df = covariate_table(tcoh, covars, mask)
    r = rank_covariates(df).set_index("covariate")
    assert r.loc["wetness", "spearman_rho"] < -0.7          # moteur dominant
    assert abs(r.loc["junk", "spearman_rho"]) < 0.2         # bruit ignoré
    # l'ordre d'importance est correct
    assert abs(r.loc["wetness", "spearman_rho"]) > abs(r.loc["rvi", "spearman_rho"])


def test_model_predicts_out_of_sample():
    tcoh, covars, mask = _synthetic()
    df = covariate_table(tcoh, covars, mask)
    res = fit_failure_model(df)
    assert res["r2_cv_mean"] > 0.8, res                     # vrai pouvoir prédictif
    c = res["coefficients"].set_index("covariate")
    assert c.loc["wetness", "std_coef"] < 0                 # bon signe
    assert abs(c.loc["junk", "std_coef"]) < 0.1             # bruit ~ sans effet


def test_model_r2_near_zero_when_no_relation():
    """Garde-fou : sans relation réelle, le R² validé croisé ne doit PAS être élevé."""
    rng = np.random.default_rng(1)
    tcoh = _da(rng.normal(0, 1, (NY, NX)))
    covars = {f"x{i}": _da(rng.normal(0, 1, (NY, NX))) for i in range(3)}
    df = covariate_table(tcoh, covars, _da(np.ones((NY, NX), bool)))
    assert fit_failure_model(df)["r2_cv_mean"] < 0.15


def test_threshold_sweep_monotone():
    tcoh, _, _ = _synthetic()
    A = np.zeros((NY, NX), bool); A[:15] = True
    C = np.zeros((NY, NX), bool); C[15:] = True
    sw = threshold_sweep(tcoh, {"A": _da(A), "C": _da(C)}, zone_names=("A", "C"))
    for z, g in sw.groupby("zone"):
        f = g.sort_values("threshold")["frac_above"].values
        assert np.all(np.diff(f) <= 1e-9), f      # décroissante en seuil


def test_dual_pol_rvi_surface_vs_volume():
    """RVI ~0 pour une surface (VH<<VV) et ~1 pour un volume (VH~VV)."""
    surface = xr.Dataset({"gamma0_vv_db": _da(np.full((NY, NX), -8.0)),
                          "gamma0_vh_db": _da(np.full((NY, NX), -28.0))})
    volume = xr.Dataset({"gamma0_vv_db": _da(np.full((NY, NX), -10.0)),
                         "gamma0_vh_db": _da(np.full((NY, NX), -10.0))})
    assert float(dual_pol_rvi(surface).mean()) < 0.1
    assert abs(float(dual_pol_rvi(volume).mean()) - 2.0) < 1e-6


if __name__ == "__main__":
    test_ranking_finds_real_driver()
    test_model_predicts_out_of_sample()
    test_model_r2_near_zero_when_no_relation()
    test_threshold_sweep_monotone()
    test_dual_pol_rvi_surface_vs_volume()
    print("ALL PHASE-H TESTS PASSED")

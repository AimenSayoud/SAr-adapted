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

from insar_wetlands.predict_failure import (
    covariate_table,
    fit_failure_model,
    rank_covariates,
    threshold_sweep,
)
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


def test_collinearity_detects_monotone_duplicate():
    """Cas REEL : rvi=4r/(1+r) et vh_vv_db=10log10(r) sont deux transformations
    monotones du meme rapport -> VIF enorme, coefficients ininterpretables."""
    from insar_wetlands.predict_failure import collinearity_report, drop_redundant

    rng = np.random.default_rng(4)
    r = rng.uniform(0.05, 0.9, 400)                 # rapport VH/VV
    df = pd.DataFrame({
        "target": 0.8 - 0.3 * r + rng.normal(0, 0.05, 400),
        "rvi": 4 * r / (1 + r),
        "vh_vv_db": 10 * np.log10(r),
        "independent": rng.normal(0, 1, 400),
    })
    rep = collinearity_report(df).set_index("covariate")
    assert rep.loc["rvi", "vif"] > 10 and rep.loc["vh_vv_db", "vif"] > 10, rep
    assert not rep.loc["independent", "redundant"], rep
    # les coefficients bruts sont aberrants (gros, signes opposes)
    raw = fit_failure_model(df)["coefficients"].set_index("covariate")["std_coef"]
    assert abs(raw["rvi"]) > 1.0 and raw["rvi"] * raw["vh_vv_db"] < 0, raw
    # apres nettoyage : une seule des deux survit, RVI protegee
    red, dropped = drop_redundant(df, keep=("rvi",))
    assert "vh_vv_db" in dropped and "rvi" in red.columns, (dropped, red.columns)
    clean = fit_failure_model(red)["coefficients"].set_index("covariate")["std_coef"]
    assert abs(clean["rvi"]) < 1.0, clean          # coefficient redevenu sain


def test_zone_fraction_excludes_invalid_pixels():
    """The bug this guards against produced a wrong number in the manuscript.

    `np.nanmean(values >= 0.7)` looks like it skips invalid pixels. It does not:
    `NaN >= 0.7` is False, not NaN, so the boolean array has nothing for nanmean
    to skip and masked pixels are counted as failures. Zones with full coverage
    still agree with the correct value, so the error stays hidden until a
    fragmented zone is compared against a complete one."""
    import numpy as np
    import xarray as xr

    from insar_wetlands.predict_failure import threshold_sweep, zone_fraction_above, zone_values
    vals = np.array([[0.9, 0.9, 0.9, 0.5],
                     [np.nan, np.nan, np.nan, np.nan]])
    field = xr.DataArray(vals, dims=("y", "x"))
    full = xr.DataArray(np.array([[True] * 4, [False] * 4]), dims=("y", "x"))
    holed = xr.DataArray(np.array([[True] * 4, [True] * 4]), dims=("y", "x"))
    zones = {"full": full, "holed": holed}

    # 3 of 4 valid pixels are >= 0.7 in both zones: the answer must not depend
    # on how many masked pixels the zone happens to contain.
    assert zone_fraction_above(field, zones, "full", 0.7) == 0.75
    assert zone_fraction_above(field, zones, "holed", 0.7) == 0.75
    assert len(zone_values(field, zones, "holed")) == 4

    # the buggy formulation, kept as a witness to what it would have returned
    naive = float(np.nanmean(field.values[holed.values] >= 0.7))
    assert naive == 0.375, naive
    assert naive != zone_fraction_above(field, zones, "holed", 0.7)

    # the sweep must agree with the single-threshold helper at that threshold
    sw = threshold_sweep(field, zones, thresholds=[0.7], zone_names=("full", "holed"))
    for z in ("full", "holed"):
        assert float(sw[sw.zone == z].frac_above.iloc[0]) == \
            zone_fraction_above(field, zones, z, 0.7), z
    print("  zone fractions ignore masked pixels; sweep and helper agree")


if __name__ == "__main__":
    test_collinearity_detects_monotone_duplicate()
    test_ranking_finds_real_driver()
    test_model_predicts_out_of_sample()
    test_model_r2_near_zero_when_no_relation()
    test_threshold_sweep_monotone()
    test_dual_pol_rvi_surface_vs_volume()
    test_zone_fraction_excludes_invalid_pixels()
    print("ALL PHASE-H TESTS PASSED")

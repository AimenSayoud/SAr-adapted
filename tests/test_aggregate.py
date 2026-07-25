"""Test synthétique de l'agrégation spatiale (Phase G).

Vérifie que les trois observables se comportent comme la physique le prédit :
(1) |R| discrimine phases ALÉATOIRES (~plancher 1/sqrt(N_eff)) de phase COMMUNE ;
(2) la double différence agrégée + inversion retrouvent un déplacement connu
    ENFOUI sous un bruit par-pixel qui rend chaque pixel inexploitable ;
(3) le biais de fermeture est NUL pour un déplacement pur et NON NUL pour une
    perturbation de type diélectrique — le discriminateur du mécanisme.

Execution : python tests/test_aggregate.py
"""

import numpy as np
import pandas as pd
import xarray as xr

from insar_wetlands.aggregate import (aggregate_unwrapped, closure_bias_by_zone,
                                      double_difference, invert_aggregate,
                                      phasor_verdict, zone_phasor)
from insar_wetlands.inversion.isbas import PHASE_TO_MM

NY = NX = 20


def _grid(dates):
    pairs = [f"{a:%Y%m%d}_{b:%Y%m%d}"
             for i, a in enumerate(dates) for b in dates[i + 1:]
             if 0 < (b - a).days <= 48]
    coords = {"pair": pairs, "y": np.arange(NY) * 40.0, "x": np.arange(NX) * 40.0}
    return pairs, coords


def _zones():
    A = np.zeros((NY, NX), bool); A[2:12, 2:12] = True     # 100 px
    C = np.zeros((NY, NX), bool); C[2:12, 13:19] = True    # 60 px
    mk = lambda m: xr.DataArray(m, dims=("y", "x"),
                                coords={"y": np.arange(NY) * 40.0,
                                        "x": np.arange(NX) * 40.0})
    return {"A": mk(A), "C": mk(C)}


def test_phasor_separates_random_from_common():
    """|R| ~ plancher si phases aléatoires ; |R| >> plancher si phase commune."""
    dates = pd.date_range("2022-01-01", periods=10, freq="12D")
    pairs, coords = _grid(dates)
    zones = _zones()
    rng = np.random.default_rng(0)
    ph = rng.uniform(-np.pi, np.pi, (len(pairs), NY, NX))   # A = bruit pur
    # zone C : phase COMMUNE (constante par paire) + petit bruit
    for k in range(len(pairs)):
        ph[k][zones["C"].values] = 0.7 + rng.normal(0, 0.1, zones["C"].values.sum())
    wrapped = xr.DataArray(ph, dims=("pair", "y", "x"), coords=coords)
    corr = xr.full_like(wrapped, 0.4)
    v = phasor_verdict(zone_phasor(wrapped, corr, zones)).set_index("zone")
    assert v.loc["A", "ratio_to_floor"] < 1.0, v      # bruit : sous le plancher
    assert v.loc["C", "ratio_to_floor"] > 3.0, v      # signal commun : au-dessus


def test_aggregation_recovers_buried_signal():
    """Un déplacement noyé sous un bruit par-pixel énorme est récupéré."""
    dates = pd.date_range("2022-01-01", periods=14, freq="12D")
    pairs, coords = _grid(dates)
    zones = _zones()
    rng = np.random.default_rng(1)
    t_years = (dates - dates[0]).days.values / 365.25
    truth_mm = -20.0 * t_years                      # tapis : -20 mm/an
    di = {d: i for i, d in enumerate(dates)}
    ph = np.zeros((len(pairs), NY, NX))
    for k, p in enumerate(pairs):
        a, b = pd.Timestamp(p[:8]), pd.Timestamp(p[9:])
        dphi = (truth_mm[di[b]] - truth_mm[di[a]]) / PHASE_TO_MM
        ph[k][zones["A"].values] = dphi              # A bouge
        # bruit par-pixel ENORME (1.5 rad) -> aucun pixel exploitable seul
        ph[k] += rng.normal(0, 1.5, (NY, NX))
    unw = xr.DataArray(ph, dims=("pair", "y", "x"), coords=coords)
    corr = xr.full_like(unw, 0.4)
    dd = aggregate_unwrapped(unw, corr, zones, "A", "C")
    res = invert_aggregate(dd)
    assert res.attrs["n_obs"] > res.attrs["n_unknowns"]      # surdéterminé
    assert abs(res.attrs["velocity_mm_yr"] - (-20.0)) < 4.0, res.attrs
    # la double différence enroulée existe aussi et reste bornée
    ddw = double_difference(xr.apply_ufunc(lambda a: np.angle(np.exp(1j * a)), unw),
                            corr, zones, "A", "C")
    assert len(ddw) == len(pairs)
    assert np.all(np.abs(ddw["ddphase_rad"]) <= np.pi + 1e-9)


def test_closure_bias_zero_for_motion_nonzero_for_dielectric():
    """Déplacement pur -> fermeture nulle ; perturbation diélectrique -> biais."""
    dates = pd.date_range("2022-01-01", periods=12, freq="12D")
    pairs, coords = _grid(dates)
    zones = _zones()
    rng = np.random.default_rng(2)
    di = {d: i for i, d in enumerate(dates)}
    theta = np.cumsum(rng.normal(0, 0.2, len(dates)))     # historique de phase
    ph = np.zeros((len(pairs), NY, NX))
    for k, p in enumerate(pairs):
        a, b = pd.Timestamp(p[:8]), pd.Timestamp(p[9:])
        dphi = theta[di[b]] - theta[di[a]]                # additif -> ferme a 0
        ph[k][zones["C"].values] = dphi                   # C : mouvement pur
        # A : mouvement + biais CONSTANT par paire (non additif) = diélectrique.
        # closure = b0 + b0 - b0 = b0 != 0
        ph[k][zones["A"].values] = dphi + 0.30
        ph[k] += rng.normal(0, 0.05, (NY, NX))
    wrapped = xr.DataArray(np.angle(np.exp(1j * ph)), dims=("pair", "y", "x"),
                           coords=coords)
    cb = closure_bias_by_zone(wrapped, zones, max_triplets=200).set_index("zone")
    assert abs(cb.loc["C", "mean_closure_rad"]) < 0.05, cb   # mouvement -> ~0
    assert abs(cb.loc["A", "mean_closure_rad"]) > 0.2, cb    # diélectrique -> biais
    assert bool(cb.loc["A", "bias_significant"]) is True
    assert bool(cb.loc["C", "bias_significant"]) is False


def test_seasonal_amplitude_recovers_breathing():
    """Une respiration saisonnière connue est retrouvée — et la VITESSE est ~0,
    ce qui montre qu'un test de vitesse n'a AUCUNE puissance sur ce signal."""
    from insar_wetlands.aggregate import seasonal_amplitude

    d = pd.date_range("2022-01-01", periods=90, freq="12D")
    t = (d - d[0]).days.values / 365.25
    rng = np.random.default_rng(3)
    y = 25.0 * np.cos(2 * np.pi * (t - 0.4)) + rng.normal(0, 2.0, t.size)
    s = seasonal_amplitude(pd.DataFrame({"date": d, "disp_mm": y}))
    assert abs(s["amplitude_mm"] - 25.0) < 2.0, s
    assert s["r2_seasonal"] > 0.9, s
    # la tendance linéaire est ~nulle : tester une vitesse aurait tout raté
    assert abs(s["trend_mm_yr"]) < 3.0, s
    # bruit pur -> amplitude faible
    n = seasonal_amplitude(pd.DataFrame({"date": d,
                                         "disp_mm": rng.normal(0, 2.0, t.size)}))
    assert n["amplitude_mm"] < 2.0, n


def test_correlation_length_and_neff():
    """N_eff doit etre MESURE, pas suppose : un champ lisse a une longueur de
    correlation plus grande qu'un champ de bruit blanc, donc moins de pixels
    independants."""
    from insar_wetlands.aggregate import correlation_length, effective_looks

    rng = np.random.default_rng(11)
    ny = nx = 40
    coords = {"y": np.arange(ny) * 40.0, "x": np.arange(nx) * 40.0}
    mask = xr.DataArray(np.ones((ny, nx), bool), dims=("y", "x"), coords=coords)
    white = xr.DataArray(rng.normal(0, 1, (ny, nx)), dims=("y", "x"), coords=coords)
    # champ lisse : bruit blanc convolue par une moyenne glissante 5x5
    k = 5
    sm = np.copy(white.values)
    for _ in range(3):
        sm = np.pad(sm, k // 2, mode="edge")
        sm = np.array([[sm[i:i + k, j:j + k].mean() for j in range(nx)]
                       for i in range(ny)])
    smooth = xr.DataArray(sm, dims=("y", "x"), coords=coords)

    l_white = correlation_length(white, mask)
    l_smooth = correlation_length(smooth, mask)
    assert l_smooth > l_white, (l_white, l_smooth)
    # plus la correlation est longue, moins il y a d'echantillons independants
    assert (effective_looks(mask, field=smooth)
            < effective_looks(mask, field=white)), "N_eff doit chuter"


def test_exclude_winter_months():
    """Test de falsification : on doit pouvoir retirer les paires hivernales
    (neige/gel = explication alternative du signal saisonnier)."""
    from insar_wetlands.aggregate import filter_pairs

    dd = pd.DataFrame({
        "pair": ["20220115_20220127",   # janvier -> exclu
                 "20220601_20220613",   # ete     -> garde
                 "20221201_20221213",   # decembre-> exclu
                 "20220701_20220713"],  # ete     -> garde
        "dt_days": [12, 12, 12, 12], "weight": [1.0] * 4})
    kept = filter_pairs(dd, max_dt_days=None, exclude_months=(12, 1, 2))
    assert list(kept["pair"]) == ["20220601_20220613", "20220701_20220713"], kept


def test_matched_null_respects_sizes_and_pvalue():
    """Le nul doit avoir la MEME taille que les zones reelles (le bruit d'un
    agregat decroit en 1/sqrt(N) : un nul 4x plus grand sous-estime le plancher
    et fabrique de fausses detections)."""
    from insar_wetlands.aggregate import (empirical_pvalue, matched_null_zones)

    D = np.zeros((NY, NX), bool); D[1:19, 1:19] = True
    zones = {"D": xr.DataArray(D, dims=("y", "x"),
                               coords={"y": np.arange(NY) * 40.0,
                                       "x": np.arange(NX) * 40.0})}
    tmpl = xr.DataArray(np.zeros((NY, NX)), dims=("y", "x"), coords=zones["D"].coords)
    zn = matched_null_zones(zones, tmpl, n_target=100, n_reference=60, seed=0)
    assert int(zn["A"].sum()) == 100 and int(zn["C"].sum()) == 60
    assert int((zn["A"] & zn["C"]).sum()) == 0
    # deux graines differentes -> deux tirages differents
    zn2 = matched_null_zones(zones, tmpl, 100, 60, seed=7)
    assert not np.array_equal(zn["A"].values, zn2["A"].values)
    # p-value empirique : jamais 0, et correcte aux deux extremes
    nulls = np.linspace(0, 1, 50)
    assert empirical_pvalue(2.0, nulls)["p_value"] < 0.05     # observe >> nul
    assert empirical_pvalue(0.0, nulls)["p_value"] > 0.9      # observe <= tout


def test_seasonal_zone_scan_separates_breathing_zone():
    """Une zone qui 'respire' est detectee ; une zone plate ne l'est pas.
    C'est le design du test lac(B)-vs-tapis(A) : le lac ne peut pas respirer
    mecaniquement, donc s'il oscille aussi, le signal est dielectrique."""
    from insar_wetlands.aggregate import seasonal_zone_scan

    ny = nx = 26
    dates = pd.date_range("2022-01-01", periods=40, freq="12D")
    pairs = [f"{a:%Y%m%d}_{b:%Y%m%d}"
             for i, a in enumerate(dates) for b in dates[i + 1:]
             if 0 < (b - a).days <= 36]
    coords = {"pair": pairs, "y": np.arange(ny) * 40.0, "x": np.arange(nx) * 40.0}
    tmpl = xr.DataArray(np.zeros((ny, nx)), dims=("y", "x"),
                        coords={"y": coords["y"], "x": coords["x"]})
    mk = lambda m: xr.DataArray(m, coords=tmpl.coords, dims=tmpl.dims)
    A = np.zeros((ny, nx), bool); A[2:8, 2:8] = True        # respire
    B = np.zeros((ny, nx), bool); B[2:8, 9:15] = True       # plate
    C = np.zeros((ny, nx), bool); C[2:8, 16:22] = True      # reference
    D = np.zeros((ny, nx), bool); D[10:25, 1:25] = True     # reservoir de nuls
    zones = {"A": mk(A), "B": mk(B), "C": mk(C), "D": mk(D)}

    di = {d: i for i, d in enumerate(dates)}
    t = (dates - dates[0]).days.values / 365.25
    breath = 30.0 * np.cos(2 * np.pi * t)                   # respiration connue
    rng = np.random.default_rng(5)
    ph = rng.normal(0, 0.8, (len(pairs), ny, nx))
    for k, p in enumerate(pairs):
        a, b = pd.Timestamp(p[:8]), pd.Timestamp(p[9:])
        ph[k][A] += (breath[di[b]] - breath[di[a]]) / PHASE_TO_MM
    unw = xr.DataArray(ph, dims=("pair", "y", "x"), coords=coords)
    corr = xr.full_like(unw, 0.5)

    scan = seasonal_zone_scan(unw, corr, zones, tmpl, reference="C",
                              targets=("A", "B"), n_trials=30).set_index("zone")
    assert scan.loc["A", "amplitude_mm"] > 20, scan     # respiration retrouvee
    assert scan.loc["A", "p_value"] < 0.10, scan
    assert scan.loc["B", "amplitude_mm"] < 10, scan     # zone plate : rien
    assert scan.loc["B", "p_value"] > scan.loc["A", "p_value"], scan


if __name__ == "__main__":
    test_seasonal_zone_scan_separates_breathing_zone()
    test_correlation_length_and_neff()
    test_exclude_winter_months()
    test_matched_null_respects_sizes_and_pvalue()
    test_phasor_separates_random_from_common()
    test_aggregation_recovers_buried_signal()
    test_closure_bias_zero_for_motion_nonzero_for_dielectric()
    test_seasonal_amplitude_recovers_breathing()
    print("ALL AGGREGATE TESTS PASSED")

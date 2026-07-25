"""Tests Phase I (lien hydrologique) + visualisation/validation des zones.

Verifie que : (1) lag_scan retrouve un forcage plante ET son DECALAGE, en
ignorant un forcage sans lien ; (2) le detrend evite la correlation fallacieuse
entre deux series simplement tendancielles ; (3) driver_pvalues se comporte
correctement aux deux extremes ; (4) zone_areas convertit bien en hectares et
mesure l'ecart a la surface connue du site ; (5) les fonctions de trace
s'executent sans erreur (backend Agg).

Execution : python tests/test_phase_i.py
"""

import numpy as np
import pandas as pd
import xarray as xr

try:                        # matplotlib est present sur Colab, pas forcement en CI
    import matplotlib
    matplotlib.use("Agg")
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

from insar_wetlands.hydro_link import driver_pvalues, lag_scan
from insar_wetlands.zone_viz import (zone_areas, zone_field_table,
                                     zone_label_array)

NY = NX = 12


def _tmpl():
    return xr.DataArray(np.zeros((NY, NX)), dims=("y", "x"),
                        coords={"y": np.arange(NY) * 40.0,
                                "x": np.arange(NX) * 40.0})


def _zones(tmpl):
    mk = lambda m: xr.DataArray(m, coords=tmpl.coords, dims=tmpl.dims)
    A = np.zeros((NY, NX), bool); A[2:6, 2:6] = True       # 16 px
    B = np.zeros((NY, NX), bool); B[6:7, 2:6] = True       # 4 px
    C = np.zeros((NY, NX), bool); C[2:6, 7:10] = True      # 12 px
    D = np.zeros((NY, NX), bool); D[8:11, 1:11] = True     # 30 px
    return {"A": mk(A), "B": mk(B), "C": mk(C), "D": mk(D)}


def test_lag_scan_finds_driver_and_lag():
    """Le forcage plante est retrouve, avec son decalage ; le bruit est ignore."""
    dates = pd.date_range("2022-01-01", periods=120, freq="12D")
    days = pd.date_range("2021-10-01", "2026-01-01", freq="1D")
    rng = np.random.default_rng(0)
    tt = (days - days[0]).days.values / 365.25
    real = pd.Series(np.sin(2 * np.pi * tt), index=days)
    junk = pd.Series(rng.normal(0, 1, len(days)), index=days)
    drivers = pd.DataFrame({"real": real, "junk": junk})
    LAG = 24
    obs = real.shift(LAG).reindex(dates, method="nearest").values
    series = pd.DataFrame({"date": dates,
                           "disp_mm": 5 * obs + rng.normal(0, 0.5, len(dates))})
    sc = lag_scan(series, drivers, max_lag_days=60, step=6).set_index("driver")
    assert abs(sc.loc["real", "r"]) > 0.85, sc
    assert abs(sc.loc["real", "lag_days"] - LAG) <= 6, sc     # decalage retrouve
    assert abs(sc.loc["junk", "r"]) < 0.5, sc


def test_detrend_kills_spurious_trend_correlation():
    """Deux series purement tendancielles ne doivent PAS sembler correlees."""
    dates = pd.date_range("2022-01-01", periods=100, freq="12D")
    days = pd.date_range("2021-12-01", "2026-01-01", freq="1D")
    drivers = pd.DataFrame({"trend": pd.Series(
        np.linspace(0, 10, len(days)), index=days)})
    series = pd.DataFrame({"date": dates,
                           "disp_mm": np.linspace(0, 50, len(dates))})
    with_dt = lag_scan(series, drivers, max_lag_days=12, step=6, detrend=True)
    no_dt = lag_scan(series, drivers, max_lag_days=12, step=6, detrend=False)
    assert abs(no_dt.iloc[0]["r"]) > 0.95                     # piege sans detrend
    assert abs(with_dt.iloc[0]["r"]) < 0.5, with_dt           # neutralise


def test_deseasonalize_kills_shared_annual_cycle():
    """Le piege REEL rencontre en Phase I : deux signaux a cycle annuel
    INDEPENDANTS correlent fortement a un certain decalage (le balayage aligne
    les phases). Seule la desaisonnalisation le revele."""
    dates = pd.date_range("2022-01-01", periods=90, freq="12D")
    days = pd.date_range("2021-10-01", "2026-01-01", freq="1D")
    rng = np.random.default_rng(3)
    td = (days - days[0]).days.values / 365.25
    ts = (dates - dates[0]).days.values / 365.25
    # forcage saisonnier pur + bruit PROPRE
    drivers = pd.DataFrame({"seasonal": pd.Series(
        np.cos(2 * np.pi * td) + rng.normal(0, .1, len(days)), index=days)})
    # serie saisonniere avec un bruit INDEPENDANT : aucun lien causal
    series = pd.DataFrame({"date": dates,
                           "disp_mm": np.cos(2 * np.pi * (ts - 0.15))
                           + rng.normal(0, .1, len(dates))})
    naive = lag_scan(series, drivers, max_lag_days=90, step=6)
    deseas = lag_scan(series, drivers, max_lag_days=90, step=6,
                      deseasonalize=True)
    # mesure : |r| passe de 0.84 (illusoire) a 0.24 -> l'illusion s'effondre
    assert abs(naive.iloc[0]["r"]) > 0.80, naive
    assert abs(deseas.iloc[0]["r"]) < 0.40, deseas
    assert abs(deseas.iloc[0]["r"]) < 0.5 * abs(naive.iloc[0]["r"]), (naive, deseas)


def test_deseasonalize_keeps_real_anomaly_coupling():
    """A l'inverse, un couplage REEL sur les ANOMALIES doit survivre."""
    dates = pd.date_range("2022-01-01", periods=90, freq="12D")
    days = pd.date_range("2021-10-01", "2026-01-01", freq="1D")
    rng = np.random.default_rng(4)
    td = (days - days[0]).days.values / 365.25
    anom = pd.Series(rng.normal(0, 1, len(days)), index=days).rolling(
        30, min_periods=1).mean()
    drv = pd.Series(np.cos(2 * np.pi * td), index=days) + anom
    drivers = pd.DataFrame({"driver": drv})
    obs = drv.reindex(dates, method="nearest").values
    series = pd.DataFrame({"date": dates,
                           "disp_mm": obs + rng.normal(0, .05, len(dates))})
    deseas = lag_scan(series, drivers, max_lag_days=30, step=6,
                      deseasonalize=True)
    assert abs(deseas.iloc[0]["r"]) > 0.7, deseas      # le vrai lien survit


def test_driver_pvalues_extremes():
    obs = pd.DataFrame([{"driver": "d", "r": 0.9, "lag_days": 0}])
    weak = pd.DataFrame({"trial": range(40), "driver": "d",
                         "r": np.linspace(-0.2, 0.2, 40)})
    assert driver_pvalues(obs, weak).iloc[0]["p_value"] < 0.05
    strong = pd.DataFrame({"trial": range(40), "driver": "d",
                           "r": np.linspace(0.9, 0.99, 40)})
    assert driver_pvalues(obs, strong).iloc[0]["p_value"] > 0.9


def test_zone_areas_and_site_check():
    """Surfaces en hectares + ecart a la surface connue (validation objective)."""
    tmpl = _tmpl(); zones = _zones(tmpl)
    df = zone_areas(zones, tmpl, expected_site_ha=3.2)
    a = df.set_index("zone")
    assert a.loc["A", "n_px"] == 16
    # 16 px de 40x40 m = 25600 m2 = 2.56 ha
    assert abs(a.loc["A", "area_ha"] - 2.56) < 1e-6, a
    assert abs(df.attrs["inside_ha"] - 3.20) < 1e-6           # A+B
    assert abs(df.attrs["area_error_pct"]) < 1e-6             # coincide


def test_label_array_and_field_table():
    tmpl = _tmpl(); zones = _zones(tmpl)
    lab = zone_label_array(zones, tmpl)
    assert set(np.unique(lab.values[np.isfinite(lab.values)])) == {1., 2., 3., 4.}
    fld = xr.DataArray(np.random.default_rng(1).normal(0, 1, (NY, NX)),
                       coords=tmpl.coords, dims=tmpl.dims)
    assert len(zone_field_table(fld, zones)) == 4


def test_save_figure_writes_png():
    """L'export doit ecrire un PNG non vide, au nom exact demande."""
    if not HAS_MPL:
        print("  (matplotlib absent -> test d'export saute)")
        return
    import tempfile
    from pathlib import Path

    import matplotlib.pyplot as plt

    from insar_wetlands.zone_viz import save_figure

    d = Path(tempfile.mkdtemp()) / "figs"       # dossier inexistant : doit etre cree
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])
    paths = save_figure(fig, "F99_test", d, dpi=72)
    assert len(paths) == 1 and paths[0].name == "F99_test.png", paths
    assert paths[0].stat().st_size > 1000, "PNG suspicieusement petit"


def test_plots_run():
    """Fumee : les traces s'executent sans erreur (saute si matplotlib absent)."""
    if not HAS_MPL:
        print("  (matplotlib absent -> test de trace saute)")
        return
    from insar_wetlands.zone_viz import (plot_zone_distributions, plot_zone_map,
                                         plot_zones_over_field)
    tmpl = _tmpl(); zones = _zones(tmpl)
    fld = xr.DataArray(np.random.default_rng(1).normal(0, 1, (NY, NX)),
                       coords=tmpl.coords, dims=tmpl.dims)
    plot_zone_map(zones, tmpl)
    plot_zones_over_field(fld, zones)
    plot_zone_distributions(fld, zones)


if __name__ == "__main__":
    test_lag_scan_finds_driver_and_lag()
    test_detrend_kills_spurious_trend_correlation()
    test_deseasonalize_kills_shared_annual_cycle()
    test_deseasonalize_keeps_real_anomaly_coupling()
    test_driver_pvalues_extremes()
    test_zone_areas_and_site_check()
    test_label_array_and_field_table()
    test_save_figure_writes_png()
    test_plots_run()
    print("ALL PHASE-I TESTS PASSED")

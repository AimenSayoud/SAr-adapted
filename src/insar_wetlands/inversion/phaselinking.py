"""Phase-linking / DS-InSAR LÉGER (EVD) — sans ISCE, sans SLC, sans corégistration.

Le phase-linking estime, par pixel, l'historique de phase optimal a partir de la
MATRICE DE COHÉRENCE COMPLEXE N×N (N = nombre de dates), dont l'entrée (i,j) est
la cohérence × exp(i·phase interférométrique) de la paire (i,j). Or nos produits
HyP3 (phase enroulée + cohérence) SONT ces entrées, et HyP3 a déjà corégistré les
paires -> on peut faire du phase-linking en pur numpy, sans ISCE/MiaplPy.

Estimateur EVD (Fornaro/Ansari) : l'historique de phase est la phase du vecteur
propre de plus grande valeur propre de la matrice de cohérence. C'est le maximum
de vraisemblance sous modèle gaussien circulaire ; il utilise la PONDÉRATION par
cohérence de TOUTES les paires simultanément, la ou notre ISBAS résout paire par
paire. Qualité = « temporal coherence » (ajustement de l'historique estimé aux
interférogrammes observés) — l'analogue de MiaplPy.temporalCoherence.

Limite : notre réseau est SPARSE (sous-ensemble de paires) ; les entrées
manquantes sont mises a 0 (elles ne contraignent pas). C'est une version
allégée du phase-linking plein (qui utiliserait toutes les paires), mais elle
teste directement H1 : la pondération par matrice de cohérence récupère-t-elle
un signal la ou le WLS échoue ?
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import xarray as xr

from ..stack import dates_from_pairs
from .isbas import PHASE_TO_MM


def wrap_unw(unw: xr.DataArray) -> xr.DataArray:
    """Ré-enroule une phase HyP3 DÉROULÉE en phase enroulée [-π, π].

    Le phase-linking opère sur la phase ENROULÉE (l'observation brute). Comme
    nos crops ne contiennent que unw_phase, on la ré-enroule : angle(exp(i·φ))
    annule les cycles ajoutés par le déroulement et redonne l'observation
    d'origine (les erreurs de déroulement, multiples de 2π, disparaissent)."""
    return xr.apply_ufunc(lambda a: np.angle(np.exp(1j * a)), unw,
                          keep_attrs=True)


def zone_tcoh_summary(ds: xr.Dataset, zones: dict) -> "pd.DataFrame":
    """Résumé par zone de la temporal_coherence du phase-linking EVD.

    Colonnes : zone, n_px, tcoh_median, tcoh_p25, tcoh_p75, frac_ge_0p7.
    frac_ge_0p7 = part de pixels dépassant le seuil usuel de fiabilité
    MiaplPy (~0.7) — l'indicateur direct du test H1 par zone."""
    tc = ds.temporal_coherence
    rows = []
    for z in ("A", "B", "C", "D"):
        if z not in zones:
            continue
        v = tc.values[zones[z].values]
        v = v[np.isfinite(v)]
        if not v.size:
            continue
        rows.append({"zone": z, "n_px": int(v.size),
                     "tcoh_median": float(np.median(v)),
                     "tcoh_p25": float(np.percentile(v, 25)),
                     "tcoh_p75": float(np.percentile(v, 75)),
                     "frac_ge_0p7": float((v >= 0.7).mean())})
    return pd.DataFrame(rows)


def _pair_date_index(pairs: list[str]):
    dates = dates_from_pairs(pairs)
    di = {d.strftime("%Y%m%d"): i for i, d in enumerate(dates)}
    idx = [(di[p.split("_")[0]], di[p.split("_")[1]]) for p in pairs]
    return dates, idx


def evd_pixel(phi: np.ndarray, coh: np.ndarray, idx, n: int,
              coh_floor: float = 0.0) -> tuple[np.ndarray, float]:
    """Phase-linking EVD pour UN pixel.

    phi, coh : phase enroulée et cohérence des paires (longueur = n_pairs).
    Convention : idx[k] = (i, j) avec i = date la plus ancienne, j la plus
    récente, et phi[k] ≈ θ_j − θ_i (déplacement récent−ancien, comme l'ISBAS
    où phase = (disp_b − disp_a) / PHASE_TO_MM).

    Sous le modèle rang-1 Γ = c cᴴ avec c_k = exp(i·θ_k), on a
    Γ[i, j] = exp(i(θ_i − θ_j)) = exp(−i·φ). On remplit donc Γ[i, j] avec le
    CONJUGUÉ de exp(i·φ) pour que la phase du vecteur propre dominant redonne
    directement θ (même signe que l'historique ISBAS).

    Retourne (theta[n], temporal_coherence)."""
    gamma = np.zeros((n, n), dtype=np.complex128)
    np.fill_diagonal(gamma, 1.0)
    ok = np.isfinite(phi) & np.isfinite(coh) & (coh > coh_floor)
    n_ok = 0
    for k, (i, j) in enumerate(idx):
        if not ok[k]:
            continue
        g = coh[k] * np.exp(-1j * phi[k])   # Γ[i,j] = coh·exp(i(θ_i−θ_j))
        gamma[i, j] = g
        gamma[j, i] = np.conj(g)
        n_ok += 1
    if n_ok < n - 1:                      # sous-déterminé
        return np.full(n, np.nan), np.nan
    # vecteur propre de plus grande valeur propre
    w, v = np.linalg.eigh(gamma)
    vec = v[:, -1]
    theta = np.angle(vec)
    theta = theta - theta[0]             # référence a la 1re date
    # temporal coherence : ajustement de theta aux interférogrammes observés
    num = 0.0 + 0.0j
    den = 0.0
    for k, (i, j) in enumerate(idx):
        if not ok[k]:
            continue
        num += coh[k] * np.exp(1j * (phi[k] - (theta[j] - theta[i])))
        den += coh[k]
    tcoh = float(abs(num) / den) if den > 0 else np.nan
    return theta, tcoh


def evd_phase_linking(wrapped: xr.DataArray, corr: xr.DataArray,
                      ref_yx: tuple[float, float] | None = None,
                      coh_floor: float = 0.0, aoi: xr.DataArray | None = None) -> xr.Dataset:
    """Phase-linking EVD sur tout le stack (ou l'AOI si fournie).

    wrapped, corr : (pair, y, x) — phase ENROULÉE et cohérence des paires.
    Retourne un Dataset : temporal_coherence(y,x), velocity_mm_yr(y,x),
    displacement_mm(time,y,x). La qualité (temporal_coherence) est LE résultat :
    élevée => le phase-linking a récupéré un historique cohérent.
    """
    pairs = [str(p) for p in wrapped.pair.values]
    dates, idx = _pair_date_index(pairs)
    n = len(dates)
    ph = wrapped.values.astype("float32")     # (n_pairs, ny, nx)
    ch = corr.values.astype("float32")
    ny, nx = ph.shape[1:]

    tcoh = np.full((ny, nx), np.nan, "float32")
    disp = np.full((n, ny, nx), np.nan, "float32")
    sel = aoi.values if aoi is not None else np.ones((ny, nx), bool)
    for iy in range(ny):
        for ix in range(nx):
            if not sel[iy, ix]:
                continue
            theta, tc = evd_pixel(ph[:, iy, ix], ch[:, iy, ix], idx, n, coh_floor)
            tcoh[iy, ix] = tc
            if np.isfinite(tc):
                disp[:, iy, ix] = theta * PHASE_TO_MM

    # référencement optionnel (soustraire l'historique d'un pixel stable)
    if ref_yx is not None:
        ry = int(np.argmin(np.abs(wrapped.y.values - ref_yx[0])))
        rx = int(np.argmin(np.abs(wrapped.x.values - ref_yx[1])))
        ref_hist = disp[:, ry, rx]
        if np.all(np.isfinite(ref_hist)):
            disp = disp - ref_hist[:, None, None]

    # vitesse (régression linéaire de l'historique). Les pixels hors `aoi`
    # sont tout-NaN sur l'axe temps -> nanmean déclenche un RuntimeWarning
    # inoffensif (ils sont masqués juste après) : on le tait.
    t_years = np.array([(d - dates[0]).days / 365.25 for d in dates])
    tc_mean = t_years - t_years.mean()
    denom = float((tc_mean ** 2).sum())
    with np.errstate(invalid="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        disp_c = disp - np.nanmean(disp, axis=0)
    vel = np.tensordot(tc_mean, np.nan_to_num(disp_c), axes=(0, 0)) / denom
    vel[np.isnan(tcoh)] = np.nan

    return xr.Dataset(
        {
            "temporal_coherence": (("y", "x"), tcoh),
            "velocity_mm_yr": (("y", "x"), vel.astype("float32")),
            "displacement_mm": (("time", "y", "x"), disp),
        },
        coords={"time": dates, "y": wrapped.y, "x": wrapped.x},
    )

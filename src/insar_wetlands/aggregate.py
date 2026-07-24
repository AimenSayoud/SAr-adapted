"""Phase G — Extraire le signal du bruit par AGRÉGATION SPATIALE (super-pixel).

Changement d'OBSERVABLE, pas énième inversion. Les Phases 8-E2 ont toutes
cherché une carte **pixel par pixel** ; six estimateurs différents échouent. Mais
l'écart-type de phase d'un pixel a γ=0.4 est ~1.5 rad, alors que la moyenne
complexe sur N pixels le divise par ~sqrt(N_eff) : sur les 499 px du tapis,
~0.07 rad (~0.3 mm), et même avec N_eff=50 on est a ~1 mm — **très en-dessous de
la respiration de 10-40 mm recherchée**. Le signal n'est pas sous le plancher de
bruit ; il est sous le plancher de bruit *par pixel*.

Hypothèse physique qui l'autorise : le tapis flottant est **une unité
hydrologique** — il respire en bloc. Moyenner ne détruit donc pas le signal
(contrairement a un champ de déformation hétérogène), ça tue la composante
aléatoire et garde le mode commun.

Deux observables, deux rôles :

1. **|R| (module du vecteur résultant)** — `zone_phasor`. Test FALSIFIABLE de
   l'existence d'une phase commune : phases aléatoires => |R| ~ 1/sqrt(N_eff)
   (~0.045 pour N=499) ; phase commune => |R| >> ça. Calculé sur la phase
   ENROULÉE (le phaseur est insensible aux erreurs de déroulement).

2. **Double différence agrégée A−C** — `double_difference`. Les deux zones sont
   a ~1 km : même écran atmosphérique, même baseline (même paire) -> l'atmosphère
   et l'orbite s'annulent, il ne reste que le mouvement DIFFÉRENTIEL du tapis par
   rapport au sol stable. C'est exactement le signal de respiration.

Et un discriminateur physique du mécanisme :

3. **Biais de phase de fermeture** — `closure_bias_by_zone`. Un déplacement (même
   non rigide) est cohérent entre dates -> le triplet ferme a ZÉRO. Une variation
   diélectrique / de profondeur de pénétration (tourbe saturée) introduit un biais
   systématiquement NON NUL (De Zan et al. 2015 ; Ansari et al. 2021). C'est le
   test qui sépare (a) diélectrique de (b) micro-mouvement non rigide, question
   restée ouverte en Phase E2. NB : distinct de `isbas.phase_closure`, qui ne
   cherche que les erreurs de déroulement (multiples de 2*pi).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import xarray as xr

from .inversion.isbas import PHASE_TO_MM
from .stack import dates_from_pairs


def _wrap(x):
    return np.angle(np.exp(1j * np.asarray(x)))


def effective_looks(mask: xr.DataArray, corr_len_px: float = 2.0) -> float:
    """N_eff = N / (facteur de corrélation spatiale).

    Les pixels voisins ne sont PAS indépendants (multilooking + corrélation
    naturelle). On divise par l'aire d'un patch de corrélation (corr_len_px^2),
    estimation CONSERVATRICE : sous-estimer N_eff sur-estime le plancher de
    bruit, donc rend le test plus sévère (jamais optimiste)."""
    n = float(mask.sum())
    return max(1.0, n / max(1.0, corr_len_px ** 2))


def zone_phasor(wrapped: xr.DataArray, corr: xr.DataArray, zones: dict,
                zone_names=("A", "B", "C", "D")) -> pd.DataFrame:
    """Phaseur complexe agrégé par zone et par paire (le « super-pixel »).

    Pour chaque paire k et zone Z : R = Σ w·exp(i·phi) / Σ w, avec w = cohérence
    (pondération par qualité). Retourne un DataFrame long :
    zone, pair, dt_days, R_abs, R_phase_rad, n_px, n_eff, noise_floor.

    **Lecture** : `R_abs` >> `noise_floor` (=1/sqrt(n_eff)) => il existe une
    phase COMMUNE a la zone, même si aucun pixel n'est individuellement
    exploitable. C'est le test central de la Phase G.
    """
    pairs = [str(p) for p in wrapped.pair.values]
    ph = wrapped.values
    ch = corr.values
    rows = []
    masks = {z: zones[z].values for z in zone_names if z in zones}
    n_eff = {z: effective_looks(zones[z]) for z in masks}
    for k, p in enumerate(pairs):
        a, b = p.split("_")
        dt = (pd.Timestamp(b) - pd.Timestamp(a)).days
        pk, ck = ph[k], ch[k]
        for z, m in masks.items():
            ok = m & np.isfinite(pk) & np.isfinite(ck)
            if ok.sum() < 5:
                continue
            w = ck[ok]
            R = np.sum(w * np.exp(1j * pk[ok])) / np.sum(w)
            rows.append({
                "zone": z, "pair": p, "dt_days": dt,
                "R_abs": float(abs(R)), "R_phase_rad": float(np.angle(R)),
                "n_px": int(ok.sum()), "n_eff": round(n_eff[z], 1),
                "noise_floor": round(1.0 / np.sqrt(n_eff[z]), 4),
            })
    return pd.DataFrame(rows)


def phasor_verdict(df: pd.DataFrame) -> pd.DataFrame:
    """Résumé par zone : |R| médian vs plancher de bruit 1/sqrt(N_eff).

    `ratio` > 1 => phase commune détectée au-dessus du hasard ; `frac_above` =
    part des paires dépassant le plancher."""
    rows = []
    for z, g in df.groupby("zone"):
        floor = float(g["noise_floor"].iloc[0])
        med = float(g["R_abs"].median())
        rows.append({"zone": z, "n_pairs": len(g), "n_eff": float(g["n_eff"].iloc[0]),
                     "noise_floor": floor, "R_abs_median": round(med, 4),
                     "ratio_to_floor": round(med / floor, 2),
                     "frac_above_floor": round(float((g["R_abs"] > floor).mean()), 3)})
    return pd.DataFrame(rows).sort_values("ratio_to_floor", ascending=False)


def double_difference(wrapped: xr.DataArray, corr: xr.DataArray, zones: dict,
                      target: str = "A", reference: str = "C") -> pd.DataFrame:
    """Phase agrégée DOUBLE-DIFFÉRENCIÉE target − reference, par paire.

    Annule l'atmosphère et l'orbite (mêmes paire et échelle spatiale ~1 km).
    Retourne : pair, dt_days, ddphase_rad, ddisp_mm, weight, R_abs_target,
    R_abs_reference. `ddisp_mm` = déplacement LOS différentiel du tapis par
    rapport au sol stable, pour CETTE paire.

    ATTENTION déroulement : la phase agrégée est ENROULÉE dans (-pi, pi]. Une
    respiration de 40 mm vaut ~9 rad -> elle s'enroulerait. Cette sortie est donc
    fiable pour les DIFFÉRENCES faibles (|ddisp| << 13.9 mm = pi*|PHASE_TO_MM|) ;
    au-dela il faut passer par `invert_aggregate` sur la phase déroulée.
    """
    ph = zone_phasor(wrapped, corr, zones, zone_names=(target, reference))
    t = ph[ph.zone == target].set_index("pair")
    r = ph[ph.zone == reference].set_index("pair")
    common = t.index.intersection(r.index)
    dd = _wrap(t.loc[common, "R_phase_rad"].values
               - r.loc[common, "R_phase_rad"].values)
    return pd.DataFrame({
        "pair": common,
        "dt_days": t.loc[common, "dt_days"].values,
        "ddphase_rad": dd,
        "ddisp_mm": dd * PHASE_TO_MM,
        # poids = qualité conjointe des deux agrégats
        "weight": t.loc[common, "R_abs"].values * r.loc[common, "R_abs"].values,
        "R_abs_target": t.loc[common, "R_abs"].values,
        "R_abs_reference": r.loc[common, "R_abs"].values,
    }).reset_index(drop=True)


def aggregate_unwrapped(unw: xr.DataArray, corr: xr.DataArray, zones: dict,
                        target: str = "A", reference: str = "C") -> pd.DataFrame:
    """Double différence sur la phase DÉROULÉE (pas d'ambiguïté 2*pi).

    Moyenne pondérée par cohérence de la phase déroulée HyP3 dans chaque zone,
    puis différence. Complémentaire de `double_difference` : pas de repliement,
    mais sensible aux erreurs de déroulement des pixels peu cohérents (a croiser
    avec |R| et le biais de fermeture)."""
    pairs = [str(p) for p in unw.pair.values]
    u, c = unw.values, corr.values
    mt, mr = zones[target].values, zones[reference].values
    rows = []
    for k, p in enumerate(pairs):
        a, b = p.split("_")
        uk, ck = u[k], c[k]

        def wmean(m):
            ok = m & np.isfinite(uk) & np.isfinite(ck)
            if ok.sum() < 5:
                return np.nan, 0.0
            w = ck[ok]
            return float(np.sum(w * uk[ok]) / np.sum(w)), float(w.mean())

        pt, wt = wmean(mt)
        pr, wr = wmean(mr)
        if not (np.isfinite(pt) and np.isfinite(pr)):
            continue
        rows.append({"pair": p, "dt_days": (pd.Timestamp(b) - pd.Timestamp(a)).days,
                     "ddphase_rad": pt - pr, "ddisp_mm": (pt - pr) * PHASE_TO_MM,
                     "weight": wt * wr})
    return pd.DataFrame(rows)


def invert_aggregate(dd: pd.DataFrame) -> pd.DataFrame:
    """SBAS sur le super-pixel : série temporelle depuis les doubles différences.

    UN seul « pixel » (le tapis entier) : ~356 observations pour ~90 inconnues,
    massivement surdéterminé — la ou l'inversion par pixel était au plancher de
    bruit. Moindres carrés pondérés (poids = qualité conjointe), pseudo-inverse
    (réseau possiblement non connexe). Retourne : date, disp_mm (relatif a la 1re
    date), et la vitesse en attribut `.attrs['velocity_mm_yr']`.
    """
    from .network import design_matrix

    pairs = list(dd["pair"])
    dates = dates_from_pairs(pairs)
    A = design_matrix(pairs, dates)
    y = dd["ddphase_rad"].values * PHASE_TO_MM      # -> mm LOS
    w = np.sqrt(np.clip(dd["weight"].values, 1e-6, None))
    x = np.linalg.pinv(A * w[:, None]) @ (y * w)    # incréments (n_dates-1)
    disp = np.concatenate([[0.0], np.cumsum(x)])
    out = pd.DataFrame({"date": dates, "disp_mm": disp})
    t = (dates - dates[0]).days.values / 365.25
    tc = t - t.mean()
    out.attrs["velocity_mm_yr"] = float((tc * (disp - disp.mean())).sum()
                                        / (tc ** 2).sum())
    out.attrs["n_obs"] = len(pairs)
    out.attrs["n_unknowns"] = len(dates) - 1
    return out


def closure_bias_by_zone(wrapped: xr.DataArray, zones: dict,
                         max_triplets: int = 300,
                         zone_names=("A", "B", "C", "D")) -> pd.DataFrame:
    """Biais de phase de FERMETURE par zone — discriminateur du mécanisme.

    Pour chaque triplet (i<j<k) présent dans le réseau, la fermeture agrégée est
    arg(<exp(i*phi_ij)> * <exp(i*phi_jk)> * conj(<exp(i*phi_ik)>)) sur la zone.

    **Physique** : un DÉPLACEMENT (même non rigide) est cohérent entre dates ->
    fermeture nulle en moyenne. Une variation DIÉLECTRIQUE / de profondeur de
    pénétration (tourbe saturée, Sphagnum humide) brise la réciprocité et
    produit un biais SYSTÉMATIQUEMENT non nul (De Zan 2015 ; Ansari 2021).

    Retourne par zone : n_triplets, mean_closure_rad (biais signé),
    median_abs_rad, et `bias_significant` (test t sur la moyenne != 0).
    Calculé sur la phase ENROULÉE : insensible aux erreurs de déroulement,
    contrairement a `isbas.phase_closure` qui, lui, ne cherche QUE ces erreurs.
    """
    pairs = [str(p) for p in wrapped.pair.values]
    pset = {p: i for i, p in enumerate(pairs)}
    dates = dates_from_pairs(pairs)
    triplets = []
    for i, d1 in enumerate(dates):
        for d2 in dates[i + 1:]:
            for d3 in dates[i + 2:]:
                if d3 <= d2:
                    continue
                p12 = f"{d1:%Y%m%d}_{d2:%Y%m%d}"
                p23 = f"{d2:%Y%m%d}_{d3:%Y%m%d}"
                p13 = f"{d1:%Y%m%d}_{d3:%Y%m%d}"
                if p12 in pset and p23 in pset and p13 in pset:
                    triplets.append((pset[p12], pset[p23], pset[p13]))
        if len(triplets) >= max_triplets:
            break
    if not triplets:
        raise ValueError("aucun triplet ferme dans le reseau")

    ph = wrapped.values
    rows = []
    for z in zone_names:
        if z not in zones:
            continue
        m = zones[z].values
        vals = []
        for a, b, c in triplets:
            def agg(idx):
                v = ph[idx][m]
                v = v[np.isfinite(v)]
                return np.mean(np.exp(1j * v)) if v.size >= 5 else np.nan

            ra, rb, rc = agg(a), agg(b), agg(c)
            if not (np.isfinite(ra) and np.isfinite(rb) and np.isfinite(rc)):
                continue
            vals.append(np.angle(ra * rb * np.conj(rc)))
        if not vals:
            continue
        v = np.asarray(vals)
        se = v.std(ddof=1) / np.sqrt(v.size) if v.size > 1 else np.nan
        rows.append({"zone": z, "n_triplets": v.size,
                     "mean_closure_rad": float(v.mean()),
                     "median_abs_rad": float(np.median(np.abs(v))),
                     "se": float(se) if np.isfinite(se) else np.nan,
                     "bias_significant": bool(np.isfinite(se) and se > 0
                                              and abs(v.mean()) > 2 * se)})
    return pd.DataFrame(rows)

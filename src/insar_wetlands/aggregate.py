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


def correlation_length(field: xr.DataArray, mask: xr.DataArray,
                       max_lag_px: int = 12) -> float:
    """Longueur de corrélation spatiale EMPIRIQUE (en pixels), par autocorrélation.

    Estime la portée réelle de la corrélation spatiale au lieu de la supposer.
    On calcule l'autocorrélation du champ (démoyenné sur le masque) en fonction
    du décalage, et on retient le décalage ou elle passe sous 1/e ≈ 0.368.

    Motivation : le facteur 1/sqrt(N) suppose des pixels INDÉPENDANTS. Les
    pixels voisins ne le sont pas (multilooking, corrélation naturelle du
    couvert). Sans mesure, N_eff est une hypothèse — c'est le point le plus
    attaquable de l'argument d'agrégation."""
    v = np.where(mask.values, field.values, np.nan).astype(float)
    v = v - np.nanmean(v)
    ok = np.isfinite(v)
    v0 = np.where(ok, v, 0.0)
    denom = float((v0[ok] ** 2).sum())
    if denom <= 0:
        return 1.0
    for lag in range(1, max_lag_px + 1):
        # moyenne des autocorrélations selon x et y a ce décalage
        acs = []
        for a, b, m1, m2 in ((v0[:, lag:], v0[:, :-lag], ok[:, lag:], ok[:, :-lag]),
                             (v0[lag:, :], v0[:-lag, :], ok[lag:, :], ok[:-lag, :])):
            both = m1 & m2
            if both.sum() > 10:
                acs.append(float((a[both] * b[both]).sum())
                           / max(denom * both.sum() / max(ok.sum(), 1), 1e-12))
        if acs and np.mean(acs) < np.exp(-1.0):
            return float(lag)
    return float(max_lag_px)


def effective_looks(mask: xr.DataArray, corr_len_px: float = 2.0,
                    field: xr.DataArray | None = None) -> float:
    """N_eff = N / (aire d'un patch de corrélation).

    Si `field` est fourni, la longueur de corrélation est **mesurée** sur les
    données (`correlation_length`) au lieu d'être supposée. Sinon on retombe sur
    `corr_len_px` (défaut 2.0), hypothèse CONSERVATRICE : sous-estimer N_eff
    sur-estime le plancher de bruit, donc durcit le test.

    ⚠️ Cette quantité n'intervient QUE dans le plancher indicatif de |R|. Les
    tests de significativité (amplitude saisonnière, corrélations) reposent sur
    des **nuls empiriques appariés en taille**, qui intègrent la corrélation
    spatiale réelle **sans aucune hypothèse sur N_eff**."""
    n = float(mask.sum())
    if field is not None:
        corr_len_px = correlation_length(field, mask)
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

        def wmean(m, uk=uk, ck=ck):
            # uk/ck bound as defaults: wmean is called inside this iteration,
            # so behaviour is unchanged, but the late-binding trap is removed.
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


def filter_pairs(dd: pd.DataFrame, max_dt_days: int | None = 60,
                 min_weight: float | None = None,
                 exclude_months: tuple | None = None) -> pd.DataFrame:
    """Retire les paires selon la baseline, le poids, ou la SAISON.

    - `max_dt_days` : les paires annuelles (~370 j) et bi-annuelles (~740 j)
      sont quasi certainement mal DÉROULÉES sur le tapis (±25 mm de dispersion
      observée) et injectent des sauts de 2*pi dans l'inversion agrégée.
    - `exclude_months` : **test de falsification**. Un manteau NEIGEUX ou le GEL
      affectent différemment tourbière et prairie et possèdent un cycle annuel :
      ils constituent donc une explication alternative du signal saisonnier.
      Retirer les paires dont l'une des dates tombe en hiver
      (ex. `exclude_months=(12, 1, 2)`) permet de vérifier que le signal ne
      provient PAS de l'hiver. S'il survit, neige et gel sont écartés.
    """
    out = dd
    if max_dt_days is not None:
        out = out[out["dt_days"] <= max_dt_days]
    if min_weight is not None:
        out = out[out["weight"] >= min_weight]
    if exclude_months:
        def _keep(p):
            a, b = str(p).split("_")
            return (int(a[4:6]) not in exclude_months
                    and int(b[4:6]) not in exclude_months)
        out = out[out["pair"].map(_keep)]
    return out.reset_index(drop=True)


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


def adjacent_null_zones(zones: xr.DataArray | dict, template: xr.DataArray,
                        ref: str = "D") -> dict:
    """Contrôle nul APPARIÉ SPATIALEMENT : deux moitiés ADJACENTES de `ref`.

    Un `argwhere` naïf coupe la zone en haut/bas de l'image — deux régions
    ÉLOIGNÉES, dont l'atmosphère ne s'annule PAS, alors que A et C sont
    voisines. Le nul était donc trop sévère (il mesurait le gradient
    atmosphérique, pas le bruit de la méthode).

    Ici on coupe par la MÉDIANE DE X (bande gauche / bande droite), et on ne
    garde qu'une bande centrale de part et d'autre de la coupure -> deux
    demi-zones **contiguës**, à la même échelle spatiale que A vs C.
    """
    m = (zones[ref] if isinstance(zones, dict) else zones).values
    yy, xx = np.where(m)
    xmid = np.median(xx)
    # largeur de bande = celle qui donne ~autant de px que la plus petite zone
    span = np.percentile(np.abs(xx - xmid), 40)
    left = m & (template.x.values[None, :] <= template.x.values[int(xmid)])
    right = m & (template.x.values[None, :] > template.x.values[int(xmid)])
    near = np.abs(np.arange(m.shape[1])[None, :] - xmid) <= max(2.0, span)
    mk = lambda a: xr.DataArray(a & near, coords=template.coords,
                                dims=template.dims)
    return {"A": mk(left), "C": mk(right)}


def seasonal_amplitude(series: pd.DataFrame, date_col: str = "date",
                       value_col: str = "disp_mm") -> dict:
    """Amplitude du cycle ANNUEL d'une série agrégée (+ tendance).

    **Pourquoi c'est LA bonne observable.** La « respiration » d'une tourbière
    est un gonflement/dégonflement SAISONNIER de 10-40 mm (Hrysiewicz 2024), pas
    une dérive linéaire. Tester une VITESSE (mm/an) sur un signal purement
    périodique donne ~0 par construction, même si la respiration est parfaitement
    mesurée : la vitesse n'a aucune puissance sur ce signal.

    Ajuste y = c + d*t + a*cos(2*pi*t) + b*sin(2*pi*t) (t en années) et retourne
    amplitude = sqrt(a^2+b^2), sa phase (jour du maximum), la tendance d, et le
    R^2 du terme saisonnier. À comparer a la même quantité sur le CONTRÔLE NUL :
    l'amplitude n'est un signal que si elle dépasse celle du nul.
    """
    d = pd.to_datetime(series[date_col])
    t = (d - d.iloc[0]).dt.days.values / 365.25
    y = series[value_col].values.astype(float)
    ok = np.isfinite(y)
    t, y = t[ok], y[ok]
    if t.size < 6:
        return {"amplitude_mm": np.nan, "phase_doy": np.nan,
                "trend_mm_yr": np.nan, "r2_seasonal": np.nan, "n": int(t.size)}
    M = np.column_stack([np.ones_like(t), t,
                         np.cos(2 * np.pi * t), np.sin(2 * np.pi * t)])
    beta, *_ = np.linalg.lstsq(M, y, rcond=None)
    resid = y - M @ beta
    # R^2 du SEUL terme saisonnier : gain par rapport au modèle constante+tendance
    M0 = M[:, :2]
    r0 = y - M0 @ np.linalg.lstsq(M0, y, rcond=None)[0]
    ss0 = float((r0 ** 2).sum())
    amp = float(np.hypot(beta[2], beta[3]))
    # jour de l'année du maximum du cosinus ajusté
    phase = float((np.arctan2(beta[3], beta[2]) / (2 * np.pi)) % 1.0 * 365.25)
    doy0 = int(pd.Timestamp(d.iloc[0]).dayofyear)
    return {"amplitude_mm": round(amp, 3),
            "phase_doy": round((doy0 + phase) % 365.25, 1),
            "trend_mm_yr": round(float(beta[1]), 3),
            "r2_seasonal": round(1.0 - float((resid ** 2).sum()) / ss0, 4)
            if ss0 > 0 else np.nan,
            "n": int(t.size)}


def _compact_blob(cand_yx: np.ndarray, seed_i: int, n: int) -> np.ndarray:
    """Les `n` pixels candidats les plus proches d'une graine -> tache COMPACTE."""
    d = np.hypot(cand_yx[:, 0] - cand_yx[seed_i, 0],
                 cand_yx[:, 1] - cand_yx[seed_i, 1])
    return np.argsort(d)[:n]


def matched_null_zones(zones: dict, template: xr.DataArray, n_target: int,
                       n_reference: int, ref: str = "D",
                       seed: int = 0) -> dict | None:
    """Contrôle nul APPARIÉ EN TAILLE : deux taches compactes adjacentes de `ref`,
    de **mêmes effectifs** que les zones réelles.

    Pourquoi c'est indispensable : le bruit d'un agrégat décroît en 1/sqrt(N).
    Un nul construit sur 2200 px alors que A n'en a que 499 a ~2x moins de bruit
    et **sous-estime donc le plancher** — ce qui fabrique de fausses détections.
    Le nul doit avoir exactement la taille des zones qu'il imite.
    """
    m = zones[ref].values
    yx = np.argwhere(m)
    need = n_target + n_reference
    if len(yx) < need:
        return None
    rng = np.random.default_rng(seed)
    seed_i = int(rng.integers(len(yx)))
    both = _compact_blob(yx, seed_i, need)          # tache compacte de taille need
    sub = yx[both]
    # coupe la tache en deux selon une direction ALÉATOIRE -> deux moitiés adjacentes
    th = rng.uniform(0, np.pi)
    proj = sub[:, 0] * np.cos(th) + sub[:, 1] * np.sin(th)
    order = np.argsort(proj)
    t_idx, r_idx = order[:n_target], order[n_target:n_target + n_reference]
    out = {}
    for name, idx in (("A", t_idx), ("C", r_idx)):
        a = np.zeros_like(m)
        a[sub[idx, 0], sub[idx, 1]] = True
        out[name] = xr.DataArray(a, coords=template.coords, dims=template.dims)
    return out


def null_distribution(unw: xr.DataArray, corr: xr.DataArray, zones: dict,
                      template: xr.DataArray, n_target: int, n_reference: int,
                      n_trials: int = 50, ref: str = "D",
                      max_dt_days: int | None = None) -> pd.DataFrame:
    """Distribution NULLE de l'amplitude saisonnière (et de la vitesse).

    Répète `n_trials` fois : deux taches compactes adjacentes de sol stable, de
    tailles identiques a celles des zones réelles -> inversion agrégée ->
    amplitude saisonnière. On obtient ainsi une **distribution** du plancher, et
    donc une **p-value empirique**, au lieu d'une comparaison a une seule
    réalisation (qui n'est pas un test).
    """
    rows = []
    for t in range(n_trials):
        zn = matched_null_zones(zones, template, n_target, n_reference,
                                ref=ref, seed=t)
        if zn is None:
            break
        try:
            dd = aggregate_unwrapped(unw, corr, zn, "A", "C")
            if max_dt_days is not None:
                dd = filter_pairs(dd, max_dt_days=max_dt_days)
            if len(dd) < 10:
                continue
            r = invert_aggregate(dd)
            s = seasonal_amplitude(r)
            rows.append({"trial": t, "amplitude_mm": s["amplitude_mm"],
                         "r2_seasonal": s["r2_seasonal"],
                         "velocity_mm_yr": r.attrs["velocity_mm_yr"]})
        except Exception:
            continue
    # colonnes garanties même si AUCUN tirage n'a abouti (ex. zone trop grande :
    # matched_null_zones exige n_target+n_reference pixels dans `ref`, ce que D
    # ne peut pas fournir pour elle-même) -> évite un KeyError en aval.
    return pd.DataFrame(rows, columns=["trial", "amplitude_mm", "r2_seasonal",
                                       "velocity_mm_yr"])


def empirical_pvalue(observed: float, null_values, tail: str = "greater") -> dict:
    """p-value empirique de `observed` contre une distribution nulle.

    p = (1 + #{nul >= observé}) / (1 + n) — l'ajout de 1 évite p=0, qui n'est
    jamais justifiable avec un nombre fini de tirages."""
    v = np.asarray([x for x in np.asarray(null_values, float) if np.isfinite(x)])
    if v.size == 0 or not np.isfinite(observed):
        return {"p_value": np.nan, "n_null": 0}
    k = int((v >= observed).sum()) if tail == "greater" else int((v <= observed).sum())
    return {"p_value": round((1 + k) / (1 + v.size), 4), "n_null": int(v.size),
            "null_median": round(float(np.median(v)), 4),
            "null_p95": round(float(np.percentile(v, 95)), 4),
            "observed": round(float(observed), 4)}


def seasonal_zone_scan(unw: xr.DataArray, corr: xr.DataArray, zones: dict,
                       template: xr.DataArray, reference: str = "C",
                       targets=("A", "B", "D"), n_trials: int = 100,
                       max_dt_days: int | None = None) -> pd.DataFrame:
    """Amplitude saisonnière de CHAQUE zone vs la référence, + p-value appariée.

    **Le test qui discriminerait mouvement vs humidité.** Une amplitude
    saisonnière n'est pas une preuve de MOUVEMENT : un cycle annuel d'humidité
    (profondeur de pénétration variable) produit le même signal de phase sans
    aucun déplacement de surface — 3 mm ne valent que 0.75 rad en bande C.

    Or le **lac (B) ne peut pas respirer mécaniquement**. S'il présente lui
    aussi une amplitude saisonnière comparable a celle du tapis, le signal est
    d'origine **diélectrique**. S'il est plat alors que le tapis oscille, le
    mouvement redevient l'explication la plus simple.

    Chaque zone est testée contre son PROPRE nul apparié en taille (le plancher
    dépend de N : comparer B, qui a peu de pixels, au nul de A serait trompeur).
    """
    rows = []
    for z in targets:
        if z not in zones or int(zones[z].sum()) < 20:
            continue
        n_t, n_r = int(zones[z].sum()), int(zones[reference].sum())
        # le réservoir de nuls doit contenir n_t + n_r pixels : une zone aussi
        # grande que le réservoir lui-même (cas de D vs D) est intestable.
        if int(zones.get("D", zones[reference]).sum()) < n_t + n_r:
            print(f"  ! zone {z}: réservoir de nuls trop petit "
                  f"({n_t}+{n_r} px requis) -> non testable")
            continue
        try:
            dd = aggregate_unwrapped(unw, corr, zones, target=z, reference=reference)
            if max_dt_days is not None:
                dd = filter_pairs(dd, max_dt_days=max_dt_days)
            s = seasonal_amplitude(invert_aggregate(dd))
            nulls = null_distribution(unw, corr, zones, template, n_t, n_r,
                                      n_trials=n_trials, max_dt_days=max_dt_days)
            pv = empirical_pvalue(s["amplitude_mm"], nulls["amplitude_mm"])
        except Exception as e:
            print(f"  ! zone {z}: {e}")
            continue
        rows.append({"zone": z, "n_px": n_t, "amplitude_mm": s["amplitude_mm"],
                     "phase_doy": s["phase_doy"], "r2_seasonal": s["r2_seasonal"],
                     "null_median": pv["null_median"], "null_p95": pv["null_p95"],
                     "p_value": pv["p_value"], "n_null": pv["n_null"]})
    return pd.DataFrame(rows)


def closure_bias_by_zone(wrapped: xr.DataArray, zones: dict,
                         max_triplets: int = 3000,
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
    # On énumère TOUS les triplets fermés puis on sous-échantillonne AU HASARD
    # si besoin. (Un `break` sur le compteur ne garderait que les premières
    # dates -> biais temporel.) Plus de triplets = SE plus petit = plus de
    # puissance pour détecter le biais : c'est le paramètre critique du test.
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
    if not triplets:
        raise ValueError("aucun triplet ferme dans le reseau")
    if len(triplets) > max_triplets:
        idx = np.random.default_rng(0).choice(len(triplets), max_triplets,
                                              replace=False)
        triplets = [triplets[i] for i in idx]

    ph = wrapped.values
    rows = []
    for z in zone_names:
        if z not in zones:
            continue
        m = zones[z].values
        vals = []
        for a, b, c in triplets:
            def agg(idx, m=m):  # bound as a default: called within this iteration, so behaviour is unchanged
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

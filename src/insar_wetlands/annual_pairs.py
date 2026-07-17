"""Approche 'seasonal-annual' (Ghezelayagh et al. 2024, Ecological Indicators
166:112305) : quelques paires a ~1 an d'ecart, MEME mois chaque annee, plutot
qu'un reseau SBAS dense multi-saisons.

Rationale (leur article, section 2.3.1, et notre propre experience Phase 8-9) :
  - Sur tourbiere en C-band, un ecart de 12-48 j accumule un signal de
    subsidence de l'ordre de 0.05-0.2 mm -> noye dans le bruit de phase,
    d'ou la necessite d'un reseau dense + inversion ponderee (SBAS/ISBAS).
  - Sur un ecart de ~365 j, le signal attendu (~10-20 mm/an pour une
    tourbiere degradee) est largement au-dessus du bruit -> UNE SEULE paire
    bien choisie suffit, sans reseau ni inversion.
  - Comparer le meme mois d'une annee sur l'autre (phenologie/etat de
    surface comparables) evite la decorrelation saisonniere : pas besoin de
    "reseau connexe", pas de bridging/phase_closure MintPy.

Cette approche est INDEPENDANTE du reste du pipeline SBAS/ISBAS (Phases
2-10) : elle reutilise seulement l'inventaire S1 (Phase 1) et les fonctions
de soumission/telechargement HyP3 (Phase 2), sur un tout petit nombre de
paires nouvelles.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import xarray as xr


def pick_seasonal_dates(inv: pd.DataFrame, target_month: int = 4,
                        target_day: int = 15) -> pd.DataFrame:
    """Pour chaque annee de l'inventaire, la date S1 la plus proche du
    (target_month, target_day) — ex: 15 avril, vegetation minimale.

    Retourne un DataFrame indexe par annee : date, granule, ecart_jours.
    """
    inv = inv.copy()
    inv["date"] = pd.to_datetime(inv["date"])
    inv["year"] = inv["date"].dt.year
    rows = []
    for year, g in inv.groupby("year"):
        target = pd.Timestamp(year=year, month=target_month, day=target_day)
        g = g.assign(gap=(g["date"] - target).abs())
        best = g.sort_values("gap").iloc[0]
        rows.append({"year": year, "date": best["date"],
                     "granule": best["granule"],
                     "gap_days": int(best["gap"].days)})
    return pd.DataFrame(rows).sort_values("year").reset_index(drop=True)


def gather_candidate_dates(inv: pd.DataFrame, target_month: int = 4,
                           target_day: int = 15, window_days: int = 15) -> pd.DataFrame:
    """TOUTES les dates S1 (pas seulement la plus proche) dans une fenetre de
    +/- window_days autour du (target_month, target_day), pour chaque annee.

    C'est le "pool" de ~5 candidats/mois de Ghezelayagh et al. 2024 (section
    2.3.1, cycle de 6 j -> jusqu'a 25 paires candidates par transition) —
    contrairement a `pick_seasonal_dates` qui ne retient que la date la plus
    proche sans verification de sa qualite atmospherique.
    """
    inv = inv.copy()
    inv["date"] = pd.to_datetime(inv["date"])
    inv["year"] = inv["date"].dt.year
    rows = []
    for year, g in inv.groupby("year"):
        target = pd.Timestamp(year=year, month=target_month, day=target_day)
        g = g.assign(gap=(g["date"] - target).dt.days)
        cand = g[g["gap"].abs() <= window_days]
        for _, r in cand.iterrows():
            rows.append({"year": year, "date": r["date"], "granule": r["granule"],
                        "gap_days": int(r["gap"])})
    return pd.DataFrame(rows).sort_values(["year", "date"]).reset_index(drop=True)


def score_candidates_by_era5(candidates: pd.DataFrame, era5_path,
                             lon: float, lat: float) -> pd.DataFrame:
    """Score de 'bruit atmospherique' par date candidate, a partir d'ERA5 deja
    telecharge (Phase 1) — remplace les releves manuels sur wunderground.com
    de l'article par une mesure quantitative et reproductible :
      - tcwv (vapeur d'eau colonne totale) a l'heure la plus proche de
        l'acquisition S1 : le retard de phase InSAR est ~ proportionnel au
        TCWV, donc c'est le proxy le plus direct du bruit qui nous interesse.
      - tp (precipitation cumulee sur les 24h precedentes) : penalise les
        episodes pluvieux, source de turbulence atmospherique a petite
        echelle spatiale, mal corrigee meme par difference de paire.
    Retourne `candidates` + colonnes tcwv, tp_24h, atmo_score (z-score
    combine ; plus bas = conditions plus favorables).
    """
    import xarray as xr

    with xr.open_dataset(era5_path) as ds:
        ds = ds.load()
    time_dim = "valid_time" if "valid_time" in ds.dims else "time"
    ds = ds.sortby(time_dim)
    pt = ds.sel(longitude=lon, latitude=lat, method="nearest")

    out = candidates.copy()
    tcwv_vals, tp_vals = [], []
    for _, r in out.iterrows():
        near = pt.sel({time_dim: r["date"]}, method="nearest")
        tcwv_vals.append(float(near["tcwv"].values))
        win = pt.sel({time_dim: slice(r["date"] - pd.Timedelta("1D"), r["date"])})
        tp_vals.append(float(win["tp"].sum().values) if win.sizes[time_dim] else np.nan)
    out["tcwv"] = tcwv_vals
    out["tp_24h"] = tp_vals

    def _z(s: pd.Series) -> pd.Series:
        std = s.std()
        return (s - s.mean()) / std if std > 0 else s * 0.0

    out["atmo_score"] = (_z(out["tcwv"])
                         + _z(out["tp_24h"].fillna(out["tp_24h"].median())))
    return out


def select_optimal_annual_pairs(scored_candidates: pd.DataFrame,
                                max_gap_years: int = 2) -> pd.DataFrame:
    """Pour chaque transition (annee_i, annee_j), choisit — parmi TOUS les
    couples du pool ERA5-note — celui qui minimise le bruit atmospherique
    combine des deux dates, au lieu de prendre la seule date la plus proche
    du 15 avril sans verification (`pick_seasonal_dates`).

    Score de paire = |Δtcwv| (retard differentiel, qui contamine
    directement l'interferogramme puisque seule la DIFFERENCE de vapeur
    d'eau entre les deux dates survit a la formation de l'interferogramme)
    + 0.5*(atmo_score_ref + atmo_score_sec) (penalise en plus les episodes
    individuellement pluvieux/turbulents). Pas de terme de baseline
    perpendiculaire : structurellement < 150 m sur ces bursts S1 (orbite
    tube, cf. hyp3/jobs.py), donc non discriminant ici — contrairement au
    Sentinel-1 pleine-scene de l'article, qui doit encore l'optimiser.

    C'est l'implementation quantitative de la section 2.3.1 de Ghezelayagh
    et al. 2024 (pool de candidats + optimisation meteo explicite), que
    `pick_seasonal_dates`/`build_annual_pair_list` ne faisaient pas.
    """
    rows = []
    years = sorted(scored_candidates["year"].unique())
    for i, yi in enumerate(years):
        for yj in years[i + 1:]:
            if yj - yi > max_gap_years:
                continue
            pool_i = scored_candidates[scored_candidates.year == yi]
            pool_j = scored_candidates[scored_candidates.year == yj]
            best = None
            for _, ri in pool_i.iterrows():
                for _, rj in pool_j.iterrows():
                    d_tcwv = abs(ri["tcwv"] - rj["tcwv"])
                    combined = d_tcwv + 0.5 * (ri["atmo_score"] + rj["atmo_score"])
                    if best is None or combined < best["combined_score"]:
                        best = {"ref_date": ri["date"], "sec_date": rj["date"],
                               "d_tcwv": d_tcwv, "combined_score": combined}
            if best is not None:
                a, b = best["ref_date"], best["sec_date"]
                rows.append({
                    "ref_date": a, "sec_date": b, "dt_days": (b - a).days,
                    "pair": f"{a:%Y%m%d}_{b:%Y%m%d}",
                    "d_tcwv_kg_m2": round(best["d_tcwv"], 3),
                    "combined_atmo_score": round(best["combined_score"], 3),
                })
    return pd.DataFrame(rows)


def select_topk_annual_pairs(scored_candidates: pd.DataFrame,
                             max_gap_years: int = 2, k: int = 3) -> pd.DataFrame:
    """Comme `select_optimal_annual_pairs`, mais retient les K meilleures
    paires PAR TRANSITION au lieu d'une seule — notre amelioration de
    redondance par rapport a l'article (qui n'a qu'UNE paire par transition,
    Table 1) : la mediane d'ensemble sur k paires ecrase le bruit
    atmospherique residuel et les sauts de deroulement isoles, la ou une
    paire unique porte tout son bruit sans recours.
    """
    rows = []
    years = sorted(scored_candidates["year"].unique())
    for i, yi in enumerate(years):
        for yj in years[i + 1:]:
            if yj - yi > max_gap_years:
                continue
            pool_i = scored_candidates[scored_candidates.year == yi]
            pool_j = scored_candidates[scored_candidates.year == yj]
            combos = []
            for _, ri in pool_i.iterrows():
                for _, rj in pool_j.iterrows():
                    d_tcwv = abs(ri["tcwv"] - rj["tcwv"])
                    combined = d_tcwv + 0.5 * (ri["atmo_score"] + rj["atmo_score"])
                    combos.append({"ref_date": ri["date"], "sec_date": rj["date"],
                                  "d_tcwv": d_tcwv, "combined_score": combined})
            combos.sort(key=lambda c: c["combined_score"])
            seen_dates = set()
            kept = 0
            for c in combos:
                if kept >= k:
                    break
                # diversite : eviter de reutiliser exactement le meme couple ;
                # autoriser le partage d'UNE date mais pas des deux
                key = (c["ref_date"], c["sec_date"])
                if key in seen_dates:
                    continue
                seen_dates.add(key)
                a, b = c["ref_date"], c["sec_date"]
                rows.append({
                    "transition": f"{yi}-{yj}", "rank": kept + 1,
                    "ref_date": a, "sec_date": b, "dt_days": (b - a).days,
                    "pair": f"{a:%Y%m%d}_{b:%Y%m%d}",
                    "d_tcwv_kg_m2": round(c["d_tcwv"], 3),
                    "combined_atmo_score": round(c["combined_score"], 3),
                })
                kept += 1
    return pd.DataFrame(rows)


def ensemble_rate(rates: dict[str, dict], transitions: pd.DataFrame,
                  corr_min: float = 0.20) -> xr.Dataset:
    """Taux vertical d'ensemble (mm/an) : mediane pixel a pixel sur toutes
    les paires, ponderee implicitement par le masquage de coherence, avec
    incertitude = ecart interquartile inter-paires — la ou l'article ne
    fournit qu'un RMSE global contre le terrain (qu'on n'a pas ici).

    `rates` : sortie de compute_pair_rate par paire ; `transitions` : le
    DataFrame de select_topk_annual_pairs (colonnes pair, transition).
    """
    per_pair = []
    for _, row in transitions.iterrows():
        p = row["pair"]
        if p not in rates:
            continue
        r = rates[p]
        rate = r["rate_mm_yr"].where(r["corr"] >= corr_min)
        per_pair.append(rate.assign_coords(pair=p))
    if not per_pair:
        raise ValueError("aucune paire disponible pour l'ensemble")
    stack = xr.concat(per_pair, dim="pair")
    med = stack.median("pair", skipna=True)
    q75 = stack.quantile(0.75, dim="pair", skipna=True).drop_vars("quantile")
    q25 = stack.quantile(0.25, dim="pair", skipna=True).drop_vars("quantile")
    n = stack.notnull().sum("pair")
    med.attrs = {"long_name": "Taux vertical median d'ensemble", "units": "mm/an",
                 "n_pairs_total": int(stack.sizes["pair"]), "corr_min": corr_min}
    iqr = (q75 - q25)
    iqr.attrs = {"long_name": "Ecart interquartile inter-paires", "units": "mm/an"}
    return xr.Dataset({"rate_mm_yr": med, "iqr_mm_yr": iqr,
                       "n_pairs_used": n.astype("int16")})


def build_annual_pair_list(seasonal_dates: pd.DataFrame,
                           max_gap_years: int = 2) -> pd.DataFrame:
    """Toutes les paires (annee_i, annee_j) avec 1 <= j-i <= max_gap_years.

    Retourne un DataFrame compatible avec hyp3.jobs.submit_pairs
    (colonnes: ref_date, sec_date, dt_days, pair).
    """
    rows = []
    years = seasonal_dates.sort_values("year").reset_index(drop=True)
    for i in range(len(years)):
        for j in range(i + 1, len(years)):
            gap_y = years.loc[j, "year"] - years.loc[i, "year"]
            if gap_y > max_gap_years:
                continue
            a, b = years.loc[i, "date"], years.loc[j, "date"]
            rows.append({
                "ref_date": a, "sec_date": b,
                "dt_days": (b - a).days,
                "pair": f"{a:%Y%m%d}_{b:%Y%m%d}",
            })
    return pd.DataFrame(rows)


def compute_pair_rate(cropped_root, pair: str, ref_yx: tuple[float, float],
                      lv_theta: xr.DataArray) -> dict:
    """Taux de deplacement vertical (mm/an) pour UNE paire annuelle,
    reference au meme pixel Classe A que le reste du pipeline.

    Pas de reseau, pas d'inversion : phase relative au pixel de reference,
    conversion directe en mm, projection LOS->vertical, division par dt.
    """
    from .geometry import incidence_angle
    from .inversion.isbas import PHASE_TO_MM
    from .stack import load_layer

    unw = load_layer(cropped_root, "unw_phase", [pair]).isel(pair=0)
    corr = load_layer(cropped_root, "corr", [pair]).isel(pair=0)

    ref_phase = float(unw.sel(y=ref_yx[0], x=ref_yx[1], method="nearest"))
    rel_phase = unw - ref_phase
    los_mm = rel_phase * PHASE_TO_MM

    inc = incidence_angle(lv_theta)
    if inc.sizes.get("y") != los_mm.sizes.get("y"):
        inc = inc.interp(y=los_mm.y, x=los_mm.x, method="nearest")
    vert_mm = los_mm / np.cos(inc)

    a, b = pd.Timestamp(pair[:8]), pd.Timestamp(pair[9:])
    dt_years = (b - a).days / 365.25

    rate_mm_yr = (vert_mm / dt_years).rename("rate_mm_yr")
    rate_mm_yr.attrs = {"long_name": "Taux de deplacement vertical",
                        "units": "mm/an", "pair": pair,
                        "dt_years": round(dt_years, 3)}
    return {"pair": pair, "dt_years": dt_years, "rate_mm_yr": rate_mm_yr,
           "vert_mm": vert_mm, "corr": corr}


def joint_annual_trend(cropped_root, annual_pairs: list[str],
                       ref_yx: tuple[float, float], lv_theta: xr.DataArray,
                       gamma_min: float = 0.15) -> xr.Dataset:
    """Tendance UNIQUE ajustee conjointement sur toutes les paires
    annuelles disponibles, au lieu de faire confiance a une seule paire.

    Avec 2-3 dates seulement, une paire isolee porte tout son bruit
    atmospherique et son risque de saut de deroulement sans rien pour se
    corriger (contrairement a l'article source qui moyenne sur 7 paires et
    ~96 000 ha). Reutilise l'inversion WLS de l'ISBAS (design_matrix +
    moindres carres) sur ce tout petit reseau : min_pairs bas (les dates
    manquent, pas les paires) mais la logique de ponderation par coherence
    est identique. Retourne une vitesse (mm/an) + erreur standard, projetee
    en vertical.
    """
    from .geometry import los_to_vertical
    from .inversion.isbas import invert_stack
    from .compare import fit_velocity
    from .stack import load_layer

    unw = load_layer(cropped_root, "unw_phase", annual_pairs)
    corr = load_layer(cropped_root, "corr", annual_pairs)
    dry = xr.ones_like(unw, dtype=bool)  # pas de masque d'eau ici (dates rares)

    isbas = invert_stack(unw, corr, dry, ref_yx,
                         min_pairs=1, gamma_min=gamma_min)
    d_vert = los_to_vertical(isbas.los_displacement_mm, lv_theta)
    vel = fit_velocity(d_vert)
    vel["n_valid_pairs"] = isbas.n_valid_pairs
    return vel


def deramp(field: xr.DataArray, fit_mask: xr.DataArray) -> xr.Dataset:
    """Retire une rampe planaire (a*x + b*y + c) ajustee par moindres
    carres sur les pixels de `fit_mask` (typiquement classes A/B, loin de
    la zone d'interet -> pas de vraie deformation attendue la-bas).

    Une rampe orbitale/atmospherique residuelle non corrigee produit un
    degrade lisse et symetrique sur toute l'image, sans rapport avec la
    geometrie du site (contrairement a un vrai signal de subsidence,
    concentre autour de la zone humide) -> reconnaissable a l'oeil et
    corrigeable par un simple plan ajuste hors de la zone d'interet.
    """
    yy, xx = np.meshgrid(field.y.values, field.x.values, indexing="ij")
    valid = fit_mask.values & np.isfinite(field.values)
    if valid.sum() < 10:
        raise ValueError("pas assez de pixels de reference pour ajuster une rampe")
    A = np.column_stack([xx[valid], yy[valid], np.ones(valid.sum())])
    coef, *_ = np.linalg.lstsq(A, field.values[valid], rcond=None)
    ramp = coef[0] * xx + coef[1] * yy + coef[2]
    corrected = field - xr.DataArray(ramp, coords=field.coords, dims=field.dims)
    return xr.Dataset({"corrected": corrected,
                       "ramp": xr.DataArray(ramp, coords=field.coords, dims=field.dims)},
                      attrs={"ramp_coef_x_per_m": float(coef[0]),
                             "ramp_coef_y_per_m": float(coef[1]),
                             "ramp_offset": float(coef[2])})


def summarize_rates(rates: dict[str, dict], aoi: xr.DataArray,
                    cls: xr.DataArray | None = None) -> pd.DataFrame:
    """Statistiques (mediane, p10/p90) du taux annuel sur l'AOI, par paire
    et par classe si `cls` est fourni — pour comparer aux ~1-2 cm/an de la
    litterature (Ghezelayagh et al. 2024 : -1.44 cm/an en moyenne)."""
    rows = []
    classes = {"AOI_total": aoi}
    if cls is not None:
        classes.update({f"class_{k}": (cls == k) for k in [1, 2, 3, 4, 5]})
    for pair, r in rates.items():
        rate = r["rate_mm_yr"]
        for name, sel in classes.items():
            v = rate.where(sel)
            finite = v.values[np.isfinite(v.values)]
            if finite.size == 0:
                continue
            rows.append({
                "pair": pair, "dt_years": round(r["dt_years"], 2),
                "zone": name, "n_pixels": finite.size,
                "median_mm_yr": float(np.median(finite)),
                "p10_mm_yr": float(np.percentile(finite, 10)),
                "p90_mm_yr": float(np.percentile(finite, 90)),
            })
    return pd.DataFrame(rows)

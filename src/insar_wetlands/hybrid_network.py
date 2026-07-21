"""Phase A — Test InSAR decisif : reseau HYBRIDE (lignes de base courtes +
annuelles) inverse conjointement, facon Hrysiewicz et al. (2024, RSE 291).

Motivation (cf. revue critique, doc revue_critique_litterature.md) : notre
echec SBAS n'utilisait QUE des paires <= 48 j. Les etudes tourbieres qui
REUSSISSENT en bande C combinent des lignes de base LONGUES (annuelles, qui
portent le signal cumule >> bruit et court-circuitent la decorrelation
estivale) ET courtes (qui densifient le reseau temporel et contraignent la
respiration saisonniere). Ce module construit ce reseau combine et l'inverse
avec le moteur ISBAS deja valide (WLS ponctuel, tolerant aux paires
manquantes), apres correction troposphérique optionnelle.

Ce test tranche la question ouverte : notre echec vient-il du SITE (physique)
ou de la METHODE (design de reseau) ? Si le reseau hybride recupere la
respiration au pixel du laser, c'etait la methode.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import xarray as xr


def build_hybrid_pairs(dates: pd.Series, short_max_days: int = 48,
                       max_pairs_per_date: int = 4,
                       annual_target_month: int = 4, annual_target_day: int = 15,
                       annual_window_days: int = 19,
                       max_gap_years: int = 2) -> pd.DataFrame:
    """Reseau hybride = union {paires courtes <= short_max_days} U {paires
    annuelles meme-mois}. Colonnes : ref_date, sec_date, dt_days, pair, kind
    ('short'|'annual'). Deduplique (une paire ne peut etre que d'un type).
    """
    from .hyp3.jobs import build_sbas_pairs
    from .annual_pairs import gather_candidate_dates, build_annual_pair_list

    short = build_sbas_pairs(dates, short_max_days, max_pairs_per_date).copy()
    short["kind"] = "short"

    inv = pd.DataFrame({"date": pd.to_datetime(dates.values),
                        "granule": [f"d{d:%Y%m%d}" for d in pd.to_datetime(dates.values)]})
    cand = gather_candidate_dates(inv, annual_target_month, annual_target_day,
                                  annual_window_days)
    # une date par annee (la plus proche de la cible) pour ancrer le reseau annuel
    cand = cand.assign(absgap=cand["gap_days"].abs()).sort_values(
        ["year", "absgap"]).groupby("year", as_index=False).first()
    seasonal = cand[["year", "date", "granule"]].copy()
    annual = build_annual_pair_list(seasonal, max_gap_years).copy()
    annual["kind"] = "annual"

    both = pd.concat([short, annual], ignore_index=True)
    both = both.drop_duplicates("pair", keep="first").reset_index(drop=True)
    return both


def invert_hybrid(cropped_root, pairs: list[str], ref_yx: tuple[float, float],
                  lv_theta: xr.DataArray, dry_pairs: xr.DataArray | None = None,
                  tropo_delay_mm: pd.Series | None = None,
                  min_pairs: int = 3, gamma_min: float = 0.20,
                  deramp_mask: xr.DataArray | None = None) -> xr.Dataset:
    """Inversion du reseau hybride -> serie temporelle verticale (mm).

    Etapes : (1) charge unw+corr ; (2) correction troposphérique optionnelle
    (delai differentiel par paire, cf. atmosphere.apply_isbas_tropo) ;
    (3) inversion ISBAS WLS ; (4) projection LOS->vertical ; (5) deramp
    optionnel par pas de temps sur pixels stables.

    Retourne un Dataset {vertical_mm(time,y,x), rms_residual_rad, n_valid_pairs}.
    """
    from .stack import load_layer
    from .geometry import incidence_angle, los_to_vertical
    from .inversion.isbas import invert_stack
    from .annual_pairs import deramp

    unw = load_layer(cropped_root, "unw_phase", pairs)
    corr = load_layer(cropped_root, "corr", pairs)
    if dry_pairs is None:
        dry_pairs = xr.ones_like(unw, dtype=bool)

    if tropo_delay_mm is not None:
        from .atmosphere import apply_isbas_tropo
        inc_mean = float(incidence_angle(lv_theta).mean())
        # aligne la serie de delai sur l'ordre des paires du stack
        dz = tropo_delay_mm.reindex(list(unw.pair.values))
        unw = apply_isbas_tropo(unw, dz, inc_mean)

    isbas = invert_stack(unw, corr, dry_pairs, ref_yx,
                         min_pairs=min_pairs, gamma_min=gamma_min)
    vert = los_to_vertical(isbas.los_displacement_mm, lv_theta)

    if deramp_mask is not None:
        corrected = []
        for ti in range(vert.sizes["time"]):
            layer = vert.isel(time=ti)
            try:
                dr = deramp(layer, deramp_mask)["corrected"]
            except (ValueError, np.linalg.LinAlgError):
                dr = layer
            corrected.append(dr)
        vert = xr.concat(corrected, dim="time").assign_coords(time=vert.time)

    return xr.Dataset({
        "vertical_mm": vert,
        "rms_residual_rad": isbas.rms_residual_rad,
        "n_valid_pairs": isbas.n_valid_pairs,
    })


def network_summary(pairs: pd.DataFrame) -> dict:
    """Statistiques descriptives du reseau hybride (pour la publication)."""
    from .stack import dates_from_pairs

    dates = dates_from_pairs(list(pairs["pair"]))
    n_short = int((pairs["kind"] == "short").sum())
    n_annual = int((pairs["kind"] == "annual").sum())
    return {
        "n_dates": len(dates),
        "n_pairs_total": len(pairs),
        "n_short": n_short,
        "n_annual": n_annual,
        "dt_days_min": int(pairs["dt_days"].min()),
        "dt_days_max": int(pairs["dt_days"].max()),
        "span_days": int((dates.max() - dates.min()).days),
    }

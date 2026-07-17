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

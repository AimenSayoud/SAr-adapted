"""Phase A (sans laser) — validations internes et hydrologiques.

Faute de vérité terrain (laser) immédiate, on peut déjà falsifier ou soutenir
plusieurs de nos conclusions par des tests qui ne demandent AUCUNE donnée
externe :

  1. Décroissance de cohérence vs ligne de base temporelle, par saison
     -> teste directement l'affirmation « décorrélation en quelques jours »
     (Reviewer R1.1). Résultat quantitatif attendu par RSE.
  2. Gain de pixels réseau hybride vs court seul -> teste « échec = méthode ».
  3. Cohérence interne annuelle (fermeture de chaîne : taux 2022->2024 doit
     égaler taux 2022->2023 + 2023->2024) -> teste la fiabilité du déroulement
     sans vérité terrain.
  4. Corrélation du signal InSAR (respiration) avec un moteur hydrologique
     ERA5 (pluie/tcwv) -> si le mouvement de surface suit la physique de la
     nappe, c'est du signal, pas du bruit. Validation SANS laser.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import xarray as xr


def coherence_vs_baseline(corr_stack: xr.DataArray, aoi: xr.DataArray,
                          season_by_month: dict[int, str] | None = None) -> pd.DataFrame:
    """Cohérence moyenne (sur l'AOI) par paire, avec ligne de base temporelle
    et saison de la date de référence. Colonnes : pair, dt_days, month, season,
    mean_coh. Sert à tracer la courbe de décorrélation et à comparer été/hiver.
    """
    season_by_month = season_by_month or {
        12: "hiver", 1: "hiver", 2: "hiver", 3: "printemps", 4: "printemps",
        5: "printemps", 6: "été", 7: "été", 8: "été", 9: "automne",
        10: "automne", 11: "automne"}
    rows = []
    aoi_v = aoi.values
    n = float(aoi_v.sum())
    for p in corr_stack.pair.values:
        ref, sec = str(p).split("_")
        dt = (pd.Timestamp(sec) - pd.Timestamp(ref)).days
        c = corr_stack.sel(pair=p).values
        m = float(np.nansum(np.where(aoi_v, c, np.nan)) / n) if n else np.nan
        month = pd.Timestamp(ref).month
        rows.append({"pair": str(p), "dt_days": dt, "month": month,
                     "season": season_by_month[month], "mean_coh": m})
    return pd.DataFrame(rows)


def decorrelation_summary(coh_df: pd.DataFrame,
                          gamma_usable: float = 0.30) -> pd.DataFrame:
    """Par saison : cohérence médiane, fraction de paires exploitables
    (mean_coh >= gamma_usable), dt médian. Verdict lisible pour l'article."""
    rows = []
    for season, g in coh_df.groupby("season"):
        rows.append({
            "season": season, "n_pairs": len(g),
            "median_coh": float(g["mean_coh"].median()),
            "frac_usable": float((g["mean_coh"] >= gamma_usable).mean()),
            "median_dt_days": float(g["dt_days"].median()),
        })
    return pd.DataFrame(rows).sort_values("median_coh")


def annual_chain_closure(rates_mm_yr: dict[str, float]) -> dict:
    """Cohérence interne des taux annuels : le taux 2 ans doit égaler la
    moyenne des deux taux 1 an consécutifs. `rates_mm_yr` : {"YYYY-YYYY": taux}.
    Retourne les triplets (i->j->k) avec l'écart de fermeture. Un écart faible
    (<< bruit) soutient la fiabilité ; un écart fort révèle des sauts de
    déroulement — SANS vérité terrain.
    """
    def parse(k): a, b = k.split("-"); return int(a), int(b)
    have = {parse(k): v for k, v in rates_mm_yr.items()}
    out = []
    for (i, j), r_ij in have.items():
        if j - i != 1:
            continue
        for (j2, k), r_jk in have.items():
            if j2 != j or k - j != 1:
                continue
            if (i, k) in have:  # doit être un pas de 2 ans
                r_ik = have[(i, k)]
                # taux 2 ans vs moyenne pondérée des deux taux 1 an
                pred = (r_ij + r_jk) / 2.0
                out.append({"chain": f"{i}->{j}->{k}", "r_2yr": r_ik,
                            "r_mean_1yr": pred, "closure_mm_yr": r_ik - pred})
    df = pd.DataFrame(out)
    return {"triplets": df,
            "max_abs_closure": float(df["closure_mm_yr"].abs().max()) if len(df) else np.nan}


def quality_filter(ds: xr.Dataset, aoi: xr.DataArray,
                   rms_max_rad: float = 1.0, min_pairs: int = 10) -> dict:
    """Distingue les pixels VRAIMENT resolus des pixels 'resolus par du bruit'.

    Un pixel avec seulement min_pairs=3 paires satisfait le critere ISBAS mais
    peut n'etre que du bruit. Ce filtre exige un residu WLS faible (rms) ET
    assez de paires. Retourne les comptes brut vs filtre, sur l'AOI et hors.
    """
    solved = np.isfinite(ds["rms_residual_rad"].values)
    good = solved & (ds["rms_residual_rad"].values <= rms_max_rad) \
        & (ds["n_valid_pairs"].values >= min_pairs)
    aoi_v = aoi.values
    return {
        "n_solved_total": int(solved.sum()),
        "n_good_total": int(good.sum()),
        "n_solved_aoi": int((solved & aoi_v).sum()),
        "n_good_aoi": int((good & aoi_v).sum()),
        "frac_good_of_solved": float(good.sum() / max(1, solved.sum())),
        "good_mask": xr.DataArray(good, coords=aoi.coords, dims=aoi.dims),
    }


def clip_to_polygon(da: xr.DataArray, cfg: dict | None = None) -> xr.DataArray:
    """Decoupe une carte au polygone EXACT de la tourbiere (geojson), en
    mettant NaN hors-polygone — pour l'affichage et les stats finales sur la
    tourbiere seule (le traitement, lui, garde le contour pour la reference)."""
    from .stack import aoi_mask
    aoi = aoi_mask(da if da.ndim == 2 else da.isel({d: 0 for d in da.dims if d not in ("y", "x")}), cfg)
    return da.where(aoi)


def coherence_vs_water(coh_mean: xr.DataArray, flooded_frac: xr.DataArray,
                       aoi: xr.DataArray) -> dict:
    """Validation S2 (sans laser) : la cohérence InSAR chute-t-elle la ou S2
    voit de l'eau ? Correle la cohérence moyenne au 'flooded_fraction' (Phase
    5-6) sur l'AOI. Correlation NEGATIVE forte = interpretation physique
    confirmee (l'eau decorrele) par un capteur independant."""
    a = coh_mean.where(aoi).values.ravel()
    b = flooded_frac.reindex_like(coh_mean, method="nearest").where(aoi).values.ravel()
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 10:
        return {"n": int(ok.sum()), "r": np.nan}
    r = float(np.corrcoef(a[ok], b[ok])[0, 1])
    return {"n": int(ok.sum()), "r": r,
            "note": "r negatif attendu : plus d'eau -> moins de coherence"}


def hydrology_proxy(era5: xr.Dataset, lon: float, lat: float,
                    tau_days: int = 30) -> pd.Series:
    """Proxy de niveau de nappe à partir d'ERA5 : anomalie de précipitation
    cumulée sur une fenêtre glissante (tau_days). Une nappe monte avec la
    pluie récente ; c'est un proxy grossier mais physique et INDÉPENDANT de
    l'InSAR — donc utilisable comme quasi-vérité pour tester le signal.
    """
    tp = era5["tp"].sel(latitude=lat, longitude=lon, method="nearest")
    tc = "valid_time" if "valid_time" in tp.coords else "time"
    s = tp.to_series()
    s.index = pd.to_datetime(s.index if tc == "time" else tp[tc].values)
    daily = s.resample("1D").sum()
    roll = daily.rolling(f"{tau_days}D").sum()
    return (roll - roll.mean()).rename("wtd_proxy")


def correlate_insar_hydrology(insar_series: pd.Series, proxy: pd.Series,
                              tolerance_days: int = 6) -> dict:
    """Corrèle la série InSAR (détendance recommandée) au proxy hydrologique
    ERA5 aux dates SAR. r élevé => le mouvement de surface suit la physique de
    la nappe => signal réel (respiration), validé SANS laser."""
    proxy = proxy.copy()
    proxy.index = pd.to_datetime(proxy.index)
    rows = []
    for t, v in insar_series.items():
        dt = np.abs((proxy.index - t).total_seconds().values) / 86400
        j = int(np.argmin(dt))
        if dt[j] <= tolerance_days:
            rows.append((v, proxy.iloc[j]))
    if len(rows) < 3:
        return {"n": len(rows), "r": np.nan}
    a, b = map(np.array, zip(*rows))
    a = a - a.mean(); b = b - b.mean()
    r = float(np.corrcoef(a, b)[0, 1])
    return {"n": len(rows), "r": r, "r2": float(r ** 2)}

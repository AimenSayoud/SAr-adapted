"""Phase A — Validation par le laser in situ (le test decisif).

Le laser mesure l'elevation de surface vraie en UN point, en continu. On
extrait la serie temporelle InSAR au pixel contenant le laser et on la
correle a la serie laser. Deux verdicts possibles :
  - forte correlation (r eleve, faible RMSE) meme si le trend net reste
    incertain  ->  la bande C FONCTIONNE ici (elle recupere la respiration,
    comme les cameras de Hrysiewicz et al. 2024) ; l'echec anterieur venait
    de la METHODE (reseau court seul), pas du site ;
  - correlation faible/nulle  ->  le coeur du fen est genuinement decorrele ;
    l'echec est bien lie au SITE.

Ce module reste agnostique au format exact du laser : `load_laser` accepte un
CSV avec colonnes date + valeur (noms configurables) et une convention de
signe (distance capteur->sol vs elevation).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import xarray as xr


def load_laser(path, date_col: str = "date", value_col: str = "elevation_mm",
               sign: float = 1.0, tz_naive: bool = True) -> pd.Series:
    """Charge la serie laser en mm d'elevation (positif = montee de surface).

    `sign=-1.0` si la valeur brute est une DISTANCE capteur->sol (elle
    diminue quand la surface monte). Reechantillonne en moyenne journaliere.
    """
    df = pd.read_csv(path)
    if date_col not in df or value_col not in df:
        raise KeyError(f"colonnes {date_col}/{value_col} absentes ; "
                       f"disponibles : {list(df.columns)}")
    t = pd.to_datetime(df[date_col])
    if tz_naive:
        t = t.dt.tz_localize(None)
    s = pd.Series((df[value_col].astype(float) * sign).values, index=t,
                  name="laser_mm").sort_index()
    return s.resample("1D").mean().dropna()


def decompose(series: pd.Series, n_harmonics: int = 1) -> dict:
    """Ajuste CONJOINTEMENT tendance lineaire + cycle(s) saisonnier(s).

    L'ajustement joint est indispensable ici : sur seulement ~3 ans, un simple
    detrend lineaire ALIASE le signal saisonnier dans la pente (biais observe :
    -7 mm/an au lieu de -3). On ajuste donc y = a*t + b + Sigma_k
    [c_k sin(2pi k t) + d_k cos(2pi k t)] par moindres carres, ce qui separe
    proprement tendance et respiration.

    Retourne {trend_mm_yr, trend_se_mm_yr, amplitude_mm (crete-a-crete du
    terme saisonnier ajuste), seasonal, trend, detrended}.
    """
    t = series.index
    ty = (t - t[0]).total_seconds().values / (365.25 * 86400)
    cols = [ty, np.ones_like(ty)]
    for k in range(1, n_harmonics + 1):
        cols += [np.sin(2 * np.pi * k * ty), np.cos(2 * np.pi * k * ty)]
    A = np.column_stack(cols)
    coef, *_ = np.linalg.lstsq(A, series.values, rcond=None)
    pred = A @ coef
    resid = series.values - pred
    seasonal = pred - (coef[0] * ty + coef[1])          # part periodique seule
    # erreur-type de la pente
    dof = max(1, len(ty) - A.shape[1])
    sigma2 = float((resid ** 2).sum() / dof)
    cov = sigma2 * np.linalg.inv(A.T @ A)
    return {
        "trend_mm_yr": float(coef[0]),
        "trend_se_mm_yr": float(np.sqrt(cov[0, 0])),
        "amplitude_mm": float(seasonal.max() - seasonal.min()),
        "seasonal": pd.Series(seasonal, index=t, name="seasonal"),
        "trend": pd.Series(coef[0] * ty + coef[1], index=t, name="trend"),
        "detrended": pd.Series(series.values - (coef[0] * ty + coef[1]), index=t,
                               name="detrended"),
    }


def insar_series_at_point(ts: xr.DataArray, y: float, x: float,
                          window_px: int = 1) -> pd.Series:
    """Serie InSAR (mm) au pixel (y,x), moyennee sur une fenetre optionnelle
    (window_px de rayon) pour reduire le bruit ponctuel. `ts` : (time,y,x)."""
    if window_px <= 0:
        pt = ts.sel(y=y, x=x, method="nearest")
    else:
        yi = int(np.argmin(np.abs(ts.y.values - y)))
        xi = int(np.argmin(np.abs(ts.x.values - x)))
        sl = dict(y=slice(max(0, yi - window_px), yi + window_px + 1),
                  x=slice(max(0, xi - window_px), xi + window_px + 1))
        pt = ts.isel(**sl).mean(("y", "x"))
    s = pd.Series(pt.values, index=pd.to_datetime(ts.time.values), name="insar_mm")
    return s.dropna()


def validate_against_laser(insar: pd.Series, laser: pd.Series,
                           tolerance_days: int = 6) -> dict:
    """Correle la serie InSAR aux mesures laser aux dates SAR (plus proche
    voisin, tolerance +/- tolerance_days). Retourne r de Pearson, RMSE,
    biais, n, et les series appariees + detendance des deux signaux.

    On compare les ANOMALIES detendance (respiration) ET les valeurs brutes :
    le test cle est la respiration, la ou l'InSAR est le plus credible.
    """
    laser_d = laser.copy()
    laser_d.index = pd.to_datetime(laser_d.index)
    rows = []
    for t, v in insar.items():
        dt = np.abs((laser_d.index - t).total_seconds().values) / 86400
        j = int(np.argmin(dt))
        if dt[j] <= tolerance_days:
            rows.append((t, v, laser_d.iloc[j]))
    if len(rows) < 3:
        return {"n": len(rows), "r": np.nan, "rmse_mm": np.nan,
                "bias_mm": np.nan, "note": "trop peu de coincidences"}
    tt, iv, lv = zip(*rows)
    iv = np.array(iv); lv = np.array(lv)
    # anomalies (retrait de la moyenne : la reference InSAR est arbitraire)
    ia = iv - iv.mean(); la = lv - lv.mean()
    r = float(np.corrcoef(ia, la)[0, 1]) if len(iv) > 2 else np.nan
    rmse = float(np.sqrt(np.mean((ia - la) ** 2)))
    return {
        "n": len(rows), "r": r, "r2": float(r ** 2) if np.isfinite(r) else np.nan,
        "rmse_mm": rmse, "bias_mm": float((iv - lv).mean()),
        "paired": pd.DataFrame({"date": tt, "insar_mm": iv, "laser_mm": lv}),
    }

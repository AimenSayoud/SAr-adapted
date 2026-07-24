"""Phase I — Le signal diélectrique est-il un CAPTEUR hydrologique ?

La Phase G a établi que le signal saisonnier agrégé (3.3 mm, p=0.011) est
**diélectrique** : le lac oscille identiquement (2.63 mm, même phase) et A−B
s'annule (0.90 mm, p=0.45). Question naturelle : **ce signal suit-il
l'hydrologie ?** Si oui, on passe de « S1 ne mesure pas le mouvement » a
**« S1 mesure l'état hydrique de surface de la tourbière »** — une capacité de
télédétection positive, et la réponse complète a la question « que mesure
réellement Sentinel-1 ici ? ».

**Le piège statistique central, et comment on l'évite.** Série InSAR et forçage
hydrologique sont tous deux fortement **autocorrélés** dans le temps. Une
p-value de corrélation naïve est alors massivement trop optimiste (le nombre
d'observations effectivement indépendantes est bien inférieur a n). On n'utilise
donc PAS de p-value paramétrique : on refait la corrélation sur les **séries
nulles appariées en taille** de la Phase G, qui possèdent **la même structure
temporelle** — la distribution nulle absorbe ainsi l'autocorrélation par
construction. C'est le même protocole que G, appliqué a une autre statistique.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import xarray as xr

from .aggregate import (aggregate_unwrapped, filter_pairs, invert_aggregate,
                        matched_null_zones)


def _detrend(y: np.ndarray, t: np.ndarray) -> np.ndarray:
    """Retire constante + tendance (une dérive commune créerait une corrélation
    fallacieuse avec tout forçage lui-même tendanciel)."""
    M = np.column_stack([np.ones_like(t), t])
    b, *_ = np.linalg.lstsq(M, y, rcond=None)
    return y - M @ b


def build_drivers(era5: xr.Dataset, lon: float, lat: float,
                  s2: xr.Dataset | None = None,
                  zone_mask: xr.DataArray | None = None,
                  api_k: float = 0.9) -> pd.DataFrame:
    """Table des forçages hydrologiques journaliers.

    - `api_mm` : indice de précipitations antécédentes (mémoire hydrologique du
      sol) — le proxy de nappe le plus pertinent en l'absence de WTD in situ.
    - `precip_mm`, `t2m_c` : ERA5 brut.
    - `s2_wetness` : NDWI moyen sur la zone (si S2 et masque fournis) — mesure
      **optique directe** de l'humidité de surface, indépendante du radar.
      C'est le forçage le plus probant : deux capteurs indépendants.
    """
    from .hydro import antecedent_precipitation_index, daily_era5_point

    df = daily_era5_point(era5, lon, lat)
    if "precip_mm" in df:
        df["api_mm"] = antecedent_precipitation_index(df["precip_mm"], k=api_k)
    if s2 is not None and zone_mask is not None:
        var = "ndwi" if "ndwi" in s2 else list(s2.data_vars)[0]
        w = s2[var].where(zone_mask).mean(("y", "x")).to_series()
        w.index = pd.to_datetime(w.index).normalize()
        df = df.join(w.rename("s2_wetness").groupby(level=0).mean(), how="outer")
    return df.sort_index()


def lag_scan(series: pd.DataFrame, drivers: pd.DataFrame,
             date_col: str = "date", value_col: str = "disp_mm",
             max_lag_days: int = 90, step: int = 6,
             detrend: bool = True) -> pd.DataFrame:
    """Corrélation série agrégée ↔ chaque forçage, balayée en DÉCALAGE.

    Le décalage est physiquement informatif : une réponse diélectrique a
    l'humidité est quasi **instantanée** (lag ~0), alors qu'un tassement
    mécanique de la tourbe suivrait la nappe avec un **retard** de plusieurs
    semaines. Retourne, par forçage, le meilleur |r| et le décalage associé.
    """
    s = series[[date_col, value_col]].dropna().copy()
    s[date_col] = pd.to_datetime(s[date_col])
    t = (s[date_col] - s[date_col].iloc[0]).dt.days.values / 365.25
    y = s[value_col].values.astype(float)
    if detrend:
        y = _detrend(y, t)
    rows = []
    for name in drivers.columns:
        d = drivers[name].dropna()
        if d.size < 10:
            continue
        d = d.reindex(pd.date_range(d.index.min(), d.index.max(), freq="1D")
                      ).interpolate()
        best = {"driver": name, "r": np.nan, "lag_days": np.nan, "n": 0}
        for lag in range(0, max_lag_days + 1, step):
            al = d.shift(lag).reindex(s[date_col], method="nearest",
                                      tolerance=pd.Timedelta("3D")).values
            ok = np.isfinite(al) & np.isfinite(y)
            if ok.sum() < 10:
                continue
            x = _detrend(al[ok], t[ok]) if detrend else al[ok]
            if x.std() == 0:
                continue
            r = float(np.corrcoef(x, y[ok])[0, 1])
            if not np.isfinite(best["r"]) or abs(r) > abs(best["r"]):
                best = {"driver": name, "r": round(r, 4), "lag_days": lag,
                        "n": int(ok.sum())}
        rows.append(best)
    return pd.DataFrame(rows).sort_values("r", key=lambda c: c.abs(),
                                          ascending=False).reset_index(drop=True)


def null_lag_scan(unw: xr.DataArray, corr: xr.DataArray, zones: dict,
                  template: xr.DataArray, drivers: pd.DataFrame,
                  n_target: int, n_reference: int, n_trials: int = 50,
                  max_lag_days: int = 90, step: int = 6,
                  max_dt_days: int | None = None) -> pd.DataFrame:
    """Distribution NULLE du |r| maximal, forçage par forçage.

    Chaque tirage rejoue toute la chaîne (deux taches de sol stable appariées en
    taille → agrégation → inversion → balayage de décalage) sur **les mêmes
    forçages**. Comme le balayage retient le MEILLEUR |r| sur ~16 décalages, il
    y a un biais de sélection : le nul doit donc subir **exactement le même
    balayage** pour que la comparaison soit honnête."""
    rows = []
    for t in range(n_trials):
        zn = matched_null_zones(zones, template, n_target, n_reference, seed=t)
        if zn is None:
            break
        try:
            dd = aggregate_unwrapped(unw, corr, zn, "A", "C")
            if max_dt_days is not None:
                dd = filter_pairs(dd, max_dt_days=max_dt_days)
            if len(dd) < 10:
                continue
            sc = lag_scan(invert_aggregate(dd), drivers,
                          max_lag_days=max_lag_days, step=step)
            for _, r in sc.iterrows():
                rows.append({"trial": t, "driver": r["driver"], "r": r["r"]})
        except Exception:
            continue
    return pd.DataFrame(rows, columns=["trial", "driver", "r"])


def driver_pvalues(observed: pd.DataFrame, nulls: pd.DataFrame) -> pd.DataFrame:
    """p-value empirique par forçage : P(|r_nul| >= |r_observé|).

    Test sur |r| : on ne présuppose pas le signe de la relation."""
    rows = []
    for _, o in observed.iterrows():
        v = nulls[nulls.driver == o["driver"]]["r"].abs().dropna().values
        if v.size == 0 or not np.isfinite(o["r"]):
            continue
        k = int((v >= abs(o["r"])).sum())
        rows.append({"driver": o["driver"], "r": o["r"], "lag_days": o["lag_days"],
                     "null_median_absr": round(float(np.median(v)), 4),
                     "null_p95_absr": round(float(np.percentile(v, 95)), 4),
                     "p_value": round((1 + k) / (1 + v.size), 4),
                     "n_null": int(v.size)})
    return (pd.DataFrame(rows).sort_values("p_value").reset_index(drop=True))

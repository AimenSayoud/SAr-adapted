"""Phase 14 — Quantification des incertitudes et produits finaux.

Sorties : vitesse annuelle moyenne (mm/an) + erreur standard + RMSE,
amplitude saisonniere de 'respiration', carte respiration vs subsidence
irreversible. GeoTIFF + PNG.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr


def annual_trend(ts_mm: xr.DataArray) -> xr.DataArray:
    """Tendance lineaire (mm/an) estimee CONJOINTEMENT au cycle annuel.

    Pourquoi pas `fit_velocity`. Une regression lineaire seule sur un signal
    periodique ne rend pas ~0 : sa pente depend de la PHASE du cycle dans la
    fenetre d'observation. Sur deux periodes entieres d'un sinus pur de
    tendance nulle, la pente OLS vaut -0.477 par unite d'amplitude et par an
    (elle est nulle seulement en phase cosinus). Pour un cycle de 7 mm cela
    fait -3.3 mm/an de tendance entierement fictive.

    On ajuste donc y = c + d*t + a*cos(2*pi*t) + b*sin(2*pi*t) et on ne retient
    que d : le terme harmonique absorbe le cycle, ce qui laisse la tendance
    non contaminee. C'est le meme modele que `aggregate.seasonal_amplitude`.
    """
    t = ((ts_mm.time - ts_mm.time[0]) / pd.Timedelta("365.25D")).values
    y = ts_mm.values                                        # (t, y, x)
    shp = y.shape[1:]
    flat = y.reshape(len(t), -1)
    valid = np.isfinite(flat).all(axis=0)
    M = np.column_stack([np.ones_like(t), t,
                         np.cos(2 * np.pi * t), np.sin(2 * np.pi * t)])
    beta, *_ = np.linalg.lstsq(M, np.nan_to_num(flat), rcond=None)
    d = beta[1].astype(float)
    d[~valid] = np.nan
    return xr.DataArray(d.reshape(shp), dims=("y", "x"),
                        coords={k: ts_mm.coords[k] for k in ("y", "x")
                                if k in ts_mm.coords},
                        name="trend_mm_yr")


def seasonal_amplitude(ts_mm: xr.DataArray) -> xr.DataArray:
    """Amplitude mediane pic-a-pic annuelle apres retrait de la tendance.

    L'amplitude reste EMPIRIQUE (max - min par annee civile, puis mediane), et
    non l'amplitude ajustee sqrt(a^2+b^2) : la respiration d'une tourbiere
    n'est pas exactement sinusoidale, et le max-min ne suppose rien de sa
    forme. Seule la tendance retiree est desormais ajustee conjointement au
    cycle -- sinon la pente fictive decrite dans `annual_trend` inclinait la
    serie et rabotait l'amplitude d'environ 12 %.
    """
    t_years = ((ts_mm.time - ts_mm.time[0]) / pd.Timedelta("365.25D")).values
    trend = annual_trend(ts_mm) * xr.DataArray(t_years, dims="time",
                                               coords={"time": ts_mm.time})
    detr = ts_mm - trend
    years = detr.time.dt.year
    amps = []
    for y in np.unique(years):
        sel = detr.sel(time=detr.time.dt.year == y)
        if sel.sizes["time"] >= 8:
            amps.append(sel.max("time") - sel.min("time"))
    amp = xr.concat(amps, dim="year").median("year")
    return amp.rename("seasonal_amplitude_mm")


def breathing_classification(vel: xr.Dataset, amp: xr.DataArray,
                             subsidence_thr_mm_yr: float = -2.0,
                             breathing_thr_mm: float = 10.0) -> xr.DataArray:
    """1=stable, 2=respiration active, 3=subsidence irreversible, 4=les deux.

    Subsidence 'irreversible' si la tendance lineaire depasse le seuil ET
    est significative (|v| > 2 x erreur standard).
    """
    v = vel.velocity_mm_yr
    sig_subs = (v < subsidence_thr_mm_yr) & (np.abs(v) > 2 * vel.velocity_se_mm_yr)
    breathing = amp > breathing_thr_mm
    out = xr.full_like(v, np.nan)
    out = out.where(~(v.notnull() & ~sig_subs & ~breathing), 1)
    out = out.where(~(breathing & ~sig_subs), 2)
    out = out.where(~(sig_subs & ~breathing), 3)
    out = out.where(~(sig_subs & breathing), 4)
    return out.rename("breathing_class")


def export_geotiff(da: xr.DataArray, template: xr.DataArray,
                   path: str | Path) -> Path:
    """Ecrit un produit final en GeoTIFF georeference (CRS de la grille HyP3)."""
    import rioxarray  # noqa: F401

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    out = da.rio.write_crs(template.rio.crs)
    out = out.rio.write_transform(template.rio.transform())
    out.rio.to_raster(path)
    return path


def summary_table(vel: xr.Dataset, amp: xr.DataArray, cls: xr.DataArray,
                  aoi: xr.DataArray) -> pd.DataFrame:
    """Statistiques finales par classe comportementale (pour la these)."""
    rows = []
    for code, name in [(1, "A"), (2, "B"), (3, "C_peat_core"),
                       (4, "D_transition"), (5, "E_water")]:
        sel = cls == code
        v = vel.velocity_mm_yr.where(sel)
        rows.append({
            "class": name,
            "n_solved": int(v.notnull().sum()),
            "vel_mean_mm_yr": float(v.mean()),
            "vel_p10": float(v.quantile(0.10)),
            "vel_p90": float(v.quantile(0.90)),
            "se_mean_mm_yr": float(vel.velocity_se_mm_yr.where(sel).mean()),
            "rmse_mean_mm": float(vel.rmse_mm.where(sel).mean()),
            "amp_median_mm": float(amp.where(sel).median()),
        })
    return pd.DataFrame(rows)

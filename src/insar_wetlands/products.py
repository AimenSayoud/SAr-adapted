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

from .compare import fit_velocity


def seasonal_amplitude(ts_mm: xr.DataArray) -> xr.DataArray:
    """Amplitude mediane pic-a-pic annuelle apres retrait de la tendance."""
    vel = fit_velocity(ts_mm)
    t_years = ((ts_mm.time - ts_mm.time[0]) / pd.Timedelta("365.25D")).values
    trend = vel.velocity_mm_yr * xr.DataArray(t_years, dims="time",
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

"""Phase 10 — Comparaison quantitative SBAS (MintPy) vs ISBAS (custom)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import xarray as xr


def fit_velocity(ts_mm: xr.DataArray) -> xr.Dataset:
    """Vitesse lineaire (mm/an) + erreur standard, pixel a pixel."""
    t_years = ((ts_mm.time - ts_mm.time[0]) / pd.Timedelta("365.25D")).values
    y = ts_mm.values  # (t, y, x)
    valid = np.isfinite(y).all(axis=0)
    t_mean = t_years.mean()
    t_c = t_years - t_mean
    denom = (t_c ** 2).sum()
    y0 = np.nan_to_num(y)
    vel = np.tensordot(t_c, y0 - y0.mean(axis=0), axes=(0, 0)) / denom
    pred = y0.mean(axis=0)[None] + vel[None] * t_c[:, None, None]
    resid = np.where(np.isfinite(y), y - pred, np.nan)
    rmse = np.sqrt(np.nanmean(resid ** 2, axis=0))
    se = rmse / np.sqrt(denom)
    vel[~valid] = np.nan
    rmse[~valid] = np.nan
    se[~valid] = np.nan
    return xr.Dataset({
        "velocity_mm_yr": (("y", "x"), vel),
        "velocity_se_mm_yr": (("y", "x"), se),
        "rmse_mm": (("y", "x"), rmse),
    }, coords={"y": ts_mm.y, "x": ts_mm.x})


def density_report(sbas_vel: xr.Dataset, isbas_vel: xr.Dataset,
                   aoi: xr.DataArray, cls: xr.DataArray) -> pd.DataFrame:
    """Gain en densite de points par classe comportementale."""
    rows = []
    for code, name in [(1, "A"), (2, "B"), (3, "C"), (4, "D"), (5, "E")]:
        in_cls = cls == code
        n_tot = int(in_cls.sum())
        n_sbas = int((in_cls & sbas_vel.velocity_mm_yr.notnull()).sum())
        n_isbas = int((in_cls & isbas_vel.velocity_mm_yr.notnull()).sum())
        rows.append({"class": name, "n_pixels": n_tot,
                     "sbas_solved": n_sbas, "isbas_solved": n_isbas,
                     "sbas_pct": 100 * n_sbas / max(n_tot, 1),
                     "isbas_pct": 100 * n_isbas / max(n_tot, 1),
                     "gain_pixels": n_isbas - n_sbas})
    return pd.DataFrame(rows)


def agreement(sbas_vel: xr.Dataset, isbas_vel: xr.Dataset) -> dict:
    """Accord sur les pixels communs (biais, RMSE, correlation)."""
    a = sbas_vel.velocity_mm_yr
    b = isbas_vel.velocity_mm_yr.interp_like(a, method="nearest")
    both = a.notnull() & b.notnull()
    da, db = a.where(both), b.where(both)
    diff = (db - da)
    return {
        "n_common": int(both.sum()),
        "bias_mm_yr": float(diff.mean()),
        "rmse_mm_yr": float(np.sqrt((diff ** 2).mean())),
        "pearson_r": float(xr.corr(da, db)),
    }


def spatial_continuity(vel: xr.DataArray) -> float:
    """Rugosite du champ de vitesse : ecart median aux 4 voisins (mm/an).

    Des valeurs elevees sur les pixels 'sauves' par l'ISBAS = sauts de phase
    artificiels probables (point de vigilance Phase 10).
    """
    v = vel.values
    diffs = []
    for shift in [(0, 1), (1, 0)]:
        d = np.abs(v - np.roll(v, shift, axis=(0, 1)))
        diffs.append(d)
    all_d = np.concatenate([d[np.isfinite(d)] for d in diffs])
    return float(np.median(all_d)) if all_d.size else np.nan

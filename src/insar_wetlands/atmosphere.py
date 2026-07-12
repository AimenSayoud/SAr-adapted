"""Phase 11 — Correction atmospherique comparative.

Pipeline 1 : calibration locale relative uniquement (le point de reference
Classe A a < 1 km absorbe l'atmosphere commune — hypothese : a l'echelle de
60-90 ha, le delai troposphérique est spatialement quasi uniforme).

Pipeline 2 : correction troposphérique modele AVANT calibration :
  - SBAS/MintPy : tropo 'pyaps' (ERA5) integre a smallbaselineApp ;
  - ISBAS custom : delai zenithal humide estime du TCWV ERA5
    (delai [mm] ~ 6.2 x PWV [mm]), projete en LOS, differencie par paire.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import xarray as xr

PWV_TO_DELAY = 6.2  # facteur de conversion PWV -> delai humide zenithal (~6.1-6.5)


def zenith_wet_delay_mm(era5: xr.Dataset, lon: float, lat: float) -> xr.DataArray:
    """Serie temporelle du delai humide zenithal (mm) au point du site."""
    tcwv = era5["tcwv"].sel(latitude=lat, longitude=lon, method="nearest")
    time_coord = "valid_time" if "valid_time" in tcwv.coords else "time"
    zwd = (tcwv * PWV_TO_DELAY).rename("zwd_mm")
    return zwd.rename({time_coord: "time"}) if time_coord != "time" else zwd


def pair_delays_mm(zwd: xr.DataArray, pairs: list[str],
                   acq_time_utc: str = "05:00") -> pd.Series:
    """Delai differentiel LOS-zenithal par paire (interpole a l'heure S1)."""
    zt = zwd.to_series()
    zt.index = pd.to_datetime(zt.index)

    def at(date_str):
        t = pd.Timestamp(f"{date_str[:4]}-{date_str[4:6]}-{date_str[6:]} {acq_time_utc}")
        i = zt.index.get_indexer([t], method="nearest")[0]
        return zt.iloc[i]

    return pd.Series({p: at(p.split("_")[1]) - at(p.split("_")[0])
                      for p in pairs}, name="dzwd_mm")


def apply_isbas_tropo(unw: xr.DataArray, dzwd_mm: pd.Series,
                      inc_angle_rad: float) -> xr.DataArray:
    """Retire le delai modele (constant spatialement) de chaque interferogramme.

    NB : comme le delai applique est uniforme sur la scene, la calibration par
    point de reference l'annulerait de toute facon — ce pipeline ne differe du
    pipeline 1 que par le bruit du modele. C'est exactement le test voulu.
    """
    from .inversion.isbas import PHASE_TO_MM

    corr_rad = xr.DataArray(
        [dzwd_mm[p] for p in unw.pair.values], dims="pair",
        coords={"pair": unw.pair},
    ) / np.cos(inc_angle_rad) / PHASE_TO_MM
    return unw - corr_rad


def residual_noise_report(ts_by_pipeline: dict[str, xr.DataArray],
                          aoi: xr.DataArray) -> pd.DataFrame:
    """Ecart-type temporel du signal detendance, par pipeline (test aveugle).

    Le 'meilleur' pipeline est celui qui minimise le bruit residuel sur les
    pixels stables — si le modele global degrade, on le documente (Phase 11).
    """
    from .compare import fit_velocity

    rows = []
    for name, ts in ts_by_pipeline.items():
        vel = fit_velocity(ts)
        t_years = ((ts.time - ts.time[0]) / pd.Timedelta("365.25D")).values
        trend = vel.velocity_mm_yr * xr.DataArray(t_years, dims="time",
                                                  coords={"time": ts.time})
        resid = ts - trend
        rows.append({
            "pipeline": name,
            "std_all_mm": float(resid.std("time").mean()),
            "std_aoi_mm": float(resid.std("time").where(aoi).mean()),
            "rmse_mm": float(vel.rmse_mm.mean()),
        })
    return pd.DataFrame(rows)

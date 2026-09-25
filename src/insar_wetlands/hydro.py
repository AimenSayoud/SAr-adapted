"""Phase 13 — Analyse hydrologique et separation du signal dielectrique.

Confrontation de la serie de deplacement vertical avec :
  - precipitations accumulees (API, antecedent precipitation index) ;
  - cycles gel/degel (T2m ERA5) ;
  - NDWI moyen du site (etat hydrique de surface).

Drapeau 'dielectrique suspect' : affaissement InSAR synchrone d'une chute
drastique du NDWI -> le signal peut etre une penetration plus profonde du
radar dans la tourbe seche, pas un tassement physique.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr


def era5_samples_per_day(era5: xr.Dataset) -> float:
    """Median number of time steps per day in an ERA5 dataset (24 = hourly)."""
    tname = "valid_time" if "valid_time" in era5.coords else "time"
    t = pd.to_datetime(era5[tname].values)
    if len(t) < 2:
        return float("nan")
    step_h = float(np.median(np.diff(t).astype("timedelta64[m]").astype(float))) / 60.0
    return 24.0 / step_h


def open_era5(drive: str | Path, hourly_dir: str | Path | None = None) -> xr.Dataset:
    """ERA5 at the site for anything that sums precipitation (X-043).

    ``tp`` is an hourly accumulation, so a daily total needs all 24 hours. The original
    download (``era5_rzecin.nc``) kept 00/06/12/18 UTC only: daily sums covered 4 of 24 h.
    Prefer the hourly re-download in ``<drive>/era5_hourly/`` (``era5_<year>_tp_hourly.nc`` +
    ``era5_<year>_inst_hourly.nc``); fall back to ``era5_rzecin.nc`` with a warning.

    Phases that only use instantaneous fields (t2m for freeze flags, tcwv for delays and
    pair scoring) still open ``era5_rzecin.nc`` directly, so their results are unchanged.
    """
    drive = Path(drive)
    h = Path(hourly_dir) if hourly_dir else drive / "era5_hourly"
    tp_files = sorted(h.glob("era5_*_tp_hourly.nc")) if h.is_dir() else []
    if tp_files:
        parts = []
        for tp_f in tp_files:
            inst_f = tp_f.with_name(tp_f.name.replace("_tp_hourly", "_inst_hourly"))
            tp = xr.open_dataset(tp_f)[["tp"]]
            inst = xr.open_dataset(inst_f)[["t2m", "tcwv"]] if inst_f.exists() else None
            parts.append(xr.merge([tp, inst], compat="override", join="outer") if inst is not None else tp)
        tname = "valid_time" if "valid_time" in parts[0].coords else "time"
        ds = xr.concat(parts, dim=tname).sortby(tname)
        return ds.drop_vars([v for v in ("number", "expver") if v in ds.coords or v in ds.data_vars])
    legacy = drive / "era5_rzecin.nc"
    ds = xr.open_dataset(legacy)
    n = era5_samples_per_day(ds)
    if "tp" in ds and n < 23.5:
        warnings.warn(f"{legacy.name}: {n:.0f} samples/day — daily precipitation covers {n:.0f}/24 h "
                      f"(X-043). Put the hourly download in {h}.", stacklevel=2)
    return ds


def daily_era5_point(era5: xr.Dataset, lon: float, lat: float) -> pd.DataFrame:
    """Precip journaliere (mm) et T2m moyenne (degC) au point du site."""
    pt = era5.sel(latitude=lat, longitude=lon, method="nearest")
    tname = "valid_time" if "valid_time" in pt.coords else "time"
    df = pd.DataFrame(index=pd.to_datetime(pt[tname].values))
    if "tp" in pt:
        df["precip_mm"] = pt["tp"].values * 1000.0
    if "t2m" in pt:
        df["t2m_c"] = pt["t2m"].values - 273.15
    daily = df.resample("1D").agg({"precip_mm": "sum", "t2m_c": "mean"})
    return daily


def antecedent_precipitation_index(precip_mm: pd.Series,
                                   k: float = 0.9) -> pd.Series:
    """API_t = k * API_{t-1} + P_t (memoire hydrologique du sol)."""
    api = np.zeros(len(precip_mm))
    p = precip_mm.fillna(0).values
    for i in range(1, len(p)):
        api[i] = k * api[i - 1] + p[i]
    return pd.Series(api, index=precip_mm.index, name="api_mm")


def freeze_flags(t2m_c: pd.Series, window_days: int = 5) -> pd.Series:
    """True si la moyenne glissante T2m < 0 degC (coherence hivernale douteuse)."""
    return (t2m_c.rolling(f"{window_days}D").mean() < 0).rename("frozen")


def site_series(d_vert: xr.DataArray, aoi: xr.DataArray,
                cls: xr.DataArray | None = None,
                class_code: int = 3) -> pd.Series:
    """Serie moyenne du deplacement vertical sur le coeur de la tourbiere."""
    sel = aoi if cls is None else (aoi & (cls == class_code))
    ts = d_vert.where(sel).mean(("y", "x")).to_series()
    return ts.rename("d_vert_mm")


def lag_correlation(insar: pd.Series, driver: pd.Series,
                    max_lag_days: int = 60, step: int = 6) -> pd.DataFrame:
    """Correlation croisee InSAR vs forcage hydrologique, par decalage."""
    drv = driver.reindex(
        pd.date_range(driver.index.min(), driver.index.max(), freq="1D")
    ).interpolate()
    rows = []
    for lag in range(0, max_lag_days + 1, step):
        shifted = drv.shift(lag)
        aligned = shifted.reindex(insar.index, method="nearest",
                                  tolerance=pd.Timedelta("3D"))
        ok = insar.notna() & aligned.notna()
        if ok.sum() > 5:
            r = float(np.corrcoef(insar[ok], aligned[ok])[0, 1])
            rows.append({"lag_days": lag, "pearson_r": r, "n": int(ok.sum())})
    return pd.DataFrame(rows)


def dielectric_suspect_epochs(insar: pd.Series, ndwi_site: pd.Series,
                              subsidence_mm: float = -5.0,
                              ndwi_drop: float = -0.10,
                              window: str = "24D") -> pd.DataFrame:
    """Epoques ou affaissement et chute de NDWI sont synchrones (ambiguite).

    Ces dates doivent etre discutees (Phase 13) : tassement poroelastique reel
    vs artefact dielectrique (penetration accrue dans la tourbe seche).
    """
    d_insar = insar.diff()
    ndwi = ndwi_site.reindex(insar.index, method="nearest",
                             tolerance=pd.Timedelta(window))
    d_ndwi = ndwi.diff()
    suspect = (d_insar <= subsidence_mm) & (d_ndwi <= ndwi_drop)
    return pd.DataFrame({"d_insar_mm": d_insar, "d_ndwi": d_ndwi,
                         "dielectric_suspect": suspect}).dropna()

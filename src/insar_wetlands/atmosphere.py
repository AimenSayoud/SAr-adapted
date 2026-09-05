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
                   acq_time_utc: str = "17:00") -> pd.Series:
    """Delai differentiel LOS-zenithal par paire (interpole a l'heure S1).

    NB : la trace S1 de Rzecin (orbite 175 ascendante) passe a ~16:36 UTC
    (cf. noms de granules T1636xx) ; l'heure ERA5 la plus proche est 17:00.
    L'ancien defaut 05:00 echantillonnait l'atmosphere du MATIN (biais).
    Pour une requete GACOS, utiliser aussi 16:36 UTC comme heure d'acquisition.
    """
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


def load_gacos_ztd(gacos_dir, date: str, template: xr.DataArray) -> xr.DataArray:
    """Charge un delai troposphérique GACOS (.ztd + .ztd.rsc) pour une date et
    le reechantillonne sur la grille du crop.

    GACOS (http://www.gacos.net) fournit un delai zenithal total SPATIALEMENT
    VARIABLE (~0.125 deg, composante stratifiee + turbulente), contrairement
    au TCWV ERA5 scalaire qui, uniforme sur une petite AOI, est annule par la
    calibration au point de reference. GACOS est donc la seule correction
    susceptible d'ameliorer un petit site plat — a condition que l'utilisateur
    ait telecharge les fichiers (requete manuelle sur gacos.net).

    Le .ztd est un binaire float32 (row-major), geo-reference par le .rsc
    (WIDTH, FILE_LENGTH, X_FIRST, Y_FIRST, X_STEP, Y_STEP). Delai en metres.
    """
    from pathlib import Path

    import rioxarray  # noqa: F401

    gacos_dir = Path(gacos_dir)
    # GACOS GeoTIFF (format demande ici) : deja geo-reference -> reproject_match
    # aligne CRS + grille sur le crop en une passe (robuste).
    tif = next(gacos_dir.glob(f"*{date}*.tif"), None)
    if tif is not None:
        da = rioxarray.open_rasterio(tif, masked=True).squeeze("band", drop=True)
        if da.rio.crs is None:
            da = da.rio.write_crs("EPSG:4326")
        return da.rio.reproject_match(template)
    # Fallback binaire .ztd + .rsc (format historique GACOS)
    ztd = next(gacos_dir.glob(f"*{date}*.ztd"), None)
    if ztd is None:
        raise FileNotFoundError(f"GACOS absent pour {date} dans {gacos_dir} "
                                "(.tif ou .ztd)")
    rsc = Path(str(ztd) + ".rsc")
    meta = {}
    for line in rsc.read_text().splitlines():
        parts = line.split()
        if len(parts) >= 2:
            meta[parts[0]] = parts[1]
    w, h = int(meta["WIDTH"]), int(meta["FILE_LENGTH"])
    x0, y0 = float(meta["X_FIRST"]), float(meta["Y_FIRST"])
    dx, dy = float(meta["X_STEP"]), float(meta["Y_STEP"])
    arr = np.fromfile(ztd, dtype="float32").reshape(h, w)
    lons = x0 + dx * np.arange(w)
    lats = y0 + dy * np.arange(h)
    da = xr.DataArray(arr, coords={"y": lats, "x": lons}, dims=("y", "x"),
                      name="ztd_m").rio.write_crs("EPSG:4326")
    return da.rio.reproject_match(template)


def apply_gacos_tropo(unw: xr.DataArray, gacos_dir, template: xr.DataArray,
                      inc: xr.DataArray) -> xr.DataArray:
    """Retire le delai GACOS differentiel (SPATIALEMENT VARIABLE) par paire.

    Pour chaque paire ref_sec : delta_ztd = ztd(sec) - ztd(ref) [m] -> mm ->
    LOS (/cos inc) -> radians. Contrairement a apply_isbas_tropo (scalaire),
    ceci n'est PAS annule par la reference et corrige la structure spatiale.
    """
    from .inversion.isbas import PHASE_TO_MM

    cache: dict[str, xr.DataArray] = {}

    def ztd(date):
        if date not in cache:
            cache[date] = load_gacos_ztd(gacos_dir, date, template)
        return cache[date]

    inc_grid = inc.interp(y=unw.y, x=unw.x, method="nearest") \
        if inc.sizes.get("y") != unw.sizes.get("y") else inc
    out = []
    for p in unw.pair.values:
        ref, sec = p.split("_")
        d_ztd_mm = (ztd(sec) - ztd(ref)) * 1000.0
        d_los_rad = (d_ztd_mm / np.cos(inc_grid)) / PHASE_TO_MM
        out.append(unw.sel(pair=p) - d_los_rad)
    return xr.concat(out, dim="pair").assign_coords(pair=unw.pair)


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

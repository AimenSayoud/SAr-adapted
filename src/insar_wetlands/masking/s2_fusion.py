"""Phase 4 — Fusion Sentinel-2 -> grille UTM HyP3 (nearest, jamais bilineaire).

Les produits HyP3 sont deja geocodes UTM : la 'fusion' est un crop + resample
nearest de S2 (10 m) vers la grille radar (~40 m). Le nearest preserve les
discontinuites eau/tourbe du fen transitionnel.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

# Classes SCL Sentinel-2 considerees invalides (nuage/ombre/neige/no-data)
SCL_INVALID = [0, 1, 3, 8, 9, 10, 11]


def load_s2_date(item, bands: list[str], template: xr.DataArray) -> xr.Dataset:
    """Charge les bandes d'un item STAC, croppe et reprojette sur la grille HyP3."""
    import rioxarray

    out = {}
    bounds = template.rio.transform_bounds("EPSG:4326")
    for band in bands:
        href = item.assets[band].href
        da = rioxarray.open_rasterio(href, masked=True).squeeze("band", drop=True)
        da = da.rio.clip_box(*bounds, crs="EPSG:4326")
        out[band] = da.rio.reproject_match(
            template, resampling=_nearest())
    ds = xr.Dataset(out)
    valid = ~ds["scl"].isin(SCL_INVALID) if "scl" in ds else None
    g, n = ds["green"], ds["nir"]
    ds["ndwi"] = (g - n) / (g + n)
    if "swir16" in ds:
        s = ds["swir16"]
        ds["mndwi"] = (g - s) / (g + s)
    if valid is not None:
        for v in ("ndwi", "mndwi"):
            if v in ds:
                ds[v] = ds[v].where(valid)
    return ds


def _nearest():
    from rasterio.enums import Resampling

    return Resampling.nearest


def build_s2_stack(items: list, cfg: dict, template: xr.DataArray,
                   out_nc: str | Path) -> Path:
    """Construit le stack temporel NDWI/MNDWI et l'ecrit sur Drive (netCDF).

    Idempotent : reprend ou le fichier s'est arrete si des dates manquent.
    """
    out_nc = Path(out_nc)
    out_nc.parent.mkdir(parents=True, exist_ok=True)
    existing = None
    if out_nc.exists():
        existing = xr.open_dataset(out_nc)
        done = set(pd.to_datetime(existing.time.values))
    else:
        done = set()
    new = []
    for item in items:
        t = pd.to_datetime(item.datetime).tz_localize(None).normalize()
        if t in done:
            continue
        try:
            ds = load_s2_date(item, cfg["sentinel2"]["bands"], template)
        except Exception as e:  # scene corrompue/href mort : on saute
            print(f"  ! {item.id}: {e}")
            continue
        keep = ds[[v for v in ("ndwi", "mndwi") if v in ds]]
        new.append(keep.expand_dims(time=[t]))
        print(f"  + {t.date()} ({item.id})")
    if new:
        stack = xr.concat(new, dim="time").sortby("time")
        if existing is not None:
            stack = xr.concat([existing, stack], dim="time").sortby("time")
            existing.close()
        tmp = out_nc.with_suffix(".tmp.nc")
        stack.to_netcdf(tmp)
        tmp.replace(out_nc)
    return out_nc

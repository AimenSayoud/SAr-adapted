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
    else:  # garantit un jeu de variables constant pour les concat
        ds["mndwi"] = xr.full_like(ds["ndwi"], np.nan)
    if valid is not None:
        for v in ("ndwi", "mndwi"):
            ds[v] = ds[v].where(valid)
    return ds


def _nearest():
    from rasterio.enums import Resampling

    return Resampling.nearest


def _item_sort_key(item):
    """Trie par (date, nuages) : en cas de doublon de date (2 tuiles MGRS),
    la scene la moins nuageuse est traitee en premier et gagne."""
    t = pd.to_datetime(item.datetime).tz_localize(None).normalize()
    return (t, item.properties.get("eo:cloud_cover", 100.0))


def _flush(existing: xr.Dataset | None, new: list, out_nc: Path) -> xr.Dataset:
    """Ecrit existant+nouveau sur disque (atomique) et retourne le cumul."""
    stack = xr.concat(new, dim="time")
    if existing is not None:
        stack = xr.concat([existing, stack], dim="time")
    stack = stack.sortby("time").load()
    tmp = out_nc.with_suffix(".tmp.nc")
    stack.to_netcdf(tmp)
    tmp.replace(out_nc)
    return stack


def build_s2_stack(items: list, cfg: dict, template: xr.DataArray,
                   out_nc: str | Path, checkpoint_every: int = 20) -> Path:
    """Construit le stack temporel NDWI/MNDWI et l'ecrit sur Drive (netCDF).

    Idempotent et resistant aux coupures Colab : une date deja presente est
    sautee, un checkpoint est ecrit toutes les `checkpoint_every` dates.
    Les doublons de date (2 tuiles MGRS le meme jour) sont dedupliques —
    la scene la moins nuageuse gagne.
    """
    out_nc = Path(out_nc)
    out_nc.parent.mkdir(parents=True, exist_ok=True)
    # .load() : tout en memoire (stack ~20 Mo) pour pouvoir reecrire le
    # fichier sans garder de references lazy dessus.
    existing = xr.load_dataset(out_nc) if out_nc.exists() else None
    done = (set(pd.to_datetime(existing.time.values))
            if existing is not None else set())

    new = []
    for item in sorted(items, key=_item_sort_key):
        t = pd.to_datetime(item.datetime).tz_localize(None).normalize()
        if t in done:
            continue
        try:
            ds = load_s2_date(item, cfg["sentinel2"]["bands"], template)
        except Exception as e:  # scene corrompue/href mort : on saute
            print(f"  ! {item.id}: {e}")
            continue
        new.append(ds[["ndwi", "mndwi"]].expand_dims(time=[t]))
        done.add(t)
        print(f"  + {t.date()} ({item.id})")
        if len(new) >= checkpoint_every:
            existing = _flush(existing, new, out_nc)
            new = []
            print(f"  ... checkpoint ({existing.time.size} dates sur disque)")
    if new:
        _flush(existing, new, out_nc)
    return out_nc

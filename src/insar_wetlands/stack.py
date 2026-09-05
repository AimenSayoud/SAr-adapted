"""Chargement des produits HyP3 croppes en stacks xarray (dim 'pair')."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr


def align_grid(da: xr.DataArray | xr.Dataset,
              template: xr.DataArray) -> xr.DataArray | xr.Dataset:
    """Reassigne les coordonnees y/x de `da` sur celles de `template`.

    Necessaire quand deux produits representant EXACTEMENT le meme crop
    (meme shape) ont des coordonnees legerement differentes en valeur
    (ex: grille MintPy reconstruite depuis X_FIRST/Y_STEP du HDF5, vs
    grille rioxarray lue directement du GeoTIFF) : xarray aligne par
    egalite EXACTE de coordonnees, pas par position -> toute comparaison
    booleenne/where entre les deux donne une intersection vide en
    silence (aucune erreur), symptome typique : des comptages a 0 alors
    que les deux cartes ont clairement des donnees.
    """
    y_dim = "y" if "y" in da.dims else None
    x_dim = "x" if "x" in da.dims else None
    if y_dim and da.sizes["y"] != template.sizes["y"]:
        raise ValueError(f"tailles y incompatibles: {da.sizes['y']} vs "
                         f"{template.sizes['y']} -> pas le meme crop")
    if x_dim and da.sizes["x"] != template.sizes["x"]:
        raise ValueError(f"tailles x incompatibles: {da.sizes['x']} vs "
                         f"{template.sizes['x']} -> pas le meme crop")
    coords = {}
    if y_dim:
        coords["y"] = template.y.values
    if x_dim:
        coords["x"] = template.x.values
    return da.assign_coords(coords)


def list_pairs(cropped_root: str | Path) -> list[str]:
    root = Path(cropped_root)
    return sorted(d.name for d in root.iterdir()
                  if d.is_dir() and (d / "ok").exists())


def _find_layer(pair_dir: Path, layer: str) -> Path | None:
    hits = list(pair_dir.glob(f"*_{layer}.tif"))
    return hits[0] if hits else None


def load_layer(cropped_root: str | Path, layer: str,
               pairs: list[str] | None = None) -> xr.DataArray:
    """Stack (pair, y, x) d'une couche HyP3 ('unw_phase', 'corr', ...).

    Toutes les paires d'un meme burst partagent la meme grille UTM ; on aligne
    par securite sur la premiere (nearest, tolerance demi-pixel).
    """
    import rioxarray  # noqa: F401  (active l'accessor .rio)

    root = Path(cropped_root)
    pairs = pairs or list_pairs(root)
    das, kept = [], []
    template = None
    for pair in pairs:
        tif = _find_layer(root / pair, layer)
        if tif is None:
            continue
        da = rioxarray.open_rasterio(tif, masked=True).squeeze("band", drop=True)
        if template is None:
            template = da
        elif da.shape != template.shape or not np.allclose(
                da.x, template.x, atol=abs(float(template.x[1] - template.x[0])) / 2):
            da = da.rio.reproject_match(template)
        das.append(da)
        kept.append(pair)
    stack = xr.concat(das, dim="pair")
    stack = stack.assign_coords(pair=kept)
    ref = pd.to_datetime([p.split("_")[0] for p in kept])
    sec = pd.to_datetime([p.split("_")[1] for p in kept])
    stack = stack.assign_coords(ref_date=("pair", ref), sec_date=("pair", sec))
    stack.name = layer
    return stack


def load_static_layer(cropped_root: str | Path, layer: str) -> xr.DataArray:
    """Couche statique (dem, lv_theta...) : lue depuis la premiere paire."""
    import rioxarray

    root = Path(cropped_root)
    for pair in list_pairs(root):
        tif = _find_layer(root / pair, layer)
        if tif is not None:
            return rioxarray.open_rasterio(tif, masked=True).squeeze(
                "band", drop=True)
    raise FileNotFoundError(f"couche '{layer}' absente de {root}")


def aoi_mask(template: xr.DataArray, cfg: dict | None = None) -> xr.DataArray:
    """Masque booleen (y, x) : True a l'interieur du polygone Rzecin."""
    from pyproj import Transformer
    from rasterio.features import geometry_mask
    from shapely.ops import transform as shp_transform

    from .aoi import load_aoi

    geom = load_aoi(cfg)
    tr = Transformer.from_crs("EPSG:4326", template.rio.crs, always_xy=True)
    geom_utm = shp_transform(tr.transform, geom)
    mask = ~geometry_mask([geom_utm], out_shape=template.shape,
                          transform=template.rio.transform(), invert=False)
    return xr.DataArray(mask, coords={"y": template.y, "x": template.x})


def dates_from_pairs(pairs: list[str]) -> pd.DatetimeIndex:
    ds = set()
    for p in pairs:
        a, b = p.split("_")
        ds.add(a)
        ds.add(b)
    return pd.DatetimeIndex(sorted(pd.to_datetime(list(ds))))

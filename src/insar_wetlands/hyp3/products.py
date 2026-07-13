"""Phase 2 — Telechargement, extraction et crop des produits HyP3 vers Drive.

Les zips HyP3 (~100 Mo) sont extraits, croppes autour de l'AOI (+ marge pour
le point de reference Classe A), puis supprimes : seuls les GeoTIFF croppes
(quelques Mo) et le .txt de metadonnees sont conserves.
"""

from __future__ import annotations

import shutil
import zipfile
from pathlib import Path

# Couches conservees apres crop (suffixes des GeoTIFF HyP3 burst InSAR).
# NB: le produit ISCE burst nomme les composantes connexes '_conncomp.tif'
# (PAS '_conn_comp.tif' comme les produits GAMMA) — les deux formes sont
# gardees. 'amp' = amplitude de retrodiffusion.
KEEP_LAYERS = ["unw_phase", "wrapped_phase", "corr", "lv_theta", "lv_phi",
               "dem", "water_mask", "conn_comp", "conncomp", "amp"]


def utm_crop_bounds(cfg: dict, template_crs, margin_m: float | None = None):
    """BBox de crop en coordonnees du produit (UTM), AOI + marge config."""
    from pyproj import Transformer

    from ..aoi import load_aoi

    margin = float(margin_m if margin_m is not None
                   else cfg["hyp3"].get("crop_buffer_m", 2000))
    geom = load_aoi(cfg)
    tr = Transformer.from_crs("EPSG:4326", template_crs, always_xy=True)
    minx, miny, maxx, maxy = geom.bounds
    x0, y0 = tr.transform(minx, miny)
    x1, y1 = tr.transform(maxx, maxy)
    return (min(x0, x1) - margin, min(y0, y1) - margin,
            max(x0, x1) + margin, max(y0, y1) + margin)


def crop_product(product_dir: Path, cfg: dict, out_root: Path,
                 pair: str) -> Path:
    """Croppe toutes les couches utiles d'un produit HyP3 extrait."""
    import rioxarray

    out_dir = Path(out_root) / pair
    out_dir.mkdir(parents=True, exist_ok=True)
    tifs = list(Path(product_dir).glob("*.tif"))
    bounds = None
    for tif in tifs:
        layer = next((l for l in KEEP_LAYERS if tif.stem.endswith(l)), None)
        if layer is None:
            continue
        da = rioxarray.open_rasterio(tif, masked=True)
        if bounds is None:
            bounds = utm_crop_bounds(cfg, da.rio.crs)
        da.rio.clip_box(*bounds).rio.to_raster(out_dir / tif.name)
        da.close()
    for txt in Path(product_dir).glob("*.txt"):
        shutil.copy(txt, out_dir / txt.name)
    return out_dir


def _crop_is_complete(pair_dir: Path, required_layers=("conncomp",)) -> bool:
    """Vrai si le crop existe ET contient les couches requises.

    Sert au rattrapage : les crops historiques ont perdu _conncomp.tif
    (mauvais motif de nom) -> ils sont consideres incomplets et la paire est
    re-telechargee (gratuit tant que le job HyP3 n'a pas expire, 14 jours).
    """
    if not (pair_dir / "ok").exists():
        return False
    for layer in required_layers:
        if not list(pair_dir.glob(f"*_{layer}.tif")):
            return False
    return True


def download_and_crop(jobs, cfg: dict, drive_root: str | Path,
                      pair_of_job: dict | None = None,
                      required_layers=("conncomp",)) -> list[str]:
    """Telecharge chaque job reussi, extrait, croppe, nettoie. Idempotent.

    pair_of_job : mapping job_id -> nom de paire (sinon deduit du nom de zip).
    required_layers : couches dont l'absence declenche un re-telechargement
    du produit (rattrapage des crops incomplets).
    Retourne la liste des paires traitees.
    """
    drive_root = Path(drive_root)
    tmp = drive_root / "_tmp_zips"
    cropped_root = drive_root / "hyp3_cropped"
    tmp.mkdir(parents=True, exist_ok=True)
    done = []
    for job in jobs:
        if not job.succeeded():
            continue
        pair = (pair_of_job or {}).get(job.job_id)
        if pair is None:
            # nom de paire deductible des parametres du job sans telecharger
            g = (getattr(job, "job_parameters", None) or {}).get("granules") or []
            if len(g) >= 2:
                pair = f"{g[0].split('_')[3][:8]}_{g[1].split('_')[3][:8]}"
        if pair and _crop_is_complete(cropped_root / pair, required_layers):
            done.append(pair)
            continue
        files = job.download_files(location=str(tmp))
        for z in files:
            z = Path(z)
            extract_dir = tmp / z.stem
            with zipfile.ZipFile(z) as zf:
                zf.extractall(extract_dir)
            product_dir = next(d for d in extract_dir.iterdir() if d.is_dir()) \
                if not list(extract_dir.glob("*.tif")) else extract_dir
            if pair is None:
                pair = _pair_from_product_name(product_dir.name)
            out = crop_product(product_dir, cfg, cropped_root, pair)
            (out / "ok").touch()
            shutil.rmtree(extract_dir, ignore_errors=True)
            z.unlink(missing_ok=True)
            done.append(pair)
    return done


def _pair_from_product_name(name: str) -> str:
    """Extrait 'YYYYMMDD_YYYYMMDD' d'un nom de produit HyP3 burst."""
    import re

    m = re.findall(r"(20\d{6})T?", name)
    if len(m) >= 2:
        return f"{m[0]}_{m[1]}"
    return name

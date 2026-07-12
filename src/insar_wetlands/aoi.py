"""Chargement et manipulation de la zone d'interet (AOI)."""

from __future__ import annotations

import json
import math
from pathlib import Path

from shapely.geometry import shape
from shapely.geometry.base import BaseGeometry

from .config import load_config, repo_root


def load_aoi(cfg: dict | None = None) -> BaseGeometry:
    cfg = cfg or load_config()
    path = repo_root() / cfg["site"]["aoi_geojson"]
    with open(path) as f:
        gj = json.load(f)
    geom = shape(gj["features"][0]["geometry"])
    # Les KML exportent souvent des coordonnees 3D (lon, lat, 0) — on aplatit.
    if geom.has_z:
        from shapely.ops import transform

        geom = transform(lambda x, y, z=None: (x, y), geom)
    # Les traces KML contiennent frequemment des auto-intersections : on repare
    # et on garde le plus grand polygone.
    if not geom.is_valid:
        from shapely.validation import make_valid

        geom = make_valid(geom)
        if geom.geom_type in ("MultiPolygon", "GeometryCollection"):
            polys = [g for g in geom.geoms if g.geom_type == "Polygon"]
            geom = max(polys, key=lambda g: g.area)
    return geom


def aoi_wkt(cfg: dict | None = None, simplify_deg: float = 1e-4) -> str:
    """WKT simplifie pour les requetes API (3380 sommets = URL trop longue)."""
    geom = load_aoi(cfg)
    simple = geom.simplify(simplify_deg, preserve_topology=True)
    if not simple.is_valid or simple.is_empty:
        simple = geom.convex_hull
    return simple.wkt


def buffered_bbox(cfg: dict | None = None) -> tuple[float, float, float, float]:
    """(minx, miny, maxx, maxy) en degres, avec la marge buffer_m du config."""
    cfg = cfg or load_config()
    geom = load_aoi(cfg)
    buffer_m = float(cfg["site"]["buffer_m"])
    lat0 = geom.centroid.y
    dlat = buffer_m / 111_132.0
    dlon = buffer_m / (111_320.0 * math.cos(math.radians(lat0)))
    minx, miny, maxx, maxy = geom.bounds
    return (minx - dlon, miny - dlat, maxx + dlon, maxy + dlat)


def area_ha(cfg: dict | None = None) -> float:
    """Surface du polygone en hectares (projection equirectangulaire locale)."""
    geom = load_aoi(cfg)
    lat0 = geom.centroid.y
    mlat = 111_132.0
    mlon = 111_320.0 * math.cos(math.radians(lat0))
    from shapely.ops import transform

    proj = transform(lambda x, y: (x * mlon, y * mlat), geom)
    return proj.area / 1e4

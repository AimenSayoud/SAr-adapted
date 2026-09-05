"""Contrôle externe indépendant — European Ground Motion Service (EGMS).

EGMS (Copernicus) fournit, pour toute l'Europe, une vitesse de déformation du
sol dérivée de Sentinel-1 par un traitement PSI/DS opérationnel de tiers
(produit Ortho : grille 100 m, EPSG:3035, composantes verticale + E-O, en
mm/an). C'est un arbitre INDÉPENDANT de notre pipeline.

Réserve de rigueur importante : EGMS repose principalement sur du PSI
(diffuseurs PERSISTANTS = bâti, roche, routes). Sur une tourbière nue et
végétalisée, il y a peu ou pas de diffuseurs persistants -> un résultat VIDE
sur le tapis est ATTENDU et ne prouve PAS à lui seul « site difficile » (il
confirmerait seulement que le PSI ne convient pas aux cibles distribuées
végétalisées, ce qui est connu). Ce que le contrôle apporte réellement :
  - le CONTEXTE régional (y a-t-il de la subsidence détectée autour du site ?);
  - un test faible mais gratuit : si EGMS montre malgré tout des points sur/en
    bordure de tourbière avec un signal cohérent, c'est informatif ;
  - une comparaison des marges stables (routes, digues) avec notre référence.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import xarray as xr


def load_egms_ortho(tif_path, aoi_geom_wgs84=None) -> xr.DataArray:
    """Charge une tuile EGMS Ortho (vitesse verticale, mm/an) en DataArray.

    Le GeoTIFF EGMS est en EPSG:3035. On le laisse dans son CRS natif ; le
    filtrage AOI se fait en reprojetant l'AOI vers 3035 (moins de rééchant.).
    """
    import rioxarray  # noqa: F401

    da = rioxarray.open_rasterio(tif_path, masked=True).squeeze("band", drop=True)
    da.name = "egms_vel_mm_yr"
    return da


def points_in_aoi(egms: xr.DataArray, cfg: dict | None = None,
                  buffer_m: float = 500.0) -> dict:
    """Extrait les valeurs EGMS valides dans l'AOI et dans un anneau tampon.

    Retourne un dict de statistiques : nombre de cellules valides (non-NaN)
    dans l'AOI et dans le tampon, médiane/p10/p90 de la vitesse, et le masque
    AOI reprojeté. Un `n_aoi = 0` = aucun point EGMS sur la tourbière (attendu
    en PSI) ; regarder alors `n_buffer` pour le contexte de bordure.
    """
    from pyproj import Transformer
    from rasterio.features import geometry_mask
    from shapely.ops import transform as shp_transform

    from .aoi import load_aoi

    geom = load_aoi(cfg)
    crs = egms.rio.crs
    tr = Transformer.from_crs("EPSG:4326", crs, always_xy=True)
    geom_p = shp_transform(tr.transform, geom)
    geom_buf = geom_p.buffer(buffer_m)

    def mask_for(g):
        return ~geometry_mask([g], out_shape=egms.shape,
                              transform=egms.rio.transform(), invert=False)

    aoi_mask = mask_for(geom_p)
    buf_mask = mask_for(geom_buf) & ~aoi_mask

    vals = egms.values
    finite = np.isfinite(vals)

    def stats(m):
        v = vals[m & finite]
        if v.size == 0:
            return {"n": 0}
        return {"n": int(v.size), "median": float(np.median(v)),
                "p10": float(np.percentile(v, 10)),
                "p90": float(np.percentile(v, 90)),
                "min": float(v.min()), "max": float(v.max())}

    return {
        "aoi": stats(aoi_mask),
        "buffer": stats(buf_mask),
        "aoi_mask": xr.DataArray(aoi_mask, coords=egms.coords, dims=egms.dims),
        "buffer_mask": xr.DataArray(buf_mask, coords=egms.coords, dims=egms.dims),
    }


def verdict(stats: dict) -> str:
    """Interprétation prête à lire du contrôle EGMS."""
    a, b = stats["aoi"], stats["buffer"]
    if a["n"] == 0 and b["n"] == 0:
        return ("Aucun point EGMS sur l'AOI NI en bordure. Attendu en PSI sur "
                "tourbière — non discriminant seul, mais montre qu'aucun "
                "traitement PSI opérationnel ne trouve de cible ici.")
    if a["n"] == 0 and b["n"] > 0:
        return (f"0 point sur la tourbière, {b['n']} en bordure "
                f"(médiane {b.get('median', float('nan')):.1f} mm/an). "
                "Cohérent avec PSI = bâti/routes seulement ; donne le contexte "
                "régional mais ne mesure pas le tapis.")
    return (f"{a['n']} points EGMS sur la tourbière "
            f"(médiane {a.get('median', float('nan')):.1f} mm/an) — un tiers "
            "détecte un signal ici : à comparer directement à notre résultat.")


def load_egms_timeseries_csv(csv_path, cfg: dict | None = None,
                             buffer_m: float = 0.0) -> pd.DataFrame:
    """Charge un export CSV EGMS (points + séries temporelles) et filtre ceux
    dans l'AOI (+ tampon). Le CSV EGMS a des colonnes easting/northing (3035)
    ou lat/lon selon l'export, puis des colonnes de dates (YYYYMMDD) en mm."""
    from shapely.geometry import Point

    from .aoi import load_aoi

    df = pd.read_csv(csv_path)
    geom = load_aoi(cfg).buffer(buffer_m / 111_000.0)  # approx deg
    latc = next((c for c in df.columns if c.lower() in ("latitude", "lat")), None)
    lonc = next((c for c in df.columns if c.lower() in ("longitude", "lon")), None)
    if latc and lonc:
        inside = df.apply(lambda r: geom.contains(Point(r[lonc], r[latc])), axis=1)
        return df[inside].reset_index(drop=True)
    raise KeyError(f"colonnes lat/lon EGMS introuvables ; colonnes: {list(df.columns)[:8]}")

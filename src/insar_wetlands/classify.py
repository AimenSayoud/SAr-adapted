"""Phase 6 — Classification spatio-temporelle en 5 classes comportementales.

Regles explicites (defendables en these, contrairement a un clustering opaque) :
  A  sol ferme/stable peripherique : hors AOI, jamais inonde, coherence haute
  B  vegetation stable             : jamais inonde, coherence moyenne
  C  coeur de la tourbiere         : dans l'AOI, inondation rare
  D  transition temporairement inondee : inondation intermittente
  E  eau libre permanente          : inonde quasi en permanence
"""

from __future__ import annotations

import numpy as np
import xarray as xr

CLASS_NAMES = {1: "A_stable_ground", 2: "B_stable_vegetation", 3: "C_peat_core",
               4: "D_transition_flooded", 5: "E_permanent_water"}


def classify(flooded_frac: xr.DataArray, mean_coh: xr.DataArray,
             aoi: xr.DataArray, cfg: dict) -> xr.DataArray:
    c = cfg.get("classification", {})
    permanent = float(c.get("permanent_water_frac", 0.80))
    intermittent = float(c.get("intermittent_frac", 0.15))
    coh_stable = float(c.get("stable_coherence", 0.45))

    # Les classes hydro (C/D/E) n'ont de sens physique que DANS l'AOI : hors
    # AOI, un pixel intermittent est un fosse/champ agricole, pas une zone de
    # transition de la tourbiere. On borne D et E a l'AOI (C l'est deja), et
    # hors AOI tout devient A/B (sol/vegetation de reference).
    never = flooded_frac <= intermittent
    flooded_in = flooded_frac.where(aoi)

    cls = xr.full_like(flooded_frac, np.nan)
    cls = cls.where(~(aoi & (flooded_in >= permanent)), 5)                # E
    cls = cls.where(~(aoi & (flooded_frac > intermittent)
                      & (flooded_frac < permanent)), 4)                   # D
    cls = cls.where(~(aoi & never), 3)                                    # C
    cls = cls.where(~(~aoi & (mean_coh >= coh_stable)), 1)                # A
    cls = cls.where(~(~aoi & (mean_coh < coh_stable)), 2)                 # B
    cls.attrs = {"long_name": "Classe comportementale (1=A..5=E)"}
    return cls.rename("behavior_class")


def pick_reference_pixel(cls: xr.DataArray, mean_coh: xr.DataArray,
                         centroid_xy: tuple[float, float],
                         max_dist_m: float = 1000.0) -> dict:
    """Point de reference Classe A : coherence max a moins de max_dist_m du site.

    (Phase 11, pipeline 'calibration locale relative'.)
    """
    xx, yy = np.meshgrid(cls.x.values, cls.y.values)
    dist = np.hypot(xx - centroid_xy[0], yy - centroid_xy[1])
    cand = mean_coh.where((cls == 1) & (xr.DataArray(
        dist, coords={"y": cls.y, "x": cls.x}) <= max_dist_m))
    if cand.isnull().all():
        cand = mean_coh.where(cls == 1)  # fallback : n'importe quel Classe A
    flat = cand.stack(z=("y", "x"))
    best = flat.idxmax("z").item()
    y, x = best
    return {"y": float(y), "x": float(x),
            "coherence": float(cand.sel(y=y, x=x))}

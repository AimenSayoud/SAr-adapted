"""Phase 12 — Decomposition geometrique LOS -> vertical.

d_vert = d_LOS / cos(theta_inc), justifie par la poroelasticite du fen confine
(deplacement quasi exclusivement vertical : gonflement/tassement de la tourbe).

HyP3 fournit lv_theta : angle d'ELEVATION du vecteur de visee par rapport a
l'horizontale (radians). L'angle d'incidence est donc :
    theta_inc = pi/2 - lv_theta
"""

from __future__ import annotations

import numpy as np
import xarray as xr


def incidence_angle(lv_theta: xr.DataArray) -> xr.DataArray:
    inc = (np.pi / 2 - lv_theta).rename("incidence_rad")
    return inc


def los_to_vertical(ts_los_mm: xr.DataArray,
                    lv_theta: xr.DataArray) -> xr.DataArray:
    """Projette la serie LOS (mm) en deplacement vertical (mm), pixel a pixel."""
    inc = incidence_angle(lv_theta)
    if ts_los_mm.sizes.get("y") != inc.sizes.get("y"):
        inc = inc.interp(y=ts_los_mm.y, x=ts_los_mm.x, method="nearest")
    d_vert = (ts_los_mm / np.cos(inc)).rename("vertical_displacement_mm")
    d_vert.attrs["units"] = "mm"
    d_vert.attrs["note"] = ("hypothese: composante horizontale negligeable "
                            "(poroelasticite du fen confine)")
    return d_vert

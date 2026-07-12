"""Phase 5 — Masque dynamique de l'eau : hybride optique (NDWI/MNDWI) + radar.

Masque(t) = f(NDWI, MNDWI, sigma0_VV) pour chaque date radar t :
  - eau libre : NDWI/MNDWI eleves ET/OU sigma0 tres faible (reflexion speculaire) ;
  - eau sous vegetation : le radar 'triche' (double-bounce brillant) -> on
    garde un drapeau separe pour ces pixels ambigus.
L'optique etant asynchrone (nuages), chaque date S1 recoit l'etat optique le
plus proche (nearest) si distant de moins de `max_gap_days`.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import xarray as xr


def optical_to_s1_dates(s2_stack: xr.Dataset, s1_dates: pd.DatetimeIndex,
                        max_gap_days: int = 12) -> xr.Dataset:
    """Interpole (nearest) le stack optique sur les dates radar."""
    tol = pd.Timedelta(days=max_gap_days)
    out = s2_stack.reindex(time=s1_dates, method="nearest", tolerance=tol)
    gap = np.abs(pd.to_datetime(out.time.values)
                 - s2_stack.time.sel(time=s1_dates, method="nearest").values)
    out["optical_gap_days"] = ("time", (gap / pd.Timedelta("1D")).astype(float))
    return out


def water_mask(s2_on_s1: xr.Dataset, rtc: xr.Dataset, cfg: dict) -> xr.Dataset:
    """Masque d'eau par date S1 + drapeau double-bounce (vegetation inondee)."""
    m = cfg["masking"]
    ndwi_w = s2_on_s1["ndwi"] > m["ndwi_water_threshold"]
    mndwi_w = (s2_on_s1["mndwi"] > m["mndwi_water_threshold"]
               if "mndwi" in s2_on_s1 else xr.zeros_like(ndwi_w))
    optical_water = ndwi_w | mndwi_w

    g0 = rtc["gamma0_vv_db"].reindex(time=optical_water.time, method="nearest",
                                     tolerance=pd.Timedelta(days=1))
    radar_dark = g0 < m["sigma0_vv_water_db"]

    # Double-bounce : fort rehaussement par rapport a la mediane temporelle du
    # pixel => vegetation emergente probablement inondee (eau cachee).
    g0_med = g0.median("time")
    double_bounce = (g0 - g0_med) > m.get("double_bounce_delta_db", 4.0)

    water = optical_water.fillna(False) | radar_dark.fillna(False)
    hidden = double_bounce.fillna(False) & ~water
    return xr.Dataset({
        "water": water,
        "hidden_water": hidden,
        "water_or_hidden": water | hidden,
        "optical_valid": s2_on_s1["ndwi"].notnull(),
    })


def flooded_fraction(mask: xr.Dataset, var: str = "water_or_hidden") -> xr.DataArray:
    """Fraction du temps ou chaque pixel est inonde (base de la Phase 6/7)."""
    return mask[var].mean("time").rename("flooded_fraction")

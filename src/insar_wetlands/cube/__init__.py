"""The data cube: every dataset of the project on one space, one time base and one set of units, for fusion v3.

Three layers, all under the hub's ``06_data/cube`` (private: field values enter from the silver layer on):

* ``raw/``    public downloads as received (Sentinel-2, Landsat, WorldCover, ERA5-Land, IMGW, LiDAR, NISAR); never edited.
* ``silver/`` harmonised: interferogram and backscatter stacks on the 40 m template (and 2026 at 20 m), optical indices,
              static pixel layers, the hourly site table with its flags, the cross-validation folds.
* ``gold/``   model-ready tables: per plot (acquisitions, pairs, truth) and per mat/reference pixel (pairs, acquisitions).

Rules the code enforces, so no model re-solves them:

* time is UTC; an overpass is an instant (ascending ≈ 16:36, descending ≈ 05:09); joins across sources are as-of joins
  with a tolerance and an explicit lag column, never a silent nearest;
* phase is stored in radians, unreferenced (wrapped phasor, unwrapped, connected component); the reference is chosen when
  a model reads it; conversion to millimetres goes through one constant (``core.WAVELENGTH_MM``);
* flags are columns with their source, never filtered-out rows.

This package is public; it contains no field values. ``scripts/cube_build.py`` runs it.
"""
from __future__ import annotations

import os
from pathlib import Path


def cube_root(root: str | Path | None = None) -> Path:
    """The hub's ``06_data/cube``: explicit argument, ``RZECIN_CUBE_ROOT``, or the hub containing this repository."""
    if root is not None:
        return Path(root)
    if os.environ.get("RZECIN_CUBE_ROOT"):
        return Path(os.environ["RZECIN_CUBE_ROOT"])
    return Path(__file__).resolve().parents[5] / "06_data" / "cube"

"""Phase 7 — Indice de qualite probabiliste W_i et matrice de ponderation.

Version corrigee (evite la double penalisation des pixels inondes) :

    W_i = coherence_conditionnelle_hors_eau_i x (1 - fraction_inondee_i)

ou la coherence conditionnelle est moyennee UNIQUEMENT sur les paires dont les
deux dates sont hors-eau pour ce pixel. Un pixel du fen tres coherent en
saison seche garde ainsi un poids honnete au lieu d'etre ecrase par ses
mauvaises paires d'hiver/inondation.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import xarray as xr


def pair_dry_mask(water: xr.DataArray, pairs: list[str]) -> xr.DataArray:
    """(pair, y, x) True si le pixel est HORS-eau aux deux dates de la paire."""
    dry = ~water
    ref = pd.to_datetime([p.split("_")[0] for p in pairs])
    sec = pd.to_datetime([p.split("_")[1] for p in pairs])
    d_ref = dry.reindex(time=ref, method="nearest").rename({"time": "pair"})
    d_sec = dry.reindex(time=sec, method="nearest").rename({"time": "pair"})
    out = (d_ref.values & d_sec.values)
    return xr.DataArray(out, dims=("pair", "y", "x"),
                        coords={"pair": pairs, "y": water.y, "x": water.x})


def quality_index(corr_stack: xr.DataArray, water: xr.DataArray,
                  flooded_frac: xr.DataArray) -> xr.Dataset:
    """W_i + composantes. corr_stack: (pair,y,x) ; water: (time,y,x)."""
    dry_pairs = pair_dry_mask(water, list(corr_stack.pair.values))
    coh_cond = corr_stack.where(dry_pairs).mean("pair")
    coh_all = corr_stack.mean("pair")
    n_dry = dry_pairs.sum("pair").astype(float)
    w = (coh_cond.fillna(0) * (1.0 - flooded_frac)).clip(0, 1)
    return xr.Dataset({
        "W": w.rename("W"),
        "coh_conditional_dry": coh_cond,
        "coh_all_pairs": coh_all,
        "n_dry_pairs": n_dry,
        "flooded_fraction": flooded_frac,
    })


def weight_matrix_for_pixel(corr_pix: np.ndarray, dry_pix: np.ndarray,
                            gamma_min: float = 0.30) -> np.ndarray:
    """Poids par paire pour l'inversion ISBAS d'un pixel donne.

    poids = coherence si (paire seche ET coherence >= gamma_min), sinon 0.
    """
    w = np.where(dry_pix & (corr_pix >= gamma_min), corr_pix, 0.0)
    return w

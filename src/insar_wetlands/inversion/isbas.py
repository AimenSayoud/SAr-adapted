"""Phase 9 — Inversion ISBAS customisee (pixels intermittents).

Contrairement au SBAS (qui exige un pixel coherent sur TOUTES les paires),
l'ISBAS resout, pixel par pixel, le sous-systeme forme des seules paires
'valides' pour ce pixel (coherentes ET hors-eau aux deux dates), des lors que
N >= min_pairs et que le sous-graphe temporel reste connexe.

Moindres carres ponderes : W = diag(poids Phase 7) ; increments de phase entre
dates consecutives, rattaches au point de reference.
"""

from __future__ import annotations

import numpy as np
import xarray as xr

from ..network import design_matrix
from ..stack import dates_from_pairs

WAVELENGTH_M = 0.0555465763  # Sentinel-1 C-band
PHASE_TO_MM = -WAVELENGTH_M / (4 * np.pi) * 1000.0  # rad -> mm LOS (convention MintPy)


def invert_pixel(phase: np.ndarray, A: np.ndarray, w: np.ndarray,
                 min_pairs: int = 20) -> tuple[np.ndarray | None, float, int]:
    """Inversion WLS d'un pixel. Retourne (increments, rms_residu, n_pairs)."""
    ok = (w > 0) & np.isfinite(phase)
    n = int(ok.sum())
    if n < min_pairs:
        return None, np.nan, n
    Ai, pi, wi = A[ok], phase[ok], w[ok]
    # connexite du sous-graphe : chaque increment doit etre couvert
    covered = (Ai != 0).any(axis=0)
    if not covered.all():
        # increments non couverts -> systeme singulier ; on resout au sens
        # des moindres carres regularise minimal (les trous restent NaN)
        pass
    sw = np.sqrt(wi)
    sol, *_ = np.linalg.lstsq(Ai * sw[:, None], pi * sw, rcond=None)
    sol = np.where(covered, sol, np.nan)
    resid = pi - Ai @ np.nan_to_num(sol)
    rms = float(np.sqrt(np.nanmean(resid ** 2)))
    return sol, rms, n


def invert_stack(unw: xr.DataArray, corr: xr.DataArray,
                 dry_pairs: xr.DataArray, ref_yx: tuple[float, float],
                 min_pairs: int = 20, gamma_min: float = 0.30) -> xr.Dataset:
    """Inversion ISBAS sur tout le stack. Grille ~60x50 px : boucle OK.

    unw, corr, dry_pairs : (pair, y, x) alignes. ref_yx : point de reference
    (coordonnees UTM y, x) soustrait de chaque interferogramme.
    """
    pairs = list(unw.pair.values)
    dates = dates_from_pairs(pairs)
    A = design_matrix(pairs, dates)

    # Rattachement au point de reference (Classe A)
    ref = unw.sel(y=ref_yx[0], x=ref_yx[1], method="nearest")
    phase = (unw - ref).values          # (pair, y, x)
    coh = corr.values
    dry = dry_pairs.values
    ny, nx = phase.shape[1:]
    n_inc = len(dates) - 1

    incs = np.full((n_inc, ny, nx), np.nan, dtype="float32")
    rms = np.full((ny, nx), np.nan, dtype="float32")
    npairs = np.zeros((ny, nx), dtype="int16")
    for iy in range(ny):
        for ix in range(nx):
            w = np.where(dry[:, iy, ix] & (coh[:, iy, ix] >= gamma_min),
                         coh[:, iy, ix], 0.0)
            sol, r, n = invert_pixel(phase[:, iy, ix], A, w, min_pairs)
            npairs[iy, ix] = n
            if sol is not None:
                incs[:, iy, ix] = sol
                rms[iy, ix] = r

    cum = np.concatenate([np.zeros((1, ny, nx), dtype="float32"),
                          np.nancumsum(incs, axis=0)], axis=0)
    cum[:, np.isnan(rms)] = np.nan
    ts_mm = cum * PHASE_TO_MM
    return xr.Dataset(
        {
            "los_displacement_mm": (("time", "y", "x"), ts_mm),
            "rms_residual_rad": (("y", "x"), rms),
            "n_valid_pairs": (("y", "x"), npairs),
        },
        coords={"time": dates, "y": unw.y, "x": unw.x},
    )


def phase_closure(unw: xr.DataArray, max_triplets: int = 500) -> xr.DataArray:
    """Controle des sauts de phase : fermeture des triplets sur phase deroulee.

    c = phi_ij + phi_jk - phi_ik ; toute valeur ~ 2*pi*k (k!=0) signale une
    erreur de deroulement SNAPHU. Retourne la fraction de triplets fautifs
    par pixel (a croiser avec les pixels 'sauves' par l'ISBAS, Phase 10).
    """
    pairs = list(unw.pair.values)
    pset = {p: i for i, p in enumerate(pairs)}
    dates = dates_from_pairs(pairs)
    triplets = []
    for i, d1 in enumerate(dates):
        for d2 in dates[i + 1:]:
            for d3 in dates[i + 2:]:
                if d3 <= d2:
                    continue
                p12 = f"{d1:%Y%m%d}_{d2:%Y%m%d}"
                p23 = f"{d2:%Y%m%d}_{d3:%Y%m%d}"
                p13 = f"{d1:%Y%m%d}_{d3:%Y%m%d}"
                if p12 in pset and p23 in pset and p13 in pset:
                    triplets.append((pset[p12], pset[p23], pset[p13]))
        if len(triplets) >= max_triplets:
            break
    if not triplets:
        raise ValueError("aucun triplet ferme dans le reseau")
    data = unw.values
    bad = np.zeros(data.shape[1:], dtype="float32")
    tot = np.zeros(data.shape[1:], dtype="float32")
    for a, b, c in triplets:
        closure = data[a] + data[b] - data[c]
        valid = np.isfinite(closure)
        bad += (np.abs(closure) > np.pi).astype("float32") * valid
        tot += valid
    with np.errstate(invalid="ignore", divide="ignore"):
        # tot==0 (pixel hors emprise) -> 0/0 ; np.where choisit nan ensuite,
        # mais bad/tot est evalue AVANT le where -> avertissement inoffensif
        # a silencer explicitement plutot que de laisser echapper.
        frac = np.where(tot > 0, bad / tot, np.nan)
    return xr.DataArray(frac, coords={"y": unw.y, "x": unw.x},
                        name="closure_error_fraction",
                        attrs={"n_triplets": len(triplets)})

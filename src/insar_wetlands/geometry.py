"""Phase 12 — Decomposition geometrique LOS -> vertical. Phase 16 — 2-LOS.

d_vert = d_LOS / cos(theta_inc), justifie par la poroelasticite du fen confine
(deplacement quasi exclusivement vertical : gonflement/tassement de la tourbe).

HyP3 fournit lv_theta : angle d'ELEVATION du vecteur de visee par rapport a
l'horizontale (radians). L'angle d'incidence est donc :
    theta_inc = pi/2 - lv_theta

Phase 16 removes the pure-vertical assumption above by combining two
geometries (ascending + descending) instead of asserting it. A single LOS
measurement is one equation in two unknowns (vertical, east — north is not
observable from any single- or dual-track C-band SAR combination and is
assumed negligible, as elsewhere in this project); two geometries make it
solvable. `two_los_decompose` is deliberately generic over what "LOS
measurement" means: called on raw per-date series it decomposes displacement
directly, called on a seasonal harmonic fit's (a, b) cosine/sine coefficients
(`aggregate.seasonal_amplitude`) it decomposes the seasonal cycle itself — the
system is linear, so the same solve applies to either.
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


def los_unit_vector(incidence_deg: float, heading_deg: float) -> tuple[float, float]:
    """(e_east, e_up): a right-looking SAR track's LOS unit vector, decomposed
    onto East and Up (North is not resolved by this pair and is assumed
    negligible — see module docstring).

    `incidence_deg` from vertical, `heading_deg` clockwise from North (the
    satellite's flight direction) — both as already recorded per track in
    `config.yaml`'s `sentinel1.tracks`. With the geometry measured from the HyP3
    look vectors (X-042): ascending (32.26°, 348.6°) gives (-0.5232, 0.8456) and
    descending (39.17°, 190.1°) gives (+0.6218, 0.7753) — matching the unit vectors
    read directly from lv_theta/lv_phi over zone A to 1e-4.
    """
    theta = np.radians(incidence_deg)
    heading = np.radians(heading_deg)
    e_up = float(np.cos(theta))
    e_east = float(-np.sin(theta) * np.cos(heading))
    return e_east, e_up


def two_los_design_matrix(geom_asc: tuple[float, float],
                          geom_desc: tuple[float, float]) -> np.ndarray:
    """The 2x2 matrix mapping (d_east, d_vert) -> (los_asc, los_desc).

    `geom_asc`/`geom_desc` are (incidence_deg, heading_deg) pairs. Well
    conditioned here (~1.57) because the two tracks' East sensitivities have
    opposite sign (ascending looks ENE, descending looks WNW) — this is what
    makes the system solvable without imposing `d_east = 0`.
    """
    e_e_asc, e_u_asc = los_unit_vector(*geom_asc)
    e_e_desc, e_u_desc = los_unit_vector(*geom_desc)
    return np.array([[e_e_asc, e_u_asc], [e_e_desc, e_u_desc]])


def two_los_decompose(los_asc, los_desc,
                      geom_asc: tuple[float, float],
                      geom_desc: tuple[float, float]):
    """Solve for (d_east, d_vert) from LOS measurements in two geometries.

    `los_asc`/`los_desc` are scalars or equal-length arrays — the same LOS
    quantity from each track (raw displacement, or a seasonal fit's harmonic
    coefficient; the solve is linear so either is valid, see module
    docstring). Returns `(d_east, d_vert)` in the same shape as the inputs.
    """
    A = two_los_design_matrix(geom_asc, geom_desc)
    stacked = np.stack([np.asarray(los_asc, dtype=float),
                        np.asarray(los_desc, dtype=float)], axis=0)
    d_east, d_vert = np.linalg.solve(A, stacked)
    if np.ndim(los_asc) == 0:
        return float(d_east), float(d_vert)
    return d_east, d_vert


def two_los_amplitude_uncertainty(fit_asc: dict, fit_desc: dict,
                                  geom_asc: tuple[float, float],
                                  geom_desc: tuple[float, float],
                                  n_trials: int = 5000,
                                  rng: np.random.Generator | None = None) -> dict:
    """Monte Carlo propagation of each track's seasonal-fit uncertainty
    through the 2-LOS decomposition, onto vertical and east amplitude.

    A deterministic `two_los_decompose` call on two point estimates answers
    "what does the geometry imply", nothing about how much to trust it — a
    decomposition anchored on a fit through mostly noise (low
    `r2_seasonal`, wide `cov_a_b`) looks exactly as precise as one anchored
    on a strong fit unless this is done. Draws each track's (a_cos_mm,
    b_sin_mm) from `N(mean, cov_a_b)` (`aggregate.seasonal_amplitude`'s OLS
    parameter covariance), decomposes every draw, and summarizes the
    resulting vertical/east amplitude distributions.

    `fit_asc`/`fit_desc` are `seasonal_amplitude(...)` result dicts (must
    include `a_cos_mm`, `b_sin_mm`, `cov_a_b`). Returns a dict with
    `vertical_amplitude_mm` and `east_amplitude_mm`, each
    `{median, mean, std, ci95: [lo, hi]}`.
    """
    rng = rng or np.random.default_rng(0)
    mean_asc = [fit_asc["a_cos_mm"], fit_asc["b_sin_mm"]]
    mean_desc = [fit_desc["a_cos_mm"], fit_desc["b_sin_mm"]]
    draws_asc = rng.multivariate_normal(mean_asc, fit_asc["cov_a_b"], size=n_trials)
    draws_desc = rng.multivariate_normal(mean_desc, fit_desc["cov_a_b"], size=n_trials)

    a_east, a_vert = two_los_decompose(draws_asc[:, 0], draws_desc[:, 0], geom_asc, geom_desc)
    b_east, b_vert = two_los_decompose(draws_asc[:, 1], draws_desc[:, 1], geom_asc, geom_desc)
    vert_amplitude = np.hypot(a_vert, b_vert)
    east_amplitude = np.hypot(a_east, b_east)

    def summarize(x: np.ndarray) -> dict:
        return {"median": round(float(np.median(x)), 3),
                "mean": round(float(np.mean(x)), 3),
                "std": round(float(np.std(x)), 3),
                "ci95": [round(float(np.percentile(x, 2.5)), 3),
                        round(float(np.percentile(x, 97.5)), 3)]}

    return {"vertical_amplitude_mm": summarize(vert_amplitude),
            "east_amplitude_mm": summarize(east_amplitude),
            "n_trials": n_trials}

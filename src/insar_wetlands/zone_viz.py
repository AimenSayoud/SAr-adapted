"""Visualisation et VALIDATION VISUELLE des zones A/B/C/D.

Objectif : montrer que les zones ne sont pas un découpage arbitraire mais des
**unités physiques réelles**, par plusieurs voies indépendantes — carte
catégorielle, contours sur fonds physiques (cohérence, σ0), distributions par
zone, profil radial au bord, et **contrôle de surface** (A+B doit retrouver les
~89.7 ha connus de la tourbière : une vérification OBJECTIVE, non visuelle,
que le masque correspond bien au site réel).

Codes couleur constants dans tout le projet :
A=rouge (tapis) · B=bleu (lac) · C=vert (prairie appariée) · D=gris (autres).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import xarray as xr

ZONE_ORDER = ("A", "B", "C", "D")
ZONE_COLORS = {"A": "#d62728", "B": "#1f77b4", "C": "#2ca02c", "D": "#9e9e9e"}
ZONE_LABELS = {"A": "A — tapis flottant", "B": "B — lac résiduel",
               "C": "C — prairie stable appariée", "D": "D — autres couverts"}


def zone_label_array(zones: dict, template: xr.DataArray) -> xr.DataArray:
    """Carte catégorielle : 1=A, 2=B, 3=C, 4=D, NaN ailleurs."""
    lab = np.full(template.shape, np.nan, "float32")
    for i, z in enumerate(ZONE_ORDER, start=1):
        if z in zones:
            lab[zones[z].values] = i
    return xr.DataArray(lab, coords=template.coords, dims=template.dims,
                        name="zone")


def pixel_area_m2(template: xr.DataArray) -> float:
    dx = abs(float(template.x[1] - template.x[0]))
    dy = abs(float(template.y[1] - template.y[0]))
    return dx * dy


def zone_areas(zones: dict, template: xr.DataArray,
               expected_site_ha: float | None = 89.7) -> pd.DataFrame:
    """Effectifs et SURFACES par zone (+ contrôle vs la surface connue du site).

    **Validation objective du masque** : A+B (intérieur du polygone, végétalisé
    + en eau) doit retrouver la surface documentée de la tourbière. Un écart
    important signalerait un masque mal géoréférencé — ce qu'aucune inspection
    visuelle ne garantit."""
    a = pixel_area_m2(template)
    rows = []
    for z in ZONE_ORDER:
        if z in zones:
            n = int(zones[z].sum())
            rows.append({"zone": z, "label": ZONE_LABELS[z], "n_px": n,
                         "area_ha": round(n * a / 1e4, 2)})
    df = pd.DataFrame(rows)
    if expected_site_ha and {"A", "B"} <= set(df.zone):
        inside = float(df[df.zone.isin(["A", "B"])].area_ha.sum())
        df.attrs["inside_ha"] = round(inside, 2)
        df.attrs["expected_ha"] = expected_site_ha
        df.attrs["area_error_pct"] = round(
            100 * (inside - expected_site_ha) / expected_site_ha, 1)
    return df


def zone_field_table(field: xr.DataArray, zones: dict,
                     name: str = "field") -> pd.DataFrame:
    """Médiane / quartiles d'un champ par zone (comparaison chiffrée)."""
    rows = []
    for z in ZONE_ORDER:
        if z not in zones:
            continue
        v = field.values[zones[z].values]
        v = v[np.isfinite(v)]
        if not v.size:
            continue
        rows.append({"zone": z, "field": name, "n": int(v.size),
                     "median": round(float(np.median(v)), 4),
                     "p25": round(float(np.percentile(v, 25)), 4),
                     "p75": round(float(np.percentile(v, 75)), 4)})
    return pd.DataFrame(rows)


def plot_zone_map(zones: dict, template: xr.DataArray, ax=None,
                  title: str = "Zones A/B/C/D"):
    """Vue 1 — carte catégorielle des quatre zones."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import BoundaryNorm, ListedColormap
    from matplotlib.patches import Patch

    if ax is None:
        _, ax = plt.subplots(figsize=(6, 5))
    lab = zone_label_array(zones, template)
    cmap = ListedColormap([ZONE_COLORS[z] for z in ZONE_ORDER])
    ax.imshow(lab.values, cmap=cmap, norm=BoundaryNorm([.5, 1.5, 2.5, 3.5, 4.5],
                                                       cmap.N),
              interpolation="nearest")
    ax.set_title(title)
    ax.legend(handles=[Patch(color=ZONE_COLORS[z], label=ZONE_LABELS[z])
                       for z in ZONE_ORDER if z in zones],
              fontsize=7, loc="upper right", framealpha=.9)
    return ax


def plot_zones_over_field(field: xr.DataArray, zones: dict, ax=None,
                          title: str = "", cmap: str = "viridis",
                          which=("A", "B", "C"), robust: bool = True):
    """Vue 2 — CONTOURS des zones sur un fond physique.

    La démonstration la plus convaincante : si les contours (tracés a partir du
    geojson + WorldCover + S2) épousent une structure visible dans un champ
    **radar indépendant** (cohérence, σ0), les zones sont des unités réelles et
    non un découpage plaqué."""
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots(figsize=(6, 5))
    v = field.values.astype(float)
    lo, hi = (np.nanpercentile(v, [2, 98]) if robust
              else (np.nanmin(v), np.nanmax(v)))
    im = ax.imshow(v, cmap=cmap, vmin=lo, vmax=hi, interpolation="nearest")
    for z in which:
        if z in zones:
            ax.contour(zones[z].values.astype(float), levels=[0.5],
                       colors=[ZONE_COLORS[z]], linewidths=1.6)
    ax.set_title(title)
    plt.colorbar(im, ax=ax, fraction=0.046)
    return ax


def plot_zone_distributions(field: xr.DataArray, zones: dict, ax=None,
                            title: str = "", xlabel: str = ""):
    """Vue 3 — distributions par zone (preuve STATISTIQUE de la séparation).

    Une carte peut tromper l'œil ; des distributions séparées, non."""
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots(figsize=(6, 4))
    data, labels, colors = [], [], []
    for z in ZONE_ORDER:
        if z not in zones:
            continue
        v = field.values[zones[z].values]
        v = v[np.isfinite(v)]
        if v.size:
            data.append(v); labels.append(z); colors.append(ZONE_COLORS[z])
    bp = ax.boxplot(data, labels=labels, patch_artist=True, showfliers=False)
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c); patch.set_alpha(.6)
    ax.set_ylabel(xlabel); ax.set_title(title); ax.grid(alpha=.3, axis="y")
    return ax


def plot_radial(field: xr.DataArray, signed_dist: xr.DataArray, ax=None,
                edges=None, title: str = "", ylabel: str = ""):
    """Vue 4 — profil radial : le champ selon la distance signée au bord.

    Une **marche** au passage du bord (distance 0) prouve une frontière
    physique ; une variation lisse indiquerait au contraire un gradient diffus
    et donc un découpage arbitraire."""
    import matplotlib.pyplot as plt

    from .stratify import radial_profile

    if ax is None:
        _, ax = plt.subplots(figsize=(6, 4))
    if edges is None:
        edges = np.arange(-400, 1401, 100.0)
    prof = radial_profile(field, signed_dist, edges)
    col = "median" if "median" in prof else prof.columns[-1]
    ax.plot(prof["dist_center_m"], prof[col], "o-", color="k", ms=4)
    ax.axvline(0, ls="--", c="r", lw=1.5)
    ax.text(0, ax.get_ylim()[1], " bord", color="r", va="top", fontsize=8)
    ax.set_xlabel("distance signée au bord (m) — négatif = DEDANS")
    ax.set_ylabel(ylabel); ax.set_title(title); ax.grid(alpha=.3)
    return ax

"""Stratification intérieur/extérieur de la tourbière — le tapis flottant
décorrèle-t-il PLUS que la même végétation sur sol stable ?

Compare le comportement InSAR (cohérence, décroissance temporelle, résidu
d'inversion) entre 4 zones du crop rectangulaire :
  A  tapis végétalisé (intérieur polygone, hors eau)   -> le sujet
  B  lac résiduel (intérieur, eau)                     -> plancher (contrôle -)
  C  végétation extérieure APPARIÉE à A (même couvert) -> le contrôle clé
  D  autres couverts extérieurs (forêt, sol nu)        -> contexte / plafond

Idée statistique centrale : comparer A et C **par interférogramme** (chaque
paire vue dans les deux zones) neutralise les confusions au niveau paire
(baseline perpendiculaire, atmosphère du jour) — il ne reste que la
différence de surface/couvert. Si, à couvert égal (WorldCover + features S2),
A décorrèle systématiquement plus que C, l'instabilité du tapis flottant est
un facteur causal propre, distinct de « la végétation en bande C ».
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import xarray as xr

# ESA WorldCover : classes utiles ici
WC_TREE, WC_SHRUB, WC_GRASS, WC_CROP = 10, 20, 30, 40
WC_BUILT, WC_BARE, WC_WATER, WC_WETLAND, WC_MOSS = 50, 60, 80, 90, 100


def s2_landcover_features(s2: xr.Dataset) -> xr.Dataset:
    """Features de couvert par pixel, robustes au contenu réel du stack S2.

    Deux formats gérés :
      - bandes brutes (green, nir, swir16) -> GNDVI et NDMI directs ;
      - indices déjà calculés (ndwi, mndwi) : notre s2_stack.nc ne stocke que
        ceux-là. Or GNDVI = (nir-green)/(nir+green) = -NDWI, donc la verdure se
        reconstruit exactement depuis -NDWI ; l'humidité depuis MNDWI (swir).

    - greenness_mean : verdure moyenne (GNDVI ≈ -NDWI ; végétation -> positif)
    - greenness_amp  : amplitude saisonnière (p90-p10) -> phénologie
    - wetness_mean   : humidité moyenne (MNDWI, ou NDMI si bandes dispo)
    """
    if "green" in s2 and "nir" in s2:
        g = s2["green"].astype("float32"); n = s2["nir"].astype("float32")
        gndvi = (n - g) / (n + g)
    elif "ndwi" in s2:
        gndvi = -s2["ndwi"].astype("float32")          # GNDVI = -NDWI
    else:
        raise KeyError(f"stack S2 sans green/nir ni ndwi ; variables: {list(s2.data_vars)}")

    if "swir16" in s2 and "nir" in s2:
        n2 = s2["nir"].astype("float32"); sw = s2["swir16"].astype("float32")
        wet = (n2 - sw) / (n2 + sw)                    # NDMI
    elif "mndwi" in s2:
        wet = s2["mndwi"].astype("float32")
    else:
        wet = -gndvi                                   # repli : NDWI comme proxy humidité

    q90 = gndvi.quantile(0.9, "time").drop_vars("quantile")
    q10 = gndvi.quantile(0.1, "time").drop_vars("quantile")
    feat = xr.Dataset({
        "greenness_mean": gndvi.mean("time", skipna=True),
        "greenness_amp": q90 - q10,
        "wetness_mean": wet.mean("time", skipna=True),
    })
    for v in feat.data_vars:
        feat[v].attrs = {}
    return feat


def worldcover_tile_name(lat: float, lon: float, year: int = 2021) -> str:
    """Nom de la tuile ESA WorldCover 3°x3° contenant (lat, lon)."""
    t_lat = int(np.floor(lat / 3.0) * 3)
    t_lon = int(np.floor(lon / 3.0) * 3)
    ns = f"N{t_lat:02d}" if t_lat >= 0 else f"S{-t_lat:02d}"
    ew = f"E{t_lon:03d}" if t_lon >= 0 else f"W{-t_lon:03d}"
    return f"ESA_WorldCover_10m_{year}_v200_{ns}{ew}_Map.tif"


def load_worldcover(template: xr.DataArray, cfg: dict | None = None,
                    year: int = 2021, cache_dir="/content/worldcover") -> xr.DataArray:
    """Charge ESA WorldCover 10 m (label de couvert indépendant) sur la grille
    du crop.

    ATTENTION MÉMOIRE : une tuile WorldCover fait 3°x3° = ~33000x33000 px
    (~1 Go). Il ne faut JAMAIS la reprojeter entière (satur RAM). On DÉCOUPE
    d'abord à la bbox du site (lecture fenêtrée COG -> ~500x500 px) puis on
    reprojette ce petit extrait. Lecture directe via /vsicurl (pas de
    téléchargement du Go entier) ; repli sur téléchargement si /vsicurl échoue.
    """
    from pathlib import Path
    import rioxarray  # noqa: F401
    from .config import load_config
    from .aoi import buffered_bbox

    cfg = cfg or load_config()
    lon, lat = cfg["site"]["centroid"]
    tile = worldcover_tile_name(lat, lon, year)
    base = (f"https://esa-worldcover.s3.eu-central-1.amazonaws.com/"
            f"v200/{year}/map/{tile}")
    minx, miny, maxx, maxy = buffered_bbox(cfg)   # degrés (CRS de la tuile)

    def _open_clip(src):
        da = rioxarray.open_rasterio(src, masked=True, chunks=True).squeeze("band", drop=True)
        return da.rio.clip_box(minx, miny, maxx, maxy)   # fenêtre COG -> minuscule

    try:
        wc_small = _open_clip("/vsicurl/" + base)
    except Exception:
        cache = Path(cache_dir); cache.mkdir(parents=True, exist_ok=True)
        local = cache / tile
        if not local.exists():
            import urllib.request
            urllib.request.urlretrieve(base, local)
        wc_small = _open_clip(local)

    wc_small = wc_small.load()                    # petit extrait : OK en RAM
    return wc_small.rio.reproject_match(template).rename("worldcover")


def dominant_class(worldcover: xr.DataArray, mask: xr.DataArray) -> int:
    """Classe WorldCover dominante sous `mask` (ex: le couvert du tapis)."""
    v = worldcover.where(mask).values
    v = v[np.isfinite(v)].astype(int)
    if v.size == 0:
        return -1
    return int(np.bincount(v).argmax())


def _feature_range(feat: xr.Dataset, mask: xr.DataArray,
                   q=(0.1, 0.9)) -> dict:
    """Intervalle [q10, q90] de chaque feature sous `mask` (gabarit de A)."""
    out = {}
    for v in feat.data_vars:
        vals = feat[v].where(mask).values
        vals = vals[np.isfinite(vals)]
        if vals.size:
            out[v] = (float(np.quantile(vals, q[0])), float(np.quantile(vals, q[1])))
    return out


def define_zones(template: xr.DataArray, cfg: dict, flooded_frac: xr.DataArray,
                 worldcover: xr.DataArray | None = None,
                 s2_feat: xr.Dataset | None = None,
                 dem: xr.DataArray | None = None,
                 water_thresh: float = 0.30, slope_max_deg: float = 5.0,
                 match_features: bool = True) -> dict:
    """Construit les masques A/B/C/D.

    C (végétation extérieure appariée) combine DEUX critères, comme demandé :
      - même classe WorldCover dominante que A (si worldcover fourni) ;
      - features S2 (verdure/amplitude/humidité) dans l'intervalle [p10,p90]
        de A (si s2_feat fourni et match_features).
    Contrôle topographique : exclusion des pentes > slope_max_deg (DEM), pour
    éviter la décorrélation géométrique de la forêt en pente.
    """
    from .stack import aoi_mask

    inside = aoi_mask(template, cfg)
    ff = flooded_frac.reindex_like(template, method="nearest") \
        if flooded_frac.shape != template.shape else flooded_frac
    water = ff > water_thresh
    outside = ~inside

    A = inside & ~water
    B = inside & water

    # contrôle de pente
    flat = xr.ones_like(template, dtype=bool)
    if dem is not None:
        gy, gx = np.gradient(dem.values)
        px = abs(float(template.x[1] - template.x[0]))
        slope = np.degrees(np.arctan(np.sqrt(gx**2 + gy**2) / px))
        flat = xr.DataArray(slope <= slope_max_deg, coords=template.coords, dims=template.dims)

    cand = outside & ~water & flat
    C = cand
    info = {}
    if worldcover is not None:
        dom = dominant_class(worldcover, A)
        info["dominant_class"] = dom
        C = C & (worldcover == dom)
    if s2_feat is not None and match_features:
        rng = _feature_range(s2_feat, A)
        info["feature_range_A"] = rng
        for v, (lo, hi) in rng.items():
            C = C & (s2_feat[v] >= lo) & (s2_feat[v] <= hi)

    D = outside & ~water & flat & ~C
    for name, m in [("A", A), ("B", B), ("C", C), ("D", D)]:
        info[f"n_{name}"] = int(m.sum())
    return {"A": A, "B": B, "C": C, "D": D, "info": info}


def coherence_by_zone_perpair(corr: xr.DataArray, zones: dict) -> pd.DataFrame:
    """Cohérence moyenne PAR PAIRE dans chaque zone (+ dt, saison).

    Colonnes : zone, pair, dt_days, season, mean_coh, n_px. La forme 'longue'
    permet l'analyse de décroissance ET la comparaison appariée A vs C.
    """
    season_by_month = {12: "hiver", 1: "hiver", 2: "hiver", 3: "printemps",
                       4: "printemps", 5: "printemps", 6: "été", 7: "été",
                       8: "été", 9: "automne", 10: "automne", 11: "automne"}
    rows = []
    zmasks = {k: zones[k].values for k in ("A", "B", "C", "D")}
    for p in corr.pair.values:
        ref, sec = str(p).split("_")
        dt = (pd.Timestamp(sec) - pd.Timestamp(ref)).days
        c = corr.sel(pair=p).values
        for z, m in zmasks.items():
            vals = c[m & np.isfinite(c)]
            if vals.size:
                rows.append({"zone": z, "pair": str(p), "dt_days": dt,
                             "season": season_by_month[pd.Timestamp(ref).month],
                             "mean_coh": float(vals.mean()), "n_px": int(vals.size)})
    return pd.DataFrame(rows)


def coherence_by_zone_stream(cropped_root, pairs, zones: dict,
                             template: xr.DataArray) -> tuple[pd.DataFrame, xr.DataArray]:
    """Version ÉCONOME EN MÉMOIRE de coherence_by_zone_perpair.

    Charge la cohérence UNE PAIRE À LA FOIS depuis le disque (jamais les ~349
    couches en RAM simultanément — ce qui faisait sauter Colab), accumule les
    stats par zone et la carte de cohérence moyenne (somme/compteur courants),
    puis libère chaque couche. Retourne (df_long, coh_mean).
    """
    from .stack import load_layer

    season_by_month = {12: "hiver", 1: "hiver", 2: "hiver", 3: "printemps",
                       4: "printemps", 5: "printemps", 6: "été", 7: "été",
                       8: "été", 9: "automne", 10: "automne", 11: "automne"}
    zmasks = {k: zones[k].values for k in ("A", "B", "C", "D")}
    tshape = template.shape
    rows = []
    ssum = np.zeros(tshape, dtype="float64")
    scnt = np.zeros(tshape, dtype="int32")
    for p in pairs:
        try:
            da = load_layer(cropped_root, "corr", [p]).isel(pair=0)
        except Exception:
            continue
        if da.shape != tshape:                      # grilles legerement differentes
            da = da.rio.reproject_match(template)
        c = da.values.astype("float32")
        del da
        finite = np.isfinite(c)
        ssum[finite] += c[finite]
        scnt[finite] += 1
        ref, sec = str(p).split("_")
        dt = (pd.Timestamp(sec) - pd.Timestamp(ref)).days
        season = season_by_month[pd.Timestamp(ref).month]
        for z, m in zmasks.items():
            vals = c[m & finite]
            if vals.size:
                rows.append({"zone": z, "pair": str(p), "dt_days": dt,
                             "season": season, "mean_coh": float(vals.mean()),
                             "n_px": int(vals.size)})
        del c, finite
    with np.errstate(invalid="ignore"):
        coh_mean = np.where(scnt > 0, ssum / np.maximum(scnt, 1), np.nan).astype("float32")
    coh_mean_da = xr.DataArray(coh_mean, coords=template.coords, dims=template.dims,
                               name="coh_mean")
    return pd.DataFrame(rows), coh_mean_da


def _decay(dt, g0, ginf, tau):
    return ginf + (g0 - ginf) * np.exp(-dt / tau)


def fit_decay(dt_days, coh) -> dict:
    """Ajuste gamma(t) = ginf + (g0-ginf)*exp(-t/tau). Retourne g0, ginf, tau.

    ginf = plancher de cohérence long-terme ; tau = temps de décorrélation.
    Un tapis flottant instable -> ginf plus bas et/ou tau plus court que la
    même végétation sur sol stable.
    """
    from scipy.optimize import curve_fit

    dt = np.asarray(dt_days, float); y = np.asarray(coh, float)
    ok = np.isfinite(dt) & np.isfinite(y)
    dt, y = dt[ok], y[ok]
    if dt.size < 4:
        return {"g0": np.nan, "ginf": np.nan, "tau": np.nan, "n": int(dt.size)}
    try:
        popt, _ = curve_fit(_decay, dt, y, p0=[0.6, 0.2, 30.0],
                            bounds=([0, 0, 1], [1, 1, 2000]), maxfev=10000)
        g0, ginf, tau = popt
        pred = _decay(dt, *popt)
        ss = 1 - np.sum((y - pred)**2) / max(1e-9, np.sum((y - y.mean())**2))
    except Exception:
        return {"g0": np.nan, "ginf": np.nan, "tau": np.nan, "n": int(dt.size)}
    return {"g0": float(g0), "ginf": float(ginf), "tau": float(tau),
            "r2": float(ss), "n": int(dt.size)}


def paired_zone_diff(df: pd.DataFrame, a: str = "A", b: str = "C",
                     n_boot: int = 2000, seed: int = 0) -> dict:
    """Comparaison APPARIÉE par interférogramme : pour chaque paire présente
    dans A ET C, delta = coh_a - coh_b. Bootstrap sur les paires (gère
    l'autocorrélation temporelle mieux qu'un N brut). Retourne le delta médian,
    son IC95, la fraction de paires où a<b, et un test de signe.
    """
    pa = df[df.zone == a].set_index("pair")["mean_coh"]
    pb = df[df.zone == b].set_index("pair")["mean_coh"]
    common = pa.index.intersection(pb.index)
    if len(common) < 5:
        return {"n_pairs": len(common), "note": "trop peu de paires communes"}
    delta = (pa.loc[common] - pb.loc[common]).values
    rng = np.random.RandomState(seed)
    boots = [rng.choice(delta, len(delta), replace=True).mean() for _ in range(n_boot)]
    return {
        "n_pairs": int(len(common)),
        "delta_median_mm": float(np.median(delta)),   # coh units, name kept generic
        "delta_mean": float(delta.mean()),
        "ci95": [float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))],
        "frac_a_lower": float((delta < 0).mean()),
        "significant": bool(np.percentile(boots, 2.5) > 0 or np.percentile(boots, 97.5) < 0),
    }


# ============================================================================
# Phase D-bis : mécanisme (mécanique/flottaison vs diélectrique) + spatial
# ============================================================================

def pair_hydro_change(pairs: list[str], era5: xr.Dataset, lon: float, lat: float,
                      tau_days: int = 30) -> pd.DataFrame:
    """Par paire : variation absolue du proxy de nappe (|Δ pluie cumulée
    tau_days|) et température minimale des deux dates (K). Base des tests
    couplage-nappe et gel."""
    from .validation import hydrology_proxy

    proxy = hydrology_proxy(era5, lon, lat, tau_days)
    proxy.index = pd.to_datetime(proxy.index)
    t2 = era5["t2m"].sel(latitude=lat, longitude=lon, method="nearest")
    tc = "valid_time" if "valid_time" in t2.coords else "time"
    ts = t2.to_series(); ts.index = pd.to_datetime(t2[tc].values)
    tdaily = ts.resample("1D").mean()

    def at(series, d):
        try:
            return float(series.asof(pd.Timestamp(d)))
        except Exception:
            return np.nan

    rows = []
    for p in pairs:
        r, s = str(p).split("_")
        wr, ws = at(proxy, r), at(proxy, s)
        rows.append({"pair": str(p), "dwtd": abs(ws - wr),
                     "tmin": min(at(tdaily, r), at(tdaily, s))})
    return pd.DataFrame(rows).set_index("pair")


def coherence_vs_hydro(df_perpair: pd.DataFrame, pair_hydro: pd.DataFrame) -> pd.DataFrame:
    """Par zone : régression cohérence ~ |Δ nappe|. Une pente PLUS négative
    pour A que pour C = la cohérence du tapis est plus sensible à la variation
    de nappe -> mécanisme de flottaison (mécanique/hydrologique) propre au tapis."""
    out = []
    for z, g in df_perpair.groupby("zone"):
        m = g.merge(pair_hydro, left_on="pair", right_index=True)
        m = m[np.isfinite(m["dwtd"]) & np.isfinite(m["mean_coh"])]
        if len(m) < 5 or m["dwtd"].std() == 0:
            continue
        slope, _ = np.polyfit(m["dwtd"], m["mean_coh"], 1)
        r = float(np.corrcoef(m["dwtd"], m["mean_coh"])[0, 1])
        out.append({"zone": z, "slope_coh_per_wtd": float(slope), "r": r, "n": int(len(m))})
    return pd.DataFrame(out)


def freeze_coherence_gain(df_perpair: pd.DataFrame, pair_hydro: pd.DataFrame,
                          t_freeze_k: float = 273.15) -> pd.DataFrame:
    """Par zone : cohérence des paires 'froides' (tmin<=0°C, surface figée) vs
    'chaudes'. Si A gagne PLUS que C au gel -> le tapis se stabilise en gelant
    = signature mécanique (flottaison stoppée par le gel)."""
    out = []
    for z, g in df_perpair.groupby("zone"):
        m = g.merge(pair_hydro, left_on="pair", right_index=True)
        cold = m[m["tmin"] <= t_freeze_k]["mean_coh"]
        warm = m[m["tmin"] > t_freeze_k]["mean_coh"]
        if len(cold) >= 3 and len(warm) >= 3:
            out.append({"zone": z, "coh_cold": float(cold.mean()),
                        "coh_warm": float(warm.mean()),
                        "freeze_gain": float(cold.mean() - warm.mean()),
                        "n_cold": int(len(cold)), "n_warm": int(len(warm))})
    return pd.DataFrame(out)


def signed_distance_to_aoi(template: xr.DataArray, cfg: dict) -> xr.DataArray:
    """Distance signée (m) au bord du polygone : négative DEDANS, positive
    DEHORS. Base du profil radial tapis-centre -> bord -> extérieur."""
    from scipy.ndimage import distance_transform_edt
    from .stack import aoi_mask

    aoi = aoi_mask(template, cfg).values
    px = abs(float(template.x[1] - template.x[0]))
    d_out = distance_transform_edt(~aoi) * px
    d_in = distance_transform_edt(aoi) * px
    signed = d_out - d_in
    return xr.DataArray(signed, coords=template.coords, dims=template.dims,
                        name="signed_dist_m")


def radial_profile(field: xr.DataArray, signed_dist: xr.DataArray,
                   edges: np.ndarray) -> pd.DataFrame:
    """Moyenne/médiane de `field` par tranche de distance signée (`edges` en m).
    Montre le gradient centre-tapis -> bord -> extérieur en une courbe."""
    d = signed_dist.values.ravel(); f = field.values.ravel()
    ok = np.isfinite(d) & np.isfinite(f)
    d, f = d[ok], f[ok]
    idx = np.digitize(d, edges)
    rows = []
    for b in range(1, len(edges)):
        m = idx == b
        if m.sum() > 3:
            rows.append({"dist_center_m": 0.5 * (edges[b - 1] + edges[b]),
                         "mean": float(f[m].mean()), "median": float(np.median(f[m])),
                         "n": int(m.sum())})
    return pd.DataFrame(rows)


def zone_field_stats(field: xr.DataArray, zones: dict) -> pd.DataFrame:
    """Stats (médiane, p10/p90, n) d'un champ quelconque par zone — réutilisable
    pour le résidu d'inversion, la variabilité de rétrodiffusion, etc."""
    rows = []
    for z in ("A", "B", "C", "D"):
        v = field.values[zones[z].values & np.isfinite(field.values)]
        if v.size:
            rows.append({"zone": z, "median": float(np.median(v)),
                         "p10": float(np.percentile(v, 10)),
                         "p90": float(np.percentile(v, 90)), "n": int(v.size)})
    return pd.DataFrame(rows)

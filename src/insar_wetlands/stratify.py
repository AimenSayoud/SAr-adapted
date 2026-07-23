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
    """Features de couvert par pixel depuis le stack S2 (green, nir, swir16).

    - greenness_mean : GNDVI moyen = (nir-green)/(nir+green) (biomasse/verdure)
    - greenness_amp  : amplitude saisonnière (p90-p10) de GNDVI (phénologie :
      forêt permanente ~ faible amp, herbe/culture ~ forte amp)
    - wetness_mean   : NDMI moyen = (nir-swir16)/(nir+swir16) (humidité couvert)
    """
    g = s2["green"].astype("float32"); n = s2["nir"].astype("float32")
    s = s2["swir16"].astype("float32")
    gndvi = (n - g) / (n + g)
    ndmi = (n - s) / (n + s)
    q90 = gndvi.quantile(0.9, "time").drop_vars("quantile")
    q10 = gndvi.quantile(0.1, "time").drop_vars("quantile")
    feat = xr.Dataset({
        "greenness_mean": gndvi.mean("time", skipna=True),
        "greenness_amp": q90 - q10,
        "wetness_mean": ndmi.mean("time", skipna=True),
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
        da = rioxarray.open_rasterio(src, masked=False, chunks=True).squeeze("band", drop=True)
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

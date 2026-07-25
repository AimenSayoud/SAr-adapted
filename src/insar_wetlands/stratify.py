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
    # Découper sur l'emprise COMPLÈTE du crop (pas seulement le buffer 500 m) :
    # sinon WorldCover ne couvre qu'un petit carré central et la zone C
    # (prairie extérieure appariée) est confinée à ~500 m + du nodata 'classe 0'
    # partout ailleurs. On prend les bornes du template en degrés.
    try:
        minx, miny, maxx, maxy = template.rio.transform_bounds("EPSG:4326")
    except Exception:
        minx, miny, maxx, maxy = buffered_bbox(cfg)   # repli

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
    """Comparaison APPARIÉE par interférogramme (delta = coh_a - coh_b par
    paire), avec la rigueur exigée en review :

    - `wilcoxon_p` : test de Wilcoxon signé sur les differences appariees (le
      VRAI test ; `frac_a_lower` n'est que la fraction de meme signe, pas un
      test) ;
    - `ci95_pairs` : bootstrap NAÏF sur les paires — OPTIMISTE car les ~356
      paires ne sont pas independantes (elles partagent ~90 dates) ;
    - `date_jackknife` : on retire CHAQUE date (et toutes ses paires) tour a
      tour et on recalcule le delta moyen. Repond directement a « une
      acquisition anormale porte-t-elle le resultat ? » : si le delta reste du
      meme signe pour TOUTE date retiree, il est robuste. C'est le controle qui
      compte face a la non-independance.
    """
    from scipy.stats import wilcoxon

    pa = df[df.zone == a].set_index("pair")["mean_coh"]
    pb = df[df.zone == b].set_index("pair")["mean_coh"]
    common = list(pa.index.intersection(pb.index))
    if len(common) < 5:
        return {"n_pairs": len(common), "note": "trop peu de paires communes"}
    delta = pd.Series((pa.loc[common] - pb.loc[common]).values, index=common)

    # bootstrap naif sur les paires (optimiste)
    rng = np.random.RandomState(seed)
    dv = delta.values
    boots = [rng.choice(dv, len(dv), replace=True).mean() for _ in range(n_boot)]

    # test de Wilcoxon signe (differences appariees)
    try:
        wstat, wp = wilcoxon(dv)
    except Exception:
        wstat, wp = np.nan, np.nan

    # jackknife par DATE : retirer une date -> toutes ses paires
    dates_of = {p: str(p).split("_") for p in common}
    all_dates = sorted({d for p in common for d in dates_of[p]})
    jack = []
    for dd in all_dates:
        keep = [p for p in common if dd not in dates_of[p]]
        if len(keep) >= 5:
            jack.append(float(delta.loc[keep].mean()))
    jack = np.array(jack)

    return {
        "n_pairs": int(len(common)),
        "n_dates": int(len(all_dates)),
        "delta_mean": float(delta.mean()),
        "delta_median": float(np.median(dv)),
        "frac_a_lower": float((dv < 0).mean()),          # descriptif, PAS un test
        "wilcoxon_stat": float(wstat), "wilcoxon_p": float(wp),
        "ci95_pairs": [float(np.percentile(boots, 2.5)),
                       float(np.percentile(boots, 97.5))],  # optimiste
        "date_jackknife_min": float(jack.min()) if jack.size else np.nan,
        "date_jackknife_max": float(jack.max()) if jack.size else np.nan,
        "date_jackknife_se": float(jack.std(ddof=1) * np.sqrt(len(jack) - 1)) if jack.size > 1 else np.nan,
        "robust_same_sign": bool(jack.size and (np.all(jack < 0) or np.all(jack > 0))),
        # 'significant' conservateur : Wilcoxon ET robustesse au jackknife
        "significant": bool(np.isfinite(wp) and wp < 0.05 and jack.size
                            and (np.all(jack < 0) or np.all(jack > 0))),
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


def _pair_indexed(pair_hydro: pd.DataFrame) -> pd.DataFrame:
    """Return `pair_hydro` keyed by `pair`, wherever the key currently lives.

    The key survives a CSV round trip either as the index (written with
    index=True) or as an ordinary column (written after reset_index), and the
    caller cannot always know which. Accepting both is what stops a caching
    detail from silently breaking the merge: without this, a frame whose key
    was dropped merges a date string against a RangeIndex, matches nothing, and
    surfaces much later as an obscure dtype error or an empty result.

    Raises if the key is genuinely absent — that is unrecoverable, and saying so
    plainly beats returning an empty table that looks like a real no-result."""
    if pair_hydro.index.name == "pair":
        return pair_hydro
    if "pair" in pair_hydro.columns:
        return pair_hydro.set_index("pair")
    raise ValueError(
        "pair_hydro must be keyed by `pair`, as index or column; got index "
        f"name {pair_hydro.index.name!r} and columns {list(pair_hydro.columns)}. "
        "A cache written with index=False drops the key — delete it and re-run.")


def coherence_vs_hydro(df_perpair: pd.DataFrame, pair_hydro: pd.DataFrame) -> pd.DataFrame:
    """Par zone : régression cohérence ~ |Δ nappe|. Une pente PLUS négative
    pour A que pour C = la cohérence du tapis est plus sensible à la variation
    de nappe -> mécanisme de flottaison (mécanique/hydrologique) propre au tapis."""
    pair_hydro = _pair_indexed(pair_hydro)
    out = []
    for z, g in df_perpair.groupby("zone"):
        m = g.merge(pair_hydro, left_on="pair", right_index=True)
        m = m[np.isfinite(m["dwtd"]) & np.isfinite(m["mean_coh"])]
        if len(m) < 5 or m["dwtd"].std() == 0:
            continue
        slope, _ = np.polyfit(m["dwtd"], m["mean_coh"], 1)
        r = float(np.corrcoef(m["dwtd"], m["mean_coh"])[0, 1])
        out.append({"zone": z, "slope_coh_per_wtd": float(slope), "r": r, "n": int(len(m))})
    # Guaranteed columns even when no zone qualifies: a column-less DataFrame
    # turns a legitimate "no result" into a KeyError several frames away.
    return pd.DataFrame(out, columns=["zone", "slope_coh_per_wtd", "r", "n"])


def freeze_coherence_gain(df_perpair: pd.DataFrame, pair_hydro: pd.DataFrame,
                          t_freeze_k: float = 273.15) -> pd.DataFrame:
    """Par zone : cohérence des paires 'froides' (tmin<=0°C, surface figée) vs
    'chaudes'. Si A gagne PLUS que C au gel -> le tapis se stabilise en gelant
    = signature mécanique (flottaison stoppée par le gel)."""
    pair_hydro = _pair_indexed(pair_hydro)
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
    return pd.DataFrame(out, columns=["zone", "coh_cold", "coh_warm",
                                      "freeze_gain", "n_cold", "n_warm"])


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


# ============================================================================
# Phase D-ter : mécanisme de diffusion & faisabilité DS-InSAR
# ============================================================================

def amplitude_dispersion_stream(cropped_root, pairs: list[str], template: xr.DataArray):
    """Indice de dispersion d'amplitude D_A = std/mean sur le stack d'amplitude
    (streaming, RAM maîtrisée). D_A BAS = diffuseur stable (candidat PS) ;
    D_A HAUT = pas de cible persistante.

    Répond à H1 (le DS-InSAR/PSI aurait-il la moindre cible ici ?) SANS ISCE :
    si le tapis est uniformément à D_A élevé, aucun PSI ne trouvera de point.
    Retourne la carte D_A (xr.DataArray)."""
    from .stack import load_layer

    tshape = template.shape
    ssum = np.zeros(tshape, "float64"); ssq = np.zeros(tshape, "float64")
    scnt = np.zeros(tshape, "int32")
    for p in pairs:
        try:
            a = load_layer(cropped_root, "amp", [p]).isel(pair=0)
        except Exception:
            continue
        if a.shape != tshape:
            a = a.rio.reproject_match(template)
        v = a.values.astype("float32")
        ok = np.isfinite(v) & (v > 0)
        ssum[ok] += v[ok]; ssq[ok] += v[ok] ** 2; scnt[ok] += 1
        del a, v, ok
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = np.where(scnt > 0, ssum / np.maximum(scnt, 1), np.nan)
        var = np.where(scnt > 1, ssq / np.maximum(scnt, 1) - mean ** 2, np.nan)
        da = np.sqrt(np.maximum(var, 0)) / mean
    return xr.DataArray(da.astype("float32"), coords=template.coords,
                        dims=template.dims, name="amplitude_dispersion")


def backscatter_by_zone(rtc: xr.Dataset, zones: dict) -> pd.DataFrame:
    """Par zone : σ0 VV moyen, sa variabilité temporelle (std), et — si VH
    présent — le ratio VH/VV moyen (discriminant volume vs double-bounce)."""
    rows = []
    vv = rtc["gamma0_vv_db"]
    vv_mean = vv.mean("time", skipna=True); vv_std = vv.std("time", skipna=True)
    ratio = rtc["ratio_vh_vv_db"].mean("time", skipna=True) if "ratio_vh_vv_db" in rtc else None
    for z in ("A", "B", "C", "D"):
        m = zones[z].values
        def stat(field):
            v = field.values[m & np.isfinite(field.values)]
            return float(np.median(v)) if v.size else np.nan
        row = {"zone": z, "sigma0_vv_db": stat(vv_mean),
               "sigma0_vv_temporal_std": stat(vv_std)}
        if ratio is not None:
            row["ratio_vh_vv_db"] = stat(ratio)
        rows.append(row)
    return pd.DataFrame(rows)


def clean_lake_mask(template: xr.DataArray, cfg: dict,
                    worldcover: xr.DataArray | None = None,
                    s2: xr.Dataset | None = None, ndwi_persist: float = 0.2,
                    persist_frac: float = 0.6) -> xr.DataArray:
    """Lac propre = eau PERSISTANTE dans le polygone : WorldCover=eau OU NDWI>seuil
    une fraction persist_frac du temps. Contrôle négatif « vraie eau libre »
    (γ≈0, σ0 spéculaire très bas) — plus propre que flooded_fraction>0.3."""
    from .stack import aoi_mask
    inside = aoi_mask(template, cfg)
    lake = xr.zeros_like(template, dtype=bool)
    if worldcover is not None:
        lake = lake | (worldcover == WC_WATER)
    if s2 is not None and "ndwi" in s2:
        wet_frac = (s2["ndwi"] > ndwi_persist).mean("time")
        lake = lake | (wet_frac.reindex_like(template, method="nearest") > persist_frac
                       if wet_frac.shape != template.shape else wet_frac > persist_frac)
    return (lake & inside).rename("clean_lake")


def s2_phenology_by_zone(s2: xr.Dataset, zones: dict) -> pd.DataFrame:
    """Cycle saisonnier moyen (verdure = -NDWI, humidité = MNDWI ou NDWI) par
    zone et par mois — caractérise l'écologie (fen vs prairie vs lac) et la
    'respiration' (dynamique d'humidité)."""
    green = -s2["ndwi"] if "ndwi" in s2 else None
    wet = s2["mndwi"] if "mndwi" in s2 else s2.get("ndwi")
    tc = "time"
    months = pd.to_datetime(s2[tc].values).month
    rows = []
    for z in ("A", "B", "C", "D"):
        m = zones[z].values
        for var, name in [(green, "greenness"), (wet, "wetness")]:
            if var is None:
                continue
            series = var.where(xr.DataArray(m, coords=zones[z].coords, dims=zones[z].dims)).mean(("y", "x"))
            for mth in range(1, 13):
                sel = series.values[months == mth]
                sel = sel[np.isfinite(sel)]
                if sel.size:
                    rows.append({"zone": z, "var": name, "month": mth,
                                 "value": float(sel.mean())})
    return pd.DataFrame(rows)


def dual_pol_rvi(rtc: xr.Dataset, reduce: str = "median") -> xr.DataArray:
    """Radar Vegetation Index dual-pol : RVI = 4*VH / (VV + VH), en PUISSANCE.

    Quantifie la diffusion de VOLUME de façon NORMALISÉE : ~0 pour une surface
    lisse (VH << VV, réflexion spéculaire / double-bounce), ~1+ pour un volume
    dépolarisant dense (canopée). Plus rigoureux que le simple ratio VH/VV
    utilisé en Phase D-ter : le RVI est insensible a un biais de calibration
    commun aux deux polarisations, car il normalise par la puissance totale.

    C'est le descripteur direct de notre hypothèse dominante (« volume diffusant
    humide »), et un prédicteur de premier plan pour la Phase H.
    """
    vv = 10 ** (rtc["gamma0_vv_db"] / 10.0)      # dB -> puissance linéaire
    vh = 10 ** (rtc["gamma0_vh_db"] / 10.0)
    rvi = (4.0 * vh / (vv + vh)).rename("rvi")
    if reduce and "time" in rvi.dims:
        rvi = getattr(rvi, reduce)("time")
    return rvi


def amplitude_dispersion_from_rtc(rtc: xr.Dataset, pol: str = "vv") -> xr.DataArray:
    """Indice de dispersion d'amplitude D_A calcule depuis la serie RTC σ0 PAR
    DATE (meilleure source que l'amplitude par interferogramme, absente des
    crops) : σ0[dB] -> puissance lineaire -> amplitude = sqrt -> D_A = std/mean
    sur le temps.

    NB : le RTC est multi-vu (~80 m), donc D_A est structurellement plus BAS
    qu'en pleine resolution SLC — le seuil PS 0.25 est indicatif ; c'est la
    comparaison RELATIVE entre zones qui compte (A a-t-il plus de diffuseurs
    stables que le lac ? moins que le sol nu ?)."""
    db = rtc[f"gamma0_{pol}_db"]
    amp = np.sqrt(10.0 ** (db / 10.0))
    mean = amp.mean("time", skipna=True)
    std = amp.std("time", skipna=True)
    with np.errstate(invalid="ignore", divide="ignore"):
        da = (std / mean)
    return da.rename("amplitude_dispersion")

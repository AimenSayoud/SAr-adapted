"""Phase H — PRÉDIRE où Sentinel-1 échoue (et pas seulement constater qu'il échoue).

Passage de « le tapis décorrèle » a « la décorrélation survient sous TELLES
conditions physiques identifiables ». C'est le saut qui transforme un constat
d'échec en modèle prédictif transposable a d'autres tourbières.

Trois apports méthodologiques :

1. **Cible CONTINUE, pas binaire.** Le seuil 0.7 (convention MiaplPy) est
   arbitraire et jette l'essentiel de l'information : réduire 499 pixels a
   « 5.4 % au-dessus de 0.7 » perd toute la structure. On modélise la
   temporal_coherence en CONTINU, et `threshold_sweep` montre la courbe
   complète au lieu d'un seuil unique.

2. **Analyse INTRA-tapis.** Les Phases D-G comparaient A vs C (entre zones).
   Ici on exploite la VARIABILITÉ INTERNE de A : pourquoi certains pixels du
   tapis tiennent-ils et d'autres non ? C'est la question qui a un pouvoir
   prédictif.

3. **Modèles interprétables d'abord.** Corrélations de Spearman (monotones, sans
   hypothèse de forme) + régression multiple standardisée (effets partiels, pour
   distinguer une covariable causale d'une simple corrélée) + R² validé
   croisé. Une forêt aléatoire n'est proposée qu'en complément : pour une thèse,
   un coefficient interprétable vaut mieux qu'une boîte noire.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import xarray as xr


def covariate_table(target: xr.DataArray, covars: dict[str, xr.DataArray],
                    mask: xr.DataArray, min_px: int = 30) -> pd.DataFrame:
    """Table pixel-par-pixel : cible + covariables, restreinte a `mask`.

    Toutes les entrées doivent partager la grille. Les pixels a NaN (sur la
    cible ou une covariable) sont retirés — on garde ainsi un jeu commun a
    toutes les analyses (indispensable pour comparer des modèles)."""
    m = mask.values.astype(bool)
    data = {"target": target.values[m]}
    for k, v in covars.items():
        arr = v.values
        if arr.shape != m.shape:
            raise ValueError(f"covariable '{k}' de forme {arr.shape} != {m.shape}")
        data[k] = arr[m]
    df = pd.DataFrame(data).replace([np.inf, -np.inf], np.nan).dropna()
    if len(df) < min_px:
        raise ValueError(f"trop peu de pixels valides ({len(df)}) — "
                         f"vérifier l'alignement des covariables")
    return df


def rank_covariates(df: pd.DataFrame, target_col: str = "target") -> pd.DataFrame:
    """Corrélation de Spearman de chaque covariable avec la cible.

    Spearman (rang) plutôt que Pearson : capte toute relation MONOTONE sans
    supposer la linéarité, et résiste aux valeurs extrêmes. Retourne rho, p,
    et |rho| trié — la première lecture de « qu'est-ce qui prédit l'échec »."""
    from scipy.stats import spearmanr

    y = df[target_col].values
    rows = []
    for c in df.columns:
        if c == target_col:
            continue
        rho, p = spearmanr(df[c].values, y)
        rows.append({"covariate": c, "spearman_rho": round(float(rho), 4),
                     "p_value": float(p), "abs_rho": abs(float(rho))})
    return (pd.DataFrame(rows).sort_values("abs_rho", ascending=False)
            .drop(columns="abs_rho").reset_index(drop=True))


def collinearity_report(df: pd.DataFrame, target_col: str = "target",
                        vif_max: float = 10.0) -> pd.DataFrame:
    """VIF (facteur d'inflation de la variance) de chaque covariable.

    VIF_j = 1/(1-R²_j), R²_j étant obtenu en régressant la covariable j sur
    TOUTES les autres. VIF > 10 = redondance sévère : le coefficient de la
    variable n'est plus interprétable, il ajuste le bruit dans l'écart entre
    variables quasi identiques (symptôme typique : deux gros coefficients de
    signes opposés sur un couple corrélé).

    Cas réel rencontré ici : `rvi` = 4r/(1+r) et `vh_vv_db` = 10log10(r) sont
    deux transformations MONOTONES du même rapport r = VH/VV — donc rang
    identique, Spearman identique, et colinéarité quasi parfaite."""
    cols = [c for c in df.columns if c != target_col]
    X = df[cols].values.astype(float)
    rows = []
    for j, c in enumerate(cols):
        y = X[:, j]
        others = np.delete(X, j, axis=1)
        A = np.column_stack([np.ones(len(others)), others])
        b, *_ = np.linalg.lstsq(A, y, rcond=None)
        ss_res = float(((y - A @ b) ** 2).sum())
        ss_tot = float(((y - y.mean()) ** 2).sum())
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
        vif = np.inf if r2 >= 1 - 1e-12 else 1.0 / (1.0 - r2)
        rows.append({"covariate": c, "vif": round(float(vif), 2),
                     "redundant": bool(vif > vif_max)})
    return (pd.DataFrame(rows).sort_values("vif", ascending=False)
            .reset_index(drop=True))


def drop_redundant(df: pd.DataFrame, target_col: str = "target",
                   vif_max: float = 10.0, keep: tuple = ()) -> tuple:
    """Retire itérativement la covariable de VIF le plus élevé jusqu'a VIF<=max.

    Itératif et non en une passe : retirer une variable change les VIF des
    autres. `keep` protège les variables a conserver de préférence (on préfère
    garder RVI, normalisé, plutôt que le ratio VH/VV brut).
    Retourne (df_réduit, liste_des_variables_retirées)."""
    out, dropped = df.copy(), []
    while True:
        rep = collinearity_report(out, target_col, vif_max)
        bad = rep[rep.redundant & ~rep.covariate.isin(keep)]
        if bad.empty or len([c for c in out.columns if c != target_col]) <= 2:
            return out, dropped
        worst = bad.iloc[0]["covariate"]
        out = out.drop(columns=[worst])
        dropped.append(worst)


def fit_failure_model(df: pd.DataFrame, target_col: str = "target",
                      n_folds: int = 5, seed: int = 0) -> dict:
    """Régression multiple STANDARDISÉE + R² validé croisé (k-fold).

    Standardiser (z-score) rend les coefficients COMPARABLES entre eux : le
    coefficient est l'effet en écarts-types de cible par écart-type de
    covariable. C'est ce qui distingue une variable a effet PROPRE d'une
    variable simplement corrélée aux autres (contrairement au Spearman
    marginal).

    R² validé croisé : mesure le pouvoir PRÉDICTIF hors échantillon. Un R²
    d'ajustement seul serait optimiste et non publiable.
    Retourne {coefficients, r2_in_sample, r2_cv_mean, r2_cv_std, n}.
    """
    cols = [c for c in df.columns if c != target_col]
    X = df[cols].values.astype(float)
    y = df[target_col].values.astype(float)
    mu, sd = X.mean(0), X.std(0)
    sd[sd == 0] = 1.0
    Xs = (X - mu) / sd
    ys = (y - y.mean()) / (y.std() or 1.0)
    A = np.column_stack([np.ones(len(Xs)), Xs])

    beta, *_ = np.linalg.lstsq(A, ys, rcond=None)
    pred = A @ beta
    r2_in = 1.0 - ((ys - pred) ** 2).sum() / ((ys - ys.mean()) ** 2).sum()

    rng = np.random.default_rng(seed)
    fold = rng.permutation(len(Xs)) % n_folds
    scores = []
    for f in range(n_folds):
        tr, te = fold != f, fold == f
        if te.sum() < 3 or tr.sum() < len(cols) + 2:
            continue
        b, *_ = np.linalg.lstsq(A[tr], ys[tr], rcond=None)
        p = A[te] @ b
        ss_res = ((ys[te] - p) ** 2).sum()
        ss_tot = ((ys[te] - ys[tr].mean()) ** 2).sum()
        if ss_tot > 0:
            scores.append(1.0 - ss_res / ss_tot)
    coefs = (pd.DataFrame({"covariate": cols, "std_coef": beta[1:].round(4)})
             .assign(abs_coef=lambda d: d.std_coef.abs())
             .sort_values("abs_coef", ascending=False)
             .drop(columns="abs_coef").reset_index(drop=True))
    return {"coefficients": coefs, "r2_in_sample": round(float(r2_in), 4),
            "r2_cv_mean": round(float(np.mean(scores)), 4) if scores else np.nan,
            "r2_cv_std": round(float(np.std(scores)), 4) if scores else np.nan,
            "n": int(len(df))}


def spatial_block_cv(df: pd.DataFrame, target_col: str = "target",
                     coords: np.ndarray | None = None,
                     n_blocks: int = 5, seed: int = 0) -> dict:
    """Spatial block cross-validation for multiple regression.

    Instead of random pixel k-fold CV (which leaks autocorrelated spatial
    information between training and validation folds and inflates apparent skill),
    partitions pixels into contiguous spatial blocks larger than the ~160 m
    spatial correlation length.

    Parameters:
        df: DataFrame containing features and target column.
        target_col: Name of target column.
        coords: Optional (N, 2) array of spatial coordinates (e.g. y, x).
                If None, regular contiguous spatial index blocks are used.
        n_blocks: Number of spatial evaluation blocks.
        seed: Random seed for block assignment if spatial clustering is used.

    Returns:
        dict with:
            r2_spatial_cv: Overall out-of-block cross-validated R²
            r2_blocks_mean: Mean of individual block R² scores
            r2_blocks: List of individual block R² scores
            n_blocks: Number of evaluated blocks
            n: Total pixel count
            coefficients: Standardized model coefficients
    """
    cols = [c for c in df.columns if c not in (target_col, "y", "x", "row", "col")]
    X = df[cols].values.astype(float)
    y = df[target_col].values.astype(float)
    mu, sd = X.mean(0), X.std(0)
    sd[sd == 0] = 1.0
    Xs = (X - mu) / sd
    ys = (y - y.mean()) / (y.std() or 1.0)
    A = np.column_stack([np.ones(len(Xs)), Xs])

    beta, *_ = np.linalg.lstsq(A, ys, rcond=None)
    pred_full = A @ beta
    r2_in = 1.0 - ((ys - pred_full) ** 2).sum() / ((ys - ys.mean()) ** 2).sum()

    if coords is not None and len(coords) == len(df):
        cy, cx = coords[:, 0], coords[:, 1]
        grid_dim = int(np.ceil(np.sqrt(n_blocks)))
        y_bins = np.linspace(cy.min(), cy.max() + 1e-6, grid_dim + 1)
        x_bins = np.linspace(cx.min(), cx.max() + 1e-6, grid_dim + 1)
        y_idx = np.digitize(cy, y_bins) - 1
        x_idx = np.digitize(cx, x_bins) - 1
        raw_blocks = y_idx * grid_dim + x_idx
        unique_blocks = np.unique(raw_blocks)
        block_map = {b: i for i, b in enumerate(unique_blocks)}
        block_ids = np.array([block_map[b] for b in raw_blocks])
        actual_n_blocks = len(unique_blocks)
    else:
        chunk_size = int(np.ceil(len(df) / n_blocks))
        block_ids = np.clip(np.arange(len(df)) // chunk_size, 0, n_blocks - 1)
        actual_n_blocks = len(np.unique(block_ids))

    scores = []
    y_pred_cv = np.zeros_like(ys)
    for b_id in np.unique(block_ids):
        te = block_ids == b_id
        tr = ~te
        if te.sum() < 3 or tr.sum() < len(cols) + 2:
            continue
        b, *_ = np.linalg.lstsq(A[tr], ys[tr], rcond=None)
        p = A[te] @ b
        y_pred_cv[te] = p
        ss_res = ((ys[te] - p) ** 2).sum()
        ss_tot = ((ys[te] - ys[tr].mean()) ** 2).sum()
        if ss_tot > 0:
            scores.append(1.0 - ss_res / ss_tot)

    ss_res_tot = ((ys - y_pred_cv) ** 2).sum()
    ss_tot_tot = ((ys - ys.mean()) ** 2).sum()
    r2_spatial_cv = 1.0 - ss_res_tot / ss_tot_tot if ss_tot_tot > 0 else 0.0

    coefs = (pd.DataFrame({"covariate": cols, "std_coef": beta[1:].round(4)})
             .assign(abs_coef=lambda d: d.std_coef.abs())
             .sort_values("abs_coef", ascending=False)
             .drop(columns="abs_coef").reset_index(drop=True))

    return {
        "coefficients": coefs,
        "r2_in_sample": round(float(r2_in), 4),
        "r2_spatial_cv": round(float(r2_spatial_cv), 4),
        "r2_blocks_mean": round(float(np.mean(scores)), 4) if scores else np.nan,
        "r2_blocks": [round(float(s), 4) for s in scores],
        "n_blocks": int(actual_n_blocks),
        "n": int(len(df))
    }


def random_forest_importance(df: pd.DataFrame, target_col: str = "target",
                             seed: int = 0) -> pd.DataFrame | None:
    """Importance de permutation (forêt aléatoire) — COMPLÉMENT non linéaire.

    Capte interactions et seuils que la régression rate. Retourne None si
    scikit-learn est absent (dépendance optionnelle) : la conclusion de la
    Phase H ne doit jamais dépendre d'une boîte noire."""
    try:
        from sklearn.ensemble import RandomForestRegressor
        from sklearn.inspection import permutation_importance
        from sklearn.model_selection import train_test_split
    except ImportError:
        return None
    cols = [c for c in df.columns if c != target_col]
    X, y = df[cols].values, df[target_col].values
    Xtr, Xte, ytr, yte = train_test_split(X, y, test_size=0.3, random_state=seed)
    rf = RandomForestRegressor(n_estimators=300, random_state=seed, n_jobs=-1)
    rf.fit(Xtr, ytr)
    imp = permutation_importance(rf, Xte, yte, n_repeats=10, random_state=seed)
    return (pd.DataFrame({"covariate": cols,
                          "perm_importance": imp.importances_mean.round(4),
                          "std": imp.importances_std.round(4)})
            .sort_values("perm_importance", ascending=False)
            .assign(r2_test=round(float(rf.score(Xte, yte)), 4))
            .reset_index(drop=True))


def zone_values(field: xr.DataArray, zones: dict, zone: str) -> np.ndarray:
    """Finite pixel values of `field` inside `zone`.

    The single place invalid pixels are dropped, so every per-zone statistic
    shares one denominator. Written out because the obvious inline alternative
    is silently wrong::

        np.nanmean(field.values[mask] >= 0.7)   # WRONG

    ``NaN >= 0.7`` evaluates to ``False``, not ``NaN``, so the boolean array
    holds no NaN for ``nanmean`` to skip: invalid pixels are counted as
    failures and the fraction is diluted by whatever share of the zone is
    masked. Zones without invalid pixels agree with the correct value, so the
    error hides until a fragmented or edge-clipped zone is compared against
    one that is fully covered."""
    v = field.values[zones[zone].values]
    return v[np.isfinite(v)]


def zone_fraction_above(field: xr.DataArray, zones: dict, zone: str,
                        threshold: float = 0.7) -> float:
    """Fraction of VALID pixels in `zone` at or above `threshold`."""
    v = zone_values(field, zones, zone)
    return float((v >= threshold).mean()) if v.size else float("nan")


def threshold_sweep(field: xr.DataArray, zones: dict,
                    thresholds=np.arange(0.4, 0.96, 0.05),
                    zone_names=("A", "B", "C", "D")) -> pd.DataFrame:
    """Fraction de pixels au-dessus de CHAQUE seuil, par zone (multi-seuil).

    Remplace le verdict binaire « % >= 0.7 » par la COURBE complète : on voit si
    les zones diffèrent par un décalage global ou seulement dans une queue de
    distribution — information perdue par un seuil unique. Le croisement des
    courbes (s'il existe) identifie le régime ou les zones se ressemblent."""
    rows = []
    for z in zone_names:
        if z not in zones:
            continue
        v = zone_values(field, zones, z)
        if not v.size:
            continue
        for t in thresholds:
            rows.append({"zone": z, "threshold": round(float(t), 3),
                         "frac_above": float((v >= t).mean()),  # rounded once, at export (C-039)
                         "n_px": int(v.size)})
    return pd.DataFrame(rows)

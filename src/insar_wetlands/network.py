"""Phase 3 — Controle qualite du reseau interferometrique et connectivite."""

from __future__ import annotations

import numpy as np
import pandas as pd
import xarray as xr


def pair_stats(corr_stack: xr.DataArray, aoi: xr.DataArray) -> pd.DataFrame:
    """Coherence moyenne/mediane par paire, sur l'AOI et hors AOI."""
    inside = corr_stack.where(aoi)
    outside = corr_stack.where(~aoi)
    # NB: (inside > 0.3).mean() serait faux — la comparaison transforme les
    # NaN hors-AOI en False, donc la moyenne se ferait sur TOUTE la grille
    # croppee (AOI ~3% de la grille) au lieu des seuls pixels AOI.
    n_aoi = float(aoi.sum())
    df = pd.DataFrame({
        "pair": corr_stack.pair.values,
        "ref_date": pd.to_datetime(corr_stack.ref_date.values),
        "sec_date": pd.to_datetime(corr_stack.sec_date.values),
        "coh_aoi_mean": inside.mean(("y", "x")).values,
        "coh_aoi_median": inside.median(("y", "x")).values,
        "coh_out_mean": outside.mean(("y", "x")).values,
        "frac_aoi_coh_gt_0p3": (inside > 0.3).sum(("y", "x")).values / n_aoi,
    })
    df["dt_days"] = (df.sec_date - df.ref_date).dt.days
    df["season"] = df.ref_date.dt.month.map(
        lambda m: "DJF" if m in (12, 1, 2) else
                  "MAM" if m in (3, 4, 5) else
                  "JJA" if m in (6, 7, 8) else "SON")
    return df


def select_pairs(stats: pd.DataFrame, min_coherence: float = 0.30) -> pd.DataFrame:
    """Marque les paires conservees (keep=True) selon la coherence AOI."""
    out = stats.copy()
    out["keep"] = out["coh_aoi_mean"] >= min_coherence
    return out


def connectivity(stats: pd.DataFrame, keep_col: str = "keep"):
    """Composantes connexes du graphe temporel avec les paires conservees.

    Retourne (n_components, labels_par_date, dates). n_components > 1 =>
    l'inversion SBAS sera singuliere : il faut des paires-ponts.
    """
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components

    kept = stats[stats[keep_col]]
    dates = sorted(set(stats.ref_date) | set(stats.sec_date))
    idx = {d: i for i, d in enumerate(dates)}
    i = [idx[d] for d in kept.ref_date]
    j = [idx[d] for d in kept.sec_date]
    n = len(dates)
    adj = coo_matrix((np.ones(len(i)), (i, j)), shape=(n, n))
    n_comp, labels = connected_components(adj, directed=False)
    return n_comp, labels, pd.DatetimeIndex(dates)


def component_report(stats: pd.DataFrame, keep_col: str = "keep") -> dict:
    """Diagnostic de connectivite exploitable pour le go/no-go SBAS.

    Distingue deux situations tres differentes derriere un meme n_comp :
      - des dates ORPHELINES (toutes leurs paires rejetees) : MintPy les
        ignore simplement, ce n'est PAS bloquant ;
      - une vraie SCISSION du reseau en gros blocs (ex: un ete decorrele qui
        coupe annee N / annee N+1) : la, il faut des ponts.
    """
    n_comp, labels, dates = connectivity(stats, keep_col)
    sizes = pd.Series(labels).value_counts().sort_values(ascending=False)
    comp_dates = {int(c): sorted(dates[labels == c]) for c in sizes.index}
    orphans = [d for c in sizes.index if sizes[c] == 1
               for d in comp_dates[int(c)]]
    big = sizes[sizes > 1]
    largest = int(sizes.iloc[0])
    verdict = ("connexe" if n_comp == 1
               else "orphelins_seulement" if len(big) <= 1
               else "scission_reelle")
    return {
        "n_components": int(n_comp),
        "n_dates": len(dates),
        "largest_component": largest,
        "coverage_pct": round(100 * largest / len(dates), 1),
        "n_orphan_dates": len(orphans),
        "orphan_dates": [d.date().isoformat() for d in orphans],
        "big_component_sizes": [int(s) for s in big.tolist()],
        "verdict": verdict,
    }


def suggest_bridges(stats: pd.DataFrame, max_bridge_days: int = 120,
                    per_boundary: int = 3) -> pd.DataFrame:
    """Paires-ponts NOUVELLES a soumettre pour reconnecter les composantes.

    Exclut les paires deja presentes dans le reseau (rejetees pour coherence
    basse : les resoumettre a HyP3 donnerait le meme resultat). Pour chaque
    frontiere entre composantes, propose les paires inter-composantes les
    plus courtes qui n'existent pas encore.
    """
    n_comp, labels, dates = connectivity(stats)
    cols = ["ref_date", "sec_date", "dt_days", "pair"]
    if n_comp == 1:
        return pd.DataFrame(columns=cols)
    existing = set(stats["pair"])
    rows = []
    for ca in range(n_comp):
        for cb in range(ca + 1, n_comp):
            da = dates[labels == ca]
            db = dates[labels == cb]
            cands = sorted({(a, b) if a < b else (b, a)
                            for a in da for b in db
                            if abs((b - a).days) <= max_bridge_days},
                           key=lambda p: (p[1] - p[0]).days)
            kept = 0
            for a, b in cands:
                pair = f"{a:%Y%m%d}_{b:%Y%m%d}"
                if pair in existing:
                    continue
                rows.append({"ref_date": a, "sec_date": b,
                             "dt_days": (b - a).days, "pair": pair})
                kept += 1
                if kept >= per_boundary:
                    break
    return (pd.DataFrame(rows, columns=cols)
            .drop_duplicates("pair")
            .sort_values("dt_days")
            .reset_index(drop=True))


def design_matrix(pairs: list[str], dates: pd.DatetimeIndex) -> np.ndarray:
    """Matrice de design SBAS A (n_pairs x n_dates-1), phase incrementale.

    phi_pair = somme des increments entre ref_date et sec_date.
    """
    idx = {d: i for i, d in enumerate(dates)}
    A = np.zeros((len(pairs), len(dates) - 1))
    for k, p in enumerate(pairs):
        a, b = p.split("_")
        i, j = idx[pd.Timestamp(a)], idx[pd.Timestamp(b)]
        A[k, i:j] = 1.0
    return A

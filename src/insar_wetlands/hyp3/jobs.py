"""Phase 2 — Construction du reseau SBAS et soumission des jobs HyP3 burst InSAR."""

from __future__ import annotations

import pandas as pd


def build_sbas_pairs(dates: pd.Series, max_temporal_days: int = 48,
                     max_pairs_per_date: int = 4) -> pd.DataFrame:
    """Reseau small-baseline : toutes les paires <= max_temporal_days.

    max_pairs_per_date limite le nombre de connexions vers l'avant par date
    (les 6/12/24/36-48 jours les plus proches) pour maitriser le nombre de jobs.
    NB: la baseline spatiale des bursts S1 est toujours < 150 m (orbite tube),
    donc seul le critere temporel filtre ici.
    """
    dates = sorted(pd.to_datetime(dates).drop_duplicates())
    rows = []
    for i, d1 in enumerate(dates):
        forward = 0
        for d2 in dates[i + 1:]:
            dt = (d2 - d1).days
            if dt > max_temporal_days or forward >= max_pairs_per_date:
                break
            rows.append({"ref_date": d1, "sec_date": d2, "dt_days": dt})
            forward += 1
    df = pd.DataFrame(rows)
    df["pair"] = (df.ref_date.dt.strftime("%Y%m%d") + "_"
                  + df.sec_date.dt.strftime("%Y%m%d"))
    return df


def select_test_batch(pairs: pd.DataFrame, n: int = 12) -> pd.DataFrame:
    """Lot test go/no-go : la paire la plus courte de chaque trimestre.

    Permet de verifier la coherence C-band sur la tourbiere dans toutes les
    saisons AVANT de soumettre les centaines de paires du reseau complet.
    """
    p = pairs.copy()
    p["quarter"] = p.ref_date.dt.to_period("Q")
    test = (p.sort_values("dt_days")
             .groupby("quarter", as_index=False)
             .first()
             .sort_values("ref_date"))
    if len(test) > n:
        import numpy as np

        test = test.iloc[np.linspace(0, len(test) - 1, n).round().astype(int)]
    return test.drop(columns="quarter").reset_index(drop=True)


def date_to_granule(inventory: pd.DataFrame) -> dict:
    """Mapping date -> nom de granule burst (depuis l'inventaire Phase 1)."""
    inv = inventory.copy()
    inv["date"] = pd.to_datetime(inv["date"])
    return dict(zip(inv["date"], inv["granule"]))


def submit_pairs(pairs: pd.DataFrame, granules: dict, name: str,
                 looks: str = "10x2") -> pd.DataFrame:
    """Soumet les paires en jobs INSAR_ISCE_BURST. Retourne le suivi (job_id)."""
    import hyp3_sdk as sdk

    hyp3 = sdk.HyP3()  # lit ~/.netrc (setup_earthdata)
    rows = []
    for _, p in pairs.iterrows():
        ref_g = granules.get(pd.Timestamp(p.ref_date))
        sec_g = granules.get(pd.Timestamp(p.sec_date))
        if not ref_g or not sec_g:
            rows.append({"pair": p.pair, "job_id": None,
                         "status": "MISSING_GRANULE"})
            continue
        job = hyp3.submit_insar_isce_burst_job(
            granule1=ref_g, granule2=sec_g, name=name, looks=looks,
            apply_water_mask=False,
        )
        j = job.jobs[0] if hasattr(job, "jobs") else job
        rows.append({"pair": p.pair, "job_id": j.job_id, "status": "PENDING"})
    return pairs.merge(pd.DataFrame(rows), on="pair")


def fetch_jobs(name: str) -> "object":
    """Recupere le batch de jobs HyP3 portant ce nom (suivi/telechargement)."""
    import hyp3_sdk as sdk

    hyp3 = sdk.HyP3()
    return hyp3.find_jobs(name=name)


def check_credits() -> float:
    import hyp3_sdk as sdk

    return sdk.HyP3().check_credits()

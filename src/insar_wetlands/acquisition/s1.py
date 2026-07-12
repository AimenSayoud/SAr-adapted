"""Inventaire Sentinel-1 : bursts SLC couvrant l'AOI, via asf_search.

Strategie HyP3 burst InSAR : on ne telecharge PAS les SLC. On construit le
catalogue (trace, burst, dates) qui servira a soumettre les paires en Phase 2.
"""

from __future__ import annotations

import pandas as pd


def search_s1_bursts(aoi_wkt: str, start: str, end: str,
                     polarization: str = "VV") -> pd.DataFrame:
    """Recherche tous les bursts S1 SLC-BURST intersectant l'AOI.

    Retourne un DataFrame : une ligne par acquisition de burst, avec
    trace relative, direction, burst ID, date.
    """
    import asf_search as asf

    results = asf.search(
        platform=[asf.PLATFORM.SENTINEL1],
        processingLevel=asf.PRODUCT_TYPE.BURST,
        polarization=polarization,
        intersectsWith=aoi_wkt,
        start=start,
        end=end,
    )
    rows = []
    for r in results:
        p = r.properties
        burst = p.get("burst", {}) or {}
        rows.append({
            "granule": p.get("fileID") or p.get("sceneName"),
            "date": pd.to_datetime(p["startTime"]).tz_localize(None).normalize(),
            "datetime": pd.to_datetime(p["startTime"]).tz_localize(None),
            "flight_direction": p.get("flightDirection"),
            "relative_orbit": p.get("pathNumber"),
            "full_burst_id": burst.get("fullBurstID"),
            "subswath": burst.get("subswath"),
            "polarization": p.get("polarization"),
            "platform": p.get("platform"),
        })
    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values("datetime").reset_index(drop=True)
    return df


def summarize_tracks(df: pd.DataFrame) -> pd.DataFrame:
    """Resume par (direction, trace, burst) : nb d'acquisitions, intervalle median.

    Sert a choisir LA trace/burst a figer dans config.yaml.
    """
    groups = []
    for (fd, orbit, burst), g in df.groupby(
            ["flight_direction", "relative_orbit", "full_burst_id"]):
        dates = g["date"].drop_duplicates().sort_values()
        gaps = dates.diff().dt.days.dropna()
        groups.append({
            "flight_direction": fd,
            "relative_orbit": orbit,
            "full_burst_id": burst,
            "n_acquisitions": len(dates),
            "first": dates.iloc[0].date(),
            "last": dates.iloc[-1].date(),
            "median_gap_days": gaps.median() if len(gaps) else None,
            "max_gap_days": gaps.max() if len(gaps) else None,
            "n_gaps_gt_24d": int((gaps > 24).sum()) if len(gaps) else None,
        })
    out = pd.DataFrame(groups).sort_values("n_acquisitions", ascending=False)
    return out.reset_index(drop=True)

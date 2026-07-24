"""Inventaire Sentinel-2 L2A via STAC (Earth-Search AWS, sans cle API)."""

from __future__ import annotations

import pandas as pd


def search_s2(cfg: dict, bbox: tuple[float, float, float, float]) -> pd.DataFrame:
    """Inventaire des scenes S2 L2A sous le seuil de nuages, sur l'AOI."""
    from pystac_client import Client

    s2 = cfg["sentinel2"]
    client = Client.open(s2["stac_url"])
    search = client.search(
        collections=[s2["collection"]],
        bbox=list(bbox),
        datetime=f"{cfg['time']['start']}/{cfg['time']['end']}",
        query={"eo:cloud_cover": {"lt": s2["max_cloud_cover"]}},
    )
    rows = []
    for item in search.items():
        rows.append({
            "item_id": item.id,
            "datetime": pd.to_datetime(item.datetime).tz_localize(None),
            "date": pd.to_datetime(item.datetime).tz_localize(None).normalize(),
            "cloud_cover": item.properties.get("eo:cloud_cover"),
            "mgrs_tile": item.properties.get("grid:code")
                         or item.properties.get("s2:mgrs_tile"),
            "platform": item.properties.get("platform"),
        })
    df = pd.DataFrame(rows)
    if not df.empty:
        # Une meme date peut couvrir 2 tuiles MGRS -> on garde la moins nuageuse
        df = (df.sort_values(["date", "cloud_cover"])
                .drop_duplicates("date")
                .sort_values("date")
                .reset_index(drop=True))
    return df

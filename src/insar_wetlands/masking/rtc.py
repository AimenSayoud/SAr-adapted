"""Phase 5 (support) — Retrodiffusion sigma0/gamma0 VV par date S1.

Source : collection 'sentinel-1-rtc' du Microsoft Planetary Computer
(RTC deja calcule, gratuit) — evite de payer des jobs RTC HyP3.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

MPC_STAC = "https://planetarycomputer.microsoft.com/api/stac/v1"


def search_rtc(cfg: dict, bbox, relative_orbit: int | None = None) -> list:
    import planetary_computer as pc
    from pystac_client import Client

    client = Client.open(MPC_STAC, modifier=pc.sign_inplace)
    query = {"sat:relative_orbit": {"eq": relative_orbit}} if relative_orbit else None
    search = client.search(
        collections=["sentinel-1-rtc"],
        bbox=list(bbox),
        datetime=f"{cfg['time']['start']}/{cfg['time']['end']}",
        query=query,
    )
    return list(search.items())


def build_rtc_stack(items: list, template: xr.DataArray,
                    out_nc: str | Path, checkpoint_every: int = 20) -> Path:
    """Stack gamma0 VV (dB) reprojete nearest sur la grille HyP3.

    Idempotent, dedup par date, checkpoint regulier (memes protections que
    build_s2_stack : coupures Colab, doublons de frames adjacentes).
    """
    import planetary_computer as pc
    import rioxarray
    from rasterio.enums import Resampling

    from .s2_fusion import _flush

    out_nc = Path(out_nc)
    out_nc.parent.mkdir(parents=True, exist_ok=True)
    existing = xr.load_dataset(out_nc) if out_nc.exists() else None
    done = (set(pd.to_datetime(existing.time.values))
            if existing is not None else set())
    bounds = template.rio.transform_bounds("EPSG:4326")
    new = []
    for item in sorted(items, key=lambda it: pd.to_datetime(it.datetime)):
        t = pd.to_datetime(item.datetime).tz_localize(None).normalize()
        if t in done or "vv" not in item.assets:
            continue
        try:
            href = pc.sign(item.assets["vv"].href)
            da = rioxarray.open_rasterio(href, masked=True).squeeze("band", drop=True)
            da = da.rio.clip_box(*bounds, crs="EPSG:4326")
            da = da.rio.reproject_match(template, resampling=Resampling.nearest)
            db = 10 * np.log10(da.where(da > 0))
            new.append(db.rename("gamma0_vv_db").expand_dims(time=[t]).to_dataset())
            done.add(t)
            print(f"  + {t.date()}")
        except Exception as e:
            print(f"  ! {item.id}: {e}")
        if len(new) >= checkpoint_every:
            existing = _flush(existing, new, out_nc)
            new = []
            print(f"  ... checkpoint ({existing.time.size} dates sur disque)")
    if new:
        _flush(existing, new, out_nc)
    return out_nc

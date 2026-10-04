"""The data cube (C-053…C-057, chain data-cube): ingest → harmonise → gold → check, every step resumable.

    INSAR_DRIVE_ROOT=../local/drive_pristine PYTHONPATH=src python scripts/cube_build.py ingest s2
    ... ingest {s2,landsat,worldcover,meteo,imgw,lidar,nisar} | harmonise {ifg,rtc,optical,static,site,folds}
    ... gold {plots,mat} | check | all

Mac-safe: reads crops already on disk, windows of public rasters and small tables; nothing here inverts a network.
Field values are read from the hub and written only to the hub (``06_data/cube``).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from insar_wetlands.cube import cube_root
from insar_wetlands.cube import sources as src

HUB = Path(__file__).resolve().parents[3]
CUBE = cube_root()
RAW, SILVER, GOLD = CUBE / "raw", CUBE / "silver", CUBE / "gold"
START, END = "2020-01-01", f"{pd.Timestamp.now(tz='UTC'):%Y-%m-%d}"
SITE_LATLON = (52.761413, 16.3099)        # P6/CR, from the field team's plot table (readme sheet)


def log(msg: str) -> None:
    print(f"[{pd.Timestamp.now():%H:%M:%S}] {msg}", flush=True)


_CTX = None


def ctx():
    global _CTX
    if _CTX is None:
        from insar_wetlands.bootstrap import start
        _CTX = start("cube_build", mount=False, git=False)
    return _CTX


def template_xy():
    t = ctx().template
    return t.x.values.astype(float), t.y.values.astype(float)


def peat_bounds_utm(buffer_m: float = 300.0) -> tuple[float, float, float, float]:
    """Zones A–C plus a buffer: where fine data (LiDAR) are worth fetching."""
    from insar_wetlands.cube.build import peat_bounds_utm as pb
    return pb(ctx(), buffer_m)


# ------------------------------------------------------------------------------ M1 ingest

def ingest(what: str) -> None:
    tx, ty = template_xy()
    bounds = src.template_bounds(tx, ty)
    if what == "s2":
        tab = src.fetch_s2(RAW / "s2", bounds, START, END, log=log)
        log(f"S2 scenes {len(tab)}, with valid pixels {int((tab.valid_share > 0).sum())}")
    elif what == "landsat":
        tab = src.fetch_landsat_lst(RAW / "landsat", bounds, START, END, log=log)
        log(f"Landsat scenes {len(tab)}, with clear pixels {int((tab.clear_share > 0).sum())}")
    elif what == "worldcover":
        src.fetch_worldcover(RAW / "worldcover", bounds, log=log)
    elif what == "meteo":
        df = src.fetch_open_meteo(RAW / "open_meteo", *SITE_LATLON, START, END, log=log)
        df.to_csv(RAW / "open_meteo" / "hourly.csv.gz")
        log(f"Open-Meteo hourly rows {len(df)}")
    elif what == "imgw":
        df = src.fetch_imgw(RAW / "imgw", int(START[:4]), int(END[:4]), log=log)
        df.to_csv(RAW / "imgw" / "daily.csv", index=False)
        log(f"IMGW station-days {len(df)}")
    elif what == "lidar":
        from pyproj import Transformer
        b = peat_bounds_utm()
        t = Transformer.from_crs(32633, 2180, always_xy=True)
        xs, ys = t.transform([b[0], b[2], b[0], b[2]], [b[1], b[1], b[3], b[3]])
        b2180 = (np.floor(min(xs) / 500) * 500, np.floor(min(ys) / 500) * 500, np.ceil(max(xs) / 500) * 500,
                 np.ceil(max(ys) / 500) * 500)
        tab = src.fetch_gugik_tiles(RAW / "lidar", b2180, log=log)
        tab.to_csv(RAW / "lidar" / "tiles.csv", index=False)
    elif what == "ecostress":
        tab = src.fetch_ecostress(RAW / "ecostress", bounds, START, END, log=log)
        log(f"ECOSTRESS acquisitions {len(tab)}, with clear pixels {int((tab.valid_share > 0).sum())}")
    elif what == "nisar":
        from insar_wetlands.cube import nisar_crops
        nisar_crops.fetch(RAW / "nisar", ctx(), log=log)
    else:
        raise SystemExit(f"unknown source {what}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("step", choices=["ingest", "harmonise", "gold", "check", "all"])
    ap.add_argument("what", nargs="*")
    a = ap.parse_args(argv)
    for d in (RAW, SILVER, GOLD):
        d.mkdir(parents=True, exist_ok=True)
    if a.step == "ingest":
        for w in a.what:
            ingest(w)
    else:
        from insar_wetlands.cube import build
        build.run(a.step, a.what, ctx, RAW, SILVER, GOLD, log=log)


if __name__ == "__main__":
    sys.exit(main())

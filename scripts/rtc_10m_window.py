"""Sentinel-1 RTC γ⁰ at its native 10 m over the site, summer 2022–2024, both tracks (X-060, X-064).

The analysis grid is 40 m; the RTC product is 10 m. Two questions need the finer grid: whether
a narrow linear structure (a boardwalk) shows in backscatter, and which 40 m cells are pure
open water. This writes a 10 m grid nested in the 40 m one (each 40 m cell = 4×4 cells of
10 m, same origin) for the snow-free months only (May–September), every acquisition, using the
same per-day slice mosaicking as the 40 m stacks (``build_rtc_track_stack``, nearest).

    PYTHONPATH=src python scripts/rtc_10m_window.py --track ascending
    PYTHONPATH=src python scripts/rtc_10m_window.py --track descending

Idempotent: rerun to retry failed days. Output: ``<out-root>/rtc10m_summer_2022_2024_<track>.nc``.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from insar_wetlands.bootstrap import start
from insar_wetlands.masking import rtc

MONTHS = range(5, 10)


def nested_grid(template: xr.DataArray, factor: int = 4) -> xr.DataArray:
    """A grid ``factor`` times finer than ``template``, each coarse cell split into factor×factor
    cells with the same outer edges (cell centres offset by ±(k+½)·res/factor from the coarse edge)."""
    import rioxarray  # noqa: F401

    dx = float(template.x[1] - template.x[0])
    dy = float(template.y[1] - template.y[0])
    off = (np.arange(factor) + 0.5) / factor - 0.5
    x = (template.x.values[:, None] + off[None, :] * dx).ravel()
    y = (template.y.values[:, None] + off[None, :] * dy).ravel()
    da = xr.DataArray(np.zeros((y.size, x.size), "float32"), coords={"y": y, "x": x}, dims=("y", "x"))
    return da.rio.write_crs(template.rio.crs).rio.write_transform(
        rtc_transform(x, y))


def rtc_transform(x: np.ndarray, y: np.ndarray):
    from affine import Affine

    rx, ry = x[1] - x[0], y[1] - y[0]
    return Affine(rx, 0, x[0] - rx / 2, 0, ry, y[0] - ry / 2)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--track", choices=["ascending", "descending"], required=True)
    ap.add_argument("--out-root", default="../local/drive_mirror")
    a = ap.parse_args(argv)
    if "drive_pristine" in str(Path(a.out_root).resolve()):
        raise SystemExit("refusing to write into the read-only snapshot")
    ctx = start("rtc_10m_window", mount=False, git=False)
    tpl10 = nested_grid(ctx.template)
    bbox = ctx.template.rio.transform_bounds("EPSG:4326")
    orbit = rtc.track_relative_orbit(ctx.cfg, a.track)
    items = rtc.search_rtc(ctx.cfg, rtc.widen_bbox(bbox), orbit, start="2022-05-01", end="2024-09-30")
    items = [it for it in items if rtc.item_day(it).month in MONTHS]
    print(f"{a.track}: {len({rtc.item_day(it) for it in items})} summer days, grid {tpl10.shape} at 10 m")
    out = Path(a.out_root) / f"rtc10m_summer_2022_2024_{a.track}.nc"
    rtc.build_rtc_track_stack(items, tpl10, out, track=a.track, relative_orbit=orbit, bbox=bbox,
                              resampling="nearest", checkpoint_every=10)
    ds = xr.open_dataset(out)
    print(f"wrote {out}: {ds.time.size} dates, {pd.Timestamp(ds.time.values[0]).date()} → {pd.Timestamp(ds.time.values[-1]).date()}")


if __name__ == "__main__":
    main()

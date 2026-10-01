"""Sentinel-1 RTC backscatter (γ⁰ VV, VH, VH−VV) for one track, on the HyP3 crop grid.

The ascending stack (``rtc_dualpol_stack.nc``, 2022–2024, built in phaseDter) is the only
backscatter the project had, so the descending track had coherence but no VV/VH/RVI. This script
builds the same stack for any track and period from the Microsoft Planetary Computer collection
``sentinel-1-rtc`` (free, no HyP3 credits): VV and VH read over the AOI only, put on the track's
interferogram grid, converted to dB.

- Track geometry comes from ``config.yaml`` (``sentinel1.tracks``); the grid is the track's first
  cropped coherence layer; the output is ``rtc_dualpol_<y0>_<y1>.nc`` suffixed by track
  (``Paths.track_file``: ``rtc_dualpol_2020_2024.nc``, ``rtc_dualpol_2020_2024_descending.nc``) — never
  the 2022–2024 ascending ``rtc_dualpol_stack.nc``, so nothing that reads it changes.
- Resampling defaults to ``nearest``, as the ascending stack was built, so the tracks compare.
- A day is a mosaic of its slices (best first, the others only where it has no data).
- Idempotent: rerun to resume or to retry a failed day. A coverage table says, date by date,
  whether the RTC catalogue has it, the stack has it and an interferogram uses it.

On this Mac — read the crops from the read-only snapshot, write to the local mirror (copying to
Drive is a separate step):

    INSAR_DRIVE_ROOT=../local/drive_pristine PYTHONPATH=src python scripts/rtc_download.py \\
        --track descending --out-root ../local/drive_mirror \\
        --extra-cropped ../local/s1_2020_2021/hyp3_cropped_descending
    ... --check          # coverage table only, downloads nothing

On Colab (Drive mounted), ``--out-root`` defaults to the data root: the stack lands beside the crops.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import xarray as xr

from insar_wetlands.aoi import buffered_bbox
from insar_wetlands.config import load_config
from insar_wetlands.masking.rtc import (
    aoi_cover,
    build_rtc_track_stack,
    items_by_day,
    rtc_coverage,
    search_rtc,
    track_relative_orbit,
    widen_bbox,
)
from insar_wetlands.paths import make_paths
from insar_wetlands.stack import list_pairs, load_layer


def stack_names(start: str, end: str) -> tuple[str, str]:
    """``rtc_dualpol_<y0>_<y1>.nc`` and its coverage table. Never the name of the 2022–2024
    ascending stack (``rtc_dualpol_stack.nc``), whose consumers must not change under them."""
    period = f"{pd.Timestamp(start).year}_{pd.Timestamp(end).year}"
    return f"rtc_dualpol_{period}.nc", f"rtc_coverage_{period}.csv"


def interferogram_days(roots: list[Path]) -> set[pd.Timestamp]:
    """Every date used by a cropped pair (``YYYYMMDD_YYYYMMDD`` directories with ``ok``)."""
    days = set()
    for root in roots:
        for pair in list_pairs(root):
            days |= set(pd.to_datetime(pair.split("_")[:2]))
    return days


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--track", choices=("ascending", "descending"), default="descending")
    ap.add_argument("--start", default="2020-01-01")
    ap.add_argument("--end", default="2024-12-31")
    ap.add_argument("--out-root", type=Path, help="where the stack goes (default: the data root)")
    ap.add_argument("--extra-cropped", type=Path, nargs="*", default=[],
                    help="more crop folders whose dates count as interferogram dates (e.g. 2020–2021)")
    ap.add_argument("--resampling", choices=("nearest", "average"), default="nearest")
    ap.add_argument("--check", action="store_true", help="coverage table only, no download")
    args = ap.parse_args()

    cfg = load_config()
    src = make_paths(cfg=cfg, track=args.track)
    out = make_paths(cfg=cfg, root=args.out_root or src.drive, track=args.track)
    if "drive_pristine" in out.drive.resolve().parts:
        raise SystemExit(f"{out.drive} is the read-only snapshot: pass --out-root")
    stack_nc, coverage_csv = (out.track_file(n) for n in stack_names(args.start, args.end))

    orbit = track_relative_orbit(cfg, args.track)
    bbox = buffered_bbox(cfg)
    # Search wider than the AOI: the slice holding the rest of a pass can have a footprint that
    # misses it. Keep the days on which at least one slice reaches the AOI, with all their slices.
    found = search_rtc(cfg, widen_bbox(bbox), relative_orbit=orbit, start=args.start, end=args.end)
    days = {d: v for d, v in items_by_day(found, bbox).items() if aoi_cover(v[0], bbox) > 0}
    items = [it for v in days.values() for it in v]
    print(f"{args.track} (relative orbit {orbit}), {args.start} → {args.end}: "
          f"{len(items)} RTC items, {len(days)} days")

    if not args.check:
        template = load_layer(src.cropped, "corr", list_pairs(src.cropped)[:1]).isel(pair=0)
        build_rtc_track_stack(items, template, stack_nc, track=args.track, relative_orbit=orbit,
                              bbox=bbox, resampling=args.resampling)

    on_disk = pd.to_datetime(xr.open_dataset(stack_nc).time.values) if stack_nc.exists() else []
    ifg = {d for d in interferogram_days([src.cropped, *args.extra_cropped])
           if pd.Timestamp(args.start) <= d <= pd.Timestamp(args.end)}
    cov = rtc_coverage(list(days), on_disk, ifg)
    cov.insert(1, "source_candidates", [",".join(it.id for it in days.get(d, [])) for d in cov.date])
    coverage_csv.parent.mkdir(parents=True, exist_ok=True)
    cov.to_csv(coverage_csv, index=False)

    gaps = cov[cov.in_interferograms & ~cov.in_stack]
    print(f"stack: {len(on_disk)} dates · catalogue: {int(cov.in_catalogue.sum())} · "
          f"interferogram dates: {int(cov.in_interferograms.sum())}, without RTC: {len(gaps)}")
    not_in_catalogue = gaps[~gaps.in_catalogue].date.dt.date.astype(str).tolist()
    if not_in_catalogue:
        print(f"  not in the RTC catalogue: {not_in_catalogue}")
    pending = cov[cov.in_catalogue & ~cov.in_stack].date.dt.date.astype(str).tolist()
    if pending:
        shown = pending if len(pending) <= 6 else [*pending[:3], "…", *pending[-2:]]
        print(f"  in the catalogue, not yet in the stack (rerun): {len(pending)} {shown}")
    print(f"→ {stack_nc}\n→ {coverage_csv}")


if __name__ == "__main__":
    main()

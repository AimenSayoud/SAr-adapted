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


def search_rtc(cfg: dict, bbox, relative_orbit: int | None = None,
               start: str | None = None, end: str | None = None) -> list:
    """Items of 'sentinel-1-rtc' over bbox; period from the config unless start/end given."""
    import planetary_computer as pc
    from pystac_client import Client

    client = Client.open(MPC_STAC, modifier=pc.sign_inplace)
    query = {"sat:relative_orbit": {"eq": relative_orbit}} if relative_orbit else None
    search = client.search(
        collections=["sentinel-1-rtc"],
        bbox=list(bbox),
        datetime=f"{start or cfg['time']['start']}/{end or cfg['time']['end']}",
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


def build_rtc_dualpol_stack(items: list, template: xr.DataArray,
                            out_nc: str | Path, checkpoint_every: int = 20) -> Path:
    """Comme build_rtc_stack mais garde VV ET VH (dB) + le ratio VH-VV (dB).

    Le ratio de polarisation croisee VH/VV separe la diffusion de VOLUME
    (vegetation -> VH eleve) du double-bounce/surface (eau, humidite -> VV
    domine). C'est le discriminant polarimetrique pour trancher mecanique
    (mouvement) vs dielectrique (humidite) sur le tapis.
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

    def _load(href):
        da = rioxarray.open_rasterio(pc.sign(href), masked=True).squeeze("band", drop=True)
        da = da.rio.clip_box(*bounds, crs="EPSG:4326")
        return da.rio.reproject_match(template, resampling=Resampling.nearest)

    new = []
    for item in sorted(items, key=lambda it: pd.to_datetime(it.datetime)):
        t = pd.to_datetime(item.datetime).tz_localize(None).normalize()
        if t in done or "vv" not in item.assets or "vh" not in item.assets:
            continue
        try:
            vv = _load(item.assets["vv"].href)
            vh = _load(item.assets["vh"].href)
            vv_db = 10 * np.log10(vv.where(vv > 0))
            vh_db = 10 * np.log10(vh.where(vh > 0))
            ds = xr.Dataset({"gamma0_vv_db": vv_db, "gamma0_vh_db": vh_db,
                             "ratio_vh_vv_db": vh_db - vv_db})
            new.append(ds.expand_dims(time=[t]))
            done.add(t)
            print(f"  + {t.date()}")
        except Exception as e:
            print(f"  ! {item.id}: {e}")
        if len(new) >= checkpoint_every:
            existing = _flush(existing, new, out_nc); new = []
            print(f"  ... checkpoint ({existing.time.size} dates)")
    if new:
        _flush(existing, new, out_nc)
    return out_nc


# ------------------------------------------------------------------ any track (descending, 2020–2024)

RTC_POLS = ("vv", "vh")


def track_relative_orbit(cfg: dict, track: str) -> int:
    """Relative orbit of a track, from the config's ``sentinel1.tracks`` block."""
    return int(cfg["sentinel1"]["tracks"][track]["relative_orbit"])


def item_day(item) -> pd.Timestamp:
    """Acquisition day (UTC) of a STAC item, as the stacks index it."""
    t = pd.Timestamp(item.datetime)
    return (t.tz_localize(None) if t.tzinfo else t).normalize()


def widen_bbox(bbox, margin_deg: float = 0.25) -> tuple[float, float, float, float]:
    """The bbox grown by a margin, for searching: STAC footprints are approximate, and the slice
    that holds the rest of a pass can have a footprint that misses the AOI (descending, spring 2020)."""
    minx, miny, maxx, maxy = bbox
    return (minx - margin_deg, miny - margin_deg, maxx + margin_deg, maxy + margin_deg)


def aoi_cover(item, bbox) -> float:
    """Share of ``bbox`` (lon/lat) inside the item's footprint (1 when either is unknown)."""
    from shapely.geometry import box, shape

    if bbox is None or not getattr(item, "geometry", None):
        return 1.0
    aoi = box(*bbox)
    return shape(item.geometry).intersection(aoi).area / aoi.area


def items_by_day(items: list, bbox=None) -> dict[pd.Timestamp, list]:
    """Dual-pol items grouped by day, each day's slices best first.

    Adjacent GRD slices of one pass are separate items of the same day. Best = largest share of
    ``bbox`` inside the footprint, then item id, so the order is deterministic. The others stay
    as fill: a footprint can claim the whole AOI while its data stops inside it.
    """
    days: dict[pd.Timestamp, list] = {}
    for it in items:
        if all(p in it.assets for p in RTC_POLS):
            days.setdefault(item_day(it), []).append(it)
    return {d: sorted(v, key=lambda it: (-aoi_cover(it, bbox), it.id)) for d, v in sorted(days.items())}


def dualpol_db(vv: xr.DataArray, vh: xr.DataArray) -> xr.Dataset:
    """Linear γ⁰ VV and VH → dB, plus the VH−VV ratio (dB); non-positive power → NaN."""
    vv_db = 10 * np.log10(vv.where(vv > 0))
    vh_db = 10 * np.log10(vh.where(vh > 0))
    return xr.Dataset({"gamma0_vv_db": vv_db, "gamma0_vh_db": vh_db,
                       "ratio_vh_vv_db": vh_db - vv_db})


def mpc_loader(template: xr.DataArray, resampling: str = "nearest",
               retries: int = 3, wait_s: float = 10.0):
    """``load(href)`` → linear γ⁰ on the template grid: signed URL, AOI window only, retried.

    ``nearest`` is what the ascending 2022–2024 stack used (one 10 m sample per 40 m cell);
    ``average`` averages the 10 m power over the cell (less speckle). Resampling happens in
    linear power, before the dB conversion.
    """
    import time

    import planetary_computer as pc
    import rioxarray  # noqa: F401 — registers .rio
    from rasterio.enums import Resampling

    bounds = template.rio.transform_bounds("EPSG:4326")
    how = Resampling[resampling]

    def load(href: str) -> xr.DataArray:
        for attempt in range(retries):
            try:
                da = rioxarray.open_rasterio(pc.sign(href), masked=True).squeeze("band", drop=True)
                da = da.rio.clip_box(*bounds, crs="EPSG:4326")
                return da.rio.reproject_match(template, resampling=how).load()
            except Exception:
                if attempt == retries - 1:
                    raise
                time.sleep(wait_s * (attempt + 1))

    return load


def _flush_items(existing: xr.Dataset | None, new: list, out_nc: Path) -> xr.Dataset:
    """``_flush`` with ``source_item`` as variable-length text: read back, it is fixed-width
    (``<U66``) and that encoding would cut a mosaic's longer ``+``-joined ids on the next write."""
    from .s2_fusion import _flush

    for ds in ([existing] if existing is not None else []) + new:
        ds["source_item"] = ds["source_item"].astype(object)
        ds["source_item"].encoding = {}
    return _flush(existing, new, out_nc)


def build_rtc_track_stack(items: list, template: xr.DataArray, out_nc: str | Path, *,
                          track: str, relative_orbit: int, bbox=None,
                          resampling: str = "nearest", load=None,
                          checkpoint_every: int = 20) -> Path:
    """Dual-pol γ⁰ stack (VV, VH, VH−VV in dB) of one track on the HyP3 crop grid.

    Same variables as ``build_rtc_dualpol_stack`` (the ascending 2022–2024 stack), plus
    ``source_item`` — the RTC item(s) used for each date, ``+``-joined — and attributes naming the
    track, orbit and resampling. Idempotent: dates already on disk are skipped; a file made for
    another track or with another resampling is refused, not extended.

    Each day is a mosaic of its slices (``items_by_day`` order): the best slice first, the next
    ones only where it has no data, stopping once the grid is full. Slices of one pass share orbit
    and processing, so the fill is radiometrically the same acquisition. A day left incomplete
    because a slice failed to load is not written (rerun to retry); one left incomplete by the
    data itself is written and reported. ``load`` defaults to ``mpc_loader``.
    """
    out_nc = Path(out_nc)
    out_nc.parent.mkdir(parents=True, exist_ok=True)
    attrs = {"source": f"Microsoft Planetary Computer, collection sentinel-1-rtc ({MPC_STAC})",
             "track": track, "relative_orbit": int(relative_orbit), "resampling": resampling,
             "grid": f"{template.rio.crs} {template.sizes['y']}x{template.sizes['x']}"}
    existing = xr.load_dataset(out_nc) if out_nc.exists() else None
    if existing is not None:
        clash = {k: (existing.attrs.get(k), attrs[k]) for k in ("track", "relative_orbit", "resampling")
                 if str(existing.attrs.get(k)) != str(attrs[k])}
        if clash:
            raise ValueError(f"{out_nc} was built with other settings {clash}; write to a new file")
    done = set(pd.to_datetime(existing.time.values)) if existing is not None else set()
    load = load or mpc_loader(template, resampling)

    new, failed, partial = [], [], []
    for day, candidates in items_by_day(items, bbox).items():
        if day in done:
            continue
        ds, used, errored = None, [], False
        for it in candidates:
            try:
                part = dualpol_db(load(it.assets["vv"].href), load(it.assets["vh"].href))
            except Exception as e:  # noqa: BLE001 — try the day's next slice; retry the day on rerun
                print(f"  ! {it.id}: {type(e).__name__}: {e}")
                errored = True
                continue
            if not bool(part.gamma0_vv_db.notnull().any()):
                continue
            ds = part if ds is None else ds.combine_first(part)
            used.append(it.id)
            if bool(ds.gamma0_vv_db.notnull().all()):
                break
        complete = ds is not None and bool(ds.gamma0_vv_db.notnull().all())
        if ds is None or (errored and not complete):
            failed.append(day)
        else:
            if not complete:
                partial.append(day)
            ds = ds.expand_dims(time=[day]).assign(source_item=("time", ["+".join(used)]))
            ds.attrs = attrs
            new.append(ds)
            print(f"  + {day.date()} {'+'.join(used)}")
        if len(new) >= checkpoint_every:
            existing = _flush_items(existing, new, out_nc)
            new = []
            print(f"  ... checkpoint ({existing.time.size} dates)")
    if new:
        _flush_items(existing, new, out_nc)
    if partial:
        print(f"  {len(partial)} day(s) written with gaps (no slice has data there): "
              f"{[d.date().isoformat() for d in partial]}")
    if failed:
        print(f"  {len(failed)} day(s) not written, rerun to retry: {[d.date().isoformat() for d in failed]}")
    return out_nc


def rtc_coverage(catalogue_days, stack_days, ifg_days) -> pd.DataFrame:
    """One row per date seen anywhere: in the RTC catalogue, in the stack, in an interferogram."""
    sets = {k: set(pd.to_datetime(list(v))) for k, v in
            (("in_catalogue", catalogue_days), ("in_stack", stack_days), ("in_interferograms", ifg_days))}
    days = sorted(set().union(*sets.values()))
    return pd.DataFrame({"date": days, **{k: [d in s for d in days] for k, s in sets.items()}})

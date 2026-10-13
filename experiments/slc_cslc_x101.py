"""X-101 — SLC pilot step 2 (Colab): the 2021 bursts geocoded onto one grid over Rzecin (coregistered CSLCs).

Geocoding every date onto the same map grid with precise orbits and a DEM is what coregisters them (the OPERA CSLC-S1
design, COMPASS). Per track, in batches of dates, so a disconnected session loses at most one batch:

1. **SAFE**: ``burst2safe <granule> --all-anns`` builds a one-burst SAFE per date from the X-100 inventory
   (``--all-anns`` is what s1-reader, COMPASS's reader, needs). Earthdata credentials from Colab's key manager.
2. **Orbits**: ``eof --force-asf`` (ASF's public orbit bucket; the default source, Copernicus Data Space, would need
   another account and prompt for it).
3. **DEM** (once): ``sardem`` Copernicus GLO-30 over the AOI plus a margin, the EGM2008 undulation added here, and the
   result checked against the cube's own heights + N on the mat (sardem's geoid handling cannot be trusted: it skips
   the shift silently, or with --keep-egm applies it the wrong way).
4. **CSLC**: COMPASS ``s1_geocode_stack.py`` writes one run file per date; each is executed. COMPASS takes its output
   grid from a burst database; a one-row-per-burst SQLite database holding the AOI box (the cube's 40 m grid extent in
   EPSG:32633) is written here, so only the AOI is geocoded, at 5 m × 10 m.
5. **Keep**: each CSLC HDF5 (a few MB over the AOI) and the static layers go to Drive
   ``insar_rzecin/slc_pilot/2021/cslc/<track>/`` as soon as the batch is done; then the batch's SAFEs are deleted.
   Dates already on Drive are skipped, so re-running after a disconnect resumes.

The QC against HyP3 on hard targets is X-101's verification and runs next (``slc_qc_x101``) on the files kept here.

    (Colab) python -m insar_wetlands.register run slc_cslc_x101 -- --tracks ascending descending
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import re
import shutil
import sqlite3
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

REPO = Path(__file__).resolve().parents[1]
INVENTORY = REPO / "data" / "slc_pilot" / "inventory_2021.csv"
WORK = Path("/content/slc_pilot")
EPSG, XPOST, YPOST = 32633, 5, 10
MARGIN_M = 200.0          # AOI margin around the cube grid (m)
DEM_MARGIN_DEG = 0.05
GEOID_N_RANGE = (25.0, 45.0)   # EGM2008 undulation around Rzecin is about 34 m
DEM_TOL_M = 3.0                # ellipsoid DEM vs cube orthometric + N on the mat
BATCH = 8                 # dates per batch: disk stays small, a disconnect loses at most one batch
DATE_RE = re.compile(r"(20\d{6})")


def sh(cmd: list[str], **kw) -> None:
    print(">>", " ".join(map(str, cmd)), flush=True)
    subprocess.run([str(c) for c in cmd], check=True, **kw)


def aoi() -> tuple[tuple[float, float, float, float], tuple[float, float, float, float]]:
    """The cube's 40 m grid extent plus a margin: (UTM 33N box, lon/lat box), each as (xmin, ymin, xmax, ymax)."""
    import xarray as xr  # noqa: PLC0415
    from pyproj import Transformer  # noqa: PLC0415

    from insar_wetlands.cube import cube_root  # noqa: PLC0415
    with xr.open_dataset(cube_root() / "silver" / "static_40m.nc") as st:
        x, y = st.x.values, st.y.values
    utm = (float(x.min()) - MARGIN_M, float(y.min()) - MARGIN_M, float(x.max()) + MARGIN_M, float(y.max()) + MARGIN_M)
    tr = Transformer.from_crs(EPSG, 4326, always_xy=True)
    lon, lat = tr.transform([utm[0], utm[2], utm[0], utm[2]], [utm[1], utm[1], utm[3], utm[3]])
    return utm, (min(lon), min(lat), max(lon), max(lat))


def credentials() -> None:
    """Earthdata login for burst2safe: the environment variables, also written to ~/.netrc."""
    # the launcher's secrets cell sets them; google.colab.userdata cannot be read from this subprocess (it hangs)
    vals = {k: os.environ.get(k) for k in ("EARTHDATA_USERNAME", "EARTHDATA_PASSWORD")}
    netrc = Path.home() / ".netrc"
    if not all(vals.values()):
        if netrc.exists() and "urs.earthdata.nasa.gov" in netrc.read_text():
            return                  # uploaded by colab_cli netrc (headless route); burst2safe reads it
        raise SystemExit("Earthdata credentials missing: run the launcher's secrets cell (browser route) or "
                         "`python -m insar_wetlands.colab_cli netrc` (terminal route)")
    netrc.write_text(f"machine urs.earthdata.nasa.gov login {vals['EARTHDATA_USERNAME']} "
                     f"password {vals['EARTHDATA_PASSWORD']}\n")
    netrc.chmod(0o600)


def jpl_id(burst_id: str) -> str:
    """'175_374052_IW1' → 't175_374052_iw1', the form s1-reader prints and COMPASS looks up."""
    t, n, sw = burst_id.split("_")
    return f"t{int(t):03d}_{n}_{sw.lower()}"


def burst_db(path: Path, ids: list[str], utm_box) -> Path:
    """The burst database COMPASS reads its output grid from, holding only the AOI box for our bursts."""
    path.unlink(missing_ok=True)
    con = sqlite3.connect(path)
    con.execute("CREATE TABLE burst_id_map (burst_id_jpl TEXT PRIMARY KEY, epsg INTEGER, "
                "xmin REAL, ymin REAL, xmax REAL, ymax REAL)")
    for i in ids:
        con.execute("INSERT INTO burst_id_map VALUES (?, ?, ?, ?, ?, ?)", (i, EPSG, *utm_box))
    con.commit()
    con.close()
    return path


def geoid_undulation(lon, lat) -> np.ndarray:
    """EGM2008 geoid height above the WGS84 ellipsoid (m), from PROJ's grid (us_nga_egm08_25.tif)."""
    from pyproj import Transformer  # noqa: PLC0415
    t = Transformer.from_crs("EPSG:4326+3855", "EPSG:4979", always_xy=True)
    _, _, h = t.transform(np.asarray(lon, float), np.asarray(lat, float), np.zeros(np.shape(lon)))
    return np.asarray(h, float)


def to_ellipsoid(src: Path, dst: Path) -> Path:
    """Copernicus heights are above EGM2008; COMPASS needs WGS84 ellipsoid heights. sardem's own conversion targets a 2-D
    CRS, so GDAL silently drops the vertical shift (checked 2026-10-13: the DEM came out orthometric); the undulation
    is added here, and refused if PROJ returns an implausible value (no grid gives 0 without an error)."""
    import rasterio  # noqa: PLC0415
    with rasterio.open(src) as r:
        z = r.read(1).astype("float64")
        prof = r.profile
        rows, cols = np.mgrid[0:r.height, 0:r.width]
        lon, lat = rasterio.transform.xy(r.transform, rows.ravel(), cols.ravel(), offset="center")
    n = geoid_undulation(lon, lat).reshape(z.shape)
    if not (GEOID_N_RANGE[0] <= np.nanmin(n) and np.nanmax(n) <= GEOID_N_RANGE[1]):
        raise RuntimeError(f"geoid undulation {np.nanmin(n):.2f}–{np.nanmax(n):.2f} m outside {GEOID_N_RANGE}: "
                           "PROJ has no EGM2008 grid (projsync --file us_nga_egm08_25.tif)")
    prof.update(dtype="float32", crs="EPSG:4326")   # ellipsoidal now; sardem's compound CRS (EPSG:9518) is unknown to isce3
    with rasterio.open(dst, "w", **prof) as w:
        w.write((z + n).astype("float32"), 1)
    print(f"DEM to ellipsoid: + {np.nanmean(n):.2f} m (EGM2008)", flush=True)
    return dst


def check_ellipsoid(dem_file: Path) -> float:
    """Independent check: on the mat, the DEM must equal the cube's own Copernicus heights (orthometric, EGM2008) plus
    the geoid undulation. sardem silently skips its geoid shift, and with --keep-egm shifts the wrong way (both seen on
    the VM 2026-10-13); only this comparison, not the tool's word, decides that the heights are ellipsoidal."""
    import rasterio  # noqa: PLC0415
    import xarray as xr  # noqa: PLC0415
    from pyproj import Transformer  # noqa: PLC0415

    from insar_wetlands.cube import cube_root  # noqa: PLC0415
    with xr.open_dataset(cube_root() / "silver" / "static_40m.nc") as st:
        h = st.dem_glo30_m.where(st.zone == 1)
        stacked = h.stack(p=("y", "x")).dropna("p")
        sel = stacked.isel(p=np.linspace(0, stacked.sizes["p"] - 1, 40).astype(int))
        xs, ys, ortho = sel.x.values, sel.y.values, sel.values
    lon, lat = Transformer.from_crs(EPSG, 4326, always_xy=True).transform(xs, ys)
    with rasterio.open(dem_file) as r:
        got = np.array([v[0] for v in r.sample(list(zip(lon, lat)))], float)
    diff = float(np.median(got - (ortho + geoid_undulation(lon, lat))))
    if abs(diff) > DEM_TOL_M:
        raise RuntimeError(f"DEM is not ellipsoidal: median (DEM - (cube orthometric + N)) = {diff:.2f} m on the mat")
    print(f"DEM check on the mat: median difference to cube orthometric + N = {diff:.2f} m", flush=True)
    return diff


def dem(lonlat_box) -> Path:
    f = WORK / "dem_ellipsoid.tif"
    if not f.exists():
        w, s, e, n = lonlat_box
        raw = WORK / "dem_sardem.tif"             # sardem's default output: Copernicus heights, orthometric in practice
        sh(["sardem", "--bbox", w - DEM_MARGIN_DEG, s - DEM_MARGIN_DEG, e + DEM_MARGIN_DEG, n + DEM_MARGIN_DEG,
            "--data-source", "COP", "-o", raw])
        to_ellipsoid(raw, f)
    check_ellipsoid(f)
    return f


def drive_dest(track: str) -> Path:
    from insar_wetlands.paths import resolve_drive  # noqa: PLC0415
    d = resolve_drive() / "slc_pilot" / "2021" / "cslc" / track
    d.mkdir(parents=True, exist_ok=True)
    return d


def done_dates(dest: Path) -> set[str]:
    return {m.group(1) for f in dest.glob("*.h5") if "static" not in f.name.lower()
            for m in [DATE_RE.search(f.name)] if m}


def build_safes(granules: list[str], safe_dir: Path) -> None:
    safe_dir.mkdir(parents=True, exist_ok=True)
    for g in granules:
        date = DATE_RE.search(g).group(1)
        if glob.glob(str(safe_dir / f"*{date}*.SAFE")):
            continue
        sh(["burst2safe", g, "--all-anns", "--output-dir", safe_dir], cwd=safe_dir)
    for f in list(safe_dir.glob("*.tiff")) + list(safe_dir.glob("*.xml")):     # burst downloads left beside the SAFEs
        f.unlink()


def fix_runconfig(path: Path) -> None:
    """s1_geocode_stack.py takes -p as a list and writes ``polarization: [co-pol]``, which COMPASS's own schema
    (s1_cslc_geo.yaml) rejects — it wants the string. Checked 2026-10-13 on COMPASS 0.5.6: every date failed."""
    cfg = yaml.safe_load(path.read_text())
    proc = cfg["runconfig"]["groups"]["processing"]
    if isinstance(proc.get("polarization"), list):
        proc["polarization"] = proc["polarization"][0]
        path.write_text(yaml.safe_dump(cfg, sort_keys=False))


def geocode(bid: str, safe_dir: Path, dem_file: Path, db: Path, lonlat_box, pol: str,
            corrections: bool) -> tuple[Path, list[str], bool]:
    """COMPASS on one batch. Returns (work dir, dates done, corrections used). If a date fails with COMPASS's LUT
    corrections (tides, ionosphere, … — optional inputs), that date is redone without them and recorded so."""
    work = safe_dir.parent / "compass"
    orbit_dir = safe_dir.parent / "orbits"
    orbit_dir.mkdir(exist_ok=True)
    sh(["eof", "--search-path", safe_dir, "--save-dir", orbit_dir, "--force-asf"])
    w, s, e, n = lonlat_box

    def stack(flag: bool) -> list[str]:
        shutil.rmtree(work, ignore_errors=True)
        sh(["s1_geocode_stack.py", "-s", safe_dir, "--unzipped", "-d", dem_file, "-o", orbit_dir, "-w", work,
            "-b", bid, "-p", "co-pol" if pol in ("VV", "HH") else "cross-pol", "-dx", XPOST, "-dy", YPOST,
            "--bbox", w, s, e, n, "-e", EPSG, "--burst-db-file", db, *([] if flag else ["-nc"])])
        runs = sorted(glob.glob(str(work / "run_files" / "run_*.sh")))
        if not runs:
            raise RuntimeError(f"COMPASS wrote no run files for {bid} in {safe_dir}")
        for rc in glob.glob(str(work / "runconfigs" / "*.yaml")):
            fix_runconfig(Path(rc))
        return runs

    runs, failed = stack(corrections), []
    for rf in runs:
        try:
            sh(["bash", rf])
        except subprocess.CalledProcessError:
            if not corrections:
                raise
            failed.append(rf)
    if failed:                      # redo the whole batch without corrections, so a batch is made one way
        print(f"{len(failed)} run(s) failed with corrections; redoing this batch with -nc", flush=True)
        for rf in stack(False):
            sh(["bash", rf])
        corrections = False
    return work, [DATE_RE.search(Path(r).name).group(1) for r in runs], corrections


def keep(work: Path, dest: Path, track: str, corrections: bool) -> list[dict]:
    """Copy every CSLC (and the static layers) to Drive, the date always in the file name so a re-run can skip it."""
    rows = []
    for h5 in sorted(Path(work).rglob("*.h5")):
        rel = h5.relative_to(work)
        if "scratch" in rel.parts:
            continue
        m = DATE_RE.search(str(rel))
        name = h5.name if (DATE_RE.search(h5.name) or "static" in h5.name.lower() or not m) else f"{m.group(1)}_{h5.name}"
        out = dest / name
        if not out.exists():
            shutil.copy2(h5, out)
        rows.append({"track": track, "date": m.group(1) if m else "", "file": out.name, "bytes": out.stat().st_size,
                     "lut_corrections": corrections})
    return rows


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--tracks", nargs="+", default=["ascending", "descending"])
    ap.add_argument("--pol", default="VV", choices=["VV", "VH"])
    ap.add_argument("--max-dates", type=int, default=0, help="process at most this many new dates per track (test)")
    a = ap.parse_args(argv)
    credentials()
    utm_box, lonlat_box = aoi()
    print("AOI (EPSG:32633)", utm_box, "lon/lat", lonlat_box)
    tracks = yaml.safe_load((REPO / "config" / "config.yaml").read_text())["sentinel1"]["tracks"]
    WORK.mkdir(parents=True, exist_ok=True)
    db = burst_db(WORK / "burst_aoi.sqlite3", [jpl_id(tracks[t]["burst_id"]) for t in a.tracks], utm_box)
    dem_file = dem(lonlat_box)
    inv = pd.read_csv(INVENTORY)
    inv["date"] = pd.to_datetime(inv.time_utc).dt.strftime("%Y%m%d")
    summary = {}
    for t in a.tracks:
        dest = drive_dest(t)
        g = inv[(inv.track == t) & (inv.polarization == a.pol)].drop_duplicates("date").sort_values("date")
        todo = g[~g.date.isin(done_dates(dest))]
        if a.max_dates:
            todo = todo.head(a.max_dates)
        print(f"{t}: {len(g)} dates in the inventory, {len(done_dates(dest))} already on Drive, {len(todo)} to do now",
              flush=True)
        rows = []
        for i in range(0, len(todo), BATCH):
            batch = todo.iloc[i:i + BATCH]
            safe_dir = WORK / t / f"batch_{i // BATCH:02d}" / "safe"
            build_safes(batch.granule.tolist(), safe_dir)
            work, _, used = geocode(jpl_id(tracks[t]["burst_id"]), safe_dir, dem_file, db, lonlat_box, a.pol,
                                    corrections=True)
            rows += keep(work, dest, t, used)
            shutil.rmtree(safe_dir.parent, ignore_errors=True)
            print(f"{t}: batch {i // BATCH + 1} kept; {len(done_dates(dest))} dates on Drive", flush=True)
        man = pd.DataFrame(rows)
        if len(man):
            old = dest / "manifest.csv"
            if old.exists():
                man = pd.concat([pd.read_csv(old), man]).drop_duplicates("file", keep="last")
            man.to_csv(old, index=False)
        summary[t] = len(done_dates(dest))
    print(json.dumps({"dates_on_drive": summary}))


if __name__ == "__main__":
    main()

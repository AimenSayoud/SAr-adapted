"""X-101 — SLC pilot step 2 (Colab): the 2021 bursts geocoded onto one grid over Rzecin (coregistered CSLCs).

One Colab session, because the VM disk is wiped between sessions and the raw bursts (tens of GB) never go to Drive:

1. **SAFE**: ``burst2safe`` builds a minimal SAFE per date from the X-100 inventory's burst granules (VV by default).
2. **Orbits**: ``eof`` (sentineleof) fetches the precise orbit of each SAFE.
3. **DEM**: ``sardem`` makes a Copernicus GLO-30 DEM over the box plus a margin.
4. **CSLC**: COMPASS ``s1_geocode_stack.py`` writes one run-config per date and burst on a common geocoded grid
   (EPSG:32633, 5 m × 10 m) over the AOI box; the run files are executed. Geocoding every date onto the same grid with
   precise orbits and the DEM is what coregisters them (OPERA CSLC-S1 design).
5. **Keep**: the CSLC HDF5 files (one per date, a few MB each over the box) are copied to Drive
   ``insar_rzecin/slc_pilot/2021/cslc/<track>/``, with a manifest; the raw bursts and SAFEs stay on the VM.

The QC against HyP3 (multilooked, filtered pairs of the same dates on hard targets) is X-101's verification and runs
next, in ``slc_qc_x101`` (Mac or Colab), on the files kept here.

**First run = the test of this chain**: these tools were not run on this project before. A failure is recorded in the
run record; fix and re-run from the step that failed (each step skips what it already produced).

    (Colab) python -m insar_wetlands.register run slc_cslc_x101 -- --tracks ascending descending
"""
from __future__ import annotations

import argparse
import glob
import json
import shutil
import subprocess
from pathlib import Path

import pandas as pd
import yaml

REPO = Path(__file__).resolve().parents[1]
INVENTORY = REPO / "data" / "slc_pilot" / "inventory_2021.csv"
WORK = Path("/content/slc_pilot")
EPSG, XPOST, YPOST = 32633, 5, 10
MARGIN_DEG = 0.05


def sh(cmd: list[str] | str, **kw) -> None:
    print(">>", cmd if isinstance(cmd, str) else " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True, shell=isinstance(cmd, str), **kw)


def aoi_bbox() -> tuple[float, float, float, float]:
    """(west, south, east, north) in degrees: the cube's 40 m grid (EPSG:32633), so the CSLCs cover the mat, its
    references (grassland, hard targets) and every cell the HyP3 comparison uses."""
    import xarray as xr  # noqa: PLC0415
    from pyproj import Transformer  # noqa: PLC0415

    from insar_wetlands.cube import cube_root  # noqa: PLC0415
    with xr.open_dataset(cube_root() / "silver" / "static_40m.nc") as st:
        x, y = st.x.values, st.y.values
    tr = Transformer.from_crs(EPSG, 4326, always_xy=True)
    lon, lat = tr.transform([x.min() - 20, x.max() + 20, x.min() - 20, x.max() + 20],
                            [y.min() - 20, y.min() - 20, y.max() + 20, y.max() + 20])
    return min(lon), min(lat), max(lon), max(lat)


def credentials() -> dict:
    out = {}
    try:
        from google.colab import userdata  # noqa: PLC0415
        for k in ("EARTHDATA_USERNAME", "EARTHDATA_PASSWORD"):
            out[k] = userdata.get(k)
    except Exception:
        import os
        out = {k: os.environ.get(k) for k in ("EARTHDATA_USERNAME", "EARTHDATA_PASSWORD")}
    if not all(out.values()):
        raise SystemExit("EARTHDATA_USERNAME / EARTHDATA_PASSWORD missing from Colab's key manager")
    netrc = Path.home() / ".netrc"
    netrc.write_text(f"machine urs.earthdata.nasa.gov login {out['EARTHDATA_USERNAME']} "
                     f"password {out['EARTHDATA_PASSWORD']}\n")
    netrc.chmod(0o600)
    return out


def safes(track: str, pols: list[str]) -> Path:
    inv = pd.read_csv(INVENTORY)
    inv = inv[(inv.track == track) & inv.polarization.isin(pols)]
    d = WORK / track / "safe"
    d.mkdir(parents=True, exist_ok=True)
    inv["date"] = pd.to_datetime(inv.time_utc).dt.strftime("%Y%m%d")
    for date, g in inv.groupby("date"):
        if glob.glob(str(d / f"*{date}*.SAFE")):
            continue
        sh(["burst2safe", *g.granule.tolist(), "--output-dir", str(d)])
    return d


def orbits(safe_dir: Path) -> Path:
    o = safe_dir.parent / "orbits"
    o.mkdir(exist_ok=True)
    sh(["eof", "--search-path", str(safe_dir), "--save-dir", str(o)])
    return o


def dem(bbox) -> Path:
    f = WORK / "dem.tif"
    if not f.exists():
        w, s, e, n = bbox
        sh(["sardem", "--bbox", str(w - MARGIN_DEG), str(s - MARGIN_DEG), str(e + MARGIN_DEG), str(n + MARGIN_DEG),
            "--data-source", "COP", "-o", str(f)])
    return f


def cslc(track: str, burst_id: str, safe_dir: Path, orbit_dir: Path, dem_file: Path, bbox, pol: str) -> Path:
    work = WORK / track / "compass"
    work.mkdir(parents=True, exist_ok=True)
    bid = "t" + burst_id.lower()
    w, s, e, n = bbox
    sh(["s1_geocode_stack.py", "-s", str(safe_dir), "-o", str(orbit_dir), "-d", str(dem_file), "-w", str(work),
        "-b", bid, "-p", pol, "--bbox", str(w), str(s), str(e), str(n), "-e", str(EPSG),
        "-x", str(XPOST), "-y", str(YPOST)])
    for rf in sorted(glob.glob(str(work / "run_files" / "*.sh"))):
        sh(["bash", rf])
    return work


def keep(track: str, work: Path) -> pd.DataFrame:
    from insar_wetlands.paths import resolve_drive  # noqa: PLC0415
    dest = resolve_drive() / "slc_pilot" / "2021" / "cslc" / track
    dest.mkdir(parents=True, exist_ok=True)
    rows = []
    for h5 in sorted(glob.glob(str(work / "**" / "*.h5"), recursive=True)):
        if "static" in Path(h5).name.lower():
            continue
        out = dest / Path(h5).name
        if not out.exists():
            shutil.copy2(h5, out)
        rows.append({"track": track, "file": str(out), "bytes": out.stat().st_size})
    man = pd.DataFrame(rows)
    man.to_csv(dest / "manifest.csv", index=False)
    print(f"{track}: {len(man)} CSLC files, {man.bytes.sum() / 1e6:.0f} MB → {dest}")
    return man


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--tracks", nargs="+", default=["ascending", "descending"])
    ap.add_argument("--pols", nargs="+", default=["VV"])
    a = ap.parse_args(argv)
    credentials()
    bbox = aoi_bbox()
    print("AOI bbox", bbox)
    bursts = yaml.safe_load((REPO / "config" / "config.yaml").read_text())["sentinel1"]["tracks"]
    dem_file = dem(bbox)
    summary = {}
    for t in a.tracks:
        sd = safes(t, a.pols)
        od = orbits(sd)
        work = cslc(t, bursts[t]["burst_id"], sd, od, dem_file, bbox, a.pols[0])
        summary[t] = len(keep(t, work))
    print(json.dumps(summary))


if __name__ == "__main__":
    main()

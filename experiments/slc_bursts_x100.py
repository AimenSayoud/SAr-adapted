"""X-100 — SLC pilot step 1: the Sentinel-1 single-burst SLCs over Rzecin for one laser season.

Two modes:

- ``inventory`` (Mac or Colab; search only, no credentials): queries ASF's SLC-BURST catalogue for the two bursts in
  ``config/config.yaml`` (ascending 175_374052_IW1, descending 022_045867_IW2), VV and VH, from 2021-01-01 to
  2022-02-28 (the season plus early-2022 dates for the long pairs of X-104), and writes the list — public ASF
  metadata only — to ``data/slc_pilot/inventory_2021.csv`` in this repository, so a Colab session has it with the clone.
- ``download`` (Colab; Earthdata credentials from Colab's key manager): fetches the listed bursts to the VM's local
  disk (``/content/slc_pilot``), never to Drive. One season is tens of GB; only the geocoded crops made from it by
  X-101 go to Drive. Each burst comes as its SLC (``.tiff``) and its annotation (``.xml``), both needed to geocode it;
  an empty or missing file is fetched again. A manifest of what was fetched goes to
  Drive ``insar_rzecin/slc_pilot/2021/download_manifest.csv``.

    make run E=slc_bursts_x100 ARGS="inventory"
    (Colab) python -m insar_wetlands.register run slc_bursts_x100 -- download
"""
from __future__ import annotations

import argparse
import io
import os
import sys
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd
import yaml

REPO = Path(__file__).resolve().parents[1]
ASF = "https://api.daac.asf.alaska.edu/services/search/param"
INVENTORY = REPO / "data" / "slc_pilot" / "inventory_2021.csv"
START, END = "2021-01-01", "2022-02-28"
POLS = ("VV", "VH")
LOCAL = Path("/content/slc_pilot")


def bursts() -> dict[str, str]:
    cfg = yaml.safe_load((REPO / "config" / "config.yaml").read_text(encoding="utf-8"))
    return {t: v["burst_id"] for t, v in cfg["sentinel1"]["tracks"].items()}


def inventory() -> pd.DataFrame:
    rows = []
    for track, bid in bursts().items():
        for pol in POLS:
            q = urllib.parse.urlencode({"dataset": "SLC-BURST", "fullBurstID": bid, "polarization": pol,
                                        "start": START, "end": END, "output": "csv"})
            with urllib.request.urlopen(f"{ASF}?{q}", timeout=300) as r:
                d = pd.read_csv(io.BytesIO(r.read()))
            if d.empty:
                continue
            rows.append(pd.DataFrame({
                "track": track, "burst_id": bid, "polarization": pol,
                "time_utc": pd.to_datetime(d["Start Time"], utc=True),
                "platform": d["Platform"], "granule": d["Granule Name"], "url": d["URL"],
                "samples": pd.to_numeric(d["samplesPerBurst"], errors="coerce"),
                # bursts are cut on demand and the catalogue gives no size: about 1 500 lines x samples x 4 bytes (CInt16)
                "est_bytes": pd.to_numeric(d["samplesPerBurst"], errors="coerce") * 1500 * 4}))
    inv = pd.concat(rows, ignore_index=True).drop_duplicates("granule").sort_values(["track", "polarization", "time_utc"])
    INVENTORY.parent.mkdir(parents=True, exist_ok=True)
    inv.to_csv(INVENTORY, index=False)
    print(inv.groupby(["track", "polarization"]).agg(n=("granule", "size"), first=("time_utc", "min"),
                                                      last=("time_utc", "max"), est_gb=("est_bytes", lambda b: b.sum() / 1e9)))
    return inv


def _session():
    """An authenticated session from Colab secrets or the environment (EARTHDATA_TOKEN, or USERNAME/PASSWORD)."""
    import requests
    token, user, pwd = (os.environ.get(k) for k in ("EARTHDATA_TOKEN", "EARTHDATA_USERNAME", "EARTHDATA_PASSWORD"))
    try:
        from google.colab import userdata  # noqa: PLC0415 — only on Colab
        token = token or userdata.get("EARTHDATA_TOKEN")
        user = user or userdata.get("EARTHDATA_USERNAME")
        pwd = pwd or userdata.get("EARTHDATA_PASSWORD")
    except Exception:
        pass
    s = requests.Session()
    if token:
        s.headers["Authorization"] = f"Bearer {token}"
    elif user and pwd:
        s.auth = (user, pwd)
    else:
        sys.exit("no Earthdata credentials: set EARTHDATA_TOKEN (or USERNAME/PASSWORD) in Colab's key manager")
    return s


def download(dest: Path = LOCAL, tracks: tuple[str, ...] = ("ascending", "descending"),
             pols: tuple[str, ...] = ("VV",)) -> pd.DataFrame:
    """VV by default (the 6-day phase, PS, SHP, long pairs); add VH only for X-105 (the scattering regime)."""
    inv = pd.read_csv(INVENTORY)
    inv = inv[inv.track.isin(tracks) & inv.polarization.isin(pols)]
    s = _session()
    dest.mkdir(parents=True, exist_ok=True)
    out = []
    for r in inv.itertuples():
        base = dest / r.track / r.polarization / r.granule
        base.parent.mkdir(parents=True, exist_ok=True)
        files = {}
        for ext, url in ((".tiff", r.url), (".xml", r.url.rsplit(".", 1)[0] + ".xml")):
            f = base.with_suffix(ext)
            if not f.exists() or f.stat().st_size == 0:
                with s.get(url, stream=True, timeout=900) as resp:
                    resp.raise_for_status()
                    with open(f, "wb") as fh:
                        for chunk in resp.iter_content(1 << 20):
                            fh.write(chunk)
            files[ext] = f
        out.append({"granule": r.granule, "track": r.track, "polarization": r.polarization, "time_utc": r.time_utc,
                    "slc": str(files[".tiff"]), "annotation": str(files[".xml"]),
                    "bytes": files[".tiff"].stat().st_size + files[".xml"].stat().st_size})
        print(f"{len(out)}/{len(inv)} {r.granule}", flush=True)
    man = pd.DataFrame(out)
    from insar_wetlands.paths import resolve_drive
    drive_dir = resolve_drive() / "slc_pilot" / "2021"
    drive_dir.mkdir(parents=True, exist_ok=True)
    man.to_csv(drive_dir / "download_manifest.csv", index=False)
    print(f"{len(man)} bursts, {man.bytes.sum() / 1e9:.1f} GB on {dest}; manifest on {drive_dir}")
    return man


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("mode", choices=["inventory", "download"])
    ap.add_argument("--dest", default=str(LOCAL))
    ap.add_argument("--pols", nargs="+", default=["VV"], choices=POLS)
    ap.add_argument("--tracks", nargs="+", default=["ascending", "descending"])
    a = ap.parse_args(argv)
    inventory() if a.mode == "inventory" else download(Path(a.dest), tuple(a.tracks), tuple(a.pols))


if __name__ == "__main__":
    main()

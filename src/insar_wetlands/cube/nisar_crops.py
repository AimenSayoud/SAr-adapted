"""NISAR L-band GUNW crops over Rzecin for the cube: every beta and provisional granule, the template window only (HTTP
range requests with the Earthdata login in ~/.netrc, as X-068/X-069), stored per granule in ``raw/nisar`` and gathered on
the 40 m template in ``silver/nisar_lband.nc``.

Layers: unwrapped phase, coherence, ionosphere screen, connected components (80 m grid → the 40 m cell's parent value);
wrapped phase and its coherence (20 m grid → 40 m phasor / mean). Phase in radians, unreferenced; the wavelength per pair.
"""
from __future__ import annotations

import io
import json
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from . import core

CMR = "https://cmr.earthdata.nasa.gov/search/granules.json"
COLLECTIONS = {"beta": "C2850261892-ASF", "provisional": "C2854335566-ASF"}
GRID = "science/LSAR/GUNW/grids/frequencyA"
C_LIGHT = 299_792_458.0
UNW_KEYS = ("unwrappedPhase", "coherenceMagnitude", "ionospherePhaseScreen", "connectedComponents")


class HTTPRange(io.RawIOBase):
    """A read-only file over HTTP range requests, 1 MB blocks cached (Earthdata auth from ~/.netrc on redirect)."""

    def __init__(self, url: str, block: int = 1 << 20):
        import requests
        self.s = requests.Session()
        r = self.s.get(url, headers={"Range": "bytes=0-0"}, allow_redirects=True, timeout=120)
        r.raise_for_status()
        self.url, self.size = r.url, int(r.headers["Content-Range"].split("/")[-1])
        self.pos, self.block, self.cache, self.requests = 0, block, {}, 0

    def readable(self):
        return True

    def seekable(self):
        return True

    def tell(self):
        return self.pos

    def seek(self, off, whence=0):
        self.pos = off if whence == 0 else self.pos + off if whence == 1 else self.size + off
        return self.pos

    def _blk(self, i):
        if i not in self.cache:
            a = i * self.block
            r = self.s.get(self.url, headers={"Range": f"bytes={a}-{min(a + self.block, self.size) - 1}"}, timeout=120)
            r.raise_for_status()
            self.cache[i] = r.content
            self.requests += 1
        return self.cache[i]

    def readinto(self, buf):
        n = min(len(buf), self.size - self.pos)
        if n <= 0:
            return 0
        out = bytearray()
        while len(out) < n:
            i, o = divmod(self.pos + len(out), self.block)
            out += self._blk(i)[o:o + n - len(out)]
        buf[:n] = out
        self.pos += n
        return n


def granules(bbox) -> list[dict]:
    out = []
    for kind, coll in COLLECTIONS.items():
        page = 1
        while True:
            q = urllib.parse.urlencode({"collection_concept_id": coll, "bounding_box": ",".join(map(str, bbox)),
                                        "page_size": 200, "page_num": page})
            with urllib.request.urlopen(f"{CMR}?{q}", timeout=120) as r:
                feed = json.load(r)["feed"]["entry"]
            for e in feed:
                h5 = [lk["href"] for lk in e.get("links", []) if lk.get("href", "").endswith(".h5")]
                if h5:
                    out.append({"title": e["title"], "url": h5[0], "kind": kind})
            if len(feed) < 200:
                break
            page += 1
    return sorted({g["title"]: g for g in out}.values(), key=lambda g: g["title"])


class OutsideGrid(ValueError):
    """The granule's footprint touches the site's box but its product grid does not reach the template."""


def _window(c: np.ndarray, lo: float, hi: float) -> slice:
    k = np.flatnonzero((c >= lo) & (c <= hi))
    if not len(k):
        raise OutsideGrid("product grid does not reach the template window")
    return slice(int(k.min()), int(k.max()) + 1)


def crop(g: dict, bounds, out: Path) -> dict:
    import h5py
    f = HTTPRange(g["url"])
    h = h5py.File(io.BufferedReader(f, buffer_size=1 << 20), "r")
    lam = C_LIGHT / float(h[f"{GRID}/centerFrequency"][()])
    idn = h["science/LSAR/identification"]
    t1 = idn["referenceZeroDopplerStartTime"][()].decode()[:19]
    t2 = idn["secondaryZeroDopplerStartTime"][()].decode()[:19]
    left, bottom, right, top = bounds
    arrays = {"wavelength_m": lam, "t1": t1, "t2": t2}
    for grid, keys, tag in (("unwrappedInterferogram", UNW_KEYS, "u"),
                            ("wrappedInterferogram", ("wrappedInterferogram", "coherenceMagnitude"), "w")):
        gx, gy = h[f"{GRID}/{grid}/xCoordinates"][:], h[f"{GRID}/{grid}/yCoordinates"][:]
        sx, sy = _window(gx, left - 100, right + 100), _window(gy, bottom - 100, top + 100)
        arrays[f"{tag}_x"], arrays[f"{tag}_y"] = gx[sx], gy[sy]
        for key in keys:
            path = f"{GRID}/{grid}/HH/{key}"
            if path not in h:
                continue
            a = h[path][sy, sx]
            if np.iscomplexobj(a):
                a = np.angle(a)
            arrays[f"{tag}_{key}"] = a.astype("float32")
    np.savez_compressed(out, **{k: np.asarray(v) for k, v in arrays.items()})
    return {"title": g["title"], "kind": g["kind"], "t1": t1, "t2": t2, "requests": f.requests, "file": out.name}


def fetch(raw_dir: Path, ctx, log=print) -> pd.DataFrame:
    raw_dir.mkdir(parents=True, exist_ok=True)
    tpl = ctx.template
    tx, ty = tpl.x.values.astype(float), tpl.y.values.astype(float)
    bounds = (tx.min() - 20, ty.min() - 20, tx.max() + 20, ty.max() + 20)
    gs = granules(tpl.rio.transform_bounds("EPSG:4326"))
    log(f"NISAR: {len(gs)} granules over the site")
    rows = []
    for i, g in enumerate(gs, 1):
        out = raw_dir / f"{g['title']}.npz"
        if out.exists():
            z = np.load(out)
            rows.append({"title": g["title"], "kind": g["kind"], "t1": str(z["t1"]), "t2": str(z["t2"]), "file": out.name})
            continue
        try:
            rows.append(crop(g, bounds, out))
            log(f"  {i}/{len(gs)} {g['title'][:60]}")
        except OutsideGrid:
            log(f"  - {g['title'][:60]}: product grid does not reach the site (skipped)")
        except Exception as e:  # noqa: BLE001 — an unreadable granule is reported, not fatal
            log(f"  ! {g['title']}: {type(e).__name__}: {e}")
    tab = pd.DataFrame(rows)
    tab.to_csv(raw_dir / "granules.csv", index=False)
    return tab


def to_silver(raw_dir: Path, tx, ty, out: Path, log=print) -> None:
    """Every cropped granule on the 40 m template: unwrapped-grid layers by parent cell, wrapped grid as a phasor."""
    from ..regrid import to_template
    tab = pd.read_csv(raw_dir / "granules.csv")
    layers = {k: [] for k in ("unw", "coh", "iono", "conncomp", "wrapped", "coh_w")}
    keep = []
    for r in tab.itertuples():
        z = np.load(raw_dir / r.file)

        def parent(a, x, y):
            col = np.clip(np.rint((tx - x[0]) / (x[1] - x[0])).astype(int), 0, len(x) - 1)
            row = np.clip(np.rint((ty - y[0]) / (y[1] - y[0])).astype(int), 0, len(y) - 1)
            return a[np.ix_(row, col)].astype("float32")

        if "u_unwrappedPhase" not in z:
            continue
        ux, uy = z["u_x"], z["u_y"]
        layers["unw"].append(parent(z["u_unwrappedPhase"], ux, uy))
        layers["coh"].append(parent(z["u_coherenceMagnitude"], ux, uy))
        layers["iono"].append(parent(z["u_ionospherePhaseScreen"], ux, uy) if "u_ionospherePhaseScreen" in z
                              else np.full((len(ty), len(tx)), np.nan, "float32"))
        layers["conncomp"].append(parent(z["u_connectedComponents"], ux, uy) if "u_connectedComponents" in z
                                  else np.full((len(ty), len(tx)), np.nan, "float32"))
        if "w_wrappedInterferogram" in z:
            layers["wrapped"].append(to_template(z["w_wrappedInterferogram"], z["w_x"], z["w_y"], tx, ty, kind="phase"))
            layers["coh_w"].append(to_template(z["w_coherenceMagnitude"], z["w_x"], z["w_y"], tx, ty, kind="mean"))
        else:
            layers["wrapped"].append(np.full((len(ty), len(tx)), np.nan))
            layers["coh_w"].append(np.full((len(ty), len(tx)), np.nan))
        keep.append(r)
    if not keep:
        log("NISAR: nothing to harmonise")
        return
    k = pd.DataFrame(keep)
    t1, t2 = pd.to_datetime(k.t1, utc=True), pd.to_datetime(k.t2, utc=True)
    parts = k.title.str.split("_")
    def txt(v):
        return np.asarray([str(x) for x in v], dtype=object)

    ds = xr.Dataset({n: (("pair", "y", "x"), np.stack(v).astype("float32")) for n, v in layers.items()},
                    coords={"pair": txt(k.title), "x": tx, "y": ty,
                            "t1": ("pair", t1.dt.tz_convert(None).to_numpy()),
                            "t2": ("pair", t2.dt.tz_convert(None).to_numpy()),
                            "revisit_days": ("pair", ((t2 - t1).dt.total_seconds() / 86400).round().astype(int).to_numpy()),
                            "track_no": ("pair", txt(parts.str[5])), "direction": ("pair", txt(parts.str[6])),
                            "product": ("pair", txt(k.kind)),
                            "fold_year": ("pair", core.pair_year_fold(t1, t2)),
                            "valid_share": ("pair", np.isfinite(np.stack(layers["unw"])).mean((1, 2)))})
    ds.attrs = {"source": "NISAR L2 GUNW (beta, provisional), HH, frequency A; 80 m unwrapped grid by parent cell, "
                          "20 m wrapped grid as a phasor", "phase_units": "rad, unreferenced",
                "wavelength_m": float(np.load(raw_dir / k.file.iloc[0])["wavelength_m"])}
    ds.to_netcdf(out, encoding={v: {"dtype": "float32", "zlib": True, "complevel": 4} for v in ds.data_vars})
    log(f"NISAR: {ds.sizes['pair']} pairs → {out.name}")

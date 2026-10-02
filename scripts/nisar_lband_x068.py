"""X-068: NISAR L-band over Rzecin, Oct 2025 – Jan 2026 — does the 24-cm radar see the mat better than C-band?

X-035 (2026-09-13) found no NISAR product over Rzecin; NASA CMR now lists beta GUNW interferograms (track 158,
ascending, 03:39 UTC, 12-day consecutive pairs, 2025-10-28 → 2026-01-20, full AOI cover). This script reads only the
AOI window of each 2.5 GB HDF5 file by HTTP range requests (Earthdata login from ~/.netrc; no download) and reports,
per pair: coherence by zone (80 m unwrapped grid and 20 m wrapped grid), the P6 line-of-sight change (P6 − grassland
median, with and without the product's ionosphere screen) and the laser change at the overpasses (laser to
2025-12-31; snow mask as fusion v2), next to the Sentinel-1 C-band pairs of the same weeks (X-066 cache).

Beta product, a handful of pairs, one laser: descriptive. Sign: −λ/4π·φ as for Sentinel-1 (MintPy convention), checked
against the laser, not assumed. Writes only to the hub (laser = field data).

    INSAR_DRIVE_ROOT=../local/drive_pristine PYTHONPATH=src python scripts/nisar_lband_x068.py
"""
from __future__ import annotations

import io
import json
import sys
import urllib.parse
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import h5py  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import requests  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import field_fusion_v2_x066 as V  # noqa: E402

from insar_wetlands.bootstrap import start  # noqa: E402
from insar_wetlands.stack import list_pairs, load_layer  # noqa: E402

OUT = V.DLV / "nisar_x068"
CMR = "https://cmr.earthdata.nasa.gov/search/granules.json"
COLLECTION = "C2850261892-ASF"                 # NISAR_L2_GUNW_BETA_V1
C_LIGHT = 299_792_458.0
GRID = "science/LSAR/GUNW/grids/frequencyA"


class HTTPRange(io.RawIOBase):
    """A read-only file over HTTP range requests, 1 MB blocks cached (Earthdata auth from ~/.netrc on redirect)."""

    def __init__(self, url: str, block: int = 1 << 20):
        self.s = requests.Session()
        r = self.s.get(url, headers={"Range": "bytes=0-0"}, allow_redirects=True, timeout=120)
        r.raise_for_status()
        self.url, self.size, self.pos, self.block, self.cache, self.requests = r.url, int(r.headers["Content-Range"].split("/")[-1]), 0, block, {}, 0

    def readable(self): return True
    def seekable(self): return True
    def tell(self): return self.pos

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
    q = urllib.parse.urlencode({"collection_concept_id": COLLECTION, "bounding_box": ",".join(map(str, bbox)), "page_size": 50})
    with urllib.request.urlopen(f"{CMR}?{q}", timeout=120) as r:
        feed = json.load(r)["feed"]["entry"]
    out = []
    for e in feed:
        h5 = [lk["href"] for lk in e.get("links", []) if lk.get("href", "").endswith(".h5")]
        if h5:
            out.append({"title": e["title"], "url": h5[0]})
    return sorted(out, key=lambda g: g["title"])


def window(coords: np.ndarray, lo: float, hi: float) -> slice:
    k = np.flatnonzero((coords >= lo) & (coords <= hi))
    return slice(int(k.min()), int(k.max()) + 1)


def read_pair(g: dict, tpl, zones, p6_xy) -> dict:
    f = HTTPRange(g["url"])
    h = h5py.File(io.BufferedReader(f, buffer_size=1 << 20), "r")
    lam = C_LIGHT / float(h[f"{GRID}/centerFrequency"][()])
    idn = h["science/LSAR/identification"]
    t1 = pd.Timestamp(idn["referenceZeroDopplerStartTime"][()].decode()[:19], tz="UTC")
    t2 = pd.Timestamp(idn["secondaryZeroDopplerStartTime"][()].decode()[:19], tz="UTC")
    xs, ys = tpl.x.values, tpl.y.values
    lo_x, hi_x, lo_y, hi_y = xs.min() - 100, xs.max() + 100, ys.min() - 100, ys.max() + 100
    out = {"granule": g["title"], "t1": t1, "t2": t2, "wavelength_m": lam}
    for grid, keys in (("unwrappedInterferogram", ("unwrappedPhase", "coherenceMagnitude", "ionospherePhaseScreen", "connectedComponents")),
                       ("wrappedInterferogram", ("coherenceMagnitude",))):
        gx, gy = h[f"{GRID}/{grid}/xCoordinates"][:], h[f"{GRID}/{grid}/yCoordinates"][:]
        sx, sy = window(gx, lo_x, hi_x), window(gy, lo_y, hi_y)
        cx, cy = np.meshgrid(gx[sx], gy[sy])
        # each product pixel → the atlas grid cell it falls in → its zone
        col = np.rint((cx - xs[0]) / (xs[1] - xs[0])).astype(int)
        row = np.rint((cy - ys[0]) / (ys[1] - ys[0])).astype(int)
        inside = (row >= 0) & (row < len(ys)) & (col >= 0) & (col < len(xs))
        zlab = np.full(cx.shape, "", dtype=object)
        for z in "ABCD":
            m = np.zeros(cx.shape, bool)
            m[inside] = zones[z][row[inside], col[inside]]
            zlab[m] = z
        k6 = np.unravel_index(np.argmin((cx - p6_xy[0]) ** 2 + (cy - p6_xy[1]) ** 2), cx.shape)
        for key in keys:
            a = h[f"{GRID}/{grid}/HH/{key}"][sy, sx].astype(float)
            tag = f"{'w20' if grid.startswith('wrapped') else 'u80'}_{key}"
            if key == "coherenceMagnitude":
                for z in "ABC":
                    out[f"{tag}_{z}"] = float(np.nanmedian(a[zlab == z]))
                out[f"{tag}_p6"] = float(a[k6])
            elif key == "connectedComponents":
                out["p6_component"] = int(a[k6])
                out["p6_component_same_as_grassland"] = bool(a[k6] != 0 and a[k6] in set(a[zlab == "C"].astype(int).tolist()))
            else:
                out[f"{key}_p6_minus_C"] = float(a[k6] - np.nanmedian(a[zlab == "C"]))
    rad = -lam / (4 * np.pi) * 1000.0
    out["los_mm"] = out["unwrappedPhase_p6_minus_C"] * rad
    out["los_mm_iono_corrected"] = (out["unwrappedPhase_p6_minus_C"] - out["ionospherePhaseScreen_p6_minus_C"]) * rad
    # incidence at P6 from the radar-grid cube (nearest height layer to the site, ~60 m)
    rg = h["science/LSAR/GUNW/metadata/radarGrid"]
    rx, ry, hz = rg["xCoordinates"][:], rg["yCoordinates"][:], rg["heightAboveEllipsoid"][:]
    out["incidence_deg"] = float(rg["incidenceAngle"][int(np.argmin(np.abs(hz - 60))), int(np.argmin(np.abs(ry - p6_xy[1]))),
                                                      int(np.argmin(np.abs(rx - p6_xy[0])))])
    out["http_requests"] = f.requests
    return out


def weather():
    """Open-Meteo hourly (X-066) and its calibration on the station (mean air-temperature bias for 'frozen at a time')."""
    om = pd.read_csv(V.HUB / "06_data" / "local_small" / "open_meteo" / "open_meteo_rzecin_hourly_2020_2025.csv", parse_dates=["time"]).set_index("time")
    om.index = pd.to_datetime(om.index, utc=True)
    meta = json.loads((V.DLV / "fusion_v2" / "flags_2025" / "open_meteo_vs_station.json").read_text())
    return om, meta


def frozen_at(om, meta, t) -> bool:
    """Air temperature (reanalysis minus its mean bias against the station) at or below 0 °C at time t."""
    a = om.Air_2m.iloc[om.index.get_indexer([pd.Timestamp(t)], method="nearest")[0]] - meta["air_bias_c"]
    return bool(a <= 0)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    ctx = start("nisar_lband_x068", mount=False, git=False)
    inp = V.Inputs(ctx, V.HUB / "05_code" / "local" / "s1_2020_2021")
    tpl, zones = ctx.template, inp.zones
    px = pd.read_csv(V.DLV / "field_first" / "plot_pixels.csv").set_index("plot")
    p6_xy = (float(px.loc["P6", "E"]), float(px.loc["P6", "N"]))
    from insar_wetlands.aoi import buffered_bbox
    cached = OUT / "nisar_pairs.csv"
    reuse = "--reuse" in sys.argv and cached.exists()
    gs = [] if reuse else granules(buffered_bbox(ctx.cfg))
    rows = [] if not reuse else pd.read_csv(cached, parse_dates=["t1", "t2"]).to_dict("records")
    om, om_meta = weather()
    for g in gs:
        r = read_pair(g, tpl, zones, p6_xy)
        la, lb = inp.laser_at(r["t1"]), inp.laser_at(r["t2"])
        across = inp.across_fill(r["t1"], r["t2"]) or inp.segment(r["t1"]) != inp.segment(r["t2"])
        r["laser_dh_mm"] = np.nan if across else lb - la
        r["laser_los_mm"] = r["laser_dh_mm"] * np.cos(np.radians(r["incidence_deg"]))
        rows.append(r)
        print(f"  {r['t1']:%Y-%m-%d}→{r['t2']:%Y-%m-%d}: coh mat {r['u80_coherenceMagnitude_A']:.2f} grass {r['u80_coherenceMagnitude_C']:.2f}"
              f" | P6 LOS {r['los_mm']:+.1f} mm (iono-corr {r['los_mm_iono_corrected']:+.1f}) | laser LOS {r['laser_los_mm']:+.1f} | {r['http_requests']} requests", flush=True)
    d = pd.DataFrame(rows)
    # surface state at 03:39 UTC (station-calibrated reanalysis, X-066 flags) and the laser without the snow screen
    d["frozen"] = [bool(frozen_at(om, om_meta, a) or frozen_at(om, om_meta, b)) for a, b in zip(d.t1, d.t2)]
    raw = inp.laser.surface_cm * 10
    near = lambda t: (float(raw.iloc[raw.index.get_indexer([t], method="nearest")[0]])  # noqa: E731
                      if abs((raw.index[raw.index.get_indexer([t], method="nearest")[0]] - t).total_seconds()) <= 3600 else np.nan)
    d["laser_dh_mm_unscreened_NOT_USED"] = [near(b) - near(a) for a, b in zip(d.t1, d.t2)]
    d.to_csv(OUT / "nisar_pairs.csv", index=False)

    # Sentinel-1 C-band, same weeks (X-066 2025 cache): coherence by zone and the P6 change against the laser
    lo, hi = pd.Timestamp(d.t1.min()).tz_convert(None), pd.Timestamp(d.t2.max()).tz_convert(None)
    cb = []
    fl25 = pd.read_csv(V.DLV / "fusion_v2" / "flags_2025" / "wetness_at_overpasses_2025_open_meteo.csv", parse_dates=["date"]).set_index(["date", "track"])
    for track in ("ascending", "descending"):
        root = V.HUB / "05_code" / "local" / "s1_2025" / f"hyp3_cropped_{track}"
        for pid in list_pairs(root):
            a, b = (pd.Timestamp(x) for x in pid.split("_"))
            if (b - a).days != 12 or a < lo.normalize() or b > hi.normalize():
                continue
            co = load_layer(root, "corr", [pid]).isel(pair=0).values
            fz = any(bool(fl25.frozen.get((x, track), True)) for x in (a, b))
            cb.append({"track": track, "t1": a, "t2": b, "dt": 12, "frozen": fz, **{f"coh_{z}": float(np.nanmedian(co[zones[z]])) for z in "ABC"},
                       "coh_p6": float(co[inp.p6]), "los_mm": np.nan, "laser_dh_mm": np.nan, "source": "2025 crops (12-day)"})
    for track in ("ascending", "descending"):
        cache = V.CACHE / f"R3_2025_6day_both_{track}.pkl"
        if not cache.exists():
            continue
        t, dl, co = pd.read_pickle(cache)
        for k, r in t.iterrows():
            if r.t1 >= lo.normalize() and r.t2 <= hi.normalize():
                cb.append({"track": track, "t1": r.t1, "t2": r.t2, "dt": int(r["dt"]), "frozen": bool(r.frozen), "wet": bool(r.wet), "source": "X-066 cache (consecutive)",
                           **{f"coh_{z}": float(np.nanmedian(co[k][zones[z]])) for z in "ABC"}, "coh_p6": float(r.coh_p6),
                           "los_mm": float(r.s1_p6), "laser_dh_mm": float(r.dh_laser) if np.isfinite(r.dh_laser) else np.nan})
    c = pd.DataFrame(cb)
    c.to_csv(OUT / "sentinel1_same_weeks.csv", index=False)

    lz = d.dropna(subset=["laser_los_mm"])
    summ = {"pairs": len(d), "pairs_with_laser": len(lz),
            "lband_coherence_median": {z: float(d[f"u80_coherenceMagnitude_{z}"].median()) for z in "ABC"},
            "lband_coherence20m_median": {z: float(d[f"w20_coherenceMagnitude_{z}"].median()) for z in "ABC"},
            "lband_12d_by_state": {("frozen" if fz else "not frozen"): {"pairs": len(g), **{z: float(g[f"u80_coherenceMagnitude_{z}"].median()) for z in "ABC"}}
                                   for fz, g in d.groupby("frozen")},
            "cband_by_track_dt_state": {f"{tr} {dt}d {'frozen' if fz else 'not frozen'}": {"pairs": len(g), **{z: float(g[f"coh_{z}"].median()) for z in "ABC"}}
                                        for (tr, dt, fz), g in c.groupby(["track", "dt", "frozen"])} if len(c) else {},
            "p6_same_sign_as_laser": int((np.sign(lz.los_mm_iono_corrected) == np.sign(lz.laser_los_mm)).sum()) if len(lz) else 0,
            "p6_r_with_laser": float(np.corrcoef(lz.los_mm_iono_corrected, lz.laser_los_mm)[0, 1]) if len(lz) >= 3 else None}
    (OUT / "summary.json").write_text(json.dumps(summ, indent=2))

    fig, ax = plt.subplots(1, 2, figsize=(13, 4.2))
    mid = d.t1 + (d.t2 - d.t1) / 2
    for z, col in (("A", "C3"), ("C", "C2"), ("B", "C0")):
        ax[0].plot(mid, d[f"u80_coherenceMagnitude_{z}"], "o-", color=col, label=f"L-band 12 d, zone {z}")
    if len(c):
        for (tr, dt), g in c[c.track == "ascending"].groupby(["track", "dt"]):
            ax[0].plot(g.t1 + (g.t2 - g.t1) / 2, g.coh_A, "x:" if dt == 6 else "s--", color="C3", alpha=0.45 if dt == 6 else 0.8,
                       mfc="none", label=f"C-band {dt} d asc, mat")
    ax[0].set_ylim(0, 1); ax[0].set_ylabel("median coherence"); ax[0].legend(fontsize=7); ax[0].set_title("Coherence: NISAR L-band vs Sentinel-1 C-band")
    ax[1].bar(np.arange(len(d)) - 0.2, d.los_mm_iono_corrected, 0.4, label="NISAR P6 − grassland (iono-corrected)")
    ax[1].bar(np.arange(len(d)) + 0.2, d.laser_los_mm, 0.4, label="laser (LOS)")
    ax[1].set_xticks(range(len(d)), [f"{a:%m-%d}" for a in d.t1], fontsize=7); ax[1].axhline(0, c="k", lw=0.5)
    ax[1].set_ylabel("mm, + toward the satellite"); ax[1].legend(fontsize=7); ax[1].set_title("P6, 12-day change")
    fig.tight_layout(); fig.savefig(OUT / "fig_nisar.png", dpi=130); plt.close(fig)
    readme(d, c, summ)
    print(json.dumps(summ, indent=2))


def readme(d: pd.DataFrame, c: pd.DataFrame, summ: dict) -> None:
    md = lambda df: "\n".join(["| " + " | ".join(map(str, df.columns)) + " |", "|" + "---|" * len(df.columns)]  # noqa: E731
                              + ["| " + " | ".join(f"{v:.2f}" if isinstance(v, float) else str(v) for v in r) + " |" for r in df.itertuples(index=False)])
    pairs = d.assign(pair=[f"{pd.Timestamp(a):%Y-%m-%d} → {pd.Timestamp(b):%Y-%m-%d}" for a, b in zip(d.t1, d.t2)])[
        ["pair", "frozen", "u80_coherenceMagnitude_A", "u80_coherenceMagnitude_C", "u80_coherenceMagnitude_B", "w20_coherenceMagnitude_A",
         "w20_coherenceMagnitude_C", "los_mm", "los_mm_iono_corrected", "laser_los_mm", "laser_dh_mm_unscreened_NOT_USED"]].rename(columns={
        "u80_coherenceMagnitude_A": "coh mat (80 m)", "u80_coherenceMagnitude_C": "coh grass (80 m)", "u80_coherenceMagnitude_B": "coh lake (80 m)",
        "w20_coherenceMagnitude_A": "coh mat (20 m)", "w20_coherenceMagnitude_C": "coh grass (20 m)", "los_mm": "P6 − grassland, LOS mm",
        "los_mm_iono_corrected": "same, iono-corrected", "laser_los_mm": "laser LOS mm (screened)", "laser_dh_mm_unscreened_NOT_USED": "laser Δh mm UNSCREENED (not used)"})
    cmp = pd.DataFrame([{"sensor": "NISAR L-band 12 d", "state": k, "pairs": v["pairs"], "mat": v["A"], "grassland": v["C"], "lake": v["B"]}
                        for k, v in summ["lband_12d_by_state"].items()]
                       + [{"sensor": f"Sentinel-1 C-band {k.replace(' not frozen', '').removesuffix(' frozen')}",
                           "state": "not frozen" if "not frozen" in k else "frozen", "pairs": v["pairs"], "mat": v["A"], "grassland": v["C"], "lake": v["B"]}
                          for k, v in summ["cband_by_track_dt_state"].items()])
    L = ["# X-068 — NISAR L-band over Rzecin, Oct 2025 – Jan 2026 (exploratory)", "",
         "Generated by `05_code/SAr-adapted/scripts/nisar_lband_x068.py`; every number from its tables. NISAR beta GUNW, track 158 "
         "ascending (03:39 UTC), consecutive 12-day pairs, read by HTTP range requests (only the AOI window of each ~2.5 GB file). "
         "X-035 (2026-09-13) found no product over the site; the archive filled in afterwards.", "",
         "## Coherence by zone: L-band vs C-band over the same weeks, split by surface state", "", md(cmp), "",
         "Frozen = air temperature ≤ 0 °C at either acquisition (station-calibrated reanalysis for NISAR's 03:39 UTC; the 2025 overpass "
         "flags for Sentinel-1).", "",
         "## Per NISAR pair", "", md(pairs), "",
         "## Reading", "",
         "- In this late-autumn/winter window L-band is **not** more coherent than C-band on the mat; C-band 12-day coherence is "
         "already high (vegetation senescent). The window cannot test L-band where C-band fails — summer.",
         "- No P6 motion test: every laser value at the 03:39 UTC overpasses falls under the snow screen (any sub-zero hour in the "
         "72 h before); the unscreened change is listed for transparency only and is not used.",
         "- Beta product; its coherence is estimated with different multilooking than HyP3's, so absolute levels compare loosely.",
         "- The real test is NISAR summer 2026 (tracks 057 A, 158 A, 008 D exist, provisional) against Sentinel-1 2026 — without laser or "
         "water table unless the field team has 2026 data.", "",
         "## Files", "", "| File | Content |", "|---|---|",
         "| `nisar_pairs.csv` | every NISAR pair: coherence by zone (80 m, 20 m), P6 change, ionosphere screen, incidence, laser |",
         "| `sentinel1_same_weeks.csv` | Sentinel-1 pairs over the same weeks (12-day from the 2025 crops; consecutive from the X-066 cache) |",
         "| `summary.json`, `fig_nisar.png` | the summary and its figure |"]
    (OUT / "README.md").write_text("\n".join(L) + "\n")


if __name__ == "__main__":
    main()

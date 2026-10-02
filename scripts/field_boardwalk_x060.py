"""X-060 — the boardwalk: where it is on the radar grids, whether it shows, and what it does to P6.

The supervisor (2026-09-28): the plots lie along a boardwalk; use the finest reliable resolution,
the smallest extraction area with enough valid pixels, and flag or exclude pixels the boardwalk
influences. Three steps:

1. **Outline.** The national orthophoto (GUGiK Geoportal WMS, ~0.25 m, public) around the plot
   transect, reprojected to the analysis CRS. Wood and metal read as light grey against the
   local background; a centreline is traced through those pixels bin by bin along the transect
   (dynamic programming with a continuity limit), buffered to the walkway width, plus the
   structures detected inside the plot frames (platforms, instrument housings). Its share of
   every 40 m, 20 m and 10 m cell is computed.
2. **Does it show?** Each boardwalk cell against the mat cells 2–4 cells away in the same row,
   for coherence and phase scatter (40 m, consecutive pairs 2022–2024) and backscatter (VV, VH:
   40 m 2020–2024; 10 m summer 2022–2024, ``rtc_10m_window.py``); the same lateral anomaly for
   every other mat cell is the null.
3. **P6 extraction.** X-059's per-pair statistics for every window from 1×1 to 3×3, with and
   without the boardwalk cells — whether any conclusion depends on the window or on the boardwalk.

Field data (plot coordinates, laser, WTD) are read from the hub at run time; outputs go only to
the hub. No coordinate is written in this file.

    INSAR_DRIVE_ROOT=../local/drive_pristine PYTHONPATH=src python scripts/field_boardwalk_x060.py
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import xarray as xr  # noqa: E402
from scipy import ndimage as ndi  # noqa: E402
from scipy.stats import mannwhitneyu  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import field_p6_validation as V  # noqa: E402

from insar_wetlands.bootstrap import start  # noqa: E402
from insar_wetlands.inversion.isbas import PHASE_TO_MM  # noqa: E402
from insar_wetlands.stack import load_layer  # noqa: E402

HUB = V.HUB
DLV = V.DLV
WMS = "https://mapy.geoportal.gov.pl/wss/service/PZGIK/ORTO/WMS/HighResolution"
RES = 0.25                 # orthophoto pixel (m)
MARGIN = 120.0             # around the plot transect (m)
WALK_HALF_WIDTH = 0.6      # m — the walkway is about 1.2 m wide on the orthophoto
BIN = 4.0                  # m along the transect for the centreline
MAX_SHIFT = 2.0            # m of sideways shift allowed between bins
SHIFT_COST = 4.0           # light pixels given up per 0.5 m of sideways shift (keeps the line straight unless the data insist)
STRUCT_RADIUS = 20.0       # m — structures kept around each plot
LATERAL = (2, 3, 4)        # cells on each side used as the local reference


# ------------------------------------------------------------------------------ 1. outline

def orthophoto(out: Path, px: pd.DataFrame) -> tuple[np.ndarray, tuple]:
    """The orthophoto around the transect on a 0.25 m grid of the analysis CRS (cached GeoTIFF)."""
    import rasterio
    from rasterio.transform import from_bounds
    from rasterio.warp import Resampling, reproject

    E0, E1 = np.floor(px.E.min() - MARGIN), np.ceil(px.E.max() + MARGIN)
    N0, N1 = np.floor(px.N.min() - MARGIN), np.ceil(px.N.max() + MARGIN)
    tif = out / "orthophoto_transect.tif"
    if not tif.exists():
        import io
        import urllib.request

        from PIL import Image
        from pyproj import Transformer

        t = Transformer.from_crs(32633, 2180, always_xy=True)
        corners = [t.transform(x, y) for x in (E0, E1) for y in (N0, N1)]
        bb = (min(c[0] for c in corners), min(c[1] for c in corners), max(c[0] for c in corners), max(c[1] for c in corners))
        w, h = int((bb[2] - bb[0]) / RES), int((bb[3] - bb[1]) / RES)
        url = (f"{WMS}?SERVICE=WMS&VERSION=1.1.1&REQUEST=GetMap&LAYERS=Raster&STYLES=&SRS=EPSG:2180"
               f"&BBOX={bb[0]},{bb[1]},{bb[2]},{bb[3]}&WIDTH={w}&HEIGHT={h}&FORMAT=image/png")
        img = np.array(Image.open(io.BytesIO(urllib.request.urlopen(url, timeout=120).read())).convert("RGB"))
        W, H = int((E1 - E0) / RES), int((N1 - N0) / RES)
        dst = np.zeros((3, H, W), np.uint8)
        for b in range(3):
            reproject(img[..., b], dst[b], src_transform=from_bounds(*bb, w, h), src_crs="EPSG:2180",
                      dst_transform=from_bounds(E0, N0, E1, N1, W, H), dst_crs="EPSG:32633",
                      resampling=Resampling.bilinear)
        with rasterio.open(tif, "w", driver="GTiff", width=W, height=H, count=3, dtype="uint8", crs="EPSG:32633",
                           transform=from_bounds(E0, N0, E1, N1, W, H), compress="jpeg", photometric="ycbcr",
                           tiled=True) as f:
            f.write(dst)
            f.update_tags(source=f"GUGiK Geoportal orthophoto WMS ({WMS}), layer Raster, EPSG:2180 → 32633 bilinear")
    with rasterio.open(tif) as f:
        img = f.read()
    return img.transpose(1, 2, 0).astype(float), (E0, E1, N0, N1)


def light_pixels(img: np.ndarray) -> np.ndarray:
    """Light-grey pixels (wood, metal) standing out of the local background: brightness above a
    5 m median by > 18, low colour saturation, absolute brightness > 110; opened, small blobs out."""
    br = img.mean(2)
    sat = img.max(2) - img.min(2)
    m = ((br - ndi.median_filter(br, size=21)) > 18) & (sat < 45) & (br > 110)
    m = ndi.binary_opening(m, np.ones((2, 2)))
    lab, n = ndi.label(m)
    size = ndi.sum(m, lab, range(1, n + 1))
    return np.isin(lab, 1 + np.flatnonzero(size >= 40))


def polyline_distance(E, N, pts) -> np.ndarray:
    d = np.full(E.shape, np.inf)
    for a, b in zip(pts[:-1], pts[1:]):
        ab = b - a
        t = np.clip(((E - a[0]) * ab[0] + (N - a[1]) * ab[1]) / (ab @ ab), 0, 1)
        d = np.minimum(d, np.hypot(E - a[0] - t * ab[0], N - a[1] - t * ab[1]))
    return d


def centreline(mask, E, N, guide, half_corridor=25.0) -> pd.DataFrame:
    """The walkway's centreline, one point per BIN m of northing: the easting whose 1.5 m-wide strip
    holds the most light pixels, chosen jointly over all bins (dynamic programming) with at most
    MAX_SHIFT m of sideways move between bins, each move costing SHIFT_COST, inside a corridor around the plot line. Bins
    without support are marked; the line is trimmed to the supported stretch."""
    n_lo, n_hi = N.min(), N.max()
    edges = np.arange(n_lo, n_hi + BIN, BIN)
    gE = np.interp((edges[:-1] + edges[1:]) / 2, guide[:, 1], guide[:, 0])
    cand = np.arange(-half_corridor, half_corridor + 0.5, 0.5)
    rows_E = E[0]
    score = np.zeros((len(edges) - 1, cand.size))
    col_n = N[:, 0]
    for b in range(len(edges) - 1):
        rs = (col_n >= edges[b]) & (col_n < edges[b + 1])
        hit = mask[rs].sum(0)
        for k, c in enumerate(cand):
            score[b, k] = hit[np.abs(rows_E - (gE[b] + c)) <= 0.75].sum()
    step = int(MAX_SHIFT / 0.5)
    acc, back = score.copy(), np.zeros(score.shape, int)
    for b in range(1, len(score)):
        for k in range(cand.size):
            lo, hi = max(0, k - step), min(cand.size, k + step + 1)
            prev = acc[b - 1, lo:hi] - SHIFT_COST * np.abs(np.arange(lo, hi) - k)
            j = lo + int(np.argmax(prev))
            acc[b, k] += prev[j - lo]
            back[b, k] = j
    path = [int(np.argmax(acc[-1]))]
    for b in range(len(score) - 1, 0, -1):
        path.append(back[b, path[-1]])
    path = path[::-1]
    cl = pd.DataFrame({"N": (edges[:-1] + edges[1:]) / 2, "E": [gE[b] + cand[k] for b, k in enumerate(path)],
                       "support_px": [score[b, k] for b, k in enumerate(path)]})
    sup = cl.support_px >= 8
    # gaps up to 3 bins inside the supported stretch are walkway hidden by vegetation; trim the ends
    first, last = np.flatnonzero(sup)[[0, -1]]
    cl = cl.iloc[first:last + 1].reset_index(drop=True)
    run = (~(cl.support_px >= 8)).astype(int)
    gap = run.groupby((run == 0).cumsum()).transform("sum")
    cl["supported"] = (cl.support_px >= 8) | (gap <= 3)
    return cl


def outline(img, ext, px) -> tuple[np.ndarray, pd.DataFrame, pd.DataFrame]:
    E0, E1, N0, N1 = ext
    H, W, _ = img.shape
    E = E0 + (np.arange(W) + 0.5) * RES
    N = N1 - (np.arange(H) + 0.5) * RES
    EE, NN = np.meshgrid(E, N)
    light = light_pixels(img)
    guide = px.sort_values("N")[["E", "N"]].to_numpy()
    guide = np.vstack([[guide[0, 0], N0], guide, [guide[-1, 0], N1]])
    near = polyline_distance(EE, NN, guide) < 30
    cl = centreline(light & near, EE, NN, guide)
    keep = cl[cl.supported]
    pts = keep[["E", "N"]].to_numpy()
    walk = np.zeros(light.shape, bool)
    # one strip per run of supported bins (a long unsupported gap is not drawn)
    for _, g in keep.groupby((~cl.supported).cumsum().loc[keep.index]):
        if len(g) >= 2:
            walk |= polyline_distance(EE, NN, g[["E", "N"]].to_numpy()) <= WALK_HALF_WIDTH
    lab, n = ndi.label(light & near & ~ndi.binary_dilation(walk, iterations=8))
    comps = []
    for i, sl in enumerate(ndi.find_objects(lab), start=1):
        m = lab[sl] == i
        ce, cn = EE[sl][m].mean(), NN[sl][m].mean()
        dplot = np.hypot(px.E - ce, px.N - cn)
        comps.append({"id": i, "E": ce, "N": cn, "area_m2": m.sum() * RES ** 2,
                      "nearest_plot": px["plot"].iloc[int(np.argmin(dplot))], "dist_plot_m": float(dplot.min()),
                      "kept": bool(dplot.min() <= STRUCT_RADIUS)})
    comps = pd.DataFrame(comps)
    kept_ids = comps.loc[comps.kept, "id"].to_numpy() if len(comps) else []
    struct = ndi.binary_dilation(np.isin(lab, kept_ids), iterations=1)
    bw = walk | struct
    cl["length_m"] = BIN
    _ = pts
    return np.stack([bw, walk, struct]), cl, comps


def cell_fraction(mask: np.ndarray, ext, x: np.ndarray, y: np.ndarray, res: float) -> np.ndarray:
    """Share of each grid cell (centres x, y; size res) covered by the 0.25 m mask (0 outside it)."""
    E0, E1, N0, N1 = ext
    H, W = mask.shape
    f = np.zeros((y.size, x.size))
    k = int(round(res / RES))
    for i, yc in enumerate(y):
        r0 = int(round((N1 - (yc + res / 2)) / RES))
        if r0 + k <= 0 or r0 >= H:
            continue
        for j, xc in enumerate(x):
            c0 = int(round((xc - res / 2 - E0) / RES))
            if c0 + k <= 0 or c0 >= W:
                continue
            blk = mask[max(r0, 0):min(r0 + k, H), max(c0, 0):min(c0 + k, W)]
            f[i, j] = blk.sum() / (k * k)
    return f


# ------------------------------------------------------------------------------ 2. does it show?

def lateral_anomaly(val: np.ndarray, mat: np.ndarray) -> np.ndarray:
    """Value minus the median of the mat cells LATERAL cells to each side in the same row."""
    out = np.full(val.shape, np.nan)
    H, W = val.shape
    for r in range(H):
        for c in range(W):
            if not mat[r, c] or not np.isfinite(val[r, c]):
                continue
            ref = [val[r, c + s * d] for d in LATERAL for s in (-1, 1)
                   if 0 <= c + s * d < W and mat[r, c + s * d] and np.isfinite(val[r, c + s * d])]
            if len(ref) >= 3:
                out[r, c] = val[r, c] - np.median(ref)
    return out


def test_cells(name: str, track: str, res_m: int, val: np.ndarray, mat: np.ndarray, frac: np.ndarray,
               labels: dict) -> tuple[list[dict], list[dict]]:
    """Boardwalk cells' lateral anomalies against every other mat cell's (the null)."""
    an = lateral_anomaly(val, mat)
    bw = (frac > 0) & mat & np.isfinite(an)
    null = an[~(ndi.binary_dilation(frac > 0, iterations=1)) & mat & np.isfinite(an)]
    cells = []
    for r, c in zip(*np.nonzero(bw)):
        cells.append({"variable": name, "track": track, "res_m": res_m, "row": int(r), "col": int(c),
                      "boardwalk_share": float(frac[r, c]), "plot": labels.get((int(r), int(c)), ""),
                      "value": float(val[r, c]), "lateral_anomaly": float(an[r, c]),
                      "null_percentile": float((null < an[r, c]).mean() * 100)})
    a = an[bw]
    p = mannwhitneyu(a, null).pvalue if len(a) >= 3 else np.nan
    summ = [{"variable": name, "track": track, "res_m": res_m, "n_boardwalk_cells": int(bw.sum()),
             "n_null_cells": int(null.size), "median_anomaly_boardwalk": float(np.median(a)) if len(a) else np.nan,
             "median_anomaly_null": float(np.median(null)), "null_q05": float(np.percentile(null, 5)),
             "null_q95": float(np.percentile(null, 95)), "share_beyond_null_90": float(np.mean(
                 (a < np.percentile(null, 5)) | (a > np.percentile(null, 95)))) if len(a) else np.nan,
             "p_mannwhitney": p}]
    return cells, summ


def radar_layers(ctx, df: pd.DataFrame, rtc_root: Path) -> dict:
    """Per-pixel 40 m summaries, both tracks: mean coherence and phase SD of the consecutive pairs
    (2022–2024), VV and VH mean (dB of mean power, 2020–2024)."""
    zones = {k: ctx.zones[k].values.astype(bool) for k in "ABCD"}
    crop = {"ascending": ctx.paths.cropped,
            "descending": ctx.paths.for_phase("field_boardwalk_x060", track="descending").cropped}
    names = {"ascending": "rtc_dualpol_2020_2024.nc", "descending": "rtc_dualpol_2020_2024_descending.nc"}
    lay = {}
    for track in crop:
        d = df[(df.window == V.PRIMARY) & (df.track == track) & df.in_network]
        corr = load_layer(crop[track], "corr", list(d.pair)).values
        unw = load_layer(crop[track], "unw_phase", list(d.pair)).values
        refC = np.array([np.nanmedian(u[zones["C"] & np.isfinite(u)]) for u in unw])
        lay[("coherence, consecutive pairs (mean)", track)] = np.nanmean(corr, 0)
        lay[("phase vs grassland, SD over consecutive pairs (mm)", track)] = np.nanstd(
            (unw - refC[:, None, None]) * PHASE_TO_MM, 0)
        rtc = xr.open_dataset(rtc_root / names[track])
        for v, lab in (("gamma0_vv_db", "VV mean (dB)"), ("gamma0_vh_db", "VH mean (dB)")):
            g = ctx.to_grid(rtc[v]).values
            lay[(lab, track)] = 10 * np.log10(np.nanmean(10 ** (g / 10), 0))
    return lay


def rtc10_layers(rtc_root: Path, x10, y10) -> dict:
    lay = {}
    for track in ("ascending", "descending"):
        ds = xr.open_dataset(rtc_root / f"rtc10m_summer_2022_2024_{track}.nc").sel(x=x10, y=y10, method="nearest")
        for v, lab in (("gamma0_vv_db", "VV mean (dB)"), ("gamma0_vh_db", "VH mean (dB)")):
            lin = 10 ** (ds[v].values / 10)
            lay[(lab, track)] = 10 * np.log10(np.nanmean(lin, 0))
            lay[(lab.replace("mean", "temporal SD"), track)] = np.nanstd(ds[v].values, 0)
    return lay


# ------------------------------------------------------------------------------ 3. P6 windows

def p6_windows(r0: int, c0: int, flagged: np.ndarray) -> dict[str, list[tuple[int, int]]]:
    """Every extraction window around P6's cell, from 1×1 to 3×3, with and without boardwalk cells."""
    w = {"1x1 (P6 cell)": [(r0, c0)],
         "2x1 north": [(r0, c0), (r0 - 1, c0)], "2x1 south": [(r0, c0), (r0 + 1, c0)],
         "1x2 west": [(r0, c0), (r0, c0 - 1)], "1x2 east": [(r0, c0), (r0, c0 + 1)]}
    for dr, dc, lab in ((-1, -1, "NW"), (-1, 0, "NE"), (0, -1, "SW"), (0, 0, "SE")):
        w[f"2x2 {lab}"] = [(r0 + dr + i, c0 + dc + j) for i in (0, 1) for j in (0, 1)]
    w["3x3"] = [(r0 + i, c0 + j) for i in (-1, 0, 1) for j in (-1, 0, 1)]
    w["3x3 without boardwalk cells"] = [rc for rc in w["3x3"] if not flagged[rc]]
    w["west neighbour only"] = [(r0, c0 - 1)]
    w["east neighbour only"] = [(r0, c0 + 1)]
    return {k: v for k, v in w.items() if v}


def window_stats_table(ctx, df, floor_mm, r0, c0, flagged, mat) -> pd.DataFrame:
    zones = {k: ctx.zones[k].values.astype(bool) for k in "ABCD"}
    crop = {"ascending": ctx.paths.cropped,
            "descending": ctx.paths.for_phase("field_boardwalk_x060", track="descending").cropped}
    rows = []
    wins = p6_windows(r0, c0, flagged)
    for track in crop:
        base = df[(df.window == V.PRIMARY) & (df.track == track)].sort_values("mid").reset_index(drop=True)
        d = base[base.in_network]
        unw = load_layer(crop[track], "unw_phase", list(d.pair)).values
        corr = load_layer(crop[track], "corr", list(d.pair)).values
        refC = np.array([np.nanmedian(u[zones["C"] & np.isfinite(u)]) for u in unw])
        for name, cells in wins.items():
            cells = [rc for rc in cells if mat[rc]]
            rr, cc = np.array([c[0] for c in cells]), np.array([c[1] for c in cells])
            s1 = (np.nanmedian(unw[:, rr, cc], 1) - refC) * PHASE_TO_MM
            coh = np.nanmedian(corr[:, rr, cc], 1)
            w = base.copy()
            w.loc[w.in_network, "s1_dlos_mm"] = s1
            w.loc[w.in_network, "coh"] = coh
            for subset, sel in (("all", np.ones(len(w), bool)), ("dry", (~w.wet_either & ~w.frozen_either).to_numpy())):
                st = V.pair_stats(w[sel], floor_mm)
                rows.append({"track": track, "window": name, "subset": subset, "n_cells": len(cells),
                             "boardwalk_cells": int(sum(flagged[rc] for rc in cells)), **st})
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------------ figures and README

def fig_outline(out, img, ext, masks, x40, y40, frac40, px):
    E0, E1, N0, N1 = ext
    fig, ax = plt.subplots(1, 2, figsize=(15, 9))
    for a in ax:
        a.imshow(img.astype(np.uint8), extent=(E0, E1, N0, N1))
        a.imshow(np.ma.masked_where(~masks[0], masks[0]), extent=(E0, E1, N0, N1), cmap="autumn", alpha=0.8, interpolation="none")
        a.scatter(px.E, px.N, s=14, c="cyan", zorder=3)
        for r in px.itertuples():
            a.annotate(r.plot, (r.E, r.N), xytext=(5, 3), textcoords="offset points", color="cyan", fontsize=8)
    for xe in np.r_[x40 - 20, x40[-1] + 20]:
        ax[0].axvline(xe, color="w", lw=0.4, alpha=0.6)
    for yn in np.r_[y40 + 20, y40[-1] - 20]:
        ax[0].axhline(yn, color="w", lw=0.4, alpha=0.6)
    for i, yc in enumerate(y40):
        for j, xc in enumerate(x40):
            if frac40[i, j] > 0:
                ax[0].add_patch(plt.Rectangle((xc - 20, yc - 20), 40, 40, fill=False, ec="yellow", lw=1.2))
                ax[0].text(xc, yc - 17, f"{frac40[i, j] * 100:.1f}%", color="yellow", fontsize=6, ha="center")
    ax[0].set_xlim(E0, E1)
    ax[0].set_ylim(N0, N1)
    ax[0].set_title("Boardwalk traced on the orthophoto (red); 40 m cells it touches (yellow, % covered)", fontsize=9)
    sel = px.set_index("plot").loc[["P5", "P6", "P8"]]
    ax[1].set_xlim(sel.E.min() - 35, sel.E.max() + 35)
    ax[1].set_ylim(sel.N.min() - 25, sel.N.max() + 25)
    ax[1].set_title("Zoom P5–P8: walkway and platforms", fontsize=9)
    fig.tight_layout()
    fig.savefig(out / "fig_boardwalk_outline.png", dpi=110)
    plt.close(fig)


def fig_tests(out, cells: pd.DataFrame):
    v = cells.groupby(["variable", "res_m"]).size().index.tolist()
    fig, axes = plt.subplots(1, len(v), figsize=(3.2 * len(v), 4), sharey=True)
    for a, (name, res) in zip(np.atleast_1d(axes), v):
        g = cells[(cells.variable == name) & (cells.res_m == res)]
        for k, (t, col) in enumerate((("ascending", "tab:orange"), ("descending", "tab:blue"))):
            gg = g[g.track == t]
            a.scatter(np.full(len(gg), k) + np.random.default_rng(0).uniform(-0.15, 0.15, len(gg)), gg.null_percentile,
                      c=col, s=12)
        a.axhspan(5, 95, color="0.9", zorder=0)
        a.set_xticks([0, 1], ["asc", "desc"])
        a.set_title(f"{name}\n{res} m", fontsize=7)
    np.atleast_1d(axes)[0].set_ylabel("boardwalk cell's lateral anomaly, percentile among mat cells")
    fig.tight_layout()
    fig.savefig(out / "fig_boardwalk_tests.png", dpi=110)
    plt.close(fig)


def md(df: pd.DataFrame, fmt="{:.2f}") -> str:
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for r in df.itertuples(index=False):
        lines.append("| " + " | ".join(fmt.format(v) if isinstance(v, float) else str(v) for v in r) + " |")
    return "\n".join(lines)


def readme(out, cl, comps, cells40, frac_tab, summ, win, p6_cell):
    s = summ.copy()
    s["p"] = s.p_mannwhitney.map(V.fp)
    s = s[["variable", "track", "res_m", "n_boardwalk_cells", "median_anomaly_boardwalk", "median_anomaly_null",
           "null_q05", "null_q95", "share_beyond_null_90", "p"]]
    w = win[win.subset == "all"][["track", "window", "n_cells", "boardwalk_cells", "n_laser", "r_s1_vs_laser", "slope_tls",
                                    "ratio_median", "n_same_sign", "n_above_floor", "r_s1_vs_dwtd", "r_coh_vs_abs_dwtd"]].copy()
    w["same sign"] = [f"{a}/{b}" for a, b in zip(w.n_same_sign, w.n_above_floor)]
    w = w.drop(columns=["n_same_sign", "n_above_floor"])
    p6c = cells40[cells40["plot"].str.contains("P6") & (cells40.res_m == 40)]
    ref = pd.read_csv(DLV / "field_p6" / "p6_stats.csv")
    ref = ref[(ref.window == V.PRIMARY) & (ref.subset == "all")].set_index("track")
    one = win[(win.window == "1x1 (P6 cell)") & (win.subset == "all")].set_index("track")
    chk = max(abs(one.loc[t, k] - ref.loc[t, k]) for t in one.index for k in ("r_s1_vs_laser", "slope_tls", "ratio_median"))
    find = []
    for t in ("ascending", "descending"):
        c = summ[(summ.track == t) & (summ.res_m == 40) & summ.variable.str.startswith("coherence")].iloc[0]
        find.append(f"- **Coherence, {t}, 40 m**: boardwalk cells {c.median_anomaly_boardwalk:+.3f} against their row neighbours "
                    f"(other mat cells {c.median_anomaly_null:+.3f}; p {V.fp(c.p_mannwhitney)}); "
                    f"{c.share_beyond_null_90 * 100:.0f} % of them outside the null's 5–95 % range.")
    sig = summ[summ.p_mannwhitney < 0.01]
    find.append("- **Every test with p < 0.01**: " + "; ".join(
        f"{r.variable}, {r.track}, {r.res_m} m ({r.median_anomaly_boardwalk:+.3f} vs {r.median_anomaly_null:+.3f})"
        for r in sig.itertuples()) + ".")
    for t in ("ascending", "descending"):
        g = win[(win.track == t) & (win.subset == "all")]
        find.append(f"- **P6, {t}, {len(g)} windows**: r with the laser {g.r_s1_vs_laser.min():.2f}–{g.r_s1_vs_laser.max():.2f}, "
                    f"TLS slope {g.slope_tls.min():.2f}–{g.slope_tls.max():.2f}, median per-pair ratio "
                    f"{g.ratio_median.min():.2f}–{g.ratio_median.max():.2f}, same sign {int(g.n_same_sign.min())}–"
                    f"{int(g.n_same_sign.max())} of {int(g.n_above_floor.max())}; coherence vs |ΔWTD| "
                    f"{g.r_coh_vs_abs_dwtd.min():.2f} to {g.r_coh_vs_abs_dwtd.max():.2f}.")
    L = ["# X-060 — the boardwalk on the radar grids, and P6's extraction window (exploratory)", "",
         "Generated by `05_code/SAr-adapted/scripts/field_boardwalk_x060.py` (outputs here only; plot positions are field data). "
         "The supervisor's request of 2026-09-28 (`01_admin/supervisor/2026-09-28_email_p6_priority.md`): finest reliable "
         "resolution, smallest extraction area, flag pixels the boardwalk influences. Every number is computed by the script.", "",
         "## Findings", "", *find, "",
         f"Check: the 1×1 window reproduces X-059 (`field_p6/p6_stats.csv`) to {chk:.1e}.", "",
         "## 1. Where the boardwalk is", "",
         f"Traced on the national orthophoto (GUGiK, ~{RES} m; `orthophoto_transect.tif`): light-grey pixels against the local "
         f"background, a centreline chosen jointly over {BIN:.0f} m bins (sideways shift ≤ {MAX_SHIFT:.0f} m per bin) buffered to "
         f"±{WALK_HALF_WIDTH} m, plus structures within {STRUCT_RADIUS:.0f} m of a plot (platforms, instrument housings). "
         f"Centreline: {int(cl.supported.sum())} of {len(cl)} bins supported (`boardwalk_centreline.csv`); "
         f"structures kept: {int(comps.kept.sum()) if len(comps) else 0} of {len(comps)} light components (`boardwalk_components.csv`). "
         "Check by eye in `fig_boardwalk_outline.png` — an automatic trace, not a survey.", "",
         "Share of each plot's cell covered by boardwalk, by grid size (`boardwalk_cells.csv` has every cell):", "",
         md(frac_tab, "{:.3f}"), "",
         f"P6's 40 m cell (row {p6_cell[0]}, col {p6_cell[1]}) contains the whole P6 platform. Any 40 m extraction at a plot "
         "contains boardwalk; the question is whether it changes the signal.", "",
         "## 2. Does the boardwalk show in the radar data?", "",
         f"Each boardwalk cell on the mat against the mat cells {LATERAL[0]}–{LATERAL[-1]} cells to each side in the same row "
         "(lateral anomaly); the null is the same anomaly for every mat cell away from the boardwalk. A boardwalk effect = "
         "boardwalk cells outside the null's 5–95 % range more often than 10 % and a small p.", "",
         md(s, "{:.3f}"), "",
         "`share_beyond_null_90`: share of boardwalk cells outside the null's 5–95 % range (10 % expected by chance). "
         "Per cell: `boardwalk_cell_tests.csv`.", "",
         "P6's own 40 m cell:", "",
         md(p6c[["variable", "track", "value", "lateral_anomaly", "null_percentile"]], "{:.3f}"), "",
         "## 3. P6: does any X-059 result depend on the window or on the boardwalk?", "",
         "X-059's statistics on consecutive pairs 2022–2024 for every window around P6 (median of the window's mat pixels "
         "minus the grassland median), all pairs. `ratio_median` = median per-pair S1/laser between the laser noise floor "
         "and λ/4. Dry subset and every column in `p6_windows.csv`.", "",
         md(w, "{:.2f}"), "",
         "## How to read this", "",
         "- Section 1 answers *where*; section 2 *whether it matters on the radar grids*; section 3 *whether P6's "
         "conclusions depend on it*. A 1.2 m walkway covers a few per cent of a 40 m cell; it matters only if it scatters "
         "strongly (dihedral with the mat or water, metal) or stays still while the mat moves.",
         "- 10 m backscatter is the finest radar grid available; interferometric products are 40 m (D-021: 20 m reprocessing).",
         "- The 10 m RTC is terrain-corrected GRD; its effective resolution is coarser than its pixel (multilooked), so a "
         "narrow structure is diluted even at 10 m.", "",
         "## Files", "", "| File | Content |", "|---|---|",
         "| `orthophoto_transect.tif` | GUGiK orthophoto around the transect, EPSG:32633, 0.25 m |",
         "| `boardwalk_mask.tif` | the traced boardwalk (1 = walkway or structure), same grid |",
         "| `boardwalk_centreline.csv`, `boardwalk_components.csv` | the trace |",
         "| `boardwalk_cells.csv` | boardwalk share of every 40/20/10 m cell it touches |",
         "| `boardwalk_cell_tests.csv`, `boardwalk_tests_summary.csv` | section 2 |",
         "| `p6_windows.csv` | section 3, all statistics, all and dry pairs |",
         "| `fig_boardwalk_outline.png`, `fig_boardwalk_tests.png`, `fig_p6_windows.png` | figures |"]
    (out / "README.md").write_text("\n".join(L) + "\n")


def fig_windows(out, win):
    w = win[win.subset == "all"]
    stats = [("r_s1_vs_laser", "r, S1 vs laser"), ("slope_tls", "TLS slope S1/laser"), ("ratio_median", "median per-pair ratio"),
             ("r_coh_vs_abs_dwtd", "r, coherence vs |ΔWTD|")]
    fig, axes = plt.subplots(1, len(stats), figsize=(16, 5), sharey=True)
    names = list(dict.fromkeys(w.window))
    for a, (k, lab) in zip(axes, stats):
        for t, col, off in (("ascending", "tab:orange", -0.15), ("descending", "tab:blue", 0.15)):
            g = w[w.track == t].set_index("window").reindex(names)
            a.scatter(g[k], np.arange(len(names)) + off, c=col, s=18, label=t)
        a.axvline(1 if "slope" in k or "ratio" in k else 0, color="0.6", lw=0.8)
        a.set_title(lab, fontsize=9)
    axes[0].set_yticks(range(len(names)), names, fontsize=8)
    axes[0].invert_yaxis()
    axes[0].legend(fontsize=7)
    fig.suptitle("P6: X-059 statistics by extraction window (consecutive pairs 2022–2024)", fontsize=10)
    fig.tight_layout()
    fig.savefig(out / "fig_p6_windows.png", dpi=110)
    plt.close(fig)


# ------------------------------------------------------------------------------ main

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default=str(DLV / "field_boardwalk_x060"))
    ap.add_argument("--rtc", default=str(V.REPO.parent / "local" / "drive_mirror"))
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    rtc_root = Path(a.rtc)
    px = pd.read_csv(DLV / "field_first" / "plot_pixels.csv")

    img, ext = orthophoto(out, px)
    masks, cl, comps = outline(img, ext, px)
    cl.to_csv(out / "boardwalk_centreline.csv", index=False)
    comps.to_csv(out / "boardwalk_components.csv", index=False)
    import rasterio
    from rasterio.transform import from_bounds
    E0, E1, N0, N1 = ext
    with rasterio.open(out / "boardwalk_mask.tif", "w", driver="GTiff", width=masks.shape[2], height=masks.shape[1], count=1,
                       dtype="uint8", crs="EPSG:32633", transform=from_bounds(E0, N0, E1, N1, masks.shape[2], masks.shape[1]),
                       compress="deflate") as f:
        f.write(masks[0].astype("uint8")[None])

    ctx = start("field_boardwalk_x060", mount=False, git=False)
    tpl = ctx.template
    x40, y40 = tpl.x.values, tpl.y.values
    mat40 = ctx.zones["A"].values.astype(bool)
    frac = {40: cell_fraction(masks[0], ext, x40, y40, 40.0)}
    off = {20: (-0.25, 0.25), 10: (-0.375, -0.125, 0.125, 0.375)}
    grids = {40: (x40, y40)}
    for res, o in off.items():
        gx = (x40[:, None] + np.array(o)[None, :] * 40).ravel()
        gy = (y40[:, None] - np.array(o)[None, :] * 40).ravel()
        grids[res] = (gx, gy)
        frac[res] = cell_fraction(masks[0], ext, gx, gy, float(res))
    rows = []
    for res, f in frac.items():
        gx, gy = grids[res]
        for i, j in zip(*np.nonzero(f)):
            rows.append({"res_m": res, "row": int(i), "col": int(j), "E": float(gx[j]), "N": float(gy[i]),
                         "boardwalk_share": float(f[i, j])})
    cells_bw = pd.DataFrame(rows)
    cells_bw.to_csv(out / "boardwalk_cells.csv", index=False)

    def plot_cell(res, e, n):
        gx, gy = grids[res]
        return int(np.argmin(np.abs(gy - n))), int(np.argmin(np.abs(gx - e)))
    ft = []
    for p in px.itertuples():
        ft.append({"plot": p.plot, **{f"{res} m": float(frac[res][plot_cell(res, p.E, p.N)]) for res in (40, 20, 10)}})
    frac_tab = pd.DataFrame(ft)
    labels40 = {}
    for p in px.itertuples():
        rc = (int(p.row), int(p.col))
        labels40[rc] = (labels40.get(rc, "") + " " + p.plot).strip()

    # 2. does it show?
    acq, p6, ov, qc, interp, wtd_s1, wet, zc = V.load_inputs()
    runs = V.filled_runs(interp)
    from insar_wetlands import field
    df = V.pair_table(acq, p6, ov, wtd_s1, wet, field.load_wtd_hourly(), runs)
    lay40 = radar_layers(ctx, df, rtc_root)
    cells, summ = [], []
    for (name, track), val in lay40.items():
        c, s = test_cells(name, track, 40, val, mat40, frac[40], labels40)
        cells += c
        summ += s
    gx10, gy10 = grids[10]
    lay10 = rtc10_layers(rtc_root, gx10, gy10)
    mat10 = np.kron(mat40, np.ones((4, 4), bool))
    labels10 = {plot_cell(10, p.E, p.N): p.plot for p in px.itertuples()}
    for (name, track), val in lay10.items():
        c, s = test_cells(name, track, 10, val, mat10, frac[10], labels10)
        cells += c
        summ += s
    cells = pd.DataFrame(cells)
    summ = pd.DataFrame(summ)
    cells.to_csv(out / "boardwalk_cell_tests.csv", index=False)
    summ.to_csv(out / "boardwalk_tests_summary.csv", index=False)

    # 3. P6 windows
    r0, c0 = int(px.set_index("plot").loc["P6", "row"]), int(px.set_index("plot").loc["P6", "col"])
    win = window_stats_table(ctx, df, float(qc["noise_floor_mm"]), r0, c0, frac[40] > 0, mat40)
    win.to_csv(out / "p6_windows.csv", index=False)

    fig_outline(out, img, ext, masks, x40, y40, frac[40], px)
    fig_tests(out, cells)
    fig_windows(out, win)
    readme(out, cl, comps, cells, frac_tab, summ, win, (r0, c0))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()

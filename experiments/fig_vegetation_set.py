"""Six figures on vegetation and the radar at Rzecin, each with one message, from registered outputs only (X-132, X-134,
C-074, C-065, the dated laser scans and the 10 cm orthophoto); nothing is refitted.

V1  What a radar cell holds: the nine cells around P6 as the 10 cm photo, the UAV true colour and NDVI, and the laser
    canopy height in September 2024 and April 2011 (leaf-off), and its change.
V2  Where the trees are and what the radar does there: canopy and shrub cover beside six-day coherence and the
    agreement with the P6 laser, every 40 m cell within 1 km.
V3  Coherence against cover: binned medians and interquartile ranges per track and pair set — tree cover on and off
    the mat, shrub cover on the tree-free mat.
V4  How large each effect is: the matched and dose contrasts of X-132 and X-134 on coherence, VH backscatter and the
    agreement with the laser, with 95 % spatial-block intervals.
V5  Is it the leaves? The same contrasts in summer (May–September) and winter (December–February).
V6  The vegetation through the seasons (UAV NDVI of the strip's cells, 2022–2024) beside the coherence of P6's cell
    over the same years.

Writes only to the hub (``08_deliverables/vegetation_and_radar/``). No field value in this file.

    make run E=fig_vegetation_set
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import fig_vegetation_radar as V  # noqa: E402  (its readers: effects(), ortho_rgb(), grid_image())
import landcover_coherence_x132 as X132  # noqa: E402  (its matching, to draw the pairs it scored)
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import shrub_coherence_x134 as X134  # noqa: E402
import xarray as xr  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib_scalebar.scalebar import ScaleBar  # noqa: E402

from insar_wetlands import contrasts as C  # noqa: E402
from insar_wetlands import figstyle as fs  # noqa: E402
from insar_wetlands.paths import cube_root, hub_root  # noqa: E402

HUB, CUBE = hub_root(), cube_root()
FINAL = HUB / "08_deliverables" / "article_plan" / "figures" / "final"
OUT = V.OUT
SETS = V.SETS
GROUPS = ["trees on the mat\nwooded − open", "trees off the mat\nwooded − open", "dense shrub − open\non the mat",
          "shrub 0 → 15 %\non the mat"]
GROUP_COL = {GROUPS[0]: "#2f6b3a", GROUPS[1]: "#8c6d31", GROUPS[2]: "#7a3fb0", GROUPS[3]: "#c27fd6"}
P6_COL = "#d55e00"


def nine_box(nine, pad=25):
    return (nine.x.min() - pad, nine.y.min() - pad, nine.x.max() + pad, nine.y.max() + pad)


def draw_grid(ax, nine, color="white"):
    for v in sorted(nine.x.unique()):
        for e in (v - 20, v + 20):
            ax.axvline(e, color=color, lw=0.6)
    for v in sorted(nine.y.unique()):
        for e in (v - 20, v + 20):
            ax.axhline(e, color=color, lw=0.6)
    for r in nine.itertuples():
        ax.text(r.x - 18, r.y + 18, r.cell, fontsize=7, color=color, va="top", ha="left", fontweight="bold")


def crop(da, box):
    sub = da.sel(x=slice(box[0], box[2]), y=slice(box[3], box[1]))
    return sub.values, (float(sub.x.min()), float(sub.x.max()), float(sub.y.min()), float(sub.y.max()))


def finish(ax, box):
    ax.set_xlim(box[0], box[2])
    ax.set_ylim(box[1], box[3])
    ax.set_xticks([])
    ax.set_yticks([])


# ------------------------------------------------------------------------------------------- V1

def v1(nine, uav):
    box = nine_box(nine)
    rgb, ext = V.ortho_rgb(box, res=0.1)
    d = V.UAV_DATE
    b, g, r = (uav[f"remx_{k}"].sel(date_remx=d) for k in ("Blue_475", "Green_560", "Red_668"))
    tc = np.dstack([crop(x, box)[0] for x in (r, g, b)])
    _, ext_u = crop(r, box)
    tc = np.clip(tc / np.nanpercentile(tc, 99), 0, 1)
    tc = np.where(np.isfinite(tc), tc, 0.87)
    nd = (uav.remx_NIR_842 - uav.remx_Red_668) / (uav.remx_NIR_842 + uav.remx_Red_668)
    lid = xr.open_dataset(CUBE / "silver" / "lidar_dated_1m.nc")
    with fs.style():
        fig, axes = plt.subplots(2, 3, figsize=fs.size("full", 140),
                                 gridspec_kw={"hspace": 0.28, "wspace": 0.06, "left": 0.01, "right": 0.99, "top": 0.93,
                                              "bottom": 0.07})
        ax = axes[0, 0]
        ax.imshow(rgb, extent=ext)
        draw_grid(ax, nine)
        ax.add_artist(ScaleBar(1, "m", location="lower right", length_fraction=0.25, box_alpha=0.85,
                               font_properties={"size": 7}, fixed_value=40))
        ax.set_title("(a) aerial photo, 10 cm, 19 Sep 2024", fontsize=7, loc="left")
        ax = axes[0, 1]
        ax.imshow(tc, extent=ext_u)
        draw_grid(ax, nine)
        ax.set_facecolor("#dddddd")
        ax.set_title(f"(b) UAV true colour, {pd.Timestamp(d):%-d %b %Y} (grey: no flight)", fontsize=7, loc="left")
        ax = axes[0, 2]
        a, e = crop(nd.sel(date_remx=d), box)
        im = ax.imshow(a, extent=e, cmap="YlGn", vmin=0.2, vmax=0.85)
        ax.set_facecolor("#dddddd")
        draw_grid(ax, nine, fs.INK)
        cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02, orientation="horizontal")
        cb.set_label("NDVI", fontsize=7)
        cb.ax.tick_params(labelsize=6.5)
        ax.set_title(f"(c) UAV NDVI, {pd.Timestamp(d):%-d %b %Y}", fontsize=7, loc="left")
        for k, (var, title) in enumerate((("chm_2024_m", "(d) canopy height, laser 19 Sep 2024 (leaf-on)"),
                                         ("chm_2011_m", "(e) canopy height, laser 24 Apr 2011 (leaf-off)"))):
            ax = axes[1, k]
            a, e = crop(lid[var], box)
            im = ax.imshow(a, extent=e, cmap="YlGn", vmin=0, vmax=2.5)
            draw_grid(ax, nine, fs.INK)
            cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02, orientation="horizontal")
            cb.set_label("m", fontsize=7)
            cb.ax.tick_params(labelsize=6.5)
            ax.set_title(title, fontsize=7, loc="left")
        ax = axes[1, 2]
        a, e = crop(lid.chm_2024_m - lid.chm_2011_m, box)
        im = ax.imshow(a, extent=e, cmap="BrBG", vmin=-1.5, vmax=1.5)
        draw_grid(ax, nine, fs.INK)
        cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02, orientation="horizontal")
        cb.set_label("m (2024 − 2011)", fontsize=7)
        cb.ax.tick_params(labelsize=6.5)
        ax.set_title("(f) canopy change 2011 → 2024", fontsize=7, loc="left")
        for ax in axes.ravel():
            finish(ax, box)
            ax.plot(nine[nine.cell == "P6"].x, nine[nine.cell == "P6"].y, "*", ms=7, color=P6_COL, mec=fs.INK,
                    mew=0.3)
        fs.export(fig, FINAL / "S16_nine_cells_contents",
                  provenance={"source": "fig_vegetation_set v1", "role": "Supplement S16"})
        fig.suptitle("V1 — what the nine 40 m radar cells around P6 hold", fontsize=9)
        fs.export(fig, OUT / "V1_what_the_cells_hold", provenance={"source": "fig_vegetation_set v1"})
        plt.close(fig)


# ------------------------------------------------------------------------------------------- V2

def cell_cover(cells):
    """Shrub share (0.5–2 m) and tree share (≥ 2 m) of each cell's 1 600 1 m laser pixels, 2024."""
    chm = xr.open_dataset(CUBE / "silver" / "lidar_dated_1m.nc").chm_2024_m
    out = []
    for r in cells.itertuples():
        a = chm.sel(x=slice(r.x - 20, r.x + 20), y=slice(r.y + 20, r.y - 20)).values
        a = a[np.isfinite(a)]
        out.append((((a >= 0.5) & (a < 2)).mean() if a.size else np.nan, (a >= 2).mean() if a.size else np.nan))
    return pd.DataFrame(out, columns=["shrub", "tree"], index=cells.index)


def v2(cells, outc, p6):
    c = cells.join(cell_cover(cells))
    o = outc.pivot_table(index="cell", columns=["outcome", "track", "set"], values="value")
    panels = [("tree", None, "(a) tree cover (≥ 2 m), share of cell", "Greens", (0, 1)),
              ("shrub", None, "(b) shrub cover (0.5–2 m), share of cell", "Purples", (0, 0.6)),
              (("O1_coh", "ascending", "2020_2021"), None, "(c) six-day coherence, ascending 2020–21", "viridis", (0.3, 0.9)),
              (("O1_coh", "descending", "2020_2021"), None, "(d) six-day coherence, descending 2020–21", "viridis", (0.3, 0.9)),
              (("O3_r_laser", "ascending", "2021+2025"), None, "(e) r with the P6 laser, ascending (6-day)", "RdBu_r", (-1, 1)),
              (("O3_r_laser", "descending", "2021+2025"), None, "(f) r with the P6 laser, descending (6-day)", "RdBu_r", (-1, 1))]
    with fs.style():
        fig, axes = plt.subplots(2, 3, figsize=fs.size("full", 145),
                                 gridspec_kw={"hspace": 0.22, "wspace": 0.05, "left": 0.01, "right": 0.99, "top": 0.92,
                                              "bottom": 0.08})
        mat, ext = V.grid_image(c, c.on_mat.to_numpy(float))
        wood, _ = V.grid_image(c, (c.cls == "wooded").to_numpy(float))
        gx = np.linspace(ext[0] + 20, ext[1] - 20, mat.shape[1])
        gy = np.linspace(ext[3] - 20, ext[2] + 20, mat.shape[0])
        for ax, (key, _, title, cmap, lim) in zip(axes.ravel(), panels):
            vals = c[key].to_numpy(float) if isinstance(key, str) else c.cell.map(o[key]).to_numpy(float)
            img, ext = V.grid_image(c, vals)
            im = ax.imshow(img, extent=ext, cmap=cmap, vmin=lim[0], vmax=lim[1], interpolation="nearest")
            ax.contour(gx, gy, np.where(np.isnan(mat), 0, mat), levels=[0.5], colors="black", linewidths=0.9)
            if not isinstance(key, str):
                ax.contour(gx, gy, np.where(np.isnan(wood), 0, wood), levels=[0.5], colors="#2f6b3a",
                           linewidths=0.5)
            ax.plot(p6.x, p6.y, "*", ms=8, color=P6_COL, mec="black", mew=0.4)
            ax.set_aspect("equal")
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(title, fontsize=7, loc="left")
            fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02, orientation="horizontal").ax.tick_params(labelsize=6)
        fig.legend(handles=[Line2D([], [], color="black", lw=0.9), Line2D([], [], color="#2f6b3a", lw=0.5),
                            Line2D([], [], ls="none", marker="*", color=P6_COL, mec="black")],
                   labels=["floating mat", "wooded cells (canopy ≥ 2 m)", "P6"], loc="lower center", ncol=3,
                   fontsize=6.5, frameon=False, bbox_to_anchor=(0.5, 0.0))
        fig.suptitle("V2 — where the trees are, and what the radar does there (40 m cells within 1 km of P6)",
                     fontsize=9, x=0.01, ha="left")
        fs.export(fig, OUT / "V2_maps", provenance={"source": "fig_vegetation_set v2"})
        plt.close(fig)


# ------------------------------------------------------------------------------------------- V3

BINS = [-0.001, 0.01, 0.05, 0.2, 0.5, 0.8, 1.0]
BIN_LAB = ["< 1", "1–5", "5–20", "20–50", "50–80", "> 80"]


def binned(x, y):
    b = pd.cut(x, BINS, labels=False)
    g = pd.DataFrame({"b": b, "y": y}).dropna().groupby("b").y
    q = g.quantile([0.25, 0.5, 0.75]).unstack()
    q["n"] = g.size()
    return q[q.n >= 5]


def v3(cells, outc, shrub_cells, outc_s):
    with fs.style():
        fig, axes = plt.subplots(2, 3, figsize=fs.size("full", 120), sharex=True,
                                 gridspec_kw={"hspace": 0.3, "wspace": 0.22, "left": 0.08, "right": 0.99, "top": 0.88,
                                              "bottom": 0.1})
        series = [("trees, on the mat", cells[cells.on_mat == 1], "canopy2024_cover", outc, GROUP_COL[GROUPS[0]], -0.12),
                  ("trees, off the mat", cells[cells.on_mat == 0], "canopy2024_cover", outc, GROUP_COL[GROUPS[1]], 0.0),
                  ("shrubs, tree-free mat", shrub_cells, "shrub_2024", outc_s, GROUP_COL[GROUPS[2]], 0.12)]
        for i, track in enumerate(("ascending", "descending")):
            for j, s in enumerate(SETS):
                ax = axes[i, j]
                for lab, cc, xcol, oo, col, dx in series:
                    y = oo[(oo.outcome == "O1_coh") & (oo.track == track) & (oo.set == s)].set_index("cell").value
                    q = binned(cc.set_index("cell")[xcol], y)
                    xs = q.index.to_numpy() + dx
                    ax.vlines(xs, q[0.25], q[0.75], color=col, lw=1.2)
                    ax.plot(xs, q[0.5], "-o", color=col, ms=3.5, lw=0.9, label=lab if (i, j) == (0, 0) else None)
                ax.set_xticks(range(len(BIN_LAB)), BIN_LAB)
                ax.set_title(f"{track}, {SETS[s][0]}", fontsize=7, loc="left")
                if j == 0:
                    ax.set_ylabel("median coherence of\nconsecutive pairs (per cell)")
                if i == 1:
                    ax.set_xlabel("cover of the 40 m cell (%)")
        fig.legend(loc="upper left", ncol=3, fontsize=7, frameon=False, bbox_to_anchor=(0.08, 0.97))
        fig.suptitle("V3 — coherence against tree and shrub cover: median and interquartile range of the cells in each "
                     "cover class (classes with ≥ 5 cells)", fontsize=8, x=0.01, ha="left", y=0.995)
        fs.export(fig, OUT / "V3_coherence_against_cover", provenance={"source": "fig_vegetation_set v3"})
        plt.close(fig)


# ------------------------------------------------------------------------------------------- V4, V5

def forest_cols(df_by_col, titles, path, suptitle, sets_by_col, ref=None, sharex=False):
    with fs.style():
        fig, axes = plt.subplots(1, len(df_by_col), figsize=fs.size("full", 95), sharey=True, sharex=sharex,
                                 gridspec_kw={"wspace": 0.08, "left": 0.16, "right": 0.99, "top": 0.82, "bottom": 0.14})
        for ax, df, title, sets in zip(axes, df_by_col, titles, sets_by_col):
            for gi, g in enumerate(GROUPS):
                if gi % 2 == 0:
                    ax.axhspan(gi - 0.5, gi + 0.5, color=fs.GRID, alpha=0.35, lw=0)
                rows = [(t, s) for t in ("ascending", "descending") for s in sets]
                for k, (t, s) in enumerate(rows):
                    r = df[(df.group == g) & (df.track == t) & (df.set == s)]
                    if not len(r):
                        continue
                    r = r.iloc[0]
                    y = gi - 0.36 + 0.72 * k / max(len(rows) - 1, 1)
                    ax.hlines(y, r.lo, r.hi, color=fs.TRACK[t], lw=1.1)
                    ax.plot(r.est, y, SETS[s][1], ms=4, color=fs.TRACK[t], mec="white", mew=0.4)
            ax.axvline(0, color=fs.INK, lw=0.7)
            if ref and title.startswith("coherence"):
                ax.axvspan(-ref, ref, color="#f3e3c3", alpha=0.6, lw=0, zorder=0)
            ax.set_title(title, fontsize=7.5, loc="left")
            ax.grid(axis="y", visible=False)
        axes[0].set_yticks(range(len(GROUPS)), GROUPS, fontsize=7)
        axes[0].set_ylim(len(GROUPS) - 0.5, -0.5)
        handles = [Line2D([], [], color=fs.TRACK[t], lw=1.1, label=f"{t} track") for t in ("ascending", "descending")]
        handles += [Line2D([], [], ls="none", marker=SETS[s][1], color=fs.MUTED, label=SETS[s][0]) for s in SETS]
        fig.legend(handles=handles, loc="upper left", ncol=5, fontsize=6.5, frameon=False, bbox_to_anchor=(0.16, 0.93))
        fig.suptitle(suptitle, fontsize=8.5, x=0.01, ha="left", y=0.99)
        fs.export(fig, path, provenance={"source": "fig_vegetation_set"})
        plt.close(fig)


def v4():
    o3 = V.effects("O3_r_laser")
    cols = [V.effects("O1_coh"), V.effects("O4_vh_median_db"), o3[o3.group != GROUPS[1]]]
    titles = ["coherence of consecutive pairs (Δ, 0–1)", "VH backscatter (Δ dB)",
              "r with the P6 laser (Δ, six-day), mat only"]
    forest_cols(cols, titles, OUT / "V4_effect_sizes",
                "V4 — how large each vegetation effect is: contrast with 95 % interval (200 m spatial-block bootstrap); "
                "shaded: X-134's ±0.03 margin", [list(SETS), list(SETS), ["2020_2021", "2025"]], ref=0.03)


def v5():
    cols = [V.effects("O1_coh_summer"), V.effects("O1_coh_winter")]
    forest_cols(cols, ["coherence, summer (May–September, leaf-on)", "coherence, winter (December–February, leaf-off)"],
                OUT / "V5_summer_winter",
                "V5 — is it the leaves? The same contrasts on the coherence of consecutive pairs, by season",
                [list(SETS), list(SETS)], sharex=True)


# ------------------------------------------------------------------------------------------- V6

def v6(nine):
    u = pd.read_csv(HUB / "08_deliverables" / "field_uav_c065" / "cells_uav.csv")
    meta = pd.read_csv(HUB / "08_deliverables" / "field_uav_c065" / "cells_meta.csv")
    u = u[(u.sensor == "RedEdge_MX_Dual") & u.qc_ok].merge(meta[["x", "y", "zone", "depth_in_mat_m"]], on=["x", "y"])
    u["date"] = pd.to_datetime(u.date)
    p = pd.read_csv(CUBE / "gold" / "plots_pair.csv.gz", low_memory=False,
                    usecols=["plot", "grid", "track", "t1", "t2", "coh_1x1", "platforms", "period"])
    p = p[(p["plot"] == "P6") & (p.grid == "40m") & ~p.platforms.astype(str).str.contains("S1D")]
    p["t1"] = pd.to_datetime(p.t1, format="ISO8601", utc=True).dt.tz_convert(None)
    p["t2"] = pd.to_datetime(p.t2, format="ISO8601", utc=True).dt.tz_convert(None)
    p = p[((p.t2 - p.t1).dt.days.between(11, 13)) & p.t1.dt.year.between(2022, 2024)]
    p["mid"] = p.t1 + (p.t2 - p.t1) / 2
    ninexy = set(zip(nine.x, nine.y))
    with fs.style():
        fig, axes = plt.subplots(2, 1, figsize=fs.size("full", 110), sharex=True,
                                 gridspec_kw={"hspace": 0.28, "left": 0.09, "right": 0.98, "top": 0.88, "bottom": 0.12})
        ax = axes[0]
        for (x, y, yr), g in u.assign(yr=u.date.dt.year).groupby(["x", "y", "yr"]):
            g = g.sort_values("date")
            hit = (x, y) in ninexy
            ax.plot(g.date, g.ndvi_median, "-o" if hit else "-o", ms=3 if hit else 1.2, lw=1.2 if hit else 0.5,
                    color=P6_COL if hit else fs.MUTED, alpha=1 if hit else 0.6, zorder=3 if hit else 1)
        ax.set_ylabel("UAV NDVI, median of the cell")
        ax.legend(handles=[Line2D([], [], color=P6_COL, marker="o", ms=3, lw=1.2),
                           Line2D([], [], color=fs.MUTED, lw=0.5)],
                  labels=["the nine cells around P6 that the strip covers", "other strip cells on the mat"],
                  fontsize=6.5, loc="upper left", bbox_to_anchor=(0.0, -0.02), ncol=2, frameon=False)
        ax.set_title("(a) how green the cells are, flight by flight (usable flights only; lines join flights of one "
                     "year)", fontsize=7.5, loc="left")
        ax = axes[1]
        for t in ("ascending", "descending"):
            g = p[p.track == t].sort_values("mid")
            ax.plot(g.mid, g.coh_1x1, "-", lw=0.7, color=fs.TRACK[t], marker="o", ms=1.8, label=f"{t} track")
        ax.set_ylabel("coherence of P6's cell\n(consecutive 12-day pairs)")
        ax.legend(fontsize=6.5, loc="upper left", bbox_to_anchor=(0.0, -0.08), ncol=2, frameon=False)
        ax.set_title("(b) the radar over the same years: lowest in summer, when the cells are greenest (all pairs, frost and snow included)", fontsize=7.5, loc="left")
        for a in axes:
            for y in (2022, 2023, 2024):
                a.axvspan(pd.Timestamp(f"{y}-05-01"), pd.Timestamp(f"{y}-09-30"), color="#e6f2dc", lw=0, zorder=0)
        fig.suptitle("V6 — the vegetation through the seasons and the coherence of P6's cell, 2022–2024 "
                     "(green bands: May–September)", fontsize=8.5, x=0.01, ha="left")
        fs.export(fig, OUT / "V6_seasons", provenance={"source": "fig_vegetation_set v6"})
        plt.close(fig)


CAPTIONS = """# Vegetation and the radar — six figures (generated)

Generated by `05_code/SAr-adapted/experiments/fig_vegetation_set.py`; do not edit. Every value is read from
registered outputs (X-132 trees, X-134 shrubs, C-074 the nine cells, C-065 the UAV stack, the dated laser scans and the
10 cm orthophoto); nothing is refitted and no new test is made. The single summary figure is `vegetation_and_radar.png`
(README.md).

## V1 — what the nine radar cells around P6 hold
![V1](V1_what_the_cells_hold.png)
The nine 40 m cells as the 10 cm photo of 19 September 2024, the UAV true colour and NDVI of the next day (the strip
covers the middle column), and the laser canopy height in September 2024 (leaf-on) and April 2011 (leaf-off), with
the change. The cells are open fen with scattered low shrubs; the boardwalk and the instrument platforms cross P6's
cell and the one south of it. Part of the 2011 → 2024 increase is season (leaf-off against leaf-on), not growth alone.

## V2 — where the trees are, and what the radar does there
![V2](V2_maps.png)
Every 40 m cell within 1 km of P6: tree and shrub cover from the 2024 laser scan beside the six-day coherence
(2020–21) and the agreement of each cell's six-day phase with the P6 laser, both tracks; black: the floating mat;
green outline: wooded cells. On the ascending track coherence is highest on the open mat around P6 and low in the
woods off the mat; on the descending track it is lower over much of the mat too. Agreement with the laser covers the
mat and fades at its edge, through wooded and open cells alike.

## V3 — coherence against cover (raw, not adjusted)
![V3](V3_coherence_against_cover.png)
Median and interquartile range of the cells' coherence in each cover class (classes with at least five cells), per
track and pair set: tree cover on and off the mat, shrub cover on the tree-free mat. These are raw, not adjusted for
depth into the mat or position, which move together with cover: the matched and adjusted contrasts are V4.

## V4 — how large each effect is
![V4](V4_effect_sizes.png)
The pre-registered contrasts with their 95 % intervals (200 m spatial-block bootstrap): trees (wooded − open, matched
on depth) on and off the mat, dense shrub − open (matched), and 0 → 15 % shrub cover (dose–response); on coherence,
VH backscatter and the agreement of the six-day phase with the P6 laser. Shaded: X-134's ±0.03 margin for "no material
difference". Trees lower coherence and raise VH, on the mat as off it; shrubs do the same, less; the agreement with the
laser is not weakened on the mat by the registered rules (one exception, descending 2021, the laser's part of the 2020–21 set; X-132).

## V5 — is it the leaves?
![V5](V5_summer_winter.png)
The same contrasts on coherence in summer (May–September) and winter (December–February), on one scale. The tree
effect stays in winter: the loss does not need leaves, though stems, branches and evergreen crowns cannot be told apart
here; the shrub effect is weaker in winter in most sets, as leaves would make it. Descriptive reading of registered secondaries; not a new test.

## V6 — the vegetation through the seasons
![V6](V6_seasons.png)
(a) The UAV NDVI of the strip's 40 m cells on the mat, flight by flight (usable flights only); orange: the cells of
the nine that the strip covers. (b) The coherence of P6's cell on consecutive twelve-day pairs over the same years.
The cells green up from May and the coherence is lowest in summer; the two are shown together, not tested together.
"""


# ------------------------------------------------------------------------------------------- Supplement S15

def matched_pairs(cells, shrub_cells):
    """The pairs X-132 and X-134 scored, rebuilt with their own matching (checked against their registered counts)."""
    c = cells.set_index("cell")
    m = X132.matched(c)
    sc = shrub_cells.set_index("cell")
    ms = C.match_nearest(sc[sc.cls == "shrubby"], sc[sc.cls == "open"], on="depth_in_mat_m",
                         caliper=X134.CALIPER_M, tie="dist_p6_m", order="dist_p6_m")
    s1 = pd.read_csv(V.X132 / "stage1.csv")
    s1 = s1[(s1.variant == "primary") & (s1.outcome == "O1_coh")].iloc[0]
    mt = pd.read_csv(V.X134 / "matched.csv")
    mt = mt[(mt.variant == "primary") & (mt.outcome == "O1_coh")].iloc[0]
    if (len(m["on_mat"]), len(m["off_mat"]), len(ms)) != (int(s1.pairs_on_mat), int(s1.pairs_off_mat), int(mt.pairs)):
        raise SystemExit("S15: the rebuilt pairs do not match X-132 / X-134's registered counts")
    return m["on_mat"], m["off_mat"], ms, c, sc


def pair_distances(on, off, ms, c, sc) -> pd.DataFrame:
    """How far apart each contrast's paired cells are, and how their distances to P6 compare (pairs are matched on
    depth into the mat, not on position)."""
    rows = []
    for name, pairs, tab in (("trees on the mat", on, c), ("trees off the mat", off, c), ("dense shrub − open", ms, sc)):
        a, b = tab.loc[pairs.treated], tab.loc[pairs.control]
        sep = np.hypot(a.x.to_numpy() - b.x.to_numpy(), a.y.to_numpy() - b.y.to_numpy())
        dd = a.dist_p6_m.to_numpy() - b.dist_p6_m.to_numpy()
        rows.append({"contrast": name, "pairs": len(pairs), "separation_median_m": float(np.median(sep)),
                     "separation_max_m": float(sep.max()), "dist_p6_difference_median_m": float(np.median(dd)),
                     "dist_p6_abs_difference_median_m": float(np.median(np.abs(dd)))})
    return pd.DataFrame(rows)


def s15(cells, shrub_cells, nine, p6):
    on, off, ms, c, sc = matched_pairs(cells, shrub_cells)
    pair_distances(on, off, ms, c, sc).to_csv(OUT / "S15_pair_distances.csv", index=False)
    with fs.style():
        fig = plt.figure(figsize=fs.size("full", 215))
        gs = fig.add_gridspec(2, 3, height_ratios=[1.15, 1.0], hspace=0.32, wspace=0.08, left=0.16, right=0.99,
                              top=0.94, bottom=0.06)
        top_gs = gs[0, :].subgridspec(1, 2, width_ratios=[1.05, 1.0], wspace=0.10)
        mat, ext = V.grid_image(cells, cells.on_mat.to_numpy(float))
        gx = np.linspace(ext[0] + 20, ext[1] - 20, mat.shape[1])
        gy = np.linspace(ext[3] - 20, ext[2] + 20, mat.shape[0])
        # (a) trees
        ax = fig.add_subplot(top_gs[0])
        ax.set_anchor("NW")
        col = cells.cls.map({"open": 0, "mid": 1, "wooded": 2}).to_numpy(float)
        img, e = V.grid_image(cells, col)
        from matplotlib.colors import ListedColormap
        ax.imshow(img, extent=e, cmap=ListedColormap(["#efe6cf", "#e4e4e4", "#2f6b3a"]), vmin=-0.5, vmax=2.5,
                  interpolation="nearest")
        for pairs, colr in ((on, "#1b9e77"), (off, "#8c6d31")):
            for r in pairs.itertuples():
                a, b = c.loc[r.treated], c.loc[r.control]
                ax.plot([a.x, b.x], [a.y, b.y], color=colr, lw=0.35, alpha=0.3)
            t, k = c.loc[pairs.treated], c.loc[pairs.control]
            ax.plot(t.x, t.y, "s", ms=2.2, mfc="none", mec=colr, mew=0.6)
            ax.plot(k.x, k.y, "o", ms=2.0, color=colr)
        ax.contour(gx, gy, np.where(np.isnan(mat), 0, mat), levels=[0.5], colors="black", linewidths=0.9)
        ax.plot(p6.x, p6.y, "*", ms=9, color=P6_COL, mec="black", mew=0.4, zorder=6)
        ax.add_artist(ScaleBar(1, "m", location="lower right", length_fraction=0.22, box_alpha=0.85,
                               font_properties={"size": 7}, fixed_value=500))
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(f"(a) trees: {len(on)} wooded–open pairs on the mat, {len(off)} off it — matched on depth into (or "
                     "distance to) the mat,\nnot on position, so a pair can be far apart (faint lines)", fontsize=7,
                     loc="left")
        ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, color="#2f6b3a"), plt.Rectangle((0, 0), 1, 1, color="#efe6cf"),
                           plt.Rectangle((0, 0), 1, 1, color="#e4e4e4"), Line2D([], [], color="#1b9e77", lw=0.8),
                           Line2D([], [], color="#8c6d31", lw=0.8), Line2D([], [], color="black", lw=0.9),
                           Line2D([], [], ls="none", marker="s", mfc="none", mec=fs.INK, ms=3),
                           Line2D([], [], ls="none", marker="o", color=fs.INK, ms=3)],
                  labels=["wooded (canopy ≥ 2 m)", "open (canopy < 0.5 m)", "in between (not used)",
                          "pair on the mat", "pair off the mat", "floating mat", "wooded cell used",
                          "open cell used"],
                  fontsize=6.5, loc="upper left", bbox_to_anchor=(-0.02, -0.01), ncol=2, frameon=False,
                  columnspacing=1.1, handletextpad=0.45, labelspacing=0.22)
        # (b) shrubs
        ax = fig.add_subplot(top_gs[1])
        ax.set_anchor("N")
        allc = cells.set_index("cell")
        base, e0 = V.grid_image(cells, np.where(cells.on_mat == 1, 0.0, np.nan))
        ax.imshow(base, extent=e0, cmap=ListedColormap(["#f2f2f2"]), interpolation="nearest")
        img, e = V.grid_image(shrub_cells, shrub_cells.shrub_2024.to_numpy(float))
        im = ax.imshow(img, extent=e, cmap="Purples", vmin=0, vmax=0.6, interpolation="nearest")
        for r in ms.itertuples():
            a, b = sc.loc[r.treated], sc.loc[r.control]
            ax.plot([a.x, b.x], [a.y, b.y], color="#d95f02", lw=0.45, alpha=0.38, zorder=2)
            ax.plot(b.x, b.y, "o", ms=1.8, color="#d95f02", alpha=0.75, zorder=3)
        ax.contour(gx, gy, np.where(np.isnan(mat), 0, mat), levels=[0.5], colors="black", linewidths=0.9, zorder=4)
        ax.add_patch(plt.Rectangle((nine.x.min() - 20, nine.y.min() - 20), 120, 120, fill=False, ec=P6_COL, lw=1.1,
                                   zorder=5))
        ax.plot(p6.x, p6.y, "*", ms=8.5, color=P6_COL, mec="black", mew=0.45, zorder=6)
        ax.add_artist(ScaleBar(1, "m", location="lower right", length_fraction=0.22, box_alpha=0.85,
                               font_properties={"size": 7}, fixed_value=250))
        xs = [allc.x.min(), allc.x.max()]
        ax.set_xlim(min(shrub_cells.x) - 60, max(shrub_cells.x) + 60)
        ax.set_ylim(min(shrub_cells.y) - 60, max(shrub_cells.y) + 60)
        _ = xs
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.02, orientation="horizontal")
        cb.set_label("shrub cover 0.5–2 m (share of cell)", fontsize=7)
        cb.ax.tick_params(labelsize=6.5)
        ax.set_title(f"(b) shrubs: {len(sc)} tree-free mat cells,\n{len(ms)} dense–open pairs (orange)",
                     fontsize=7, loc="left")
        # (c–e) effect sizes
        o3 = V.effects("O3_r_laser")
        o3 = o3[o3.group != GROUPS[1]]                       # off the mat the ground does not move with the mat
        near = V.effects("O3_r_laser", "radius 500 m")
        near = near[near.group != GROUPS[1]]
        cols = [V.effects("O1_coh"), V.effects("O4_vh_median_db"), o3]
        titles = ["(c) coherence of consecutive pairs (Δ)", "(d) VH backscatter (Δ dB)",
                  "(e) r with the P6 laser (Δ), 2021 and 2025"]
        sets_by = [list(SETS), list(SETS), ["2020_2021", "2025"]]
        axes = [fig.add_subplot(gs[1, k]) for k in range(3)]
        for k, (ax, df, title, sets) in enumerate(zip(axes, cols, titles, sets_by)):
            for gi, g in enumerate(GROUPS):
                if gi % 2 == 0:
                    ax.axhspan(gi - 0.5, gi + 0.5, color=fs.GRID, alpha=0.35, lw=0)
                rows = [(t, st) for t in ("ascending", "descending") for st in sets]
                for j, (t, st) in enumerate(rows):
                    r = df[(df.group == g) & (df.track == t) & (df.set == st)]
                    if not len(r):
                        continue
                    r = r.iloc[0]
                    y = gi - 0.36 + 0.72 * j / max(len(rows) - 1, 1)
                    ax.hlines(y, r.lo, r.hi, color=fs.TRACK[t], lw=1.1)
                    ax.plot(r.est, y, SETS[st][1], ms=3.8, color=fs.TRACK[t], mec="white", mew=0.4)
                    if k == 2:
                        q = near[(near.group == g) & (near.track == t) & (near.set == st)]
                        if len(q):
                            ax.plot(q.est.iloc[0], y + 0.07, SETS[st][1], ms=3.6, mfc="none", mec=fs.TRACK[t],
                                    mew=0.7)
            ax.axvline(0, color=fs.INK, lw=0.7)
            if k == 2:
                ax.text(0.5, 1 / len(GROUPS) * 2.5, "not shown: solid ground\ndoes not move with the mat,\nso its "
                        "agreement with the\nP6 laser has no meaning", transform=ax.transAxes, ha="center",
                        va="center", fontsize=6.5, color=fs.MUTED,
                        bbox={"boxstyle": "round,pad=0.3", "fc": "white", "ec": "none"})
            if k == 0:
                ax.axvspan(-0.03, 0.03, color="#f3e3c3", alpha=0.6, lw=0, zorder=0)
                ax.set_yticks(range(len(GROUPS)), GROUPS, fontsize=7)
            else:
                ax.set_yticks(range(len(GROUPS)), [""] * len(GROUPS))
            ax.set_ylim(len(GROUPS) - 0.5, -0.5)
            ax.set_title(title, fontsize=7, loc="left")
            ax.grid(axis="y", visible=False)
            ax.tick_params(axis="x", labelsize=7)
        handles = [Line2D([], [], color=fs.TRACK[t], lw=1.1, label=f"{t} track") for t in ("ascending", "descending")]
        handles += [Line2D([], [], ls="none", marker=SETS[st][1], color=fs.MUTED, label=SETS[st][0]) for st in SETS]
        handles += [plt.Rectangle((0, 0), 1, 1, color="#f3e3c3", label="±0.03 (X-134's margin)"),
                     Line2D([], [], ls="none", marker="o", mfc="none", mec=fs.MUTED,
                            label="(e) only cells within 500 m of P6")]
        fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=6.8, frameon=False, bbox_to_anchor=(0.57, 0.0))
        for stem in (OUT / "S15_vegetation_effects", FINAL / "S15_vegetation_effects"):
            fs.export(fig, stem, provenance={"source": "fig_vegetation_set s15", "pairs_checked": True})
        plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    cells = pd.read_csv(V.X132 / "cells.csv")
    outc = pd.read_csv(V.X132 / "outcomes.csv.gz")
    shrub_cells = pd.read_csv(V.X134 / "cells.csv")
    outc_s = pd.read_csv(V.X134 / "outcomes.csv.gz")
    nine = pd.read_csv(V.C074 / "nine_cells.csv")
    p6 = nine[nine.cell == "P6"].iloc[0]
    uav = xr.open_dataset(V.UAV)
    v1(nine, uav)
    v2(cells, outc, p6)
    v3(cells, outc, shrub_cells, outc_s)
    v4()
    v5()
    v6(nine)
    s15(cells, shrub_cells, nine, p6)
    (OUT / "FIGURES.md").write_text(CAPTIONS, encoding="utf-8")
    print("V1–V6 →", OUT)


if __name__ == "__main__":
    main()

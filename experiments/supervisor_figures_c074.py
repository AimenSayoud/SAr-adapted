"""C-074 — the figures of the supervisor pack: the nine cells around P6 (site, tree cover, backscatter, coherence and
phase against the laser), the phase against the laser at P6 and its neighbours, P6 through time; and the C-073
figures of the wider analysis copied in beside them.

Drawn from C-074's nine-cell tables, the cube and the committed C-073 figures; nothing is re-estimated. House style
(``insar_wetlands.figstyle``); colour is never the only cue. Writes only to the hub
(``08_deliverables/supervisor_pack_c074/``).

    make run E=supervisor_figures_c074
"""
from __future__ import annotations

import shutil
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import xarray as xr  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import figures_landcover_footprint_c073 as F73  # noqa: E402  (map axes, export, the C-073 figures)
from insar_wetlands import figstyle as fs  # noqa: E402
from insar_wetlands.cube.core import WAVELENGTH_MM  # noqa: E402
from insar_wetlands.paths import cube_root, hub_root  # noqa: E402

HUB, CUBE = hub_root(), cube_root()
PACK = HUB / "08_deliverables" / "supervisor_pack_c074"
TAB = PACK / "tables"
C073 = HUB / "08_deliverables" / "figures_landcover_footprint"
TRACKS = ("ascending", "descending")
GROUPS = ["2021", "2022–2024", "2025"]
GLABEL = {"2021": "2021, 6-day", "2022–2024": "2022–24, 12-day", "2025": "2025, 6-day"}
LAYOUT = [["NW", "N", "NE"], ["W", "P6", "E"], ["SW", "S", "SE"]]
P6C, NBC = fs.CLASS["laser"], fs.CLASS["mat"]          # P6 cell; mean of its 8 neighbours
CELL = 40.0
QW = WAVELENGTH_MM / 4


def export(fig, name, sources, notes, data):
    (PACK / "data").mkdir(parents=True, exist_ok=True)
    for k, d in data.items():
        d.to_csv(PACK / "data" / f"{name}_{k}.csv", index=False)
    with fs.style():
        fs.export(fig, PACK / name, provenance={"sources": sources, "data": [f"data/{name}_{k}.csv" for k in data],
                                                **notes})
    plt.close(fig)


def p6_station() -> tuple[float, float]:
    px = pd.read_csv(HUB / "08_deliverables" / "field_first" / "plot_pixels.csv").set_index("plot")
    return float(px.loc["P6", "E"]), float(px.loc["P6", "N"])


# ------------------------------------------------------------------------------ S1 — the nine cells, 19 Sep 2024

def ortho_rgb(xlim, ylim, res=0.2) -> np.ndarray:
    """The 10 cm true-colour orthophoto of 19 Sep 2024 over the box, in UTM 33N at ``res`` m (band means), as delivered."""
    import rasterio
    from pyproj import Transformer
    from rasterio.merge import merge
    from rasterio.transform import from_origin
    from rasterio.warp import Resampling, reproject
    raw = CUBE / "raw" / "gugik_dated"
    tab = pd.read_csv(raw / "sheets.csv")
    files = [raw / f for f in tab[(tab.kind == "rgb") & (tab.year == 2024)].file]
    t = Transformer.from_crs(32633, 2180, always_xy=True)
    ex, ey = t.transform([xlim[0], xlim[1], xlim[0], xlim[1]], [ylim[0], ylim[0], ylim[1], ylim[1]])
    b = (min(ex) - 5, min(ey) - 5, max(ex) + 5, max(ey) + 5)
    srcs = [rasterio.open(f) for f in files]
    img, tr = merge(srcs, bounds=b)
    crs = srcs[0].crs
    for s in srcs:
        s.close()
    W, H = int(round((xlim[1] - xlim[0]) / res)), int(round((ylim[1] - ylim[0]) / res))
    out = np.zeros((3, H, W), "uint8")
    for k in range(3):
        reproject(img[k], out[k], src_transform=tr, src_crs=crs, dst_transform=from_origin(xlim[0], ylim[1], res, res),
                  dst_crs="EPSG:32633", resampling=Resampling.average)
    return np.moveaxis(out, 0, -1)


def grid_and_labels(ax, cells, label_fn, fc="white", ec="white", bottom=False):
    for name, c in cells.iterrows():
        ax.add_patch(Rectangle((c.x - CELL / 2, c.y - CELL / 2), CELL, CELL, fc="none",
                               ec=P6C if name == "P6" else ec, lw=1.4 if name == "P6" else 0.8, zorder=4))
        y = c.y - CELL / 2 + 2 if bottom else c.y + CELL / 2 - 3
        ax.text(c.x, y, label_fn(name, c), ha="center", va="bottom" if bottom else "top", fontsize=6.3, zorder=5,
                bbox={"boxstyle": "round,pad=0.25", "fc": fc, "ec": "none", "alpha": 0.85})


def s1_site(cells: pd.DataFrame) -> None:
    pad = 25
    xlim = (cells.x.min() - CELL / 2 - pad, cells.x.max() + CELL / 2 + pad)
    ylim = (cells.y.min() - CELL / 2 - pad, cells.y.max() + CELL / 2 + pad)
    img = ortho_rgb(xlim, ylim)
    li = xr.open_dataset(CUBE / "silver" / "lidar_dated_1m.nc")
    chm = li.chm_2024_m.sel(x=slice(*xlim), y=slice(ylim[1], ylim[0]))
    sx, sy = p6_station()
    lab = lambda n, c: f"{n}\nshrub {c.shrub_share_2024:.0%}\ntree {c.tree_share_2024:.0%}"  # noqa: E731
    with fs.style():
        fig, axs = plt.subplots(1, 2, figsize=fs.size("full", 102), layout="constrained")
        ax = axs[0]
        ax.imshow(img, extent=[*xlim, *ylim], origin="upper", interpolation="none", zorder=0)
        grid_and_labels(ax, cells, lambda n, c: f"{n}\nboardwalk {c.boardwalk_share:.0%}")
        ax.plot(sx, sy, "*", ms=10, mfc=P6C, mec="white", mew=0.7, zorder=6)
        ax.set_title("True-colour orthophoto, 10 cm, 19 Sep 2024", loc="left")
        F73.map_axes(ax, xlim, ylim, scalebar=True, north=True)
        ticks = np.arange(np.ceil(xlim[0] / 50) * 50, xlim[1], 50)
        ax.set_xticks(ticks)
        ax.set_yticks(np.arange(np.ceil(ylim[0] / 50) * 50, ylim[1], 50))
        ax.xaxis.set_major_formatter(lambda v, _: f"{v:.0f}")
        ax.yaxis.set_major_formatter(lambda v, _: f"{v:.0f}")
        ax.tick_params(axis="x", labelrotation=30)
        ax.set_xlabel("Easting (m, UTM 33N)")
        ax.set_ylabel("Northing (m, UTM 33N)")
        fs.panel_label(ax, "a")
        ax = axs[1]
        im = ax.imshow(chm.values, extent=[float(chm.x[0]) - 0.5, float(chm.x[-1]) + 0.5, float(chm.y[-1]) - 0.5,
                                           float(chm.y[0]) + 0.5], origin="upper", interpolation="none",
                       cmap=plt.get_cmap("YlGn"), vmin=0, vmax=3, zorder=0)
        grid_and_labels(ax, cells, lab, ec=fs.MUTED, bottom=True)
        ax.plot(sx, sy, "*", ms=8, mfc=P6C, mec="white", mew=0.6, zorder=4.5)        # under the cell labels
        ax.set_title("Vegetation height, laser scan the same day (1 m)", loc="left")
        F73.map_axes(ax, xlim, ylim, scalebar=False, north=False, ylabel=False)
        ax.set_xticks(ticks)
        ax.xaxis.set_major_formatter(lambda v, _: f"{v:.0f}")
        ax.tick_params(axis="x", labelrotation=30)
        ax.set_xlabel("Easting (m, UTM 33N)")
        fig.colorbar(im, ax=ax, shrink=0.8, extend="max", label="Canopy height (m); shrub 0.5–2 m, tree ≥ 2 m")
        fs.panel_label(ax, "b")
    export(fig, "S1_nine_cells_site", ["cube:raw/gugik_dated (RGB 2024, 10 cm)", "cube:silver/lidar_dated_1m.nc",
                                       "hub:supervisor_pack_c074/tables/nine_cells.csv",
                                       "hub:field_first/plot_pixels.csv (P6 station)"],
           {"display": "orthophoto as delivered, resampled to 0.2 m; white squares: the 40 m radar cells (corrected grid)"},
           {"cells": cells.reset_index()})


# ------------------------------------------------------------------------------ S2 — what differs between the nine cells

def tiles(ax, values: dict, fmt: str, cmap: str, vmin=None, vmax=None, title=""):
    a = np.array([[values[c] for c in row] for row in LAYOUT], float)
    vmin = np.nanmin(a) if vmin is None else vmin
    vmax = np.nanmax(a) if vmax is None else vmax
    im = ax.imshow(a, cmap=cmap, vmin=vmin, vmax=vmax, interpolation="nearest")
    cm = plt.get_cmap(cmap)
    for i in range(3):
        for j in range(3):
            v = a[i, j]
            rgba = cm((v - vmin) / (vmax - vmin) if vmax > vmin else 0.5)
            lum = 0.2126 * rgba[0] + 0.7152 * rgba[1] + 0.0722 * rgba[2]
            ax.text(j, i, f"{LAYOUT[i][j]}\n{format(v, fmt)}", ha="center", va="center", fontsize=6.2,
                    color="white" if lum < 0.5 else fs.INK)
    ax.add_patch(Rectangle((0.5, 0.5), 1, 1, fc="none", ec=P6C, lw=1.6))
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title(title, fontsize=7, loc="left")
    return im


READ = ("Each 3 × 3 block is the nine 40 m radar cells around P6 as they lie on the ground (north up); "
        "P6's own cell is framed in orange.")
ROWLAB = {"2021": "2021 · 6-day pairs\n(Sentinel-1A + 1B)", "2022–2024": "2022–24 · 12-day pairs\n(Sentinel-1A only)",
          "2025": "2025 · 6-day pairs\n(Sentinel-1A + 1C)"}


def s2_cells(cells, bs):
    """S2 — what the nine cells are: land cover and position, then backscatter over all the laser years."""
    def bsp(track, col):
        return bs[(bs.track == track) & (bs.years == "2021–2025")].set_index("cell")[col].to_dict()
    top = [("Shrub cover 0.5–2 m (%)", (cells.shrub_share_2024 * 100).to_dict(), ".0f", "YlGn", 0, 20),
           ("Boardwalk share (%)", (cells.boardwalk_share * 100).to_dict(), ".0f", "Greys", 0, 20),
           ("Depth into the mat (m)", cells.depth_in_mat_m.to_dict(), ".0f", "Blues", 150, 280),
           ("Canopy height, 90th pct (m)", cells.canopy2024_p90_m.to_dict(), ".2f", "YlGn", 0, 1.5)]
    bot = [("VV (dB)", bsp("ascending", "vv_median_db"), ".1f", "magma", -14, -9),
           ("VH (dB)", bsp("ascending", "vh_median_db"), ".1f", "magma", -20, -15),
           ("VV (dB)", bsp("descending", "vv_median_db"), ".1f", "magma", -14, -9),
           ("VH (dB)", bsp("descending", "vh_median_db"), ".1f", "magma", -20, -15)]
    with fs.style():
        fig = plt.figure(figsize=fs.size("full", 105), layout="constrained")
        fig.suptitle(READ, x=0.01, ha="left", fontsize=7, color=fs.MUTED)
        sf = fig.subfigures(2, 1)
        heads = ["The cells: laser scan and orthophoto of 19 Sep 2024 (no cell holds trees ≥ 2 m: at most 1 %)",
                 "Backscatter, median of all scenes 2021–2025 · left pair: ascending (dusk, ≈ 16:36 UTC) · right "
                 "pair: descending (dawn, ≈ 05:09 UTC)"]
        for k, (row, head) in enumerate(zip([top, bot], heads)):
            axs = sf[k].subplots(1, 4)
            sf[k].suptitle(head, x=0.01, ha="left", fontsize=7.5, fontweight="bold")
            for j, (ax, (title, vals, fmt, cmap, lo, hi)) in enumerate(zip(axs, row)):
                if k == 1:
                    title = f"{'Ascending' if j < 2 else 'Descending'} · {title}"
                tiles(ax, vals, fmt, cmap, lo, hi, title)
    data = [{"panel": t, "cell": c, "value": v} for row in (top, bot) for (t, vals, *_r) in row for c, v in vals.items()]
    export(fig, "S2_nine_cells_what_they_are", ["hub:supervisor_pack_c074/tables/nine_cells.csv",
                                                "nine_cells_backscatter.csv (years 2021–2025)"],
           {"layout": READ, "scales": "fixed: VV −14…−9 dB, VH −20…−15 dB; shrub and boardwalk 0–20 %"},
           {"tiles": pd.DataFrame(data)})


def s2c_coherence_laser(by):
    """S2c — per pair set (rows) and track: coherence of the pairs and how the cell's phase follows the laser."""
    cols = [("Coherence", "ascending", "median_coh", "viridis", 0.5, 1),
            ("Coherence", "descending", "median_coh", "viridis", 0.5, 1),
            ("r with the P6 laser", "ascending", "pearson_increments", "cividis", 0, 1),
            ("r with the P6 laser", "descending", "pearson_increments", "cividis", 0, 1)]
    with fs.style():
        fig = plt.figure(figsize=fs.size("full", 150), layout="constrained")
        fig.suptitle(READ + " Rows: the three pair sets of the laser years. Columns: median coherence of those "
                     "pairs, and Pearson r of each cell's phase change with the P6 laser's change, per track.",
                     x=0.01, ha="left", fontsize=7, color=fs.MUTED, wrap=True)
        gs = fig.add_gridspec(3, 4)
        for i, g in enumerate(GROUPS):
            n = {tr: int(by[(by.cell == "P6") & (by.track == tr) & (by.group == g)].pairs.iloc[0]) for tr in TRACKS}
            for j, (q, tr, col, cmap, lo, hi) in enumerate(cols):
                ax = fig.add_subplot(gs[i, j])
                d = by[(by.track == tr) & (by.group == g)].set_index("cell")[col].to_dict()
                tiles(ax, d, ".2f", cmap, lo, hi, f"{q} · {tr[:3]}." if i == 0 else "")
                if j == 0:
                    ax.set_ylabel(f"{ROWLAB[g]}\nn = {n['ascending']} asc. / {n['descending']} desc.", fontsize=7.5,
                                  fontweight="bold", rotation=0, ha="right", va="center", labelpad=8)
    export(fig, "S2c_nine_cells_coherence_laser", ["hub:supervisor_pack_c074/tables/nine_cells_by_group.csv"],
           {"layout": "rows: 2021 6-day, 2022–24 12-day, 2025 6-day consecutive pairs; columns: coherence asc/desc, "
                      "r with the laser asc/desc", "scales": "coherence 0.5…1, r 0…1"},
           {"tiles": by[by.cell.isin(LAYOUT[0] + LAYOUT[1] + LAYOUT[2])]})


def s2b_forest(cells, by):
    order = cells.sort_values("depth_in_mat_m").index.tolist() + ["mean of the 8 neighbours"]
    lab = {c: (f"{c} · {cells.loc[c, 'depth_in_mat_m']:.0f} m · shrub {cells.loc[c, 'shrub_share_2024']:.0%}"
               if c in cells.index else "mean of 8 neighbours") for c in order}
    with fs.style():
        fig, axs = plt.subplots(2, 3, figsize=fs.size("full", 124), layout="constrained", sharey=True, sharex="col")
        fig.suptitle("How strongly each cell's phase follows the P6 laser (slope; 1 = the laser's own amplitude). Rows "
                     "of each panel: the nine cells, shallow to deep in the mat. Where the intervals overlap, the cells "
                     "cannot be ranked.", x=0.01, ha="left", fontsize=7, color=fs.MUTED)
        for i, tr in enumerate(TRACKS):
            for j, g in enumerate(GROUPS):
                ax = axs[i, j]
                d = by[(by.track == tr) & (by.group == g)].set_index("cell").loc[order]
                y = np.arange(len(order))
                p6 = float(d.loc["P6", "slope_through_origin"])
                ax.axvline(p6, color=P6C, lw=0.7, ls=":")
                ax.axvline(1, color=fs.MUTED, lw=0.5, ls="--")
                for yy, (c, r) in zip(y, d.iterrows()):
                    col = P6C if c == "P6" else (NBC if c.startswith("mean") else fs.INK)
                    mk = "*" if c == "P6" else ("D" if c.startswith("mean") else "o")
                    ax.plot([r.slope_lo95, r.slope_hi95], [yy, yy], color=col, lw=1.0)
                    ax.plot(r.slope_through_origin, yy, mk, ms=6 if mk == "*" else 4, color=col)
                ax.set_yticks(y, [lab[c] for c in order])
                ax.invert_yaxis()
                if i == 1:
                    ax.set_xlabel("Slope through origin vs laser (95 % CI)")
                ax.set_title(f"{tr.capitalize()} · {GLABEL[g]} (n = {int(d.pairs.iloc[0])} pairs)", loc="left",
                             fontsize=7)
        fig.legend(handles=[Line2D([], [], color=P6C, marker="*", ms=6, label="P6 cell"),
                            Line2D([], [], color=fs.INK, marker="o", ms=4, label="neighbour cell"),
                            Line2D([], [], color=NBC, marker="D", ms=4, label="mean of the 8 neighbours"),
                            Line2D([], [], color=P6C, lw=0.7, ls=":", label="P6's slope"),
                            Line2D([], [], color=fs.MUTED, lw=0.5, ls="--", label="slope 1 (the laser's amplitude)")],
                   loc="outside lower center", ncol=5, fontsize=6.5)
    export(fig, "S2b_nine_cells_slopes", ["hub:supervisor_pack_c074/tables/nine_cells_by_group.csv"],
           {"rows": "cells ordered by depth into the mat (shallow on top); label: depth · shrub cover 2024",
            "intervals": "2000 resamples of the pairs, seed 0 (C-074)"},
           {"slopes": by[by.cell.isin(order)]})


# ------------------------------------------------------------------------------ S3 — phase against the laser

def s3_scatter(pairs, by):
    nb = pairs[pairs.cell != "P6"].groupby(["track", "group", "pair"]).s1_los_mm.mean().rename("nb")
    p6 = pairs[pairs.cell == "P6"].set_index(["track", "group", "pair"])
    d = p6.join(nb).reset_index()
    d = d[pd.to_datetime(d.t2, format="ISO8601", utc=True) < pd.Timestamp("2026-01-01", tz="UTC")]  # design years
    lim = np.ceil(np.nanmax(np.abs(d.dropna(subset=["laser_los_mm"])[["laser_los_mm", "s1_los_mm", "nb"]]
                                  .to_numpy())) / 5) * 5
    with fs.style():
        fig, axs = plt.subplots(2, 3, figsize=fs.size("full", 132), layout="constrained", sharex=True, sharey=True)
        fig.suptitle("Each point is one consecutive pair: the laser's change at P6 (x) against the radar's change (y). "
                     "On the dashed line the radar reads the laser exactly; outside the blue band the laser moved more "
                     "than one pair can read without ambiguity.", x=0.01, ha="left", fontsize=7, color=fs.MUTED)
        for i, tr in enumerate(TRACKS):
            for j, g in enumerate(GROUPS):
                ax = axs[i, j]
                e = d[(d.track == tr) & (d.group == g)].dropna(subset=["laser_los_mm"])
                ax.axhspan(-QW, QW, color="#dbe9f6", zorder=0, lw=0)
                ax.plot([-lim, lim], [-lim, lim], color=fs.INK, lw=0.7, ls="--", zorder=1)
                ax.axhline(0, color=fs.GRID, lw=0.5)
                ax.axvline(0, color=fs.GRID, lw=0.5)
                for col, cname, mk, colr, lab in (("s1_los_mm", "P6", "o", P6C, "P6 cell"),
                                                  ("nb", "mean of the 8 neighbours", "^", NBC, "mean of 8 neighbours")):
                    st = by[(by.track == tr) & (by.group == g) & (by.cell == cname)].iloc[0]
                    ax.scatter(e.laser_los_mm, e[col], s=14, marker=mk, facecolors=colr if mk == "o" else "white",
                               edgecolors=colr, linewidths=0.8, zorder=3,
                               label=f"{lab}\nr {st.pearson_increments:.2f} [{st.r_lo95:.2f}, {st.r_hi95:.2f}]; "
                                     f"slope {st.slope_through_origin:.2f} [{st.slope_lo95:.2f}, {st.slope_hi95:.2f}]")
                    xs = np.array([-lim, lim])
                    ax.plot(xs, st.slope_through_origin * xs, color=colr, lw=1.0, ls="-" if mk == "o" else "-.",
                            zorder=2)
                ax.set_xlim(-lim, lim)
                ax.set_ylim(-lim, lim)
                ax.set_aspect("equal")
                ax.set_title(f"{tr.capitalize()} · {GLABEL[g]} · n = {len(e)}", loc="left", fontsize=7)
                ax.legend(loc="upper left", fontsize=5.6, handletextpad=0.3, borderaxespad=0.2, labelspacing=0.3)
                if i == 1:
                    ax.set_xlabel("Laser change along the line of sight (mm)")
                if j == 0:
                    ax.set_ylabel("Sentinel-1 change (mm, LOS)")
        fig.legend(handles=[Line2D([], [], color=fs.INK, lw=0.7, ls="--", label="1 : 1"),
                            Rectangle((0, 0), 1, 1, fc="#dbe9f6", ec="none",
                                      label=f"± λ/4 = ± {QW:.1f} mm (one pair's unambiguous range)"),
                            Line2D([], [], color=P6C, lw=1.0, label="slope through origin, P6"),
                            Line2D([], [], color=NBC, lw=1.0, ls="-.", label="slope through origin, neighbours")],
                   loc="outside lower center", ncol=4, fontsize=6.5)
    export(fig, "S3_phase_vs_laser", ["hub:supervisor_pack_c074/tables/nine_cells_pairs.csv", "nine_cells_by_group.csv"],
           {"pairs": "consecutive pairs as X-126 (hard-target reference); laser on the snow-evidence rule",
            "statistics": "r and slope through the origin with 2000-resample pair-bootstrap intervals"},
           {"pairs": d})


# ------------------------------------------------------------------------------ S4 — P6 through time

def s4_time(pairs):
    s = pd.read_csv(CUBE / "silver" / "site_hourly.csv.gz", usecols=["time_utc", "laser_surface_cm", "laser_ok",
                                                                     "wtd_P6_cm"], index_col=0, low_memory=False)
    s.index = pd.to_datetime(s.index, format="ISO8601", utc=True)
    ok = s.laser_ok.map({True: True, False: False, "True": True, "False": False}).fillna(False)
    las = (s.laser_surface_cm.where(ok) * 10).resample("1D").mean()
    wtd = s.wtd_P6_cm.resample("1D").mean()
    t0, t1 = pd.Timestamp("2021-05-01", tz="UTC"), pd.Timestamp("2025-12-31", tz="UTC")
    las, wtd = las[t0:t1], wtd[t0:t1]
    gap = (pd.Timestamp("2022-07-01", tz="UTC"), pd.Timestamp("2022-09-01", tz="UTC"))
    p = pairs[pairs.cell == "P6"].copy()
    for c in ("t1", "t2"):
        p[c] = pd.to_datetime(p[c], format="ISO8601", utc=True)
    p = p[p.t2 < pd.Timestamp("2026-01-01", tz="UTC")]          # pairs ending in 2026 (held-out year) left out
    p["mid"] = p.t1 + (p.t2 - p.t1) / 2
    with fs.style():
        fig, axs = plt.subplots(5, 1, figsize=fs.size("full", 175), layout="constrained", sharex=True,
                                gridspec_kw={"height_ratios": [1, 0.8, 1.1, 1.1, 0.8]})
        ax = axs[0]
        before, after = las[:gap[0]], las[gap[1]:]
        ax.plot(before.index, before - np.nanmedian(before), color=fs.CLASS["laser"], lw=0.9)
        ax.plot(after.index, after - np.nanmedian(after), color=fs.CLASS["laser"], lw=0.9, ls="-")
        ax.axvspan(*gap, color=fs.GRID, alpha=0.6, lw=0)
        ax.text(gap[0], 0.97, " laser re-levelled in this gap:\n levels before and after on different datums",
                transform=ax.get_xaxis_transform(), fontsize=6.2, va="top", color=fs.MUTED)
        ax.set_ylabel("Laser surface\n(mm, each datum\nabout its median)")
        fig.suptitle("P6: the mat surface (laser), the water table, and Sentinel-1", x=0.01, ha="left", fontsize=8)
        fs.panel_label(ax, "a")
        ax = axs[1]
        ax.plot(wtd.index, wtd, color=fs.CLASS["open water"], lw=0.9)
        ax.set_ylabel("Water-table\ndepth P6 (cm)")
        ax.text(0.995, 0.06, "2025: not yet delivered", transform=ax.transAxes, ha="right", fontsize=6.3,
                color=fs.MUTED)
        fs.panel_label(ax, "b")
        for k, tr in enumerate(TRACKS):
            ax = axs[2 + k]
            e = p[p.track == tr]
            ax.axhline(0, color=fs.GRID, lw=0.5)
            ax.scatter(e.mid, e.laser_los_mm, s=9, marker="*", color=fs.INK, label="laser change", zorder=3)
            ax.scatter(e.mid, e.s1_los_mm, s=9, marker="o" if tr == "ascending" else "s", facecolors="none",
                       edgecolors=fs.TRACK[tr], linewidths=0.8, label=f"Sentinel-1 change, {tr}", zorder=2)
            ax.set_ylabel(f"Pair change,\n{tr} (mm LOS)")
            ax.text(1.0, 1.01, "radar point without a star: no laser value (snow evidence or a gap)",
                    transform=ax.transAxes, ha="right", va="bottom", fontsize=6.2, color=fs.MUTED)
            ax.legend(loc="lower left", bbox_to_anchor=(0.0, 1.0), fontsize=6.2, ncol=2, borderaxespad=0.1)
            fs.panel_label(ax, "cd"[k], x=-0.06)
        ax = axs[4]
        for tr in TRACKS:
            e = p[p.track == tr]
            ax.scatter(e.mid, e.coherence, s=8, marker="o" if tr == "ascending" else "s", color=fs.TRACK[tr],
                       label=tr, alpha=0.8, linewidths=0)
        ax.set_ylim(0, 1)
        ax.set_ylabel("Coherence\n(P6 cell)")
        ax.legend(loc="lower left", bbox_to_anchor=(0.0, 1.0), fontsize=6.2, ncol=2, borderaxespad=0.1)
        fs.panel_label(ax, "e", x=-0.06)
        for g, lab in (("2021", "6-day pairs"), ("2022–2024", "12-day pairs (Sentinel-1B lost)"), ("2025", "6-day pairs (S1A + S1C)")):
            e = p[p.group == g]
            if len(e):
                axs[0].annotate(lab, (e.mid.min() + (e.mid.max() - e.mid.min()) / 2, 1.0), xycoords=("data", "axes fraction"),
                                ha="center", va="bottom", fontsize=6.2, color=fs.MUTED, annotation_clip=False)
    export(fig, "S4_p6_through_time", ["cube:silver/site_hourly.csv.gz (laser_surface_cm, laser_ok, wtd_P6_cm)",
                                       "hub:supervisor_pack_c074/tables/nine_cells_pairs.csv (P6)"],
           {"display": "daily means; gaps left open; pair changes plotted at the pair's middle date, never joined",
            "laser": "snow-evidence rule; the July–August 2022 re-levelling splits the level into two datums",
            "left_out": "pairs ending in 2026 (held-out year), among them 26 Dec 2025 – 1 Jan 2026, which shifts all nine "
                        "cells by about +218 mm against the hard targets (no laser value)"},
           {"p6_pairs": p.drop(columns=["mid"])})


# ------------------------------------------------------------------------------ the C-073 figures, beside

COPIES = {"F2_trees_vs_open": "S5_trees_vs_open_1km", "F4_footprint": "S6_footprint",
          "F5_gain_vs_depth": "S7_gain_vs_depth", "F6_season_maps": "S8_season_maps",
          "F3_drivers_of_coherence": "S9_drivers_of_coherence", "F1_dated_land_cover": "S10_dated_land_cover"}


def copy_c073() -> None:
    for a, b in COPIES.items():
        for ext in ("pdf", "png", "json"):
            shutil.copy2(C073 / f"{a}.{ext}", PACK / f"{b}.{ext}")


def readme(cells, by, bs) -> None:
    def g(cell, tr, grp, col):
        return float(by[(by.cell == cell) & (by.track == tr) & (by.group == grp)][col].iloc[0])
    sh = cells.shrub_share_2024
    tr_max = cells.tree_share_2024.max()
    a21 = {c: g(c, "ascending", "2021", "slope_through_origin") for c in cells.index}
    lo_c, hi_c = min(a21, key=a21.get), max(a21, key=a21.get)
    vv = bs[(bs.track == "ascending") & (bs.years == "2022–2024")].vv_median_db
    figs = [
        ("S1_nine_cells_site", "The nine cells around P6, 19 September 2024",
         f"(a) The 10 cm true-colour orthophoto of 19 Sep 2024 (GUGiK) with the nine 40 m radar cells around P6 on the "
         f"corrected grid (white; P6's cell orange) and the boardwalk share of each cell; star: the P6 station. (b) "
         f"Vegetation height from the national laser scan flown the same day (1 m), with each cell's shrub (0.5–2 m) and "
         f"tree (≥ 2 m) cover. The nine cells hold no trees (at most {tr_max:.0%} of a cell); shrub cover runs from "
         f"{sh.min():.0%} to {sh.max():.0%}, highest to the north and east. P6's cell is the one the boardwalk and the "
         f"plot platforms cross ({cells.loc['P6', 'boardwalk_share']:.0%} of it).",
         "Two square images of the same 160 m area. Left, an aerial photograph of the open fen with the boardwalk "
         "running north–south and the instrument platforms, overlaid by a three-by-three grid of radar cells. Right, a "
         "pale vegetation-height map with scattered green shrub patches, densest in the north-east cells."),
        ("S2_nine_cells_what_they_are", "What the nine cells are",
         f"Each 3 × 3 block is the nine 40 m cells as on the ground (north up; P6 framed). Top: shrub cover "
         f"(0.5–2 m) and canopy height from the laser scan of 19 Sep 2024, the boardwalk's share of each cell and "
         f"its depth into the mat. No cell holds trees (≥ 2 m on at most {tr_max:.0%} of a cell). Bottom: median "
         f"backscatter of all scenes 2021–2025, VV and VH, ascending (dusk) and descending (dawn); the nine cells "
         f"differ by about a decibel. Fixed colour scales.",
         "Eight small three-by-three grids. Shrub cover and depth rise towards the north-east; the boardwalk "
         "crosses P6's cell and the one south of it; backscatter barely changes from cell to cell on either track."),
        ("S2c_nine_cells_coherence_laser", "Coherence and the laser, cell by cell, in each pair set",
         "Rows: the three pair sets of the laser years — six-day pairs of 2021 (Sentinel-1A + 1B), twelve-day pairs of "
         "2022–24 (Sentinel-1A alone), six-day pairs of 2025 (Sentinel-1A + 1C). Columns: the median coherence of "
         "those pairs, then Pearson r of each cell's phase change with the P6 laser's change, ascending and "
         "descending (n = pairs). The rows differ far more than the cells: the pair set decides how well the phase "
         "follows the laser, not the cell.",
         "Twelve three-by-three grids in three rows. In the 2021 and 2025 rows coherence and correlation are high in "
         "every cell; in the 2022–24 row they drop, most on the descending track, again in every cell alike."),
        ("S2b_nine_cells_slopes", "The slope against the laser, cell by cell",
         f"Slope through the origin of each cell's Sentinel-1 change against the P6 laser's, with 95 % intervals (2000 "
         f"resamples of the pairs), per track and pair set; cells ordered by depth into the mat, labelled with depth and "
         f"shrub cover; P6 (star), the mean of its eight neighbours (diamond). In 2021 ascending the slope runs from "
         f"{a21[lo_c]:.2f} ({lo_c}) to {a21[hi_c]:.2f} ({hi_c}); the intervals overlap, and the cells that read lower "
         f"are the shallow ones with little shrub, the higher the deep ones with more: shrub and depth rise together "
         f"across nine cells and cannot be separated here; S5 (trees, matched over 1 km) and S7 (depth, over the whole mat) take each apart.",
         "Six forest plots of nine cells each. In every panel the intervals overlap widely; slopes are close to one "
         "on the six-day pairs and about half on the twelve-day pairs, lowest of all on the descending twelve-day pairs."),
        ("S3_phase_vs_laser", "Sentinel-1 against the laser at P6 and around it",
         f"Each point is one consecutive pair: the laser's change at P6 along the line of sight against the Sentinel-1 "
         f"change of P6's cell (filled) and of the mean of its eight neighbours (open), hard-target reference. Lines: "
         f"slope through the origin; dashed: 1 : 1; band: ± λ/4 = ± {QW:.1f} mm, the largest change one pair can read "
         f"without ambiguity. r and slope with 95 % intervals in each legend. Six-day pairs (2021, 2025) follow the "
         f"laser on both tracks; twelve-day pairs (2022–24, Sentinel-1A alone) follow it partly on the ascending track "
         f"and weakly on the descending.",
         "Six scatter plots, ascending on top and descending below, by pair set. In the four six-day panels the points "
         "lie along the one-to-one line for P6 and its neighbours alike; in the two twelve-day panels they scatter "
         "widely and the fitted slopes are well below one."),
        ("S4_p6_through_time", "P6 through time: laser, water table, Sentinel-1, coherence",
         "(a) The laser's surface level, daily (the snow-evidence rule; re-levelled in the July–August 2022 gap, so the "
         "levels before and after are on different datums and are drawn each about its own median). (b) The P6 "
         "water-table depth (logger, daily; delivered to the end of 2024). (c, d) Each consecutive pair's change: laser "
         "(stars) and Sentinel-1 at P6's cell (open markers), ascending and descending, at the pair's middle date, never "
         "joined. (e) The coherence of those pairs. Pairs ending in 2026 are left out (a held-out year).",
         "Five stacked time panels from mid-2021 to the end of 2025. The laser and water table rise and fall together "
         "each year. In the six-day years the radar's pair changes sit on the laser's; in 2022–24 they scatter more and "
         "coherence varies more."),
        ("S5_trees_vs_open_1km", "Trees against open ground within 1 km (matched)",
         "From C-073 F2 (X-132): wooded − open coherence, on and off the mat, matched on depth and distance; the "
         "difference between the two; the effect on how the phase follows the laser.", "See C-073 F2."),
        ("S6_footprint", "How far P6's six-day motion reaches", "From C-073 F4 (X-133).", "See C-073 F4."),
        ("S7_gain_vs_depth", "The motion grows towards the mat's interior", "From C-073 F5 (X-133).", "See C-073 F5."),
        ("S8_season_maps", "A first season map of the mat", "From C-073 F6 (X-133).", "See C-073 F6."),
        ("S9_drivers_of_coherence", "What drives coherence", "From C-073 F3 (X-132).", "See C-073 F3."),
        ("S10_dated_land_cover", "Dated land cover around P6", "From C-073 F1 (C-072).", "See C-073 F1."),
    ]
    lines = ["# Supervisor pack, October 2026 — figures and tables (generated)", "",
             "Generated by SAr-adapted `experiments/supervisor_figures_c074.py` from the tables of "
             "`experiments/p6_nine_cells_c074.py` (`tables/`, pre-registered: `../insar_full_test/preregistration/"
             "C-074.yaml`) and the C-073 figures (`../figures_landcover_footprint/`, whose README holds the full "
             "captions of S5–S10). Do not edit.", "",
             "Inputs: the corrected radar grids (HyP3 half-cell, RTC placement), the P6 laser on the snow-evidence rule, "
             "the hard-target reference, the dated GUGiK laser scans and orthophotos. Design years only. Each figure: "
             "vector PDF (Arial embedded), PNG, provenance JSON, and `data/` CSVs of what it plots.", "",
             "These replace the figures shown in early October (wooden-boardwalk / wood-free neighbour set), which were "
             "withdrawn on 8 October (X-079): their per-cell numbers were real, but the conclusion — that wood-free "
             "neighbours follow the laser 1 : 1 and woody cover or the boardwalk changes the slope — did not hold, one "
             "series labelled 'wood-free neighbour' was the mean of the west and east cells, and the slopes were TLS "
             "slopes, which noise pushes up.", ""]
    for name, title, cap, alt in figs:
        lines += [f"## {name} — {title}", "", f"![{title}]({name}.png)", "", f"**Caption.** {cap}", "",
                  f"**Alt text.** {alt}", ""]
    (PACK / "README.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    PACK.mkdir(parents=True, exist_ok=True)
    cells = pd.read_csv(TAB / "nine_cells.csv").set_index("cell")
    by = pd.read_csv(TAB / "nine_cells_by_group.csv")
    bs = pd.read_csv(TAB / "nine_cells_backscatter.csv")
    pairs = pd.read_csv(TAB / "nine_cells_pairs.csv")
    s1_site(cells)
    s2_cells(cells, bs)
    s2c_coherence_laser(by)
    s2b_forest(cells, by)
    s3_scatter(pairs, by)
    s4_time(pairs)
    copy_c073()
    readme(cells, by, bs)
    print("pack figures written to", PACK)


if __name__ == "__main__":
    main()

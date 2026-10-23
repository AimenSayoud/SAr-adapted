"""C-073 — the figures of the dated land cover (C-072), trees against open ground (X-132) and the footprint and season
map of P6's six-day motion (X-133).

Drawn only from the committed outputs of X-132 and X-133 and the cube's static layers and dated orthophotos; nothing
is re-estimated except descriptive summaries that are written beside each figure (binned medians, the P6 uncertainty
over the laser-covered pairs). House style (``insar_wetlands.figstyle``): Elsevier widths, ≥ 7 pt text, one colour
per class and per track, and every colour paired with a marker, line style or fill so no reading rests on colour alone.
Each figure is written as a vector PDF, a PNG preview, a provenance record and the CSV of what it plots.

Writes only to the hub (``08_deliverables/figures_landcover_footprint/``).

    make run E=figures_landcover_footprint_c073
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import xarray as xr  # noqa: E402
from matplotlib.colors import ListedColormap, TwoSlopeNorm  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import FancyArrowPatch, Patch  # noqa: E402
from matplotlib_scalebar.scalebar import ScaleBar  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from article_fig_f1 import cell_outline  # noqa: E402  (the 40 m outline along cell edges, as in the article's F1)
from insar_wetlands import figstyle as fs  # noqa: E402
from insar_wetlands.paths import cube_root, hub_root  # noqa: E402

HUB, CUBE = hub_root(), cube_root()
FT = HUB / "08_deliverables" / "insar_full_test"
X132, X133 = FT / "landcover_x132", FT / "footprint_x133"
OUT = HUB / "08_deliverables" / "figures_landcover_footprint"
CELL = 40.0
TRACKS = ("ascending", "descending")
TMARK = {"ascending": "o", "descending": "s"}
TLW = {"ascending": 1.5, "descending": 0.8}         # the two track hues are close in grey: line weight carries track too
SETS = {"2020_2021": "2020–21, 6-day", "2022_2024": "2022–24, 12-day", "2025": "2025, 6-day"}
SIDE = {"on_mat": ("on the mat", fs.CLASS["mat"], "o", "-"), "off_mat": ("off the mat", fs.MUTED, "s", "--")}
BAD = "#eeeeee"
FACTOR_LABEL = {"ndvi_season": "Season (cell NDVI of the month)", "descending": "Descending track (vs ascending)",
                "pair_days": "Pair length", "wind_ms": "Wind at the overpasses", "frost_snow": "Frost or snow (either date)",
                "abs_dwtd_cm": "|Water-table change|", "on_mat": "On the mat (vs off)", "rain_mm": "Rain between dates",
                "canopy2024_p90_m": "Canopy height (2024, p90)", "boardwalk_share": "Boardwalk share",
                "dwtd_x_on_mat": "|Water-table change| × on mat", "wc_tree": "WorldCover tree share",
                "canopy_x_on_mat": "Canopy height × on mat", "canopy2024_cover": "Canopy cover (2024, ≥ 2 m)",
                "depth_in_mat_m": "Depth into the mat"}


# ------------------------------------------------------------------------------ shared map pieces

def static() -> xr.Dataset:
    return xr.open_dataset(CUBE / "silver" / "static_40m.nc")


def grid_of(st: xr.Dataset, rows, cols, values) -> np.ndarray:
    a = np.full(st.zone.shape, np.nan)
    a[np.asarray(rows), np.asarray(cols)] = np.asarray(values, float)
    return a


def extent(st) -> list[float]:
    x, y = st.x.values, st.y.values
    return [x[0] - CELL / 2, x[-1] + CELL / 2, y[-1] - CELL / 2, y[0] + CELL / 2]


def map_axes(ax, xlim, ylim, scalebar: bool = True, north: bool = True, ylabel: bool = True) -> None:
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal")
    ticks = lambda lo, hi: np.arange(np.ceil(lo / 500) * 500, hi, 500)  # noqa: E731
    ax.set_xticks(ticks(*xlim))
    ax.set_yticks(ticks(*ylim))
    ax.xaxis.set_major_formatter(lambda v, _: f"{v / 1000:.1f}")
    ax.yaxis.set_major_formatter(lambda v, _: f"{v / 1000:.1f}")
    ax.set_xlabel("Easting (km, UTM 33N)")
    if ylabel:
        ax.set_ylabel("Northing (km, UTM 33N)")
    else:
        ax.set_yticklabels([])
    for s in ("top", "right"):
        ax.spines[s].set_visible(True)
    if scalebar:
        ax.add_artist(ScaleBar(1, "m", location="lower left", length_fraction=0.25, box_alpha=0.85,
                               font_properties={"size": 7}, border_pad=0.5))
    if north:
        x0, x1 = ax.get_xlim()
        y0, y1 = ax.get_ylim()
        nx, ny = x1 - 0.07 * (x1 - x0), y1 - 0.17 * (y1 - y0)
        ax.add_patch(FancyArrowPatch((nx, ny), (nx, ny + 0.08 * (y1 - y0)), arrowstyle="-|>", mutation_scale=8,
                                     color=fs.INK, lw=0.8, zorder=7))
        ax.text(nx, ny + 0.09 * (y1 - y0), "N", ha="center", va="bottom", fontsize=7, zorder=7)


def mat_outline(ax, st, **kw) -> None:
    kw = {"colors": fs.CLASS["mat"], "linewidths": 1.3, "zorder": 5, **kw}
    cell_outline(ax, st.x.values, st.y.values, st.zone.values == 1, **kw)


def p6_star(ax, x, y) -> None:
    ax.plot(x, y, "*", ms=9, mfc=fs.CLASS["laser"], mec="white", mew=0.6, zorder=8)


def box_around(cells: pd.DataFrame, pad: float = 60.0):
    return (cells.x.min() - CELL / 2 - pad, cells.x.max() + CELL / 2 + pad), \
           (cells.y.min() - CELL / 2 - pad, cells.y.max() + CELL / 2 + pad)


def cir_image(xlim, ylim, res: float = 2.0) -> np.ndarray:
    """The 2022-07-19 colour-infrared orthophoto over the box, reprojected to UTM 33N at ``res`` m (band means),
    shown as delivered (near-infrared, red, green as R, G, B; no contrast stretch)."""
    import rasterio
    from pyproj import Transformer
    from rasterio.merge import merge
    from rasterio.transform import from_origin
    from rasterio.warp import Resampling, reproject
    raw = CUBE / "raw" / "gugik_dated"
    tab = pd.read_csv(raw / "sheets.csv")
    fs_ = [raw / f for f in tab[(tab.kind == "cir") & (tab.year == 2022)].file]
    t = Transformer.from_crs(32633, 2180, always_xy=True)
    ex, ey = t.transform([xlim[0], xlim[1], xlim[0], xlim[1]], [ylim[0], ylim[0], ylim[1], ylim[1]])
    b = (min(ex) - 20, min(ey) - 20, max(ex) + 20, max(ey) + 20)
    srcs = [rasterio.open(f) for f in fs_]
    img, tr = merge(srcs, bounds=b, res=1.0, resampling=Resampling.average)
    crs = srcs[0].crs
    for s in srcs:
        s.close()
    W, H = int((xlim[1] - xlim[0]) / res), int((ylim[1] - ylim[0]) / res)
    out = np.zeros((3, H, W), "uint8")
    for k in range(3):
        reproject(img[k], out[k], src_transform=tr, src_crs=crs, dst_transform=from_origin(xlim[0], ylim[1], res, res),
                  dst_crs="EPSG:32633", resampling=Resampling.average)
    return np.moveaxis(out, 0, -1)


def export(fig, name: str, sources: list[str], notes: dict, data: dict[str, pd.DataFrame]) -> None:
    for k, d in data.items():
        d.to_csv(OUT / "data" / f"{name}_{k}.csv", index=False)
    with fs.style():                         # text is laid out at save time: fonts and Type 42 must be in force here
        fs.export(fig, OUT / name, provenance={"sources": sources, "data": [f"data/{name}_{k}.csv" for k in data],
                                               **notes})
    plt.close(fig)


# ------------------------------------------------------------------------------ F1 — dated land cover

def f1(st, cells, p6):
    xlim, ylim = box_around(cells)
    img = cir_image(xlim, ylim)
    c24 = st.canopy2024_p90_m.values
    dch = st.canopy2024_p90_m.values - st.canopy2011_p90_m.values
    inside = grid_of(st, cells.row, cells["col"], np.ones(len(cells)))
    with fs.style():
        fig, axs = plt.subplots(1, 3, figsize=fs.size("full", 92), layout="constrained")
        ax = axs[0]
        ax.imshow(img, extent=[xlim[0], xlim[1], ylim[0], ylim[1]], origin="upper", interpolation="none", zorder=0)
        w, o = cells[cells.cls == "wooded"], cells[cells.cls == "open"]
        ax.scatter(w.x, w.y, s=4, marker="s", c=fs.INK, linewidths=0, zorder=3, label="wooded (p90 ≥ 2 m)")
        ax.scatter(o.x, o.y, s=4, marker="s", facecolors="white", edgecolors=fs.INK, linewidths=0.3, zorder=3,
                   label="open (p90 < 0.5 m)")
        mat_outline(ax, st)
        p6_star(ax, *p6)
        map_axes(ax, xlim, ylim)
        ax.set_title("Orthophoto (CIR), 19 Jul 2022", loc="left")
        ax.legend(handles=[Line2D([], [], ls="", marker="s", ms=3.5, mfc=fs.INK, mec=fs.INK, label="wooded (2024 scan)"),
                           Line2D([], [], ls="", marker="s", ms=3.5, mfc="white", mec=fs.INK, label="open (2024 scan)"),
                           Line2D([], [], color=fs.CLASS["mat"], lw=1.3, label="floating mat"),
                           Line2D([], [], ls="", marker="*", ms=7, mfc=fs.CLASS["laser"], mec="white", label="P6 (laser)")],
                  loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=2, fontsize=6.5, handletextpad=0.3,
                  columnspacing=0.8)
        fs.panel_label(ax, "a")
        ax = axs[1]
        im = ax.imshow(np.where(inside == 1, c24, np.nan), extent=extent(st), origin="upper", interpolation="nearest",
                       cmap=plt.get_cmap("YlGn").with_extremes(bad=BAD), vmin=0, vmax=25, zorder=1)
        mat_outline(ax, st)
        p6_star(ax, *p6)
        map_axes(ax, xlim, ylim, scalebar=False, ylabel=False)
        ax.set_title("Canopy height, 19 Sep 2024", loc="left")
        fig.colorbar(im, cax=ax.inset_axes([0.08, -0.3, 0.84, 0.045]), orientation="horizontal", extend="max",
                     label="Canopy height, cell 90th percentile (m)")
        fs.panel_label(ax, "b")
        ax = axs[2]
        im = ax.imshow(np.where(inside == 1, dch, np.nan), extent=extent(st), origin="upper", interpolation="nearest",
                       cmap=plt.get_cmap("BrBG").with_extremes(bad=BAD), norm=TwoSlopeNorm(0, -10, 10), zorder=1)
        mat_outline(ax, st, colors=fs.INK, linewidths=1.0)
        p6_star(ax, *p6)
        map_axes(ax, xlim, ylim, scalebar=False, north=False, ylabel=False)
        ax.set_title("Change, Apr 2011 → Sep 2024", loc="left")
        fig.colorbar(im, cax=ax.inset_axes([0.08, -0.3, 0.84, 0.045]), orientation="horizontal", extend="both",
                     label="Δ canopy height p90 (m)")
        fs.panel_label(ax, "c")
    tr = pd.crosstab(cells[cells.on_mat == 1].cls_2011, cells[cells.on_mat == 1].cls)
    export(fig, "F1_dated_land_cover",
           ["cube:raw/gugik_dated (CIR 2022)", "cube:silver/static_40m.nc canopy2024_p90_m, canopy2011_p90_m",
            "hub:landcover_x132/cells.csv"],
           {"display": "orthophoto as delivered (NIR, R, G as RGB), no stretch; cells outside the X-132 sample grey"},
           {"cells": cells.reset_index()[["cell", "x", "y", "on_mat", "cls", "cls_2011", "canopy2024_p90_m",
                                          "canopy2011_p90_m"]],
            "mat_class_2011_vs_2024": tr.reset_index()})


# ------------------------------------------------------------------------------ F2 — trees against open ground

def f2(s1):
    p = s1[s1.variant == "primary"]
    o1 = p[p.outcome == "O1_coh"]
    rows = [(st, tr) for st in SETS for tr in TRACKS]
    with fs.style():
        fig, axs = plt.subplots(1, 3, figsize=fs.size("full", 80), layout="constrained",
                                gridspec_kw={"width_ratios": [1.2, 0.95, 1.0]})
        ax = axs[0]
        for i, (st, tr) in enumerate(rows):
            for k, (side, (lab, col, mk, _)) in enumerate(SIDE.items()):
                r = o1[(o1.set == st) & (o1.track == tr) & (o1.contrast == side)].iloc[0]
                yy = i + (-0.16 if k == 0 else 0.16)
                ax.plot([r.lo95, r.hi95], [yy, yy], color=col, lw=1.1, ls=SIDE[side][3])
                ax.plot(r.wooded_minus_open, yy, mk, ms=4, mfc=col if k == 0 else "white", mec=col, mew=0.9)
        ax.axvline(0, color=fs.MUTED, lw=0.6)
        ax.set_yticks(range(len(rows)), [f"{SETS[s]} · {t[:3]}" for s, t in rows])
        ax.invert_yaxis()
        ax.set_xlabel("Δ coherence, wooded − open\n← trees lower coherence")
        ax.set_title("Coherence (both tracks, three pair sets)", loc="left")
        fs.panel_label(ax, "a")
        ax = axs[1]
        for i, (st, tr) in enumerate(rows):
            r = o1[(o1.set == st) & (o1.track == tr) & (o1.contrast == "interaction")].iloc[0]
            ax.plot([r.lo95, r.hi95], [i, i], color=fs.TRACK[tr], lw=1.1)
            ax.plot(r.wooded_minus_open, i, TMARK[tr], ms=4, color=fs.TRACK[tr])
        ax.axvline(0, color=fs.MUTED, lw=0.6)
        ax.set_yticks(range(len(rows)), [])
        ax.invert_yaxis()
        ax.set_xlabel("Δ on the mat − Δ off the mat\n← trees cost more on the mat · less →")
        ax.set_title("Mat against solid ground", loc="left")
        fs.panel_label(ax, "b")
        ax = axs[2]
        o3 = p[(p.outcome == "O3_r_laser") & (p.contrast == "on_mat")]
        lab = {"2020_2021": "2021", "2025": "2025", "2021+2025": "2021 + 2025"}
        yrows = [(y, t) for y in ("2021+2025", "2020_2021", "2025") for t in TRACKS]
        for i, (y, t) in enumerate(yrows):
            r = o3[(o3.set == y) & (o3.track == t)].iloc[0]
            ax.plot([r.lo95, r.hi95], [i, i], color=fs.TRACK[t], lw=1.3 if y == "2021+2025" else 0.9)
            ax.plot(r.wooded_minus_open, i, TMARK[t], ms=4.5 if y == "2021+2025" else 3.5, color=fs.TRACK[t],
                    mfc=fs.TRACK[t] if y == "2021+2025" else "white")
        ax.axvline(0, color=fs.MUTED, lw=0.6)
        ax.set_yticks(range(len(yrows)), [f"{lab[y]} · {t[:3]}" for y, t in yrows])
        ax.invert_yaxis()
        ax.set_xlabel("Δ r, wooded − open\n← trees weaken the laser read")
        ax.set_title("r with the P6 laser, mat", loc="left")
        fig.legend(handles=[Line2D([], [], color=c, ls=ls, marker=m, ms=4, mfc=c if s == "on_mat" else "white", label=l)
                            for s, (l, c, m, ls) in SIDE.items()] +
                   [Line2D([], [], color=fs.TRACK[t], marker=TMARK[t], ms=4, label=t) for t in TRACKS] +
                   [Line2D([], [], color=fs.MUTED, marker="o", ms=3.5, mfc="white", lw=0.9, label="single year (c)")],
                   loc="outside lower center", ncol=6, fontsize=6.5)
        fs.panel_label(ax, "c")
    export(fig, "F2_trees_vs_open", ["hub:landcover_x132/stage1.csv (variant primary)"],
           {"intervals": "95 % block bootstrap, 2000 resamples of 200 m blocks, seed 0 (X-132)"},
           {"coherence": o1, "laser_read": o3})


# ------------------------------------------------------------------------------ F3 — what drives coherence

def f3(s2, r2_note: str):
    t = s2.sort_values("share_of_r2", ascending=True).reset_index(drop=True)
    y = np.arange(len(t))
    with fs.style():
        fig, axs = plt.subplots(1, 2, figsize=fs.size("full", 92), layout="constrained", sharey=True,
                                gridspec_kw={"width_ratios": [1.2, 1]})
        ax = axs[0]
        for i, r in t.iterrows():
            col = fs.INK if r.matters else fs.MUTED
            ax.plot([r.lo95, r.hi95], [i, i], color=col, lw=1.0)
            mk = {"0/1": "s", "per SD": "o", "product": "D"}[r.scale]
            ax.plot(r.effect_fisher_z, i, mk, ms=4, mfc=col if r.matters else "white", mec=col, mew=0.9)
        ax.axvline(0, color=fs.MUTED, lw=0.6)
        ax.set_yticks(y, [FACTOR_LABEL.get(f, f) for f in t.factor])
        ax.set_xlabel("Effect on coherence (Fisher z = atanh γ)\n← lowers coherence · raises →")
        ax.set_title("Effects (95 % block-bootstrap intervals)", loc="left")
        fig.legend(handles=[Line2D([], [], ls="", marker="o", ms=4, color=fs.INK, label="per SD of the factor"),
                           Line2D([], [], ls="", marker="s", ms=4, color=fs.INK, label="0/1 factor"),
                           Line2D([], [], ls="", marker="D", ms=4, color=fs.INK, label="interaction"),
                           Line2D([], [], ls="", marker="o", ms=4, mfc="white", mec=fs.MUTED,
                                  label="below the pre-registered rule")],
                   loc="outside lower center", ncol=4, fontsize=6.3)
        fs.panel_label(ax, "a")
        ax = axs[1]
        ax.barh(y, t.r2_drop, color=[fs.INK if m else "#c9ced4" for m in t.matters], height=0.62,
                hatch=[None if m else "////" for m in t.matters], edgecolor=fs.MUTED, linewidth=0.3)
        ax.hlines(y, t.r2_drop_lo95, t.r2_drop_hi95, color=fs.MUTED, lw=0.8)       # percentile interval as drawn
        ax.set_xlim(left=0)
        ax.set_xlabel("R² lost when the factor alone is left out")
        ax.set_title("Share of explained variance", loc="left")
        ax.text(0.98, 0.02, r2_note, transform=ax.transAxes, ha="right", va="bottom", fontsize=6.5, color=fs.MUTED)
        fs.panel_label(ax, "b")
    export(fig, "F3_drivers_of_coherence", ["hub:landcover_x132/stage2.csv", "hub:landcover_x132/README.md (R²)"],
           {"model": "one least-squares fit on cell × consecutive pair, both tracks, 2020–2025; intervals 2000 block "
                     "resamples (X-132); filled / dark = 'matters' by the pre-registered rule"},
           {"stage2": t})


# ------------------------------------------------------------------------------ F4 — footprint

def f4(st, cells, fp, share, sec, p6):
    xlim, ylim = box_around(cells[cells.role == "mat"], pad=500)
    with fs.style():
        fig = plt.figure(figsize=fs.size("full", 120), layout="constrained")
        sf = fig.subfigures(2, 1, height_ratios=[1.45, 1])
        axm = sf[0].subplots(1, 2)
        axb = sf[1].subplots(1, 2)
        norm = TwoSlopeNorm(0, -1, 1)
        cmap = plt.get_cmap("RdBu_r").with_extremes(bad=BAD)
        for i, tr in enumerate(TRACKS):
            d = fp[(fp.track == tr) & (fp.years == "2021+2025")].set_index("cell").join(cells)
            thr = float(share[share.track == tr].threshold_r.iloc[0])
            ax = axm[i]
            im = ax.imshow(grid_of(st, d.row, d["col"], d.r), extent=extent(st), origin="upper", cmap=cmap, norm=norm,
                           interpolation="nearest", zorder=1)
            fol = (d.role == "mat") & (d.r > thr)
            cell_outline(ax, st.x.values, st.y.values, grid_of(st, d.row[fol], d["col"][fol], np.ones(fol.sum())) == 1,
                         colors=fs.INK, linewidths=0.9, zorder=6)
            mat_outline(ax, st, linewidths=1.0, linestyles=":")
            p6_star(ax, *p6)
            map_axes(ax, xlim, ylim, scalebar=(i == 0), north=(i == 1), ylabel=(i == 0))
            sr = share[share.track == tr].iloc[0]
            ax.set_title(f"{tr.capitalize()}: {sr.share:.0%} of mat cells follow P6 "
                         f"[{sr.lo95:.0%}, {sr.hi95:.0%}]", loc="left")
            fs.panel_label(ax, "ab"[i])
        sf[0].colorbar(im, ax=list(axm), shrink=0.8, label="r with the P6 laser\n(six-day increments)")
        sf[0].legend(handles=[Line2D([], [], color=fs.INK, lw=0.9, label="mat cells above the solid-ground 95th percentile"),
                              Line2D([], [], color=fs.CLASS["mat"], lw=1.0, ls=":", label="floating mat"),
                              Line2D([], [], ls="", marker="*", ms=7, mfc=fs.CLASS["laser"], mec="white", label="P6")],
                     loc="outside lower center", ncol=3, fontsize=6.5)
        ax = axb[0]
        bins = np.linspace(-0.6, 1, 33)
        for tr in TRACKS:
            d = fp[(fp.track == tr) & (fp.years == "2021+2025")].set_index("cell").join(cells)
            for role, ls in (("mat", "-"), ("null", "--")):
                v = d[d.role == role].r.dropna()
                ax.hist(v, bins=bins, density=True, histtype="step", color=fs.TRACK[tr], ls=ls, lw=TLW[tr])
            ax.axvline(float(share[share.track == tr].threshold_r.iloc[0]), color=fs.TRACK[tr], lw=0.6, ls=":")
        ax.set_xlabel("r with the P6 laser")
        ax.set_ylabel("Density of cells")
        ax.set_title("Mat cells against the solid-ground null", loc="left")
        ax.legend(handles=[Line2D([], [], color=fs.TRACK[t], lw=TLW[t], label=t) for t in TRACKS] +
                  [Line2D([], [], color=fs.INK, lw=1, ls="-", label="mat"),
                   Line2D([], [], color=fs.INK, lw=1, ls="--", label="solid ground (null)"),
                   Line2D([], [], color=fs.INK, lw=0.6, ls=":", label="null 95th percentile")],
                  loc="upper center", bbox_to_anchor=(0.5, -0.3), ncol=3, fontsize=6.3)
        ax.set_xlim(-0.6, 1)
        fs.panel_label(ax, "c")
        ax = axb[1]
        dist = sec[sec.subset.str.startswith("dist")].copy()
        labels = list(dict.fromkeys(dist.subset))
        for k, tr in enumerate(TRACKS):
            d = dist[dist.track == tr].set_index("subset").loc[labels]
            xx = np.arange(len(labels)) + (k - 0.5) * 0.18
            ax.vlines(xx, d.lo95, d.hi95, color=fs.TRACK[tr], lw=0.9)
            ax.plot(xx, d.share, TMARK[tr], ms=4, color=fs.TRACK[tr], label=tr)
        ax.axhline(0.05, color=fs.MUTED, lw=0.6, ls=":")
        ax.text(len(labels) - 0.6, 0.065, "rate expected by chance (5 %)", ha="right", va="bottom", fontsize=6.3,
                color=fs.MUTED)
        ncell = dist[dist.track == TRACKS[0]].set_index("subset").cells
        ax.set_xticks(range(len(labels)), [f"{s.replace('dist ', '')}\n(n = {int(ncell[s])})" for s in labels])
        ax.set_xlim(-0.5, len(labels) - 0.5)
        ax.set_ylim(0, 1.02)
        ax.set_xlabel("Distance from P6")
        ax.set_ylabel("Share of mat cells following P6")
        ax.set_title("How far the motion reaches", loc="left")
        ax.legend(loc="lower left", bbox_to_anchor=(0.0, 0.12), fontsize=6.5)
        fs.panel_label(ax, "d")
    export(fig, "F4_footprint", ["hub:footprint_x133/footprint.csv.gz", "hub:footprint_x133/share.csv",
                                 "hub:footprint_x133/secondary.csv", "hub:footprint_x133/cells.csv"],
           {"pairs": "X-126's consecutive six-day pairs, 2021 (from 27 May) + 2025; hard-target reference",
            "intervals": "95 % block bootstrap over the mat's 200 m blocks, 2000 resamples, seed 0"},
           {"share": share, "by_distance": dist})


# ------------------------------------------------------------------------------ F5 — gain against depth

def f5(cells, fp, share):
    edges = np.arange(0, 400, 40)
    rows = []
    with fs.style():
        fig, ax = plt.subplots(figsize=fs.size("one_half", 75), layout="constrained")
        for tr in TRACKS:
            d = fp[(fp.track == tr) & (fp.years == "2021+2025")].set_index("cell").join(cells)
            thr = float(share[share.track == tr].threshold_r.iloc[0])
            m = d[(d.role == "mat") & (d.r > thr)]
            ax.scatter(m.depth_in_mat_m, m.gain, s=4, marker=TMARK[tr], color=fs.TRACK[tr], alpha=0.25, linewidths=0,
                       rasterized=True)
            for a, b in zip(edges[:-1], edges[1:]):
                g = m[(m.depth_in_mat_m >= a) & (m.depth_in_mat_m < b)].gain.dropna()
                if len(g) >= 5:
                    rows.append({"track": tr, "depth_lo_m": a, "depth_hi_m": b, "cells": len(g),
                                 "median": g.median(), "q25": g.quantile(0.25), "q75": g.quantile(0.75)})
            bt = pd.DataFrame([r for r in rows if r["track"] == tr])
            mid = (bt.depth_lo_m + bt.depth_hi_m) / 2
            ax.fill_between(mid, bt.q25, bt.q75, color=fs.TRACK[tr], alpha=0.18, linewidth=0)
            sr = share[share.track == tr].iloc[0]
            ax.plot(mid, bt["median"], marker=TMARK[tr], ms=4, color=fs.TRACK[tr], lw=1.2,
                    label=f"{tr}: Spearman ρ {sr.gain_depth_rho:.2f} [{sr.rho_lo95:.2f}, {sr.rho_hi95:.2f}]")
        ax.axhline(1, color=fs.MUTED, lw=0.6, ls=":")
        ax.text(0.99, 1.0, "gain 1: the P6 laser's amplitude", transform=ax.get_yaxis_transform(), ha="right",
                va="bottom", fontsize=6.3, color=fs.MUTED, bbox={"fc": "white", "ec": "none", "pad": 0.5})
        ax.set_xlabel("Depth into the mat (m from the edge)")
        ax.set_ylabel("Gain against the P6 laser (slope through origin)")
        ax.set_title("The motion grows towards the mat's interior (cells following P6)", loc="left")
        ax.legend(loc="lower right", fontsize=6.5)
    export(fig, "F5_gain_vs_depth", ["hub:footprint_x133/footprint.csv.gz", "hub:footprint_x133/share.csv"],
           {"summary": "points = cells; lines = median per 40 m of depth, band = interquartile range (descriptive); "
                       "ρ and its interval from X-133 (block bootstrap)",
            "caveat": "cells within 40–80 m of the edge mix mat and solid ground, which lowers their gain"},
           {"binned": pd.DataFrame(rows)})


# ------------------------------------------------------------------------------ F6 — season maps

def f6(st, cells, mp, agree, runs, p6):
    mat = cells[cells.role == "mat"]
    xlim, ylim = box_around(mat, pad=60)
    v = mp[mp.following]
    lim = np.ceil(float(np.nanpercentile(-v.vertical_change_mm, 98)) / 25) * 25           # sinking
    up = max(25.0, np.ceil(float(np.nanpercentile(v.vertical_change_mm, 99.5)) / 25) * 25)  # rising, if any
    norm = TwoSlopeNorm(0, -lim, up)
    cmap = plt.get_cmap("RdBu_r").with_extremes(bad=BAD)
    with fs.style():
        fig = plt.figure(figsize=fs.size("full", 146), layout="constrained")
        fig.suptitle("May–October sinking of the mat, from chains of six-day pairs. Only 2021 and 2025: the test was "
                     "pre-registered on six-day pairs,\nand 2022–24 had twelve-day pairs only (one satellite).",
                     x=0.01, ha="left", fontsize=7, color=fs.MUTED)
        gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 0.95])
        for i, yr in enumerate(("2021", "2025")):
            for j, tr in enumerate(TRACKS):
                ax = fig.add_subplot(gs[i, j])
                d = mp[(mp.year.astype(str) == yr) & (mp.track == tr)].set_index("cell").join(mat, how="inner")
                d = d[d.following]
                im = ax.imshow(grid_of(st, d.row, d["col"], d.vertical_change_mm), extent=extent(st), origin="upper",
                               cmap=cmap, norm=norm, interpolation="nearest", zorder=1)
                mat_outline(ax, st, linewidths=1.0)
                p6_star(ax, *p6)
                map_axes(ax, xlim, ylim, scalebar=(i == 1 and j == 0), north=(i == 0 and j == 1), ylabel=(j == 0))
                rr = runs[(runs.year.astype(str) == yr) & (runs.track == tr)].iloc[0]
                ax.set_title(f"{yr}, {tr}\n{pd.Timestamp(rr.start):%d %b}–{pd.Timestamp(rr.end):%d %b}, "
                             f"{int(rr.pairs)} pairs, 1σ ≈ {float(d.sigma_mm.median()):.0f} mm", loc="left",
                             fontsize=7)
                fs.panel_label(ax, "abde"[2 * i + j])
        cb = fig.colorbar(im, ax=fig.axes, shrink=0.6, location="bottom", aspect=40,
                          label="Vertical change over the season's run (mm; cells following P6 on that track; grey: not following)")
        cb.ax.set_xscale("linear")
        ax = fig.add_subplot(gs[0, 2])
        for yr, mk, col in (("2021", "o", fs.INK), ("2025", "^", fs.MUTED)):
            w = mp[mp.year.astype(str) == yr].pivot_table(index="cell", columns="track", values="vertical_change_mm")
            both = mp[mp.year.astype(str) == yr].groupby("cell").following.all()
            w = w[both.reindex(w.index).fillna(False).astype(bool)].dropna()
            a = agree[agree.year.astype(str) == yr].iloc[0]
            ax.scatter(w.ascending, w.descending, s=5, marker=mk, color=col, alpha=0.5, linewidths=0, rasterized=True,
                       label=f"{yr}: r {a.r_asc_desc:.2f} [{a.lo95:.2f}, {a.hi95:.2f}], n={int(a.cells)}")
        ax.plot([-lim * 1.6, lim * 0.5], [-lim * 1.6, lim * 0.5], color=fs.MUTED, lw=0.6, ls=":")
        ax.set_xlabel("Ascending (mm)")
        ax.set_ylabel("Descending (mm)")
        ax.set_aspect("equal", adjustable="datalim")
        ax.set_title("The two tracks, cell by cell", loc="left")
        ax.legend(loc="lower right", fontsize=6.2, markerscale=1.6)
        fs.panel_label(ax, "c")
        ax = fig.add_subplot(gs[1, 2])
        xt, labels = [], []
        for k, (yr, tr) in enumerate([(y, t) for y in ("2021", "2025") for t in TRACKS]):
            rr = runs[(runs.year.astype(str) == yr) & (runs.track == tr)].iloc[0]
            sig = float(rr.p6_sigma_mm) * np.sqrt(rr.pairs_with_laser / rr.pairs)
            ax.errorbar(k - 0.12, rr.p6_vertical_on_laser_pairs_mm, yerr=sig, fmt=TMARK[tr], ms=4, color=fs.TRACK[tr],
                        capsize=2, elinewidth=0.9)
            ax.plot(k + 0.12, rr.laser_vertical_on_laser_pairs_mm, "*", ms=7, color=fs.INK)
            xt.append(k)
            labels.append(f"{yr}\n{tr[:3]}")
        ax.axhline(0, color=fs.MUTED, lw=0.6)
        ax.set_xticks(xt, labels)
        ax.set_ylabel("Season change at P6 (mm, vertical)")
        ax.set_title("P6: radar (±1σ) against the laser", loc="left")
        ax.legend(handles=[Line2D([], [], ls="", marker="o", color=fs.TRACK["ascending"], label="radar, ascending"),
                           Line2D([], [], ls="", marker="s", color=fs.TRACK["descending"], label="radar, descending"),
                           Line2D([], [], ls="", marker="*", ms=7, color=fs.INK, label="laser")],
                  loc="lower right", fontsize=6.3)
        fs.panel_label(ax, "f")
    runs = runs.assign(p6_sigma_on_laser_pairs_mm=runs.p6_sigma_mm * np.sqrt(runs.pairs_with_laser / runs.pairs))
    export(fig, "F6_season_maps", ["hub:footprint_x133/season_map.csv.gz", "hub:footprint_x133/tracks_agree.csv",
                                   "hub:footprint_x133/season_runs.csv"],
           {"colour_scale": f"diverging, white at 0; −{lim:.0f} mm (98th pct of sinking) to +{up:.0f} mm",
            "p6_panel": "radar and laser summed over the run's pairs that have a laser change; σ = the run's σ × "
                        "sqrt(pairs with laser / pairs)"},
           {"runs": runs, "agree": agree})


# ------------------------------------------------------------------------------ README (captions, alt text)

def readme(c132, s1, s2, r2, share, agree, runs) -> None:
    mat = c132[c132.on_mat == 1]
    newly = int(((mat.cls == "wooded") & (mat.cls_2011 != "wooded")).sum())
    o1 = s1[(s1.variant == "primary") & (s1.outcome == "O1_coh")]
    rng = lambda side: (o1[o1.contrast == side].wooded_minus_open.min(), o1[o1.contrast == side].wooded_minus_open.max())  # noqa: E731
    on, off = rng("on_mat"), rng("off_mat")
    sa, sd = (share[share.track == t].iloc[0] for t in TRACKS)
    top = s2.sort_values("share_of_r2", ascending=False)
    lead = ", ".join(FACTOR_LABEL[f].lower() for f in top[top.matters].factor[:4])
    a21, a25 = (agree[agree.year.astype(str) == y].iloc[0] for y in ("2021", "2025"))
    rr = runs.set_index(["year", "track"])
    p = lambda y, t, c: float(rr.loc[(int(y), t), c])  # noqa: E731
    figs = [
        ("F1_dated_land_cover", "Dated land cover around P6",
         f"(a) Colour-infrared orthophoto of 19 July 2022 (GUGiK; near-infrared, red and green shown as red, green and "
         f"blue, as delivered, resampled to 2 m) with the 40 m cells of X-132 within 1 km of P6: wooded (90th percentile "
         f"of 1 m canopy height ≥ 2 m in the 19 September 2024 laser scan; filled) and open (< 0.5 m; hollow). Green "
         f"line: the floating mat; star: P6, the laser point. (b) Canopy height per cell, 2024 scan (its own surface "
         f"minus terrain model). (c) Change since the 24 April 2011 scan; the 2011 scan is leaf-off, so part of the "
         f"increase in deciduous stands may be phenology rather than growth. On the mat, {newly} of the "
         f"{int((mat.cls == 'wooded').sum())} cells wooded in 2024 were not wooded in 2011. Grey: outside the analysed cells.",
         "Three maps of the same 2 km square around the fen. Left, a false-colour aerial photo with small black squares "
         "on forest and white squares on open ground; the open squares cover the floating mat and the grassland round "
         "it. Middle, canopy height: tall forest to the west and south, low vegetation on the mat. Right, height change "
         "since 2011: mostly small gains, a few felled stands in the south."),
        ("F2_trees_vs_open", "Trees against open ground, on and off the mat",
         f"Median over matched pairs of wooded − open cells (matched on depth into the mat, or distance to its edge, and "
         f"on distance to P6), with 95 % block-bootstrap intervals (2000 resamples of 200 m blocks). (a) Coherence of "
         f"consecutive pairs per set and track: on the mat {on[0]:.2f} to {on[1]:.2f}, off it {off[0]:.2f} to "
         f"{off[1]:.2f}. (b) The difference of the two effects (on − off the mat); no interval excludes zero on both "
         f"tracks of any set. (c) The tree effect on how a mat cell's six-day phase follows the P6 laser (r), 2021 and "
         f"2025 pooled (filled) and apart (hollow).",
         "Three forest plots. In the first, every wooded-minus-open coherence difference is negative, on the mat and "
         "off it, with off-mat effects slightly larger. In the second, the differences between the two effects sit "
         "around zero to slightly positive with intervals touching zero. In the third, the effects on the laser "
         "correlation scatter around zero."),
        ("F3_drivers_of_coherence", "What drives coherence",
         f"One least-squares regression of Fisher-z coherence on every factor over all cell × consecutive pair rows "
         f"(both tracks, 2020–2025; {r2.split('. Effects')[0]}). (a) Effects per standard deviation of the factor "
         f"(circles), between the two states of a 0/1 factor (squares) or of an interaction (diamonds), with 95 % "
         f"block-bootstrap intervals; hollow: below the pre-registered rule (interval excluding zero and at least 1 % of "
         f"R²). (b) R² lost when the factor alone is left out, with its interval. Leading factors: {lead}.",
         "A dot-and-whisker chart of fifteen factors next to a bar chart of the variance each explains. Season and "
         "track have by far the longest bars; land-cover factors such as canopy height and canopy cover have almost "
         "none, and depth into the mat none."),
        ("F4_footprint", "How far P6's six-day motion reaches",
         f"(a, b) Pearson r of every cell's consecutive six-day phase increments (hard-target reference) with the P6 "
         f"laser's line-of-sight increments, the pairs of 2021 (from 27 May) and 2025. Black outline: mat cells whose r "
         f"exceeds the 95th percentile of solid-ground cells read the same way (r > {sa.threshold_r:.2f} ascending, "
         f"> {sd.threshold_r:.2f} descending). (c) Distributions of r on the mat and on solid ground. (d) The share of "
         f"mat cells that follow, by distance from P6, with block-bootstrap intervals. Overall {sa.share:.0%} "
         f"[{sa.lo95:.0%}, {sa.hi95:.0%}] ascending and {sd.share:.0%} [{sd.lo95:.0%}, {sd.hi95:.0%}] descending. An "
         f"off-mat patch to the south-east also follows P6; it lies in the null and only raises the threshold.",
         "Two maps where the floating mat stands out in red against pale surroundings, meaning its cells move in step "
         "with the laser at P6 while solid ground does not. Below, the mat's correlations sit well to the right of the "
         "solid-ground distribution, and the share of mat cells following P6 stays high out to 1 km on the ascending "
         "track and falls to under half beyond 500 m on the descending track."),
        ("F5_gain_vs_depth", "The motion grows towards the mat's interior",
         f"Gain (slope through the origin of a cell's six-day increments against the P6 laser's) of the mat cells that "
         f"follow P6, against depth into the mat. Points: cells; lines: median per 40 m of depth; bands: interquartile "
         f"range (descriptive). Spearman ρ with its block-bootstrap interval in the legend. Cells within 40–80 m of the "
         f"edge mix mat with solid ground, which lowers their gain.",
         "A scatter plot with two rising curves: from about half the laser's amplitude near the mat edge to about the "
         "laser's amplitude, or slightly above on the ascending track, beyond 200 m into the mat."),
        ("F6_season_maps", "A first season map of the mat",
         f"(a, b, d, e) Vertical change from May to October summed over the longest unbroken run of consecutive six-day "
         f"pairs (dates and number of pairs in each title), for the cells that follow P6 on that track; 1σ = P6's "
         f"increment error against the laser × √pairs. (c) The two tracks cell by cell (cells following on both): r "
         f"{a21.r_asc_desc:.2f} [{a21.lo95:.2f}, {a21.hi95:.2f}] in 2021, {a25.r_asc_desc:.2f} [{a25.lo95:.2f}, "
         f"{a25.hi95:.2f}] in 2025. (f) At P6, the radar's season change against the laser's over the run's pairs that "
         f"have a laser change: 2021 ascending {p(2021, 'ascending', 'p6_vertical_on_laser_pairs_mm'):.0f} mm against "
         f"{p(2021, 'ascending', 'laser_vertical_on_laser_pairs_mm'):.0f} mm; 2025 ascending "
         f"{p(2025, 'ascending', 'p6_vertical_on_laser_pairs_mm'):.0f} mm against "
         f"{p(2025, 'ascending', 'laser_vertical_on_laser_pairs_mm'):.0f} mm, so the 2025 map overstates the sinking.",
         "Four maps of the mat, all blue, meaning it sank over each summer, deepest along the north and around the "
         "inner lake. A scatter shows the two tracks agreeing loosely. A small chart at P6 shows radar and laser close "
         "in 2021 and the radar sinking about twice as far as the laser in 2025."),
    ]
    lines = ["# Figures — dated land cover, trees against open ground, and the footprint of P6's motion (generated)", "",
             "Generated by SAr-adapted `experiments/figures_landcover_footprint_c073.py` "
             "(`make run E=figures_landcover_footprint_c073`, ticket C-073); do not edit. Drawn from the committed "
             "outputs of C-072 (dated land cover), X-132 (`../insar_full_test/landcover_x132/`) and X-133 "
             "(`../insar_full_test/footprint_x133/`); nothing is re-estimated. Both analyses were pre-registered.", "",
             "Each figure: `<name>.pdf` (vector, Arial embedded as TrueType, rasters at 600 dpi), `<name>.png` "
             "(preview, 300 dpi), `<name>.json` (provenance: sources, code commit, size), and `data/<name>_*.csv` "
             "(what is plotted). Widths are Elsevier's (190 mm full, 140 mm one and a half); a journal is not chosen, so "
             "its own rules are still to check. Colour is never the only cue: tracks also differ by marker or line "
             "weight, classes by fill or line style. The two track colours pass 3:1 against white but are close in "
             "grey, which is why they carry a second cue.", ""]
    for name, title, cap, alt in figs:
        lines += [f"## {name} — {title}", "", f"![{title}]({name}.png)", "", f"**Caption.** {cap}", "",
                  f"**Alt text.** {alt}", ""]
    lines += ["## Limits shared by every figure", "",
              "- One laser point (P6): the footprint measures whether a cell follows P6's timing, not its own motion.",
              "- Design years only (2020–2025; the laser from 27 May 2021); the held-out years are untouched.",
              "- The season map's uncertainty assumes P6's per-pair error elsewhere and adds no unwrapping error; its "
              "test against the UAV is fixed in X-133's pre-registration and waits for the rasters (C-065).", ""]
    (OUT / "README.md").write_text("\n".join(lines), encoding="utf-8")


# ------------------------------------------------------------------------------ main

def main() -> None:
    (OUT / "data").mkdir(parents=True, exist_ok=True)
    st = static()
    pu = pd.read_csv(CUBE / "silver" / "plot_units.csv").set_index("plot")
    r6, c6 = int(pu.loc["P6", "row"]), int(pu.loc["P6", "col"])
    p6 = (float(st.x.values[c6]), float(st.y.values[r6]))
    c132 = pd.read_csv(X132 / "cells.csv", index_col=0)
    c133 = pd.read_csv(X133 / "cells.csv", index_col=0, keep_default_na=False, na_values=[""])   # role "null" is a label
    s1, s2 = pd.read_csv(X132 / "stage1.csv"), pd.read_csv(X132 / "stage2.csv")
    r2 = next(l for l in (X132 / "README.md").read_text().splitlines() if l.startswith("R² "))
    fp = pd.read_csv(X133 / "footprint.csv.gz")
    share, sec = pd.read_csv(X133 / "share.csv"), pd.read_csv(X133 / "secondary.csv")
    mp, agree, runs = (pd.read_csv(X133 / f) for f in ("season_map.csv.gz", "tracks_agree.csv", "season_runs.csv"))
    f1(st, c132, p6)
    f2(s1)
    f3(s2, "full model " + r2.split(". Effects")[0])
    f4(st, c133, fp, share, sec, p6)
    f5(c133, fp, share)
    f6(st, c133, mp, agree, runs, p6)
    readme(c132, s1, s2, r2, share, agree, runs)
    print("figures written to", OUT)


if __name__ == "__main__":
    main()

"""F1, final version — the site: where it is, the floating mat on a colour-infrared orthophoto, its wells, the laser at
P6, the references and the radar's two looks.

Round 2 (2026-10-31, after the figure review in the hub NOTES): (a) where Rzecin lies in Poland and Europe (Natural
Earth outlines); (b) the peatland on the GUGiK colour-infrared orthophoto of 19 July 2022 (0.15 m, vegetation in red),
with the mat's 40 m outline, the residual lake, the hard targets, the matched grassland, the wells and the laser; (c) the
transect of wells on the same orthophoto at full detail, with the 40 m cells around it; (d) the two tracks' flight and
look directions with their incidence at the mat, computed from the cube's look vectors. Round 1 drew (b) and (c) on a
grey LiDAR hillshade. Every value drawn is read from the cube or the field table; no number is typed. Reads the cube and
the hub; writes only to the hub: ``08_deliverables/article_plan/figures/final/``.

    make run E=article_fig_f1
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import xarray as xr  # noqa: E402
import geopandas as gpd  # noqa: E402
import rasterio  # noqa: E402
from rasterio.enums import Resampling  # noqa: E402
from rasterio.transform import from_origin  # noqa: E402
from rasterio.warp import reproject, transform_bounds  # noqa: E402
from rasterio.windows import from_bounds  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import FancyArrowPatch, Patch, Rectangle  # noqa: E402
from matplotlib_scalebar.scalebar import ScaleBar  # noqa: E402

from insar_wetlands import figstyle as fs  # noqa: E402
from insar_wetlands.paths import cube_root, hub_root  # noqa: E402

HUB, CUBE = hub_root(), cube_root()
OUT = HUB / "08_deliverables" / "article_plan" / "figures" / "final"
MAT_WELLS = ("P6", "P7", "P8", "P9")          # field_plots_x061: P6–P9 on the mat, P1–P5 the water-table reference
CELL = 40.0
GUGIK = CUBE / "raw" / "gugik_dated"
ORTHO_YEAR = 2022                                  # 19 July 2022, 0.15 m: the summer flight (sheets.csv)
NE = CUBE / "raw" / "naturalearth" / "ne_50m_admin_0_countries.zip"
GRID_CRS = "EPSG:32633"
SITE_LAT, SITE_LON = 52 + 45 / 60, 16 + 18 / 60    # §2.1: 52°45′ N, 16°18′ E
LAKE = "#7cc7ff"


BOARDWALK = "#cf222e"                              # as in F5


def ortho_sheets() -> tuple[list, str]:
    """The colour-infrared sheets of the chosen flight and its date, from the GUGiK index the cube keeps."""
    sh = pd.read_csv(GUGIK / "sheets.csv")
    sh = sh[(sh.kind == "cir") & (sh.year == ORTHO_YEAR)]
    return [GUGIK / f for f in sh.file], str(sh.date.iloc[0])


def ortho(bounds: tuple[float, float, float, float], res: float) -> np.ndarray:
    """The orthophoto on the radar grid's projection over ``bounds`` (x0, x1, y0, y1) at ``res`` metres: each sheet read
    decimated toward ``res`` (averaging), then reprojected; bands as delivered (NIR, red, green: vegetation red),
    stretched 1–99 % per band for display."""
    x0, x1, y0, y1 = bounds
    w, h = int(round((x1 - x0) / res)), int(round((y1 - y0) / res))
    dst = np.zeros((3, h, w), dtype=np.float32)
    dt = from_origin(x0, y1, res, res)
    for f in ortho_sheets()[0]:
        with rasterio.open(f) as src:
            l, b_, r, t = transform_bounds(GRID_CRS, src.crs, x0, y0, x1, y1, densify_pts=21)
            try:
                win = from_bounds(l - 20, b_ - 20, r + 20, t + 20, src.transform).intersection(
                    rasterio.windows.Window(0, 0, src.width, src.height))
            except rasterio.errors.WindowError:                          # the sheet lies outside the bounds
                continue
            k = max(1, int(res / src.res[0] / 1.5))
            out = (3, max(1, int(win.height // k)), max(1, int(win.width // k)))
            arr = src.read(window=win, out_shape=out, resampling=Resampling.average).astype(np.float32)
            st = src.window_transform(win) * rasterio.Affine.scale(win.width / out[2], win.height / out[1])
            tmp = np.zeros_like(dst)
            reproject(arr, tmp, src_transform=st, src_crs=src.crs, dst_transform=dt, dst_crs=GRID_CRS,
                      resampling=Resampling.bilinear, src_nodata=0, dst_nodata=0)
            dst = np.where(tmp > 0, tmp, dst)
    rgb = np.zeros((h, w, 3), dtype=np.float32)
    for i in range(3):
        v = dst[i][dst[i] > 0]
        lo, hi = np.percentile(v, [1, 99]) if v.size else (0, 255)
        rgb[..., i] = np.clip((dst[i] - lo) / (hi - lo), 0, 1)
    rgb[(dst == 0).all(axis=0)] = 1.0
    return rgb


def cell_outline(ax, x, y, mask, **kw) -> None:
    """The outline of a set of 40 m cells, drawn along cell edges (the grid the analysis uses, not a smoothed shape)."""
    m = np.pad(mask.astype(bool), 1)
    xe = np.concatenate([x - CELL / 2, [x[-1] + CELL / 2]])
    ye = np.concatenate([y - np.sign(y[1] - y[0]) * CELL / 2, [y[-1] + np.sign(y[1] - y[0]) * CELL / 2]])
    segs = []
    for r in range(m.shape[0] - 1):
        for c in range(m.shape[1] - 1):
            if m[r + 1, c] != m[r + 1, c + 1] and r < len(y):        # vertical edge between (r, c-1) and (r, c)
                segs.append(((xe[c], ye[r]), (xe[c], ye[r + 1])))
            if m[r, c + 1] != m[r + 1, c + 1] and c < len(x):        # horizontal edge between (r-1, c) and (r, c)
                segs.append(((xe[c], ye[r]), (xe[c + 1], ye[r])))
    from matplotlib.collections import LineCollection
    ax.add_collection(LineCollection(segs, **kw))


def look_geometry(st: xr.Dataset) -> dict:
    """Per track, over the mat: look azimuth (ground direction the radar looks toward), flight heading (right-looking),
    and incidence — medians of the cube's look vectors (ground → satellite)."""
    m = st.zone.values == 1
    out = {}
    for tr, k in (("ascending", "asc"), ("descending", "des")):
        e, n = (float(np.nanmedian(st[f"los_{c}_{k}"].values[m])) for c in "en")
        look = np.degrees(np.arctan2(-e, -n)) % 360
        out[tr] = {"look_az": look, "heading": (look - 90) % 360,
                   "incidence": float(np.nanmedian(st[f"incidence_{k}_deg"].values[m]))}
    return out


def panel_locator(ax) -> None:
    """Where the site lies: Poland among its neighbours (Natural Earth 1:50m, Lambert azimuthal equal-area for
    Europe), the site from the coordinates of §2.1."""
    w = gpd.read_file(f"zip://{NE}").to_crs("EPSG:3035")
    site = gpd.GeoSeries(gpd.points_from_xy([SITE_LON], [SITE_LAT]), crs="EPSG:4326").to_crs("EPSG:3035").iloc[0]
    pl = w[w.ADMIN == "Poland"]
    w.plot(ax=ax, fc="#eef1f4", ec=fs.MUTED, lw=0.3, zorder=0)
    pl.plot(ax=ax, fc="#d9dde2", ec=fs.INK, lw=0.6, zorder=1)
    ax.plot(site.x, site.y, "*", ms=9, mfc=fs.CLASS["laser"], mec="white", mew=0.5, zorder=3)
    c = pl.geometry.iloc[0].representative_point()
    ax.text(c.x + 60e3, c.y - 60e3, "Poland", fontsize=7, ha="center", color=fs.INK, zorder=2)
    ax.text(site.x, site.y + 55e3, "Rzecin", fontsize=7, ha="center", va="bottom", fontweight="bold", zorder=3,
            bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.8))
    b = pl.total_bounds
    pad = 0.55 * (b[2] - b[0])
    ax.set_xlim(b[0] - pad, b[2] + pad), ax.set_ylim(b[1] - 0.45 * pad, b[3] + 0.45 * pad)
    ax.set_aspect("equal")
    ax.set_xticks([]), ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(True), sp.set_linewidth(0.6)


def panel_site(ax, st, li, px, box) -> None:
    x, y = st.x.values, st.y.values
    rm, cm = np.where(st.zone.values == 1)
    pad = 14
    xlim = (x[max(cm.min() - pad, 0)] - CELL / 2, x[min(cm.max() + pad, len(x) - 1)] + CELL / 2)
    ylim = (y[min(rm.max() + pad, len(y) - 1)] - CELL / 2, y[max(rm.min() - pad, 0)] + CELL / 2)
    ax.imshow(ortho((xlim[0], xlim[1], ylim[0], ylim[1]), 2.0), extent=[xlim[0], xlim[1], ylim[0], ylim[1]],
              origin="upper", interpolation="lanczos", zorder=0)
    grass = st.zone.values == 3
    cell_outline(ax, x, y, grass, colors=fs.CLASS["grassland"], linewidths=0.7, zorder=2)
    cell_outline(ax, x, y, st.zone.values == 2, colors=LAKE, linewidths=1.2, zorder=3)
    cell_outline(ax, x, y, st.zone.values == 1, colors="white", linewidths=2.6, zorder=3)
    cell_outline(ax, x, y, st.zone.values == 1, colors=fs.CLASS["mat"], linewidths=1.4, zorder=3)
    hr, hc = np.where(st.hard_target.values.astype(bool))
    for r, c in zip(hr, hc):
        ax.add_patch(Rectangle((x[c] - CELL / 2, y[r] - CELL / 2), CELL, CELL, fc="none", ec="white", lw=2.0,
                               zorder=4))
        ax.add_patch(Rectangle((x[c] - CELL / 2, y[r] - CELL / 2), CELL, CELL, fc="none", ec=fs.CLASS["hard targets"],
                               lw=1.0, zorder=4))
    ref = px[~px.index.isin(MAT_WELLS)]
    mat = px[px.index.isin(MAT_WELLS)]
    ax.plot(ref.E, ref.N, "o", ms=2.6, mfc="white", mec=fs.INK, mew=0.6, zorder=5)
    ax.plot(mat.E, mat.N, "o", ms=2.6, mfc=fs.INK, mec="white", mew=0.4, zorder=5)
    ax.plot(px.loc["P6", "E"], px.loc["P6", "N"], "*", ms=8, mfc=fs.CLASS["laser"], mec="white", mew=0.5, zorder=6)
    ax.add_patch(Rectangle((box[0], box[2]), box[1] - box[0], box[3] - box[2], fc="none", ec="white", lw=1.0,
                           ls="--", zorder=6))
    ax.text(box[1] + 8, box[3], "(c)", fontsize=7, va="top", color="white", fontweight="bold", zorder=6)
    ax.set_xlim(*xlim), ax.set_ylim(*ylim)
    ax.set_aspect("equal")
    ticks = lambda lo, hi: np.arange(np.ceil(lo / 500) * 500, hi, 500)
    ax.set_xticks(ticks(*xlim)), ax.set_yticks(ticks(*ylim))
    ax.xaxis.set_major_formatter(lambda v, _: f"{v / 1000:.1f}")
    ax.yaxis.set_major_formatter(lambda v, _: f"{v / 1000:.1f}")
    ax.set_xlabel("Easting (km, UTM 33N)")
    ax.set_ylabel("Northing (km, UTM 33N)")
    for s in ("top", "right"):
        ax.spines[s].set_visible(True)
    ax.add_artist(ScaleBar(1, "m", location="lower left", length_fraction=0.18, box_alpha=0.85, font_properties={
        "size": 7}, border_pad=0.6))
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    nx, ny = x1 - 0.05 * (x1 - x0), y1 - 0.16 * (y1 - y0)
    ax.add_patch(FancyArrowPatch((nx, ny), (nx, ny + 0.09 * (y1 - y0)), arrowstyle="-|>", mutation_scale=9,
                                 color="white", lw=0.9, zorder=7))
    ax.text(nx, ny + 0.10 * (y1 - y0), "N", ha="center", va="bottom", fontsize=7, color="white", zorder=7)


def panel_transect(ax, st, li, px, box) -> None:
    ax.imshow(ortho(box, 0.3), extent=list(box), origin="upper", interpolation="lanczos", zorder=0)
    x, y = st.x.values, st.y.values
    for xv in np.concatenate([x - CELL / 2, [x[-1] + CELL / 2]]):
        if box[0] <= xv <= box[1]:
            ax.axvline(xv, color="white", lw=0.4, alpha=0.7, zorder=1)
    for yv in np.concatenate([y + CELL / 2, [y[-1] - CELL / 2]]):
        if box[2] <= yv <= box[3]:
            ax.axhline(yv, color="white", lw=0.4, alpha=0.7, zorder=1)
    cell_outline(ax, x, y, st.zone.values == 1, colors=fs.CLASS["mat"], linewidths=1.6, zorder=2)
    bw = pd.read_csv(HUB / "08_deliverables" / "field_boardwalk_x060" / "boardwalk_centreline.csv").sort_values("N")
    bw = bw[bw.supported]                          # X-060's traced centreline (orthophoto), supported bins only
    ax.plot(bw.E, bw.N, color=BOARDWALK, lw=1.4, solid_capstyle="round", zorder=3)
    left = {"P9", "P7", "P5"}                      # wells west of their neighbours: label on the west side
    for p, row in px.iterrows():
        on = p in MAT_WELLS
        if p == "P6":
            ax.plot(row.E, row.N, "*", ms=10, mfc=fs.CLASS["laser"], mec="white", mew=0.5, zorder=5)
        else:
            ax.plot(row.E, row.N, "o", ms=4, mfc=fs.INK if on else "white", mec="white" if on else fs.INK,
                    mew=0.8, zorder=4)
        dx = -9 if p in left else 9
        ax.text(row.E + dx, row.N, p + (" laser" if p == "P6" else ""), fontsize=7, va="center",
                ha="right" if p in left else "left", zorder=6,
                bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.8))
    ax.set_xlim(box[0], box[1]), ax.set_ylim(box[2], box[3])
    ax.set_aspect("equal")
    ax.set_xticks([]), ax.set_yticks([])
    for s in ("top", "right"):
        ax.spines[s].set_visible(True)
    ax.add_artist(ScaleBar(1, "m", location="lower right", length_fraction=0.5, box_alpha=0.85,
                           font_properties={"size": 7}, fixed_value=40, border_pad=0.4))


def panel_geometry(ax, geo: dict) -> None:
    """Map view of each track: the flight direction (dashed) and the look direction (solid), both from the cube."""
    ax.set_xlim(-1.45, 1.45), ax.set_ylim(-4.75, 1.05)
    ax.set_aspect("equal")
    ax.axis("off")
    for k, (tr, g) in enumerate(geo.items()):
        cy = 0 if k == 0 else -2.75
        ax.add_patch(matplotlib.patches.Circle((0, cy), 0.75, fc="none", ec=fs.GRID, lw=0.6))
        ax.text(0, cy + 0.8, "N", ha="center", va="bottom", fontsize=6.5, color=fs.MUTED)
        col = fs.TRACK[tr]
        for az, ls, lab in ((g["heading"], (0, (3, 2)), "flight"), (g["look_az"], "-", "look")):
            dx, dy = np.sin(np.radians(az)) * 0.7, np.cos(np.radians(az)) * 0.7
            ax.add_patch(FancyArrowPatch((0, cy), (dx, cy + dy), arrowstyle="-|>", mutation_scale=8, color=col,
                                         lw=1.1, linestyle=ls))
        ax.text(0, cy - 0.85, f"{tr}\nheading {g['heading']:.0f}°, look {g['look_az']:.0f}°\n"
                              f"incidence {g['incidence']:.0f}°", ha="center", va="top", fontsize=6.5, linespacing=1.2)


def main() -> None:
    st = xr.open_dataset(CUBE / "silver" / "static_40m.nc")
    px = pd.read_csv(HUB / "08_deliverables" / "field_first" / "plot_pixels.csv").set_index("plot")
    box = (px.E.min() - 60, px.E.max() + 60, px.N.min() - 40, px.N.max() + 40)
    geo = look_geometry(st)
    with fs.style():
        fig = plt.figure(figsize=fs.size("full", 128))
        gs = fig.add_gridspec(2, 3, width_ratios=[1.0, 3.0, 0.95], height_ratios=[0.85, 1.55], left=0.025,
                              right=0.995, bottom=0.17, top=0.965, wspace=0.26, hspace=0.22)
        a = fig.add_subplot(gs[0, 0])
        d = fig.add_subplot(gs[1, 0])
        b = fig.add_subplot(gs[:, 1])
        c = fig.add_subplot(gs[:, 2])
        panel_locator(a)
        panel_site(b, st, None, px, box)
        panel_transect(c, st, None, px, box)
        panel_geometry(d, geo)
        d.text(0.5, 1.0, "radar geometry at the mat", transform=d.transAxes, ha="center", va="bottom", fontsize=7.5)
        d.legend(handles=[Line2D([], [], color=fs.INK, lw=1.1, ls=(0, (3, 2)), label="flight"),
                          Line2D([], [], color=fs.INK, lw=1.1, label="look")], loc="upper center", ncol=2,
                 bbox_to_anchor=(0.5, -0.02), handlelength=2.0, columnspacing=1.0)
        for ax, l in ((a, "a"), (b, "b"), (c, "c")):
            fs.panel_label(ax, l, x=0.0, y=1.005)
        fs.panel_label(d, "d", x=0.0, y=1.07)
        handles = [Line2D([], [], color=fs.CLASS["mat"], lw=1.6, label="the mat (40 m cells)"),
                   Line2D([], [], color=LAKE, lw=1.2, label="residual lake"),
                   Line2D([], [], color=fs.CLASS["grassland"], lw=0.8, label="grassland (first study's reference)"),
                   Patch(fc="none", ec=fs.CLASS["hard targets"], lw=1.0, label="hard targets (phase reference)"),
                   Line2D([], [], ls="", marker="o", ms=4, mfc=fs.INK, mec="white", label="well, mat interior (P6–P9)"),
                   Line2D([], [], ls="", marker="o", ms=4, mfc="white", mec=fs.INK, label="well, mat edge: reference (P1–P5)"),
                   Line2D([], [], ls="", marker="*", ms=8, mfc=fs.CLASS["laser"], mec="white", label="laser (P6)"),
                   Line2D([], [], color=fs.MUTED, lw=0.4, label="40 m radar grid"),
                   Line2D([], [], color=BOARDWALK, lw=1.4, label="boardwalk (traced on the orthophoto)")]
        fig.legend(handles=handles, loc="lower center", ncol=3, bbox_to_anchor=(0.5, 0.0), columnspacing=1.6,
                   handlelength=1.8)
        sheets, date = ortho_sheets()
        rec = fs.export(fig, OUT / "F1_site", provenance={
            "inputs": ["cube silver/static_40m.nc (zone, hard_target, los_*, incidence_*)",
                       *[f"cube raw/gugik_dated/{f.name} (colour-infrared orthophoto)" for f in sheets],
                       "cube raw/naturalearth/ne_50m_admin_0_countries.zip (Natural Earth, public domain)",
                       "hub 08_deliverables/field_first/plot_pixels.csv (well positions)",
                       "hub 08_deliverables/field_boardwalk_x060/boardwalk_centreline.csv"],
            "orthophoto_date": date,
            "transformations": ["orthophoto sheets read decimated (average), reprojected EPSG:2180 → EPSG:32633 "
                                "(bilinear), 2 m in (b), 0.3 m in (c); bands as delivered (NIR, red, green); 1–99 % "
                                "stretch per band", "locator in EPSG:3035", "look geometry: medians over mat cells"],
            "look_geometry": geo})
        plt.close(fig)
    print("F1 →", OUT, rec["width_mm"], "×", rec["height_mm"], "mm")


if __name__ == "__main__":
    main()

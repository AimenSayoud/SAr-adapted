"""F1, final version — the site: the floating mat, its wells, the laser at P6, the references and the radar's two looks.

Built in the figure workbook (hub NOTES, 2026-10-16): (a) the mat in its forest on a grey hillshade of the 1 m LiDAR
surface model, with the matched grassland, the hard targets and the 40 m radar grid's mat outline; (b) the transect of
wells at 1 m with the 40 m cells around it; (c) the two tracks' flight and look directions with their incidence at the
mat, computed from the cube's look vectors. Every value drawn is read from the cube or the field table; no number is
typed. Reads the cube and the hub; writes only to the hub: ``08_deliverables/article_plan/figures/final/``.

    make run E=article_fig_f1
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import xarray as xr  # noqa: E402
from matplotlib.colors import LightSource  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import FancyArrowPatch, Patch, Rectangle  # noqa: E402
from matplotlib_scalebar.scalebar import ScaleBar  # noqa: E402

from insar_wetlands import figstyle as fs  # noqa: E402
from insar_wetlands.paths import cube_root, hub_root  # noqa: E402

HUB, CUBE = hub_root(), cube_root()
OUT = HUB / "08_deliverables" / "article_plan" / "figures" / "final"
MAT_WELLS = ("P6", "P7", "P8", "P9")          # field_plots_x061: P6–P9 on the mat, P1–P5 the water-table reference
CELL = 40.0


def hillshade(z: np.ndarray, res: float) -> np.ndarray:
    """Grey relief from the surface model, lit from the north-west (the cartographic convention)."""
    return LightSource(azdeg=315, altdeg=45).hillshade(np.nan_to_num(z, nan=np.nanmedian(z)), vert_exag=2, dx=res,
                                                       dy=res)


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


def panel_site(ax, st, li, px, box) -> None:
    x, y = st.x.values, st.y.values
    rm, cm = np.where(st.zone.values == 1)
    pad = 14
    xlim = (x[max(cm.min() - pad, 0)] - CELL / 2, x[min(cm.max() + pad, len(x) - 1)] + CELL / 2)
    ylim = (y[min(rm.max() + pad, len(y) - 1)] - CELL / 2, y[max(rm.min() - pad, 0)] + CELL / 2)
    sub = li.sel(x=slice(*xlim), y=slice(ylim[1], ylim[0]))
    hs = hillshade(sub.dsm_m.values, 1.0)
    ax.imshow(hs, cmap="Greys_r", vmin=0, vmax=1, alpha=0.6,
              extent=[float(sub.x[0]) - 0.5, float(sub.x[-1]) + 0.5, float(sub.y[-1]) - 0.5, float(sub.y[0]) + 0.5],
              origin="upper", interpolation="none", zorder=0)
    ext = [x[0] - CELL / 2, x[-1] + CELL / 2, y[-1] - CELL / 2, y[0] + CELL / 2]
    water = np.where(st.zone.values == 2, 1.0, np.nan)
    ax.imshow(water, extent=ext, cmap=matplotlib.colors.ListedColormap([fs.CLASS["open water"]]), alpha=0.85,
              origin="upper", interpolation="nearest", zorder=1)
    grass = st.zone.values == 3
    gy, gx = np.where(grass)
    for r, c in zip(gy, gx):
        ax.add_patch(Rectangle((x[c] - CELL / 2, y[r] - CELL / 2), CELL, CELL, fc="none", ec=fs.CLASS["grassland"],
                               hatch="////", lw=0, zorder=2))
    cell_outline(ax, x, y, grass, colors=fs.CLASS["grassland"], linewidths=0.6, zorder=2)
    cell_outline(ax, x, y, st.zone.values == 1, colors=fs.CLASS["mat"], linewidths=1.6, zorder=3)
    hr, hc = np.where(st.hard_target.values.astype(bool))
    for r, c in zip(hr, hc):
        ax.add_patch(Rectangle((x[c] - CELL / 2, y[r] - CELL / 2), CELL, CELL, fc="none", ec=fs.CLASS["hard targets"],
                               lw=1.0, zorder=4))
    ref = px[~px.index.isin(MAT_WELLS)]
    mat = px[px.index.isin(MAT_WELLS)]
    ax.plot(ref.E, ref.N, "o", ms=2.6, mfc="white", mec=fs.INK, mew=0.6, zorder=5)
    ax.plot(mat.E, mat.N, "o", ms=2.6, mfc=fs.INK, mec="white", mew=0.4, zorder=5)
    ax.plot(px.loc["P6", "E"], px.loc["P6", "N"], "*", ms=8, mfc=fs.CLASS["laser"], mec="white", mew=0.5, zorder=6)
    ax.add_patch(Rectangle((box[0], box[2]), box[1] - box[0], box[3] - box[2], fc="none", ec=fs.INK, lw=0.8,
                           ls="--", zorder=6))
    ax.text(box[1] + 8, box[3], "(b)", fontsize=7, va="top", zorder=6)
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
                                 color=fs.INK, lw=0.9, zorder=7))
    ax.text(nx, ny + 0.10 * (y1 - y0), "N", ha="center", va="bottom", fontsize=7, zorder=7)


def panel_transect(ax, st, li, px, box) -> None:
    sub = li.sel(x=slice(box[0], box[1]), y=slice(box[3], box[2]))
    ax.imshow(hillshade(sub.dsm_m.values, 1.0), cmap="Greys_r", vmin=0, vmax=1, alpha=0.75,
              extent=[float(sub.x[0]) - 0.5, float(sub.x[-1]) + 0.5, float(sub.y[-1]) - 0.5, float(sub.y[0]) + 0.5],
              origin="upper", interpolation="none", zorder=0)
    x, y = st.x.values, st.y.values
    for xv in np.concatenate([x - CELL / 2, [x[-1] + CELL / 2]]):
        if box[0] <= xv <= box[1]:
            ax.axvline(xv, color=fs.MUTED, lw=0.4, zorder=1)
    for yv in np.concatenate([y + CELL / 2, [y[-1] - CELL / 2]]):
        if box[2] <= yv <= box[3]:
            ax.axhline(yv, color=fs.MUTED, lw=0.4, zorder=1)
    cell_outline(ax, x, y, st.zone.values == 1, colors=fs.CLASS["mat"], linewidths=1.6, zorder=2)
    left = {"P9", "P7", "P5"}                      # wells west of their neighbours: label on the west side
    for p, row in px.iterrows():
        on = p in MAT_WELLS
        if p == "P6":
            ax.plot(row.E, row.N, "*", ms=10, mfc=fs.CLASS["laser"], mec="white", mew=0.5, zorder=5)
        else:
            ax.plot(row.E, row.N, "o", ms=4, mfc=fs.INK if on else "white", mec=fs.INK if not on else "white",
                    mew=0.7, zorder=4)
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
    ax.set_xlim(-1.3, 1.3), ax.set_ylim(-4.4, 1.1)
    ax.set_aspect("equal")
    ax.axis("off")
    for k, (tr, g) in enumerate(geo.items()):
        cy = 0 if k == 0 else -2.55
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
    li = xr.open_dataset(CUBE / "silver" / "lidar_1m.nc")
    px = pd.read_csv(HUB / "08_deliverables" / "field_first" / "plot_pixels.csv").set_index("plot")
    box = (px.E.min() - 60, px.E.max() + 60, px.N.min() - 40, px.N.max() + 40)
    geo = look_geometry(st)
    with fs.style(), matplotlib.rc_context({"hatch.linewidth": 0.5}):
        fig = plt.figure(figsize=fs.size("full", 112))
        gs = fig.add_gridspec(1, 3, width_ratios=[3.3, 0.95, 1.05], left=0.075, right=0.995, bottom=0.19, top=0.975,
                              wspace=0.1)
        a, b, c = fig.add_subplot(gs[0]), fig.add_subplot(gs[1]), fig.add_subplot(gs[2])
        panel_site(a, st, li, px, box)
        panel_transect(b, st, li, px, box)
        panel_geometry(c, geo)
        c.set_title("radar geometry at the mat", fontsize=7.5, pad=2)
        c.legend(handles=[Line2D([], [], color=fs.INK, lw=1.1, ls=(0, (3, 2)), label="flight"),
                          Line2D([], [], color=fs.INK, lw=1.1, label="look")], loc="lower center", ncol=2,
                 bbox_to_anchor=(0.5, -0.02), handlelength=2.0, columnspacing=1.0)
        for ax, l in ((a, "a"), (b, "b"), (c, "c")):
            fs.panel_label(ax, l, x=0.0, y=1.005)
        handles = [Line2D([], [], color=fs.CLASS["mat"], lw=1.6, label="floating mat (40 m cells)"),
                   Patch(fc=fs.CLASS["open water"], label="residual lake"),
                   Patch(fc="none", ec=fs.CLASS["grassland"], hatch="////", label="matched grassland (reference)"),
                   Patch(fc="none", ec=fs.CLASS["hard targets"], lw=1.0, label="hard targets (phase reference)"),
                   Line2D([], [], ls="", marker="o", ms=4, mfc=fs.INK, mec="white", label="well, mat interior (P6–P9)"),
                   Line2D([], [], ls="", marker="o", ms=4, mfc="white", mec=fs.INK, label="well, mat edge: reference (P1–P5)"),
                   Line2D([], [], ls="", marker="*", ms=8, mfc=fs.CLASS["laser"], mec="white", label="laser (P6)"),
                   Line2D([], [], color=fs.MUTED, lw=0.4, label="40 m radar grid")]
        fig.legend(handles=handles, loc="lower center", ncol=4, bbox_to_anchor=(0.5, 0.0), columnspacing=1.4,
                   handlelength=1.8)
        rec = fs.export(fig, OUT / "F1_site", provenance={
            "inputs": ["cube silver/static_40m.nc (zone, hard_target, los_*, incidence_*)",
                       "cube silver/lidar_1m.nc (dsm_m: GUGiK NMPT, 1 m)",
                       "hub 08_deliverables/field_first/plot_pixels.csv (well positions)"],
            "transformations": ["hillshade of the DSM, azimuth 315°, altitude 45°, vertical exaggeration 2",
                                "both panels at the DSM's own 1 m (no resampling)", "look geometry: medians over mat cells"],
            "look_geometry": geo})
        plt.close(fig)
    print("F1 →", OUT, rec["width_mm"], "×", rec["height_mm"], "mm")


if __name__ == "__main__":
    main()

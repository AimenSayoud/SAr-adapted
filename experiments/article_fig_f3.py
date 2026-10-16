"""F3, final version — a schematic of what C-band sees on a floating mat and which pairs keep its motion. No data and
no numbers: the quantities it names are measured in Figs. 2, 4 and 5 and in §§4–5.

(a) A cross-section at two dates. The water rises; the mat, anchored at its edge, rises by a fraction of that
(buoyancy). C-band rays at the tracks' incidence scatter from the mat surface, so a pair's phase follows the surface;
only where stems stand in open water (the lake's reed margin) can a double bounce carry the water level. Coherence
changes with the wetness of the surface. (b) A year of acquisitions: consecutive short pairs follow the motion, pairs
that skip an acquisition and long pairs lose it, and across the winter the chain's level must come from outside the
radar (§5.2; the estimator of §6 takes it from the weather).
Writes only to the hub ``08_deliverables/article_plan/figures/final/``.

    make run E=article_fig_f3
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Arc, FancyArrowPatch, Polygon, Rectangle  # noqa: E402

from insar_wetlands import figstyle as fs  # noqa: E402
from insar_wetlands.paths import hub_root  # noqa: E402

OUT = hub_root() / "08_deliverables" / "article_plan" / "figures" / "final"
WATER, WATER_FILL, PEAT, STEM, BANK = "#2C7FB8", "#dcecf7", "#8a6d3b", "#3f7d3a", "#c9c3b6"
INC = 35.0                                            # schematic incidence, between the two tracks'


def arrow(ax, p0, p1, color=fs.INK, lw=0.9, style="-|>", ls="-", ms=8):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle=style, mutation_scale=ms, color=color, lw=lw, linestyle=ls,
                                 shrinkA=0, shrinkB=0))


def stems(ax, x0, x1, base, n, h=0.32, color=STEM, lw=0.8):
    for xv in np.linspace(x0, x1, n):
        ax.plot([xv, xv + 0.03], [base, base + h], color=color, lw=lw, solid_capstyle="round")


def panel_section(ax) -> None:
    ax.set_xlim(0, 20), ax.set_ylim(0.6, 7.75)
    ax.set_aspect("equal")
    ax.axis("off")
    w1, w2 = 2.2, 3.3                                  # water level, date 1 and date 2
    rise = 0.4 * (w2 - w1)                             # the mat rises by a fraction of the water's rise (drawn large)
    ax.add_patch(Rectangle((2.6, 0), 15.2, w2, fc=WATER_FILL, ec="none", zorder=0))
    ax.plot([2.6, 17.8], [w2, w2], color=WATER, lw=1.0, zorder=1)
    ax.plot([2.6, 17.8], [w1, w1], color=WATER, lw=0.8, ls=(0, (3, 2)), zorder=1)
    ax.add_patch(Polygon([(0, 0), (0, 3.9), (1.6, 3.9), (3.0, 2.4), (3.4, 0)], fc=BANK, ec="none", zorder=2))
    ax.text(0.8, 0.6, "firm\nground", ha="center", fontsize=6.5, color=fs.MUTED, zorder=3)
    top2 = w2 + 0.35                                   # a floating mat sits low: small freeboard, base below both
    base = top2 - 1.95                                 # water levels at either date (it floats on date 1 too)
    mat = Polygon([(3.05, top2), (12.6, top2), (12.6, base), (3.35, base)], fc=PEAT, ec="none", zorder=3)
    ax.add_patch(mat)
    ax.plot([3.05, 13.2], [top2 - rise, top2 - rise], color=fs.INK, lw=0.7, ls=(0, (2, 1.5)), zorder=4)
    ax.plot([12.6, 13.2], [top2, top2], color=fs.INK, lw=0.7, zorder=4)
    stems(ax, 3.3, 12.3, top2, 30, h=0.3)
    ax.text(7.9, base + 0.45, "floating mat: roots and peat", ha="center", va="center", fontsize=6.5, color="white",
            zorder=5)
    stems(ax, 14.0, 15.2, w2 - 0.9, 8, h=1.6, lw=0.9)                       # reeds in the lake margin
    ax.text(16.2, 1.0, "residual lake", ha="center", fontsize=6.5, color=WATER)
    # changes between the dates
    arrow(ax, (17.3, w1), (17.3, w2), color=WATER, style="<|-|>", ms=6)
    ax.text(17.55, (w1 + w2) / 2, "water\nrises", ha="left", va="center", fontsize=6.5, color=WATER)
    arrow(ax, (13.05, top2 - rise), (13.05, top2), color=fs.INK, style="<|-|>", ms=6)
    ax.text(12.5, top2 + 0.2, "mat surface rises\nby a fraction of it", ha="right", va="bottom", fontsize=6.5)
    ax.text(3.3, top2 - rise - 0.05, "mat surface, date 1", ha="left", va="top", fontsize=6, color="white", zorder=5)
    # radar rays at the incidence angle, from the upper left
    dx, dy = np.sin(np.radians(INC)), -np.cos(np.radians(INC))
    for xe in (5.2, 7.0, 8.8):
        L = 3.3
        arrow(ax, (xe - dx * L, top2 - dy * L * -1 if False else top2 + abs(dy) * L), (xe, top2 + 0.05), lw=0.8)
    ax.text(2.4, top2 + 3.05, "C-band", ha="left", va="bottom", fontsize=7, fontweight="bold")
    ax.annotate("surface scattering: the pair's phase\nfollows the mat surface", (7.0, top2 + 0.1),
                xytext=(9.3, top2 + 2.75), fontsize=6.5, ha="left", va="center",
                arrowprops=dict(arrowstyle="-", color=fs.INK, lw=0.5))
    # double bounce at the reed margin
    xb = 14.3
    arrow(ax, (xb - dx * 3.2, w2 + abs(dy) * 3.2), (xb, w2 + 0.02), lw=0.8)
    arrow(ax, (xb, w2 + 0.02), (xb + 0.45, w2 + 0.65), lw=0.6, ms=5)
    ax.annotate("double bounce, water and stems:\nwould follow the water level\n(too rare at the plots to read)",
                (xb + 0.3, w2 + 0.5), xytext=(15.0, top2 + 2.2), fontsize=6.5, ha="left", va="center",
                arrowprops=dict(arrowstyle="-", color=fs.INK, lw=0.5))
    # coherence and wetness
    for xv in np.linspace(3.6, 6.4, 7):
        ax.plot(xv, top2 + 0.08, "o", ms=1.8, color=WATER, zorder=6)
    ax.annotate("coherence: changes with the\nwetness of the surface", (4.4, top2 + 0.1), xytext=(0.2, top2 + 1.7),
                fontsize=6.5, ha="left", va="center", arrowprops=dict(arrowstyle="-", color=fs.INK, lw=0.5))
    ax.plot([], [], color=WATER, lw=0.8, ls=(0, (3, 2)), label="date 1")
    ax.plot([], [], color=WATER, lw=1.0, label="date 2")
    ax.legend(loc="lower right", bbox_to_anchor=(1.0, -0.02), fontsize=6.5, handlelength=1.6, ncol=2,
              title="water level", title_fontsize=6.5)


def panel_pairs(ax) -> None:
    n = 31
    t = np.arange(n)
    winter = (17, 23)
    ax.set_xlim(-0.8, n - 0.2), ax.set_ylim(-1.6, 6.2)
    ax.axis("off")
    ax.add_patch(Rectangle((winter[0] - 0.5, -0.4), winter[1] - winter[0] + 1, 6.3, fc="#eef1f4", ec="none", zorder=0))
    ax.text(winter[1] + 0.3, 5.4, "winter: snow and frost", ha="right", fontsize=6.5, color=fs.MUTED)
    col_in = [fs.MUTED if winter[0] <= k <= winter[1] else fs.INK for k in t]
    ax.scatter(t, np.zeros(n), s=9, c=col_in, zorder=3)
    ax.text(-0.6, -0.65, "acquisitions of one track", ha="left", fontsize=6.5, color=fs.MUTED)

    def arc(a, b, color, ls="-", lw=1.2):
        w = b - a
        ax.add_patch(Arc(((a + b) / 2, 0), w, w * 0.9, theta1=0, theta2=180, color=color, lw=lw, ls=ls, zorder=2))

    for k in range(1, 5):
        arc(k, k + 1, fs.CLASS["mat"])
    ax.text(3.0, 1.0, "consecutive pairs:\nfollow the motion", ha="center", va="bottom", fontsize=6.5,
            color=fs.CLASS["mat"])
    arc(7, 9, fs.CLASS["laser"], ls=(0, (2.5, 1.5)), lw=1.0)
    arc(8, 10, fs.CLASS["laser"], ls=(0, (2.5, 1.5)), lw=1.0)
    ax.text(8.5, 1.35, "skipping an\nacquisition", ha="center", va="bottom", fontsize=6.5, color=fs.CLASS["laser"])
    arc(11, 16, fs.CLASS["laser"], ls=(0, (2.5, 1.5)), lw=1.0)
    arc(12, 21, fs.CLASS["laser"], ls=(0, (2.5, 1.5)), lw=1.0)
    ax.text(16.5, 4.15, "long pairs: lose the motion", ha="center", va="bottom", fontsize=6.5,
            color=fs.CLASS["laser"])
    for k in range(winter[0] - 1, winter[1] + 1):
        arc(k, k + 1, fs.MUTED, lw=0.8)
    arrow(ax, (winter[0] - 1, -1.15), (winter[1] + 1, -1.15), color="#8250df", lw=1.0)
    ax.text((winter[0] + winter[1]) / 2, -1.45, "across the winter the level must come from outside the radar",
            ha="center", va="top", fontsize=6.5, color="#8250df")
    for k in range(winter[1] + 1, n - 1):
        arc(k, k + 1, fs.CLASS["mat"])
    ax.text(27.0, 0.75, "the chain\nresumes", ha="center", va="bottom", fontsize=6.5, color=fs.CLASS["mat"])


def main() -> None:
    with fs.style():
        fig = plt.figure(figsize=fs.size("full", 118))
        a = fig.add_axes([0.01, 0.37, 0.98, 0.6])
        b = fig.add_axes([0.01, 0.02, 0.98, 0.31])
        panel_section(a)
        panel_pairs(b)
        fs.panel_label(a, "a", x=0.02, y=0.97)
        fs.panel_label(b, "b", x=0.02, y=0.97)
        rec = fs.export(fig, OUT / "F3_schematic", provenance={"inputs": [], "kind": "schematic, no data",
                                                                  "incidence_drawn_deg": INC})
        plt.close(fig)
    print("F3 →", OUT, rec["width_mm"], "×", rec["height_mm"], "mm")


if __name__ == "__main__":
    main()

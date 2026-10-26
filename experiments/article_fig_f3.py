"""F3, final version (round 3) — a schematic of what C-band sees on a floating mat and which pairs keep its motion.
No data: the quantities it names are measured in Figs. 2, 4 and 5 and in §§4–5; the only value printed is the
wavelength's quarter, from the library's constant.

(a) A cross-section at two dates, dashed = date 1, solid = date 2 for both the mat and the water. The water rises; the
mat, held at its firm edge, rises by a fraction of it and more toward its interior (§4.4), so its surface ends nearer
the water. Numbered markers and a key replace leader lines: ① surface scattering — the pair's phase follows the mat
surface; ② double bounce on stems standing in open water — would follow the water level; ③ coherence changes with the
surface's wetness. An inset projects the vertical rise d on the line of sight at incidence θ (d·cos θ).
(b) One track over a year. Top: the pairs — consecutive pairs follow the motion, pairs that skip an acquisition and
long pairs lose it (why: c), winter pairs are still made. Bottom: a sketch of the level — the chain of consecutive
pairs keeps close to it through the season; across the winter its noise adds up (the band widens), so the level has to
be re-set from outside the radar (§6: from the weather), itself with an error. (c) Why the longer pairs lose it: a
pair's phase gives the line-of-sight change only within ±λ/4 and wraps beyond it (X-127 found that ambiguity enough to
explain the loss); shaded beyond ±λ/4 as in Fig. 4. Briefs in the hub ``article_plan/figures/briefs/F3.md``.
Writes only to the hub ``08_deliverables/article_plan/figures/final/``.

    make run E=article_fig_f3
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Arc, FancyArrowPatch, Polygon, Rectangle  # noqa: E402

from insar_wetlands import figstyle as fs  # noqa: E402
from insar_wetlands.cube.core import WAVELENGTH_MM  # noqa: E402
from insar_wetlands.paths import hub_root  # noqa: E402

OUT = hub_root() / "08_deliverables" / "article_plan" / "figures" / "final"
WATER, WATER_FILL, PEAT, STEM, BANK = "#2C7FB8", "#dcecf7", "#8a6d3b", "#3f7d3a", "#c9c3b6"
LOS, LOST, OUTSIDE, SHADE = "#0072B2", fs.CLASS["laser"], "#8250df", "#eef1f4"
INC = 35.0                                            # schematic incidence, between the two tracks'
Q = WAVELENGTH_MM / 4                                 # a quarter wavelength along the line of sight, mm


def arrow(ax, p0, p1, color=fs.INK, lw=0.9, style="-|>", ls="-", ms=8, z=6):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle=style, mutation_scale=ms, color=color, lw=lw, linestyle=ls,
                                 shrinkA=0, shrinkB=0, zorder=z))


def badge(ax, x, y, n, transform=None):
    """A numbered marker: a white disc with the number, the same in the drawing and in the key."""
    kw = {} if transform is None else {"transform": transform}
    ax.plot([x], [y], "o", ms=9, mfc="white", mec=fs.INK, mew=0.8, zorder=9, clip_on=False, **kw)
    ax.text(x, y, str(n), ha="center", va="center", fontsize=7, fontweight="bold", zorder=10, **kw)


def stems(ax, xs, base, h, color=STEM, lw=0.7, z=4):
    for xv, b in zip(xs, np.broadcast_to(base, np.shape(xs))):
        ax.plot([xv, xv + 0.03], [b, b + h], color=color, lw=lw, solid_capstyle="round", zorder=z)


def panel_section(ax) -> None:
    ax.set_xlim(0, 20), ax.set_ylim(0, 7.4)
    ax.set_aspect("equal")
    ax.axis("off")
    x0, x1 = 2.9, 13.0                                 # the mat: held at the bank (x0), free toward the lake (x1)
    t1, thick = 3.55, 1.7                              # date-1 surface and the mat's thickness (vertical exaggerated)
    w1, w2 = 2.7, 3.35                                 # water level, date 1 and date 2

    def t2(x):                                         # date-2 surface: a small rise at the edge, more inside (§4.4)
        x = np.asarray(x, float)
        return t1 + 0.15 + 0.38 * (1 - np.exp(-(x - x0) / 3.0))

    # water, bank, mat at date 2 (filled) and its outline at date 1 (dashed)
    ax.add_patch(Rectangle((x0, 0), 20 - x0, w2, fc=WATER_FILL, ec="none", zorder=0))
    ax.plot([x1, 20], [w2, w2], color=WATER, lw=1.1, zorder=1)
    ax.plot([x1, 20], [w1, w1], color=WATER, lw=1.0, ls=(0, (3, 2)), zorder=1)
    ax.add_patch(Polygon([(0, 0), (0, 4.3), (1.8, 4.3), (x0 + 0.05, 3.1), (x0 + 0.25, 0)], fc=BANK, ec="none",
                         zorder=2))
    ax.text(1.35, 0.45, "firm\nground", ha="center", fontsize=7, color=fs.INK, linespacing=1.0, zorder=3)
    xs = np.linspace(x0, x1, 200)
    top = t2(xs)
    ax.add_patch(Polygon(np.r_[np.c_[xs, top], np.c_[xs[::-1], top[::-1] - thick]], fc=PEAT, ec="none", zorder=3))
    ax.plot([x0, x1, x1, x0], [t1, t1, t1 - thick, t1 - thick], color=fs.INK, lw=0.9, ls=(0, (3, 2)), zorder=5)
    stems(ax, np.linspace(x0 + 0.3, x1 - 0.3, 26), t2(np.linspace(x0 + 0.3, x1 - 0.3, 26)), 0.22)
    ax.plot([x0], [t2(x0)], "o", ms=3.5, color=fs.INK, zorder=6)                       # where the mat is held
    ax.text(7.9, 2.95, "floating mat (roots and peat)", ha="center", va="center", fontsize=7, color="white",
            zorder=6)
    ax.text(7.9, 2.5, "held at its edge, it rises more toward its interior  →", ha="center", va="center",
            fontsize=7, color="white", zorder=6)
    # the rise: small at the edge, larger inside; the inset below names it d
    xr = 12.35
    arrow(ax, (xr, t1), (xr, float(t2(xr))), style="-|>", ms=6, lw=0.9, z=7)
    ax.text(xr - 0.12, (t1 + float(t2(xr))) / 2 + 0.02, "d", ha="right", va="center", fontsize=7,
            fontstyle="italic", color="white", zorder=7)
    # the water: direct labels and its rise
    ax.text(15.9, w2 + 0.06, "water, date 2", ha="left", va="bottom", fontsize=7, color=WATER)
    ax.text(15.9, w1 + 0.06, "water, date 1", ha="left", va="bottom", fontsize=7, color=WATER)
    arrow(ax, (19.45, w1), (19.45, w2), color=WATER, style="-|>", ms=6, lw=1.0)
    ax.text(19.3, (w1 + w2) / 2 - 0.15, "rises", ha="right", va="center", fontsize=7, color=WATER)
    ax.text(17.2, 1.9, "open water (residual lake)", ha="center", va="center", fontsize=7, color=WATER)
    stems(ax, np.linspace(14.1, 15.3, 7), w2 - 0.7, 1.55, lw=0.9)                       # reeds in open water

    # radar rays at the incidence angle, from the upper left (direction of travel: down and right)
    dx, dy = np.sin(np.radians(INC)), np.cos(np.radians(INC))
    L = 2.7
    for xe in (7.6, 9.0, 10.4):
        ye = float(t2(xe)) + 0.24
        arrow(ax, (xe - dx * L, ye + dy * L), (xe, ye), lw=0.9)
    ax.text(7.6 - dx * L - 0.1, float(t2(7.6)) + 0.24 + dy * L + 0.12, "C-band radar", ha="left", va="bottom",
            fontsize=7, fontweight="bold")
    xt, yt = 10.4, float(t2(10.4)) + 0.24                # θ at one ray, from the vertical
    ax.plot([xt, xt], [yt, yt + 1.35], color=fs.MUTED, lw=0.6, ls=(0, (2, 1.5)), zorder=5)
    ax.add_patch(Arc((xt, yt), 1.9, 1.9, theta1=90, theta2=90 + INC, color=fs.MUTED, lw=0.7, zorder=5))
    ax.text(xt - 0.3, yt + 1.15, "θ", ha="center", va="center", fontsize=7.5, color=fs.INK)
    badge(ax, 11.3, float(t2(11.3)) + 0.55, 1)
    # double bounce: down to the water, across to a stem, back toward the radar
    xb = 14.0
    arrow(ax, (xb - dx * 2.6, w2 + dy * 2.6), (xb, w2 + 0.02), lw=0.9)
    arrow(ax, (xb, w2 + 0.02), (xb + 0.42, w2 + 0.6), lw=0.7, ms=5)
    arrow(ax, (xb + 0.42, w2 + 0.6), (xb + 0.42 - dx * 1.3, w2 + 0.6 + dy * 1.3), lw=0.7, ms=5, ls=(0, (2, 1.2)))
    badge(ax, 15.75, w2 + 0.85, 2)
    # wetness near the edge, where the surface ends nearest the water
    for xv in np.linspace(3.3, 5.4, 8):
        ax.plot(xv, float(t2(xv)) + 0.05, "o", ms=2.0, color=WATER, zorder=7)
    badge(ax, 4.35, float(t2(4.35)) + 0.62, 3)

    # the key, top left (in place of leader lines)
    key = [(1, "surface scattering: the pair's phase", "follows the mat surface"),
           (2, "double bounce on stems in open water:", "would follow the water level\n(rare at the plots)"),
           (3, "coherence: changes with the wetness", "of the surface")]
    y = 6.95
    for n, l1, l2 in key:
        badge(ax, 0.45, y - 0.17, n)
        ax.text(0.8, y, l1 + "\n" + l2, ha="left", va="top", fontsize=7, linespacing=1.15)
        y -= 0.78 + 0.3 * l2.count("\n")

    # inset: the vertical rise seen along the line of sight
    ox, oy, d = 17.75, 4.75, 1.55
    u = np.array([-dx, dy])                            # toward the radar
    p = np.array([ox, oy]) + d * dy * u                # the projection of the rise on the line of sight
    ax.plot([ox, ox - 2.2 * dx], [oy, oy + 2.2 * dy], color=fs.MUTED, lw=0.6, ls=(0, (2, 1.5)), zorder=5)
    arrow(ax, (ox, oy), (ox, oy + d), lw=1.1, ms=7)
    arrow(ax, (ox, oy), tuple(p), color=LOS, lw=1.6, ms=7)
    ax.plot([ox, p[0]], [oy + d, p[1]], color=fs.MUTED, lw=0.6, ls=(0, (1, 1.2)), zorder=5)
    ax.add_patch(Arc((ox, oy), 0.9, 0.9, theta1=90, theta2=90 + INC, color=fs.MUTED, lw=0.7, zorder=5))
    ax.text(ox - 0.2, oy + 0.7, "θ", ha="center", va="center", fontsize=7.5)
    ax.text(ox + 0.12, oy + d - 0.05, "d, vertical rise", ha="left", va="top", fontsize=7)
    ax.text(p[0] - 0.12, p[1] - 0.05, "d · cos θ:\nwhat the pair\nmeasures", ha="right", va="top", fontsize=7,
            color=LOS, linespacing=1.1)
    ax.text(ox - 2.2 * dx - 0.05, oy + 2.2 * dy + 0.05, "line of sight", ha="right", va="bottom", fontsize=7,
            color=fs.MUTED)
    # what dashed and solid mean, and the scale
    hs = [Line2D([], [], color=fs.INK, lw=0.9, ls=(0, (3, 2))), Line2D([], [], color=fs.INK, lw=1.1)]
    leg = ax.legend(hs, ["date 1 (earlier)", "date 2 (later)"], loc="lower right", bbox_to_anchor=(1.0, 0.0),
                    fontsize=7, handlelength=2.2, ncol=2, title="mat and water", title_fontsize=7,
                    borderaxespad=0.2)
    leg._legend_box.align = "left"
    ax.text(x0 + 0.4, 0.15, "not to scale; vertical exaggerated", ha="left", va="bottom", fontsize=7,
            color=fs.MUTED, zorder=8)


# --- (b) one track over a year --------------------------------------------------------------------------------------
MONTHS = ["May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun"]
MONTH_START = np.cumsum([0, 31, 30, 31, 31, 30, 31, 30, 31, 31, 28, 31, 30, 31, 30])  # from 1 May
WINTER = (MONTH_START[7], MONTH_START[11])            # 1 December to 31 March
STEP = 12                                             # days between acquisitions (schematic)


def _months(ax, labels: bool) -> None:
    ax.set_xlim(-6, MONTH_START[-1])
    ax.set_xticks(MONTH_START[:-1], [""] * len(MONTHS))
    if labels:
        ax.set_xticks((MONTH_START[:-1] + MONTH_START[1:]) / 2, MONTHS, minor=True)
        ax.tick_params(axis="x", which="minor", length=0, labelsize=7)
    ax.tick_params(axis="x", which="major", length=3)
    ax.spines[["left", "top", "right"]].set_visible(False)
    ax.set_yticks([])


def panel_pairs(ax, lv) -> None:
    t = np.arange(0, MONTH_START[-1], STEP)
    n = len(t)
    in_winter = (t >= WINTER[0]) & (t < WINTER[1])
    for a in (ax, lv):
        a.axvspan(*WINTER, color=SHADE, lw=0, zorder=0)
    _months(ax, labels=False)
    _months(lv, labels=True)
    ax.set_ylim(-0.25, 2.25)

    def arc(i, j, color, ls="-", lw=1.1):
        w = t[j] - t[i]
        h = 0.32 + 0.22 * (j - i - 1)
        ax.add_patch(Arc(((t[i] + t[j]) / 2, 0), w, 2 * h, theta1=0, theta2=180, color=color, lw=lw, ls=ls,
                         zorder=2))

    for i in range(n - 1):
        winter_pair = in_winter[i] or in_winter[i + 1]
        arc(i, i + 1, fs.MUTED if winter_pair else fs.CLASS["mat"], lw=0.8 if winter_pair else 1.1)
    ax.scatter(t, np.zeros(n), s=7, c=[fs.MUTED if w else fs.INK for w in in_winter], zorder=3, lw=0)
    dash = (0, (2.5, 1.5))
    for i in (5, 6):
        arc(i, i + 2, LOST, ls=dash, lw=1.0)
    for i in (11, 12):
        arc(i, i + 5, LOST, ls=dash, lw=1.0)
    ax.text(t[1] + 6, 0.5, "consecutive pairs:\nfollow the motion", ha="center", va="bottom", fontsize=7,
            color=fs.CLASS["mat"], linespacing=1.1)
    ax.text((t[5] + t[8]) / 2, 0.86, "skip one:\nlose it", ha="center", va="bottom", fontsize=7, color=LOST,
            linespacing=1.1)
    ax.text((t[11] + t[17]) / 2, 1.5, "long pairs: lose it (why: c)", ha="center", va="bottom", fontsize=7,
            color=LOST)
    ax.text(sum(WINTER) / 2, 0.5, "winter: pairs still made,\nbut snow and frost", ha="center", va="bottom",
            fontsize=7, color=fs.MUTED, linespacing=1.1)
    ax.text(MONTH_START[12] + 30, 0.5, "consecutive\nagain", ha="center", va="bottom", fontsize=7,
            color=fs.CLASS["mat"], linespacing=1.1)
    ax.text(-6, 2.2, "the pairs of one track", ha="left", va="top", fontsize=7, color=fs.INK, fontweight="bold")

    # the level: truth, the chain of consecutive pairs and how far it can be off (±2 σ √steps)
    tt = np.linspace(0, MONTH_START[-1], 600)
    truth = lambda x: 1.0 * np.cos(2 * np.pi * (np.asarray(x) - 255) / 365)      # high in winter, low in summer
    lv.plot(tt, truth(tt), color=fs.INK, lw=1.0, zorder=4)
    rng = np.random.default_rng(5)
    sig = np.where(in_winter[1:] | in_winter[:-1], 0.12, 0.035)                    # per pair, schematic
    spring = int(np.argmax(t >= WINTER[1]))
    est, band = [truth(t[0])], [0.0]
    for k in range(1, spring):
        est.append(est[-1] + truth(t[k]) - truth(t[k - 1]) + rng.normal(0, sig[k - 1]))
        band.append(np.sqrt(band[-1] ** 2 + sig[k - 1] ** 2))
    est, band = np.array(est), 2 * np.array(band)
    season = t[:spring] < WINTER[0]
    last_in = int(np.flatnonzero(season)[-1])
    for sl, col, ls in ((slice(0, last_in + 1), fs.CLASS["mat"], "-"), (slice(last_in, spring), fs.MUTED, "--")):
        lv.fill_between(t[sl], est[sl] - band[sl], est[sl] + band[sl], color=col, alpha=0.18, lw=0, zorder=2)
        lv.plot(t[sl], est[sl], color=col, lw=1.1, ls=ls, zorder=5)
    # re-set at spring from outside the radar, with its own error, then the chain again
    r0 = truth(t[spring]) + 0.12
    arrow(lv, (t[spring - 1], est[-1]), (t[spring], r0), color=OUTSIDE, lw=1.0, ms=6)
    lv.errorbar([t[spring]], [r0], yerr=0.18, fmt="o", ms=3.5, mfc="white", mec=OUTSIDE, ecolor=OUTSIDE,
                elinewidth=0.9, capsize=2, zorder=7)
    e2, b2 = [r0], [0.09]
    for k in range(spring + 1, n):
        e2.append(e2[-1] + truth(t[k]) - truth(t[k - 1]) + rng.normal(0, sig[k - 1]))
        b2.append(np.sqrt(b2[-1] ** 2 + sig[k - 1] ** 2))
    e2, b2 = np.array(e2), 2 * np.array(b2)
    lv.fill_between(t[spring:], e2 - b2, e2 + b2, color=fs.CLASS["mat"], alpha=0.18, lw=0, zorder=2)
    lv.plot(t[spring:], e2, color=fs.CLASS["mat"], lw=1.1, zorder=5)
    lv.set_ylim(-1.3, 2.4)
    lv.set_ylabel("mat level\n(sketch)", fontsize=7, linespacing=1.1)
    lv.spines["left"].set_visible(True)
    lv.text(4, -0.2, "true level", ha="left", va="bottom", fontsize=7, color=fs.INK)
    lv.text(132, 0.55, "chain of consecutive pairs:\nclose to it within the season", ha="center", va="bottom",
            fontsize=7, color=fs.CLASS["mat"], linespacing=1.1)
    lv.text(sum(WINTER) / 2, -1.2, "noise adds up\nalong the chain", ha="center", va="bottom", fontsize=7,
            color=fs.MUTED, linespacing=1.1)
    lv.text(t[spring] - 3, 1.65, "level re-set from outside the radar\n(§6: from the weather), with an error",
            ha="right", va="bottom", fontsize=7, color=OUTSIDE, linespacing=1.1)


def panel_wrap(ax) -> None:
    """What one pair's phase returns against the true line-of-sight change (two-way phase 4πd/λ wrapped to ±π)."""
    lim = 2.6 * Q
    x = np.linspace(-lim, lim, 2001)
    y = (x + Q) % (2 * Q) - Q
    y[np.r_[False, np.abs(np.diff(y)) > Q]] = np.nan   # no line across a wrap
    for sgn in (-1, 1):
        ax.axvspan(sgn * Q, sgn * lim, color=SHADE, lw=0, zorder=0)
    ax.plot([-lim, lim], [-lim, lim], color=fs.MUTED, lw=0.7, ls=(0, (3, 2)), zorder=1)
    ax.plot(x, y, color=fs.INK, lw=1.2, zorder=2)
    xe = 1.45 * Q                                      # one change beyond λ/4: returned smaller and with the wrong sign
    ye = (xe + Q) % (2 * Q) - Q
    ax.annotate("", (xe, ye), (xe, xe), arrowprops=dict(arrowstyle="-|>", color=LOST, lw=0.9, mutation_scale=7,
                                                        shrinkA=2, shrinkB=2))
    ax.plot([xe], [xe], "o", ms=3.5, mfc="white", mec=LOST, mew=0.9, zorder=3)
    ax.plot([xe], [ye], "o", ms=3.5, color=LOST, zorder=3)
    ax.text(lim - 0.05 * Q, -1.75 * Q, "beyond λ/4:\nwrong size\nand sign", ha="right", va="center", fontsize=7,
            color=LOST, linespacing=1.05)
    ax.text(-0.3 * Q, -1.9 * Q, "within ±λ/4:\nfollows 1:1", ha="center", va="center", fontsize=7, color=fs.INK,
            linespacing=1.05)
    ax.text(-lim + 0.1 * Q, lim - 0.1 * Q, f"λ/4 ≈ {Q:.0f} mm\n(C-band)", ha="left", va="top", fontsize=7,
            color=fs.MUTED, linespacing=1.1)
    ticks = [-2 * Q, -Q, 0, Q, 2 * Q]
    labels = ["−λ/2", "−λ/4", "0", "λ/4", "λ/2"]
    ax.set_xticks(ticks, labels), ax.set_yticks(ticks, labels)
    ax.set_xlim(-lim, lim), ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.set_xlabel("true change between the dates\n(line of sight)", fontsize=7, linespacing=1.1)
    ax.set_ylabel("what the pair's phase gives", fontsize=7)
    ax.tick_params(labelsize=7)


def main() -> None:
    with fs.style():
        fig = plt.figure(figsize=fs.size("full", 150))
        a = fig.add_axes([0.0, 0.515, 1.0, 0.48])
        b = fig.add_axes([0.075, 0.27, 0.6, 0.2])
        lv = fig.add_axes([0.075, 0.045, 0.6, 0.215])
        c = fig.add_axes([0.775, 0.1, 0.21, 0.34])
        panel_section(a)
        panel_pairs(b, lv)
        panel_wrap(c)
        fig.text(0.005, 0.985, "(a)", fontsize=8, fontweight="bold", va="top")
        fig.text(0.005, 0.475, "(b)", fontsize=8, fontweight="bold", va="top")
        fig.text(0.70, 0.475, "(c)", fontsize=8, fontweight="bold", va="top")
        rec = fs.export(fig, OUT / "F3_schematic", provenance={
            "inputs": [], "kind": "schematic, no data", "incidence_drawn_deg": INC,
            "quarter_wavelength_mm": round(Q, 3), "wavelength_source": "insar_wetlands.cube.core.WAVELENGTH_MM",
            "wrap": "(c) line-of-sight change wrapped to ±λ/4 (two-way phase 4πd/λ wrapped to ±π)",
            "level_sketch": "(b) chain = truth increments + N(0, σ) per pair, σ larger in winter, seed 5; band ±2σ√n",
            "brief": "08_deliverables/article_plan/figures/briefs/F3.md"})
        plt.close(fig)
    print("F3 →", OUT, rec["width_mm"], "×", rec["height_mm"], "mm")


if __name__ == "__main__":
    main()

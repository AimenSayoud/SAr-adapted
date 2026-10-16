"""F4, final version — which pairs hold the mat's motion: the pair phase at P6 against the laser, by pair length, track
and period, in the groups of Table 2.

(a) One scatter per track and group: consecutive six-day pairs (2021: the laser starts 27 May 2021), next-but-one twelve-day pairs of the same
years, consecutive twelve-day pairs (2022–2024) and next-but-one twenty-four-day pairs of those years; the 1:1 line and,
shaded, laser changes beyond a quarter of the C-band wavelength, which one pair cannot follow unambiguously. (b) Every
group's correlation with its bootstrap 95 % interval (Table 2, ``article_plan/tables/revisit_correlations.csv``). Values
and intervals are read from the tables; no number is typed. Writes only to the hub
``08_deliverables/article_plan/figures/final/``.

    make run E=article_fig_f4
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

from insar_wetlands import figstyle as fs  # noqa: E402
from insar_wetlands.cube.core import WAVELENGTH_MM  # noqa: E402
from insar_wetlands.paths import hub_root  # noqa: E402

HUB = hub_root()
DLV = HUB / "08_deliverables"
OUT = DLV / "article_plan" / "figures" / "final"
PAIRS = DLV / "insar_full_test" / "reread_rtc" / "field_p6_short_pairs_x065" / "p6_pairs_2020_2024.csv"
TABLE2 = DLV / "article_plan" / "tables" / "revisit_correlations.csv"
GROUPS = (("2021", 6, "consecutive"), ("2021", 12, "next-but-one"),
          ("2022–2024", 12, "consecutive"), ("2022–2024", 24, "next-but-one"))
TRACKS = ("ascending", "descending")


def main() -> None:
    d = pd.read_csv(PAIRS).dropna(subset=["s1", "laser_los"])
    d["period"] = np.where(d.year <= 2021, "2021", "2022–2024")   # as Table 2; the laser starts 27 May 2021
    t2 = pd.read_csv(TABLE2).set_index(["track", "period", "revisit_days"])
    qw = WAVELENGTH_MM / 4
    shown = d[d.dt.isin({dt for _, dt, _ in GROUPS})]
    lim = float(np.nanmax(np.abs(np.r_[shown.s1, shown.laser_los]))) * 1.06     # every point inside the axes
    with fs.style():
        fig = plt.figure(figsize=fs.size("full", 96))
        outer = fig.add_gridspec(1, 2, width_ratios=[4.0, 1.25], wspace=0.42, left=0.075, right=0.985, bottom=0.2,
                                 top=0.9)
        gs = outer[0].subgridspec(2, 4, wspace=0.12, hspace=0.14)
        axs = np.array([[fig.add_subplot(gs[r, c]) for c in range(4)] for r in range(2)])
        for r, tr in enumerate(TRACKS):
            col = fs.TRACK[tr]
            for c, (per, dt, kind) in enumerate(GROUPS):
                ax = axs[r, c]
                g = d[(d.track == tr) & (d.period == per) & (d.dt == dt)]
                for s in (-1, 1):
                    ax.axvspan(s * qw, s * lim, color="#eef1f4", lw=0, zorder=0)
                ax.plot([-lim, lim], [-lim, lim], color=fs.MUTED, lw=0.6, ls=(0, (3, 2)), zorder=1)
                ax.axhline(0, color=fs.GRID, lw=0.5, zorder=1), ax.axvline(0, color=fs.GRID, lw=0.5, zorder=1)
                filled = kind == "consecutive"
                ax.scatter(g.laser_los, g.s1, s=10, color=col if filled else "white", edgecolors=col, lw=0.7,
                           zorder=3)
                row = t2.loc[(tr, per, dt)]
                ax.text(0.97, 0.03, f"r = {row.r:.2f}\n[{row.r_lo95:.2f}, {row.r_hi95:.2f}]\nn = {int(row.n)}",
                        transform=ax.transAxes, ha="right", va="bottom", fontsize=6.5, linespacing=1.15, zorder=5,
                        bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.2))
                ax.set_xlim(-lim, lim), ax.set_ylim(-lim, lim)
                ax.set_aspect("equal")
                if r == 0:
                    ax.set_title(f"{dt} days, {kind}\n{per}", fontsize=7, linespacing=1.15, pad=3)
                if c == 0:
                    ax.set_ylabel(f"{tr}\npair phase (mm, LOS)")
                else:
                    ax.tick_params(labelleft=False)
                if r == 1:
                    ax.set_xlabel("laser (mm, LOS)")
                else:
                    ax.tick_params(labelbottom=False)
        fs.panel_label(axs[0, 0], "a", x=-0.08, y=1.17)
        b = fig.add_subplot(outer[1])
        ylab = [f"{dt} d, {kind}\n{per}" for per, dt, kind in GROUPS]
        for k, tr in enumerate(TRACKS):
            for j, (per, dt, kind) in enumerate(GROUPS):
                row = t2.loc[(tr, per, dt)]
                y = j + (k - 0.5) * 0.28
                b.plot([row.r_lo95, row.r_hi95], [y, y], color=fs.TRACK[tr], lw=1.1, solid_capstyle="butt")
                b.plot(row.r, y, "o", ms=4.2, color=fs.TRACK[tr],
                       mfc=fs.TRACK[tr] if kind == "consecutive" else "white", mew=0.9)
        b.axvline(0, color=fs.MUTED, lw=0.6)
        b.set_yticks(range(len(GROUPS)), ylab, fontsize=7)
        b.set_ylim(len(GROUPS) - 0.5, -0.5)
        b.set_xlim(-0.65, 1.0)
        b.set_xlabel("correlation with the laser\n(bootstrap 95 % interval)")
        b.spines["left"].set_visible(False), b.tick_params(axis="y", length=0)
        fs.panel_label(b, "b", x=-0.02, y=1.02)
        fig.legend(handles=[Line2D([], [], ls="", marker="o", ms=4.5, color=fs.TRACK["ascending"], label="ascending"),
                            Line2D([], [], ls="", marker="o", ms=4.5, color=fs.TRACK["descending"], label="descending"),
                            Line2D([], [], ls="", marker="o", ms=4.5, mfc=fs.INK, mec=fs.INK, label="consecutive pair"),
                            Line2D([], [], ls="", marker="o", ms=4.5, mfc="white", mec=fs.INK,
                                   label="next-but-one pair"),
                            Line2D([], [], color=fs.MUTED, lw=0.6, ls=(0, (3, 2)), label="1:1"),
                            matplotlib.patches.Patch(fc="#eef1f4", label="laser beyond a quarter wavelength")],
                   loc="lower center", ncol=6, bbox_to_anchor=(0.5, 0.0), handlelength=1.4, columnspacing=1.2)
        rec = fs.export(fig, OUT / "F4_pairs_vs_laser", provenance={
            "inputs": [str(PAIRS.relative_to(HUB)), str(TABLE2.relative_to(HUB)), "insar_wetlands.cube.core.WAVELENGTH_MM"],
            "groups": [f"{tr} {per} {dt} d {kind}: n={int(t2.loc[(tr, per, dt)].n)}" for tr in TRACKS
                       for per, dt, kind in GROUPS],
            "quarter_wavelength_mm": qw})
        plt.close(fig)
    print("F4 →", OUT, rec["width_mm"], "×", rec["height_mm"], "mm")


if __name__ == "__main__":
    main()

"""F4, final version — which pairs hold the mat's motion: the pair phase at P6 against the laser, by pair length, track
and period, in the groups of Table 2. Since X-126 (2026-10-17): the hard-target reference, the corrected cube and the
common scorer's pairs; this script also writes Table 2 from X-126 (adding 2025's pairs, not drawn here).

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
X126 = DLV / "insar_full_test" / "x065_hard_targets"                 # X-126: hard targets, corrected cube, one scorer
PAIRS = X126 / "pairs.csv"
TABLE2 = DLV / "article_plan" / "tables" / "revisit_correlations.csv"
REF = "hard_targets"
TABLE2_ORDER = (("2021", 6), ("2021", 12), ("2022–2024", 12), ("2022–2024", 24), ("2025", 6), ("2025", 12))
GROUPS = (("2021", 6, "consecutive"), ("2021", 12, "next-but-one"),
          ("2022–2024", 12, "consecutive"), ("2022–2024", 24, "next-but-one"))
TRACKS = ("ascending", "descending")


def table2() -> pd.DataFrame:
    """Table 2 from X-126's by-revisit statistics, hard-target reference, in the article's row order."""
    t = pd.read_csv(X126 / "by_revisit.csv")
    t = t[t.reference == REF].rename(columns={"pairs": "n", "pearson_increments": "r"})
    t["order"] = [TABLE2_ORDER.index((p, d)) if (p, d) in TABLE2_ORDER else 99 for p, d in zip(t.period, t.revisit_days)]
    t = t[t.order < 99].sort_values(["track", "order"])
    # A property of the truth, not a test: the share of each group's pairs whose laser change exceeds a quarter
    # wavelength along the line of sight, which a single pair cannot follow unambiguously.
    pr = pd.read_csv(PAIRS).dropna(subset=["d_laser_los_mm", f"s1_{REF}"])
    beyond = (pr.assign(b=pr.d_laser_los_mm.abs() > WAVELENGTH_MM / 4)
              .groupby(["track", "group", "revisit_days"]).b.mean().rename("laser_beyond_quarter_wave"))
    t = t.merge(beyond, left_on=["track", "period", "revisit_days"], right_index=True, how="left")
    t = t[["track", "period", "revisit_days", "pairs_kind", "n", "r", "r_lo95", "r_hi95", "laser_beyond_quarter_wave"]]
    t.to_csv(TABLE2, index=False)
    return t


def main() -> None:
    d = pd.read_csv(PAIRS).drop(columns=["period"]).rename(columns={f"s1_{REF}": "s1", "d_laser_los_mm": "laser_los", "group": "period",
                                           "revisit_days": "dt"}).dropna(subset=["s1", "laser_los"])
    d["period"] = d.period.astype(str)
    t2 = table2().set_index(["track", "period", "revisit_days"])
    qw = WAVELENGTH_MM / 4
    shown = d[d.dt.isin({dt for _, dt, _ in GROUPS})]
    lim = float(np.nanmax(np.abs(np.r_[shown.s1, shown.laser_los]))) * 1.06     # every point inside the axes
    with fs.style():
        fig = plt.figure(figsize=fs.size("full", 96))
        outer = fig.add_gridspec(1, 2, width_ratios=[4.0, 1.25], wspace=0.42, left=0.075, right=0.985, bottom=0.2,
                                 top=0.86)
        gs = outer[0].subgridspec(2, 4, wspace=0.12, hspace=0.3)
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
                stats = f"r = {row.r:.2f} [{row.r_lo95:.2f}, {row.r_hi95:.2f}], n = {int(row.n)}"   # above the
                # axes, so no point is hidden behind it and it hides none
                ax.set_xlim(-lim, lim), ax.set_ylim(-lim, lim)
                ax.set_aspect("equal")
                head = f"{dt} days, {kind}\n{per}\n" if r == 0 else ""
                ax.set_title(head + stats, fontsize=6.3, linespacing=1.2, pad=3)
                if c == 0:
                    ax.set_ylabel(f"{tr}\npair phase (mm, LOS)")
                else:
                    ax.tick_params(labelleft=False)
                if r == 1:
                    ax.set_xlabel("laser (mm, LOS)")
                else:
                    ax.tick_params(labelbottom=False)
        fs.panel_label(axs[0, 0], "a", x=-0.08, y=1.25)
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

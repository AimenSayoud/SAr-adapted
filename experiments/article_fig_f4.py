"""F4, final version — which pairs hold the mat's motion, and why the others lose it.

Round 2 (C-070, brief ``article_plan/figures/briefs/F4.md``): (a) the pair phase at P6 against the laser's change, one
panel per period (2021, 2022–2024, 2025) and kind of pair (consecutive; skipping one acquisition), both tracks in each;
(b) the correlation of every group with its bootstrap 95 % interval (Table 2), and for the skip-one groups the range of
correlation ambiguity alone would leave (X-127); (c) the share of each group's pairs whose laser change exceeds a quarter
wavelength (Table 2). Since X-126 (2026-10-17): the hard-target reference, the corrected cube and the
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
PERIODS = (("2021", 6, 12), ("2022–2024", 12, 24), ("2025", 6, 12))    # (period, consecutive days, skip-one days)
X127 = DLV / "insar_full_test" / "skip_pairs" / "by_group.csv"      # X-127: r under ambiguity alone, skip-one groups
MARK = {"ascending": "o", "descending": "^"}
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


PERIOD_HEAD = {"2021": "2021 — Sentinel-1A + 1B", "2022–2024": "2022–2024 — Sentinel-1A alone",
               "2025": "2025 — Sentinel-1A + 1C"}


def main() -> None:
    d = pd.read_csv(PAIRS).drop(columns=["period"]).rename(columns={f"s1_{REF}": "s1", "d_laser_los_mm": "laser_los",
                                                                  "group": "period", "revisit_days": "dt"})
    d = d.dropna(subset=["s1", "laser_los"])
    d["period"] = d.period.astype(str)
    t2 = table2()
    t2["period"] = t2.period.astype(str)
    br = pd.read_csv(X126 / "by_revisit.csv")
    br = br[br.reference == REF].assign(period=lambda x: x.period.astype(str))
    sim = pd.read_csv(X127)
    sim["period"] = sim.period.astype(str)
    qw = WAVELENGTH_MM / 4
    shown = d[[(p, dt, b) in {(p0, c, 0) for p0, c, _ in PERIODS} | {(p0, k, 1) for p0, _, k in PERIODS}
               for p, dt, b in zip(d.period, d.dt, d.between)]]
    lim = float(np.nanmax(np.abs(np.r_[shown.s1, shown.laser_los]))) * 1.06      # every point inside the axes
    with fs.style():
        fig = plt.figure(figsize=fs.size("full", 132))
        outer = fig.add_gridspec(1, 2, width_ratios=[1.0, 0.95], wspace=0.30, left=0.075, right=0.905, bottom=0.135,
                                 top=0.93)
        ga = outer[0].subgridspec(3, 2, wspace=0.08, hspace=0.12)
        axs = np.array([[fig.add_subplot(ga[r, c]) for c in range(2)] for r in range(3)])
        for r, (per, cons, skip) in enumerate(PERIODS):
            for c, (days, between) in enumerate(((cons, 0), (skip, 1))):
                ax = axs[r, c]
                for sgn in (-1, 1):
                    ax.axvspan(sgn * qw, sgn * lim, color="#eef1f4", lw=0, zorder=0)
                ax.plot([-lim, lim], [-lim, lim], color=fs.MUTED, lw=0.6, ls=(0, (3, 2)), zorder=1)
                ax.axhline(0, color=fs.GRID, lw=0.5, zorder=1), ax.axvline(0, color=fs.GRID, lw=0.5, zorder=1)
                slopes = []
                for tr in TRACKS:
                    g = d[(d.track == tr) & (d.period == per) & (d.dt == days) & (d.between == between)]
                    ax.scatter(g.laser_los, g.s1, s=9, marker=MARK[tr], color=fs.TRACK[tr] if between == 0 else "white",
                               edgecolors=fs.TRACK[tr], lw=0.6, zorder=3)
                    k = br[(br.track == tr) & (br.period == per) & (br.revisit_days == days)].slope_through_origin
                    if len(k):
                        k = float(k.iloc[0])
                        ax.plot([-lim, lim], [-k * lim, k * lim], color=fs.TRACK[tr], lw=0.9, zorder=2)
                        slopes.append((tr, k))
                for j, (tr, k) in enumerate(slopes):                # the slope through the origin, per track
                    ax.text(0.96, 0.05 + 0.1 * (len(slopes) - 1 - j), f"slope {k:.2f}".replace("-", "−"), transform=ax.transAxes,
                            ha="right", va="bottom", fontsize=7, color=fs.TRACK[tr],
                            bbox=dict(fc="white", ec="none", pad=0.6, alpha=0.85), zorder=4)
                ax.set_xlim(-lim, lim), ax.set_ylim(-lim, lim)
                ax.set_aspect("equal")
                ax.text(0.04, 0.96, f"{per}, {days} days", transform=ax.transAxes, ha="left", va="top", fontsize=7,
                        bbox=dict(fc="white", ec="none", pad=1.0, alpha=0.85), zorder=4)
                if r == 0:
                    ax.set_title("consecutive pairs" if c == 0 else "pairs that skip one acquisition", fontsize=8)
                if c == 0:
                    ax.set_ylabel("pair phase (mm, LOS)")
                else:
                    ax.tick_params(labelleft=False)
                if r == 2:
                    ax.set_xlabel("laser change (mm, LOS)")
                else:
                    ax.tick_params(labelbottom=False)
        fs.panel_label(axs[0, 0], "a", x=-0.2, y=1.08)
        # (b): one forest plot, grouped by period; the share of laser steps beyond λ/4 as a column of numbers
        b = fig.add_subplot(outer[1])
        y, ticks, rows = 0.0, [], []
        for per, cons, skip in PERIODS:
            b.text(-0.73, y, PERIOD_HEAD[per], fontsize=7.5, fontweight="bold", va="center", ha="left", zorder=5,
                   bbox=dict(fc="white", ec="none", pad=1.5))
            y += 0.75
            for days, kind in ((cons, "consecutive"), (skip, "skips one")):
                ticks.append((y + 0.15, f"{kind}, {days} d"))
                for k, tr in enumerate(TRACKS):
                    rows.append((per, days, kind, tr, y + k * 0.3))
                y += 0.95
            y += 0.25
        labelled = False
        for per, days, kind, tr, yy in rows:
            row = t2[(t2.track == tr) & (t2.period == per) & (t2.revisit_days == days)].iloc[0]
            if kind == "skips one":
                s = sim[(sim.track == tr) & (sim.period == per)].iloc[0]
                b.add_patch(plt.Rectangle((s.r_sim_p2_5, yy - 0.12), s.r_sim_p97_5 - s.r_sim_p2_5, 0.24,
                                          color="#d9dde2", lw=0, zorder=0))
                if not labelled:
                    b.annotate("range ambiguity + indep.\n6 d noise would give (X-127)", xy=(s.r_sim_p2_5 + 0.02, yy - 0.12),
                               xytext=(-0.72, yy - 0.62), fontsize=7, color=fs.MUTED, ha="left", va="bottom",
                               arrowprops=dict(arrowstyle="-", color=fs.MUTED, lw=0.5))
                    labelled = True
            b.plot([row.r_lo95, row.r_hi95], [yy, yy], color=fs.TRACK[tr], lw=1.0, solid_capstyle="butt", zorder=2)
            b.plot(row.r, yy, MARK[tr], ms=4.2, color=fs.TRACK[tr], mfc=fs.TRACK[tr] if kind == "consecutive" else "white",
                   mew=0.9, zorder=3)
            b.text(1.06, yy, f"{row.laser_beyond_quarter_wave * 100:.0f} %", fontsize=7, color=fs.TRACK[tr],
                   va="center", ha="left", clip_on=False)
        b.text(1.06, -0.55, "laser\nbeyond λ/4", fontsize=7, va="bottom", ha="left", clip_on=False)
        b.set_ylim(y - 0.2, -0.6)
        b.set_yticks([t for t, _ in ticks], [lab for _, lab in ticks], fontsize=7)
        b.spines["left"].set_visible(False), b.tick_params(axis="y", length=0)
        b.axvline(0, color=fs.MUTED, lw=0.6)
        b.set_xlim(-0.75, 1.0)
        b.set_xticks([-0.5, 0, 0.5, 1.0])
        b.set_xlabel("correlation with the laser (bootstrap 95 % interval)")
        fs.panel_label(b, "b", x=-0.02, y=1.0)
        fig.legend(handles=[Line2D([], [], ls="", marker="o", ms=4.5, color=fs.TRACK["ascending"], label="ascending"),
                            Line2D([], [], ls="", marker="^", ms=4.5, color=fs.TRACK["descending"], label="descending"),
                            Line2D([], [], ls="", marker="o", ms=4.5, mfc=fs.INK, mec=fs.INK, label="consecutive pair"),
                            Line2D([], [], ls="", marker="o", ms=4.5, mfc="white", mec=fs.INK, label="skips one acquisition"),
                            Line2D([], [], color=fs.INK, lw=0.9, label="fit through the origin"),
                            Line2D([], [], color=fs.MUTED, lw=0.6, ls=(0, (3, 2)), label="1:1"),
                            matplotlib.patches.Patch(fc="#eef1f4", label="laser beyond λ/4")],
                   loc="lower center", ncol=7, bbox_to_anchor=(0.5, 0.0), handlelength=1.3, columnspacing=1.0,
                   fontsize=7)
        rec = fs.export(fig, OUT / "F4_pairs_vs_laser", provenance={
            "inputs": [str(PAIRS.relative_to(HUB)), str(TABLE2.relative_to(HUB)), str(X127.relative_to(HUB)),
                       str((X126 / "by_revisit.csv").relative_to(HUB)), "insar_wetlands.cube.core.WAVELENGTH_MM"],
            "brief": "08_deliverables/article_plan/figures/briefs/F4.md",
            "groups": [f"{tr} {per} {days} d {kind}" for per, days, kind, tr, _ in rows],
            "fit": "slope through the origin, insar_wetlands.scoring.increment_stats (X-126 by_revisit.csv)",
            "quarter_wavelength_mm": qw})
        plt.close(fig)
    print("F4 →", OUT, rec["width_mm"], "×", rec["height_mm"], "mm")


if __name__ == "__main__":
    main()

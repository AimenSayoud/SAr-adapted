"""Graphical abstract — the article's three answers in three pictures, read from the tables the article cites.

(1) The mat floats with its water, in part: the laser's mat level and the water level at P6 over the summer of 2021, as
changes from the first hour, on one axis in centimetres (no second axis; the mat's flatter line is the buoyancy of §4.2).
(2) Only consecutive six-day pairs follow it: pair phase against the laser's change, 2021, consecutive six-day pairs
against pairs that skip one acquisition, both tracks; the correlations are read from Table 2. (3) Across a winter the
chain loses the level: the plain chain of consecutive pairs (ascending, hard-target reference, converted to vertical by
the cube's up-component at P6) against the laser from May 2023 to April 2024, both set to zero at the chain's first
date; the window ends before the laser's April 2024 step (§2.2). No number is typed; field values stay in the hub.
Writes only to the hub ``08_deliverables/article_plan/figures/final/GA_graphical_abstract.*``.

    make run E=article_graphical_abstract
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from insar_wetlands import figstyle as fs  # noqa: E402
from insar_wetlands.paths import cube_root, hub_root  # noqa: E402

HUB, CUBE = hub_root(), cube_root()
OUT = HUB / "08_deliverables" / "article_plan" / "figures" / "final"
PAIRS = HUB / "08_deliverables" / "insar_full_test" / "x065_hard_targets" / "pairs.csv"
TABLE2 = HUB / "08_deliverables" / "article_plan" / "tables" / "revisit_correlations.csv"
SUMMER = ("2021-06-01", "2021-11-15")              # datum A, snow-free
WINTER = ("2023-05-01", "2024-04-06")              # an unbroken ascending chain across winter 2023/24, before the step
WATER = "#0072B2"


def truthy(s: pd.Series) -> pd.Series:
    return s.astype(str).str.lower().isin(["true", "1", "1.0"])


def main() -> None:
    h = pd.read_csv(CUBE / "gold" / "truth_laser_hourly.csv.gz",
                    usecols=["time_utc", "laser_surface_cm", "laser_ok", "cr_raw"])
    h["time_utc"] = pd.to_datetime(h.time_utc, utc=True)
    h = h.set_index("time_utc").sort_index()
    ok = truthy(h.laser_ok)
    p = pd.read_csv(PAIRS, dtype={"group": str})
    for c in ("t1", "t2"):
        p[c] = pd.to_datetime(p[c], utc=True)
    t2 = pd.read_csv(TABLE2, dtype={"period": str})
    los_u = pd.read_csv(CUBE / "gold" / "plots_pair.csv.gz", usecols=["plot", "track", "grid", "los_u"])
    los_u = float(los_u[(los_u["plot"] == "P6") & (los_u.grid == "40m") & (los_u.track == "ascending")].los_u.median())
    if not np.isfinite(los_u):
        raise SystemExit("graphical abstract: no up-component at P6")

    with fs.style():
        fig = plt.figure(figsize=fs.size("full", 70))
        gs = fig.add_gridspec(1, 3, width_ratios=[1.1, 0.8, 1.1], left=0.06, right=0.99, bottom=0.3, top=0.81,
                              wspace=0.36)
        a, b, c = (fig.add_subplot(gs[i]) for i in range(3))

        # (1) the mat floats with its water, in part
        s = h.loc[SUMMER[0]:SUMMER[1]]
        lv = s.laser_surface_cm.where(ok.loc[s.index])
        wt = s.cr_raw
        a.plot(wt.index, wt - wt.dropna().iloc[0], color=WATER, lw=0.8, label="water level")
        a.plot(lv.index, lv - lv.dropna().iloc[0], color=fs.INK, lw=1.0, label="mat surface (laser)")
        a.axhline(0, color=fs.GRID, lw=0.5, zorder=0)
        a.set_ylabel("change since June (cm)")
        a.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%b"))
        a.xaxis.set_major_locator(matplotlib.dates.MonthLocator(bymonth=(6, 8, 10)))
        a.legend(loc="upper right", frameon=False, handlelength=1.5)
        a.set_title("The mat floats with its water,\nrising and falling by part of it", fontsize=8, loc="left")

        # (2) only consecutive six-day pairs follow it
        q = p[(p.group == "2021") & p.d_laser_los_mm.notna() & p.s1_hard_targets.notna()]
        cons, skip = q[(q.revisit_days == 6) & (q.between == 0)], q[(q.revisit_days == 12) & (q.between == 1)]
        lim = float(np.nanmax(np.abs(np.r_[q.d_laser_los_mm, q.s1_hard_targets]))) * 1.05
        b.plot([-lim, lim], [-lim, lim], color=fs.MUTED, lw=0.6, ls=(0, (3, 2)), zorder=1)
        b.scatter(skip.d_laser_los_mm, skip.s1_hard_targets, s=9, fc="white", ec=fs.MUTED, lw=0.6, zorder=2)
        b.scatter(cons.d_laser_los_mm, cons.s1_hard_targets, s=9, fc=fs.CLASS["mat"], ec="white", lw=0.3, zorder=3)
        r = t2[t2.period == "2021"].set_index(["track", "revisit_days"]).r
        b.legend(handles=[plt.Line2D([], [], ls="", marker="o", ms=4, mfc=fs.CLASS["mat"], mec="white",
                                     label=f"6-day: {r['ascending', 6]:.2f} / {r['descending', 6]:.2f}"),
                          plt.Line2D([], [], ls="", marker="o", ms=4, mfc="white", mec=fs.MUTED,
                                     label=f"skip one: {r['ascending', 12]:.2f} / {r['descending', 12]:.2f}")],
                 loc="upper center", bbox_to_anchor=(0.5, -0.3), ncol=2, frameon=False, handletextpad=0.2,
                 columnspacing=1.0, fontsize=7, title="correlation with the laser, ascending / descending",
                 title_fontsize=7)
        b.set_xlim(-lim, lim), b.set_ylim(-lim, lim), b.set_aspect("equal")
        b.set_xlabel("laser change (mm)")
        b.set_ylabel("radar pair (mm)")
        b.set_title("Only six-day pairs\nfollow it", fontsize=8, loc="left")

        # (3) across a winter the chain loses the level
        w0, w1 = (pd.Timestamp(t, tz="UTC") for t in WINTER)
        ch = p[(p.track == "ascending") & (p.between == 0) & (p.t1 >= w0) & (p.t2 <= w1)].sort_values("t1")
        if (ch.t1.to_numpy()[1:] != ch.t2.to_numpy()[:-1]).any():
            raise SystemExit("graphical abstract: the winter chain is broken")
        level = np.r_[0.0, np.cumsum(ch.s1_hard_targets.to_numpy() / los_u)]
        dates = pd.DatetimeIndex([ch.t1.iloc[0], *ch.t2])
        w = h.loc[WINTER[0]:WINTER[1]]
        lw_ = w.laser_surface_cm.where(ok.loc[w.index]) * 10
        start = lw_.loc[dates[0] - pd.Timedelta(hours=6):dates[0] + pd.Timedelta(hours=6)]
        lw_ = lw_ - float(start.mean())
        c.axvspan(pd.Timestamp("2023-11-15", tz="UTC"), pd.Timestamp("2024-03-15", tz="UTC"), color="#eef1f4", lw=0,
                  zorder=0)
        c.plot(lw_.index, lw_, color=fs.INK, lw=1.0, label="laser")
        c.plot(dates, level, color=fs.MUTED, lw=1.0, ls=(0, (3, 2)), marker="o", ms=2.2, label="chain of pairs")
        c.text(pd.Timestamp("2024-01-15", tz="UTC"), 0.04, "winter", transform=c.get_xaxis_transform(), ha="center",
               fontsize=7, color=fs.MUTED)
        c.axhline(0, color=fs.GRID, lw=0.5, zorder=0)
        c.set_ylabel("level since May (mm)")
        c.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%b"))
        c.xaxis.set_major_locator(matplotlib.dates.MonthLocator(bymonth=(5, 8, 11, 2)))
        c.legend(loc="upper left", frameon=False, handlelength=1.8)
        c.set_title("Across a winter the radar\nloses the mat's level", fontsize=8, loc="left")

        fig.text(0.5, 0.95, "A floating fen mat against Sentinel-1: a laser scores which radar pairs keep its motion",
                 ha="center", va="center", fontsize=8.5, fontweight="bold")
        rec = fs.export(fig, OUT / "GA_graphical_abstract", provenance={
            "inputs": ["cube gold/truth_laser_hourly.csv.gz (laser, laser_ok, CR water level)",
                       "hub 08_deliverables/insar_full_test/x065_hard_targets/pairs.csv (X-126 pairs)",
                       "hub 08_deliverables/article_plan/tables/revisit_correlations.csv (Table 2)",
                       "cube gold/plots_pair.csv.gz (los_u at P6, ascending)"],
            "summer": SUMMER, "winter": WINTER, "winter_pairs": int(len(ch)), "los_u_ascending": los_u})
        plt.close(fig)
    print("graphical abstract →", OUT, rec["width_mm"], "×", rec["height_mm"], "mm")


if __name__ == "__main__":
    main()

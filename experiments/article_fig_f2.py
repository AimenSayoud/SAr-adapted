"""F2, final version — how the mat moves: its level and its water at P6, the buoyancy that ties them, and the nine
radar cells around the laser that do not move alike.

(a) The laser's mat-surface level at P6 through the design years, by datum segment (the 2022 re-levelling marked;
snow periods and the delivery's filled stretch drawn apart, never as data). (b) The water level at P6 on the same time
axis. (c) Six-day changes of the surface against those of the water table, with the pooled buoyancy of the year-blocked
fit (X-074 E1.2) and, for scale, what a raft floating free would do. (d) Each of the nine cells' slope of chained phase
on the laser against the cell's depth into the mat, with 95 % intervals (X-079 on corrected inputs). Every value drawn
is read from the cube or the results tables; no number is typed. Writes only to the hub
``08_deliverables/article_plan/figures/final/``.

    make run E=article_fig_f2
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
import yaml  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

from insar_wetlands import figstyle as fs  # noqa: E402
from insar_wetlands.cube.core import WAVELENGTH_MM  # noqa: E402
from insar_wetlands.paths import cube_root, hub_root  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import field_p6_cells_x079 as X079  # noqa: E402  (the cells' depth, as X-079 computed it)
import fusion_v3_wp1_x074 as WP1  # noqa: E402  (the daily series and six-day changes of E1.2)

HUB, CUBE = hub_root(), cube_root()
DLV = HUB / "08_deliverables"
OUT = DLV / "article_plan" / "figures" / "final"
CELLS = DLV / "insar_full_test" / "reread_corrected" / "field_p6_cells_x079"
# §4.4: the gradient on the ascending track in both periods; descending 2020–2021 as the case that does not replicate.
# Descending 2022–2024 is left out: no cell's phase follows the laser there (cell r −0.05 to 0.27), so no slope exists.
SHOWN = (("ascending", "2020–2021", "o", True), ("ascending", "2022–2024", "o", False),
         ("descending", "2020–2021", "s", False))
WATER = "#2C7FB8"
SHOW_AS = {"2020–2021": "2021"}    # X-079's six-day period; the laser record starts 27 May 2021


def panel_series(a, b, t, breaks) -> None:
    ok, snow, filled = (t[c].astype(bool) for c in ("laser_ok", "laser_snow_72h", "laser_filled"))
    a.plot(t.time_utc, t.laser_surface_raw_cm.where(snow), color="#b8bfc6", lw=0.6, zorder=1)
    a.plot(t.time_utc, t.laser_surface_cm.where(ok), color=fs.INK, lw=0.7, zorder=2)
    a.plot(t.time_utc, t.laser_surface_raw_cm.where(filled), color=fs.CLASS["laser"], lw=0.9, ls=(0, (2, 1.5)),
           zorder=3)
    b.plot(t.time_utc, t.cr_raw, color=WATER, lw=0.7)
    for x0, x1 in breaks:
        for ax in (a, b):
            ax.axvspan(pd.Timestamp(x0, tz="UTC"), pd.Timestamp(x1, tz="UTC"), fc="none", ec=fs.MUTED, hatch="////",
                       lw=0, zorder=0)
    lo, hi = a.get_ylim()
    a.text(pd.Timestamp(breaks[0][1], tz="UTC"), hi, "  ← re-levelled: levels either side are not comparable",
           ha="left", va="top", fontsize=7, color=fs.MUTED)
    first, last = t.time_utc.min(), t.time_utc.max()
    b0, b1 = pd.Timestamp(breaks[0][0], tz="UTC"), pd.Timestamp(breaks[0][1], tz="UTC")
    for x, lab in (((first + (b0 - first) / 2), "datum A"), ((b1 + (last - b1) / 2), "datum B")):
        a.text(x, lo, lab, ha="center", va="bottom", fontsize=7, color=fs.MUTED)
    a.set_ylabel("Mat surface at P6\n(cm, laser, + up)")
    b.set_ylabel("Water level at P6\n(cm, logger datum)")
    b.set_xlim(first, last)
    a.legend(handles=[Line2D([], [], color=fs.INK, lw=0.9, label="laser, quality-checked"),
                      Line2D([], [], color="#b8bfc6", lw=0.9, label="snow period (not used)"),
                      Line2D([], [], color=fs.CLASS["laser"], lw=0.9, ls=(0, (2, 1.5)),
                             label="filled in the delivery (not data)"),
                      Patch(fc="none", ec=fs.MUTED, hatch="////", lw=0, label="re-levelling gap")],
             loc="upper right", ncol=4, bbox_to_anchor=(1.0, 1.2), handlelength=1.8, columnspacing=1.2)


def panel_buoyancy(ax, inc_asc: float) -> dict:
    d = WP1.daily()
    c = WP1.changes(d)
    ok = c.d_laser_cm.notna() & c.d_wtd_cm.notna()
    dw, ds = c.d_wtd_cm[ok].to_numpy(), 10 * c.d_laser_cm[ok].to_numpy()          # cm, mm
    g = 10 * float(pd.read_csv(DLV / "fusion_v3" / "wp1" / "e12_coefficients.csv")["g·ΔW"].iloc[0])   # mm per cm
    m = pd.read_csv(DLV / "fusion_v3" / "wp1" / "e12_change_models.csv").set_index("model (6-day laser change)")
    r_yb = float(m.loc["g·ΔW", "r"])                                                 # year-blocked (X-074 E1.2)
    qv = WAVELENGTH_MM / 4 / np.cos(np.radians(inc_asc))     # a quarter wavelength of LOS, as vertical, ascending
    ax.scatter(dw, ds, s=4, color=fs.INK, alpha=0.35, lw=0, rasterized=True, zorder=2)
    lim = np.nanmax(np.abs(dw)) * 1.05
    xx = np.array([-lim, lim])
    ax.plot(xx, g * xx, color=fs.CLASS["mat"], lw=1.4, zorder=3,
            label=f"buoyancy {g:.2f} mm/cm\n(r {r_yb:.2f}, years held out)")
    ax.plot(xx, 10 * xx, color=fs.MUTED, lw=0.9, ls=(0, (3, 2)), zorder=1, label="free raft: 10 mm/cm")
    ax.axhline(0, color=fs.GRID, lw=0.6, zorder=0), ax.axvline(0, color=fs.GRID, lw=0.6, zorder=0)
    ax.set_xlim(-lim, lim)
    y0, y1 = np.nanpercentile(ds, 0.5) * 1.3, np.nanpercentile(ds, 99.5) * 1.3
    ax.set_ylim(min(y0 * 1.18, -qv * 1.35), max(y1, qv * 1.25))
    for sgn in (-1, 1):
        ax.axhspan(sgn * qv, sgn * 200, color="#eef1f4", lw=0, zorder=0)
    beyond = float(np.mean(np.abs(ds) > qv))
    ax.text(0.03, 0.97, "shaded: beyond λ/4\n(ascending LOS, as vertical)", transform=ax.transAxes,
            ha="left", va="top", fontsize=7, color=fs.MUTED)
    ax.set_xlabel("Six-day change of the water table (cm)")
    ax.set_ylabel("Six-day change of the\nmat surface (mm)")
    ax.legend(loc="lower right", handlelength=1.6, fontsize=7, frameon=True, framealpha=0.92, edgecolor="none",
              borderaxespad=0.3)
    return {"g_mm_per_cm": g, "n": int(ok.sum()), "r_year_blocked": r_yb, "quarter_wave_vertical_mm": qv,
            "share_beyond_quarter_wave": beyond}


def panel_rain(top, bot) -> dict:
    """§4.3: superposed epochs after rain (X-055) — the water table's pulse above, the mat surface's step below."""
    e = pd.read_csv(DLV / "field_events_x055" / "epoch_response.csv")
    out = {}
    for ax, series, col, lab in ((top, "WTD plot median", WATER, "water table, plots"),
                                 (bot, "laser surface P6", fs.INK, "mat surface, P6")):
        g = e[e.series == series].sort_values("lag_h")
        days = g.lag_h / 24
        ax.fill_between(days, 10 * g.q1, 10 * g.q3, color=col, alpha=0.18, lw=0)
        ax.plot(days, 10 * g["median"], color=col, lw=1.1)
        ax.axvline(0, color=fs.MUTED, lw=0.6, ls=(0, (3, 2)))
        ax.axhline(0, color=fs.GRID, lw=0.6)
        ax.set_ylabel("change (mm)")
        ax.text(0.22, 0.95, f"{lab} ({int(g.n_events.iloc[0])} events)", transform=ax.transAxes, ha="left", va="top",
                fontsize=7, color=col)
        out[series] = int(g.n_events.iloc[0])
    top.tick_params(labelbottom=False)
    bot.set_xlabel("days after rain onset (median, interquartile range)")
    bot.set_xlim(-2, 10)
    return out


def panel_cells(ax, st) -> dict:
    px = pd.read_csv(DLV / "field_first" / "plot_pixels.csv").set_index("plot")
    desc = X079.descriptors(st, int(px.loc["P6", "row"]), int(px.loc["P6", "col"]))
    cells = pd.read_csv(CELLS / "cells.csv")
    drv = pd.read_csv(CELLS / "drivers.csv")
    out = {}
    for k, (tr, per, mk, filled) in enumerate(SHOWN):
        cc = cells[(cells.track == tr) & (cells.period == per)].set_index("cell").join(desc.depth_mat_edge_m)
        row = drv[(drv.track == tr) & (drv.period == per) & (drv.slope == "slope_tls")
                  & (drv["against"] == "depth_mat_edge_m")].iloc[0]
        col = fs.TRACK[tr]
        x = cc.depth_mat_edge_m + (k - 1) * 1.6
        ax.errorbar(x, cc.slope_tls, yerr=[cc.slope_tls - cc.tls_lo95, cc.tls_hi95 - cc.slope_tls], fmt=mk, ms=3.6,
                    color=col, mfc=col if filled else "white", mec=col, elinewidth=0.6, capsize=0, zorder=3,
                    label=f"{tr[:3]}. {SHOW_AS.get(per, per)}: r {row.r:.2f}, p {row.p_perm:.3f}")
        out[f"{tr} {per}"] = {"r": float(row.r), "p_perm": float(row.p_perm), "n_cells": int(len(cc))}
    ax.axhline(1, color=fs.MUTED, lw=0.7, ls=(0, (3, 2)), zorder=1, label="slope 1: as the laser")
    ax.set_ylim(0, None)
    dep = desc.depth_mat_edge_m.sort_values()
    groups, cur = [], [dep.index[0]]
    for a_, b_ in zip(dep.index[:-1], dep.index[1:]):                  # cells within 6 m of depth share one label
        (cur.append(b_) if dep[b_] - dep[a_] < 6 else (groups.append(cur), cur := [b_]))
    groups.append(cur)
    for k, g in enumerate(groups):                                  # alternate heights so close labels do not touch
        ax.text(float(dep[g].mean()), 0.02 + 0.06 * (k % 2), "/".join(g), transform=ax.get_xaxis_transform(),
                ha="center", va="bottom", fontsize=7, color=fs.MUTED)
    ax.set_xlabel("depth of the cell into the mat\n(m, from its edge)")
    ax.set_ylabel("slope of the cell's chained phase\non the laser (TLS, 95 % CI)")
    top = float(np.nanmax(np.concatenate([cells[(cells.track == tr) & (cells.period == per)].tls_hi95
                                          for tr, per, _, _ in SHOWN])))
    ax.set_ylim(0, top * 1.22)                          # headroom for the legend above the tallest interval
    ax.legend(loc="upper left", handlelength=1.4, fontsize=7, ncol=1, frameon=True, framealpha=0.92,
              edgecolor="none", borderaxespad=0.15, labelspacing=0.22)
    return out


def main() -> None:
    t = pd.read_csv(CUBE / "gold" / "truth_laser_hourly.csv.gz", parse_dates=["time_utc"])
    breaks = yaml.safe_load((DLV / "insar_full_test" / "protocol.yaml").read_text())["datum_breaks"]
    st = xr.open_dataset(CUBE / "silver" / "static_40m.nc")
    with fs.style(), matplotlib.rc_context({"hatch.linewidth": 0.5}):
        fig = plt.figure(figsize=fs.size("full", 168))
        gs = fig.add_gridspec(3, 1, height_ratios=[1.0, 0.8, 1.3], hspace=0.30, left=0.085, right=0.985,
                              bottom=0.085, top=0.94)
        a = fig.add_subplot(gs[0])
        b = fig.add_subplot(gs[1], sharex=a)
        bottom = gs[2].subgridspec(1, 3, wspace=0.55, width_ratios=[1.0, 0.85, 1.0])
        c = fig.add_subplot(bottom[0])
        rd = bottom[1].subgridspec(2, 1, hspace=0.12)
        d1, d2 = fig.add_subplot(rd[0]), fig.add_subplot(rd[1])
        e = fig.add_subplot(bottom[2])
        panel_series(a, b, t, breaks)
        plt.setp(a.get_xticklabels(), visible=False)
        px = pd.read_csv(DLV / "field_first" / "plot_pixels.csv").set_index("plot")
        inc_asc = float(st.incidence_asc_deg.values[int(px.loc["P6", "row"]), int(px.loc["P6", "col"])])
        g = panel_buoyancy(c, inc_asc)
        rain = panel_rain(d1, d2)
        cl = panel_cells(e, st)
        for ax, l in ((a, "a"), (b, "b")):
            fs.panel_label(ax, l, x=-0.01, y=1.02)
        for ax, l in ((c, "c"), (d1, "d"), (e, "e")):
            fs.panel_label(ax, l, x=-0.03, y=1.04)
        fig.align_ylabels([a, b])
        rec = fs.export(fig, OUT / "F2_mat_motion", provenance={
            "inputs": ["cube gold/truth_laser_hourly.csv.gz (laser_surface_cm, laser_ok/snow/filled, cr_raw)",
                       "hub insar_full_test/protocol.yaml (datum_breaks)",
                       "cube silver site_hourly / plots (via fusion_v3_wp1_x074.daily, changes)",
                       "hub fusion_v3/wp1/e12_coefficients.csv (g·ΔW)",
                       "hub insar_full_test/reread_corrected/field_p6_cells_x079/cells.csv, drivers.csv",
                       "cube silver/static_40m.nc (zone → depth from the mat's edge, as X-079)"],
            "buoyancy": g, "cells": cl, "rain_events": rain,
            "rain": "hub field_events_x055/epoch_response.csv (X-055 superposed epochs)",
            "brief": "08_deliverables/article_plan/figures/briefs/F2.md"})
        plt.close(fig)
    print("F2 →", OUT, rec["width_mm"], "×", rec["height_mm"], "mm")


if __name__ == "__main__":
    main()

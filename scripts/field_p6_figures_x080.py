"""X-080 — figures of P6 and its nine 40 m cells, from verified data only.

Four figures, each with one job, every number on them computed here or read from a generated table:

1. **Site and cells.** Where the nine cells sit (zones on the 40 m grid, plot cells), what they contain
   (0.25 m orthophoto, X-060's boardwalk outline, the laser station) and how coherent each is (median over
   consecutive pairs, both tracks, six-day 2020–2021 and twelve-day 2022–2024; X-079 `cells.csv`).
2. **P6 through time.** Water table (P6, censored stretches out), laser surface (X-058 QC: snow, outliers,
   gap-filled hours out), and for every consecutive ascending pair the laser's LOS change beside the P6 cell's phase
   change (X-065 pairs), with the pair's coherence. Pair changes are drawn as separate marks, never joined.
3. **Phase vs laser at P6**, both tracks × both revisits: the 1:1 line, X-065's motion coefficient m ± 1.96 SD,
   the λ/4 limit, points coloured by coherence.
4. **The nine cells.** m per cell (diverging around 1) with its interval, and m against depth into the mat and
   against LiDAR woody cover with X-079's partial correlations.

Field data are read from the hub at run time; outputs go only to the hub. No coordinate or field value is
written in this file. Needs X-079's outputs (`field_p6_cells_x079.py`) first.

    PYTHONPATH=src python scripts/field_p6_figures_x080.py
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import rasterio  # noqa: E402
import rasterio.windows  # noqa: E402
import xarray as xr  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, TwoSlopeNorm  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402

from insar_wetlands import field  # noqa: E402

HUB = field.field_root().parents[1]
DLV = HUB / "08_deliverables"
SIL = HUB / "06_data" / "cube" / "silver"
X079 = DLV / "field_p6_cells_x079"
OUT = DLV / "field_p6_figures_x080"
QW = 55.465763 / 4
CELL_M = 40.0
GRID = (("NW", "N", "NE"), ("W", "P6", "E"), ("SW", "S", "SE"))
PERIODS = ("2020–2021", "2022–2024")
REVISIT = {"2020–2021": "6-day pairs, 2020–2021 (S1A + S1B)", "2022–2024": "12-day pairs, 2022–2024 (S1A)"}

# Reference palette (dataviz skill): categorical slots 1–3 validated all-pairs; ink and chrome tokens.
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, MUTED, GRIDC, AXIS, SURFACE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"
SEQ = LinearSegmentedColormap.from_list("seq_blue", ["#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])
DIV = LinearSegmentedColormap.from_list("div_blue_red", [BLUE, "#f0efec", "#e34948"])

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": AXIS, "axes.linewidth": 0.8, "axes.labelcolor": INK2, "axes.titlecolor": INK,
    "axes.titlesize": 10, "axes.titleweight": "bold", "axes.titlelocation": "left", "axes.labelsize": 9,
    "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "axes.grid": True, "grid.color": GRIDC, "grid.linewidth": 0.6, "grid.linestyle": "-",
    "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False, "legend.fontsize": 8,
    "font.size": 9, "text.color": INK,
})


def load():
    st = xr.open_dataset(SIL / "static_40m.nc")
    px = pd.read_csv(DLV / "field_first" / "plot_pixels.csv").set_index("plot")
    cells = pd.read_csv(X079 / "cells.csv")
    desc = pd.read_csv(X079 / "descriptors.csv").set_index("cell")
    part = pd.read_csv(X079 / "partial.csv")
    pairs = pd.read_csv(DLV / "field_p6_short_pairs_x065" / "p6_pairs_2020_2024.csv")
    return st, px, cells, desc, part, pairs


def cell_grid(values: dict) -> np.ndarray:
    return np.array([[values[c] for c in row] for row in GRID], float)


def no_grid(ax):
    ax.grid(False)
    ax.set_xticks([]), ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


# ------------------------------------------------------------------------------------------- figure 1

def fig_site(st, px, cells, out):
    r0, c0 = int(px.loc["P6", "row"]), int(px.loc["P6", "col"])
    E0, N0 = float(st.x[c0]), float(st.y[r0])
    fig = plt.figure(figsize=(12, 9.2))
    gs = fig.add_gridspec(2, 4, height_ratios=[1.35, 1], hspace=0.12, wspace=0.12)

    # (a) zones around the mat, 40 m
    ax = fig.add_subplot(gs[0, :2])
    z = st.zone.values
    rr, cc = np.where(z == 1)
    r_lo, r_hi, c_lo, c_hi = rr.min() - 4, rr.max() + 5, cc.min() - 4, cc.max() + 5
    sub = z[r_lo:r_hi, c_lo:c_hi]
    xs, ys = st.x.values[c_lo:c_hi], st.y.values[r_lo:r_hi]
    ext = [xs[0] - 20, xs[-1] + 20, ys[-1] - 20, ys[0] + 20]
    ax.imshow(sub, cmap=ListedColormap([AQUA, BLUE, ORANGE, "#ecebe6"]), vmin=0.5, vmax=4.5, extent=ext,
              interpolation="nearest")
    ax.add_patch(Rectangle((E0 - 1.5 * CELL_M, N0 - 1.5 * CELL_M), 3 * CELL_M, 3 * CELL_M, fill=False, ec=INK, lw=1.4))
    for p, r in px.iterrows():
        ax.plot(r.E, r.N, "o", ms=4, mfc=INK if p == "P6" else SURFACE, mec=INK, mew=0.8)
    lab_p = px.sort_values("N")
    ax.annotate("P1", (lab_p.iloc[0].E, lab_p.iloc[0].N), xytext=(6, -2), textcoords="offset points", fontsize=7, color=INK2)
    ax.annotate("P8/P9", (lab_p.iloc[-1].E, lab_p.iloc[-1].N), xytext=(6, 2), textcoords="offset points", fontsize=7, color=INK2)
    ax.annotate("nine cells around P6", (E0 + 1.5 * CELL_M, N0), xytext=(14, 0), textcoords="offset points",
                fontsize=8, color=INK, va="center", bbox=dict(fc=SURFACE, ec="none", pad=1.5))
    no_grid(ax)
    ax.set_title("a  Zones on the 40 m grid; plots P1–P9")
    ax.legend(handles=[Patch(fc=AQUA, label="floating mat (A)"), Patch(fc=BLUE, label="lake (B)"),
                       Patch(fc=ORANGE, label="matched grassland (C)"), Patch(fc="#ecebe6", label="other ground"),
                       Line2D([], [], marker="o", ls="", mfc=SURFACE, mec=INK, label="water-table plot"),
                       Line2D([], [], marker="o", ls="", mfc=INK, mec=INK, label="P6 (laser)")],
              loc="lower left", fontsize=7, ncol=2, frameon=True, facecolor=SURFACE, edgecolor=GRIDC)
    sb = 500.0
    ax.plot([ext[1] - sb - 60, ext[1] - 60], [ext[2] + 60] * 2, color=INK, lw=2)
    ax.text(ext[1] - sb / 2 - 60, ext[2] + 90, f"{sb:.0f} m", ha="center", fontsize=7)

    # (b) orthophoto, cells, boardwalk, station
    ax = fig.add_subplot(gs[0, 2:])
    w = 1.5 * CELL_M + 10
    bounds = (E0 - w, N0 - w, E0 + w, N0 + w)
    with rasterio.open(DLV / "field_boardwalk_x060" / "orthophoto_transect.tif") as src:
        img = np.moveaxis(src.read(window=rasterio.windows.from_bounds(*bounds, src.transform)), 0, -1)
    with rasterio.open(DLV / "field_boardwalk_x060" / "boardwalk_mask.tif") as src:
        bw = src.read(1, window=rasterio.windows.from_bounds(*bounds, src.transform)) > 0
    ext_o = [bounds[0], bounds[2], bounds[1], bounds[3]]
    ax.imshow(img[..., :3], extent=ext_o)
    ax.contour(np.linspace(bounds[0], bounds[2], bw.shape[1]), np.linspace(bounds[3], bounds[1], bw.shape[0]),
               bw.astype(float), levels=[0.5], colors=[ORANGE], linewidths=1.1)
    for i, row in enumerate(GRID):
        for j, name in enumerate(row):
            x, y = E0 + (j - 1) * CELL_M, N0 - (i - 1) * CELL_M
            ax.add_patch(Rectangle((x - 20, y - 20), 40, 40, fill=False, ec="white", lw=1.0))
            ax.text(x - 17, y + 17, name, color="white", fontsize=8, fontweight="bold", va="top")
    ax.plot(px.loc["P6", "E"], px.loc["P6", "N"], "o", ms=6, mfc=INK, mec="white", mew=1.2)
    ax.set_xlim(bounds[0], bounds[2]), ax.set_ylim(bounds[1], bounds[3])
    no_grid(ax)
    ax.plot([bounds[2] - 45, bounds[2] - 5], [bounds[1] + 6] * 2, color="white", lw=2)
    ax.text(bounds[2] - 25, bounds[1] + 2, "40 m", color="white", ha="center", va="top", fontsize=7)
    ax.set_title("b  Orthophoto: boardwalk (orange), laser (dot)")

    # (c) coherence per cell
    vals = cells.median_coh_all_pairs
    vmin, vmax = np.floor(vals.min() * 10) / 10, 1.0
    for k, (track, period) in enumerate([(t, p) for t in ("ascending", "descending") for p in PERIODS]):
        ax = fig.add_subplot(gs[1, k])
        g = cells[(cells.track == track) & (cells.period == period)].set_index("cell")
        Z = cell_grid(g.median_coh_all_pairs.to_dict())
        im = ax.imshow(Z, cmap=SEQ, vmin=vmin, vmax=vmax)
        for i, row in enumerate(GRID):
            for j, name in enumerate(row):
                v = Z[i, j]
                ax.text(j, i, f"{name}\n{v:.2f}", ha="center", va="center", fontsize=7.5,
                        color="white" if (v - vmin) / (vmax - vmin) > 0.55 else INK)
        no_grid(ax)
        ax.set_title(f"{'c' if k == 0 else ''}  {track}, {'6-day' if period == '2020–2021' else '12-day'}\n{period}",
                     fontsize=9)
    cb = fig.colorbar(im, ax=fig.axes[-4:], orientation="horizontal", fraction=0.05, pad=0.04, shrink=0.6)
    cb.set_label("median coherence of consecutive pairs (each cell; same colour scale in all four)", fontsize=8)
    cb.outline.set_visible(False)
    fig.suptitle("P6 and its nine 40 m cells", x=0.125, ha="left", fontsize=12, fontweight="bold", color=INK)
    fig.savefig(out / "fig1_site_cells.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


# ------------------------------------------------------------------------------------------- figure 2

def fig_time(pairs, out):
    wtd = field.load_wtd_hourly()["P6"]
    cens = field.censored_flag(wtd)
    wtd_d = wtd[~cens.astype(bool)].resample("D").median()
    fl = pd.read_csv(DLV / "field_p6" / "laser_qc" / "laser_flags.csv", parse_dates=["time_utc"]).set_index("time_utc")
    ok = fl.surface_cm.notna() & ~fl.snow_72h & ~fl.outlier & ~fl.filled
    las = fl.surface_cm[ok].resample("D").median()
    las = las - las.median()
    t0, t1 = las.dropna().index.min() - pd.Timedelta(days=15), las.dropna().index.max() + pd.Timedelta(days=15)
    wtd_d = wtd_d[(wtd_d.index >= t0) & (wtd_d.index <= t1)]

    p = pairs[(pairs.track == "ascending") & (pairs.step == 1)].copy()
    p["t1"] = pd.to_datetime(p.pair.str[:8], utc=True)
    p["t2"] = pd.to_datetime(p.pair.str[9:], utc=True)
    p["mid"] = p.t1 + (p.t2 - p.t1) / 2
    p = p[(p.mid >= t0) & (p.mid <= t1)]
    L = p.dropna(subset=["laser_los"])

    fig, ax = plt.subplots(4, 1, figsize=(12, 9.5), sharex=True, gridspec_kw={"height_ratios": [1, 1, 1.5, 1]})
    ax[0].plot(wtd_d.index, wtd_d.values, color=BLUE, lw=1.2)
    ax[0].set_ylabel("cm")
    ax[0].set_title("a  Water table at P6 (daily median; sensor-floor stretches removed)")
    ax[1].plot(las.index, las.values, color=INK, lw=1.0)
    ax[1].set_ylabel("cm")
    ax[1].set_title("b  Mat surface at P6, laser (daily median of QC'd hours, relative to its median; up = rising)")

    for _, r in L.iterrows():
        ax[2].plot([r.mid, r.mid], [r.laser_los, r.s1], color=GRIDC, lw=1.0, zorder=1)
    ax[2].plot(L.mid, L.laser_los, "o", ms=4.5, mfc=SURFACE, mec=INK, mew=1.0, ls="", zorder=3,
               label="laser, change over the pair (LOS)")
    ax[2].plot(L.mid, L.s1, "o", ms=4.5, color=BLUE, mec=SURFACE, mew=0.6, ls="", zorder=4,
               label="Sentinel-1, P6 cell minus grassland")
    for y in (QW, -QW):
        ax[2].axhline(y, color=MUTED, lw=0.7)
    ax[2].text(t1, QW, " ±λ/4", color=MUTED, fontsize=7, va="center")
    ax[2].axhline(0, color=AXIS, lw=0.8)
    ax[2].set_ylabel("mm, line of sight")
    ax[2].set_title(f"c  Each consecutive ascending pair with a laser value (n = {len(L)}): "
                    "laser and radar change side by side, grey = their difference")
    ax[2].legend(loc="lower left", ncol=2)

    for dt, mk in ((6, "o"), (12, "s")):
        g = p[p.dt == dt]
        ax[3].plot(g.mid, g.coh, mk, ms=4, color=BLUE if dt == 6 else ORANGE, mec=SURFACE, mew=0.5, ls="",
                   label=f"{dt}-day pairs (median {g.coh.median():.2f})")
    g = p[~p.dt.isin([6, 12])]
    if len(g):
        ax[3].plot(g.mid, g.coh, "^", ms=4, color=MUTED, ls="", label=f"longer pairs (n = {len(g)})")
    ax[3].set_ylim(0, 1)
    ax[3].set_ylabel("coherence")
    ax[3].set_title("d  Coherence of every consecutive ascending pair at the P6 cell")
    ax[3].legend(loc="lower left", ncol=3)
    ax[3].set_xlim(t0, t1)
    fig.suptitle("P6 through time: water table, laser and Sentinel-1 (ascending, 16:36 UTC)", x=0.125, ha="left",
                 fontsize=12, fontweight="bold")
    fig.tight_layout()
    fig.savefig(out / "fig2_p6_through_time.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    return {"n_pairs_laser": len(L), "t0": las.dropna().index.min().date(), "t1": las.dropna().index.max().date()}


# ------------------------------------------------------------------------------------------- figure 3

def fig_scatter(cells, out):
    ser = pd.read_csv(X079 / "cells_series.csv")
    ser = ser[(ser.cell == "P6")].dropna(subset=["laser_los"])
    fig, axs = plt.subplots(2, 2, figsize=(10.5, 10), sharex=True, sharey=True)
    lim = np.ceil(np.nanmax(np.abs(ser[["laser_los", "s1_cell"]].to_numpy())) / 5) * 5
    for i, track in enumerate(("ascending", "descending")):
        for j, period in enumerate(PERIODS):
            ax = axs[i, j]
            g = ser[(ser.track == track) & (ser.period == period)]
            s = cells[(cells.track == track) & (cells.period == period) & (cells.cell == "P6")].iloc[0]
            x = np.array([-lim, lim])
            ax.fill_between(x, (s.m - 1.96 * s.m_sd) * x, (s.m + 1.96 * s.m_sd) * x, color=ORANGE, alpha=0.15, lw=0)
            ax.plot(x, x, color=INK2, lw=1.0, label="1:1 (phase = laser)")
            ax.plot(x, s.m * x, color=ORANGE, lw=1.6, label="motion coefficient m (± 1.96 SD)")
            for v in (QW, -QW):
                ax.axvline(v, color=MUTED, lw=0.7)
            sc = ax.scatter(g.laser_los, g.s1_cell, c=g.coh_cell, cmap=SEQ, vmin=0, vmax=1, s=34,
                            edgecolors=SURFACE, linewidths=0.6, zorder=3)
            ax.set_xlim(-lim, lim), ax.set_ylim(-lim, lim)
            ax.set_aspect("equal")
            p_txt = "< 0.001" if s.p_circ < 0.001 else f"{s.p_circ:.3f}"
            ax.text(0.03, 0.97, f"n = {int(s.n)} pairs\nr = {s.r:.2f} (circular-shift p {p_txt})\n"
                                f"m = {s.m:.2f} ± {s.m_sd:.2f}\nsame sign {s.same_sign:.0%} of {int(s.n_cmp)}\n"
                                f"median coherence {s.median_coh_laser_pairs:.2f}",
                    transform=ax.transAxes, va="top", fontsize=8, color=INK,
                    bbox=dict(fc=SURFACE, ec="none", alpha=0.85, pad=2), zorder=5)
            ax.set_title(f"{'abcd'[2 * i + j]}  {track}, {REVISIT[period]}", fontsize=9)
            if i == 1:
                ax.set_xlabel("laser change over the pair, LOS (mm)")
            if j == 0:
                ax.set_ylabel("Sentinel-1 phase change, P6 cell (mm LOS)")
    axs[0, 0].text(QW, lim * 0.6, " λ/4", color=MUTED, fontsize=7)
    axs[1, 1].legend(loc="lower right", fontsize=7)
    cb = fig.colorbar(sc, ax=axs, orientation="horizontal", fraction=0.04, pad=0.07, shrink=0.5)
    cb.set_label("coherence of the pair at P6", fontsize=8)
    cb.outline.set_visible(False)
    fig.suptitle("Does the phase at P6 follow the laser? Consecutive pairs, by track and revisit",
                 x=0.125, ha="left", fontsize=12, fontweight="bold")
    fig.savefig(out / "fig3_p6_phase_vs_laser.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


# ------------------------------------------------------------------------------------------- figure 4

def fig_cells(cells, desc, part, out):
    fig = plt.figure(figsize=(12, 8.6))
    gs = fig.add_gridspec(2, 4, height_ratios=[1, 1.15], hspace=0.30, wspace=0.15)
    hi = float(np.nanmax(np.abs(cells.m - 1)))
    norm = TwoSlopeNorm(vcenter=1.0, vmin=1 - hi, vmax=1 + hi)
    for k, (track, period) in enumerate([(t, p) for t in ("ascending", "descending") for p in PERIODS]):
        ax = fig.add_subplot(gs[0, k])
        g = cells[(cells.track == track) & (cells.period == period)].set_index("cell")
        Z = cell_grid(g.m.to_dict())
        im = ax.imshow(Z, cmap=DIV, norm=norm)
        for i, row in enumerate(GRID):
            for j, name in enumerate(row):
                lo, up = g.loc[name, "m"] - 1.96 * g.loc[name, "m_sd"], g.loc[name, "m"] + 1.96 * g.loc[name, "m_sd"]
                ax.text(j, i, f"{name}\n{g.loc[name, 'm']:.2f}\n[{lo:.2f}, {up:.2f}]", ha="center", va="center",
                        fontsize=6.8, color=INK)
        no_grid(ax)
        ax.set_title(f"{'a' if k == 0 else ''}  {track}, {'6-day' if period == '2020–2021' else '12-day'}\n{period}",
                     fontsize=9)
    cb = fig.colorbar(im, ax=fig.axes[:4], orientation="horizontal", fraction=0.05, pad=0.04, shrink=0.6)
    cb.set_label("motion coefficient m against the P6 laser (1 = the cell's phase moves as the laser; "
                 "95 % interval in brackets)", fontsize=8)
    cb.outline.set_visible(False)

    style = {("ascending", "2020–2021"): ("o", BLUE, -1.5), ("ascending", "2022–2024"): ("s", ORANGE, -0.5),
             ("descending", "2020–2021"): ("^", AQUA, 0.5), ("descending", "2022–2024"): ("v", MUTED, 1.5)}
    dodge = {"depth_mat_edge_m": 1.6, "lidar_woody_pct": 0.12}
    fp = lambda v: "< 0.001" if v < 0.001 else (f"{v:.3f}" if v < 0.01 else f"{v:.2f}")  # noqa: E731
    axd = fig.add_subplot(gs[1, :2])
    axv = fig.add_subplot(gs[1, 2:], sharey=axd)
    for (track, period), (mk, col, off) in style.items():
        g = cells[(cells.track == track) & (cells.period == period)].set_index("cell").loc[desc.index]
        pr = part[(part.track == track) & (part.period == period) & (part.slope == "m")].iloc[0]
        lab = f"{track} {'6-day' if period == '2020–2021' else '12-day'}"
        for ax, xcol, key in ((axd, "depth_mat_edge_m", "depth_given_cover"), (axv, "lidar_woody_pct", "cover_given_depth")):
            ax.errorbar(desc[xcol] + off * dodge[xcol], g.m, yerr=1.96 * g.m_sd, fmt=mk, ms=5.5, color=col, mec=SURFACE, mew=0.6,
                        elinewidth=0.8, capsize=0, alpha=0.9,
                        label=f"{lab}: partial r {pr['r_' + key]:+.2f} (p {fp(pr['p_' + key])})")
    for ax in (axd, axv):
        ax.axhline(1, color=INK2, lw=0.8)
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2, fontsize=7)
    a12 = cells[(cells.track == "ascending") & (cells.period == "2022–2024")].set_index("cell")
    for name in ("W", "P6", "NE"):
        axd.annotate(name, (desc.loc[name, "depth_mat_edge_m"] - 0.5 * dodge["depth_mat_edge_m"], a12.loc[name, "m"]), xytext=(5, 4),
                     textcoords="offset points", fontsize=7, color=INK2)
    axd.set_xlabel("depth into the mat: distance to its edge (m)\n(partial r holds woody cover fixed)")
    axv.set_xlabel("woody cover in the cell: LiDAR canopy > 0.5 m (%)\n(partial r holds depth fixed)")
    axd.set_ylabel("motion coefficient m (± 1.96 SD)")
    axd.set_title("b  m vs depth into the mat")
    axv.set_title("c  m vs woody cover")
    plt.setp(axv.get_yticklabels(), visible=False)
    fig.suptitle("Nine cells around P6: how much of the laser's motion each cell's phase carries",
                 x=0.125, ha="left", fontsize=12, fontweight="bold")
    fig.savefig(out / "fig4_cells_motion.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


# ------------------------------------------------------------------------------------------- README

def readme(out, cells, desc, part, tinfo):
    c = cells.set_index(["track", "period", "cell"])
    pa = part[part.slope == "m"].set_index(["track", "period"])
    a6, a12 = c.loc[("ascending", "2020–2021", "P6")], c.loc[("ascending", "2022–2024", "P6")]
    d6, d12 = c.loc[("descending", "2020–2021", "P6")], c.loc[("descending", "2022–2024", "P6")]
    coh = cells.groupby(["track", "period"]).median_coh_all_pairs.agg(["min", "max"])
    asc = pa.loc["ascending"]
    L = ["# X-080 — P6 and its nine cells, from verified data (exploratory)", "",
         "Generated by `05_code/SAr-adapted/scripts/field_p6_figures_x080.py` from X-079's tables, X-065's pairs, the "
         "QC'd laser (X-058), the P6 water table, the data cube and X-060's orthophoto and boardwalk outline. Every "
         "number below and on the figures is computed. These replace nothing by name: the withdrawn 2026-10-06 "
         "figures are documented in `99_archive/withdrawn_2026-10-06/README.md`.", "",
         "## 1. Site and cells", "", "![site](fig1_site_cells.png)", "",
         "Median coherence of consecutive pairs per cell: "
         + "; ".join(f"{t} {p}: {coh.loc[(t, p), 'min']:.2f}–{coh.loc[(t, p), 'max']:.2f}" for t, p in coh.index)
         + ". The nine cells differ little in coherence; the revisit and the overpass time matter far more than "
           "the position of the cell.", "",
         "## 2. P6 through time", "", "![time](fig2_p6_through_time.png)", "",
         f"First to last QC'd laser day: {tinfo['t0']} to {tinfo['t1']}; {tinfo['n_pairs_laser']} consecutive ascending "
         "pairs have a laser value at both dates. Each pair's laser and radar change are drawn side by side at the "
         "pair's middle date, never joined into a line: a pair is a change, not a level.", "",
         "## 3. Phase vs laser at P6", "", "![scatter](fig3_p6_phase_vs_laser.png)", "",
         f"- Ascending, 6-day: r = {a6.r:.2f}, m = {a6.m:.2f} ± {a6.m_sd:.2f} (n = {int(a6.n)}).",
         f"- Ascending, 12-day: r = {a12.r:.2f}, m = {a12.m:.2f} ± {a12.m_sd:.2f} (n = {int(a12.n)}).",
         f"- Descending, 6-day: r = {d6.r:.2f}, m = {d6.m:.2f} ± {d6.m_sd:.2f} (n = {int(d6.n)}).",
         f"- Descending, 12-day: r = {d12.r:.2f}, m = {d12.m:.2f} ± {d12.m_sd:.2f} (n = {int(d12.n)}).", "",
         "m is X-065's motion coefficient: s1 = m·laser + e·ΔWTD, coherence-weighted, through the origin; the line "
         "drawn is the motion term alone. Six-day pairs on both tracks follow the laser; at twelve days the dusk "
         "track still does, the dawn track does not (dew, X-050).", "",
         "## 4. The nine cells", "", "![cells](fig4_cells_motion.png)", "",
         f"Ascending, holding woody cover fixed, depth into the mat keeps partial r = "
         f"{asc.r_depth_given_cover.min():+.2f} to {asc.r_depth_given_cover.max():+.2f}; holding depth fixed, woody "
         f"cover keeps {asc.r_cover_given_depth.min():+.2f} to {asc.r_cover_given_depth.max():+.2f}. On the ascending "
         "track m rises with depth into the mat at both revisits — the interior of the raft moves more than its edge "
         "(fusion v3 WP3); at six days every cell stays below the laser (m < 1), at twelve days the deeper cells exceed "
         "it. Descending shows no pattern. Nine neighbouring cells, one laser: exploratory.", "",
         "## Reading notes", "",
         "- Phase: cell minus the matched-grassland median, mm LOS, + towards the satellite (X-065 conventions).",
         "- Laser: LOS = vertical × cos(incidence); snow (72 h), outliers and gap-filled hours removed (X-058).",
         "- λ/4: beyond it a single pair's phase is ambiguous; points outside it are not usable as motion.",
         "- Colours: reference palette of the project's figure guidance, validated for colour-vision deficiency; "
         "one y-axis per panel."]
    (out / "README.md").write_text("\n".join(L) + "\n")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    st, px, cells, desc, part, pairs = load()
    fig_site(st, px, cells, OUT)
    tinfo = fig_time(pairs, OUT)
    fig_scatter(cells, OUT)
    fig_cells(cells, desc, part, OUT)
    readme(OUT, cells, desc, part, tinfo)
    print(f"X-080 written to {OUT}")


if __name__ == "__main__":
    main()

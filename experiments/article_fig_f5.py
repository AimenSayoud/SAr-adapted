"""F5, final version — coherence cannot choose the pairs that hold the motion, and where on the mat the six-day pairs
keep phase (merges the drafts' F5 and F5b).

(a, b) For each track, the share of a class's cells whose coherence is above the pair's own open-water floor (the 95th
percentile of open-water coherence in that pair), by pair length: median and interquartile range over the pairs of
each length, for the mat, the matched grassland and the hard targets. (c, d) For each track, the share of six-day pairs
above the floor in every mat and grassland cell, with the mat's outline on the 40 m grid and the depth that splits its
interior from its edge (the median depth in the mat, as ``tables/short_pairs_by_depth.csv``). The pair table and the
six-day counts are computed by the C-067 functions that wrote the article's tables, so the figure and the registered
numbers come from one computation. Writes only to the hub ``08_deliverables/article_plan/figures/final/``.

    make run E=article_fig_f5
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from cmcrameri import cm as cmc  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from matplotlib_scalebar.scalebar import ScaleBar  # noqa: E402

from insar_wetlands import figstyle as fs  # noqa: E402
from insar_wetlands.paths import hub_root  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import article_fig_f1 as F1  # noqa: E402  (cell outlines on the 40 m grid)
import article_figures_c067 as C67  # noqa: E402  (the pair table and six-day counts behind the article's tables)
from interior_edge_x128 import near  # noqa: E402  (X-128's exclusion: cells within n cells of a plot or the boardwalk)

HUB = hub_root()
OUT = HUB / "08_deliverables" / "article_plan" / "figures" / "final"
PLOTS = HUB / "08_deliverables" / "field_first" / "plot_pixels.csv"
TABLE = HUB / "08_deliverables" / "article_plan" / "tables" / "plot_phase_share.csv"
X128_CELLS = 3                                                      # X-128's primary exclusion (120 m)
BINS = ((6, 6, "6"), (12, 12, "12"), (18, 18, "18"), (24, 24, "24"), (36, 36, "36"), (48, 48, "48"),
        (100, 400, "120–396"), (700, 800, "720–744"))
CLASSES = (("mat", "mat", "o"), ("grassland", "grassland", "s"), ("hard targets", "hard targets", "^"))


def panel_lengths(ax, t: pd.DataFrame, track: str) -> list[dict]:
    rows = []
    xs = np.arange(len(BINS))
    for k, (cls, col, mk) in enumerate(CLASSES):
        med, q1, q3 = [], [], []
        for lo, hi, lab in BINS:
            v = t.loc[(t.days >= lo) & (t.days <= hi), col].dropna()
            med.append(v.median()), q1.append(v.quantile(0.25)), q3.append(v.quantile(0.75))
            rows.append({"track": track, "class": cls, "days": lab, "n_pairs": int(len(v)), "median": v.median()})
        x = xs + (k - 1) * 0.16
        c = fs.CLASS[cls]
        mec = fs.INK if cls == "grassland" else "white"
        ax.vlines(x, q1, q3, color=c, lw=1.0)
        ax.plot(x, med, "-", color=c, lw=1.0, alpha=0.6)
        ax.plot(x, med, mk, ms=3.8, color=c, mec=mec, mew=0.45, label=cls)
    ax.axhline(0.05, color=fs.MUTED, lw=0.6, ls=(0, (1, 1.5)), label="open water (5 %, by construction)")
    n = [int(((t.days >= lo) & (t.days <= hi)).sum()) for lo, hi, _ in BINS]
    ax.set_xticks(xs, [f"{b[2]}\n({k})" for b, k in zip(BINS, n)], rotation=0, fontsize=7)
    top = max(float(np.nanmax(t.loc[t.days == 6, col].quantile(0.75))) for _, col, _ in CLASSES)
    ax.annotate("six-day pairs: they follow\nthe laser at P6 (Fig. 4)", xy=(0, top + 0.02), xytext=(0.35, 0.86),
                fontsize=7, color=fs.MUTED, ha="left", va="center",
                arrowprops=dict(arrowstyle="-", color=fs.MUTED, lw=0.5))
    ax.set_xlim(-0.5, len(BINS) - 0.5)
    ax.set_ylim(0, 1)
    ax.set_xlabel("Pair length in days (number of pairs)")
    ax.set_title(track, fontsize=7.5, pad=2)
    return rows


def panel_map(ax, st, masks, extra: dict, plots: pd.DataFrame):
    x, y = st.x.values, st.y.values
    share = extra["above6"] / max(extra["n6"], 1)
    show = np.where(masks["mat"] | masks["grassland"], share, np.nan)
    ext = [x[0] - 20, x[-1] + 20, y[-1] - 20, y[0] + 20]
    im = ax.imshow(show, extent=ext, cmap=cmc.batlow, vmin=0, vmax=0.8, origin="upper", interpolation="none", zorder=1)
    lake = np.where(st.zone.values == 2, 1.0, np.nan)
    ax.imshow(lake, extent=ext, cmap=matplotlib.colors.ListedColormap([fs.CLASS["open water"]]), alpha=0.6,
              origin="upper", interpolation="none", zorder=1)
    F1.cell_outline(ax, x, y, masks["mat"], colors=fs.INK, linewidths=0.9, zorder=3)
    anchors = (st.plot_cell.values > 0) | (np.nan_to_num(st.boardwalk_share.values) > 0)
    F1.cell_outline(ax, x, y, masks["mat"] & near(anchors, X128_CELLS), colors="white", linewidths=0.9,
                    linestyles=(0, (2, 1.5)), zorder=3)
    bw = np.nan_to_num(st.boardwalk_share.values) > 0
    F1.cell_outline(ax, x, y, bw, colors="#cf222e", linewidths=0.7, zorder=4)
    ax.scatter(plots.E, plots.N, s=7, color="white", edgecolors=fs.INK, lw=0.6, zorder=5)
    p6 = plots.loc["P6"]
    ax.annotate("P6", xy=(p6.E, p6.N), xytext=(p6.E + 230, p6.N + 170), fontsize=7, fontweight="bold",
                arrowprops=dict(arrowstyle="-", color=fs.INK, lw=0.6), zorder=6,
                bbox=dict(fc="white", ec="none", pad=1.0, alpha=0.9))
    rr, cc = np.where(masks["mat"])
    pad = 10
    ax.set_xlim(x[max(cc.min() - pad, 0)] - 20, x[min(cc.max() + pad, len(x) - 1)] + 20)
    ax.set_ylim(y[min(rr.max() + pad, len(y) - 1)] - 20, y[max(rr.min() - pad, 0)] + 20)
    ax.set_aspect("equal")
    ax.set_xticks([]), ax.set_yticks([])
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(True)
    ax.add_artist(ScaleBar(1, "m", location="lower left", length_fraction=0.2, box_alpha=0.85,
                           font_properties={"size": 7}, fixed_value=500))
    return im


def plot_shares(st, masks, extras, plots) -> pd.DataFrame:
    """Each plot's six-day share and its percentile among the mat's cells, per track (the representativeness of the
    field sites, written for the text)."""
    rows = []
    for tr, e in extras.items():
        share = e["above6"] / max(e["n6"], 1)
        mat = share[masks["mat"]]
        mat = mat[np.isfinite(mat)]
        for p, r in plots.iterrows():
            v = float(share[int(r.row), int(r.col)])
            rows.append({"track": tr, "plot": p, "share": v, "percentile_in_mat": float(np.mean(mat <= v) * 100),
                         "mat_median": float(np.median(mat)), "mat_p90": float(np.percentile(mat, 90)),
                         "six_day_pairs": int(e["n6"])})
    t = pd.DataFrame(rows)
    t.to_csv(TABLE, index=False)
    return t


def panel_cdf(ax, masks, extras, shares: pd.DataFrame):
    for tr, e in extras.items():
        share = e["above6"] / max(e["n6"], 1)
        v = np.sort(share[masks["mat"]][np.isfinite(share[masks["mat"]])])
        ls = "-" if tr == "ascending" else (0, (3.5, 1.8))
        ax.step(v, np.arange(1, len(v) + 1) / len(v) * 100, where="post", color=fs.TRACK[tr], lw=1.1, ls=ls, label=tr)
        p6 = shares[(shares.track == tr) & (shares["plot"] == "P6")].iloc[0]
        ax.plot(p6.share, p6.percentile_in_mat, "*", ms=8, color=fs.TRACK[tr], mec=fs.INK, mew=0.5, zorder=4)
        ax.annotate(f"P6 ({p6.percentile_in_mat:.0f}th)", xy=(p6.share, p6.percentile_in_mat),
                    xytext=(p6.share - 0.05, p6.percentile_in_mat - (14 if tr == "ascending" else 28)), fontsize=7,
                    color=fs.TRACK[tr], ha="right", arrowprops=dict(arrowstyle="-", color=fs.TRACK[tr], lw=0.5))
    ax.set_xlim(0, 0.8), ax.set_ylim(0, 102)
    ax.set_xlabel("share of six-day pairs\nabove the open-water floor")
    ax.set_ylabel("mat cells at or below (%)")
    ax.legend(loc="lower right", fontsize=7, handlelength=1.6)


def main() -> None:
    st, masks, water = C67.static()
    plots = pd.read_csv(PLOTS).set_index("plot")
    tables, extras = {}, {}
    for tr in ("ascending", "descending"):
        tables[tr], extras[tr] = C67.pair_table(tr, masks, water)
    shares = plot_shares(st, masks, extras, plots)
    with fs.style():
        fig = plt.figure(figsize=fs.size("full", 162))
        gs = fig.add_gridspec(2, 1, height_ratios=[0.8, 1.05], hspace=0.30, left=0.075, right=0.985, bottom=0.135,
                              top=0.935)
        top = gs[0].subgridspec(1, 2, wspace=0.08)
        bot = gs[1].subgridspec(1, 5, width_ratios=[1, 1, 0.045, 0.42, 0.66], wspace=0.06)
        a, b = fig.add_subplot(top[0]), fig.add_subplot(top[1])
        c, d = fig.add_subplot(bot[0]), fig.add_subplot(bot[1])
        cax = fig.add_subplot(bot[2].subgridspec(3, 1, height_ratios=[0.14, 0.72, 0.14])[1])
        e = fig.add_subplot(bot[4])
        meds = panel_lengths(a, tables["ascending"], "ascending") + panel_lengths(b, tables["descending"],
                                                                                    "descending")
        a.set_ylabel("share of the class's cells above the pair's\nopen-water floor (median, interquartile range)")
        b.tick_params(labelleft=False)
        h, l = a.get_legend_handles_labels()
        fig.legend(h, l, loc="upper center", bbox_to_anchor=(0.49, 1.0), handlelength=1.2, ncol=4, columnspacing=1.2)
        for ax in (a, b):
            ax.set_ylim(0, 1.05)
        im = panel_map(c, st, masks, extras["ascending"], plots)
        panel_map(d, st, masks, extras["descending"], plots)
        c.set_title("ascending, six-day pairs", fontsize=7.5, pad=2)
        d.set_title("descending, six-day pairs", fontsize=7.5, pad=2)
        cb = fig.colorbar(im, cax=cax, extend="max")
        cb.set_label("share of six-day pairs above the open-water floor")
        panel_cdf(e, masks, extras, shares)
        excl_handle = (Line2D([], [], color=fs.INK, lw=2.2),
                       Line2D([], [], color="white", lw=0.9, ls=(0, (2, 1.5))))
        fig.legend(handles=[Line2D([], [], color=fs.INK, lw=0.9, label="floating mat (cells outside it: grassland)"),
                            Patch(fc=fs.CLASS["open water"], ec=fs.INK, lw=0.4, alpha=0.7, label="residual lake"),
                            Line2D([], [], ls="", marker="o", ms=4, mfc="white", mec=fs.INK, label="plots P1–P9"),
                            Line2D([], [], color="#cf222e", lw=0.9, label="boardwalk cells"),
                            excl_handle],
                   labels=["floating mat (cells outside it: grassland)",
                           "residual lake",
                           "plots P1–P9",
                           "boardwalk cells",
                           "within 120 m of a plot or the boardwalk (left out in X-128)"],
                   loc="lower center", ncol=3, bbox_to_anchor=(0.5, 0.005), fontsize=7, frameon=True,
                   facecolor="#f8f9fa", edgecolor="#d0d7de", framealpha=1.0)
        for ax, lab in ((a, "a"), (b, "b"), (c, "c"), (d, "d"), (e, "e")):
            fs.panel_label(ax, lab, x=0.0, y=1.02)
        rec = fs.export(fig, OUT / "F5_coherence_floor", provenance={
            "inputs": ["cube silver/static_40m.nc (zone, hard_target, wc_water, plot_cell, boardwalk_share)",
                       "cube silver/ifg_<set>_<track>.nc (coherence of every pair; C-067 pair_table)",
                       str(PLOTS.relative_to(HUB))],
            "brief": "08_deliverables/article_plan/figures/briefs/F5.md",
            "floor": "per pair: 95th percentile of open-water (WorldCover water ≥ 0.9) coherence",
            "six_day_pairs": {tr: int(extras[tr]["n6"]) for tr in extras},
            "x128_exclusion_cells": X128_CELLS, "medians_by_length": meds,
            "table": str(TABLE.relative_to(HUB))})
        plt.close(fig)
    print("F5 →", OUT, rec["width_mm"], "×", rec["height_mm"], "mm")


if __name__ == "__main__":
    main()

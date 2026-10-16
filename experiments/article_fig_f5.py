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

HUB = hub_root()
OUT = HUB / "08_deliverables" / "article_plan" / "figures" / "final"
BINS = ((6, 6, "6"), (12, 12, "12"), (18, 18, "18"), (24, 24, "24"), (36, 36, "36"), (48, 48, "48"),
        (100, 400, "120–396"), (700, 800, "720–744"))
CLASSES = (("mat", "mat"), ("grassland", "grassland"), ("hard targets", "hard targets"))


def panel_lengths(ax, t: pd.DataFrame, track: str) -> list[dict]:
    rows = []
    xs = np.arange(len(BINS))
    for k, (cls, col) in enumerate(CLASSES):
        med, q1, q3 = [], [], []
        for lo, hi, lab in BINS:
            v = t.loc[(t.days >= lo) & (t.days <= hi), col].dropna()
            med.append(v.median()), q1.append(v.quantile(0.25)), q3.append(v.quantile(0.75))
            rows.append({"track": track, "class": cls, "days": lab, "n_pairs": int(len(v)), "median": v.median()})
        x = xs + (k - 1) * 0.16
        c = fs.CLASS[cls]
        ax.vlines(x, q1, q3, color=c, lw=1.0)
        ax.plot(x, med, "-", color=c, lw=1.0, alpha=0.6)
        ax.plot(x, med, "o", ms=3.6, color=c, mec="white", mew=0.4, label=cls)
    ax.axhline(0.05, color=fs.MUTED, lw=0.6, ls=(0, (1, 1.5)), label="open water (5 %, by construction)")
    n = [int(((t.days >= lo) & (t.days <= hi)).sum()) for lo, hi, _ in BINS]
    ax.set_xticks(xs, [f"{b[2]}\n({k})" for b, k in zip(BINS, n)], rotation=0, fontsize=6.5)
    ax.set_xlim(-0.5, len(BINS) - 0.5)
    ax.set_ylim(0, 1)
    ax.set_xlabel("Pair length in days (number of pairs)")
    ax.set_title(track, fontsize=7.5, pad=2)
    return rows


def panel_map(ax, st, masks, extra: dict, split: float):
    x, y = st.x.values, st.y.values
    share = extra["above6"] / max(extra["n6"], 1)
    show = np.where(masks["mat"] | masks["grassland"], share, np.nan)
    ext = [x[0] - 20, x[-1] + 20, y[-1] - 20, y[0] + 20]
    im = ax.imshow(show, extent=ext, cmap=cmc.batlow, vmin=0, vmax=0.8, origin="upper", interpolation="none", zorder=1)
    lake = np.where(st.zone.values == 2, 1.0, np.nan)
    ax.imshow(lake, extent=ext, cmap=matplotlib.colors.ListedColormap([fs.CLASS["open water"]]), alpha=0.6,
              origin="upper", interpolation="none", zorder=1)
    F1.cell_outline(ax, x, y, masks["mat"], colors=fs.INK, linewidths=0.9, zorder=3)
    interior = masks["mat"] & (st.depth_in_mat_m.values > split)
    F1.cell_outline(ax, x, y, interior, colors="white", linewidths=0.9, linestyles=(0, (2, 1.5)), zorder=3)
    rr, cc = np.where(masks["mat"])
    pad = 10
    ax.set_xlim(x[max(cc.min() - pad, 0)] - 20, x[min(cc.max() + pad, len(x) - 1)] + 20)
    ax.set_ylim(y[min(rr.max() + pad, len(y) - 1)] - 20, y[max(rr.min() - pad, 0)] + 20)
    ax.set_aspect("equal")
    ax.set_xticks([]), ax.set_yticks([])
    for s in ("top", "right"):
        ax.spines[s].set_visible(True)
    ax.add_artist(ScaleBar(1, "m", location="lower left", length_fraction=0.2, box_alpha=0.85,
                           font_properties={"size": 7}, fixed_value=500))
    return im


def main() -> None:
    st, masks, water = C67.static()
    split = float(np.nanmedian(st.depth_in_mat_m.values[masks["mat"]]))
    tables, extras = {}, {}
    for tr in ("ascending", "descending"):
        tables[tr], extras[tr] = C67.pair_table(tr, masks, water)
    with fs.style():
        fig = plt.figure(figsize=fs.size("full", 140))
        gs = fig.add_gridspec(2, 2, height_ratios=[0.8, 1.15], hspace=0.28, wspace=0.08, left=0.075, right=0.9,
                              bottom=0.075, top=0.93)
        a, b = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
        c, d = fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])
        meds = panel_lengths(a, tables["ascending"], "ascending") + panel_lengths(b, tables["descending"],
                                                                                    "descending")
        a.set_ylabel("Share of the class's cells\nabove the pair's open-water floor")
        b.tick_params(labelleft=False)
        h, l = a.get_legend_handles_labels()
        fig.legend(h, l, loc="upper center", bbox_to_anchor=(0.49, 1.0), handlelength=1.2, ncol=4, columnspacing=1.2)
        for ax in (a, b):
            ax.set_ylim(0, 1.05)
        im = panel_map(c, st, masks, extras["ascending"], split)
        panel_map(d, st, masks, extras["descending"], split)
        c.set_title("ascending, six-day pairs", fontsize=7.5, pad=2)
        d.set_title("descending, six-day pairs", fontsize=7.5, pad=2)
        cax = fig.add_axes([0.91, 0.12, 0.012, 0.38])
        cb = fig.colorbar(im, cax=cax, extend="max")
        cb.set_label("Share of six-day pairs above the open-water floor")
        fig.legend(handles=[Line2D([], [], color=fs.INK, lw=0.9, label="floating mat"),
                          Line2D([], [], color="white", lw=0.9, ls=(0, (2, 1.5)), label="interior (deeper than median)"),
                          Patch(fc=fs.CLASS["open water"], alpha=0.6, label="residual lake")],
                 loc="lower center", ncol=3, bbox_to_anchor=(0.5, 0.0), fontsize=7, frameon=True,
                 facecolor="#c8ccd0", edgecolor="none", framealpha=1.0)
        for ax, l in ((a, "a"), (b, "b"), (c, "c"), (d, "d")):
            fs.panel_label(ax, l, x=0.0, y=1.02)
        rec = fs.export(fig, OUT / "F5_coherence_floor", provenance={
            "inputs": ["cube silver/static_40m.nc (zone, hard_target, wc_water, depth_in_mat_m)",
                       "cube silver/ifg_<set>_<track>.nc (coherence of every pair; C-067 pair_table)"],
            "floor": "per pair: 95th percentile of open-water (WorldCover water ≥ 0.9) coherence",
            "six_day_pairs": {tr: int(extras[tr]["n6"]) for tr in extras},
            "interior_split_depth_m": split, "medians_by_length": meds})
        plt.close(fig)
    print("F5 →", OUT, rec["width_mm"], "×", rec["height_mm"], "mm")


if __name__ == "__main__":
    main()

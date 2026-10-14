"""C-067 — draft article figures from existing data (structure: hub ``08_deliverables/article_plan/ARTICLE_PLAN.md``).

F2 inset (the laser in its 40 m cell and the cells around), F3 (schematic: what C-band sees on a floating mat), F5 (the
coherence matrix: share of class cells above each pair's open-water floor, by the pair's dates, and by pair length),
F5b (map: share of 6-day pairs above the floor per cell) and F6 (the systematic test: every try by family and
verdict). Drafts — layout and styling are agreed with Aymen before the final versions. Reads the cube and the hub;
writes only to the hub: ``08_deliverables/article_plan/figures/``. No number is written into a caption.

    make run E=article_figures_c067
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import xarray as xr  # noqa: E402
from matplotlib.patches import FancyArrowPatch, Polygon, Rectangle  # noqa: E402

from insar_wetlands.paths import cube_root, hub_root  # noqa: E402

HUB = hub_root()
CUBE = cube_root()
OUT = HUB / "08_deliverables" / "article_plan" / "figures"
SETS = ("2020_2021", "2022_2024", "matched_wt", "annual")
TRACKS = ("ascending", "descending")
INK, MUTED, GRID = "#1f2328", "#6e7781", "#d0d7de"
CLASS = {"mat": "#0b6e4f", "grassland": "#c08b2c", "hard targets": "#4a5fc1"}
VERDICT = {"worked": "#1a7f37", "partial": "#bf8700", "failed": "#cf222e", "inconclusive": "#8c959f",
           "open": "#54aeff", "not_run": "#d0d7de", "superseded": "#d0d7de"}
plt.rcParams.update({"font.size": 8, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False})


def static():
    st = xr.open_dataset(CUBE / "silver" / "static_40m.nc")
    masks = {"mat": st.zone.values == 1, "grassland": st.zone.values == 3,
             "hard targets": st.hard_target.values.astype(bool)}
    return st, masks, st.wc_water.values >= 0.9


def pair_table(track: str, masks, water) -> tuple[pd.DataFrame, dict]:
    """Per pair: dates, length, the open-water floor (95th percentile of water coherence) and each class's share above it;
    plus, for the 6-day pairs, the per-cell count above the floor."""
    rows, above6, n6 = [], None, 0
    for s in SETS:
        f = CUBE / "silver" / f"ifg_{s}_{track}.nc"
        if not f.exists():
            continue
        d = xr.open_dataset(f)
        for i in range(d.sizes["pair"]):
            coh = d.coh.isel(pair=i).values
            wv = coh[water]
            if np.isfinite(wv).sum() < 10:
                continue
            floor = float(np.nanquantile(wv, 0.95))
            rec = {"set": s, "t1": pd.Timestamp(d.t1.values[i]), "t2": pd.Timestamp(d.t2.values[i]),
                   "days": int(d.revisit_days.values[i]), "floor": floor}
            for k, m in masks.items():
                rec[k] = float(np.nanmean(coh[m] > floor))
            rows.append(rec)
            if int(d.revisit_days.values[i]) == 6:
                a = (coh > floor).astype(float)
                a[~np.isfinite(coh)] = np.nan
                above6 = a if above6 is None else np.nansum(np.stack([above6, a]), 0)
                n6 += 1
    return pd.DataFrame(rows), {"above6": above6, "n6": n6}


def f5(tabs: dict) -> None:
    fig = plt.figure(figsize=(7.2, 6.2))
    gs = fig.add_gridspec(3, 3, height_ratios=[1, 1, 0.9], hspace=0.45, wspace=0.3)
    for r, track in enumerate(TRACKS):
        t = tabs[track]
        for c, k in enumerate(CLASS):
            ax = fig.add_subplot(gs[r, c])
            sc = ax.scatter(t.t1, t.t2, c=t[k], s=6, cmap="viridis", vmin=0, vmax=1, lw=0)
            ax.set_title(f"{k} — {track}", fontsize=8, color=INK)
            ax.tick_params(labelsize=6)
            yrs = pd.to_datetime([f"{y}-01-01" for y in range(2018, 2027, 2)])
            ax.set_xticks(yrs), ax.set_xticklabels([str(y.year) for y in yrs])
            ax.set_yticks(yrs), ax.set_yticklabels([str(y.year) for y in yrs])
            if c == 0:
                ax.set_ylabel("second date")
            if r == 1:
                ax.set_xlabel("first date")
    cb = fig.colorbar(sc, ax=fig.axes, shrink=0.5, location="right")
    cb.set_label("share of class cells above the pair's open-water floor")
    ax = fig.add_subplot(gs[2, :])
    for track, ls in zip(TRACKS, ("-", "--")):
        t = tabs[track]
        for k, col in CLASS.items():
            g = t.groupby(pd.cut(t.days, [0, 7, 13, 25, 49, 100, 250, 450, 800]), observed=True)[k].median()
            x = [iv.mid for iv in g.index]
            ax.plot(x, g.values, ls, color=col, lw=1.6, marker="o", ms=3, label=f"{k}, {track}")
    ax.set_xscale("log")
    ax.set_xlabel("pair length (days)")
    ax.set_ylabel("median share above floor")
    ax.axhline(0.05, color=MUTED, lw=0.8, ls=":")
    ax.text(ax.get_xlim()[1], 0.05, " open water by construction", va="center", fontsize=6, color=MUTED)
    ax.legend(fontsize=6, ncol=3, frameon=False, loc="upper right")
    fig.suptitle("F5 (draft) — where coherence holds, against each pair's own noise floor", fontsize=9, color=INK)
    fig.savefig(OUT / "F5_coherence_matrix.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def f5b(st, masks, extra: dict) -> None:
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.4), constrained_layout=True)
    ext = [st.x.values.min() - 20, st.x.values.max() + 20, st.y.values.min() - 20, st.y.values.max() + 20]
    for ax, track in zip(axs, TRACKS):
        e = extra[track]
        share = e["above6"] / max(e["n6"], 1)
        show = np.where(masks["mat"] | masks["grassland"], share, np.nan)
        im = ax.imshow(show, extent=ext, cmap="viridis", vmin=0, vmax=1, origin="upper")
        ax.contour(st.x.values, st.y.values, masks["mat"].astype(float), levels=[0.5], colors="white", linewidths=0.8)
        rr, cc = np.where(masks["mat"])
        pad = 12
        ax.set_xlim(st.x.values[max(cc.min() - pad, 0)], st.x.values[min(cc.max() + pad, len(st.x) - 1)])
        ax.set_ylim(st.y.values[min(rr.max() + pad, len(st.y) - 1)], st.y.values[max(rr.min() - pad, 0)])
        ax.set_title(f"{track} — 6-day pairs", fontsize=8, color=INK)
        ax.set_xticks([]), ax.set_yticks([])
    cb = fig.colorbar(im, ax=axs, shrink=0.8)
    cb.set_label("share of 6-day pairs above the open-water floor")
    fig.suptitle("F5b (draft) — where on the mat the shortest pairs keep phase (mat outlined; mat and grassland shown)",
                 fontsize=9, color=INK)
    fig.savefig(OUT / "F5b_short_pair_map.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def f2_inset(st, masks) -> None:
    px = pd.read_csv(HUB / "08_deliverables" / "field_first" / "plot_pixels.csv").set_index("plot")
    r, c = int(px.loc["P6", "row"]), int(px.loc["P6", "col"])
    h = 3
    sub = st.lidar_vegh_mean_m.values[r - h:r + h + 1, c - h:c + h + 1]
    xs, ys = st.x.values[c - h:c + h + 1], st.y.values[r - h:r + h + 1]
    fig, ax = plt.subplots(figsize=(3.6, 3.4), constrained_layout=True)
    im = ax.imshow(sub, extent=[xs[0] - 20, xs[-1] + 20, ys[-1] - 20, ys[0] + 20], cmap="Greens", origin="upper")
    for xv in np.r_[xs - 20, xs[-1] + 20]:
        ax.axvline(xv, color=GRID, lw=0.5)
    for yv in np.r_[ys - 20, ys[-1] + 20]:
        ax.axhline(yv, color=GRID, lw=0.5)
    x0, y0 = st.x.values[c], st.y.values[r]
    ax.add_patch(Rectangle((x0 - 60, y0 - 60), 120, 120, fill=False, ec=INK, lw=1.2, ls="--"))
    ax.add_patch(Rectangle((x0 - 20, y0 - 20), 40, 40, fill=False, ec=INK, lw=1.6))
    for p, row in px.iterrows():
        if abs(row.E - x0) < 160 and abs(row.N - y0) < 160:
            ax.plot(row.E, row.N, "o", ms=4, color="#cf222e" if p == "P6" else INK)
            ax.text(row.E + 6, row.N + 6, p, fontsize=6, color=INK)
    ax.set_xticks([]), ax.set_yticks([])
    cb = fig.colorbar(im, ax=ax, shrink=0.8)
    cb.set_label("LiDAR mean vegetation height (m)")
    ax.set_title("F2 inset (draft) — the laser (P6) in its 40 m cell; the nine cells around (dashed)", fontsize=7, color=INK)
    fig.savefig(OUT / "F2_inset_p6_cells.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def f3() -> None:
    fig, ax = plt.subplots(figsize=(7.2, 3.2), constrained_layout=True)
    ax.set_xlim(0, 10), ax.set_ylim(0, 6), ax.axis("off")
    ax.add_patch(Rectangle((0, 0), 10, 1.6, color="#cfe3f3"))
    ax.text(0.2, 0.3, "water", fontsize=8, color=MUTED)
    ax.add_patch(Polygon([(1.5, 1.4), (8.5, 1.4), (8.5, 2.2), (1.5, 2.2)], closed=True, color="#8a6d3b"))
    ax.text(5, 1.7, "floating peat mat — rises and falls with the water table (buoyancy)", ha="center", fontsize=7, color="white")
    for xv in np.linspace(1.8, 8.2, 22):
        ax.plot([xv, xv + 0.05], [2.2, 2.55], color="#3f7d3a", lw=1)
    ax.add_patch(FancyArrowPatch((8.8, 5.6), (4.2, 2.35), arrowstyle="-|>", mutation_scale=10, color=INK, lw=1.2))
    ax.add_patch(FancyArrowPatch((4.6, 2.35), (8.9, 5.05), arrowstyle="-|>", mutation_scale=10, color=MUTED, lw=1, ls="--"))
    ax.text(6.9, 3.4, "backscatter from the mat surface", fontsize=6, color=MUTED, rotation=31)
    ax.text(8.9, 5.6, "Sentinel-1\nC-band", fontsize=7, color=INK)
    notes = ["the phase reads the mat surface, not the water level (I04)",
             "coherence follows the change of moisture between the dates (I05)",
             "a homogeneous distributed scatterer — no persistent scatterers (I07)",
             "only the shortest consecutive pairs keep the motion (I08)"]
    for i, n in enumerate(notes):
        ax.text(0.2, 5.6 - 0.45 * i, "•  " + n, fontsize=7, color=INK)
    ax.set_title("F3 (draft schematic) — what C-band sees on a floating mat", fontsize=9, color=INK)
    fig.savefig(OUT / "F3_schematic.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def f6() -> None:
    t = pd.read_csv(HUB / "08_deliverables" / "independent_review" / "01_tries" / "tries.csv")
    fams = list(dict.fromkeys(t.family))
    fig, ax = plt.subplots(figsize=(7.2, 0.32 * len(fams) + 1.2), constrained_layout=True)
    for yi, fam in enumerate(fams):
        g = t[t.family == fam]
        for xi, (_, r) in enumerate(g.iterrows()):
            ax.add_patch(Rectangle((xi, yi - 0.4), 0.92, 0.8, color=VERDICT.get(r.verdict, GRID)))
            ax.text(xi + 0.46, yi, r["id"], ha="center", va="center", fontsize=5, color="white")
    ax.set_yticks(range(len(fams)))
    ax.set_yticklabels(fams)
    ax.set_xlim(0, t.family.value_counts().max() + 0.2), ax.set_ylim(len(fams) - 0.5, -0.5)
    ax.set_xticks([])
    for k, col in VERDICT.items():
        if k in set(t.verdict):
            ax.plot([], [], "s", color=col, label=k)
    ax.legend(ncol=6, fontsize=6, frameon=False, loc="lower center", bbox_to_anchor=(0.5, 1.0))
    ax.set_title("F6 (draft) — the systematic test: every try by family and verdict (tries map)", fontsize=9, color=INK, pad=18)
    fig.savefig(OUT / "F6_systematic_test.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    st, masks, water = static()
    tabs, extra = {}, {}
    for track in TRACKS:
        tabs[track], extra[track] = pair_table(track, masks, water)
        tabs[track].to_csv(OUT / f"F5_pairs_{track}.csv", index=False)
    f5(tabs)
    f5b(st, masks, extra)
    f2_inset(st, masks)
    f3()
    f6()
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()

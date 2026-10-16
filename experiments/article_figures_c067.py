"""C-067 — draft article figures from existing data (structure: hub ``08_deliverables/article_plan/ARTICLE_PLAN.md``).

F1 (the site), F2 inset (the laser in its 40 m cell and the cells around), F3 (conceptual: what a pair measures on a
floating mat and which pairs hold the motion), F4 (pair phase against the laser by revisit), F5 (the
coherence matrix: share of class cells above each pair's open-water floor, by the pair's dates, and by pair length),
F5b (map: share of 6-day pairs above the floor per cell) and F6 (the systematic test: every route by what it varies,
with its outcome; `article_plan/f6_routes.yaml`). Drafts — layout and styling are agreed with Aymen before the final
versions. Reads the cube and the hub;
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
from matplotlib.patches import Arc, FancyArrowPatch, Rectangle  # noqa: E402

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


FT = HUB / "08_deliverables" / "insar_full_test"
OUTCOME = {"holds": ("#1a7f37", "follows the motion"), "partly": ("#bf8700", "partly"), "helps": ("#1a7f37", "helps"),
           "lost": ("#cf222e", "loses it"), "not_motion": ("#8250df", "a signal, not the motion"),
           "not_applicable": ("#8c959f", "fails its own conditions"), "no_effect": ("#57606a", "no gain / no change")}
PRODUCT = {"HyP3": "HyP3", "SLC": "own SLC", "both": "HyP3 + own SLC", "sim": "simulation"}


def latest(entry: str, rel: str) -> Path:
    """An entry's table as the article reads it: the latest re-read on corrected inputs, else the original."""
    for root in (FT / "reread_rtc", FT / "reread_corrected"):
        f = root / entry / rel
        if f.exists():
            return f
    return HUB / "08_deliverables" / entry / rel


def f1(st, masks) -> None:
    """The site: the mat, the plots, P6, the hard targets, over the LiDAR canopy; the two looks."""
    px = pd.read_csv(HUB / "08_deliverables" / "field_first" / "plot_pixels.csv").set_index("plot")
    x, y = st.x.values, st.y.values
    ext = [x.min() - 20, x.max() + 20, y.min() - 20, y.max() + 20]
    fig, ax = plt.subplots(figsize=(7.2, 6.2), constrained_layout=True)
    im = ax.imshow(st.lidar_vegh_p90_m.values, extent=ext, cmap="Greens", vmin=0, vmax=25, origin="upper")
    grass = np.where(st.zone.values == 3, 1.0, np.nan)
    ax.imshow(grass, extent=ext, cmap=matplotlib.colors.ListedColormap(["#c08b2c"]), alpha=0.55, origin="upper")
    ax.add_patch(Rectangle((0, 0), 0, 0, color="#c08b2c", alpha=0.55, label="matched grassland (reference zone)"))
    for z, col, lab in ((1, "#0b6e4f", "floating mat"), (2, "#0969da", "lake")):
        ax.contour(x, y, (st.zone.values == z).astype(float), levels=[0.5], colors=col, linewidths=1.6)
        ax.plot([], [], color=col, lw=1.6, label=lab)
    rr, cc = np.where(st.hard_target.values.astype(bool))
    ax.plot(x[cc], y[rr], "s", ms=4, mfc="none", mec="#4a5fc1", mew=1.2, label="hard targets (phase reference)")
    for p, row in px.iterrows():
        big = p == "P6"
        ax.plot(row.E, row.N, "*" if big else "o", ms=12 if big else 3.5, color="#cf222e" if big else INK,
                mec="white", mew=0.6, zorder=5)
    for p, dx in (("P1", -40), ("P9", -40), ("P6", 40)):
        row = px.loc[p]
        ax.annotate(p + (" — laser" if p == "P6" else ""), (row.E, row.N), (row.E + dx * 3, row.N + 20),
                    fontsize=7, color=INK, ha="right" if dx < 0 else "left", va="center", zorder=6,
                    arrowprops=dict(arrowstyle="-", color=INK, lw=0.5),
                    bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.8))
    rm, cm = np.where(masks["mat"])
    pad = 18
    ax.set_xlim(x[max(cm.min() - pad, 0)], x[min(cm.max() + pad, len(x) - 1)])
    ax.set_ylim(y[min(rm.max() + pad, len(y) - 1)], y[max(rm.min() - pad, 0)])
    x0, y0 = ax.get_xlim()[0], ax.get_ylim()[0]
    w, h = np.diff(ax.get_xlim())[0], np.diff(ax.get_ylim())[0]
    ax.add_patch(Rectangle((x0 + 0.01 * w, y0 + 0.01 * h), 0.5 * w, 0.17 * h, fc="white", ec=GRID, alpha=0.9, zorder=6))
    for k, (tr, e, n) in enumerate((("ascending", "los_e_asc", "los_n_asc"), ("descending", "los_e_des", "los_n_des"))):
        ve, vn = -float(np.nanmedian(st[e].values)), -float(np.nanmedian(st[n].values))   # satellite → ground
        nrm = np.hypot(ve, vn)
        cx, cy = x0 + w * (0.08 + 0.15 * k), y0 + h * 0.11
        ax.add_patch(FancyArrowPatch((cx - 0.05 * w * ve / nrm, cy - 0.05 * w * vn / nrm),
                                     (cx + 0.05 * w * ve / nrm, cy + 0.05 * w * vn / nrm), arrowstyle="-|>",
                                     mutation_scale=10, color=INK, zorder=7))
        ax.text(cx, cy - 0.035 * h, f"{tr} look", ha="center", va="top", fontsize=6, color=INK, zorder=7)
    ax.add_patch(FancyArrowPatch((x0 + 0.44 * w, y0 + 0.05 * h), (x0 + 0.44 * w, y0 + 0.14 * h), arrowstyle="-|>",
                                 mutation_scale=10, color=INK, zorder=7))
    ax.text(x0 + 0.465 * w, y0 + 0.12 * h, "N", ha="left", fontsize=7, zorder=7)
    ax.plot([x0 + 0.30 * w, x0 + 0.30 * w + 200], [y0 + 0.06 * h] * 2, color=INK, lw=2, zorder=7)
    ax.text(x0 + 0.30 * w + 100, y0 + 0.075 * h, "200 m", ha="center", fontsize=6, zorder=7)
    ax.set_xticks([]), ax.set_yticks([])
    ax.legend(fontsize=6, frameon=True, loc="upper right")
    cb = fig.colorbar(im, ax=ax, shrink=0.6)
    cb.set_label("LiDAR canopy height, 90th percentile (m)")
    ax.set_title("F1 (draft) — the Rzecin floating mat, its plots, the laser at P6 and the hard-target reference",
                 fontsize=9, color=INK)
    fig.savefig(OUT / "F1_site.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def f3() -> None:
    """Conceptual figure: (a) what a pair measures on a floating mat, (b) which pairs hold the motion over the year."""
    fig, (a, b) = plt.subplots(2, 1, figsize=(7.2, 6.0), constrained_layout=True, gridspec_kw={"height_ratios": [1.2, 1]})
    a.set_xlim(0, 10), a.set_ylim(0, 6), a.axis("off")
    for x0, lvl, lab, wet in ((0.4, 1.4, "date 1 — lower water", False), (5.4, 2.2, "date 2 — higher water", True)):
        a.add_patch(Rectangle((x0, 0), 4.2, lvl, color="#cfe3f3"))
        a.plot([x0, x0 + 4.2], [lvl, lvl], color="#0969da", lw=1)
        a.add_patch(Rectangle((x0 + 0.5, lvl - 0.35), 3.2, 0.7, color="#8a6d3b"))
        for xv in np.linspace(x0 + 0.65, x0 + 3.55, 12):
            a.plot([xv, xv + 0.04], [lvl + 0.35, lvl + 0.65], color="#3f7d3a", lw=1)
        if wet:
            for xv in np.linspace(x0 + 0.7, x0 + 3.5, 9):
                a.plot(xv, lvl + 0.42, "o", ms=2, color="#0969da")
        a.text(x0 + 2.1, 0.25, lab, ha="center", fontsize=7, color=INK)
    a.add_patch(FancyArrowPatch((5.15, 1.8), (5.15, 2.6), arrowstyle="<|-|>", mutation_scale=8, color="#cf222e"))
    a.text(5.05, 3.35, "the pair's phase:\nthe change of the\nmat surface", ha="right", fontsize=7, color="#cf222e")
    a.text(9.9, 3.9, "coherence: how much the\nsurface itself changed\n(wetness, vegetation)", ha="right", fontsize=7,
           color="#0969da")
    a.text(0.4, 5.6, "a  The mat floats: it rises with its water table (buoyancy). A short pair reads that rise at the "
                     "surface, not the water level.", fontsize=7.5, color=INK)
    b.set_xlim(0, 60), b.set_ylim(-0.2, 4.2), b.axis("off")
    t = np.arange(2, 58, 2.0)
    winter = (t > 34) & (t < 46)
    b.add_patch(Rectangle((34, -0.1), 12, 3.4, color="#eaeef2"))
    b.text(40, 3.0, "winter: snow, frost", ha="center", fontsize=7, color=MUTED)
    b.plot(t[~winter], np.zeros((~winter).sum()), "o", ms=3, color=INK)
    b.plot(t[winter], np.zeros(winter.sum()), "o", ms=3, color=MUTED)
    for t1, t2, col, ls, lw in ((4, 6, "#1a7f37", "-", 1.6), (6, 8, "#1a7f37", "-", 1.6), (8, 10, "#1a7f37", "-", 1.6),
                                (12, 16, "#bf8700", "-", 1.2), (18, 30, "#cf222e", "--", 1.0), (22, 52, "#cf222e", "--", 1.0)):
        b.add_patch(Arc(((t1 + t2) / 2, 0), t2 - t1, (t2 - t1) * 0.18 + 0.4, theta1=0, theta2=180, color=col, ls=ls, lw=lw))
    b.text(7, 0.75, "consecutive shortest pairs:\nfollow the motion", ha="center", fontsize=6.5, color="#1a7f37")
    b.text(14, 1.0, "twice the revisit:\npartly", ha="center", fontsize=6.5, color="#bf8700")
    b.text(19, 1.55, "long pairs: decorrelated,\nambiguous — lose it", ha="right", fontsize=6.5, color="#cf222e")
    b.annotate("", xy=(48, -0.05), xytext=(32, -0.05), arrowprops=dict(arrowstyle="-|>", color="#8250df", lw=1.2))
    b.text(40, -0.55 + 0.3, "the chain loses its level across the winter (summed pair noise) —\nthe level step comes "
                            "from the weather instead", ha="center", va="top", fontsize=6.5, color="#8250df")
    b.text(0, 4.0, "b  Which pairs hold the motion over a year", fontsize=7.5, color=INK)
    fig.suptitle("F3 (draft) — what C-band sees on a floating mat, and which pairs keep it", fontsize=9, color=INK)
    fig.savefig(OUT / "F3_schematic.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def f4() -> None:
    """Pair phase against the laser at P6, by revisit and track (X-065, read on the corrected inputs)."""
    d = pd.read_csv(latest("field_p6_short_pairs_x065", "p6_pairs_2020_2024.csv"))
    d = d.dropna(subset=["s1", "laser_los"])
    groups = [("6 days", d.dt == 6), ("12 days", d.dt == 12), ("24 days and longer", d.dt >= 24)]
    fig, axs = plt.subplots(2, 3, figsize=(7.2, 4.8), constrained_layout=True, sharex=True, sharey=True)
    lim = float(np.nanpercentile(np.abs(np.r_[d.s1, d.laser_los]), 99)) * 1.1
    for r, track in enumerate(TRACKS):
        for c, (lab, m) in enumerate(groups):
            ax = axs[r, c]
            g = d[(d.track == track) & m]
            ax.plot([-lim, lim], [-lim, lim], color=GRID, lw=0.8)
            ax.axhline(0, color=GRID, lw=0.5), ax.axvline(0, color=GRID, lw=0.5)
            ax.scatter(g.laser_los, g.s1, s=8, color="#0969da" if track == "ascending" else "#bc4c00", lw=0)
            if len(g) >= 5:
                rr = np.corrcoef(g.laser_los, g.s1)[0, 1]
                ax.text(0.04, 0.92, f"r = {rr:.2f}   n = {len(g)}", transform=ax.transAxes, fontsize=6.5, color=INK)
            ax.set_xlim(-lim, lim), ax.set_ylim(-lim, lim)
            if r == 0:
                ax.set_title(lab, fontsize=8, color=INK)
            if c == 0:
                ax.set_ylabel(f"{track}\npair phase (mm, line of sight)")
            if r == 1:
                ax.set_xlabel("laser change (mm, line of sight)")
    fig.suptitle("F4 (draft) — only the shortest consecutive pairs follow the laser, on both tracks (P6, 2020–2024)",
                 fontsize=9, color=INK)
    fig.savefig(OUT / "F4_pairs_vs_laser.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def f6() -> None:
    """The systematic test: every route, grouped by what it varies, with its outcome against the laser."""
    import yaml
    spec = yaml.safe_load((HUB / "08_deliverables" / "article_plan" / "f6_routes.yaml").read_text())
    known = set(pd.read_csv(HUB / "08_deliverables" / "independent_review" / "01_tries" / "tries.csv")["id"])
    rows = []
    for g in spec["groups"]:
        rows.append(("group", g["name"], None, None))
        for r in g["rows"]:
            missing = set(r["tries"]) - known
            if missing:
                raise SystemExit(f"F6: {r['route']} cites unknown tries {missing}")
            rows.append(("row", r["route"], r["product"], r["outcome"]))
    fig, ax = plt.subplots(figsize=(7.2, 0.21 * len(rows) + 0.9), constrained_layout=True)
    ax.set_xlim(0, 10), ax.set_ylim(len(rows), -1), ax.axis("off")
    ax.text(0.0, -0.6, "route", fontsize=7, color=MUTED), ax.text(5.6, -0.6, "product", fontsize=7, color=MUTED)
    ax.text(7.3, -0.6, "against the laser at P6", fontsize=7, color=MUTED)
    for k, (kind, label, prod, out) in enumerate(rows):
        if kind == "group":
            ax.text(0.0, k + 0.1, label, fontsize=7.5, color=INK, weight="bold")
            continue
        ax.text(0.25, k + 0.1, label, fontsize=7, color=INK)
        ax.text(5.6, k + 0.1, PRODUCT[prod], fontsize=6.5, color=MUTED)
        col, word = OUTCOME[out]
        ax.add_patch(Rectangle((7.25, k - 0.32), 0.22, 0.64, color=col))
        ax.text(7.6, k + 0.1, word, fontsize=6.5, color=INK)
    ax.set_title("F6 (draft) — every other route tried loses the mat's motion; only the consecutive shortest pairs "
                 "keep it", fontsize=9, color=INK)
    fig.savefig(OUT / "F6_systematic_test.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def f2() -> None:
    """F2 main: the mat's level (laser) and its water at P6 through the design years, datum segments marked."""
    import yaml
    t = pd.read_csv(CUBE / "gold" / "truth_laser_hourly.csv.gz", parse_dates=["time_utc"])
    prot = yaml.safe_load((FT / "protocol.yaml").read_text())
    fig, (a, b) = plt.subplots(2, 1, figsize=(7.2, 4.6), sharex=True, constrained_layout=True)
    ok, snow, filled = t.laser_ok.astype(bool), t.laser_snow_72h.astype(bool), t.laser_filled.astype(bool)
    a.plot(t.time_utc, t.laser_surface_cm.where(ok), color=INK, lw=0.7, label="mat surface (laser)")
    a.plot(t.time_utc, t.laser_surface_raw_cm.where(snow), color=GRID, lw=0.7, label="snow period (not used)")
    a.plot(t.time_utc, t.laser_surface_raw_cm.where(filled), color="#8250df", lw=1, ls="--",
           label="filled in the delivery (not data)")
    for x0, x1 in prot["datum_breaks"]:
        for ax in (a, b):
            ax.axvspan(pd.Timestamp(x0, tz="UTC"), pd.Timestamp(x1, tz="UTC"), color="#fff8c5", lw=0)
    a.text(pd.Timestamp(prot["datum_breaks"][0][1], tz="UTC"), a.get_ylim()[1], " re-levelled: new datum",
           fontsize=6, color=MUTED, va="top")
    a.set_ylabel("mat surface (cm, + up)")
    a.legend(fontsize=6, frameon=True, loc="upper right", ncol=1)
    b.plot(t.time_utc, t.cr_raw, color="#0969da", lw=0.7)
    b.set_ylabel("water level at P6\n(cm, raw logger)")
    b.set_xlim(t.time_utc.min(), t.time_utc.max())
    fig.suptitle("F2 (draft) — the mat rises and falls with its water, by centimetres (P6, design years; levels compared "
                 "within a datum only)", fontsize=9, color=INK)
    fig.savefig(OUT / "F2_laser_water.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def tables() -> None:
    """The article's data tables (``article_plan/tables/``): per track, the consecutive pairs at P6 the chain uses and
    the acquisitions they join; the site's climate over the archive (ERA5 at the site, X-081)."""
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import p6_validation_x082 as X
    out = HUB / "08_deliverables" / "article_plan" / "tables"
    out.mkdir(parents=True, exist_ok=True)
    g = pd.read_csv(CUBE / "gold" / "plots_pair.csv.gz", low_memory=False,
                    usecols=["track", "grid", "plot", "t1", "t2", "platforms", "revisit_days"])
    g = g[(g["plot"] == "P6") & (g.grid == "40m")].copy()
    for c in ("t1", "t2"):
        g[c] = pd.to_datetime(g[c], format="ISO8601", utc=True)
    rows = []
    for track, d in g.groupby("track"):
        dd = d[X.consecutive(d) & ~d.platforms.astype(str).str.contains("S1D")]
        rows.append({"track": track, "first": dd.t1.min().date(), "last": dd.t2.max().date(),
                     "acquisitions": int(pd.concat([dd.t1, dd.t2]).dt.normalize().nunique()),
                     "consecutive_pairs": len(dd), "pairs_6_day": int((dd.revisit_days == 6).sum()),
                     "pairs_12_day": int((dd.revisit_days == 12).sum()),
                     "pairs_longer": int((dd.revisit_days > 12).sum())})
    pd.DataFrame(rows).to_csv(out / "data_by_track.csv", index=False)
    w = pd.read_csv(HUB / "08_deliverables" / "cube_analysis_x081" / "weather_by_year.csv")
    pd.DataFrame([{"years": f"{w.year.min()}–{w.year.max()}", "precip_mm_median": w.precip_mm.median(),
                   "precip_mm_min": w.precip_mm.min(), "precip_mm_max": w.precip_mm.max(),
                   "snow_share_median": w.snow_share.median(), "frozen_share_median": w.frozen_share.median()}]
                 ).to_csv(out / "climate.csv", index=False)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    st, masks, water = static()
    tabs, extra = {}, {}
    for track in TRACKS:
        tabs[track], extra[track] = pair_table(track, masks, water)
        tabs[track].to_csv(OUT / f"F5_pairs_{track}.csv", index=False)
    f5(tabs)
    f5b(st, masks, extra)
    f2()
    f2_inset(st, masks)
    f1(st, masks)
    f3()
    f4()
    f6()
    tables()
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()

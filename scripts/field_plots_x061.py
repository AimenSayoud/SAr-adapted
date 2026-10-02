"""X-061 — plot roles: P6–P9 on the floating mat (main), P1–P5 water-table reference.

The supervisor (2026-09-28): the plots are not equivalent — P1–P5 are water-table monitoring
locations, P6–P9 the floating-mat plots (P6 with the laser). Compare the water table with the
Sentinel-1 phase/apparent displacement and coherence **for P6–P9**, with the finest reliable
extraction and the boardwalk flagged, and keep P1–P5 as hydrological reference.

Per plot, per track: the consecutive Sentinel-1 pairs (X-059's observation unit) in two
independent periods — 2020–2021 (S1A+S1B, mostly 6 d) and 2022–2024 (S1A, 12 d) — with the
plot's own water-table change at the overpass times (censored stretches out), the phase of the
plot's 40 m cell minus the grassland median (mm LOS), and its coherence. Two extractions: the
plot's cell, and a boardwalk-free one (median of the nearest mat cells to the west and east
that X-060 found free of boardwalk). Statistics as X-059: correlations with circular-shift
p-values, coherence with the season and pair length regressed out, all and dry pairs (X-050).

At 40 m, P6–P9 are three radar cells (P8 and P9 share one) and P1–P5 five; the grouping is
descriptive, not a test between independent sites.

    INSAR_DRIVE_ROOT=../local/drive_pristine PYTHONPATH=src python scripts/field_plots_x061.py
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from insar_wetlands import field  # noqa: E402
from insar_wetlands import field_link as fl  # noqa: E402
from insar_wetlands.bootstrap import start  # noqa: E402
from insar_wetlands.inversion.isbas import PHASE_TO_MM  # noqa: E402
from insar_wetlands.paths import make_paths  # noqa: E402
from insar_wetlands.stack import list_pairs, load_layer  # noqa: E402

HUB = Path(__file__).resolve().parents[3]
DLV = HUB / "08_deliverables"
HOUR = {"ascending": "16:36", "descending": "05:09"}
MAT = ["P6", "P7", "P8", "P9"]
REF = ["P1", "P2", "P3", "P4", "P5"]
PERIODS = {"2020–2021": (2020, 2021), "2022–2024": (2022, 2024)}


def plot_cells(ctx, px: pd.DataFrame, bw: pd.DataFrame) -> dict:
    """Per plot: its 40 m cell, and the boardwalk-free lateral cells (nearest mat cell to the west and
    to the east with no boardwalk, up to 3 cells away)."""
    mat = ctx.zones["A"].values.astype(bool)
    flagged = {(int(r.row), int(r.col)) for r in bw[bw.res_m == 40].itertuples()}
    out = {}
    for p in px.itertuples():
        r, c = int(p.row), int(p.col)
        lat = []
        for s in (-1, 1):
            for k in (1, 2, 3):
                rc = (r, c + s * k)
                if mat[rc] and rc not in flagged:
                    lat.append(rc)
                    break
        out[p.plot] = {"cell": [(r, c)], "off-boardwalk": lat, "cell_flagged": (r, c) in flagged}
    return out


def pair_rows(ctx, cells: dict, wtd: pd.DataFrame, wet: pd.DataFrame) -> pd.DataFrame:
    C = ctx.zones["C"].values.astype(bool)
    cens = {p: field.censored_flag(wtd[p]) for p in cells}
    rows = []
    for track in HOUR:
        roots = [make_paths(cfg=ctx.cfg, track=track).cropped, HUB / f"05_code/local/s1_2020_2021/hyp3_cropped_{track}"]
        have = {p: root for root in roots if root.exists() for p in list_pairs(root)}
        dates = sorted({pd.Timestamp(d) for p in have for d in p.split("_")})
        dates = [d for d in dates if 2020 <= d.year <= 2024]
        for a, b in zip(dates[:-1], dates[1:]):
            pid = f"{a:%Y%m%d}_{b:%Y%m%d}"
            if pid not in have:
                continue
            u = load_layer(have[pid], "unw_phase", [pid]).isel(pair=0).values
            co = load_layer(have[pid], "corr", [pid]).isel(pair=0).values
            ref = np.nanmedian(u[C & np.isfinite(u)])
            ta, tb = (pd.Timestamp(f"{d.date()} {HOUR[track]}", tz="UTC") for d in (a, b))
            wa, wb = (wet.loc[(f"{d.date()}", track)] if (f"{d.date()}", track) in wet.index else None for d in (a, b))
            flags = {"wet_either": bool((wa is not None and wa.wet) or (wb is not None and wb.wet)),
                     "frozen_either": bool((wa is not None and wa.frozen) or (wb is not None and wb.frozen)),
                     "flags_known": wa is not None and wb is not None}
            for plot, cc in cells.items():
                censored = bool(cens[plot].asof(ta.tz_convert(None)) or cens[plot].asof(tb.tz_convert(None))) \
                    if cens[plot].index.tz is None else bool(cens[plot].asof(ta) or cens[plot].asof(tb))
                dw = np.nan if censored else field._interp_at(wtd[plot], tb) - field._interp_at(wtd[plot], ta)
                for ext, rcs in (("cell", cc["cell"]), ("off-boardwalk", cc["off-boardwalk"])):
                    if not rcs:
                        continue
                    rr, cl = np.array([x[0] for x in rcs]), np.array([x[1] for x in rcs])
                    rows.append({"track": track, "pair": pid, "t1": a, "t2": b, "dt_days": (b - a).days,
                                 "mid": a + (b - a) / 2, "plot": plot, "group": "mat P6–P9" if plot in MAT else "reference P1–P5",
                                 "extraction": ext, "n_cells": len(rcs),
                                 "s1_dlos_mm": float((np.nanmedian(u[rr, cl]) - ref) * PHASE_TO_MM),
                                 "coh": float(np.nanmedian(co[rr, cl])), "dwtd_cm": dw, **flags})
    d = pd.DataFrame(rows)
    d["period"] = np.where(d.t2.dt.year <= 2021, "2020–2021", "2022–2024")
    return d


def stats(g: pd.DataFrame) -> dict:
    g = g.sort_values("mid")
    W = g[g.s1_dlos_mm.notna() & g.dwtd_cm.notna()]
    out = {"n_pairs": len(g), "n_wtd": len(W), "median_dt": float(g.dt_days.median()),
           "median_coh": float(g.coh.median()), "sd_s1_mm": float(g.s1_dlos_mm.std()),
           "sd_dwtd_cm": float(W.dwtd_cm.std()) if len(W) else np.nan}
    if len(W) >= 10:
        r, p = fl.circular_shift_p(W.dwtd_cm.to_numpy(), W.s1_dlos_mm.to_numpy())
        out.update({"r_s1_vs_dwtd": r, "p_s1_vs_dwtd": p,
                    "slope_s1_mm_per_cm": fl.slopes(W.dwtd_cm, W.s1_dlos_mm)["slope_ols"]})
        Wc = W[W.coh.notna()]
        res = fl.season_residual(Wc.coh.to_numpy(), Wc.mid, Wc.dt_days.to_numpy())
        rc, pc = fl.circular_shift_p(np.abs(Wc.dwtd_cm.to_numpy()), res)
        out.update({"r_coh_vs_abs_dwtd": rc, "p_coh_vs_abs_dwtd": pc})
    return out


def table(d: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (track, period, plot, ext), g in d.groupby(["track", "period", "plot", "extraction"]):
        for subset, sel in (("all", np.ones(len(g), bool)),
                            ("dry", (g.flags_known & ~g.wet_either & ~g.frozen_either).to_numpy())):
            rows.append({"track": track, "period": period, "plot": plot, "group": g.group.iloc[0], "extraction": ext,
                         "subset": subset, **stats(g[sel])})
    return pd.DataFrame(rows)


def group_summary(t: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (track, period, ext, subset, group), g in t.groupby(["track", "period", "extraction", "subset", "group"]):
        rows.append({"track": track, "period": period, "extraction": ext, "subset": subset, "group": group,
                     "plots": len(g), "median_r_s1_vs_dwtd": g.r_s1_vs_dwtd.median(),
                     "n_p05_s1_vs_dwtd": int((g.p_s1_vs_dwtd < 0.05).sum()),
                     "median_slope_mm_per_cm": g.slope_s1_mm_per_cm.median(),
                     "median_r_coh_vs_abs_dwtd": g.r_coh_vs_abs_dwtd.median(),
                     "n_p05_coh_vs_abs_dwtd": int((g.p_coh_vs_abs_dwtd < 0.05).sum()),
                     "median_coh": g.median_coh.median(), "median_sd_dwtd_cm": g.sd_dwtd_cm.median()})
    return pd.DataFrame(rows)


def fig(out: Path, t: pd.DataFrame):
    order = REF + MAT
    fig, axes = plt.subplots(2, 2, figsize=(13, 8), sharex=True)
    for j, (stat, lab) in enumerate((("r_s1_vs_dwtd", "r, phase vs ΔWTD"), ("r_coh_vs_abs_dwtd", "r, coherence vs |ΔWTD| (season removed)"))):
        for i, track in enumerate(HOUR):
            a = axes[i, j]
            for k, (period, mk) in enumerate((("2020–2021", "s"), ("2022–2024", "o"))):
                g = t[(t.track == track) & (t.period == period) & (t.extraction == "cell") & (t.subset == "all")].set_index("plot").reindex(order)
                pcol = np.where(g[stat.replace("r_", "p_")] < 0.05, "k", "none")
                cols = ["tab:green" if p in MAT else "tab:gray" for p in order]
                a.scatter(np.arange(len(order)) + (k - 0.5) * 0.25, g[stat], c=cols, marker=mk, edgecolors=pcol, s=45,
                          label=f"{period} (black edge: p < 0.05)")
            a.axhline(0, color="0.6", lw=0.8)
            a.axvspan(len(REF) - 0.5, len(order) - 0.5, color="tab:green", alpha=0.06)
            a.set_title(f"{track} — {lab}", fontsize=9)
            a.set_xticks(range(len(order)), order)
    axes[0, 0].legend(fontsize=7)
    fig.suptitle("Per plot, consecutive pairs, plot's own 40 m cell: reference P1–P5 (grey) · floating mat P6–P9 (green; P8 and P9 share a cell)", fontsize=10)
    fig.tight_layout()
    fig.savefig(out / "fig_plots_p1_p9.png", dpi=110)
    plt.close(fig)


def fp(p):
    return "≤ 0.001" if p <= 0.001 else f"{p:.3f}"


def md(df: pd.DataFrame) -> str:
    lines = ["| " + " | ".join(df.columns) + " |", "|" + "---|" * len(df.columns)]
    for r in df.itertuples(index=False):
        lines.append("| " + " | ".join(f"{v:.2f}" if isinstance(v, float) else str(v) for v in r) + " |")
    return "\n".join(lines)


def readme(out: Path, d, t, gs, cells):
    cl = pd.DataFrame([{"plot": p, "group": "mat" if p in MAT else "reference", "cell (row, col)": str(c["cell"][0]),
                        "cell has boardwalk": c["cell_flagged"], "off-boardwalk cells": ", ".join(map(str, c["off-boardwalk"]))}
                       for p, c in cells.items()])
    tt = t[(t.extraction == "cell") & (t.subset == "all")].copy()
    tt["phase vs ΔWTD r (p)"] = [f"{r:.2f} ({fp(p)})" for r, p in zip(tt.r_s1_vs_dwtd, tt.p_s1_vs_dwtd)]
    tt["coh vs |ΔWTD| r (p)"] = [f"{r:.2f} ({fp(p)})" for r, p in zip(tt.r_coh_vs_abs_dwtd, tt.p_coh_vs_abs_dwtd)]
    tt = tt[["track", "period", "plot", "group", "n_wtd", "median_dt", "median_coh", "sd_dwtd_cm", "slope_s1_mm_per_cm",
             "phase vs ΔWTD r (p)", "coh vs |ΔWTD| r (p)"]]
    g = gs[gs.subset == "all"].drop(columns="subset")
    find = []
    for track in HOUR:
        for period in PERIODS:
            s = gs[(gs.track == track) & (gs.period == period) & (gs.extraction == "cell") & (gs.subset == "all")].set_index("group")
            m, r = s.loc["mat P6–P9"], s.loc["reference P1–P5"]
            find.append(f"- **{track}, {period}**: phase vs ΔWTD significant at {m.n_p05_s1_vs_dwtd}/{m.plots} mat plots "
                        f"(median r {m.median_r_s1_vs_dwtd:.2f}) and {r.n_p05_s1_vs_dwtd}/{r.plots} reference plots "
                        f"({r.median_r_s1_vs_dwtd:.2f}); coherence vs |ΔWTD| at {m.n_p05_coh_vs_abs_dwtd}/{m.plots} mat "
                        f"(median r {m.median_r_coh_vs_abs_dwtd:.2f}) and {r.n_p05_coh_vs_abs_dwtd}/{r.plots} reference "
                        f"({r.median_r_coh_vs_abs_dwtd:.2f}); water-table change SD {m.median_sd_dwtd_cm:.1f} cm on the mat "
                        f"vs {r.median_sd_dwtd_cm:.1f} cm at the reference plots.")
    for track in HOUR:
        for period in PERIODS:
            s = gs[(gs.track == track) & (gs.period == period) & (gs.subset == "all")].set_index(["group", "extraction"])
            m_c, m_o = s.loc[("mat P6–P9", "cell")], s.loc[("mat P6–P9", "off-boardwalk")]
            r_c, r_o = s.loc[("reference P1–P5", "cell")], s.loc[("reference P1–P5", "off-boardwalk")]
            find.append(f"- *Boardwalk-free, {track} {period}*: phase vs ΔWTD significant at {m_o.n_p05_s1_vs_dwtd}/4 mat "
                        f"plots (cell: {m_c.n_p05_s1_vs_dwtd}/4; median r {m_o.median_r_s1_vs_dwtd:.2f} vs "
                        f"{m_c.median_r_s1_vs_dwtd:.2f}) and {r_o.n_p05_s1_vs_dwtd}/5 reference (cell: {r_c.n_p05_s1_vs_dwtd}/5).")
    cmp_ = t[t.subset == "all"].pivot_table(index=["track", "period", "plot"], columns="extraction",
                                            values=["r_s1_vs_dwtd", "r_coh_vs_abs_dwtd"]).dropna()
    dif = {k: float((cmp_[(k, "off-boardwalk")] - cmp_[(k, "cell")]).abs().max()) for k in ("r_s1_vs_dwtd", "r_coh_vs_abs_dwtd")}
    L = ["# X-061 — P6–P9 on the floating mat, P1–P5 as water-table reference (exploratory)", "",
         "Generated by `05_code/SAr-adapted/scripts/field_plots_x061.py` (outputs here only; field data). The supervisor's "
         "request of 2026-09-28: the main analysis on the mat plots P6–P9, the others as hydrological reference; finest "
         "extraction, boardwalk flagged (X-060). Every number is computed by the script.", "",
         "## Findings", "", *find,
         f"- Boardwalk-free extraction (lateral cells) against the plot's cell: largest change in r over all plot × track × "
         f"period cases {dif['r_s1_vs_dwtd']:.2f} (phase vs ΔWTD) and {dif['r_coh_vs_abs_dwtd']:.2f} (coherence vs |ΔWTD|) — "
         "full comparison in `plots_stats.csv`.", "",
         "## Settings", "",
         "- Observations: consecutive pairs per track; 2020–2021 (S1A + S1B, mostly 6 d) and 2022–2024 (S1A, 12 d) are "
         "independent periods. Phase: plot cell (or lateral cells) minus the grassland C median, mm LOS (+ toward the "
         "satellite). Water table: the plot's own hourly series interpolated at the overpass; censored stretches "
         "(sensor floor) out. Coherence: season and pair length regressed out (X-059). p from circular shifts in time "
         "(floor 1/2001). Dry = neither date wet nor frozen (X-050).",
         "- Radar units: P8 and P9 share one 40 m cell, so the mat group is three radar cells with four water tables.",
         "- No laser outside P6: at P7–P9 the phase can only be related to the water table, not to measured motion.", "",
         "## Plots on the grid", "", md(cl), "",
         "## Per plot (plot cell, all pairs)", "", md(tt), "",
         "## By group (medians over plots; counts of p < 0.05)", "", md(g), "",
         "## How to read this", "",
         "- At P6 the 12-day phase follows the laser surface (X-059) and the surface follows the water table, so a phase–ΔWTD "
         "link on the mat is expected *if* the phase records motion; at the reference plots (edge, more water-table change, "
         "less floating) a link would point to moisture or to the grassland reference instead.",
         "- Correlations at different plots share the same dates and the same reference: they are not independent tests.",
         "- A link present in the plot's cell but gone in the boardwalk-free cells may come from the boardwalk itself (X-060: "
         "boardwalk cells slightly more coherent, dawn VH brighter at 10 m) or from the few metres that separate the cells; "
         "a link present in both is not a boardwalk artefact.", "",
         "## Files", "", "| File | Content |", "|---|---|",
         "| `plots_pairs.csv` | every consecutive pair × plot × extraction: phase, coherence, ΔWTD, flags |",
         "| `plots_stats.csv` | statistics per plot × track × period × extraction × subset |",
         "| `plots_groups.csv` | the group medians |",
         "| `fig_plots_p1_p9.png` | per-plot correlations, both tracks and periods |"]
    (out / "README.md").write_text("\n".join(L) + "\n")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default=str(DLV / "field_plots_x061"))
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    ctx = start("field_plots_x061", mount=False, git=False)
    px = pd.read_csv(DLV / "field_first" / "plot_pixels.csv")
    bw = pd.read_csv(DLV / "field_boardwalk_x060" / "boardwalk_cells.csv")
    cells = plot_cells(ctx, px, bw)
    wtd = field.load_wtd_hourly()
    wet = pd.read_csv(DLV / "field_dew_x050" / "wetness_at_overpasses_2020_2024.csv").set_index(["date", "track"])
    d = pair_rows(ctx, cells, wtd, wet)
    d.to_csv(out / "plots_pairs.csv", index=False)
    t = table(d)
    t.to_csv(out / "plots_stats.csv", index=False)
    gs = group_summary(t)
    gs.to_csv(out / "plots_groups.csv", index=False)
    fig(out, t)
    readme(out, d, t, gs, cells)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()

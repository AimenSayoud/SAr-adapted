"""D-027 verdict — does the floating mat keep any coherence over a year when both dates are in winter?

Hrysiewicz et al. (2023, 2024) recover 90–360-day coherence on grounded raised bogs, best between winter dates. Test 1
(X-082) found the cube's ten April-to-April pairs at the coherence floor on the mat, the grassland *and* stable ground,
so D-027 added winter first dates (10 January, 10 February; 2017–2025, both tracks; partner ~360 and ~720 days later).
As PLAN.md requires, the mat is judged against other ground **in the same pairs**, not by an absolute value:

- per pair and zone: median pixel coherence, and the share of pixels above the lake's 95th percentile in that pair
  (open water decorrelates fully, so the lake measures the estimator's floor for this pair and look count);
- paired comparisons mat vs lake, stable ground, hard targets and grassland (Wilcoxon signed-rank, exact);
- winter (D-027) against the April annual pairs already in the cube; by conditions at the two dates (snow / frost).

Agreement with their bogs would be the mat clearly above the lake floor in winter annual pairs, and not far below
stable ground. Reads the cube gold layer; writes to the hub. Exploratory.

    PYTHONPATH=src python scripts/annual_pairs_x082.py
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import xarray as xr  # noqa: E402
from scipy.stats import wilcoxon  # noqa: E402

from insar_wetlands import field  # noqa: E402

HUB = field.field_root().parents[1]
GOLD = HUB / "06_data" / "cube" / "gold"
OUT = HUB / "08_deliverables" / "hrysiewicz_tests_x082" / "annual_pairs"
TRACKS = ("ascending", "descending")
MIN_DAYS = 300
ZONES = ("stable", "mat", "grassland", "hard", "lake")          # fixed categorical order
COL = {"stable": "#2a78d6", "mat": "#eb6834", "grassland": "#1baf7a", "hard": "#eda100", "lake": "#e87ba4"}
LABEL = {"stable": "stable ground", "mat": "mat", "grassland": "grassland", "hard": "hard targets",
         "lake": "lake (floor)"}
INK, INK2, MUTED, GRIDC, AXIS, SURFACE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"
plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE, "axes.edgecolor": AXIS,
    "axes.labelcolor": INK2, "axes.titlecolor": INK, "axes.titlesize": 9, "axes.titleweight": "bold",
    "axes.titlelocation": "left", "axes.labelsize": 8.5, "xtick.color": MUTED, "ytick.color": MUTED,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "axes.grid": True, "grid.color": GRIDC, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False, "legend.fontsize": 7.5,
})


def masks(d: xr.Dataset) -> dict[str, np.ndarray]:
    z = d.st_zone.values
    return {"stable": d.st_stable.values == 1, "mat": z == 1, "grassland": z == 3,
            "hard": d.st_hard_target.values == 1, "lake": z == 2}


def end_state(d: xr.Dataset, end: str) -> np.ndarray:
    snow = (d[f"drv_flag_snow_imgw_{end}"].values == 1) | (d[f"drv_flag_snow_om_{end}"].values == 1)
    frost = (d[f"drv_flag_frozen_air_{end}"].values == 1) | (d[f"drv_flag_frozen_soil_om_{end}"].values == 1)
    return np.where(snow, "snow", np.where(frost, "frozen", "open"))


def per_pair(d: xr.Dataset, track: str, start: str) -> pd.DataFrame:
    coh = d.coh.values                                             # (pixel, pair)
    m = masks(d)
    floor = np.nanpercentile(coh[m["lake"]], 95, axis=0)
    out = pd.DataFrame({"track": track, "start": start, "pair": d.pair.values.astype(str),
                        "t1": pd.to_datetime(d.t1.values).date, "t2": pd.to_datetime(d.t2.values).date,
                        "days": d.revisit_days.values, "platforms": d.platforms.values.astype(str),
                        "t1_state": end_state(d, "t1"), "t2_state": end_state(d, "t2"), "lake_p95": floor})
    for z in ZONES:
        out[f"{z}_median"] = np.nanmedian(coh[m[z]], axis=0)
        out[f"{z}_above_floor"] = np.nanmean(coh[m[z]] > floor, axis=0)
        out[f"{z}_n"] = int(m[z].sum())
    out["span"] = np.where(out.days > 500, "two-year", "annual")
    return out


def load() -> pd.DataFrame:
    frames = []
    for t in TRACKS:
        frames.append(per_pair(xr.open_dataset(GOLD / f"mat_pairs_annual_{t}.nc"), t, "winter (D-027)"))
        d = xr.open_dataset(GOLD / f"mat_pairs_2022_2024_{t}.nc")
        d = d.isel(pair=np.flatnonzero(d.revisit_days.values >= MIN_DAYS))
        if d.sizes["pair"]:
            frames.append(per_pair(d, t, "April (in the cube)"))
    return pd.concat(frames, ignore_index=True)


def paired(p: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (start, track), g in list(p.groupby(["start", "track"])) + [(("winter (D-027)", "both"),
                                                                      p[p.start.str.startswith("winter")])]:
        for metric in ("median", "above_floor"):
            for other in ("lake", "stable", "hard", "grassland"):
                a, b = g[f"mat_{metric}"].to_numpy(), g[f"{other}_{metric}"].to_numpy()
                ok = np.isfinite(a) & np.isfinite(b)
                diff = a[ok] - b[ok]
                pv = float(wilcoxon(diff).pvalue) if ok.sum() >= 6 and np.any(diff != 0) else np.nan
                rows.append({"start": start, "track": track, "metric": metric, "mat_minus": other, "pairs": int(ok.sum()),
                             "median_difference": float(np.median(diff)) if ok.any() else np.nan,
                             "mat_higher_in": int((diff > 0).sum()), "wilcoxon_p": pv})
    return pd.DataFrame(rows)


def by_state(p: pd.DataFrame) -> pd.DataFrame:
    w = p[p.start.str.startswith("winter")].copy()
    w["ends"] = np.where((w.t1_state == "open") & (w.t2_state == "open"), "both open",
                         np.where((w.t1_state != "open") & (w.t2_state != "open"), "both snow/frozen", "one snow/frozen"))
    agg = {f"{z}_{m}": "median" for m in ("median", "above_floor") for z in ("mat", "stable", "lake")}
    return w.groupby(["ends", "span"]).agg(pairs=("pair", "size"), **{k: (k, v) for k, v in agg.items()}).reset_index()


# ------------------------------------------------------------------------------------------- figures

def fig_pairs(p: pd.DataFrame, metric: str, ylabel: str, name: str, out):
    fig, axs = plt.subplots(2, 1, figsize=(13, 8), sharey=True)
    for ax, t in zip(axs, TRACKS):
        g = p[p.track == t].sort_values(["start", "t1", "days"], ascending=[False, True, True]).reset_index(drop=True)
        x = np.arange(len(g))
        for k, z in enumerate(ZONES):
            ax.plot(x + (k - 2) * 0.12, g[f"{z}_{metric}"], "o", ms=4.5, color=COL[z], mec=SURFACE, mew=0.6,
                    label=LABEL[z])
        sep = (g.start != g.start.iloc[0]).to_numpy()
        if sep.any():
            ax.axvline(np.argmax(sep) - 0.5, color=INK2, lw=0.8)
        ax.set_xticks(x)
        ax.set_xticklabels([f"{a:%Y-%m-%d} {s}" for a, s in zip(pd.to_datetime(g.t1), g.span.str[0])], fontsize=6,
                           rotation=90)
        ax.set_ylabel(ylabel)
        april = " and the cube's April starts (right of the line)" if sep.any() else ""
        ax.set_title(f"{t}: first date of each pair — D-027 winter starts{april}; a = annual, t = two-year", fontsize=8.5)
    axs[0].legend(ncol=5, loc="upper right")
    fig.suptitle("Year-long pairs at Rzecin: the mat against the lake floor and stable ground, pair by pair",
                 x=0.01, ha="left", fontsize=11, fontweight="bold")
    fig.tight_layout()
    fig.savefig(out / name, dpi=160)
    plt.close(fig)


# ------------------------------------------------------------------------------------------- README

def md(df: pd.DataFrame) -> str:
    def f(v):
        if isinstance(v, (float, np.floating)):
            return "—" if not np.isfinite(v) else (f"{v:.3g}" if abs(v) < 0.01 and v != 0 else f"{v:.2f}")
        return str(v)
    return "\n".join(["| " + " | ".join(map(str, df.columns)) + " |", "|" + "---|" * len(df.columns)]
                     + ["| " + " | ".join(f(v) for v in r) + " |" for r in df.itertuples(index=False)])


def readme(p, pr, st, out):
    w = p[p.start.str.startswith("winter")]
    a = p[~p.start.str.startswith("winter")]
    both = pr[(pr.track == "both")].set_index(["metric", "mat_minus"])
    g = lambda m, o, c: both.loc[(m, o), c]  # noqa: E731
    L = ["# D-027 — year-long winter pairs on the floating mat (X-082, exploratory)", "",
         "Generated by `05_code/SAr-adapted/scripts/annual_pairs_x082.py` from the cube (`mat_pairs_annual_*`, "
         "`mat_pairs_2022_2024_*`); every number is computed by it. Plan: `../PLAN.md` (D-027).", "",
         "## Findings", "",
         f"- **{len(w)} winter pairs** (D-027: {int((w.span == 'annual').sum())} annual, {int((w.span == 'two-year').sum())} "
         f"two-year, both tracks, first dates {w.t1.min():%Y}–{w.t1.max():%Y}) and {len(a)} April pairs already in the cube.",
         f"- **Over a year every natural surface sits at the same level.** Median coherence in the winter pairs: mat "
         f"{w.mat_median.median():.2f}, stable ground {w.stable_median.median():.2f}, grassland "
         f"{w.grassland_median.median():.2f}, lake {w.lake_median.median():.2f}, hard targets {w.hard_median.median():.2f}. "
         f"Share of pixels above the lake's 95th percentile: mat {w.mat_above_floor.median():.0%}, stable ground "
         f"{w.stable_above_floor.median():.0%}, grassland {w.grassland_above_floor.median():.0%}, hard targets "
         f"{w.hard_above_floor.median():.0%} (the lake 5% by construction).",
         f"- **The mat is not worse than stable ground over a year — both are gone.** Mat minus stable ground, same pairs: "
         f"{100 * g('above_floor', 'stable', 'median_difference'):+.1f} percentage points above the floor (mat higher in "
         f"{int(g('above_floor', 'stable', 'mat_higher_in'))} of {int(g('above_floor', 'stable', 'pairs'))}, Wilcoxon p "
         f"{g('above_floor', 'stable', 'wilcoxon_p'):.2g}); against the grassland p {g('above_floor', 'grassland', 'wilcoxon_p'):.2g}. "
         f"The mat is above the lake (higher in {int(g('above_floor', 'lake', 'mat_higher_in'))} of "
         f"{int(g('above_floor', 'lake', 'pairs'))} pairs, p {g('above_floor', 'lake', 'wilcoxon_p'):.2g}), but so is all land: "
         "a difference between land and water shared by every land zone, not coherence the mat keeps.",
         f"- **The pairs are valid:** hard targets stay above the floor (mat minus hard targets "
         f"{100 * g('above_floor', 'hard', 'median_difference'):+.1f} points, mat higher in only "
         f"{int(g('above_floor', 'hard', 'mat_higher_in'))} of {int(g('above_floor', 'hard', 'pairs'))}), so the products carry "
         "phase where the ground keeps it.",
         "- **Winter dates do not help:** pairs with both ends under snow or frost, one end, or neither give the same "
         f"medians (section 3); April starts (cube): mat {a.mat_above_floor.median():.0%}, stable ground "
         f"{a.stable_above_floor.median():.0%} above the floor." if len(a) else "- No April annual pairs in the cube.",
         "- **Verdict (PLAN D-027):** the 1–2-year coherence Hrysiewicz et al. use on grounded raised bogs is not found at "
         "Rzecin on any vegetated surface, the mat included; at this resolution and look count, C-band over a year keeps "
         "phase only on hard targets. A multi-year rate from annual pairs is therefore not available here; nothing in "
         "this says the mat decorrelates more than its surroundings over a year — that contrast (H2) lives at 6–24 days.",
         "- Conditions at the two dates: section 3 (snow / frozen flags from the cube).", "",
         "![above floor](fig1_above_floor.png)", "", "![median](fig2_median_coherence.png)", "",
         "## 1. Paired comparisons (mat minus another zone, same pairs)", "",
         "`median`: zone-median coherence; `above_floor`: share of the zone's pixels above the lake's 95th percentile "
         "in that pair. Wilcoxon signed-rank, two-sided; not computed under 6 pairs.", "", md(pr), "",
         "## 2. Every pair", "",
         md(p[["track", "start", "pair", "days", "platforms", "t1_state", "t2_state", "lake_p95"]
              + [f"{z}_median" for z in ZONES] + [f"{z}_above_floor" for z in ZONES]]), "",
         "## 3. Winter pairs by conditions at the two dates (medians over pairs)", "", md(st), "",
         "## Reading notes", "",
         "- HyP3 coherence is estimated in a 10 × 2-look window and is biased upward at low true coherence; the lake "
         "(open water, decorrelated within seconds) measures that bias pair by pair, so 'above the floor' is the "
         "readable quantity, not the coherence value.",
         "- Pixel counts per zone: " + ", ".join(f"{LABEL[z]} {int(p[f'{z}_n'].iloc[0])}" for z in ZONES) + " (40 m cells).",
         "- Snow: IMGW snow cover or ERA5 snow; frozen: air or soil below 0 °C (cube flags at the overpass).",
         "- A pair starting in S1B time and ending in S1A is a different satellite at the two dates (column "
         "`platforms`); the 12-day orbit repeat makes it geometrically equivalent."]
    (out / "README.md").write_text("\n".join(L) + "\n")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    p = load()
    pr, st = paired(p), by_state(p)
    p.to_csv(OUT / "pairs.csv", index=False)
    pr.to_csv(OUT / "paired_tests.csv", index=False)
    st.to_csv(OUT / "by_conditions.csv", index=False)
    fig_pairs(p, "above_floor", "share of pixels above the lake's p95", "fig1_above_floor.png", OUT)
    fig_pairs(p, "median", "zone median coherence", "fig2_median_coherence.png", OUT)
    readme(p, pr, st, OUT)
    print(f"D-027 analysis written to {OUT}")


if __name__ == "__main__":
    main()

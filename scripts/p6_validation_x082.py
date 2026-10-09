"""X-082 test 2 — the P6 phase against the laser, scored the way Hrysiewicz et al. (2024) score InSAR against peat
cameras, so the two can be compared. Exploratory; plan in ``08_deliverables/hrysiewicz_tests_x082/PLAN.md``.

One InSAR chain per track: every consecutive pair at the P6 cell from the first laser date (2021) to the last (2025),
phase referenced to the grassland (X-065) or to hard targets (WP0), cumulated. The laser is a *level* series (QC'd,
X-058), converted to LOS with P6's incidence, so it only needs a value at the epochs scored — snow gaps do not break
the chain. Views, as their §3.3:

- **full span** (their recipe): one linear trend removed over 2021–2025, then Pearson r and RMSE of the detrended
  series (their Table 2), deviations by season and shares within ±5 / ±10 mm (their Fig. 11e–f);
- **per season** (May–November of each year, restarted): the same statistics one season at a time;
- **increments**: Pearson r per pair (their Table 3, right column; the X-065 statistic);
- the **fusion v3** estimator at P6 (WP5 ``series_p6.csv``): detrended, and on demeaned levels as WP5 scored it.

Multi-year rates are not compared: the laser covers 2021–2025 with winter gaps (their own caution — less than 3 years
of in-situ data is not enough); the 2017–2026 laser archive will allow it. Their published values are read from
``02_literature/notes/hrysiewicz2024_table3.csv`` (a citation table). Reads the cube and the hub; writes only to the
hub. No field value in this file.

    PYTHONPATH=src python scripts/p6_validation_x082.py
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import xarray as xr  # noqa: E402

from insar_wetlands import field  # noqa: E402
from insar_wetlands.inversion.isbas import PHASE_TO_MM  # noqa: E402

HUB = field.field_root().parents[1]
CUBE = HUB / "06_data" / "cube"
OUT = HUB / "08_deliverables" / "hrysiewicz_tests_x082" / "p6_validation"
THEIRS = HUB / "02_literature" / "notes" / "hrysiewicz2024_table3.csv"
V3 = HUB / "08_deliverables" / "fusion_v3" / "wp5" / "series_p6.csv"
ASF = HUB / "06_data" / "s1_archive" / "asf_acquisitions.csv"
TRACKS = ("ascending", "descending")
REFS = {"grassland": "ref_grassland_unw", "hard_targets": "ref_hard_unw"}
MIN_EPOCHS = 8
# The laser was re-levelled across its July 2022 gap (laser QC, X-058 §5, §7: the filled stretch moves several times
# what the water predicts and the grid offset changes): levels before and after are on different datums.
RELEVEL = (pd.Timestamp("2022-07-05", tz="UTC"), pd.Timestamp("2022-08-08", tz="UTC"))
DATUMS = {"A": "laser datum A 2021-05 → 2022-07", "B": "laser datum B 2022-08 → 2025-12"}
SEASON = {12: "winter", 1: "winter", 2: "winter", 3: "spring", 4: "spring", 5: "spring", 6: "summer", 7: "summer",
          8: "summer", 9: "autumn", 10: "autumn", 11: "autumn"}

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, MUTED, GRIDC, AXIS, SURFACE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"
SCOL = {"winter": BLUE, "spring": AQUA, "summer": ORANGE, "autumn": "#4a3aa7"}
plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE, "axes.edgecolor": AXIS,
    "axes.labelcolor": INK2, "axes.titlecolor": INK, "axes.titlesize": 9, "axes.titleweight": "bold",
    "axes.titlelocation": "left", "axes.labelsize": 8.5, "xtick.color": MUTED, "ytick.color": MUTED,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "axes.grid": True, "grid.color": GRIDC, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False, "legend.fontsize": 7.5,
})


# ------------------------------------------------------------------------------------------- data

def consecutive(p: pd.DataFrame) -> pd.Series:
    """True where no acquisition of the same track lies strictly between the pair's dates (C-059 list)."""
    acq = pd.read_csv(ASF, parse_dates=["date"])
    out = pd.Series(False, index=p.index)
    for t in TRACKS:
        d = np.sort(acq[acq.track == t].date.values.astype("datetime64[D]"))
        s = p.track == t
        a = p.loc[s, "t1"].dt.tz_convert(None).values.astype("datetime64[D]")
        b = p.loc[s, "t2"].dt.tz_convert(None).values.astype("datetime64[D]")
        out[s] = np.searchsorted(d, b) - np.searchsorted(d, a, side="right") == 0
    return out


def load_pairs() -> pd.DataFrame:
    cols = ["track", "grid", "plot", "pair", "t1", "t2", "revisit_days", "unw_1x1", "d_laser_los_mm"] + list(REFS.values())
    p = pd.read_csv(CUBE / "gold" / "plots_pair.csv.gz", usecols=cols, low_memory=False)
    p = p[(p["plot"] == "P6") & (p.grid == "40m")].copy()
    for c in ("t1", "t2"):
        p[c] = pd.to_datetime(p[c], format="ISO8601", utc=True)
    p = p[consecutive(p)].sort_values(["track", "t1"])
    for name, col in REFS.items():
        p[f"s1_{name}"] = (p.unw_1x1 - p[col]) * PHASE_TO_MM
    return p.reset_index(drop=True)


def laser() -> tuple[dict[str, pd.Series], pd.Series]:
    """QC'd laser surface (mm, vertical) at every hour, and as LOS for each track (P6's incidence)."""
    site = pd.read_csv(CUBE / "silver" / "site_hourly.csv.gz", index_col=0, low_memory=False,
                       usecols=["time_utc", "laser_surface_cm", "laser_ok"])
    site.index = pd.to_datetime(site.index, format="ISO8601", utc=True)
    ok = site.laser_ok.map({True: True, False: False, "True": True, "False": False}).fillna(False).to_numpy()
    lev = site.laser_surface_cm.where(ok) * 10.0
    st = xr.open_dataset(CUBE / "silver" / "static_40m.nc")
    px = pd.read_csv(HUB / "08_deliverables" / "field_first" / "plot_pixels.csv").set_index("plot")
    r, c = int(px.loc["P6", "row"]), int(px.loc["P6", "col"])
    inc = {"ascending": float(st.incidence_asc_deg.values[r, c]), "descending": float(st.incidence_des_deg.values[r, c])}
    return {t: lev * np.cos(np.radians(inc[t])) for t in TRACKS}, lev


def chain(g: pd.DataFrame, ref: str) -> pd.DataFrame:
    """Cumulative phase along consecutive pairs; a missing pair starts a new segment (``seg``)."""
    rows, seg, level, prev = [], 0, 0.0, None
    for r in g.itertuples():
        if prev is None or r.t1 != prev:
            if prev is not None:
                seg += 1
            level = 0.0
            rows.append({"time": r.t1, "seg": seg, "insar": 0.0})
        level += getattr(r, f"s1_{ref}")
        rows.append({"time": r.t2, "seg": seg, "insar": level})
        prev = r.t2
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------------------------- scoring

def detrend(t_days: np.ndarray, y: np.ndarray) -> np.ndarray:
    return y - np.polyval(np.polyfit(t_days, y, 1), t_days)


def stats(d: pd.DataFrame) -> tuple[dict, np.ndarray, np.ndarray]:
    tt = (d.time - d.time.iloc[0]).dt.total_seconds().to_numpy() / 86400
    dl, ds = detrend(tt, d.laser.to_numpy(float)), detrend(tt, d.insar.to_numpy(float))
    dev = ds - dl
    return ({"epochs": len(d), "pearson_detrended": float(np.corrcoef(dl, ds)[0, 1]),
             "rmse_detrended_mm": float(np.sqrt(np.mean(dev ** 2))),
             "laser_sd_mm": float(np.std(dl)), "within_5mm": float(np.mean(np.abs(dev) < 5)),
             "within_10mm": float(np.mean(np.abs(dev) < 10)), "amplitude_ratio": float(np.std(ds) / np.std(dl))}, dl, ds)


def score(p: pd.DataFrame):
    las, _ = laser()
    rows, epochs, inc, chains, winters = [], [], [], [], []
    for t in TRACKS:
        lt = las[t].dropna()
        g = p[(p.track == t) & (p.t1 >= lt.index.min() - pd.Timedelta(days=13)) &
              (p.t2 <= lt.index.max() + pd.Timedelta(days=13))].sort_values("t1")
        for ref in REFS:
            c = chain(g, ref)
            c["laser"] = lt.reindex(c.time.dt.round("h")).to_numpy()
            chains.append(c.assign(track=t, reference=ref))
            c = c.dropna(subset=["laser"]).copy()
            for _, cs in c.groupby("seg"):              # each segment starts at its first scored epoch
                c.loc[cs.index, "laser"] -= cs.laser.iloc[0]
                c.loc[cs.index, "insar"] -= cs.insar.iloc[0]
            c["datum"] = np.where(c.time < RELEVEL[1], "A", "B")
            views = [("full span 2021–2025 (crosses the re-levelling)", c)]
            views += [(DATUMS[k], c[c.datum == k]) for k in DATUMS]
            views += [(f"season {y}{'' if y != 2022 else ' ' + k}", s_) for y in range(2021, 2026) for k in DATUMS
                      if len(s_ := c[(c.time.dt.year == y) & c.time.dt.month.between(5, 11) & (c.datum == k)])]
            winters += winter_steps(c, t, ref)
            for view, d in views:
                if len(d) < MIN_EPOCHS:
                    continue
                d = d.copy()
                if view.startswith("season"):
                    d["laser"] -= d.laser.iloc[0]
                    d["insar"] -= d.insar.iloc[0]
                st, dl, ds = stats(d)
                rows.append({"track": t, "view": view, "reference": ref, "segments": d.seg.nunique(), **st})
                if view == DATUMS["B"]:
                    epochs += [{"track": t, "reference": ref, "time": w, "season": SEASON[w.month], "insar_mm": a,
                                "laser_mm": b, "deviation_mm": a - b} for w, a, b in zip(d.time, ds, dl)]
            gi = g.dropna(subset=["d_laser_los_mm"])
            inc.append({"track": t, "reference": ref, "pairs": len(gi),
                        "pearson_increments": float(np.corrcoef(gi.d_laser_los_mm, gi[f"s1_{ref}"])[0, 1])})
    return (pd.DataFrame(rows), pd.DataFrame(epochs), pd.DataFrame(inc), pd.concat(chains, ignore_index=True),
            pd.DataFrame(winters))


def winter_steps(c: pd.DataFrame, track: str, ref: str) -> list[dict]:
    """Across each winter: last scored epoch on or before 30 November → first on or after 1 March, same datum."""
    out = []
    for y in range(2021, 2025):
        a = c[(c.time.dt.year == y) & c.time.dt.month.between(5, 11)]
        b = c[(c.time >= pd.Timestamp(f"{y + 1}-03-01", tz="UTC")) & (c.time < pd.Timestamp(f"{y + 1}-07-01", tz="UTC"))]
        if not len(a) or not len(b):
            continue
        a, b = a.iloc[-1], b.iloc[0]
        if a.datum != b.datum:
            continue
        out.append({"track": track, "reference": ref, "winter": f"{y}/{str(y + 1)[2:]}", "from": a.time.date(),
                    "to": b.time.date(), "laser_change_mm": b.laser - a.laser, "insar_change_mm": b.insar - a.insar,
                    "difference_mm": (b.insar - a.insar) - (b.laser - a.laser)})
    return out


def score_v3() -> pd.DataFrame:
    """Fusion v3 at P6 against the laser (vertical): detrended (their statistic) and demeaned levels (WP5's)."""
    if not V3.exists():
        return pd.DataFrame()
    v = pd.read_csv(V3)
    v["time_utc"] = pd.to_datetime(v.time_utc, format="ISO8601", utc=True)
    _, lev = laser()
    out = []
    for (run, variant), g in v.groupby(["run", "variant"]):
        g = g.sort_values("time_utc").drop_duplicates("time_utc")
        lv = lev.reindex(g.time_utc.dt.round("h")).to_numpy()
        m = np.isfinite(lv) & np.isfinite(g.h_v3_mm.to_numpy())
        if m.sum() < MIN_EPOCHS:
            continue
        d = pd.DataFrame({"time": g.time_utc[m].to_numpy(), "laser": lv[m], "insar": g.h_v3_mm.to_numpy()[m]})
        d["time"] = pd.to_datetime(d.time, utc=True)
        st, _, _ = stats(d)
        a, b = d.laser - d.laser.mean(), d.insar - d.insar.mean()
        out.append({"run": run, "variant": variant, **st, "pearson_levels": float(np.corrcoef(a, b)[0, 1]),
                    "rmse_levels_mm": float(np.sqrt(np.mean((b - a) ** 2)))})
    return pd.DataFrame(out)


# ------------------------------------------------------------------------------------------- figures

def fig_chain(ch: pd.DataFrame, out):
    las, _ = laser()
    fig, axs = plt.subplots(2, 1, figsize=(13, 7.5), sharex=True)
    for ax, t in zip(axs, TRACKS):
        lt = las[t].dropna()
        for ref, col in (("grassland", ORANGE), ("hard_targets", BLUE)):
            c = ch[(ch.track == t) & (ch.reference == ref)]
            for k, (_, cs) in enumerate(c.groupby("seg")):
                ax.plot(cs.time, cs.insar - cs.insar.iloc[0], "-", color=col, lw=1.0,
                        label=f"InSAR chain, ref. {ref.replace('_', ' ')}" if k == 0 else None)
        c0 = ch[(ch.track == t) & (ch.reference == "hard_targets")]
        l0 = lt.reindex(c0.time.dt.round("h"))
        first = l0.dropna().iloc[0] if l0.notna().any() else 0.0
        ax.plot(c0.time, l0.to_numpy() - first, "o", color=INK, ms=2.5, label="laser at the overpasses (LOS)")
        ax.axvspan(*RELEVEL, color=GRIDC, lw=0)
        ax.annotate("laser re-levelled\n(QC §5, §7)", (RELEVEL[1], 1), xycoords=("data", "axes fraction"),
                    xytext=(4, -4), textcoords="offset points", va="top", fontsize=7, color=INK2)
        ax.set_ylabel("mm LOS")
        ax.set_title(f"{t}: one chain of consecutive pairs 2021–2025 against the laser level (not detrended)", fontsize=8.5)
    axs[0].legend(loc="lower left")
    fig.suptitle("P6: the InSAR chain and the laser over five seasons", x=0.01, ha="left", fontsize=11, fontweight="bold")
    fig.tight_layout()
    fig.savefig(out / "fig1_chain_vs_laser.png", dpi=160)
    plt.close(fig)


def fig_scatter_box(ep: pd.DataFrame, out):
    e = ep[ep.reference == "hard_targets"]
    fig, axs = plt.subplots(1, 2, figsize=(12, 5.2))
    for s, col in SCOL.items():
        h = e[e.season == s]
        axs[0].scatter(h.laser_mm, h.insar_mm, s=10, color=col, alpha=0.7, edgecolors="none", label=f"{s} (n {len(h)})")
    lim = np.nanmax(np.abs(e[["laser_mm", "insar_mm"]].to_numpy())) * 1.05
    x = np.array([-lim, lim])
    axs[0].fill_between(x, x - 10, x + 10, color=GRIDC, alpha=0.6, lw=0)
    axs[0].fill_between(x, x - 5, x + 5, color=AXIS, alpha=0.5, lw=0)
    axs[0].plot(x, x, color=INK2, lw=1)
    axs[0].set_xlim(-lim, lim), axs[0].set_ylim(-lim, lim), axs[0].set_aspect("equal")
    axs[0].set_xlabel("laser, detrended (mm LOS)"), axs[0].set_ylabel("InSAR chain, detrended (mm LOS)")
    axs[0].set_title("their Fig. 11e at P6: laser datum B (2022-08 → 2025), both tracks, ref. hard targets", fontsize=8.5)
    axs[0].legend(loc="upper left")
    seas = list(SCOL)
    data = [e[e.season == s].deviation_mm.dropna().to_numpy() for s in seas]
    axs[1].axhspan(-10, 10, color=GRIDC, alpha=0.6, lw=0), axs[1].axhspan(-5, 5, color=AXIS, alpha=0.5, lw=0)
    b = axs[1].boxplot(data, patch_artist=True, manage_ticks=False, positions=range(4), medianprops=dict(color=INK),
                       flierprops=dict(markersize=3))
    for box, s in zip(b["boxes"], seas):
        box.set(facecolor=SCOL[s], alpha=0.55, edgecolor="none")
    axs[1].set_xticks(range(4)), axs[1].set_xticklabels([f"{s}\n(n {len(d)})" for s, d in zip(seas, data)])
    axs[1].axhline(0, color=INK2, lw=0.8)
    axs[1].set_ylabel("InSAR − laser (mm)"), axs[1].set_title("their Fig. 11f at P6: deviations by season (datum B)", fontsize=8.5)
    fig.tight_layout()
    fig.savefig(out / "fig2_scatter_and_seasons.png", dpi=160)
    plt.close(fig)


# ------------------------------------------------------------------------------------------- README

def md(df: pd.DataFrame, fmt="{:.2f}") -> str:
    f = lambda v: fmt.format(v) if isinstance(v, (float, np.floating)) and np.isfinite(v) else ("—" if isinstance(v, float) else str(v))  # noqa: E731
    return "\n".join(["| " + " | ".join(map(str, df.columns)) + " |", "|" + "---|" * len(df.columns)]
                     + ["| " + " | ".join(f(v) for v in r) + " |" for r in df.itertuples(index=False)])


def readme(sc, ep, inc, wi, v3, theirs, out):
    ts, tc = theirs[theirs.network == "short TBN"], theirs[theirs.network == "combined TBN"]
    h = sc[sc.reference == "hard_targets"]
    pick = lambda v: h[h.view == v].set_index("track")  # noqa: E731
    full, dB = pick(sc.view[sc.view.str.startswith("full")].iloc[0]), pick(DATUMS["B"])
    ses = h[h.view.str.startswith("season")]
    e = ep[ep.reference == "hard_targets"]
    hi = inc[inc.reference == "hard_targets"].set_index("track")
    w = wi[wi.reference == "hard_targets"]
    line = lambda d: "; ".join(f"{t} r {d.loc[t, 'pearson_detrended']:.2f}, RMSE {d.loc[t, 'rmse_detrended_mm']:.1f} mm, "  # noqa: E731
                               f"{d.loc[t, 'within_10mm']:.0%} within ±10 mm (laser SD {d.loc[t, 'laser_sd_mm']:.0f} mm)" for t in d.index)
    L = ["# X-082 test 2 — P6 against the laser, on Hrysiewicz et al. (2024)'s statistics (exploratory)", "",
         "Generated by `05_code/SAr-adapted/scripts/p6_validation_x082.py`; every number is computed by it, their values "
         "are read from `02_literature/notes/hrysiewicz2024_table3.csv`. Plan: `../PLAN.md`.", "",
         "## Findings", "",
         "- **The laser is not on one datum over 2021–2025.** It was re-levelled across its July 2022 gap "
         "(`08_deliverables/field_p6/laser_qc/laser_qc.md` §5, §7: the grid offset changes there) and drifts against both "
         "surface references. The size of the re-levelling is unknown — the grid offset fixes it only modulo the sensor "
         "resolution, and part of the drop across the gap is predicted by the water level (QC §5 table), so it is not all "
         "datum. Their recipe (one trend removed over the whole series) is therefore applied within each datum; the score "
         "across the re-levelling is kept for transparency only: " + line(full) + ". The 2022 season has fewer than "
         f"{MIN_EPOCHS} scored epochs on either side of the gap and is not scored on its own.",
         "- **Datum B, August 2022 → 2025, their recipe** (one chain of consecutive pairs, one trend removed, reference hard "
         "targets; the closest analogue of their 2.5-year series): " + line(dB) + ".",
         f"- **Their raised bogs:** short networks r {ts.pearson_detrended.min():.2f}–{ts.pearson_detrended.max():.2f} "
         f"(RMSE {ts.rmse_mm.min():.1f}–{ts.rmse_mm.max():.1f} mm), combined r {tc.pearson_detrended.min():.2f}–"
         f"{tc.pearson_detrended.max():.2f}.",
         f"- **One season at a time (May–November):** r {ses.pearson_detrended.min():.2f}–{ses.pearson_detrended.max():.2f} "
         f"across {len(ses)} track × season cases (median {ses.pearson_detrended.median():.2f}, RMSE median "
         f"{ses.rmse_detrended_mm.median():.1f} mm). Removing a trend from one season also removes most of its seasonal "
         "decline, so these values measure the shorter fluctuations — not comparable with their multi-year series.",
         f"- **Across the winters** (last scored epoch ≤ 30 November → first ≥ 1 March, same datum): chain minus laser "
         f"{w.difference_mm.min():+.0f} to {w.difference_mm.max():+.0f} mm over {len(w)} track × winter cases; "
         f"{int((w.difference_mm.abs() > 10).sum())} exceed ±10 mm ("
         + ", ".join(f"{r.track[:3]} {r.winter}" for r in w[w.difference_mm.abs() > 10].itertuples()) + ") — section 3. "
         "What the chain loses over a winter stays in every later epoch, which is why the multi-year score is lower than "
         "the per-season ones.",
         "- **Increments (their Table 3 right column; the X-065 statistic):** "
         + "; ".join(f"{t} r {hi.loc[t, 'pearson_increments']:.2f} over {int(hi.loc[t, 'pairs'])} pairs" for t in hi.index) + ".",
         "- **Fusion v3 at P6:** section 5 — scored both detrended (their statistic) and on levels (WP5's).", "",
         "![chain](fig1_chain_vs_laser.png)", "", "![scatter](fig2_scatter_and_seasons.png)", "",
         "## 1. Scores per track, view and reference", "", md(sc), "",
         "## 2. Increments", "", md(inc), "",
         "## 3. Across the winters", "", md(wi, "{:.1f}"), "",
         "## 4. Deviations by season (datum B, reference hard targets)", "",
         md(e.groupby("season").deviation_mm.agg(n="size", median="median",
                                                 within_5mm=lambda d: float(np.mean(np.abs(d) < 5)),
                                                 within_10mm=lambda d: float(np.mean(np.abs(d) < 10))).reset_index()), "",
         "## 5. Fusion v3 at P6 (vertical)", "", md(v3) if len(v3) else "—", "",
         "## 6. Their published values (citation table)", "", md(theirs), "",
         "## Reading notes", "",
         "- InSAR: P6 cell (40 m) minus the reference zone's median phase, `PHASE_TO_MM`, + towards the satellite, cumulated "
         "over consecutive pairs (6-day in 2021 and 2025, 12-day in 2022–2024); a cycle slip in any pair stays in the "
         "chain from then on. Laser: QC'd level (X-058) × cos(incidence at P6), at the overpass hour; epochs without a "
         "laser value (snow, gaps) are not scored but do not break the chain.",
         "- The laser drift (QC §7) is linear to first order and is removed by the detrending within a datum; the "
         "re-levelling is not, hence the split.",
         "- Their in-situ series are interpolated to SAR dates; ours are read at the overpass hour.",
         "- Multi-year rates are not compared (laser 2021–2025 with winter gaps and a re-levelling); the 2017–2026 laser "
         "archive, with its datum history from the field team (X-057), would allow it.",
         "- Their stations are grounded bogs, ~2.5 years of peat-camera data; ours is one floating-mat cell."]
    (out / "README.md").write_text("\n".join(L) + "\n")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    p = load_pairs()
    sc, ep, inc, ch, wi = score(p)
    v3 = score_v3()
    theirs = pd.read_csv(THEIRS)
    for name, df in (("scores", sc), ("epochs", ep), ("increments", inc), ("chains", ch), ("winters", wi), ("fusion_v3_scores", v3)):
        df.to_csv(OUT / f"{name}.csv", index=False)
    for stale in ("rates.csv", "fig1_cumulative_runs.png"):
        (OUT / stale).unlink(missing_ok=True)
    fig_chain(ch, OUT)
    fig_scatter_box(ep, OUT)
    readme(sc, ep, inc, wi, v3, theirs, OUT)
    print(f"X-082 test 2 written to {OUT}")


if __name__ == "__main__":
    main()

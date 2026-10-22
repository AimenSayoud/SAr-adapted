"""F7a — the estimator through the design years: the laser at P6, the plain chain of consecutive pairs and the
re-levelled short-pair chain (v4, frozen), per track, each winter shaded, with the anchor's spread as a band.

The band is X-107 amendment_3's error statement drawn: v4 rebuilt with each of the five SPAMS fits that leave one design
year out (``v4/spams_loyo_params.yaml``), everything else frozen; at each epoch the band spans the five rebuilt levels.
It is in-sample (the design years) and says how much the winter steps depend on the years the anchor is fitted on — not
how well v4 does on years it has not seen (Fig. 7b, after X-107). Levels are compared as the scorer compares them: each
series is set to the laser at the first scored epoch of each chain segment within each laser datum
(``scoring.level_views``), so nothing is compared across the 2022 re-levelling. The chains are built by
``freeze_v4_x106.build_v4`` with the frozen parameter file — the code X-107 runs — and no refit.

Writes only to the hub: ``article_plan/figures/final/F7a_estimator_design_years.*``. No field value in this file.

    make run E=article_fig_f7
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import yaml  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import freeze_v4_x106 as F  # noqa: E402  (build_v4: the frozen estimator, as X-107 runs it)
import p6_validation_x082 as X  # noqa: E402  (the QC'd laser, vertical mm)

from insar_wetlands import figstyle as fs  # noqa: E402
from insar_wetlands import v4  # noqa: E402
from insar_wetlands.paths import hub_root  # noqa: E402
from insar_wetlands.scoring import Protocol, guard_lockbox  # noqa: E402

HUB = hub_root()
FT = HUB / "08_deliverables" / "insar_full_test"
OUT = HUB / "08_deliverables" / "article_plan" / "figures" / "final"
LOYO = FT / "v4" / "spams_loyo_params.yaml"
PROTOCOL = FT / "protocol.yaml"
COL = {"laser": fs.INK, "plain": "#6e7781", "v4": fs.CLASS["mat"]}


def aligned(c: pd.DataFrame, truth: pd.Series, prot: Protocol) -> pd.DataFrame:
    """``c`` (time, seg, est) with ``truth`` at its epochs; each (segment, datum) run set to the laser at its first
    epoch with a laser reading."""
    c = c.copy()
    c["truth"] = truth.reindex(pd.DatetimeIndex(c.time).round("h")).to_numpy()
    c["datum"] = prot.datum(c.time)
    for _, g in c.groupby(["seg", "datum"]):
        first = g.dropna(subset=["truth"])
        c.loc[g.index, "est"] += (first.truth.iloc[0] - first.est.iloc[0]) if len(first) else np.nan
    return c


def series(track: str, frozen: dict, loyo: dict, truth: pd.Series, prot: Protocol) -> dict:
    """The plain chain, v4 and the band, from the laser's first reading on (earlier epochs are never scored)."""
    start = truth.first_valid_index()
    out = {"plain": aligned(F.build_v4(track, frozen, F.DESIGN, relevel=False), truth, prot),
           "v4": aligned(F.build_v4(track, frozen, F.DESIGN, relevel=True), truth, prot)}
    alts = []
    for year, p in loyo["fits"].items():
        alt = {**frozen, "spams": {**frozen["spams"], "x_p": p["x_p"], "x_e": p["x_e"], "tau_days": p["tau_days"]}}
        alts.append(aligned(F.build_v4(track, alt, F.DESIGN, relevel=True), truth, prot).set_index("time").est
                    .rename(str(year)))
    band = pd.concat(alts, axis=1)
    out["band"] = pd.DataFrame({"lo": band.min(axis=1), "hi": band.max(axis=1)})
    out["band"] = out["band"][out["band"].index >= start]
    for k in ("plain", "v4"):
        out[k] = out[k][out[k].time >= start]
    return out


def main() -> None:
    frozen = yaml.safe_load(F.PARAMS.read_text(encoding="utf-8"))
    loyo = yaml.safe_load(LOYO.read_text(encoding="utf-8"))
    prot = Protocol.load(PROTOCOL)
    truth = X.laser()[1]
    truth = truth[truth.index.year.isin(F.DESIGN)]
    guard_lockbox(truth.dropna().index, prot)
    with fs.style():
        fig, axes = plt.subplots(2, 1, figsize=fs.size("full", 130), sharex=True,
                                 gridspec_kw={"hspace": 0.16, "left": 0.08, "right": 0.98, "top": 0.86, "bottom": 0.07})
        record = {}
        fig.text(0.01, 0.875, "(a) design years (in-sample)", fontsize=8, fontweight="bold", ha="left", va="bottom")
        for ax, track in zip(axes, X.TRACKS):
            ax.set_facecolor("white")
            ax.grid(False)
            s = series(track, frozen, loyo, truth, prot)
            v = s["v4"].reset_index(drop=True)
            for a, b in v4.winter_bounds(v.time):
                ax.axvspan(v.time.iloc[a], v.time.iloc[b], color=fs.GRID, alpha=0.45, lw=0, zorder=0)
            for b0, b1 in prot.datum_breaks:
                ax.axvspan(b0, b1, facecolor="none", hatch="////", edgecolor=fs.MUTED, lw=0, zorder=0)
            ax.plot(truth.index, truth.to_numpy(), color=COL["laser"], lw=0.6, label="laser at P6 (quality-checked)")
            ax.fill_between(s["band"].index, s["band"].lo, s["band"].hi, color=COL["v4"], alpha=0.22, lw=0,
                            label="re-levelled chain, anchor fitted leaving one design year out (range of five)")
            for name, label in (("plain", "plain chain of consecutive pairs"), ("v4", "re-levelled short-pair chain (frozen)")):
                c = s[name]
                for _, g in c.groupby(["seg", "datum"]):
                    ax.plot(g.time, g.est, color=COL[name], lw=1.0 if name == "v4" else 0.9,
                            ls="-" if name == "v4" else (0, (3, 1.5)), marker="o", ms=1.6, label=label)
                    label = None
            # The article calls this whole figure "Fig. 7a" (design years) and keeps "7b" for the held-out years, so
            # the two tracks are named, not lettered.
            ax.set_ylabel("mat level (mm, vertical)")
            ax.text(0.995, 0.97, f"{track} track", transform=ax.transAxes, ha="right", va="top", fontsize=8,
                    fontweight="bold")
            record[track] = {"epochs_v4": int(v.est.notna().sum()), "winters": len(v4.winter_bounds(v.time))}
        handles, labels = axes[0].get_legend_handles_labels()
        handles += [plt.Rectangle((0, 0), 1, 1, color=fs.GRID, alpha=0.45),
                    plt.Rectangle((0, 0), 1, 1, facecolor="none", hatch="////", edgecolor=fs.MUTED, lw=0)]
        labels += ["winter re-levelled", "laser re-levelling: levels not compared across"]
        fig.legend(handles, labels, loc="upper left", bbox_to_anchor=(0.08, 0.995), ncol=2, fontsize=7)
        axes[1].set_xlim(truth.first_valid_index() - pd.Timedelta(days=10), pd.Timestamp("2026-01-05", tz="UTC"))
        axes[1].xaxis.set_major_locator(mdates.MonthLocator(bymonth=(1, 7)))
        axes[1].xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
        rec = fs.export(fig, OUT / "F7a_estimator_design_years", provenance={
            "inputs": ["hub 08_deliverables/insar_full_test/v4/v4_frozen_params.yaml (no refit)",
                       "hub 08_deliverables/insar_full_test/v4/spams_loyo_params.yaml (X-107 amendment_3)",
                       "hub 08_deliverables/insar_full_test/protocol.yaml",
                       "cube silver/site_hourly.csv.gz (laser), gold/plots_pair.csv.gz and silver/ifg_* (the chains)"],
            "estimator": "freeze_v4_x106.build_v4", "alignment": "per chain segment and laser datum, at the first "
            "epoch with a laser reading (scoring.level_views)", "design_years": list(F.DESIGN), "per_track": record,
            "brief": "08_deliverables/article_plan/PIPELINE_edits_leads_figures_2026-10-18.md (Stage 3, Fig. 7a)"})
        plt.close(fig)
    print("F7a →", OUT, rec["width_mm"], "×", rec["height_mm"], "mm", record)


if __name__ == "__main__":
    main()

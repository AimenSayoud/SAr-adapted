"""F7a — the estimator through the design years: the laser at P6, the plain chain of consecutive pairs and the
re-levelled short-pair chain (v4, frozen), per track, each winter shaded, with the anchor's spread as a band.

The band is X-107 amendment_3's error statement drawn: v4 rebuilt with each of the five SPAMS fits that leave one design
year out (``v4/spams_loyo_params.yaml``), everything else frozen; at each epoch the band spans the five rebuilt levels.
It is in-sample (the design years) and says how much the winter steps depend on the years the anchor is fitted on — not
how well v4 does on years it has not seen (Fig. 7b, after X-107). Each series is set to the laser at the first scored
epoch of each chain segment within each laser datum and, since 2026-10-31, within each hydrological year from
September: every winter's loss then starts from the laser at the season's end instead of accumulating over the years
(the display the figure review asked for; the scorer's winter steps are differences, unchanged by it). Nothing is
compared across a laser datum break (the 2022 re-levelling; the April 2024 step, X-147). The chains are built by
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
import article_fig_f7b as F7B  # noqa: E402  (panel b: nested leave-one-winter-out step errors)
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
HYEAR_START = 9                                    # September: re-anchor each series before each winter


def aligned(c: pd.DataFrame, truth: pd.Series, prot: Protocol) -> pd.DataFrame:
    """``c`` (time, seg, est) with ``truth`` at its epochs; each (segment, datum) run set to the laser at its first
    epoch with a laser reading."""
    c = c.copy()
    c["truth"] = truth.reindex(pd.DatetimeIndex(c.time).round("h")).to_numpy()
    c["datum"] = prot.datum(c.time)
    t = pd.DatetimeIndex(c.time)
    c["hyear"] = t.year - (t.month < HYEAR_START)                    # hydrological year from September
    for _, g in c.groupby(["seg", "datum", "hyear"]):
        first = g.dropna(subset=["truth"])
        c.loc[g.index, "est"] += (first.truth.iloc[0] - first.est.iloc[0]) if len(first) else np.nan
    return c


def _relevel_from_plain(plain_raw: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    R = F.spams_level(cfg)
    at = lambda t: float(R.reindex([pd.Timestamp(t).tz_convert(None).normalize()]).iloc[0])
    return v4.relevel(plain_raw.copy(), at, v4.winter_bounds(plain_raw.time))


def series(track: str, frozen: dict, loyo: dict, truth: pd.Series, prot: Protocol) -> dict:
    """The plain chain, v4 and the band, from the laser's first reading on (earlier epochs are never scored)."""
    start = truth.first_valid_index()
    plain_raw = F.build_v4(track, frozen, F.DESIGN, relevel=False)
    out = {"plain": aligned(plain_raw, truth, prot),
           "v4": aligned(_relevel_from_plain(plain_raw, frozen), truth, prot)}
    alts = []
    for year, p in loyo["fits"].items():
        alt = {**frozen, "spams": {**frozen["spams"], "x_p": p["x_p"], "x_e": p["x_e"], "tau_days": p["tau_days"]}}
        alts.append(aligned(_relevel_from_plain(plain_raw, alt), truth, prot).set_index("time").est
                    .rename(str(year)))
    band = pd.concat(alts, axis=1)
    out["band"] = pd.DataFrame({"lo": band.min(axis=1), "hi": band.max(axis=1)})
    out["band"] = out["band"][out["band"].index >= start]
    for k in ("plain", "v4"):
        out[k] = out[k][out[k].time >= start]
    return out


def draw_panel_a(axes, series_by_track: dict, truth: pd.Series, prot: Protocol) -> dict:
    record = {}
    for ax, track in zip(axes, X.TRACKS):
        ax.set_facecolor("white")
        ax.grid(False)
        s = series_by_track[track]
        v = s["v4"].reset_index(drop=True)
        for a, b in v4.winter_bounds(v.time):
            ax.axvspan(v.time.iloc[a], v.time.iloc[b], color=fs.GRID, alpha=0.45, lw=0, zorder=0)
        for b0, b1 in prot.datum_breaks:
            ax.axvspan(b0, b1, facecolor="none", hatch="////", edgecolor=fs.MUTED, lw=0, zorder=0)
        ax.plot(truth.index, truth.to_numpy(), color=COL["laser"], lw=0.6, label="laser at P6 (quality-checked)")
        ax.fill_between(s["band"].index, s["band"].lo, s["band"].hi, color=COL["v4"], alpha=0.22, lw=0,
                        label="re-levelled chain, anchor fitted leaving one design year out (range of five)")
        for name, label in (("plain", "plain chain of consecutive pairs"),
                            ("v4", "re-levelled short-pair chain (frozen)")):
            c = s[name]
            for _, g in c.groupby(["seg", "datum", "hyear"]):
                ax.plot(g.time, g.est, color=COL[name], lw=1.0 if name == "v4" else 0.9,
                        ls="-" if name == "v4" else (0, (3, 1.5)), marker="o", ms=1.6, label=label)
                label = None
        ax.set_ylabel("mat level (mm, vertical)")
        ax.text(0.995, 0.97, f"{track} track", transform=ax.transAxes, ha="right", va="top", fontsize=8,
                fontweight="bold")
        record[track] = {"epochs_v4": int(v.est.notna().sum()), "winters": len(v4.winter_bounds(v.time))}
    axes[1].set_xlim(truth.first_valid_index() - pd.Timedelta(days=10), pd.Timestamp("2025-12-28", tz="UTC"))
    axes[1].xaxis.set_major_locator(mdates.MonthLocator(bymonth=(1, 7)))
    axes[1].xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
    return record


def main() -> None:
    frozen = yaml.safe_load(F.PARAMS.read_text(encoding="utf-8"))
    loyo = yaml.safe_load(LOYO.read_text(encoding="utf-8"))
    prot = Protocol.load(PROTOCOL)
    truth = X.laser()[1]
    truth = truth[truth.index.year.isin(F.DESIGN)]
    guard_lockbox(truth.dropna().index, prot)
    series_by_track = {tr: series(tr, frozen, loyo, truth, prot) for tr in X.TRACKS}
    n_steps = pd.read_csv(F7B.SRC / "nested_steps.csv")
    verdict = pd.read_csv(F7B.SRC / "verdict.csv").iloc[0]
    with fs.style():
        fig, axes = plt.subplots(2, 1, figsize=fs.size("full", 130), sharex=True,
                                 gridspec_kw={"hspace": 0.16, "left": 0.08, "right": 0.97, "top": 0.86, "bottom": 0.07})
        fig.text(0.01, 0.875, "(a) design years (in-sample)", fontsize=8, fontweight="bold", ha="left", va="bottom")
        record = draw_panel_a(axes, series_by_track, truth, prot)
        handles, labels = axes[0].get_legend_handles_labels()
        handles += [plt.Rectangle((0, 0), 1, 1, color=fs.GRID, alpha=0.45),
                    plt.Rectangle((0, 0), 1, 1, facecolor="none", hatch="////", edgecolor=fs.MUTED, lw=0)]
        labels += ["winter re-levelled", "laser datum break: levels not compared across"]
        fig.legend(handles, labels, loc="upper left", bbox_to_anchor=(0.08, 0.995), ncol=2, fontsize=7)
        rec = fs.export(fig, OUT / "F7a_estimator_design_years", provenance={
            "inputs": ["hub 08_deliverables/insar_full_test/v4/v4_frozen_params.yaml (no refit)",
                       "hub 08_deliverables/insar_full_test/v4/spams_loyo_params.yaml (X-107 amendment_3)",
                       "hub 08_deliverables/insar_full_test/protocol.yaml",
                       "cube silver/site_hourly.csv.gz (laser), gold/plots_pair.csv.gz and silver/ifg_* (the chains)"],
            "estimator": "freeze_v4_x106.build_v4", "alignment": "per chain segment and laser datum, at the first "
            "epoch with a laser reading (scoring.level_views)", "design_years": list(F.DESIGN), "per_track": record,
            "brief": "08_deliverables/article_plan/PIPELINE_edits_leads_figures_2026-10-18.md (Stage 3, Fig. 7a)"})
        plt.close(fig)

        # Unified Fig. 7 combining (a) design years and (b) nested leave-one-winter-out step errors
        fig_u = plt.figure(figsize=fs.size("full", 195))
        gs = fig_u.add_gridspec(2, 1, height_ratios=[1.85, 1.0], hspace=0.34,
                                left=0.08, right=0.97, top=0.90, bottom=0.075)
        gs_a = gs[0].subgridspec(2, 1, hspace=0.16)
        ax_a0 = fig_u.add_subplot(gs_a[0])
        ax_a1 = fig_u.add_subplot(gs_a[1], sharex=ax_a0)
        plt.setp(ax_a0.get_xticklabels(), visible=False)
        rec_a = draw_panel_a([ax_a0, ax_a1], series_by_track, truth, prot)
        fig_u.text(0.01, 0.91, "(a) design years (in-sample)", fontsize=8, fontweight="bold", ha="left", va="bottom")
        ha, la = ax_a0.get_legend_handles_labels()
        ha += [plt.Rectangle((0, 0), 1, 1, color=fs.GRID, alpha=0.45),
               plt.Rectangle((0, 0), 1, 1, facecolor="none", hatch="////", edgecolor=fs.MUTED, lw=0)]
        la += ["winter re-levelled", "laser datum break: levels not compared across"]
        fig_u.legend(ha, la, loc="upper left", bbox_to_anchor=(0.08, 0.995), ncol=2, fontsize=7)

        gs_b = gs[1].subgridspec(1, 2, wspace=0.06)
        ax_b0 = fig_u.add_subplot(gs_b[0])
        ax_b1 = fig_u.add_subplot(gs_b[1], sharey=ax_b0)
        plt.setp(ax_b1.get_yticklabels(), visible=False)
        rec_b = F7B.draw_panel_b([ax_b0, ax_b1], n_steps)
        fig_u.text(0.01, 0.345, "(b) each laser winter left out in turn", fontsize=8, fontweight="bold",
                   ha="left", va="bottom")
        hb, lb = ax_b0.get_legend_handles_labels()
        fig_u.legend(hb, lb, loc="lower left", bbox_to_anchor=(0.27, 0.342), ncol=2, fontsize=7)
        rec_u = fs.export(fig_u, OUT / "F7_estimator", provenance={
            "inputs": ["hub 08_deliverables/insar_full_test/v4/v4_frozen_params.yaml (no refit)",
                       "hub 08_deliverables/insar_full_test/v4/spams_loyo_params.yaml (X-107 amendment_3)",
                       "hub 08_deliverables/insar_full_test/protocol.yaml",
                       "cube silver/site_hourly.csv.gz (laser), gold/plots_pair.csv.gz and silver/ifg_* (the chains)",
                       "hub 08_deliverables/insar_full_test/nested_winter_x135/nested_steps.csv (X-135)",
                       "hub 08_deliverables/insar_full_test/nested_winter_x135/verdict.csv"],
            "estimator": "freeze_v4_x106.build_v4", "design_years": list(F.DESIGN),
            "panel_a": rec_a, "panel_b": rec_b, "reading": str(verdict.reading),
            "choices": str(verdict.choices), "laser_rule": str(verdict.laser_rule)})
        plt.close(fig_u)
    print("F7a →", OUT, rec["width_mm"], "×", rec["height_mm"], "mm", record)
    print("F7_estimator →", OUT, rec_u["width_mm"], "×", rec_u["height_mm"], "mm")


if __name__ == "__main__":
    main()

"""F6, round 3 — why the routes that lost the mat's motion lost it: the share of pairs over which the laser at P6 moved
beyond a quarter wavelength along the line of sight, by pair length, with the pair lengths each length-dependent route
relied on drawn above it, coloured by its outcome.

Round 2 drew the systematic test as a verdict matrix — Table 1 as a picture. Table 1 keeps every route (written here
from the same ``article_plan/f6_routes.yaml``); the figure shows what a table cannot: that the failures follow the lag
at which the mat outruns a single pair (brief: hub ``article_plan/figures/briefs/F6_round3.md``).

Laser: the cube's QC'd level at each acquisition (``gold/plots_acq.csv.gz``, P6; snow-flagged acquisitions dropped),
converted to each track's line of sight by the cube's ``los_u`` at P6. Pairs of every length between acquisitions of
one track, never across the protocol's datum break nor across a filled stretch of the delivered record
(``gold/truth_laser_hourly.csv.gz``, ``laser_filled``). Design years only (``guard_lockbox``). Route spans: fixed where
the length defines the route, else read from the cube's pair inventory (``gold/plots_pair.csv.gz``) for the periods
the route ran on. Writes only to the hub: ``article_plan/figures/final/F6_pair_length.*``,
``article_plan/tables/f6_span_exceedance.csv`` and Table 1. No field value in this file.

    make run E=article_fig_f6
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import yaml  # noqa: E402

from insar_wetlands import figstyle as fs  # noqa: E402
from insar_wetlands.lag_exceedance import exceedance, pair_changes  # noqa: E402
from insar_wetlands.paths import cube_root, hub_root  # noqa: E402
from insar_wetlands.scoring import QUARTER_WAVE_MM, Protocol, guard_lockbox  # noqa: E402

HUB = hub_root()
PLAN = HUB / "08_deliverables" / "article_plan"
OUT = PLAN / "figures" / "final"
TRIES = HUB / "08_deliverables" / "independent_review" / "01_tries" / "tries.csv"
PROTOCOL = HUB / "08_deliverables" / "insar_full_test" / "protocol.yaml"
GOLD = cube_root() / "gold"
EDGES = [5.5, 6.5, 12.5, 18.5, 30.5, 48.5, 72.5, 120.5, 180.5, 300.5, 420.5, 780.0]   # days; 6 and 12 alone
TICKS = [6, 12, 24, 48, 96, 180, 365, 730]
# outcome → (colour, open bar, words), as in Table 1; "loses it" is grey, not the laser's vermilion, which (b) uses for
# the descending track
OUTCOME = {"holds": (fs.CLASS["mat"], False, "follows the motion"),
           "partly": (fs.CLASS["mat"], True, "partly: one track only"),
           "lost": ("#6e7781", False, "loses it")}
OUTCOME_WORDS = {"holds": "follows the motion", "helps": "helps (design years)", "partly": "partly: ascending only",
                 "mixed": "one track better, one worse: fails by its rule", "lost": "loses it",
                 "not_motion": "a signal, not the motion", "not_applicable": "fails its own conditions here",
                 "no_effect": "no gain, no change", "inconclusive": "inconclusive", "not_tried": "not tried here"}
PRODUCT_WORDS = {"HyP3": "HyP3", "SLC": "own SLC", "both": "HyP3 + own SLC", "sim": "simulation", "none": "—"}


def truthy(s: pd.Series) -> pd.Series:
    return s.astype(str).str.lower().isin(["true", "1", "1.0"])


def spec_rows() -> list[dict]:
    """The routes file's rows, each checked against the tries map."""
    spec = yaml.safe_load((PLAN / "f6_routes.yaml").read_text(encoding="utf-8"))
    known = set(pd.read_csv(TRIES)["id"])
    out = []
    for g in spec["groups"]:
        for r in g["rows"]:
            missing = set(r.get("tries", [])) - known
            if missing:
                raise SystemExit(f"F6: {r['route']} cites unknown tries {missing}")
            out.append({**r, "group": g["name"]})
    return out


def table1(rows: list[dict]) -> None:
    """Table 1: every route, its product, its outcome in words and the radar years of its tries."""
    tries = pd.read_csv(TRIES).set_index("id")
    pd.DataFrame([{"group": r["group"], "route": r["route"], "product": PRODUCT_WORDS[r["product"]],
                   "outcome": r.get("words") or OUTCOME_WORDS[r["outcome"]],
                   "years": "; ".join(sorted({str(tries.loc[i, "period"]) for i in r.get("tries", [])
                                              if i in tries.index})) if r["outcome"] != "not_tried" else "",
                   "tries": ", ".join(r.get("tries", []))}
                  for r in rows]).to_csv(PLAN / "tables" / "table1_systematic_test.csv", index=False)


def laser_exceedance() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per track: the laser's line-of-sight change for every pair of acquisitions, and its share beyond λ/4 by lag."""
    acq = pd.read_csv(GOLD / "plots_acq.csv.gz", low_memory=False)
    acq = acq[acq["plot"] == "P6"].copy()
    acq["t"] = pd.to_datetime(acq["acq_time_exact_utc"], utc=True)
    pp = pd.read_csv(GOLD / "plots_pair.csv.gz", usecols=["plot", "track", "grid", "los_u"])
    los_u = pp[(pp["plot"] == "P6") & (pp["grid"] == "40m")].groupby("track")["los_u"].median()
    hourly = pd.read_csv(GOLD / "truth_laser_hourly.csv.gz", usecols=["time_utc", "laser_filled"])
    filled = np.sort(pd.DatetimeIndex(pd.to_datetime(hourly.loc[truthy(hourly["laser_filled"]), "time_utc"],
                                                     utc=True)).asi8)
    breaks = Protocol.load(PROTOCOL).datum_breaks

    def barred(t1: pd.Timestamp, t2: pd.Timestamp) -> bool:
        if any(t1 < b1 and t2 > b0 for b0, b1 in breaks):
            return True
        i = np.searchsorted(filled, t1.value)
        return bool(i < len(filled) and filled[i] <= t2.value)

    level = acq["laser_surface_cm"].where(~truthy(acq["laser_snow"]))
    guard_lockbox(acq.loc[level.notna(), "t"], Protocol.load(PROTOCOL))
    pairs, shares = [], []
    for track, g in acq.groupby("track"):
        lv = (level[g.index] * 10 * los_u[track]).to_numpy()            # cm → mm, vertical → line of sight
        ch = pair_changes(g["t"], lv, EDGES[-1], barred).assign(track=track)
        pairs.append(ch)
        shares.append(exceedance(ch, QUARTER_WAVE_MM, EDGES).assign(track=track, los_u=los_u[track]))
    return pd.concat(pairs, ignore_index=True), pd.concat(shares, ignore_index=True)


def route_spans(rows: list[dict]) -> list[dict]:
    """[shortest, longest] pair length each length-dependent route relied on, in days."""
    inv = pd.read_csv(GOLD / "plots_pair.csv.gz", usecols=["period", "track", "pair", "revisit_days"])
    inv = inv.drop_duplicates(["period", "track", "pair"])
    out = []
    for r in rows:
        if "span_days" in r:
            lo, hi = r["span_days"]
        elif "span_network" in r:
            s = r["span_network"]
            d = inv[inv["period"].astype(str).isin([str(p) for p in s["periods"]])]["revisit_days"]
            counts = d.value_counts()
            keep = counts[counts >= s.get("min_pairs", 1)].index
            keep = keep[keep >= s.get("min_days", 0)]
            if keep.empty:
                raise SystemExit(f"F6: no pair lengths for {r['route']} in {s}")
            lo, hi = float(keep.min()), float(keep.max())
        else:
            continue
        out.append({"route": r["route"], "outcome": r["outcome"], "lo": float(lo), "hi": float(hi)})
    rank = {"holds": 0, "partly": 1, "lost": 2}
    out.sort(key=lambda s: (rank.get(s["outcome"], 3), s["lo"], s["hi"]))
    return out


def draw(spans: list[dict], shares: pd.DataFrame) -> plt.Figure:
    fig = plt.figure(figsize=fs.size("full", 150))
    gs = fig.add_gridspec(2, 1, height_ratios=[len(spans) * 0.42, 4.2], hspace=0.08,
                          left=0.33, right=0.97, top=0.93, bottom=0.09)
    ax_r, ax_e = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
    for ax in (ax_r, ax_e):
        ax.set_xscale("log")
        ax.set_xlim(5, 800)
        for d in (6, 12, 24, 365):
            ax.axvline(d, color=fs.GRID, lw=0.6, zorder=0)
    # (a) routes over the pair lengths they relied on
    for i, s in enumerate(spans):
        colour, hollow, _ = OUTCOME[s["outcome"]]
        if s["hi"] / s["lo"] < 1.05:
            ax_r.plot(s["lo"], i, "o", ms=6, mfc="white" if hollow else colour, mec=colour, mew=1.3)
        else:
            ax_r.plot([s["lo"], s["hi"]], [i, i], color=colour, lw=6, solid_capstyle="butt")
            if hollow:
                ax_r.plot([s["lo"], s["hi"]], [i, i], color="white", lw=3.4, solid_capstyle="butt")
    ax_r.set_yticks(range(len(spans)), [s["route"] for s in spans], fontsize=7)
    ax_r.set_ylim(len(spans) - 0.4, -0.6)
    ax_r.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    ax_r.tick_params(axis="y", length=0)
    ax_r.spines[["left", "bottom"]].set_visible(False)
    handles = [plt.Line2D([], [], color=c, lw=6 if not h else 0, marker=None if not h else "o", mfc="white",
                          mec=c, mew=1.3, ms=6, label=w) for c, h, w in OUTCOME.values()]
    ax_r.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.0, 1.02), ncol=3, fontsize=7,
                handlelength=1.8, columnspacing=1.6)
    fs.panel_label(ax_r, "a", x=-0.47, y=1.0)
    for d, lab in ((6, "6 d"), (12, "12 d"), (24, "24 d"), (365, "1 yr")):
        ax_r.text(d, -0.75, lab, ha="center", va="bottom", fontsize=7, color=fs.MUTED)
    # (b) share of pairs over which the laser moved beyond λ/4
    for track, g in shares.groupby("track"):
        c = fs.TRACK[track]
        ls = "-" if track == "ascending" else (0, (3.5, 1.8))
        ax_e.errorbar(g["lag_median_days"], g["share"],
                      yerr=[g["share"] - g["share_lo95"], g["share_hi95"] - g["share"]],
                      color=c, ls=ls, marker="o" if track == "ascending" else "^", ms=4, lw=1.0, elinewidth=0.8,
                      capsize=0, label=f"{track} (line of sight)")
    ax_e.set_ylim(-0.02, 1.02)
    ax_e.yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0))
    ax_e.set_xticks(TICKS, [str(t) for t in TICKS])
    ax_e.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax_e.set_xlabel("pair length (days)")
    ax_e.set_ylabel("pairs with the laser beyond λ/4\n(laser at P6, 95 % interval)")
    ax_e.legend(loc="lower right")
    fs.panel_label(ax_e, "b", x=-0.47, y=1.0)
    return fig


def main() -> None:
    rows = spec_rows()
    pairs, shares = laser_exceedance()
    spans = route_spans(rows)
    shares.to_csv(PLAN / "tables" / "f6_span_exceedance.csv", index=False)
    with fs.style():
        fig = draw(spans, shares)
        rec = fs.export(fig, OUT / "F6_pair_length", provenance={
            "inputs": ["cube gold/plots_acq.csv.gz (P6 laser at each acquisition, snow-flagged dropped)",
                       "cube gold/plots_pair.csv.gz (los_u at P6; pair inventory for the route spans)",
                       "cube gold/truth_laser_hourly.csv.gz (filled stretches)",
                       "hub 08_deliverables/insar_full_test/protocol.yaml (datum breaks, lockbox)",
                       "hub 08_deliverables/article_plan/f6_routes.yaml (routes, outcomes, spans)"],
            "quarter_wave_mm": QUARTER_WAVE_MM, "lag_edges_days": EDGES,
            "pairs_by_track": pairs.groupby("track").size().to_dict(),
            "route_spans": spans, "interval": "Wilson 95 %, pairs treated as independent (overlapping pairs: narrower "
                                               "than the record's real uncertainty)",
            "brief": "08_deliverables/article_plan/figures/briefs/F6_round3.md"})
        plt.close(fig)
    table1(rows)
    print("F6 →", OUT, rec["width_mm"], "×", rec["height_mm"], "mm;", len(pairs), "laser pairs")


if __name__ == "__main__":
    main()

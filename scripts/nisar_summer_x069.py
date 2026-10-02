"""X-069: NISAR L-band in summer 2026 — does L-band keep the floating mat coherent where C-band loses it?

X-068 compared NISAR's beta interferograms (Oct 2025 – Jan 2026) with Sentinel-1: late autumn and winter, when C-band does
not lose the mat, so L-band's expected advantage could not show. NISAR's calibrated ("provisional") GUNW products were
released on 2026-07-20 and cover Rzecin from June 2026 on four tracks (057 A and 158 A near 03:30 UTC, 008 D and 080 D
near 18:00 UTC). This reads every one of them over the site (HTTP range requests on the AOI window only, as X-068) and
reports coherence by zone, the mat − grassland coherence gap, by track, pair length and surface state.

Sentinel-1 2026 is not processed yet (decision D-024), so the C-band reference is **summer in other years**: the
2025 crops (S1A + S1C, 6- and 12-day) and 2020–2021 (S1A + S1B), May–September, same zones, same wet/frozen rule.
Labelled as such; replaced by 2026 when D-024's crops exist.

Surface state at each NISAR time: Open-Meteo hourly (2026) with the station calibration of fusion v2 (X-066): wet =
RH ≥ 90 % in the reanalysis (≙ 95 % at the station) or rain in the previous 3 h; frozen = air − its mean bias ≤ 0 °C.
No field data are read (no laser or water table exists for 2026); outputs go to the hub with the other NISAR results.

    INSAR_DRIVE_ROOT=../local/drive_pristine PYTHONPATH=src python scripts/nisar_summer_x069.py [--reuse]
"""
from __future__ import annotations

import argparse
import json
import sys
import urllib.parse
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fusion_v2_flags_2025 as FL  # noqa: E402
import nisar_lband_x068 as N  # noqa: E402

from insar_wetlands import field  # noqa: E402
from insar_wetlands.bootstrap import start  # noqa: E402
from insar_wetlands.regrid import stack_to_template  # noqa: E402
from insar_wetlands.stack import list_pairs, load_layer  # noqa: E402

HUB = Path(__file__).resolve().parents[3]
OUT = HUB / "08_deliverables" / "nisar_x069"
COLLECTION = "C2854335566-ASF"            # NISAR_L2_GUNW_PROVISIONAL_V1
PERIOD = ("2026-04-01T00:00:00Z", "2026-10-15T00:00:00Z")
SUMMER = range(5, 10)                     # May–September
S1_HOUR = {"ascending": "16:36", "descending": "05:09"}


def granules(bbox) -> list[dict]:
    out, page = [], 1
    while True:
        q = urllib.parse.urlencode({"collection_concept_id": COLLECTION, "bounding_box": ",".join(map(str, bbox)),
                                    "temporal": ",".join(PERIOD), "page_size": 200, "page_num": page})
        with urllib.request.urlopen(f"{N.CMR}?{q}", timeout=120) as r:
            feed = json.load(r)["feed"]["entry"]
        for e in feed:
            h5 = [lk["href"] for lk in e.get("links", []) if lk.get("href", "").endswith(".h5")]
            if h5:
                out.append({"title": e["title"], "url": h5[0]})
        if len(feed) < 200:
            break
        page += 1
    return sorted({g["title"]: g for g in out}.values(), key=lambda g: g["title"])


def track_of(title: str) -> tuple[str, str]:
    """NISAR_L2_PR_GUNW_<cycle>_<track>_<A|D>_… → ('057', 'A')."""
    p = title.split("_")
    return p[5], p[6]


def flags(om: pd.DataFrame, meta: dict, times) -> pd.DataFrame:
    cal = om.copy()
    cal["Air_2m"] = cal["Air_2m"] - meta["air_bias_c"]
    return field.surface_wetness_at(cal, times, rh_wet=meta["rh_wet_threshold_open_meteo"])


def open_meteo(ctx, start_date: str, end_date: str) -> pd.DataFrame:
    lat, lon = ctx.cfg["site"]["centroid"][1], ctx.cfg["site"]["centroid"][0]
    cache = OUT / f"open_meteo_{start_date}_{end_date}.csv"
    if cache.exists():
        d = pd.read_csv(cache, parse_dates=["time"]).set_index("time")
        d.index = pd.to_datetime(d.index, utc=True)
        return d
    d = FL.fetch(lat, lon, start_date, end_date)
    d.to_csv(cache)
    return d


def sentinel1_summer(ctx, om_all: pd.DataFrame, meta: dict, om26: pd.DataFrame | None = None) -> pd.DataFrame:
    """C-band reference: every May–September pair ≤ 24 d of the 2026 (same weeks as NISAR; 20 m products put on the 40 m
    template, D-024), 2025 and 2020–2021 crops, zone median coherence."""
    zones = {z: ctx.zones[z].values.astype(bool) for z in "ABCD"}
    st = field.load_wtd_hourly()
    rows = []
    for track in S1_HOUR:
        for label, root in (("2026 (S1A+S1C+S1D, same weeks)", HUB / f"05_code/local/s1_2026/hyp3_cropped_{track}"),
                            ("2025 (S1A+S1C)", HUB / f"05_code/local/s1_2025/hyp3_cropped_{track}"),
                            ("2020–2021 (S1A+S1B)", HUB / f"05_code/local/s1_2020_2021/hyp3_cropped_{track}")):
            if not root.exists():
                continue
            pairs = [p for p in list_pairs(root) if pd.Timestamp(p[:8]).month in SUMMER]
            pairs = [p for p in pairs if (pd.Timestamp(p[9:]) - pd.Timestamp(p[:8])).days <= 24]
            if not pairs:
                continue
            co = stack_to_template(load_layer(root, "corr", pairs), ctx.template.x.values, ctx.template.y.values).values
            t1 = [pd.Timestamp(f"{p[:8]} {S1_HOUR[track]}", tz="UTC") for p in pairs]
            t2 = [pd.Timestamp(f"{p[9:]} {S1_HOUR[track]}", tz="UTC") for p in pairs]
            if label.startswith("2026"):
                if om26 is None:
                    continue
                f1, f2 = flags(om26, meta, t1), flags(om26, meta, t2)
            elif label.startswith("2025"):
                f1, f2 = flags(om_all, meta, t1), flags(om_all, meta, t2)
            else:
                f1, f2 = field.surface_wetness_at(st, t1), field.surface_wetness_at(st, t2)
            for k, p in enumerate(pairs):
                rows.append({"sensor": "Sentinel-1 C-band", "period": label, "track": track, "pair": p, "t1": t1[k], "t2": t2[k],
                             "dt": (t2[k] - t1[k]).days, "wet": bool(f1.wet[k] or f2.wet[k]), "frozen": bool(f1.frozen[k] or f2.frozen[k]),
                             **{f"coh_{z}": float(np.nanmedian(co[k][zones[z]])) for z in "ABC"}})
    return pd.DataFrame(rows)


def summarise(d: pd.DataFrame, keys: list[str]) -> pd.DataFrame:
    """Medians over pairs, frozen pairs out, and frames that do not cover the site (no mat or grassland value) out."""
    g = d[~d.frozen & d.coh_A.notna() & d.coh_C.notna()].groupby(keys)
    s = g.agg(pairs=("coh_A", "size"), mat=("coh_A", "median"), grassland=("coh_C", "median"), lake=("coh_B", "median")).reset_index()
    s["mat_minus_grassland"] = g.apply(lambda x: float(np.median(x.coh_A - x.coh_C)), include_groups=False).values
    return s


def fig(d: pd.DataFrame, c: pd.DataFrame):
    d = d[d.coh_A.notna() & d.coh_C.notna()]                 # frames that do not cover the site out
    fig, ax = plt.subplots(1, 2, figsize=(15, 5.2))
    for (trk, ad), g in d.groupby(["track_no", "dir"]):
        g = g.sort_values("t1")
        for dt, gg in g.groupby("dt"):
            ax[0].plot(gg.t1 + (gg.t2 - gg.t1) / 2, gg.coh_A, "o-" if dt == 12 else "s:", lw=1, label=f"NISAR {trk} {ad}, {dt} d — mat")
    ax[0].set_ylim(0, 1); ax[0].set_ylabel("median coherence, mat (zone A)"); ax[0].legend(fontsize=7)
    ax[0].set_title("NISAR L-band, summer 2026, every pair", fontsize=9)
    labels, data = [], []
    for (sensor, period, dt), g in pd.concat([d.assign(period="2026", sensor="NISAR L-band"), c]).query("~frozen").groupby(["sensor", "period", "dt"]):
        if len(g) >= 3:
            short = period.split(" ")[0]
            labels.append(f"{'L' if sensor.startswith('NISAR') else 'C'} {short}\n{dt} d, n {len(g)}"); data.append(g.coh_A - g.coh_C)
    ax[1].boxplot(data, tick_labels=labels, showfliers=False)
    ax[1].axhline(0, color="0.6", lw=0.8)
    ax[1].set_ylabel("coherence, mat − grassland (per pair)"); ax[1].tick_params(axis="x", labelsize=7, rotation=0)
    ax[1].text(0.01, 0.02, "L = NISAR L-band 2026 · C = Sentinel-1 C-band, summers 2020–2021 and 2025", transform=ax[1].transAxes, fontsize=7)
    ax[1].set_title("Does the mat lose coherence relative to the grassland? (summer, not frozen)", fontsize=9)
    fig.tight_layout(); fig.savefig(OUT / "fig_nisar_summer.png", dpi=130); plt.close(fig)


def md(df: pd.DataFrame) -> str:
    lines = ["| " + " | ".join(map(str, df.columns)) + " |", "|" + "---|" * len(df.columns)]
    for r in df.itertuples(index=False):
        lines.append("| " + " | ".join(f"{v:.2f}" if isinstance(v, float) else str(v) for v in r) + " |")
    return "\n".join(lines)


def findings(sN: pd.DataFrame, sC: pd.DataFrame) -> list[str]:
    """Generated reading of the two summary tables (no typed number)."""
    a = sN[sN.track_no == "all"].set_index("dt")
    out = []
    same = sC[sC.period.str.startswith("2026")]
    if 12 in a.index and len(same):
        s12, s6 = same[same.dt == 12], same[same.dt == 6]
        out.append(f"- **Same weeks (Sentinel-1 2026, D-024):** mat coherence C-band 12-day {s12.mat.min():.2f}–{s12.mat.max():.2f}, "
                   f"6-day {s6.mat.min():.2f}–{s6.mat.max():.2f}, against L-band 12-day {a.loc[12, 'mat']:.2f}; mat minus grassland "
                   f"C 12-day {s12.mat_minus_grassland.min():+.2f} to {s12.mat_minus_grassland.max():+.2f}, L 12-day {a.loc[12, 'mat_minus_grassland']:+.2f}.")
    if 12 in a.index:
        c12, c6 = sC[(sC.dt == 12) & ~sC.period.str.startswith("2026")], sC[(sC.dt == 6) & ~sC.period.str.startswith("2026")]
        out.append(f"- **Summer, mat:** NISAR L-band 12-day coherence {a.loc[12, 'mat']:.2f} ({int(a.loc[12, 'pairs'])} pairs) against "
                   f"Sentinel-1 C-band 12-day {c12.mat.min():.2f}–{c12.mat.max():.2f} and 6-day {c6.mat.min():.2f}–{c6.mat.max():.2f} "
                   "(other summers, both tracks).")
        out.append(f"- **The H2 gap:** mat minus grassland {a.loc[12, 'mat_minus_grassland']:+.2f} for L-band 12-day, against "
                   f"{c12.mat_minus_grassland.min():+.2f} to {c12.mat_minus_grassland.max():+.2f} for C-band 12-day — at L-band the mat is not "
                   "the less coherent surface in summer.")
    long = sN[(sN.track_no == "all") & (sN.dt >= 24)]
    if len(long):
        out.append(f"- **Longer L-band pairs** ({', '.join(f'{int(r.dt)} d' for r in long.itertuples())}): mat {long.mat.min():.2f}–{long.mat.max():.2f} "
                   f"({int(long.pairs.sum())} pairs) — near the floor where everything has decorrelated; the advantage is a 12-day one.")
    fl_c = HUB / "08_deliverables" / "lake_x064" / "masks_summary.csv"
    if fl_c.exists():
        m = pd.read_csv(fl_c)
        fC = m[(m["mask"] == "mat (zone A)") & m.coh_floor_gt96d_median.notna()].coh_floor_gt96d_median
        if len(fC) and len(long):
            out.append(f"- **Above each sensor's floor:** C-band's floor on the mat (pairs > 96 d, X-064) is {fC.min():.2f}–{fC.max():.2f}; "
                       f"L-band's is not measured yet (its longest pairs here, {long.mat.min():.2f}–{long.mat.max():.2f}, bound it from above). "
                       "Compared above the floors, the L-band 12-day margin is the larger.")
    return out


def readme(d, c, sN, sC, winter):
    L = ["# X-069 — NISAR L-band in summer 2026 against Sentinel-1 C-band (exploratory)", "",
         "Generated by `05_code/SAr-adapted/scripts/nisar_summer_x069.py`; every number from its tables. NISAR provisional (calibrated) "
         "GUNW, every pair over Rzecin from June 2026, four tracks; read by HTTP range requests (AOI window only). The C-band reference "
         "is **summer of other years** (2025, 2020–2021) until Sentinel-1 2026 is processed (D-024).", "",
         f"NISAR pairs read: {len(d)} ({int(d.coh_A.notna().sum())} cover the site) on tracks {', '.join(sorted({f'{a} {b}' for a, b in zip(d.track_no, d.dir)}))}; "
         f"Sentinel-1 summer pairs ≤ 24 d: {len(c)}.", "",
         "## Findings", "", *findings(sN, sC), "",
         "## NISAR L-band, summer 2026 (not frozen), by track and pair length", "", md(sN), "",
         "## Sentinel-1 C-band, summer (May–September) of other years (not frozen), by period, track and pair length", "", md(sC), "",
         "## NISAR winter (X-068, beta) for contrast", "", md(winter) if len(winter) else "(X-068 table not found)", "",
         "`mat_minus_grassland`: median over pairs of the mat's coherence minus the grassland's in the same pair — the H2 gap. "
         "L-band would help where the mat keeps coherence (high `mat`) and the gap closes.", "",
         "## How to read this", "",
         "- Pair lengths differ: NISAR repeats every 12 days (its shortest pairs); Sentinel-1 has 6-day pairs in 2025 and 2020–2021. "
         "Compare L 12-day with C 12-day first, then ask whether L 12-day beats C 6-day.",
         "- Times differ: NISAR ascending passes near 03:30 UTC (night, often dew), descending near 18:00 UTC; Sentinel-1 at 16:36 and "
         "05:09 UTC. Wet pairs are kept here and flagged in `nisar_pairs.csv`; the summary excludes frozen pairs only.",
         "- Grids differ: NISAR coherence on its 80 m unwrapped grid (20 m wrapped grid in the table too); Sentinel-1 on 40 m. "
         "Coherence estimators differ; compare gaps (mat − grassland) more than absolute levels.",
         "- Different years for the C-band reference; D-024 (Sentinel-1 2026) makes it the same weeks.", "",
         "## Files", "", "| File | Content |", "|---|---|",
         "| `nisar_pairs.csv` | every NISAR pair: track, dates, coherence by zone (80 m, 20 m), surface state |",
         "| `sentinel1_summer_reference.csv` | every Sentinel-1 summer pair ≤ 24 d of 2025 and 2020–2021: coherence by zone, state |",
         "| `summary_nisar.csv`, `summary_sentinel1.csv`, `summary.json` | the tables above |",
         "| `fig_nisar_summer.png` | NISAR mat coherence per pair; mat − grassland gap by sensor, period and pair length |"]
    (OUT / "README.md").write_text("\n".join(L) + "\n")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--reuse", action="store_true", help="reuse nisar_pairs.csv instead of reading NISAR again")
    a = ap.parse_args(argv)
    OUT.mkdir(parents=True, exist_ok=True)
    ctx = start("nisar_summer_x069", mount=False, git=False)
    tpl = ctx.template
    zones = {z: ctx.zones[z].values.astype(bool) for z in "ABCD"}
    px = pd.read_csv(HUB / "08_deliverables" / "field_first" / "plot_pixels.csv").set_index("plot")
    p6_xy = (float(px.loc["P6", "E"]), float(px.loc["P6", "N"]))
    meta = json.loads((HUB / "08_deliverables" / "fusion_v2" / "flags_2025" / "open_meteo_vs_station.json").read_text())
    cache = OUT / "nisar_pairs.csv"
    if a.reuse and cache.exists():
        d = pd.read_csv(cache, parse_dates=["t1", "t2"])
    else:
        bbox = tpl.rio.transform_bounds("EPSG:4326")
        gs = granules(bbox)
        print(f"{len(gs)} provisional GUNW over the site")
        rows = []
        done = set(pd.read_csv(cache).granule) if cache.exists() else set()
        old = pd.read_csv(cache, parse_dates=["t1", "t2"]) if cache.exists() else pd.DataFrame()
        for i, g in enumerate(gs, 1):
            if g["title"] in done:
                continue
            try:
                r = N.read_pair(g, tpl, zones, p6_xy)
            except Exception as e:  # noqa: BLE001 — a granule that cannot be read is reported and skipped
                print(f"  ! {g['title']}: {type(e).__name__}: {e}", flush=True)
                continue
            r["track_no"], r["dir"] = track_of(g["title"])
            rows.append(r)
            print(f"  {i}/{len(gs)} {r['track_no']}{r['dir']} {r['t1']:%Y-%m-%d}→{r['t2']:%Y-%m-%d} mat {r['u80_coherenceMagnitude_A']:.2f} "
                  f"grass {r['u80_coherenceMagnitude_C']:.2f} ({r['http_requests']} requests)", flush=True)
            pd.concat([old, pd.DataFrame(rows)], ignore_index=True).to_csv(cache, index=False)
        d = pd.read_csv(cache, parse_dates=["t1", "t2"])
    d["t1"], d["t2"] = pd.to_datetime(d.t1, utc=True), pd.to_datetime(d.t2, utc=True)
    d["dt"] = ((d.t2 - d.t1).dt.total_seconds() / 86400).round().astype(int)
    d["track_no"] = d.track_no.astype(str).str.zfill(3)
    om26 = open_meteo(ctx, "2026-04-01", (pd.Timestamp.utcnow() - pd.Timedelta(days=6)).strftime("%Y-%m-%d"))
    f1, f2 = flags(om26, meta, d.t1), flags(om26, meta, d.t2)
    d["wet"], d["frozen"] = (f1.wet.values | f2.wet.values), (f1.frozen.values | f2.frozen.values)
    d = d.rename(columns={"u80_coherenceMagnitude_A": "coh_A", "u80_coherenceMagnitude_B": "coh_B", "u80_coherenceMagnitude_C": "coh_C"})
    d.to_csv(cache.with_name("nisar_pairs_flagged.csv"), index=False)
    om_all = pd.read_csv(HUB / "06_data" / "local_small" / "open_meteo" / "open_meteo_rzecin_hourly_2020_2025.csv", parse_dates=["time"]).set_index("time")
    om_all.index = pd.to_datetime(om_all.index, utc=True)
    c = sentinel1_summer(ctx, om_all, meta, om26)
    c.to_csv(OUT / "sentinel1_summer_reference.csv", index=False)
    sN = summarise(d.assign(sensor="NISAR L-band"), ["track_no", "dir", "dt"])
    sN_all = summarise(d.assign(sensor="NISAR L-band"), ["dt"]).assign(track_no="all", dir="")
    sN = pd.concat([sN, sN_all[sN.columns]], ignore_index=True)
    sC = summarise(c, ["period", "track", "dt"])
    sN.to_csv(OUT / "summary_nisar.csv", index=False)
    sC.to_csv(OUT / "summary_sentinel1.csv", index=False)
    wint = HUB / "08_deliverables" / "nisar_x068" / "nisar_pairs.csv"
    winter = pd.DataFrame()
    if wint.exists():
        w = pd.read_csv(wint)
        winter = pd.DataFrame([{"pairs": len(w), "mat": float(w.u80_coherenceMagnitude_A.median()),
                                "grassland": float(w.u80_coherenceMagnitude_C.median()), "lake": float(w.u80_coherenceMagnitude_B.median()),
                                "mat_minus_grassland": float((w.u80_coherenceMagnitude_A - w.u80_coherenceMagnitude_C).median())}])
    (OUT / "summary.json").write_text(json.dumps({"nisar_pairs": len(d), "sentinel1_pairs": len(c),
                                                  "nisar": sN.to_dict("records"), "sentinel1": sC.to_dict("records"),
                                                  "nisar_winter_x068": winter.to_dict("records")}, indent=1, default=str))
    fig(d, c)
    readme(d, c, sN, sC, winter)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()

"""X-063: backscatter by zone on both tracks, 2020–2024 — phaseDter's backscatter, run on this Mac.

phaseDter (Colab notebook) reads the ascending 2022–2024 RTC stack only. With the 2020–2024 stacks of
both tracks (``scripts/rtc_download.py``, C-048) the same quantities are computed here, on the same
zones (the project's A–D, defined once on the radar grid — the descending crops share it):

1. **Reproduction** — the new ascending stack equals ``rtc_dualpol_stack.nc`` on its 85 dates, and the
   2022–2024 zone statistics equal the recorded phaseDter run (σ0 VV mean and temporal SD, D_A).
2. **Zone series** — per date and zone: median σ0 VV, VH, VH−VV and dual-pol RVI, both tracks.
3. **Summary** per track × period (2020–2021, 2022–2024, 2020–2024) × zone: phaseDter's statistics
   plus VH−VV and RVI.
4. **Sentinel-1A vs 1B** — 2020–2021 mixes two satellites (6-day interleave, so the seasons are the
   same): their zone medians are compared, because a calibration offset would pass for a change
   between the periods.
5. **Clean lake** — phaseDter's negative control (persistent open water: WorldCover water or S2 NDWI
   > 0.2 most of the time, inside the outline): its backscatter per track × period, and the coherence
   of every interferogram 2020–2024 on both tracks for the mat, the clean lake, the grassland, outside.
6. **Dew and rain** — anomalies (each date minus its calendar-month median, frozen dates out) on wet
   vs dry overpasses (X-050 flags: RH ≥ 95 % or rain in the previous 3 h), per zone and for the mat
   against the grassland (A − C). Exploratory.

No inversion, no null realisations: medians of a 203 × 129 × 138 stack. Outputs to the private hub
(the wet flags come from the field station).

    INSAR_DRIVE_ROOT=../local/drive_pristine PYTHONPATH=src python scripts/backscatter_x063.py
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import xarray as xr  # noqa: E402
from scipy.stats import mannwhitneyu  # noqa: E402

from insar_wetlands import field  # noqa: E402
from insar_wetlands.bootstrap import start  # noqa: E402
from insar_wetlands.paths import make_paths  # noqa: E402
from insar_wetlands.stack import list_pairs  # noqa: E402
from insar_wetlands.stratify import (  # noqa: E402
    amplitude_dispersion_from_rtc,
    backscatter_by_zone,
    clean_lake_mask,
    coherence_by_zone_stream,
    dual_pol_rvi,
    load_worldcover,
    zone_backscatter_series,
)

TRACKS = {"ascending": "rtc_dualpol_2020_2024.nc", "descending": "rtc_dualpol_2020_2024_descending.nc"}
PERIODS = {"2020–2021": (2020, 2021), "2022–2024": (2022, 2024), "2020–2024": (2020, 2024)}
VARS = {"vv_db": "σ0 VV (dB)", "vh_db": "σ0 VH (dB)", "ratio_vh_vv_db": "VH − VV (dB)", "rvi": "RVI"}
PS_DA = 0.25                     # phaseDter's indicative PS-candidate threshold on D_A


def md(df: pd.DataFrame) -> str:
    """A markdown table without extra dependencies."""
    fmt = lambda v: f"{v:.3g}" if isinstance(v, float) else str(v)  # noqa: E731
    head = "| " + " | ".join(map(str, df.columns)) + " |"
    return "\n".join([head, "|" + "---|" * len(df.columns)]
                     + ["| " + " | ".join(fmt(v) for v in row) + " |" for row in df.itertuples(index=False)])


def years(ds: xr.Dataset, y0: int, y1: int) -> xr.Dataset:
    t = pd.DatetimeIndex(ds.time.values)
    return ds.isel(time=np.flatnonzero((t.year >= y0) & (t.year <= y1)))


def satellite(ids) -> list[str]:
    """S1A / S1B from the RTC item id (a mosaic's ids share the satellite)."""
    return [str(i)[:3] for i in ids]


def summary(ds: xr.Dataset, zones: dict) -> pd.DataFrame:
    """phaseDter's per-zone statistics on one stack, plus VH−VV and RVI."""
    bz = backscatter_by_zone(ds, zones)
    da = amplitude_dispersion_from_rtc(ds)
    rvi = dual_pol_rvi(ds, reduce="median")
    rows = []
    for z in "ABCD":
        m = zones[z].values
        d = da.values[m & np.isfinite(da.values)]
        r = rvi.values[m & np.isfinite(rvi.values)]
        rows.append({"zone": z, "n_px": int(m.sum()),
                     "da_median": float(np.median(d)), f"da_frac_below_{PS_DA}": float((d < PS_DA).mean()),
                     "rvi_median": float(np.median(r))})
    return bz.merge(pd.DataFrame(rows), on="zone")


def month_anomaly(s: pd.DataFrame, cols) -> pd.DataFrame:
    """Each value minus the median of its track × zone × calendar month (frozen dates excluded)."""
    s = s.copy()
    for c in cols:
        s[f"{c}_anom"] = np.nan
    clim = s[~s.frozen].groupby(["track", "zone", s.date.dt.month])[cols].transform("median")
    s.loc[~s.frozen, [f"{c}_anom" for c in cols]] = (s.loc[~s.frozen, cols] - clim).to_numpy()
    return s


def wet_dry(series: pd.DataFrame) -> pd.DataFrame:
    """Wet minus dry median anomaly, per track × zone (and A − C), Mann–Whitney p (two-sided)."""
    cols = ["vv_db", "vh_db", "rvi"]
    s = month_anomaly(series, cols)
    s = s[~s.frozen]
    rows = []
    for track in s.track.unique():
        st = s[s.track == track]
        groups = {z: st[st.zone == z] for z in "ABCD"}
        j = st[st.zone == "A"].set_index("date").join(st[st.zone == "C"].set_index("date"), rsuffix="_C", how="inner")
        groups["A − C"] = pd.DataFrame({"wet": j.wet, **{f"{c}_anom": j[f"{c}_anom"] - j[f"{c}_anom_C"] for c in cols}})
        for label, d in groups.items():
            for c in cols:
                wet_v = d.loc[d.wet.astype(bool), f"{c}_anom"].dropna()
                dry_v = d.loc[~d.wet.astype(bool), f"{c}_anom"].dropna()
                if len(wet_v) < 5 or len(dry_v) < 5:
                    continue
                rows.append({"track": track, "zone": label, "variable": c, "n_wet": len(wet_v), "n_dry": len(dry_v),
                             "wet_minus_dry": float(wet_v.median() - dry_v.median()),
                             "p": float(mannwhitneyu(wet_v, dry_v).pvalue)})
    return pd.DataFrame(rows)


def satellite_offset(series: pd.DataFrame) -> pd.DataFrame:
    """2020–2021: median zone value on S1A dates minus S1B dates, Mann–Whitney p."""
    s = series[(series.date.dt.year <= 2021) & ~series.frozen]
    rows = []
    for (track, z), g in s.groupby(["track", "zone"]):
        for c in ("vv_db", "vh_db"):
            a, b = g.loc[g.satellite == "S1A", c].dropna(), g.loc[g.satellite == "S1B", c].dropna()
            if len(a) >= 5 and len(b) >= 5:
                rows.append({"track": track, "zone": z, "variable": c, "n_S1A": len(a), "n_S1B": len(b),
                             "S1A_minus_S1B_db": float(a.median() - b.median()),
                             "p": float(mannwhitneyu(a, b).pvalue)})
    return pd.DataFrame(rows)


def lake_coherence(ctx, extra: Path, zones: dict, lake: xr.DataArray) -> pd.DataFrame:
    """Per pair, both tracks, 2020–2024: mean coherence of the mat, the clean lake (as B), the
    grassland and outside — phaseDter's lake check, streamed one interferogram at a time."""
    zl = {"A": zones["A"], "B": lake, "C": zones["C"], "D": zones["D"]}
    out = []
    for track in TRACKS:
        for root in (make_paths(cfg=ctx.cfg, track=track).cropped, extra / f"hyp3_cropped_{track}"):
            if root.is_dir():
                df, _ = coherence_by_zone_stream(root, list_pairs(root), zl, ctx.template)
                out.append(df.assign(track=track))
    df = pd.concat(out, ignore_index=True).drop_duplicates(["track", "pair", "zone"])
    df["year"] = pd.to_datetime(df.pair.str[:8]).dt.year
    return df


def lake_coherence_table(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (track, period), (y0, y1) in [((t, p), yy) for t in TRACKS for p, yy in PERIODS.items()]:
        d = df[(df.track == track) & df.year.between(y0, y1)]
        for subset, dd in (("pairs ≤ 24 d", d[d.dt_days <= 24]), ("all pairs", d)):
            med = dd.groupby("zone").mean_coh.median()
            rows.append({"track": track, "period": period, "pairs": subset, "n_pairs": int(dd.pair.nunique()),
                         **{f"coh_{'clean_lake' if z == 'B' else z}": float(med.get(z, np.nan)) for z in "ABCD"}})
    return pd.DataFrame(rows)


def figures(out: Path, series: pd.DataFrame, stacks: dict, zones: dict) -> None:
    colors = {"A": "C3", "B": "C0", "C": "C2", "D": "0.5"}
    fig, ax = plt.subplots(4, 1, figsize=(13, 11), sharex=True)
    for i, (track, var) in enumerate([(t, v) for t in TRACKS for v in ("vv_db", "vh_db")]):
        for z, c in colors.items():
            g = series[(series.track == track) & (series.zone == z)]
            ax[i].plot(g.date, g[var], lw=0.9, color=c, label=z)
        ax[i].axvline(pd.Timestamp("2022-01-01"), color="k", ls=":", lw=0.8)
        ax[i].set_ylabel(f"{track}\n{VARS[var]}")
    ax[0].legend(ncol=4, fontsize=8)
    fig.suptitle("X-063 — zone median backscatter, 2020–2024 (dotted: start of the 2022–2024 stacks)", fontsize=10)
    fig.tight_layout(); fig.savefig(out / "fig_zone_series.png", dpi=140); plt.close(fig)

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
    s = series[~series.frozen]
    for i, var in enumerate(("vv_db", "vh_db", "rvi")):
        for z in ("A", "C", "B"):
            for track, ls in (("ascending", "-"), ("descending", "--")):
                g = s[(s.track == track) & (s.zone == z)].groupby(s.date.dt.month)[var].median()
                ax[i].plot(g.index, g.values, ls, color=colors[z], label=f"{z} {track[:4]}")
        ax[i].set_xlabel("month"); ax[i].set_title(VARS[var])
    ax[0].legend(fontsize=7, ncol=3)
    fig.suptitle("Monthly median by zone: dusk (ascending, solid) vs dawn (descending, dashed), frozen dates out",
                 fontsize=10)
    fig.tight_layout(); fig.savefig(out / "fig_monthly_dusk_dawn.png", dpi=140); plt.close(fig)

    fig, ax = plt.subplots(2, 3, figsize=(15, 9))
    outline = zones["A"].values.astype(float)
    for i, (track, ds) in enumerate(stacks.items()):
        for j, (name, f, vr) in enumerate([("σ0 VV mean (dB)", ds.gamma0_vv_db.mean("time"), (-16, -6)),
                                           ("σ0 VH mean (dB)", ds.gamma0_vh_db.mean("time"), (-22, -12)),
                                           ("RVI median", dual_pol_rvi(ds), (0.4, 1.4))]):
            f.plot(ax=ax[i, j], cmap="gray" if j < 2 else "YlGn", vmin=vr[0], vmax=vr[1])
            ax[i, j].contour(f.x, f.y, outline, levels=[0.5], colors="r", linewidths=1)
            ax[i, j].set_title(f"{track} — {name}"); ax[i, j].set_aspect("equal")
            ax[i, j].set_xlabel(""); ax[i, j].set_ylabel("")
    fig.tight_layout(); fig.savefig(out / "fig_maps.png", dpi=120); plt.close(fig)


def main(argv=None):
    hub = field.field_root().parents[1]
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--stacks", default=str(hub / "05_code" / "local" / "drive_mirror"),
                    help="folder holding rtc_dualpol_2020_2024[_descending].nc")
    ap.add_argument("--wet", default=str(hub / "08_deliverables" / "field_dew_x050" / "wetness_at_overpasses_2020_2024.csv"))
    ap.add_argument("--extra", default=str(hub / "05_code" / "local" / "s1_2020_2021"),
                    help="the 2020–2021 crops (D-020), for the clean-lake coherence")
    ap.add_argument("--reference", default=None, help="phaseDter_summary.json to reproduce (default: latest local run)")
    ap.add_argument("--out", default=str(hub / "08_deliverables" / "backscatter_x063"))
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    ctx = start("backscatter_x063", mount=False, git=False)
    zones, tpl = ctx.zones, ctx.template
    stacks = {}
    for track, name in TRACKS.items():
        ds = xr.load_dataset(Path(a.stacks) / name)
        if not (np.array_equal(ds.x, tpl.x) and np.array_equal(ds.y, tpl.y)):
            raise SystemExit(f"{name} is not on the radar grid")
        stacks[track] = ds

    # 1. reproduction
    checks = []
    canon = xr.load_dataset(ctx.paths.drive_file("rtc_dualpol_stack.nc"))
    new = stacks["ascending"].sel(time=canon.time)
    for v in ("gamma0_vv_db", "gamma0_vh_db", "ratio_vh_vv_db"):
        checks.append({"check": f"ascending stack 2022–2024 = rtc_dualpol_stack.nc, {v}", "n": int(canon.time.size),
                       "max_abs_diff": float(np.nanmax(np.abs(new[v].values - canon[v].values)))})
    ref_path = Path(a.reference) if a.reference else max(
        (Path(a.stacks) / "runs" / "phaseDter").glob("*/products/phaseDter_summary.json"), default=None)
    if ref_path and ref_path.exists():
        ref = json.loads(ref_path.read_text())
        mine = summary(years(stacks["ascending"], 2022, 2024), zones).set_index("zone")
        for r in ref["backscatter"]:
            for k in ("sigma0_vv_db", "sigma0_vv_temporal_std"):
                checks.append({"check": f"phaseDter {k}, zone {r['zone']}", "n": 1,
                               "max_abs_diff": abs(float(mine.loc[r["zone"], k]) - r[k])})
        for r in ref["amplitude_dispersion"]:
            checks.append({"check": f"phaseDter D_A median, zone {r['zone']} (n {r['n']} vs {mine.loc[r['zone'], 'n_px']})",
                           "n": 1, "max_abs_diff": abs(float(mine.loc[r["zone"], "da_median"]) - r["median"])})
    checks = pd.DataFrame(checks)

    # 2. series, with satellite and wet / frozen flags
    wet = pd.read_csv(a.wet, parse_dates=["date"])[["date", "track", "wet", "frozen"]]
    series = []
    for track, ds in stacks.items():
        s = zone_backscatter_series(ds, zones)
        sat = dict(zip(pd.to_datetime(ds.time.values), satellite(ds.source_item.values)))
        series.append(s.assign(track=track, satellite=s.date.map(sat)))
    series = pd.concat(series, ignore_index=True).merge(wet, on=["date", "track"], how="left")
    series[["wet", "frozen"]] = series[["wet", "frozen"]].fillna(False).astype(bool)

    # 3. summaries
    summ = pd.concat([summary(years(ds, y0, y1), zones).assign(track=track, period=p, n_dates=int(years(ds, y0, y1).time.size))
                      for track, ds in stacks.items() for p, (y0, y1) in PERIODS.items()], ignore_index=True)
    summ = summ[["track", "period", "n_dates", "zone", "n_px", "sigma0_vv_db", "sigma0_vv_temporal_std",
                 "ratio_vh_vv_db", "rvi_median", "da_median", f"da_frac_below_{PS_DA}"]]

    # 5. clean lake
    wc = load_worldcover(tpl, ctx.cfg, cache_dir=ctx.paths.cache)
    lake = clean_lake_mask(tpl, ctx.cfg, worldcover=wc, s2=xr.load_dataset(ctx.paths.drive_file("s2_stack.nc")))
    if ref_path and ref_path.exists():
        checks = pd.concat([checks, pd.DataFrame([{"check": "phaseDter clean-lake pixels", "n": 1,
                                                   "max_abs_diff": abs(int(lake.sum()) - ref["clean_lake_n"])}])])
    zl = {**zones, "B": lake}
    lake_bs = pd.concat([summary(years(ds, y0, y1), zl).query("zone == 'B'").assign(track=track, period=p)
                         for track, ds in stacks.items() for p, (y0, y1) in PERIODS.items()], ignore_index=True)
    lake_bs = lake_bs.assign(zone="clean lake")[["track", "period", "zone", "n_px", "sigma0_vv_db",
                                                "sigma0_vv_temporal_std", "ratio_vh_vv_db", "rvi_median", "da_median"]]
    lake_coh = lake_coherence_table(lake_coherence(ctx, Path(a.extra), zones, lake))

    # 4, 6.
    sat = satellite_offset(series)
    wd = wet_dry(series)

    series.to_csv(out / "zone_backscatter_series.csv", index=False)
    summ.to_csv(out / "zone_backscatter_summary.csv", index=False)
    sat.to_csv(out / "satellite_offset_2020_2021.csv", index=False)
    wd.to_csv(out / "wet_dry_contrast.csv", index=False)
    lake_bs.to_csv(out / "clean_lake_backscatter.csv", index=False)
    lake_coh.to_csv(out / "clean_lake_coherence.csv", index=False)
    checks.to_csv(out / "reproduction_check.csv", index=False)
    figures(out, series, stacks, zones)

    # README
    n = {t: int(ds.time.size) for t, ds in stacks.items()}
    inc = ", ".join(f"{t} {ctx.cfg['sentinel1']['tracks'][t]['incidence_angle_deg']}°" for t in TRACKS)
    span = {t: f"{pd.Timestamp(ds.time.values[0]).date()} → {pd.Timestamp(ds.time.values[-1]).date()}" for t, ds in stacks.items()}
    L = [
        "# X-063 — backscatter by zone, both tracks, 2020–2024 (exploratory)", "",
        "Generated by `05_code/SAr-adapted/scripts/backscatter_x063.py`; every number below is computed by it. "
        "phaseDter's backscatter, for the descending track and for 2020–2021, run on this Mac.", "",
        "## Inputs", "",
        f"- RTC γ⁰ stacks (Planetary Computer `sentinel-1-rtc`, `scripts/rtc_download.py`, C-048), on the 40 m radar grid, "
        f"nearest resampling like the 2022–2024 ascending stack: ascending {n['ascending']} dates ({span['ascending']}), "
        f"descending {n['descending']} dates ({span['descending']}).",
        "- Zones A (mat), B (lake), C (matched grassland), D (outside): the project's, defined once on the radar grid "
        f"(n = {', '.join(f'{z} {int(zones[z].values.sum())}' for z in 'ABCD')}).",
        "- Wet / frozen overpass flags: X-050 (station RH ≥ 95 % or rain in the previous 3 h; frozen).", "",
        "## Checks", "", md(checks.assign(max_abs_diff=checks.max_abs_diff.map(lambda v: f"{v:.2g}"))), "",
        "## Summary per track, period and zone", "",
        "σ0 VV mean and its temporal SD (phaseDter), VH − VV and RVI (volume scattering), amplitude dispersion D_A and the "
        f"share of pixels below {PS_DA} (indicative PS candidates on multi-looked RTC). Medians over the zone's pixels. "
        f"The tracks see the site at different incidence ({inc}): backscatter falls with incidence, so absolute levels "
        "differ between tracks for that reason alone — compare zones within a track, and changes in time.", "",
        md(summ), "",
        "## Clean lake — the open-water control", "",
        f"phaseDter's clean lake ({int(lake.sum())} pixels: WorldCover water or S2 NDWI > 0.2 in most scenes, inside the "
        "outline). Open water should be dark and specular in backscatter and incoherent between dates.", "",
        md(lake_bs), "",
        "Median over pairs of the zone-mean coherence (every interferogram of the period, both stacks):", "",
        md(lake_coh), "",
        "## Sentinel-1A vs 1B (2020–2021)", "",
        "Median zone value on S1A dates minus S1B dates (interleaved, same seasons; frozen dates out). An offset here "
        "would show up as a change between 2020–2021 and 2022–2024 (S1A only).", "", md(sat), "",
        "## Wet vs dry overpasses (anomalies to the calendar month)", "",
        "Median anomaly on wet dates minus dry dates; `A − C` is the mat against the grassland on the same dates. "
        "Mann–Whitney, two-sided; exploratory, many tests.", "", md(wd), "",
        "## Files", "",
        "| File | Content |", "|---|---|",
        "| `zone_backscatter_series.csv` | per track × date × zone: median VV, VH, VH − VV, RVI, valid pixels, satellite, wet, frozen |",
        "| `zone_backscatter_summary.csv` | the summary table |",
        "| `clean_lake_backscatter.csv`, `clean_lake_coherence.csv` | the open-water control |",
        "| `satellite_offset_2020_2021.csv`, `wet_dry_contrast.csv`, `reproduction_check.csv` | the tables above |",
        "| `fig_zone_series.png` | zone series, both tracks, VV and VH |",
        "| `fig_monthly_dusk_dawn.png` | monthly medians, dusk vs dawn |",
        "| `fig_maps.png` | mean VV, VH and median RVI per track, mat outlined |",
    ]
    (out / "README.md").write_text("\n".join(L) + "\n")
    pd.set_option("display.width", 220)
    print(checks.to_string(index=False))
    print(summ[summ.period == "2020–2024"].round(3).to_string(index=False))
    print(lake_bs.round(3).to_string(index=False))
    print(lake_coh.round(3).to_string(index=False))
    print(sat.round(3).to_string(index=False))
    print(wd[wd.p < 0.05].round(3).to_string(index=False))


if __name__ == "__main__":
    main()

"""The supervisor's first field-data deliverables (Sentinel-1 at the 9 monitoring plots).

Reads the Sentinel-1 products (Drive, via INSAR_DRIVE_ROOT) and the field data (the research
hub's 06_data/field, via insar_wetlands.field) and writes ONLY into the hub (private):
``--out`` (default ``<hub>/08_deliverables/field_first``). Nothing field-derived is written into
this public repository.

1. map of the plots over the Sentinel-1 grid            map_plots_s1_grid.png
2. Sentinel-1 acquisition dates 2022–2024               s1_acquisitions.csv
3. plot-level Sentinel-1 extraction                     s1_plot_by_date.csv, s1_plot_by_pair.csv,
                                                        plot_pixels.csv
   joined with WTD at the acquisitions                  plot_s1_wtd.csv, wtd_at_s1.csv
4. preliminary figure, WTD with S1 phase/coherence      fig_wtd_s1_plots.png
5. laser coverage vs Sentinel-1                         laser_coverage.csv, laser_summary.md

Windows: 1×1, 3×3 (default) and 5×5 pixels around each plot's pixel, mat pixels (zone A) only;
every statistic carries its number of valid pixels. Per-pair phase is referenced to the
grassland zone C (median over C) — the same reference as the zone aggregation.

    INSAR_DRIVE_ROOT=../local/drive_pristine PYTHONPATH=src python scripts/field_first_deliverables.py
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import xarray as xr  # noqa: E402

from insar_wetlands import field  # noqa: E402
from insar_wetlands.bootstrap import start  # noqa: E402
from insar_wetlands.inversion.isbas import PHASE_TO_MM  # noqa: E402
from insar_wetlands.stack import dates_from_pairs, list_pairs, load_layer  # noqa: E402

OVERPASS = {"ascending": "16:36", "descending": "05:09"}
WINDOWS = {"1x1": 0, "3x3": 1, "5x5": 2}
PERIOD = ("2022-01-01", "2024-12-31")


def plot_pixels(plots, template, zone_a: np.ndarray) -> pd.DataFrame:
    """One row per WTD plot: its pixel (from the mean of its polygons' centroids), mat-pixel
    counts per window, and which other plots share its 3×3 window."""
    xs, ys = template.x.values, template.y.values
    base = plots.groupby("base_plot")[["x", "y"]].mean()
    rows = []
    for p, (x, y) in base.iterrows():
        r, c = int(np.abs(ys - y).argmin()), int(np.abs(xs - x).argmin())
        rec = {"plot": p, "E": round(x, 2), "N": round(y, 2), "row": r, "col": c}
        for w, h in WINDOWS.items():
            rec[f"mat_px_{w}"] = int(zone_a[r - h:r + h + 1, c - h:c + h + 1].sum())
        rows.append(rec)
    df = pd.DataFrame(rows).sort_values("N").reset_index(drop=True)
    df["shares_3x3_with"] = [", ".join(q for q, rr, cc in zip(df["plot"], df.row, df.col)
                                       if q != p and abs(rr - r) <= 2 and abs(cc - c) <= 2)
                             for p, r, c in zip(df["plot"], df.row, df.col)]
    return df


def window(arr2d: np.ndarray, r: int, c: int, h: int, zone_a: np.ndarray) -> dict:
    sl = (slice(r - h, r + h + 1), slice(c - h, c + h + 1))
    return field.window_stats(arr2d[sl], valid=zone_a[sl])


def date_stack_rows(name: str, da: xr.DataArray, px: pd.DataFrame, zone_a, track: str) -> list[dict]:
    rows = []
    times = pd.to_datetime(da.time.values)
    vals = da.values
    for k, t in enumerate(times):
        for _, p in px.iterrows():
            for w, h in WINDOWS.items():
                rows.append({"date": t.date(), "track": track, "plot": p["plot"], "variable": name, "window": w,
                             **window(vals[k], p.row, p.col, h, zone_a)})
    return rows


def pair_rows(cropped: Path, track: str, px: pd.DataFrame, zone_a, zone_c) -> pd.DataFrame:
    pairs = list_pairs(cropped)
    corr = load_layer(cropped, "corr", pairs).values
    unw = load_layer(cropped, "unw_phase", pairs).values
    rows = []
    for k, pr in enumerate(pairs):
        a, b = pr.split("_")
        ref = unw[k][zone_c & np.isfinite(unw[k])]
        ref_med = float(np.median(ref)) if ref.size >= 5 else np.nan
        for _, p in px.iterrows():
            for w, h in WINDOWS.items():
                cs = window(corr[k], p.row, p.col, h, zone_a)
                us = window(unw[k], p.row, p.col, h, zone_a)
                d = us["median"] - ref_med
                rows.append({"pair": pr, "ref_date": pd.Timestamp(a).date(), "sec_date": pd.Timestamp(b).date(),
                             "dt_days": (pd.Timestamp(b) - pd.Timestamp(a)).days, "track": track, "plot": p["plot"],
                             "window": w, "n_valid": cs["n_valid"], "coh_median": cs["median"], "coh_iqr": cs["iqr"],
                             "dphase_vs_C_rad": d, "dlos_vs_C_mm": d * PHASE_TO_MM, "phase_iqr_rad": us["iqr"]})
    return pd.DataFrame(rows)


def coherence_by_date(pairs: pd.DataFrame, max_dt: int = 24) -> pd.DataFrame:
    """Per acquisition: mean of the window-median coherence of the short pairs that include it."""
    short = pairs[pairs.dt_days <= max_dt]
    long = pd.concat([short.rename(columns={"ref_date": "date"}), short.rename(columns={"sec_date": "date"})])
    g = long.groupby(["date", "track", "plot", "window"])
    return g.agg(coh_short_pairs=("coh_median", "mean"), n_short_pairs=("pair", "nunique"),
                 n_valid=("n_valid", "min")).reset_index()


def _runs(flag: pd.Series):
    """(start, end) of consecutive True stretches of a daily boolean series."""
    grp = (flag != flag.shift()).cumsum()
    return [(g.index[0], g.index[-1] + pd.Timedelta(days=1)) for _, g in flag.groupby(grp) if g.iloc[0]]


def main(argv=None):
    hub = field.field_root().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(hub / "08_deliverables" / "field_first"))
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    ctx = start("field_first", mount=False, git=False)
    D = ctx.paths.drive
    template, zones = ctx.template, ctx.zones
    zone_a, zone_c = zones["A"].values.astype(bool), zones["C"].values.astype(bool)
    plots = field.load_plots()
    px = plot_pixels(plots, template, zone_a)
    px.to_csv(out / "plot_pixels.csv", index=False)

    # --- 2. acquisition dates -----------------------------------------------------------
    acq = []
    crop = {"ascending": ctx.paths.cropped, "descending": ctx.paths.for_phase("field_first", track="descending").cropped}
    for track, cr in crop.items():
        for d in dates_from_pairs(list_pairs(cr)):
            if PERIOD[0] <= str(d.date()) <= PERIOD[1]:
                acq.append({"date": d.date(), "track": track, "overpass_utc": OVERPASS[track],
                            "time_utc": pd.Timestamp(f"{d.date()} {OVERPASS[track]}", tz="UTC")})
    acq = pd.DataFrame(acq).sort_values(["date", "track"])

    # --- 3. extraction ----------------------------------------------------------------------
    rows = []
    for name, f, var in (("isbas_los_mm", "ts_isbas_ref_only.nc", "los_displacement_mm"),
                         ("sbas_los_mm", "ts_sbas_ref_only.nc", "los_displacement_mm"),
                         ("evd_los_mm", "phaseE2_evd.nc", "displacement_mm"),
                         ("hybrid_vertical_mm", "hybrid_vertical_ts.nc", "vertical_mm")):
        da = xr.open_dataset(D / f)[var]
        if f == "ts_sbas_ref_only.nc":  # MintPy writes pixel corners: shift to centres (as the atlas does)
            da = da.assign_coords(x=da.x + 20.0, y=da.y - 20.0)
        rows += date_stack_rows(name, ctx.to_grid(da), px, zone_a, "ascending")
    rtc = xr.open_dataset(D / "rtc_dualpol_stack.nc")
    vv, vh = ctx.to_grid(rtc["gamma0_vv_db"]), ctx.to_grid(rtc["gamma0_vh_db"])
    lin_vv, lin_vh = 10 ** (vv / 10), 10 ** (vh / 10)
    rows += date_stack_rows("vv_db", vv, px, zone_a, "ascending")
    rows += date_stack_rows("vh_db", vh, px, zone_a, "ascending")
    rows += date_stack_rows("rvi", 4 * lin_vh / (lin_vv + lin_vh), px, zone_a, "ascending")
    by_date = pd.DataFrame(rows)
    tcoh = ctx.to_grid(xr.open_dataset(D / "phaseE2_evd.nc")["temporal_coherence"]).values
    tc = pd.DataFrame([{"plot": p["plot"], "window": w, **window(tcoh, p.row, p.col, h, zone_a)}
                       for _, p in px.iterrows() for w, h in WINDOWS.items()])
    tc.to_csv(out / "s1_plot_temporal_coherence.csv", index=False)

    pairs = pd.concat([pair_rows(crop[t], t, px, zone_a, zone_c) for t in crop])
    pairs.to_csv(out / "s1_plot_by_pair.csv", index=False)
    coh = coherence_by_date(pairs)
    by_date = pd.concat([by_date, coh.rename(columns={"coh_short_pairs": "median"}).assign(variable="coh_pairs_le24d")
                         [["date", "track", "plot", "variable", "window", "n_valid", "median"]]])
    by_date.to_csv(out / "s1_plot_by_date.csv", index=False)

    # --- WTD at the acquisitions and the joined table ---------------------------------------
    wtd = field.load_wtd_hourly()
    w_at = pd.concat([field.wtd_at_times(wtd, acq[acq.track == t].time_utc).assign(track=t) for t in crop])
    w_at["date"] = w_at.time_utc.dt.date
    cens = {p: field.censored_flag(wtd[p]) for p in field.PLOTS}
    w_at["wtd_censored"] = [bool(cens[p].asof(t)) for p, t in zip(w_at["plot"], w_at.time_utc)]
    w_at.to_csv(out / "wtd_at_s1.csv", index=False)
    wide = by_date[by_date.window == "3x3"].pivot_table(index=["date", "track", "plot"], columns="variable",
                                                         values="median").reset_index()
    nval = by_date[(by_date.window == "3x3") & (by_date.variable == "coh_pairs_le24d")][["date", "track", "plot", "n_valid"]]
    wide = wide.merge(nval, on=["date", "track", "plot"], how="left")
    wide = wide.merge(tc[tc.window == "3x3"][["plot", "median"]].rename(columns={"median": "temporal_coherence"}), on="plot")
    joined = wide.merge(w_at.drop(columns="time_utc"), on=["date", "track", "plot"], how="left")
    joined.sort_values(["plot", "track", "date"]).to_csv(out / "plot_s1_wtd.csv", index=False)
    acq["n_plots_with_s1"] = acq.apply(lambda r: int(((joined.date == r.date) & (joined.track == r.track)).sum()), axis=1)

    # --- 5. laser coverage ---------------------------------------------------------------------
    laser = field.load_laser()
    snow = field.snow_mask(wtd["Air_2m"], laser.index)
    has = laser["surface_cm"].notna()
    clean = has & ~snow
    acq["laser_within_1h"] = False
    acq["laser_snowfree_within_1h"] = False
    for track in crop:
        tt = acq.track == track
        near = [(laser.index >= t - pd.Timedelta(hours=1)) & (laser.index <= t + pd.Timedelta(hours=1)) for t in acq.time_utc[tt]]
        acq.loc[tt, "laser_within_1h"] = [bool(has[m].any()) for m in near]
        acq.loc[tt, "laser_snowfree_within_1h"] = [bool(clean[m].any()) for m in near]
    acq.drop(columns="time_utc").to_csv(out / "s1_acquisitions.csv", index=False)
    cov = pd.DataFrame({"hours": 1, "laser": has, "laser_snowfree": clean}).groupby(laser.index.to_period("M")).sum()
    cov["share_laser"] = (cov.laser / cov.hours).round(3)
    cov["share_snowfree"] = (cov.laser_snowfree / cov.hours).round(3)
    cov.to_csv(out / "laser_coverage.csv")
    s = acq.groupby("track")[["laser_within_1h", "laser_snowfree_within_1h"]].sum().astype(int)
    n = acq.groupby("track").size()
    (out / "laser_summary.md").write_text(
        "# Laser coverage vs Sentinel-1 (2022–2024)\n\n"
        f"Laser record (P6/CR): {laser.index[0]:%Y-%m-%d} → {laser.index[-1]:%Y-%m-%d}; hours with a value "
        f"{int(has.sum())} of {len(laser)} ({has.mean():.0%}); snow-free (no frost in the previous 72 h) "
        f"{int(clean.sum())} ({clean.mean():.0%}).\n\n"
        "| track | acquisitions | laser within ±1 h | snow-free laser within ±1 h |\n|---|---|---|---|\n"
        + "".join(f"| {t} | {n[t]} | {s.loc[t, 'laser_within_1h']} | {s.loc[t, 'laser_snowfree_within_1h']} |\n" for t in n.index)
        + "\nMonthly shares: `laser_coverage.csv`. Snow rule: any sub-zero hourly air temperature in the "
          "72 h before (field.snow_mask); the laser is a snow-depth sensor and reads snow as surface.\n")

    # --- 1. map ---------------------------------------------------------------------------------
    x, y = template.x.values, template.y.values
    r0, r1 = px.row.min() - 5, px.row.max() + 6
    c0, c1 = px.col.min() - 6, px.col.max() + 7
    fig, ax = plt.subplots(figsize=(6.4, 8))
    ext = [x[c0] - 20, x[c1 - 1] + 20, y[r1 - 1] - 20, y[r0] + 20]
    im = ax.imshow(tcoh[r0:r1, c0:c1], extent=ext, cmap="viridis", vmin=0.3, vmax=1.0, origin="upper")
    ax.contour(x[c0:c1], y[r0:r1], zone_a[r0:r1, c0:c1].astype(float), levels=[0.5], colors="white", linewidths=1.2)
    for xv in np.arange(x[c0] - 20, x[c1 - 1] + 21, 40):
        ax.axvline(xv, color="k", lw=0.3, alpha=0.4)
    for yv in np.arange(y[r1 - 1] - 20, y[r0] + 21, 40):
        ax.axhline(yv, color="k", lw=0.3, alpha=0.4)
    for _, p in px.iterrows():
        ax.add_patch(plt.Rectangle((x[p.col] - 20, y[p.row] - 20), 40, 40, fill=False, ec="red", lw=1.4))
    ux0, ux1 = x[px.col.min()] - 60, x[px.col.max()] + 60
    uy0, uy1 = y[px.row.max()] - 60, y[px.row.min()] + 60
    ax.add_patch(plt.Rectangle((ux0, uy0), ux1 - ux0, uy1 - uy0, fill=False, ec="orange", lw=1.2, ls="--"))
    plots.plot(ax=ax, color="red", edgecolor="red", markersize=4)
    ax.scatter(plots.x, plots.y, s=12, c="red", zorder=5)
    for _, p in px.iterrows():
        ax.annotate(p["plot"], (p.E, p.N), xytext=(8, 0), textcoords="offset points", color="white", fontsize=8,
                    fontweight="bold")
    fig.colorbar(im, ax=ax, shrink=0.6, label="temporal coherence (EVD, ascending)")
    ax.set_title("Monitoring plots on the Sentinel-1 grid (40 m pixels)\nred: plot pixel · orange dashed: extent of the 3×3 windows · white: mat edge", fontsize=9)
    ax.set_xlabel("UTM 33N easting (m)")
    ax.set_ylabel("northing (m)")
    ax.ticklabel_format(useOffset=False, style="plain")
    fig.tight_layout()
    fig.savefig(out / "map_plots_s1_grid.png", dpi=200)
    plt.close(fig)

    # --- 4. preliminary figure ---------------------------------------------------------------
    show = [p for p in ("P6", "P1", "P4", "P9") if p in set(px["plot"])]
    fig, axes = plt.subplots(len(show), 1, figsize=(11, 3.0 * len(show) + 0.8), sharex=True)
    daily = wtd[field.PLOTS].resample("D").mean()
    for ax, p in zip(np.atleast_1d(axes), show):
        ax.plot(daily.index, daily[p], color="tab:blue", lw=1, label="WTD daily mean (cm)")
        cday = cens[p].resample("D").mean() > 0.5
        for start_, end_ in _runs(cday):
            ax.axvspan(start_, end_, color="grey", alpha=0.25, lw=0)
        if p == "P6":
            ls = laser["surface_cm"].where(~snow).resample("D").mean()
            ax.plot(ls.index, ls - np.nanmedian(ls), color="black", lw=0.8, alpha=0.8, label="laser surface, snow-masked (cm, − median)")
        ax.set_ylabel("cm")
        ax.set_xlim(pd.Timestamp(PERIOD[0], tz="UTC"), pd.Timestamp(PERIOD[1], tz="UTC"))
        a2 = ax.twinx()
        for track, col in (("ascending", "tab:red"), ("descending", "tab:orange")):
            j = joined[(joined["plot"] == p) & (joined.track == track)]
            a2.plot(pd.to_datetime(j.date).dt.tz_localize("UTC"), j.coh_pairs_le24d, "o", ms=2.5, color=col,
                    label=f"coherence, pairs ≤ 24 d ({track})")
        a2.set_ylim(0, 1)
        a2.set_ylabel("coherence")
        ax.set_title(f"{p}", loc="left", fontsize=9, fontweight="bold")
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = a2.get_legend_handles_labels()
        if p == show[0]:
            h1.append(plt.Rectangle((0, 0), 1, 1, color="grey", alpha=0.25))
            l1.append("WTD probably at the sensor floor (censored)")
            fig.legend(h1 + h2, l1 + l2, fontsize=7.5, loc="upper center", bbox_to_anchor=(0.5, 0.965), ncol=3, frameon=False)
    fig.suptitle("Measured water table and Sentinel-1 coherence at the monitoring plots (3×3 windows, mat pixels)", fontsize=10, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(out / "fig_wtd_s1_plots.png", dpi=170)
    plt.close(fig)
    print(f"wrote {sorted(p.name for p in out.iterdir())}")


if __name__ == "__main__":
    main()

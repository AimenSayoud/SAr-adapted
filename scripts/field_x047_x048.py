"""X-047 / X-048: does Sentinel-1 respond to the measured water table, and does its phase see the
laser-measured surface motion at P6/CR?

Inputs: the first deliverables (hub ``08_deliverables/field_first``: ``plot_s1_wtd.csv``,
``s1_plot_by_pair.csv``) and the field data (``insar_wetlands.field``). Outputs go ONLY to the
private hub, ``--out`` (default ``<hub>/08_deliverables/field_x047_x048``).

Per acquisition: correlations of each Sentinel-1 variable with each WTD variable, raw (trend
removed) and on anomalies (trend + annual cycle removed from both), p from circular shifts.
Per interferogram (pairs ≤ 24 d): coherence (season and baseline regressed out) vs |ΔWTD|;
LOS change vs ΔWTD; at P6, LOS change vs the laser surface change projected to the LOS.
WTD at the sensor floor (``field.censored_flag``) and snow-affected laser values are excluded.

    PYTHONPATH=src python scripts/field_x047_x048.py
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
from insar_wetlands.inversion.isbas import WAVELENGTH_M  # noqa: E402

OVERPASS = {"ascending": "16:36", "descending": "05:09"}
INCIDENCE = {"ascending": 32.26, "descending": 39.17}   # config.yaml, measured (X-042)
S1_VARS = ["coh_pairs_le24d", "vv_db", "vh_db", "rvi", "isbas_los_mm", "evd_los_mm", "hybrid_vertical_mm"]
WTD_VARS = ["wtd_at", "mean_24h", "mean_prev3d", "mean_prev7d", "change_3d", "change_7d"]
QUARTER_WAVE_MM = WAVELENGTH_M / 4 * 1000


def per_date(j: pd.DataFrame) -> pd.DataFrame:
    rows = []
    j = j[~j.wtd_censored]
    for (plot, track), d in j.groupby(["plot", "track"]):
        for y in S1_VARS:
            for x in WTD_VARS:
                if d[y].notna().sum() < 12:
                    continue
                rows.append({"plot": plot, "track": track, "s1": y, "wtd": x,
                             **fl.anomaly_correlation(d.date, d[x], d[y])})
    return pd.DataFrame(rows)


def pair_table(pr: pd.DataFrame, wtd: pd.DataFrame, surface: pd.Series) -> pd.DataFrame:
    cens = {p: field.censored_flag(wtd[p]) for p in field.PLOTS}
    rows = []
    for r in pr[pr.window == "3x3"].itertuples(index=False):
        t1 = pd.Timestamp(f"{r.ref_date} {OVERPASS[r.track]}", tz="UTC")
        t2 = pd.Timestamp(f"{r.sec_date} {OVERPASS[r.track]}", tz="UTC")
        s = wtd[r.plot]
        w1, w2 = field._interp_at(s, t1), field._interp_at(s, t2)
        rec = {"pair": r.pair, "track": r.track, "plot": r.plot, "dt_days": r.dt_days, "mid": t1 + (t2 - t1) / 2,
               "coh": r.coh_median, "dlos_mm": r.dlos_vs_C_mm, "dwtd_cm": w2 - w1, "wtd_mean_cm": (w1 + w2) / 2,
               "wtd_censored": bool(cens[r.plot].asof(t1)) or bool(cens[r.plot].asof(t2))}
        if r.plot == "P6":
            rec["dsurface_cm"] = fl.value_at(surface, t2) - fl.value_at(surface, t1)
        rows.append(rec)
    return pd.DataFrame(rows)


def per_pair(pairs: pd.DataFrame, max_dt: int = 24) -> pd.DataFrame:
    rows = []
    ok = pairs[(pairs.dt_days <= max_dt) & ~pairs.wtd_censored].dropna(subset=["coh", "dlos_mm", "dwtd_cm"])
    for (plot, track), d in ok.groupby(["plot", "track"]):
        d = d.sort_values("mid")
        coh_res = fl.season_residual(d.coh, d.mid, d.dt_days)
        r_coh, p_coh = fl.circular_shift_p(np.abs(d.dwtd_cm.to_numpy()), coh_res)
        r_ph, p_ph = fl.circular_shift_p(d.dwtd_cm.to_numpy(), d.dlos_mm.to_numpy())
        rows.append({"plot": plot, "track": track, "n_pairs": len(d),
                     "r_coh_vs_abs_dwtd": r_coh, "p_coh": p_coh,
                     "r_dlos_vs_dwtd": r_ph, "p_dlos": p_ph,
                     "slope_dlos_mm_per_cm_wtd": float(np.polyfit(d.dwtd_cm, d.dlos_mm, 1)[0]),
                     "sd_dwtd_cm": float(d.dwtd_cm.std()), "sd_dlos_mm": float(d.dlos_mm.std())})
    return pd.DataFrame(rows)


def laser_p6(pairs: pd.DataFrame, wtd: pd.DataFrame, laser: pd.DataFrame, snow: pd.Series, max_dt: int = 24):
    rows = []
    q = pairs[(pairs["plot"] == "P6") & (pairs.dt_days <= max_dt) & ~pairs.wtd_censored].dropna(
        subset=["dsurface_cm", "dlos_mm", "coh", "dwtd_cm"]).sort_values("mid")
    for track, d in q.groupby("track"):
        exp = fl.los_from_vertical(d.dsurface_cm.to_numpy() * 10, INCIDENCE[track])
        r, p = fl.circular_shift_p(exp, d.dlos_mm.to_numpy())
        rs, ps = fl.circular_shift_p(d.dwtd_cm.to_numpy(), d.dsurface_cm.to_numpy())
        coh_res = fl.season_residual(d.coh, d.mid, d.dt_days)
        rc, pc = fl.circular_shift_p(np.abs(d.dsurface_cm.to_numpy()), coh_res)
        rw, pw = fl.circular_shift_p(np.abs(d.dwtd_cm.to_numpy()), coh_res)   # same pairs: WTD vs surface
        rows.append({
            "track": track, "n_pairs": len(d),
            "median_abs_dsurface_cm": float(d.dsurface_cm.abs().median()), "max_abs_dsurface_cm": float(d.dsurface_cm.abs().max()),
            "sd_expected_los_mm": float(exp.std()), "sd_observed_los_mm": float(d.dlos_mm.std()),
            "share_pairs_expected_los_over_quarter_wave": float(np.mean(np.abs(exp) > QUARTER_WAVE_MM)),
            "r_los_vs_laser": r, "p_los_vs_laser": p,
            "slope_los_per_expected_los": float(np.polyfit(exp, d.dlos_mm, 1)[0]),
            "r_dsurface_vs_dwtd": rs, "p_dsurface_vs_dwtd": ps,
            "slope_dsurface_per_dwtd": float(np.polyfit(d.dwtd_cm, d.dsurface_cm, 1)[0]),
            "r_coh_vs_abs_dsurface": rc, "p_coh_vs_abs_dsurface": pc,
            "r_coh_vs_abs_dwtd_same_pairs": rw, "p_coh_vs_abs_dwtd_same_pairs": pw})
    daily = pd.DataFrame({"surface": laser.surface_cm.where(~snow), "cr_raw": wtd["CR_raw"].reindex(laser.index),
                          "wtd_p6": wtd["P6"].reindex(laser.index)}).resample("D").mean().dropna()
    buoy = {"n_days": len(daily), "r_surface_vs_raw_level": float(daily.surface.corr(daily.cr_raw)),
            "slope_surface_per_raw_level": float(np.polyfit(daily.cr_raw, daily.surface, 1)[0]),
            "r_surface_vs_wtd_p6": float(daily.surface.corr(daily.wtd_p6))}
    return pd.DataFrame(rows), buoy, q, daily


def figures(out: Path, pd_corr: pd.DataFrame, pairs: pd.DataFrame, q: pd.DataFrame, daily: pd.DataFrame):
    # 1. heatmap of anomaly correlations with the water table at the overpass
    sub = pd_corr[pd_corr.wtd == "wtd_at"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.2), sharey=True)
    for ax, track in zip(axes, ["ascending", "descending"]):
        m = sub[sub.track == track].pivot(index="s1", columns="plot", values="r_anom").reindex(S1_VARS)
        pm = sub[sub.track == track].pivot(index="s1", columns="plot", values="p_anom").reindex(S1_VARS)
        im = ax.imshow(m.values, cmap="RdBu_r", vmin=-0.7, vmax=0.7, aspect="auto")
        for i in range(m.shape[0]):
            for k in range(m.shape[1]):
                v = m.values[i, k]
                if np.isfinite(v):
                    ax.text(k, i, f"{v:.2f}" + ("*" if pm.values[i, k] < 0.05 else ""), ha="center", va="center", fontsize=7)
        ax.set_xticks(range(m.shape[1]), m.columns)
        ax.set_yticks(range(m.shape[0]), m.index)
        empty = m.isna().all(axis=1)
        for i in np.where(empty.to_numpy())[0]:
            ax.text(m.shape[1] / 2 - 0.5, i, "no product for this track", ha="center", va="center", fontsize=7, color="grey")
        ax.set_title(track, fontsize=10)
    fig.suptitle("Correlation of anomalies (annual cycle removed from both) with the WTD at the overpass — * p < 0.05 (circular shifts)", fontsize=10)
    fig.colorbar(im, ax=axes, shrink=0.8, label="r")
    fig.savefig(out / "fig1_anomaly_correlations.png", dpi=170, bbox_inches="tight")
    plt.close(fig)

    # 2. coherence (season/baseline removed) vs |ΔWTD| per pair
    ok = pairs[(pairs.dt_days <= 24) & ~pairs.wtd_censored].dropna(subset=["coh", "dwtd_cm"])
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
    for ax, track in zip(axes, ["ascending", "descending"]):
        for plot, d in ok[ok.track == track].groupby("plot"):
            d = d.sort_values("mid")
            ax.scatter(np.abs(d.dwtd_cm), fl.season_residual(d.coh, d.mid, d.dt_days), s=6, alpha=0.5, label=plot)
        ax.axhline(0, color="k", lw=0.5)
        ax.set_xlabel("|ΔWTD| over the pair (cm)")
        ax.set_title(f"{track}, pairs ≤ 24 d", fontsize=9)
    axes[0].set_ylabel("coherence, season & baseline removed")
    axes[1].legend(fontsize=7, ncol=3)
    fig.tight_layout()
    fig.savefig(out / "fig2_coherence_vs_wtd_change.png", dpi=170)
    plt.close(fig)

    # 3. P6: laser vs radar
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.3))
    for track, col in (("ascending", "tab:red"), ("descending", "tab:orange")):
        d = q[q.track == track]
        exp = fl.los_from_vertical(d.dsurface_cm * 10, INCIDENCE[track])
        axes[0].scatter(exp, d.dlos_mm, s=10, color=col, alpha=0.7, label=track)
        axes[1].scatter(d.dwtd_cm, d.dsurface_cm, s=10, color=col, alpha=0.7, label=track)
    lim = 45
    axes[0].plot([-lim, lim], [-lim, lim], "k--", lw=0.8, label="1:1 (radar sees the surface)")
    for s_ in (-1, 1):
        axes[0].axvline(s_ * QUARTER_WAVE_MM, color="grey", lw=0.6, ls=":")
    axes[0].set_xlabel("laser surface change projected to LOS (mm)")
    axes[0].set_ylabel("Sentinel-1 LOS change, P6 3×3 − grassland (mm)")
    axes[0].set_title("P6: does the phase follow the surface? (dotted: ±λ/4)", fontsize=9)
    axes[0].legend(fontsize=7)
    axes[1].set_xlabel("ΔWTD over the pair (cm)")
    axes[1].set_ylabel("laser surface change (cm)")
    axes[1].set_title("P6: surface change vs water-table change", fontsize=9)
    axes[2].plot(daily.index, daily.cr_raw - daily.cr_raw.median(), lw=0.9, label="raw water level CR (− median)")
    axes[2].plot(daily.index, daily.surface - daily.surface.median(), lw=0.9, color="k", label="laser surface, snow-free (− median)")
    axes[2].set_ylabel("cm")
    axes[2].set_title("P6: the mat surface rides with the water level (daily)", fontsize=9)
    axes[2].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(out / "fig3_p6_laser_vs_radar.png", dpi=170)
    plt.close(fig)


def main(argv=None):
    hub = field.field_root().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument("--first", default=str(hub / "08_deliverables" / "field_first"))
    ap.add_argument("--out", default=str(hub / "08_deliverables" / "field_x047_x048"))
    a = ap.parse_args(argv)
    first, out = Path(a.first), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    j = pd.read_csv(first / "plot_s1_wtd.csv", parse_dates=["date"])
    pr = pd.read_csv(first / "s1_plot_by_pair.csv")
    wtd = field.load_wtd_hourly()
    laser = field.load_laser()
    snow = field.snow_mask(wtd["Air_2m"], laser.index)

    pdc = per_date(j)
    pdc.to_csv(out / "per_date_correlations.csv", index=False)
    pairs = pair_table(pr, wtd, laser.surface_cm.where(~snow))
    pairs.to_csv(out / "per_pair_table.csv", index=False)
    ppc = per_pair(pairs)
    ppc.to_csv(out / "per_pair_correlations.csv", index=False)
    lz, buoy, q, daily = laser_p6(pairs, wtd, laser, snow)
    lz.to_csv(out / "p6_laser_vs_radar.csv", index=False)
    pd.Series(buoy).to_csv(out / "p6_surface_vs_water_level.csv", header=["value"])
    figures(out, pdc, pairs, q, daily)
    print(f"wrote {sorted(p.name for p in out.iterdir())}")


if __name__ == "__main__":
    main()

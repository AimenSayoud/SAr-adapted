"""X-050: does surface wetness at the overpass (dew, rain) explain why the descending (05:09 UTC)
phase carries no water-table / surface signal while the ascending (16:36 UTC) one does?

The water table and the mat surface are the same at dawn and dusk (hourly field data), so a
dusk/dawn difference must come from the surface state. The station's humidity and rain give a
``wet`` flag per acquisition (``field.surface_wetness_at``); frozen acquisitions (air ≤ 0 °C) are
a different dielectric state and are left out. Three tests, each per track:

1. Zone coherence (phaseD, every pair ≤ 48 d): season- and baseline-cleaned coherence of pairs
   with a wet date vs pairs with none — on the mat (A), the lake (B), the grassland (C) and the
   mat minus grassland (A−C).
2. Aggregated A−C series (the asc/desc consistency run): |residual from the seasonal fit| on wet
   vs dry dates; B−C and the adjacent null as controls.
3. Plot phase (X-047/X-048 pair table, pairs ≤ 24 d): phase vs ΔWTD per plot and phase vs laser
   at P6, on pairs whose two dates are both dry vs the rest.

Significance of 1–2: ``field_link.flag_shift_test`` (wet flags circularly shifted along the
dates; keeps their seasonality). Of 3: ``circular_shift_p``. Outputs go ONLY to the private hub
(``--out``, default ``<hub>/08_deliverables/field_dew_x050``): the meteo is field data.

    PYTHONPATH=src python scripts/field_dew_x050.py
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

from insar_wetlands import field  # noqa: E402
from insar_wetlands import field_link as fl  # noqa: E402
from insar_wetlands.hydro_link import _detrend  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
OVERPASS = {"ascending": "16:36", "descending": "05:09"}
INCIDENCE = {"ascending": 32.26, "descending": 39.17}   # config.yaml, measured (X-042)
COH_FILES = {"ascending": "phaseD_coh_by_zone.csv", "descending": "phaseD_coh_by_zone_descending.csv"}
RH_THRESHOLDS = (90.0, 95.0, 98.0)
MAIN_RH = 95.0


def iso(d: str) -> str:
    """``20220110`` → ``2022-01-10`` (pair names use the compact form)."""
    return f"{d[:4]}-{d[4:6]}-{d[6:]}" if "-" not in d else d


def acquisition_dates(coh: pd.DataFrame) -> list[str]:
    return sorted({iso(d) for p in coh.pair for d in p.split("_")})


def wetness_table(meteo: pd.DataFrame, dates: dict[str, list[str]]) -> pd.DataFrame:
    rows = []
    for track, ds in dates.items():
        t = [pd.Timestamp(f"{d} {OVERPASS[track]}", tz="UTC") for d in ds]
        w = field.surface_wetness_at(meteo, t, rh_wet=MAIN_RH)
        for thr in RH_THRESHOLDS:
            w[f"wet_rh{int(thr)}"] = (w.rh >= thr) | (w.rain_prev_mm > 0)
        w.insert(0, "date", ds)
        w.insert(1, "track", track)
        rows.append(w)
    return pd.concat(rows, ignore_index=True)


def _members(pairs, index):
    return np.array([[index[iso(a)], index[iso(b)]] for a, b in (p.split("_") for p in pairs)])


def zone_coherence(coh: pd.DataFrame, wt: pd.DataFrame, track: str, max_dt: int = 48) -> list[dict]:
    w = wt[wt.track == track].reset_index(drop=True)
    index = {d: i for i, d in enumerate(w.date)}
    frozen = w.frozen.to_numpy()
    wide = coh[coh.dt_days <= max_dt].pivot_table(index=["pair", "dt_days"], columns="zone", values="mean_coh").reset_index()
    m = _members(wide.pair, index)
    keep = ~frozen[m].any(axis=1)
    wide, m = wide[keep].reset_index(drop=True), m[keep]
    mid = [pd.Timestamp(a) + (pd.Timestamp(b) - pd.Timestamp(a)) / 2 for a, b in (p.split("_") for p in wide.pair)]
    wide["A-C"] = wide["A"] - wide["C"]
    rows = []
    for zone in ("A", "B", "C", "A-C"):
        res = fl.season_residual(wide[zone], mid, wide.dt_days)
        for thr in RH_THRESHOLDS:
            flags = w[f"wet_rh{int(thr)}"].to_numpy()
            r = fl.flag_shift_test(res, m, flags)
            hit = flags[m].any(axis=1)
            rows.append({"track": track, "zone": zone, "rh_threshold": thr, "n_dates_wet": int(flags[~frozen].sum()),
                         "n_dates_used": int((~frozen).sum()), **r,
                         "mean_coh_clear": float(wide[zone][~hit].mean()), "mean_coh_flagged": float(wide[zone][hit].mean())})
    return rows


def series_residuals(series: dict, wt: pd.DataFrame, track: str) -> list[dict]:
    w = wt[wt.track == track].set_index("date")
    rows = []
    for label in ("A-C", "B-C", "A-B", "null"):
        s = series[f"{track}|{label}|full"]
        d = pd.DataFrame({"date": s["date"], "disp": s["disp_mm"]}).dropna()
        d = d[d.date.map(lambda x: not w.frozen.get(x, True))].reset_index(drop=True)
        t = fl.years_since(d.date)
        res = np.abs(_detrend(d.disp.to_numpy(float), t, True))
        for thr in RH_THRESHOLDS:
            flags = w.loc[d.date, f"wet_rh{int(thr)}"].to_numpy(bool)
            r = fl.flag_shift_test(res, np.arange(len(d))[:, None], flags)
            rows.append({"track": track, "series": "adjacent null" if label == "null" else label,
                         "rh_threshold": thr, "n_dates": len(d),
                         "median_abs_res_dry_mm": float(np.median(res[~flags])),
                         "median_abs_res_wet_mm": float(np.median(res[flags])) if flags.any() else np.nan,
                         "diff_mean_dry_minus_wet_mm": r["diff"], "p": r["p"]})
    return rows


def split_stat(plots: np.ndarray, x: np.ndarray, y: np.ndarray, wet_any: np.ndarray) -> float:
    """Mean over plots of r(x, y | both dates dry) − r(x, y | a wet date); plots with < 12
    pairs in either group are skipped."""
    diffs = []
    for pl in field.PLOTS:
        a, b = (plots == pl) & ~wet_any, (plots == pl) & wet_any
        if a.sum() >= 12 and b.sum() >= 12:
            diffs.append(np.corrcoef(x[a], y[a])[0, 1] - np.corrcoef(x[b], y[b])[0, 1])
    return float(np.mean(diffs)) if diffs else np.nan


def plot_phase(pairs: pd.DataFrame, wt: pd.DataFrame, max_dt: int = 24, n_shift: int = 1000
               ) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Phase vs ΔWTD per plot, and phase vs laser at P6, on pairs with two dry dates vs pairs
    with a wet one (frozen pairs out). The split test: mean over plots of r(dry) − r(wet),
    against the same with the wet flags circularly shifted along each track's dates."""
    q = pairs[(pairs.dt_days <= max_dt) & ~pairs.wtd_censored].dropna(subset=["dlos_mm", "dwtd_cm"]).copy()
    rows, lz, split = [], [], []
    for track in OVERPASS:
        w = wt[wt.track == track].reset_index(drop=True)
        index = {d: i for i, d in enumerate(w.date)}
        d_all = q[q.track == track].sort_values("mid").reset_index(drop=True)
        m = _members(d_all.pair, index)
        keep = ~w.frozen.to_numpy()[m].any(axis=1)
        d_all, m = d_all[keep].reset_index(drop=True), m[keep]
        plots = d_all["plot"].to_numpy()
        x, y = d_all.dwtd_cm.to_numpy(), d_all.dlos_mm.to_numpy()

        for thr in RH_THRESHOLDS:
            flags = w[f"wet_rh{int(thr)}"].to_numpy(bool)
            wet_any = flags[m].any(axis=1)
            stat = split_stat(plots, x, y, wet_any)
            rng = np.random.default_rng(0)
            null = np.array([split_stat(plots, x, y, np.roll(flags, k)[m].any(axis=1)) for k in rng.integers(3, len(flags) - 3, n_shift)])
            null = null[np.isfinite(null)]
            split.append({"track": track, "rh_threshold": thr, "mean_r_dry_minus_r_wet": stat,
                          "p": float((np.sum(np.abs(null) >= abs(stat)) + 1) / (len(null) + 1)),
                          "n_pairs_dry": int((~wet_any).sum()), "n_pairs_wet": int(wet_any.sum())})
            for pl in field.PLOTS:
                for st, sel in (("both dry", ~wet_any), ("any wet", wet_any)):
                    k = (plots == pl) & sel
                    r, pv = fl.circular_shift_p(x[k], y[k]) if k.sum() >= 12 else (np.nan, np.nan)
                    rows.append({"plot": pl, "track": track, "rh_threshold": thr, "surface_state": st, "n_pairs": int(k.sum()),
                                 "r_dlos_vs_dwtd": r, "p": pv,
                                 "slope_mm_per_cm": float(np.polyfit(x[k], y[k], 1)[0]) if k.sum() >= 3 else np.nan})
            p6 = (plots == "P6") & d_all.dsurface_cm.notna().to_numpy()
            for st, sel in (("both dry", ~wet_any), ("any wet", wet_any)):
                k = p6 & sel
                exp = fl.los_from_vertical(d_all.dsurface_cm.to_numpy()[k] * 10, INCIDENCE[track])
                r, pv = fl.circular_shift_p(exp, y[k]) if k.sum() >= 12 else (np.nan, np.nan)
                lz.append({"track": track, "rh_threshold": thr, "surface_state": st, "n_pairs": int(k.sum()),
                           "r_los_vs_laser": r, "p": pv,
                           "slope_los_per_expected_los": float(np.polyfit(exp, y[k], 1)[0]) if k.sum() >= 3 else np.nan})
    return pd.DataFrame(rows), pd.DataFrame(lz), pd.DataFrame(split)


def figure(out: Path, wt: pd.DataFrame, zc: pd.DataFrame, sr: pd.DataFrame, pp: pd.DataFrame, split: pd.DataFrame):
    fig, ax = plt.subplots(1, 4, figsize=(20, 4.4))
    for track, c in (("ascending", "#d97706"), ("descending", "#2563eb")):
        w = wt[wt.track == track]
        ax[0].hist(w.rh, bins=np.arange(30, 101, 2.5), alpha=0.55, color=c, label=f"{track} ({OVERPASS[track]} UTC)")
    ax[0].axvline(MAIN_RH, ls="--", c="k", lw=0.8)
    ax[0].set_xlabel("relative humidity at the overpass (%)"); ax[0].set_ylabel("acquisitions"); ax[0].legend(fontsize=8)
    ax[0].set_title("How often the surface is likely wet", fontsize=9)
    z = zc[zc.rh_threshold == MAIN_RH]
    x = np.arange(4)
    for k, (track, c) in enumerate((("ascending", "#d97706"), ("descending", "#2563eb"))):
        d = z[z.track == track].set_index("zone").loc[["A", "B", "C", "A-C"]]
        ax[1].bar(x + (k - 0.5) * 0.38, d["diff"], 0.38, color=c, label=track)
        for xi, (v, p) in enumerate(zip(d["diff"], d["p"])):
            ax[1].text(xi + (k - 0.5) * 0.38, v, "*" if p < 0.05 else "", ha="center", va="bottom")
    ax[1].axhline(0, c="k", lw=0.6); ax[1].set_xticks(x, ["A mat", "B lake", "C grass", "A−C"])
    ax[1].set_ylabel("coherence: pairs without − with a wet date"); ax[1].legend(fontsize=8)
    ax[1].set_title(f"Zone coherence, season removed (RH ≥ {MAIN_RH:.0f} % or rain; * p<0.05)", fontsize=9)
    s = sr[sr.rh_threshold == MAIN_RH]
    x = np.arange(4)
    for k, (track, c) in enumerate((("ascending", "#d97706"), ("descending", "#2563eb"))):
        d = s[s.track == track].set_index("series").loc[["A-C", "B-C", "A-B", "adjacent null"]]
        ax[2].bar(x + (k - 0.5) * 0.38, d["diff_mean_dry_minus_wet_mm"], 0.38, color=c, label=track)
        for xi, (v, p) in enumerate(zip(d["diff_mean_dry_minus_wet_mm"], d["p"])):
            ax[2].text(xi + (k - 0.5) * 0.38, v, "*" if p < 0.05 else "", ha="center", va="bottom")
    ax[2].axhline(0, c="k", lw=0.6); ax[2].set_xticks(x, ["A−C", "B−C", "A−B", "null"])
    ax[2].set_ylabel("|residual| dry − wet dates (mm)"); ax[2].legend(fontsize=8)
    ax[2].set_title("Aggregated series: |residual| from the seasonal fit\n(negative = wet dates noisier)", fontsize=9)
    ph = pp[pp.rh_threshold == MAIN_RH]
    for k, (track, c) in enumerate((("ascending", "#d97706"), ("descending", "#2563eb"))):
        for st, mk, off in (("both dry", "o", -0.12), ("any wet", "x", 0.12)):
            d = ph[(ph.track == track) & (ph.surface_state == st)].set_index("plot").reindex(field.PLOTS)
            ax[3].scatter(np.arange(9) + off + (k - 0.5) * 0.02, d.r_dlos_vs_dwtd, marker=mk, color=c, s=28,
                          label=f"{track}, {st}")
        sp = split[(split.track == track) & (split.rh_threshold == MAIN_RH)].iloc[0]
        ax[3].text(0.01, 0.97 - 0.07 * k, f"{track}: mean r(dry) − r(wet) = {sp.mean_r_dry_minus_r_wet:.2f}, p = {sp.p:.2f}",
                   transform=ax[3].transAxes, fontsize=7.5, color=c, va="top")
    ax[3].axhline(0, c="k", lw=0.6); ax[3].set_xticks(np.arange(9), field.PLOTS)
    ax[3].set_ylabel("r (LOS change vs ΔWTD)"); ax[3].legend(fontsize=7, loc="lower right")
    ax[3].set_title("Plot phase vs water table: pairs with two dry dates vs a wet one\n(p: wet flags shifted along the dates)", fontsize=9)
    fig.tight_layout()
    fig.savefig(out / "fig_dew_test.png", dpi=150)
    plt.close(fig)


def main(argv=None):
    hub = field.field_root().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument("--series", default=str(hub / "05_code" / "local" / "results" / "asc_desc" / "asc_desc_web.json"))
    ap.add_argument("--pairs", default=str(hub / "08_deliverables" / "field_x047_x048" / "per_pair_table.csv"))
    ap.add_argument("--out", default=str(hub / "08_deliverables" / "field_dew_x050"))
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    meteo = field.load_wtd_hourly()
    coh = {t: pd.read_csv(REPO / "results" / "tables" / f) for t, f in COH_FILES.items()}
    wt = wetness_table(meteo, {t: acquisition_dates(c) for t, c in coh.items()})
    wt.to_csv(out / "wetness_at_overpasses.csv", index=False)

    zc = pd.DataFrame([r for t in coh for r in zone_coherence(coh[t], wt, t)])
    zc.to_csv(out / "zone_coherence_by_wetness.csv", index=False)
    series = json.loads(Path(a.series).read_text())["series"]
    sr = pd.DataFrame([r for t in coh for r in series_residuals(series, wt, t)])
    sr.to_csv(out / "series_residual_by_wetness.csv", index=False)
    pp, lz, split = plot_phase(pd.read_csv(a.pairs), wt)
    pp.to_csv(out / "plot_phase_by_wetness.csv", index=False)
    split.to_csv(out / "plot_phase_split_test.csv", index=False)
    lz.to_csv(out / "p6_laser_by_wetness.csv", index=False)
    figure(out, wt, zc, sr, pp, split)
    print(f"wrote {sorted(p.name for p in out.iterdir())}")


if __name__ == "__main__":
    main()

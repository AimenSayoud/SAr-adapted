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
4. Zone phase (``--drive``: a read-only Drive snapshot): the coherence-weighted A−C phase of every
   pair (``aggregate_unwrapped``, as phaseG) vs the median ΔWTD over the plots and vs the P6 laser,
   dry-dry pairs vs the rest — one test on the whole mat instead of nine overlapping plots.
   B−C, A−B and the adjacent-null zone are the controls.

Significance of 1–2: ``field_link.flag_shift_test`` (wet flags circularly shifted along the
dates; keeps their seasonality). Of 3: ``circular_shift_p``. Outputs go ONLY to the private hub
(``--out``, default ``<hub>/08_deliverables/field_dew_x050``): the meteo is field data.

    PYTHONPATH=src python scripts/field_dew_x050.py [--drive <snapshot>]
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


def zone_pair_table(drive: Path, wt: pd.DataFrame, extra: dict[str, Path] | None = None) -> pd.DataFrame:
    """Per track × zone pair × interferogram: aggregated LOS change, median ΔWTD over the plots
    (dry-well stretches out; ≥ 5 plots), the P6 snow-free laser change, wet/frozen flags.
    ``extra``: per track, another cropped stack on the same grid (the 2020–2021 extension, D-020)
    aggregated with the same zones; ``stack`` says which one a pair comes from."""
    from insar_wetlands.aggregate import adjacent_null_zones, aggregate_unwrapped
    from insar_wetlands.bootstrap import start
    from insar_wetlands.stack import list_pairs, load_layer
    ctx = start("field_dew_x050", mount=False, git=False, drive_root=drive)
    zones, tpl = ctx.zones, ctx.template
    wtd = field.load_wtd_hourly()
    cens = {p: field.censored_flag(wtd[p]) for p in field.PLOTS}
    laser = field.load_laser()
    surface = laser.surface_cm.where(~field.snow_mask(wtd["Air_2m"], laser.index))
    znull = adjacent_null_zones(zones, tpl, ref="D")
    key = wt.set_index(["track", "date"])
    rows = []
    stacks = [(track, "2022–2024", ctx.paths.for_phase("x", track=track).cropped) for track in OVERPASS]
    stacks += [(track, "2020–2021", Path(d)) for track, d in (extra or {}).items()]
    for track, stack, cropped in stacks:
        pairs = list_pairs(cropped)
        unw, corr = load_layer(cropped, "unw_phase", pairs), load_layer(cropped, "corr", pairs)
        dd = {label: aggregate_unwrapped(unw, corr, z, target=t, reference=r).set_index("pair")
              for label, z, t, r in (("A-C", zones, "A", "C"), ("B-C", zones, "B", "C"),
                                     ("A-B", zones, "A", "B"), ("adjacent null", znull, "A", "C"))}
        del unw, corr
        for pair in pairs:
            a, b = (iso(x) for x in pair.split("_"))
            t1, t2 = (pd.Timestamp(f"{x} {OVERPASS[track]}", tz="UTC") for x in (a, b))
            dw = [field._interp_at(wtd[p], t2) - field._interp_at(wtd[p], t1) for p in field.PLOTS
                  if not (bool(cens[p].asof(t1)) or bool(cens[p].asof(t2)))]
            s1, s2 = key.loc[(track, a)], key.loc[(track, b)]
            base = {"track": track, "stack": stack, "pair": pair, "date1": a, "date2": b, "dt_days": (pd.Timestamp(b) - pd.Timestamp(a)).days,
                    "dwtd_median_cm": float(np.median(dw)) if len(dw) >= 5 else np.nan, "n_plots": len(dw),
                    "dsurface_p6_cm": fl.value_at(surface, t2) - fl.value_at(surface, t1),
                    "frozen_any": bool(s1.frozen or s2.frozen)}
            for thr in RH_THRESHOLDS:
                base[f"wet_any_rh{int(thr)}"] = bool(s1[f"wet_rh{int(thr)}"] or s2[f"wet_rh{int(thr)}"])
            for label, d in dd.items():
                if pair in d.index:
                    rows.append({**base, "zones": label, "ddisp_mm": float(d.loc[pair, "ddisp_mm"]),
                                 "weight": float(d.loc[pair, "weight"])})
    return pd.DataFrame(rows)


def zone_phase_tests(zt: pd.DataFrame, wt: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (track, label), d0 in zt.groupby(["track", "zones"]):
        w = wt[wt.track == track].reset_index(drop=True)
        index = {x: i for i, x in enumerate(w.date)}
        for max_dt in (24, 48):
            d = d0[(d0.dt_days <= max_dt) & ~d0.frozen_any].sort_values("date1")
            for xvar in ("dwtd_median_cm", "dsurface_p6_cm"):
                q = d.dropna(subset=[xvar, "ddisp_mm"])
                if len(q) < 24:
                    continue
                m = np.array([[index[a], index[b]] for a, b in zip(q.date1, q.date2)])
                x = q[xvar].to_numpy()
                if xvar == "dsurface_p6_cm":
                    x = fl.los_from_vertical(x * 10, INCIDENCE[track])
                r_all, p_all = fl.circular_shift_p(x, q.ddisp_mm.to_numpy())
                for thr in RH_THRESHOLDS:
                    s = fl.flag_split_correlation(x, q.ddisp_mm.to_numpy(), m, w[f"wet_rh{int(thr)}"].to_numpy(bool))
                    rows.append({"track": track, "zones": label, "max_dt": max_dt, "x": xvar, "rh_threshold": thr,
                                 "n_pairs": len(q), "r_all": r_all, "p_all": p_all, **s})
    return pd.DataFrame(rows)


def zone_phase_by_period(zt: pd.DataFrame, wt: pd.DataFrame) -> pd.DataFrame:
    """``zone_phase_tests`` on 2020–2021, 2022–2024 and both (periods by the pair's first date;
    the wet flags are shifted along all 2020–2024 dates of the track)."""
    first = pd.to_datetime(zt.date1)
    out = []
    for period, sel in (("2020–2021", first < "2022-01-01"), ("2022–2024", first >= "2022-01-01"),
                        ("2020–2024", first.notna())):
        out.append(zone_phase_tests(zt[sel], wt).assign(period=period))
    return pd.concat(out, ignore_index=True)


def dry_pair_tests(zt: pd.DataFrame, max_dt: int = 24) -> pd.DataFrame:
    """Does the zone phase follow the water table on pairs whose two dates are both dry (RH < 95 %,
    no rain, not frozen)? r and circular-shift p per period × track × zone pair — the question the
    dry/wet contrast cannot answer when one group is small."""
    d = zt[(zt.dt_days <= max_dt) & ~zt.frozen_any & ~zt[f"wet_any_rh{int(MAIN_RH)}"]].dropna(
        subset=["dwtd_median_cm", "ddisp_mm"])
    first = pd.to_datetime(d.date1)
    rows = []
    for period, sel in (("2020–2021", first < "2022-01-01"), ("2022–2024", first >= "2022-01-01"),
                        ("2020–2024", first.notna())):
        for (track, label), q in d[sel].groupby(["track", "zones"]):
            q = q.sort_values("date1")
            r, p = fl.circular_shift_p(q.dwtd_median_cm.to_numpy(), q.ddisp_mm.to_numpy())
            rows.append({"period": period, "track": track, "zones": label, "n_dry_pairs": len(q), "r": r, "p": p,
                         "slope_mm_per_cm": float(np.polyfit(q.dwtd_median_cm, q.ddisp_mm, 1)[0])})
    return pd.DataFrame(rows)


def dry_figure(out: Path, dp: pd.DataFrame):
    labels = ["A-C", "B-C", "A-B", "adjacent null"]
    periods = ["2020–2021", "2022–2024", "2020–2024"]
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.2), sharey=True)
    for a, (track, c) in zip(ax, (("ascending", "#d97706"), ("descending", "#2563eb"))):
        for k, per in enumerate(periods):
            t = dp[(dp.track == track) & (dp.period == per)].set_index("zones").reindex(labels)
            xs = np.arange(4) + (k - 1) * 0.26
            a.bar(xs, t.r, 0.24, color=c, alpha=(0.45, 0.7, 1.0)[k], label=per)
            for x, r, p, n in zip(xs, t.r, t.p, t.n_dry_pairs):
                a.text(x, (r if r > 0 else 0) + 0.015, f"{'*' if p < 0.05 else ''}\n{n}", ha="center", fontsize=6.5)
        a.axhline(0, c="k", lw=0.6); a.set_xticks(np.arange(4), ["A−C mat−grass", "B−C lake−grass", "A−B mat−lake", "null"])
        a.set_title(f"{track}: zone phase vs ΔWTD on dry-dry pairs ≤ 24 d\n(* p < 0.05, circular shifts; n = pairs)", fontsize=9)
        a.legend(fontsize=7)
    ax[0].set_ylabel("r")
    fig.tight_layout()
    fig.savefig(out / "fig_zone_phase_dry_2020_2024.png", dpi=150)
    plt.close(fig)


def zone_figure(out: Path, zp: pd.DataFrame):
    labels = ["A-C", "B-C", "A-B", "adjacent null"]
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.2), sharey=True)
    for a, xvar, title in ((ax[0], "dwtd_median_cm", "vs median ΔWTD over the plots"),
                           (ax[1], "dsurface_p6_cm", "vs P6 laser surface change (LOS)")):
        d = zp[(zp.x == xvar) & (zp.max_dt == 24) & (zp.rh_threshold == MAIN_RH)]
        for k, (track, c) in enumerate((("ascending", "#d97706"), ("descending", "#2563eb"))):
            t = d[d.track == track].set_index("zones").reindex(labels)
            xs = np.arange(4) + (k - 0.5) * 0.36
            a.bar(xs - 0.08, t.r_clear, 0.16, color=c, label=f"{track}: both dates dry")
            a.bar(xs + 0.08, t.r_flagged, 0.16, color=c, alpha=0.35, hatch="//", label=f"{track}: a wet date")
            for xi, (hi, p) in enumerate(zip(np.fmax(t.r_clear, t.r_flagged), t.p)):
                a.text(xs[xi], hi + 0.02, f"p={p:.2f}", ha="center", fontsize=7, color=c)
        a.axhline(0, c="k", lw=0.6); a.set_xticks(np.arange(4), ["A−C mat−grass", "B−C lake−grass", "A−B mat−lake", "null"])
        a.set_title(f"Aggregated zone phase {title}\n(pairs ≤ 24 d; p: dry − wet, flags shifted)", fontsize=9)
    ax[0].set_ylabel("r"); ax[0].legend(fontsize=7, loc="lower left")
    fig.tight_layout()
    fig.savefig(out / "fig_zone_phase.png", dpi=150)
    plt.close(fig)


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
    ap.add_argument("--drive", default=None, help="read-only Drive snapshot for the zone-phase test (4)")
    ap.add_argument("--extra", nargs="*", default=[], metavar="TRACK=DIR",
                    help="extra cropped stacks for test 4 (D-020: 2020–2021), e.g. ascending=<dir>")
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
    if a.drive:
        zt = zone_pair_table(Path(a.drive), wt)
        zt.to_csv(out / "zone_pair_table.csv", index=False)
        zp = zone_phase_tests(zt, wt)
        zp.to_csv(out / "zone_phase_by_wetness.csv", index=False)
        zone_figure(out, zp)
    if a.drive and a.extra:   # 4 again on 2020–2024 (D-020); the 2022–2024 outputs above are unchanged
        from insar_wetlands.stack import list_pairs
        extra = dict(x.split("=", 1) for x in a.extra)
        dates = {t: sorted(set(acquisition_dates(coh[t])) | {iso(d) for p in list_pairs(Path(extra[t])) for d in p.split("_")})
                 for t in OVERPASS}
        wt_all = wetness_table(meteo, dates)
        wt_all.to_csv(out / "wetness_at_overpasses_2020_2024.csv", index=False)
        zt = zone_pair_table(Path(a.drive), wt_all, extra)
        zt.to_csv(out / "zone_pair_table_2020_2024.csv", index=False)
        zp = zone_phase_by_period(zt, wt_all)
        zp.to_csv(out / "zone_phase_by_wetness_2020_2024.csv", index=False)
        dp = dry_pair_tests(zt)
        dp.to_csv(out / "zone_phase_dry_pairs_2020_2024.csv", index=False)
        dry_figure(out, dp)
    print(f"wrote {sorted(p.name for p in out.iterdir())}")


if __name__ == "__main__":
    main()

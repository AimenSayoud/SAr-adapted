"""X-140 — how far the mat at P6 rises back each winter, 2017/18 → 2025/26, each year's spring-to-spring loop, and what
the loops say about the change since spring 2017.

Levels are window means (autumn 15 Oct – 15 Nov, spring 15 May – 15 Jun) of the plain chain (X-099's, hard-target
reference) at the overpasses, on one unbroken chain segment; the laser at the same overpasses is the truth where it
exists on one datum. A correction for the chain's winter under-reading (a constant, pooled over tracks) is kept only if
it beats the plain chain leaving each laser winter out. The P6 water table gives a second winter estimate (secondary).

Pre-registration: hub ``08_deliverables/insar_full_test/preregistration/X-140.yaml`` (checked before the run).
Reads the cube and the hub; writes only to the hub ``08_deliverables/drying_x140/``. No field value in this file.

    make run E=drying_x140
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import p6_validation_x082 as X  # noqa: E402  (the X-082 readers: pairs and laser, unchanged)

from insar_wetlands import drying as D  # noqa: E402
from insar_wetlands import figstyle as fs  # noqa: E402
from insar_wetlands.cube import cube_root  # noqa: E402
from insar_wetlands.paths import hub_root  # noqa: E402
from insar_wetlands.scoring import Protocol, assert_preregistered, chain  # noqa: E402

HUB = hub_root()
CUBE = cube_root()
FT = HUB / "08_deliverables" / "insar_full_test"
OUT = HUB / "08_deliverables" / "drying_x140"
PREREG = FT / "preregistration" / "X-140.yaml"
PROTOCOL = FT / "protocol.yaml"
X099_WINTERS = FT / "reread_corrected" / "baseline_x099" / "winters.csv"
X139_LASER = HUB / "08_deliverables" / "drying_x139" / "laser_multiyear.csv"
WIN = {"autumn": ((10, 15), (11, 15)), "spring": ((5, 15), (6, 15))}
YEARS = range(2017, 2027)
LAT = 52.76
PRIMARY = "hard_targets"


def window(year: int, name: str) -> tuple[pd.Timestamp, pd.Timestamp]:
    (m0, d0), (m1, d1) = WIN[name]
    return (pd.Timestamp(year=year, month=m0, day=d0, tz="UTC"),
            pd.Timestamp(year=year, month=m1, day=d1, hour=23, minute=59, tz="UTC"))


def window_level(c: pd.DataFrame, year: int, name: str) -> dict | None:
    """Mean chain level and mean laser at the window's overpasses on one segment (the segment with most epochs)."""
    a, b = window(year, name)
    w = c[(c.time >= a) & (c.time <= b)]
    if len(w) == 0:
        return None
    seg = w.seg.value_counts().idxmax()
    w = w[w.seg == seg]
    if len(w) < 2:
        return None
    return {"seg": int(seg), "n": len(w), "est": float(w.est.mean()), "truth": float(w.truth.mean())
            if w.truth.notna().sum() >= 2 else np.nan, "datum": w.datum.iloc[0] if w.datum.nunique() == 1 else "mixed"}


def steps(c: pd.DataFrame, cos: float, track: str, ref: str, kind: str) -> list[dict]:
    """Winter rises (autumn y → spring y+1) or summer drops (spring y → autumn y), vertical cm."""
    out = []
    for y in YEARS:
        a_name, a_year, b_name, b_year = (("autumn", y, "spring", y + 1) if kind == "winter" else
                                          ("spring", y, "autumn", y))
        a, b = window_level(c, a_year, a_name), window_level(c, b_year, b_name)
        label = f"{y}/{str(y + 1)[2:]}" if kind == "winter" else str(y)
        r = {"period": label, "year": y, "track": track, "reference": ref, "kind": kind}
        if a is None or b is None:
            r["status"] = "window not read"
        elif a["seg"] != b["seg"]:
            r["status"] = "chain gap"
        else:
            r.update(status="read", n_from=a["n"], n_to=b["n"], radar_cm=(b["est"] - a["est"]) / cos / 10)
            if np.isfinite(a["truth"]) and np.isfinite(b["truth"]) and a["datum"] == b["datum"] != "mixed":
                r["laser_cm"] = (b["truth"] - a["truth"]) / cos / 10
        out.append(r)
    return out


def water_table_steps(prot: Protocol) -> pd.DataFrame:
    s = pd.read_csv(CUBE / "silver" / "site_hourly.csv.gz", index_col=0, low_memory=False,
                    usecols=["time_utc", "wtd_P6_cm", "censored_P6", "wtd_mat_mean_cm"])
    s.index = pd.to_datetime(s.index, format="ISO8601", utc=True)
    cen = pd.to_numeric(s.censored_P6, errors="coerce").fillna(0)
    rows = []
    for y in YEARS:
        r = {"period": f"{y}/{str(y + 1)[2:]}", "year": y}
        vals = {}
        for name, yy in (("autumn", y), ("spring", y + 1)):
            a, b = window(yy, name)
            w = s[(s.index >= a) & (s.index <= b)]
            if w.wtd_P6_cm.notna().sum() < 24 * 10:
                vals = None
                break
            vals[name] = (float(w.wtd_P6_cm.median()), float(cen[w.index].mean()), float(w.wtd_mat_mean_cm.median()))
        if vals:
            r.update(dwtd_P6_cm=vals["spring"][0] - vals["autumn"][0], dwtd_mat_cm=vals["spring"][2] - vals["autumn"][2],
                     censored_share=max(vals["spring"][1], vals["autumn"][1]))
            r["readable"] = r["censored_share"] <= 0.2
        rows.append(r)
    return pd.DataFrame(rows)


def winter_context(prot: Protocol, pairs_flags: pd.DataFrame) -> pd.DataFrame:
    im = pd.concat([pd.read_csv(CUBE / "raw" / "imgw_normals" / "daily_1991_2016.csv", parse_dates=["date"]),
                    pd.read_csv(CUBE / "raw" / "imgw" / "daily.csv", parse_dates=["date"])]).drop_duplicates(["station", "date"])
    im["et0"] = D.hargreaves_et0(im.tmax_c, im.tmin_c, im.date.dt.dayofyear, LAT)
    day = pd.DataFrame({k: D.station_mean(im.pivot(index="date", columns="station", values=v))
                        for k, v in (("rain", "precip_mm"), ("et0", "et0"), ("snow", "snow_depth_cm"))})
    rows = []
    for y in YEARS:
        a, b = pd.Timestamp(year=y, month=11, day=1), pd.Timestamp(year=y + 1, month=4, day=30)
        w = day[(day.index >= a) & (day.index <= b)]
        r = {"period": f"{y}/{str(y + 1)[2:]}", "year": y, "days": len(w), "rain_novapr_mm": w.rain.sum(min_count=1),
             "et0_novapr_mm": w.et0.sum(min_count=1), "snow_days": int((w.snow > 0).sum())}
        r["cwb_novapr_mm"] = r["rain_novapr_mm"] - r["et0_novapr_mm"]
        pf = pairs_flags[(pairs_flags.t1 >= pd.Timestamp(a, tz="UTC")) & (pairs_flags.t2 <= pd.Timestamp(b, tz="UTC"))]
        r["pairs"] = len(pf)
        r["pairs_snow_or_frozen_share"] = float(pf.cold.mean()) if len(pf) else np.nan
        r["median_revisit_days"] = float(pf.revisit_days.median()) if len(pf) else np.nan
        rows.append(r)
    return pd.DataFrame(rows)


def rmse(v) -> float:
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    return float(np.sqrt(np.mean(v ** 2))) if len(v) else np.nan


def maxabs(v) -> float:
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    return float(np.max(np.abs(v))) if len(v) else np.nan


def md(df: pd.DataFrame, fmt: str = "{:.1f}") -> str:
    cols = list(df.columns)
    out = ["| " + " | ".join(str(c) for c in cols) + " |", "|" + "---|" * len(cols)]
    for r in df.itertuples(index=False):
        out.append("| " + " | ".join(fmt.format(v) if isinstance(v, (float, np.floating)) and np.isfinite(v)
                                     else ("" if isinstance(v, float) else str(v)) for v in r) + " |")
    return "\n".join(out)


def main() -> None:
    assert_preregistered(PREREG, [OUT / "winters.csv", OUT / "loops.csv"])
    prot = Protocol.load(PROTOCOL)
    OUT.mkdir(parents=True, exist_ok=True)
    p = X.load_pairs()
    las_los, lev = X.laser()
    x099 = pd.read_csv(X099_WINTERS)
    x099 = x099[(x099.estimate == "plain chain")]

    rows, cos_by = [], {}
    for track in X.TRACKS:
        cos = float((las_los[track] / lev).dropna().median())
        cos_by[track] = cos
        for ref in X.REFS:
            g = p[p.track == track].sort_values("t1")
            c = chain(g.assign(est_mm=g[f"s1_{ref}"]), "est_mm")
            c["time"] = pd.to_datetime(c.time, utc=True)
            c["truth"] = las_los[track].reindex(pd.DatetimeIndex(c.time).round("h")).to_numpy()
            c["datum"] = prot.datum(c.time)
            # stop rule: X-099's Nov → March steps on its own epochs (line of sight, mm)
            for r in x099[(x099.track == track) & (x099.reference == ref)].itertuples():
                a = c[c.time.dt.date == pd.Timestamp(r._2).date()]
                b = c[c.time.dt.date == pd.Timestamp(r.to).date()]
                a, b = a[a.seg.isin(b.seg)], b[b.seg.isin(a.seg)]
                if len(a) != 1 or len(b) != 1 or abs((b.est.iloc[0] - a.est.iloc[0]) - r.insar_change_mm) > 1e-3:
                    raise SystemExit(f"stop: X-099 {r.winter} {track} {ref} not reproduced")
            rows += steps(c, cos, track, ref, "winter") + steps(c, cos, track, ref, "summer")
    t = pd.DataFrame(rows)
    t["error_cm"] = t.radar_cm - t.laser_cm

    # A' — leave-one-winter-out constant correction, pooled over tracks
    scores = []
    for ref in X.REFS:
        w = t[(t.kind == "winter") & (t.reference == ref) & t.laser_cm.notna()].copy()
        lw = sorted(w.period.unique())
        w["corr_lowo_cm"] = np.nan
        for out_w in lw:
            others = w[w.period != out_w]
            cc = float((others.laser_cm - others.radar_cm).mean())
            m = w.period == out_w
            w.loc[m, "corr_lowo_cm"] = w.loc[m, "radar_cm"] + cc
        rmse_a, rmse_ac = rmse(w.radar_cm - w.laser_cm), rmse(w.corr_lowo_cm - w.laser_cm)
        c_all = float((w.laser_cm - w.radar_cm).mean())
        keep = rmse_ac < rmse_a
        scores.append({"reference": ref, "laser_winters": len(lw), "track_winters": len(w), "rmse_A_cm": rmse_a,
                       "rmse_Acorr_lowo_cm": rmse_ac, "max_A_cm": maxabs(w.radar_cm - w.laser_cm),
                       "max_Acorr_lowo_cm": maxabs(w.corr_lowo_cm - w.laser_cm), "c_all_cm": c_all,
                       "kept": "A' (corrected)" if keep else "A (plain chain)", "role": "primary" if ref == PRIMARY else "secondary"})
        t.loc[w.index, "corr_lowo_cm"] = w.corr_lowo_cm
        # per-track c (secondary)
        for track in X.TRACKS:
            wt_ = w[w.track == track]
            errs = []
            for out_w in lw:
                o = wt_[wt_.period != out_w]
                m = wt_[wt_.period == out_w]
                if len(m) and len(o):
                    errs.append(float(m.radar_cm.iloc[0] + (o.laser_cm - o.radar_cm).mean() - m.laser_cm.iloc[0]))
            scores.append({"reference": ref, "laser_winters": len(lw), "track_winters": len(wt_),
                           "rmse_A_cm": rmse(wt_.radar_cm - wt_.laser_cm), "rmse_Acorr_lowo_cm": rmse(errs),
                           "kept": f"per-track c, {track}", "role": "secondary"})
    scores = pd.DataFrame(scores)
    prim = scores[(scores.reference == PRIMARY) & (scores.role == "primary")].iloc[0]
    use_corr = prim.kept.startswith("A'")
    c_all = float(prim.c_all_cm)
    w_err_typ = float(prim.rmse_Acorr_lowo_cm if use_corr else prim.rmse_A_cm)
    w_err_max = float(prim.max_Acorr_lowo_cm if use_corr else prim.max_A_cm)
    sm = t[(t.kind == "summer") & (t.reference == PRIMARY) & t.laser_cm.notna()]
    s_err_typ, s_err_max = rmse(sm.error_cm), maxabs(sm.error_cm)

    # final estimates
    tp = t[t.reference == PRIMARY].copy()
    tp["estimate_cm"] = tp.radar_cm + np.where((tp.kind == "winter") & use_corr, c_all, 0.0)
    tp["err_typ_cm"] = np.where(tp.kind == "winter", w_err_typ, s_err_typ)
    tp["err_max_cm"] = np.where(tp.kind == "winter", w_err_max, s_err_max)
    def wsent(v, e):
        if not np.isfinite(v):
            return ""
        return "no rise distinguishable" if abs(v) <= e else (f"rose by {v:.1f} cm" if v > 0 else f"fell by {abs(v):.1f} cm")
    winters = tp[tp.kind == "winter"].copy()
    winters["sentence"] = [wsent(v, e) for v, e in zip(winters.estimate_cm, winters.err_typ_cm)]

    # B — water table (secondary)
    wtab = water_table_steps(prot)
    lw_mean = (t[(t.kind == "winter") & (t.reference == PRIMARY) & t.laser_cm.notna()]
               .groupby("period").laser_cm.mean().rename("laser_mean_cm"))
    b = wtab.merge(lw_mean, left_on="period", right_index=True, how="left")
    cal = b[b.readable.fillna(False).astype(bool) & b.laser_mean_cm.notna()]
    b["B_cm"], b["B_lowo_err_cm"] = np.nan, np.nan
    if len(cal) >= 1:
        k_all = float((cal.dwtd_P6_cm * cal.laser_mean_cm).sum() / (cal.dwtd_P6_cm ** 2).sum())
        b.loc[b.readable.fillna(False).astype(bool), "B_cm"] = k_all * b.dwtd_P6_cm
        for i, r in cal.iterrows():
            o = cal.drop(i)
            if len(o):
                k = float((o.dwtd_P6_cm * o.laser_mean_cm).sum() / (o.dwtd_P6_cm ** 2).sum())
                b.loc[i, "B_lowo_err_cm"] = k * r.dwtd_P6_cm - r.laser_mean_cm
        b["k_all"] = k_all
    b_rmse = rmse(b.B_lowo_err_cm)

    # loops (spring y → spring y+1) per track and two-track mean
    loops = []
    for y in YEARS:
        r = {"year": f"spring {y} → spring {y + 1}", "y": y}
        vals = []
        for track in X.TRACKS:
            s_ = tp[(tp.kind == "summer") & (tp.year == y) & (tp.track == track)]
            w_ = tp[(tp.kind == "winter") & (tp.year == y) & (tp.track == track)]
            sv = float(s_.estimate_cm.iloc[0]) if len(s_) else np.nan
            wv = float(w_.estimate_cm.iloc[0]) if len(w_) else np.nan
            r[f"{track[:3]}_summer_cm"], r[f"{track[:3]}_winter_cm"] = sv, wv
            r[f"{track[:3]}_loop_cm"] = sv + wv
            vals.append(sv + wv)
            ls = s_.laser_cm.iloc[0] if len(s_) else np.nan
            lw_ = w_.laser_cm.iloc[0] if len(w_) else np.nan
            r[f"{track[:3]}_laser_loop_cm"] = ls + lw_
        r["mean_loop_cm"] = float(np.mean(vals)) if all(np.isfinite(vals)) else np.nan
        r["loop_err_cm"] = float(np.sqrt(s_err_typ ** 2 + w_err_typ ** 2))
        x = r["mean_loop_cm"]
        e = r["loop_err_cm"]
        r["sentence"] = ("" if not np.isfinite(x) else f"ended the year higher by {x:.1f} cm" if x > e else
                         f"ended the year lower by {abs(x):.1f} cm" if x < -e else "back to about the same level")
        loops.append(r)
    loops = pd.DataFrame(loops)
    full = loops[loops.y.between(2017, 2025)]
    if full.mean_loop_cm.notna().all():
        tot, tot_e = float(full.mean_loop_cm.sum()), float(np.sqrt((full.loop_err_cm ** 2).sum()))
        since = (f"Sentinel-1 puts the mat {abs(tot):.1f} ± {tot_e:.1f} cm {'higher' if tot > 0 else 'lower'} in spring 2026 "
                 f"than in spring 2017" if abs(tot) > 2 * tot_e else
                 f"Sentinel-1 cannot tell the change since 2017 from zero ({tot:+.1f} ± {tot_e:.1f} cm)")
    else:
        miss = ", ".join(str(v) for v in full[full.mean_loop_cm.isna()].y)
        tot, tot_e = np.nan, np.nan
        since = f"the loops cannot be summed since 2017: not read for {miss}"
    chk = loops[loops.y.isin([2023, 2024])]
    chk_radar = float(chk.mean_loop_cm.sum()) if chk.mean_loop_cm.notna().all() else np.nan
    lz = pd.read_csv(X139_LASER)
    lsp = lz[(lz.window == "spring")].set_index("year").laser_cm
    chk_laser = float(lsp.get(2025, np.nan) - lsp.get(2023, np.nan))

    ctx_flags = pd.read_csv(CUBE / "gold" / "plots_pair.csv.gz", low_memory=False,
                            usecols=["plot", "grid", "track", "t1", "t2", "revisit_days", "flag_snow_imgw_t1",
                                     "flag_snow_imgw_t2", "flag_frozen_soil_om_t1", "flag_frozen_soil_om_t2"])
    ctx_flags = ctx_flags[(ctx_flags["plot"] == "P6") & (ctx_flags.grid == "40m")].copy()
    for col in ("t1", "t2"):
        ctx_flags[col] = pd.to_datetime(ctx_flags[col], format="ISO8601", utc=True)
    fl = [c for c in ctx_flags.columns if c.startswith("flag_")]
    ctx_flags["cold"] = ctx_flags[fl].astype(str).isin(["True", "1", "1.0"]).any(axis=1)
    ctx_flags = ctx_flags.merge(p[["track", "t1", "t2"]], on=["track", "t1", "t2"])    # consecutive pairs only
    ctx = winter_context(prot, ctx_flags)
    wl = winters.groupby("period").agg(est_mean_cm=("estimate_cm", "mean"), laser_mean_cm=("laser_cm", "mean")).reset_index()
    ctx = ctx.merge(wl, on="period", how="left")

    # outputs
    t.to_csv(OUT / "steps_all.csv", index=False, float_format="%.4f")
    winters.to_csv(OUT / "winters.csv", index=False, float_format="%.4f")
    tp[tp.kind == "summer"].to_csv(OUT / "summers.csv", index=False, float_format="%.4f")
    loops.to_csv(OUT / "loops.csv", index=False, float_format="%.4f")
    scores.to_csv(OUT / "scores.csv", index=False, float_format="%.4f")
    b.to_csv(OUT / "water_table_estimate.csv", index=False, float_format="%.4f")
    ctx.to_csv(OUT / "winter_context.csv", index=False, float_format="%.2f")
    verdict = {"kept": prim.kept, "correction_cm": c_all if use_corr else 0.0, "winter_error_typ_cm": w_err_typ,
               "winter_error_max_cm": w_err_max, "summer_error_typ_cm": s_err_typ, "summer_error_max_cm": s_err_max,
               "since_2017": since, "check_2023_2025_radar_cm": chk_radar, "check_2023_2025_laser_cm": chk_laser,
               "water_table_lowo_rmse_cm": b_rmse}
    (OUT / "verdict.json").write_text(json.dumps(verdict, indent=2, ensure_ascii=False))

    # figure
    with fs.style():
        fig, axes = plt.subplots(2, 1, figsize=fs.size("full", 140), gridspec_kw={"hspace": 0.5, "left": 0.09,
                                 "right": 0.97, "top": 0.93, "bottom": 0.1})
        ax = axes[0]
        ax.axhline(0, color=fs.INK, lw=0.6)
        for track, dx in (("ascending", -0.15), ("descending", 0.15)):
            s = winters[(winters.track == track) & winters.estimate_cm.notna()]
            ax.errorbar(s.year + dx, s.estimate_cm, yerr=s.err_typ_cm, fmt="o", ms=4, lw=0.8, capsize=2,
                        color=fs.TRACK[track], label=f"radar winter rise, {track}")
            ax.plot(s.year + dx, s.laser_cm, ls="none", marker="_", ms=9, mew=1.6, color=fs.CLASS["laser"],
                    label="laser" if track == "ascending" else None)
        bb = b[b.B_cm.notna()]
        ax.plot(bb.year, bb.B_cm, "D", ms=4, mfc="none", color=fs.CLASS["mat"], label="from the water table (secondary)")
        ax.set_xticks(list(YEARS)[:-1], [f"{y}/{str(y + 1)[2:]}" for y in list(YEARS)[:-1]])
        ax.set_ylabel("autumn → spring (cm)")
        ax.set_title("(a) How far the mat rose back each winter at P6", loc="left")
        ax.legend(frameon=False, ncol=2, fontsize=7, loc="upper left")
        ax = axes[1]
        ax.axhline(0, color=fs.INK, lw=0.6)
        lp = loops[loops.mean_loop_cm.notna()]
        ax.errorbar(lp.y, lp.mean_loop_cm, yerr=lp.loop_err_cm, fmt="s", ms=4, capsize=2, color=fs.CLASS["mat"],
                    label="radar, two-track mean (± typical)")
        ll = loops[["y", "asc_laser_loop_cm", "des_laser_loop_cm"]].set_index("y").mean(axis=1).dropna()
        ax.plot(ll.index, ll.values, ls="none", marker="_", ms=9, mew=1.6, color=fs.CLASS["laser"], label="laser")
        ax.set_xticks(list(YEARS)[:-1], [f"{y}→{str(y + 1)[2:]}" for y in list(YEARS)[:-1]])
        ax.set_ylabel("spring → next spring (cm)")
        ax.set_title("(b) The full year: did the winter make up the summer's drop?", loc="left")
        ax.legend(frameon=False, fontsize=7, loc="upper left")
        fs.export(fig, OUT / "fig_x140", provenance={"script": "experiments/drying_x140.py", "ticket": "X-140"})
        plt.close(fig)

    wshow = winters[["period", "track", "status", "n_from", "n_to", "radar_cm", "estimate_cm", "err_typ_cm",
                     "laser_cm", "sentence"]]
    sshow = tp[tp.kind == "summer"][["period", "track", "status", "radar_cm", "laser_cm", "error_cm"]]
    lshow = loops[["year", "asc_summer_cm", "asc_winter_cm", "asc_loop_cm", "des_summer_cm", "des_winter_cm",
                   "des_loop_cm", "mean_loop_cm", "loop_err_cm", "asc_laser_loop_cm", "des_laser_loop_cm", "sentence"]]
    readme = f"""# X-140 — how far the mat rises back each winter, and each year's full loop (generated)

Generated by `experiments/drying_x140.py` (`make run E=drying_x140`); do not edit. Pre-registration:
`../insar_full_test/preregistration/X-140.yaml`. Stop rule passed (X-099's winter steps reproduced on its epochs).
Windows: autumn 15 Oct – 15 Nov, spring 15 May – 15 Jun; a window's level = mean at its Sentinel-1 overpasses; P6,
hard-target reference, vertical cm, positive = up.

## Answers (pre-registered rules)

- **Winter estimator kept:** {prim.kept} — leave-one-winter-out error {prim.rmse_Acorr_lowo_cm:.1f} cm against
  {prim.rmse_A_cm:.1f} cm for the plain chain{f'; correction +{c_all:.1f} cm per winter' if use_corr else ''}.
  Winter error used for winters without a laser: typical ± {w_err_typ:.1f} cm, worst seen ± {w_err_max:.1f} cm. Summer
  error: typical ± {s_err_typ:.1f} cm, worst ± {s_err_max:.1f} cm.
- **Since spring 2017:** {since}.
- **Check against the laser:** the radar loops spring 2023 → spring 2025 sum to {chk_radar:+.1f} cm; the laser says
  {chk_laser:+.1f} cm.
- **Water table as a second estimate (secondary):** leave-one-winter-out error {b_rmse:.1f} cm over the laser winters
  it can be read on.

![summary](fig_x140.png)

## Each winter (autumn → next spring)

{md(wshow)}

## Each summer (spring → autumn) and its error against the laser

{md(sshow)}

## Each year's full loop (spring → next spring)

Loop error = √(summer error² + winter error²), typical; the two-track mean keeps one track's error (the tracks share
errors, X-136).

{md(lshow)}

## Scores (laser winters, leave one winter out)

{md(scores, "{:.2f}")}

## Water table as a second estimate (secondary)

ΔWTD = change of P6 water-table depth between the windows' medians (cm, negative below the floating surface); the mat's
rise = k · ΔWTD, k fitted on the laser winters (left out in turn for the error).

{md(b, "{:.2f}")}

## Winter context (IMGW November → April; pairs with snow or frozen ground)

{md(ctx, "{:.1f}")}
"""
    (OUT / "README.md").write_text(readme)
    print(readme)


if __name__ == "__main__":
    main()

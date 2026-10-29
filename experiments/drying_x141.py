"""X-141 — the mat's winter rise at P6 from the national stations' winter weather, scored on the four laser winters;
every winter since 1991; each year's loop; and the big year table with the winter columns.

Model W: winter rise (autumn window → spring window, X-140) = a + b · IMGW November–April water balance (rain −
Hargreaves ET0, Piła + Poznań). Kept only if, leaving each laser winter out, it beats both "an average winter" and the
radar chain. Other single predictors and a nested choice among them are secondary. Routines in
``insar_wetlands.drying`` (tested).

Pre-registration: hub ``08_deliverables/insar_full_test/preregistration/X-141.yaml`` (checked before the run).
Reads the cube and the hub; writes only to the hub ``08_deliverables/drying_x141/``. No field value in this file.

    make run E=drying_x141
"""
from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from insar_wetlands import drying as D  # noqa: E402
from insar_wetlands import figstyle as fs  # noqa: E402
from insar_wetlands.cube import cube_root  # noqa: E402
from insar_wetlands.paths import hub_root  # noqa: E402
from insar_wetlands.scoring import assert_preregistered  # noqa: E402

HUB = hub_root()
CUBE = cube_root()
FT = HUB / "08_deliverables" / "insar_full_test"
OUT = HUB / "08_deliverables" / "drying_x141"
D138, D139, D140 = (HUB / "08_deliverables" / f"drying_x{k}" for k in (138, 139, 140))
PREREG = FT / "preregistration" / "X-141.yaml"
LAT = 52.76
WINTERS = range(1991, 2026)                       # winter y/y+1, November y → April y+1
NORMAL_W = (1991, 2019)                           # winters 1991/92 → 2019/20


def label(y: int) -> str:
    return f"{y}/{str(y + 1)[2:]}"


def imgw_days() -> pd.DataFrame:
    im = pd.concat([pd.read_csv(CUBE / "raw" / "imgw_normals" / "daily_1991_2016.csv", parse_dates=["date"]),
                    pd.read_csv(CUBE / "raw" / "imgw" / "daily.csv", parse_dates=["date"])]).drop_duplicates(["station", "date"])
    im["et0"] = D.hargreaves_et0(im.tmax_c, im.tmin_c, im.date.dt.dayofyear, LAT)
    return pd.DataFrame({k: D.station_mean(im.pivot(index="date", columns="station", values=v))
                         for k, v in (("rain", "precip_mm"), ("et0", "et0"), ("snow", "snow_depth_cm"))})


def winter_predictors(day: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for y in WINTERS:
        w = day[(day.index >= f"{y}-11-01") & (day.index <= f"{y + 1}-04-30")]
        a = day[(day.index >= f"{y}-08-01") & (day.index <= f"{y}-10-31")]
        ok = w.rain.notna().sum() >= 175
        rows.append({"winter": label(y), "y": y, "complete": ok,
                     "rain_novapr_mm": w.rain.sum() if ok else np.nan,
                     "cwb_novapr_mm": (w.rain.sum() - w.et0.sum()) if ok else np.nan,
                     "snow_days": int((w.snow > 0).sum()) if ok else np.nan,
                     "snow_max_cm": float(w.snow.max()) if ok else np.nan,
                     "rain_augoct_mm": a.rain.sum() if a.rain.notna().sum() >= 85 else np.nan})
    t = pd.DataFrame(rows)
    n = t[t.y.between(*NORMAL_W)]
    t["rain_novapr_pct"] = 100 * t.rain_novapr_mm / n.rain_novapr_mm.mean()
    t["cwb_novapr_anom_mm"] = t.cwb_novapr_mm - n.cwb_novapr_mm.mean()
    return t


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
    assert_preregistered(PREREG, [OUT / "big_table.csv", OUT / "models.csv"])
    OUT.mkdir(parents=True, exist_ok=True)
    day = imgw_days()
    wp = winter_predictors(day)

    # X-140: laser and radar winters (two-track means), summers
    w140 = pd.read_csv(D140 / "winters.csv")
    wm = w140.groupby("year").agg(laser_cm=("laser_cm", "mean"), radar_cm=("radar_cm", "mean"),
                                  n_laser=("laser_cm", "count")).reset_index().rename(columns={"year": "y"})
    s140 = pd.read_csv(D140 / "summers.csv")
    sm = s140.groupby("year").agg(summer_laser_cm=("laser_cm", "mean"), summer_radar_cm=("estimate_cm", "mean"),
                                  n_slaser=("laser_cm", "count")).reset_index().rename(columns={"year": "y"})
    sm.loc[sm.n_slaser < 2, "summer_laser_cm"] = np.nan
    v140 = json.loads((D140 / "verdict.json").read_text())
    s_err = float(v140["summer_error_typ_cm"])
    wp = wp.merge(wm, on="y", how="left").merge(sm[["y", "summer_laser_cm", "summer_radar_cm"]], on="y", how="left")
    wp.loc[wp.n_laser < 2, "laser_cm"] = np.nan
    wp["prev_summer_cm"] = wp.summer_laser_cm.fillna(wp.summer_radar_cm)

    # fit and score on the laser winters
    cal = wp[wp.laser_cm.notna()].sort_values("y").reset_index(drop=True)
    y = cal.laser_cm.to_numpy()
    preds = {"climatology (baseline)": D.lowo_climatology(y), "radar chain (baseline)": cal.radar_cm.to_numpy()}
    cands = {"W: CWB Nov–Apr": "cwb_novapr_mm", "rain Nov–Apr": "rain_novapr_mm", "snow days": "snow_days",
             "max snow depth": "snow_max_cm", "rain Aug–Oct": "rain_augoct_mm", "preceding summer drop": "prev_summer_cm"}
    for name, col in cands.items():
        preds[name] = D.lowo_linear(cal[col], y)
    nest, chosen = D.nested_selection({k: cal[c].to_numpy() for k, c in cands.items()}, y)
    preds["nested choice among the six"] = nest
    models = []
    for name, pr in preds.items():
        e = pr - y
        role = ("baseline" if "baseline" in name else "primary" if name.startswith("W:") else "secondary")
        r = {"model": name, "role": role, "rmse_lowo_cm": rmse(e), "max_lowo_cm": maxabs(e)}
        if name in cands:
            b, a = np.polyfit(cal[cands[name]], y, 1)
            r.update(a=a, b=b, r_in_sample=float(np.corrcoef(cal[cands[name]], y)[0, 1]))
        for wl, pv, ev in zip(cal.winter, pr, e):
            r[f"pred_{wl}"] = pv
        models.append(r)
    models = pd.DataFrame(models)
    W = models[models.model.str.startswith("W:")].iloc[0]
    base = models[models.role == "baseline"]
    keep = bool(W.rmse_lowo_cm < base.rmse_lowo_cm.min())
    w_typ, w_max = float(W.rmse_lowo_cm), float(W.max_lowo_cm)
    lo, hi = cal.cwb_novapr_mm.min(), cal.cwb_novapr_mm.max()
    wp["model_cm"] = W.a + W.b * wp.cwb_novapr_mm
    wp["extrapolated"] = (wp.cwb_novapr_mm < lo) | (wp.cwb_novapr_mm > hi)
    W_sent = (f"the winter model W is kept: leaving each laser winter out it misses by {w_typ:.1f} cm (typical), against "
              f"{base.set_index('model').rmse_lowo_cm.to_dict()}" if keep else
              f"the winter model W is not kept: leaving each laser winter out it misses by {w_typ:.1f} cm, not better than "
              f"the baselines ({', '.join(f'{k} {v:.1f} cm' for k, v in base.set_index('model').rmse_lowo_cm.items())})")

    # summer counterpart (printed only)
    smr = []
    for yy in range(2017, 2026):
        d = day[(day.index >= f"{yy}-05-01") & (day.index <= f"{yy}-10-31")]
        smr.append({"y": yy, "cwb_mayoct_mm": d.rain.sum() - d.et0.sum()})
    smr = pd.DataFrame(smr).merge(sm[["y", "summer_laser_cm", "summer_radar_cm"]], on="y", how="left")
    sc = smr[smr.summer_laser_cm.notna()]
    s_pred = D.lowo_linear(sc.cwb_mayoct_mm, sc.summer_laser_cm)
    summer_model = {"n": len(sc), "rmse_lowo_cm": rmse(s_pred - sc.summer_laser_cm),
                    "radar_rmse_cm": rmse(sc.summer_radar_cm - sc.summer_laser_cm),
                    "climatology_rmse_cm": rmse(D.lowo_climatology(sc.summer_laser_cm) - sc.summer_laser_cm)}

    # loops and the gate
    loops = []
    for yy in range(2017, 2026):
        r = wp[wp.y == yy].iloc[0]
        wv = float(r.model_cm) if keep else float(r.radar_cm)
        we = w_typ if keep else float(v140["winter_error_typ_cm"])
        sv = float(sm.set_index("y").summer_radar_cm.get(yy, np.nan))
        lv = float(r.summer_laser_cm + r.laser_cm) if np.isfinite(r.summer_laser_cm) and np.isfinite(r.laser_cm) else np.nan
        loops.append({"y": yy, "year": f"spring {yy} → spring {yy + 1}", "summer_cm": sv, "winter_cm": wv,
                      "winter_source": "model W" + (" (extrapolated)" if keep and r.extrapolated else "") if keep else "radar chain",
                      "loop_cm": sv + wv, "loop_err_cm": float(np.sqrt(s_err ** 2 + we ** 2)), "laser_loop_cm": lv})
    loops = pd.DataFrame(loops)
    l140 = pd.read_csv(D140 / "loops.csv")
    l140["laser_mean"] = l140[["asc_laser_loop_cm", "des_laser_loop_cm"]].mean(axis=1)
    g = loops[loops.y.isin([2023, 2024])]
    g_est, g_err = float(g.loop_cm.sum()), float(np.sqrt((g.loop_err_cm ** 2).sum()))
    g_laser = float(l140[l140.y.isin([2023, 2024])].laser_mean.sum())
    gate = abs(g_est - g_laser) <= 2 * g_err
    full = loops[loops.y.between(2017, 2025)]
    tot, tot_e = float(full.loop_cm.sum()), float(np.sqrt((full.loop_err_cm ** 2).sum()))
    if not gate:
        since = "the change since 2017 cannot be given from these estimates (the 2023 → 2025 gate fails)"
    elif abs(tot) <= 2 * tot_e:
        since = f"the change since spring 2017 is not distinguishable from zero ({tot:+.1f} ± {tot_e:.1f} cm)"
    else:
        since = f"the mat stood {abs(tot):.1f} ± {tot_e:.1f} cm {'higher' if tot > 0 else 'lower'} in spring 2026 than in spring 2017"
    gate_txt = (f"gate {'passes' if gate else 'fails'}: estimated loops 2023 → 2025 {g_est:+.1f} ± {g_err:.1f} cm, "
                f"laser {g_laser:+.1f} cm")

    # the big table
    rs = pd.read_csv(D138 / "radar_seasons.csv")
    rs = rs[rs.reference == "hard_targets"]
    ym = pd.read_csv(D138 / "radar_year_means.csv").set_index("year")
    yt9 = pd.read_csv(D139 / "year_table.csv").set_index("year")
    rows = []
    for yy in range(2017, 2027):
        r = {"year": yy}
        for tr in ("ascending", "descending"):
            q = rs[(rs.year == yy) & (rs.track == tr)]
            if len(q):
                q = q.iloc[0]
                r[f"summer_{tr[:3]}"] = f"{pd.Timestamp(q.start):%d %b}→{pd.Timestamp(q.end):%d %b}: {q.net_cm:+.1f}"
                r["revisit"] = q.revisit.split()[0]
                r["err"] = q.net_typical_cm
            lasers = rs[(rs.year == yy)].laser_net_cm.dropna()
            r["summer_laser"] = (" / ".join(f"{v:+.1f}" for v in lasers) if len(lasers) else "")
        r["summer_mean_cm"] = ym.net_mean_cm.get(yy, np.nan)
        r["rain_year"] = (f"{yt9.rain_mm.get(yy):.0f} mm ({yt9.rain_pct.get(yy):.0f} %, {yt9.rain_class.get(yy)})"
                          if yy in yt9.index and np.isfinite(yt9.rain_mm.get(yy)) else "")
        r["wt_P6_cm"] = yt9.wt_P6_mean_cm.get(yy, np.nan) if yy in yt9.index else np.nan
        r["open_water_ha"] = yt9.open_water_ha.get(yy, "") if yy in yt9.index else ""
        w = wp[wp.y == yy]
        if len(w) and yy <= 2025:
            w = w.iloc[0]
            r["winter"] = w.winter
            r["winter_rain"] = f"{w.rain_novapr_mm:.0f} mm ({w.rain_novapr_pct:.0f} %)"
            r["winter_cwb_mm"] = w.cwb_novapr_mm
            r["winter_laser_cm"] = w.laser_cm
            r["winter_radar_cm"] = w.radar_cm
            r["winter_model_cm"] = (f"{w.model_cm:+.1f} (extrapolated)" if w.extrapolated else
                                    f"{w.model_cm:+.1f} ± {w_typ:.1f}")
            lp = loops[loops.y == yy].iloc[0]
            r["year_loop"] = (f"laser {lp.laser_loop_cm:+.1f}" if np.isfinite(lp.laser_loop_cm) else
                              f"est. {lp.loop_cm:+.1f} ± {lp.loop_err_cm:.1f}")
        rows.append(r)
    big = pd.DataFrame(rows)

    # outputs
    wp.to_csv(OUT / "winters_1991_2026.csv", index=False, float_format="%.3f")
    models.to_csv(OUT / "models.csv", index=False, float_format="%.4f")
    loops.to_csv(OUT / "loops.csv", index=False, float_format="%.3f")
    big.to_csv(OUT / "big_table.csv", index=False, float_format="%.2f")
    pd.DataFrame({"fold_left_out": list(cal.winter), "chosen": chosen}).to_csv(OUT / "nested_choices.csv", index=False)
    (OUT / "verdict.json").write_text(json.dumps({"model_W": W_sent, "kept": keep, "winter_error_typ_cm": w_typ,
                                                 "winter_error_max_cm": w_max, "gate": gate_txt, "since_2017": since,
                                                 "summer_counterpart": summer_model}, indent=2, ensure_ascii=False))

    with fs.style():
        fig, axes = plt.subplots(2, 1, figsize=fs.size("full", 150), gridspec_kw={"hspace": 0.5, "left": 0.09,
                                 "right": 0.97, "top": 0.93, "bottom": 0.1})
        ax = axes[0]
        xx = np.linspace(min(lo, wp.cwb_novapr_mm.min()), max(hi, wp.cwb_novapr_mm.max()), 50)
        ax.plot(xx, W.a + W.b * xx, color=fs.MUTED, lw=0.8, ls="--", label="model W (fit on the four laser winters)")
        ax.scatter(cal.cwb_novapr_mm, cal.laser_cm, color=fs.CLASS["laser"], zorder=3, label="laser winter rise")
        cal_off = {"2024/25": (-5, 4, "right")}
        for r in cal.itertuples():
            dx, dy, ha = cal_off.get(r.winter, (4, 3, "left"))
            ax.annotate(r.winter, (r.cwb_novapr_mm, r.laser_cm), xytext=(dx, dy), textcoords="offset points",
                        ha=ha, fontsize=7)
        o = wp[wp.y.between(2017, 2025) & wp.laser_cm.isna()]
        ax.scatter(o.cwb_novapr_mm, o.model_cm, marker="o", facecolors="none", edgecolors=fs.CLASS["mat"],
                   label="model, winters without laser")
        o_off = {"2017/18": (-5, 5, "right"), "2020/21": (5, -10, "left"), "2019/20": (-5, 5, "right"),
                 "2018/19": (5, -9, "left"), "2025/26": (0, 6, "center")}
        for r in o.itertuples():
            dx, dy, ha = o_off.get(r.winter, (4, -9, "left"))
            ax.annotate(r.winter, (r.cwb_novapr_mm, r.model_cm), xytext=(dx, dy), textcoords="offset points",
                        ha=ha, fontsize=7, color=fs.MUTED)
        ax.axvspan(lo, hi, color=fs.GRID, alpha=0.3, lw=0)
        ax.set_ylim(-6.8, 15.2)
        ax.set_xlabel("IMGW November–April water balance (mm)")
        ax.set_ylabel("winter rise (cm)")
        ax.set_title("(a) Winter rise against winter water balance (grey: the laser winters' range)", loc="left")
        ax.legend(frameon=False, fontsize=7, loc="upper left")
        ax = axes[1]
        ww = wp[wp.complete]
        col = [fs.MUTED if e else fs.CLASS["mat"] for e in ww.extrapolated]
        ax.bar(ww.y, ww.model_cm, color=col, width=0.7)
        ax.plot(cal.y, cal.laser_cm, ls="none", marker="_", ms=10, mew=1.8, color=fs.CLASS["laser"], label="laser")
        ax.axhline(0, color=fs.INK, lw=0.6)
        ax.set_ylabel("modelled winter rise (cm)")
        ax.set_title("(b) Model W for every winter since 1991/92 (grey bars: outside the laser winters' range)", loc="left")
        ax.legend(frameon=False, fontsize=7, loc="upper left")
        fs.export(fig, OUT / "fig_x141", provenance={"script": "experiments/drying_x141.py", "ticket": "X-141"})
        fs.export(
            fig,
            HUB / "08_deliverables" / "article_plan" / "figures" / "final" / "S17_drying_trajectory",
            provenance={"script": "experiments/drying_x141.py", "ticket": "X-141", "role": "S17"},
        )
        plt.close(fig)

    mshow = models[["model", "role", "rmse_lowo_cm", "max_lowo_cm", "a", "b", "r_in_sample"]]
    wshow = wp[wp.y >= 2016][["winter", "rain_novapr_mm", "rain_novapr_pct", "cwb_novapr_mm", "snow_days", "snow_max_cm",
                              "laser_cm", "radar_cm", "model_cm", "extrapolated"]]
    readme = f"""# X-141 — the winter model and the big table (generated)

Generated by `experiments/drying_x141.py` (`make run E=drying_x141`); do not edit. Pre-registration:
`../insar_full_test/preregistration/X-141.yaml`. Builds on X-138, X-139, X-140.

## Answers (pre-registered rules)

- **Winter model:** {W_sent}.
- **Loops:** {gate_txt}.
- **Since spring 2017:** {since}.

![summary](fig_x141.png)

## The big table

Summer = each track's first → last overpass in May–November and its net (cm; X-138), the two-track mean with the
error of its revisit, the laser's season net. Winter = the winter that follows (November → April rain from IMGW,
% of the 1991/92–2019/20 normal; water balance = rain − ET0), the mat's rise from the autumn window (15 Oct – 15 Nov) to
the spring window (15 May – 15 Jun) by the laser, by the radar (two-track mean, error about ± 5 cm, X-140) and by the
model. Year loop = spring → next spring (laser where it exists, else the estimate).

{md(big)}

## Models, scored leaving each laser winter out (cm)

{md(mshow, "{:.3g}")}

Nested choice per fold: {', '.join(f'{a}: {b}' for a, b in zip(cal.winter, chosen))}.

Summer counterpart (printed only): {json.dumps({k: round(v, 2) if isinstance(v, float) else v for k, v in summer_model.items()})}.

## Winters 2016/17 → 2025/26

{md(wshow)}

All winters since 1991/92: `winters_1991_2026.csv`.

## Loops

{md(loops.drop(columns='y'))}
"""
    (OUT / "README.md").write_text(readme)
    print(readme)


if __name__ == "__main__":
    main()

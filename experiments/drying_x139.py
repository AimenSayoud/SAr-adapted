"""X-139 — was each year dry or wet by the national stations, did the drier seasons sink further, and how far has the
mat moved across years where a measurement across years exists? One table per year with everything we hold.

Rain and reference evapotranspiration (Hargreaves–Samani from daily Tmax/Tmin) from IMGW synop Piła (230) and Poznań
(330), two-station mean, against their 1991–2020 normal (Kaczorowska's classes); the season sinking of X-138 against
the season's water balance; the laser at P6 on one datum (B) in fixed spring and autumn windows. Routines in
``insar_wetlands.drying`` (tested).

Pre-registration: hub ``08_deliverables/insar_full_test/preregistration/X-139.yaml`` (checked before the run).
Reads the cube and the hub; writes only to the hub ``08_deliverables/drying_x139/``. No field value in this file.

    make run E=drying_x139
"""
from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

from insar_wetlands import drying as D  # noqa: E402
from insar_wetlands import figstyle as fs  # noqa: E402
from insar_wetlands.cube import cube_root  # noqa: E402
from insar_wetlands.paths import hub_root  # noqa: E402
from insar_wetlands.scoring import Protocol, assert_preregistered  # noqa: E402

HUB = hub_root()
CUBE = cube_root()
FT = HUB / "08_deliverables" / "insar_full_test"
OUT = HUB / "08_deliverables" / "drying_x139"
X138 = HUB / "08_deliverables" / "drying_x138"
PREREG = FT / "preregistration" / "X-139.yaml"
PROTOCOL = FT / "protocol.yaml"
X137 = FT / "season_total_x137" / "season_totals.csv"
X099_WINTERS = FT / "reread_corrected" / "baseline_x099" / "winters.csv"
LAT = 52.76
SEASON = (5, 11)                                    # May–November, the radar season window
NORMAL = (1991, 2020)
WINDOWS = {"spring": ((5, 15), (6, 15)), "autumn": ((10, 15), (11, 15))}


# --------------------------------------------------------------------------------------------- IMGW
def imgw_daily() -> pd.DataFrame:
    cols = ["station", "date", "tmax_c", "tmin_c", "precip_mm"]
    a = pd.read_csv(CUBE / "raw" / "imgw_normals" / "daily_1991_2016.csv", usecols=cols, parse_dates=["date"])
    b = pd.read_csv(CUBE / "raw" / "imgw" / "daily.csv", usecols=cols, parse_dates=["date"])
    d = pd.concat([a, b]).drop_duplicates(["station", "date"])
    d["et0_mm"] = D.hargreaves_et0(d.tmax_c, d.tmin_c, d.date.dt.dayofyear, LAT)
    out = pd.DataFrame({"precip_mm": D.station_mean(d.pivot(index="date", columns="station", values="precip_mm")),
                        "et0_mm": D.station_mean(d.pivot(index="date", columns="station", values="et0_mm"))})
    full = pd.date_range(out.index.min(), out.index.max(), freq="D")
    return out.reindex(full).rename_axis("date")


def yearly(day: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for y, g in day.groupby(day.index.year):
        s = g[g.index.month.between(*SEASON)]
        ndays = 366 if pd.Timestamp(year=y, month=12, day=31).dayofyear == 366 else 365
        nseason = (pd.Timestamp(year=y, month=SEASON[1] + 1, day=1) - pd.Timestamp(year=y, month=SEASON[0], day=1)).days
        rows.append({"year": y, "days": len(g), "missing_days": int(ndays - g.precip_mm.notna().sum()),
                     "missing_days_season": int(nseason - s.precip_mm.notna().sum()),
                     "rain_mm": g.precip_mm.sum(), "et0_mm": g.et0_mm.sum(),
                     "rain_season_mm": s.precip_mm.sum(), "et0_season_mm": s.et0_mm.sum()})
    t = pd.DataFrame(rows)
    t["cwb_mm"] = t.rain_mm - t.et0_mm
    t["cwb_season_mm"] = t.rain_season_mm - t.et0_season_mm
    t["incomplete"] = t.missing_days > 10
    t["incomplete_season"] = t.missing_days_season > 5
    return t


# --------------------------------------------------------------------------------------------- laser across years
def laser_multiyear(prot: Protocol) -> tuple[pd.DataFrame, pd.Series]:
    s = pd.read_csv(CUBE / "silver" / "site_hourly.csv.gz", index_col=0, low_memory=False,
                    usecols=["time_utc", "laser_surface_cm", "laser_ok", "p6_level_on_laser_datum_cm"])
    s.index = pd.to_datetime(s.index, format="ISO8601", utc=True)
    ok = s.laser_ok.map({True: True, False: False, "True": True, "False": False}).fillna(False).to_numpy()
    lev = s.laser_surface_cm.where(ok)
    wat = s.p6_level_on_laser_datum_cm.where(ok)
    day = pd.DataFrame({"laser": lev.resample("D").mean(), "water": wat.resample("D").mean()})
    day = day[np.asarray(prot.datum(pd.Series(day.index))) == "B"]
    rows = []
    for name, ((m0, d0), (m1, d1)) in WINDOWS.items():
        for y in sorted(set(day.index.year)):
            a, b = pd.Timestamp(year=y, month=m0, day=d0, tz="UTC"), pd.Timestamp(year=y, month=m1, day=d1, tz="UTC")
            w = day[(day.index >= a) & (day.index <= b)]
            n = int(w.laser.notna().sum())
            rows.append({"window": name, "year": y, "valid_days": n,
                         "laser_cm": float(w.laser.median()) if n >= 10 else np.nan,
                         "water_cm": float(w.water.median()) if int(w.water.notna().sum()) >= 10 else np.nan})
    t = pd.DataFrame(rows)
    for name in WINDOWS:
        m = t.window == name
        first = t[m & t.laser_cm.notna()].sort_values("year")
        if len(first):
            t.loc[m, "laser_rel_cm"] = t.loc[m, "laser_cm"] - first.laser_cm.iloc[0]
        firstw = t[m & t.water_cm.notna()].sort_values("year")
        if len(firstw):
            t.loc[m, "water_rel_cm"] = t.loc[m, "water_cm"] - firstw.water_cm.iloc[0]
    return t, day.laser


# --------------------------------------------------------------------------------------------- figure
def figure(yr, normal, ym, lm, laser_day) -> None:
    with fs.style():
        fig, axes = plt.subplots(3, 1, figsize=fs.size("full", 190), gridspec_kw={"hspace": 0.6, "left": 0.09,
                                 "right": 0.97, "top": 0.95, "bottom": 0.07})
        ax = axes[0]
        r = yr[yr.year.between(1991, 2025) & ~yr.incomplete]
        pct = 100 * r.rain_mm / normal["rain_mm"]
        col = [fs.CLASS["grassland"] if v < 90 else (fs.CLASS["open water"] if v > 110 else fs.GRID) for v in pct]
        ax.bar(r.year, pct, color=col, width=0.75)
        ax.axhline(100, color=fs.INK, lw=0.6)
        ax.axvspan(2016.5, 2025.5, color=fs.GRID, alpha=0.3, lw=0)
        ax.set_ylabel("rain, % of 1991–2020")
        ax.set_title("(a) Yearly rain at IMGW Piła + Poznań (orange: dry < 90 %, blue: wet > 110 %; shaded: the radar years)",
                     loc="left")
        ax = axes[1]
        m = ym.merge(yr[["year", "cwb_season_mm"]], on="year")
        m["anom"] = m.cwb_season_mm - normal["cwb_season_mm"]
        ax.scatter(m.anom, m.net_mean_cm, s=18, color=fs.CLASS["mat"], zorder=3)
        for row in m.itertuples():
            ax.annotate(str(row.year), (row.anom, row.net_mean_cm), xytext=(3, 3), textcoords="offset points", fontsize=7)
        ax.axvline(0, color=fs.MUTED, lw=0.5)
        ax.set_xlabel("May–November water balance, anomaly from 1991–2020 (mm)")
        ax.set_ylabel("mat May → Nov (cm)")
        ax.set_title("(b) How far each season sank (Sentinel-1, two-track mean) against its water balance", loc="left")
        ax = axes[2]
        ax.plot(laser_day.index, laser_day - laser_day.dropna().iloc[0], lw=0.6, color=fs.CLASS["laser"])
        for row in lm[lm.laser_cm.notna()].itertuples():
            (m0, d0), (m1, d1) = WINDOWS[row.window]
            a = pd.Timestamp(year=row.year, month=m0, day=d0, tz="UTC")
            b = pd.Timestamp(year=row.year, month=m1, day=d1, tz="UTC")
            ax.plot([a, b], [row.laser_cm - laser_day.dropna().iloc[0]] * 2, color=fs.INK, lw=2)
        ax.set_ylabel("laser level, datum B (cm)")
        ax.set_title("(c) The laser at P6 across years on one datum (black: spring and autumn window medians)", loc="left")
        fs.export(fig, OUT / "fig_x139", provenance={"script": "experiments/drying_x139.py", "ticket": "X-139"})
        plt.close(fig)


def md(df: pd.DataFrame, fmt: str = "{:.1f}") -> str:
    cols = list(df.columns)
    out = ["| " + " | ".join(str(c) for c in cols) + " |", "|" + "---|" * len(cols)]
    for r in df.itertuples(index=False):
        out.append("| " + " | ".join(fmt.format(v) if isinstance(v, (float, np.floating)) and np.isfinite(v)
                                     else ("" if isinstance(v, float) else str(v)) for v in r) + " |")
    return "\n".join(out)


def main() -> None:
    assert_preregistered(PREREG, [OUT / "year_table.csv", OUT / "imgw_yearly.csv"])
    prot = Protocol.load(PROTOCOL)
    OUT.mkdir(parents=True, exist_ok=True)

    # (1) IMGW years against the normal
    day = imgw_daily()
    yr = yearly(day)
    base = yr[yr.year.between(*NORMAL)]
    normal = {k: float(base[k].mean()) for k in ("rain_mm", "et0_mm", "cwb_mm", "rain_season_mm", "et0_season_mm",
                                                  "cwb_season_mm")}
    normal["years"] = int(len(base))
    normal["incomplete_years_in_normal"] = int(base.incomplete.sum())
    yr["rain_pct"] = 100 * yr.rain_mm / normal["rain_mm"]
    yr["rain_season_pct"] = 100 * yr.rain_season_mm / normal["rain_season_mm"]
    yr["rain_class"] = [D.rain_class(p) if not inc else "incomplete" for p, inc in zip(yr.rain_pct, yr.incomplete)]
    yr["rain_season_class"] = [D.rain_class(p) if not inc else "incomplete"
                               for p, inc in zip(yr.rain_season_pct, yr.incomplete_season)]
    yr["cwb_anom_mm"] = yr.cwb_mm - normal["cwb_mm"]
    yr["cwb_season_anom_mm"] = yr.cwb_season_mm - normal["cwb_season_mm"]
    last = day.precip_mm.dropna().index.max()
    win = day[(day.index.month < last.month) | ((day.index.month == last.month) & (day.index.day <= last.day))]
    todate = win.groupby(win.index.year).agg(rain_mm=("precip_mm", "sum"), et0_mm=("et0_mm", "sum"))
    todate["cwb_mm"] = todate.rain_mm - todate.et0_mm
    n_todate = todate.loc[NORMAL[0]:NORMAL[1]].mean()
    td2026 = {"window": f"01-01 → {last:%m-%d}", "rain_mm": float(todate.loc[2026, "rain_mm"]),
              "rain_pct": float(100 * todate.loc[2026, "rain_mm"] / n_todate.rain_mm),
              "cwb_mm": float(todate.loc[2026, "cwb_mm"]), "cwb_anom_mm": float(todate.loc[2026, "cwb_mm"] - n_todate.cwb_mm),
              "rank_driest_cwb_since_1991": int((todate.cwb_mm <= todate.loc[2026, "cwb_mm"]).sum())}
    yr.loc[yr.year == 2026, "rain_class"] = f"to {last:%d %b}: " + D.rain_class(td2026["rain_pct"])

    trends = []
    c = yr[yr.year.between(1991, 2025) & ~yr.incomplete]
    tr = D.trend(c.year, c.cwb_mm)
    long_sent = D.trend_sentence(tr, 0.0, "the water balance at the stations became drier since 1991",
                                 "the water balance at the stations became wetter since 1991",
                                 "no trend in the stations' water balance since 1991")
    trends.append({"series": "IMGW CWB, yearly", "years": "1991–2025", **tr, "sentence": long_sent, "role": "primary"})
    for name, col, yrs in (("IMGW rain, yearly", "rain_mm", (1991, 2025)), ("IMGW CWB, yearly", "cwb_mm", (2017, 2025)),
                           ("IMGW ET0, yearly", "et0_mm", (1991, 2025))):
        cc = yr[yr.year.between(*yrs) & ~yr.incomplete]
        t2 = D.trend(cc.year, cc[col])
        trends.append({"series": name, "years": f"{yrs[0]}–{yrs[1]}", **t2, "role": "secondary",
                       "sentence": D.trend_sentence(t2, 0.0, "drier/lower", "wetter/higher", "no trend")})

    # (2) season sinking against the season's water balance
    ym = pd.read_csv(X138 / "radar_year_means.csv")
    m = ym.merge(yr[["year", "cwb_season_mm", "rain_season_mm"]], on="year")
    rho, p = stats.spearmanr(m.cwb_season_mm, m.net_mean_cm)
    sv_sent = ("seasons with a drier water balance sank further" if rho > 0 and p < 0.05 else
               "the opposite: seasons with a drier water balance sank less" if rho < 0 and p < 0.05 else
               "a link between how far a season sank and its water balance is not shown with nine seasons")
    corr = [{"pair": "radar two-track net vs CWB May–Nov", "n": len(m), "rho": rho, "p": p, "role": "primary"}]
    for name, x, yv in (("radar largest drop vs CWB May–Nov", m.cwb_season_mm, m.drop_mean_cm),
                        ("radar net vs rain May–Nov", m.rain_season_mm, m.net_mean_cm)):
        r2, p2 = stats.spearmanr(x, yv)
        corr.append({"pair": name, "n": len(m), "rho": r2, "p": p2, "role": "secondary"})
    six = m[m.year.isin([2017, 2018, 2019, 2020, 2021, 2025])]
    r2, p2 = stats.spearmanr(six.cwb_season_mm, six.net_mean_cm)
    corr.append({"pair": "radar net vs CWB May–Nov, six-day years", "n": len(six), "rho": r2, "p": p2, "role": "secondary"})
    x137 = pd.read_csv(X137)
    lz = x137[(x137.reference == "hard_targets") & (x137.partial.isna())].groupby("year").laser_net_cm.mean()
    lm_ = yr.set_index("year").loc[lz.index, "cwb_season_mm"]
    r2, p2 = stats.spearmanr(lm_, lz)
    corr.append({"pair": "laser season net (X-137) vs CWB May–Nov", "n": len(lz), "rho": r2, "p": p2, "role": "secondary"})
    corr = pd.DataFrame(corr)

    # (3) the laser across years, datum B
    lm, laser_day = laser_multiyear(prot)
    au = lm[(lm.window == "autumn") & lm.laser_cm.notna()].set_index("year")
    sp = lm[(lm.window == "spring") & lm.laser_cm.notna()].set_index("year")
    ch_au = float(au.laser_cm.get(2025, np.nan) - au.laser_cm.get(2022, np.nan))
    ch_sp = float(sp.laser_cm.get(2025, np.nan) - sp.laser_cm.get(2023, np.nan))
    def word(v): return "lower" if v < 0 else "higher"
    laser_sent = (f"on the laser, the mat at P6 stood {abs(ch_au):.1f} cm {word(ch_au)} in autumn 2025 than in autumn "
                  f"2022, and {abs(ch_sp):.1f} cm {word(ch_sp)} in spring 2025 than in spring 2023"
                  if np.isfinite(ch_au) and np.isfinite(ch_sp) else "the laser windows are not all readable")
    winters = pd.read_csv(X099_WINTERS)
    winters = winters[(winters.estimate == "plain chain") & (winters.reference == "hard_targets")]
    wsteps = winters.groupby("winter").laser_change_mm.first().div(10).rename("laser_winter_step_cm").reset_index()

    # (4) the year table
    rs = pd.read_csv(X138 / "radar_seasons.csv")
    rs = rs[rs.reference == "hard_targets"]
    rows = []
    wt = pd.read_csv(X138 / "water_table.csv")
    ow = pd.read_csv(X138 / "open_water.csv").set_index("year")
    era = pd.read_csv(X138 / "water_climate.csv").set_index("year")
    yri = yr.set_index("year")
    ymi = ym.set_index("year")
    for y in range(2017, 2027):
        r = {"year": y}
        for tr_ in ("ascending", "descending"):
            q = rs[(rs.year == y) & (rs.track == tr_)]
            if len(q):
                q = q.iloc[0]
                r[f"{tr_[:3]}_window"] = f"{pd.Timestamp(q.start):%d %b} → {pd.Timestamp(q.end):%d %b}"
                r[f"{tr_[:3]}_net_cm"] = q.net_cm
                r[f"{tr_[:3]}_laser_cm"] = q.laser_net_cm
                r["revisit"] = q.revisit
                r["error_typical_cm"] = q.net_typical_cm
                r["error_worst_cm"] = q.net_worst_cm
        r["mean_net_cm"] = ymi.net_mean_cm.get(y, np.nan)
        w6 = wt[(wt["plot"] == "P6") & (wt.year == y)]
        wm = wt[(wt["plot"] == "mat_mean") & (wt.year == y)]
        r["wt_P6_mean_cm"] = w6.mean_cm.iloc[0] if len(w6) else np.nan
        r["wt_P6_mayoct_cm"] = w6.mayoct_mean_cm.iloc[0] if len(w6) else np.nan
        r["wt_P6_flag"] = (w6.flag.iloc[0] if len(w6) and isinstance(w6.flag.iloc[0], str) else "")
        r["wt_mat_mean_cm"] = wm.mean_cm.iloc[0] if len(wm) else np.nan
        r["rain_mm"] = yri.rain_mm.get(y, np.nan)
        r["rain_pct"] = yri.rain_pct.get(y, np.nan)
        r["rain_class"] = yri.rain_class.get(y, "")
        r["rain_mayNov_mm"] = yri.rain_season_mm.get(y, np.nan)
        r["rain_mayNov_class"] = yri.rain_season_class.get(y, "")
        r["cwb_mayNov_anom_mm"] = yri.cwb_season_anom_mm.get(y, np.nan)
        r["era5_rain_mm"] = era.era5_p_mm.get(y, np.nan)
        n = ow.scenes_used.get(y, 0)
        r["open_water_ha"] = f"{ow.area_ha.get(y):.1f}" if n >= 5 else "no clear image"
        lsp, lau = sp.laser_rel_cm.get(y, np.nan), au.laser_rel_cm.get(y, np.nan)
        r["laser_spring_rel_cm"], r["laser_autumn_rel_cm"] = lsp, lau
        rows.append(r)
    yt = pd.DataFrame(rows)

    # outputs
    yr.to_csv(OUT / "imgw_yearly.csv", index=False, float_format="%.2f")
    pd.DataFrame([normal]).to_csv(OUT / "imgw_normal.csv", index=False, float_format="%.2f")
    pd.DataFrame(trends).to_csv(OUT / "trends.csv", index=False, float_format="%.5f")
    corr.to_csv(OUT / "season_vs_balance.csv", index=False, float_format="%.4f")
    lm.to_csv(OUT / "laser_multiyear.csv", index=False, float_format="%.3f")
    wsteps.to_csv(OUT / "laser_winter_steps_x099.csv", index=False, float_format="%.2f")
    yt.to_csv(OUT / "year_table.csv", index=False, float_format="%.2f")
    (OUT / "verdict.json").write_text(json.dumps({"long_trend": long_sent, "season_vs_balance": sv_sent,
                                                 "laser_multiyear": laser_sent, "2026_to_date": td2026},
                                                indent=2, ensure_ascii=False))
    figure(yr, normal, ym, lm, laser_day)

    # error-band uncertainty (chi-square interval on a standard deviation from n differences)
    bands = pd.read_csv(X138 / "error_bands.csv")
    bands = bands[bands.reference == "hard_targets"].copy()
    bands["typical_95_lo_cm"] = bands.net_typical_cm * np.sqrt(bands.n_year_tracks / stats.chi2.ppf(0.975, bands.n_year_tracks))
    bands["typical_95_hi_cm"] = bands.net_typical_cm * np.sqrt(bands.n_year_tracks / stats.chi2.ppf(0.025, bands.n_year_tracks))
    x = x137[x137.reference == "hard_targets"][["year", "track", "revisit", "laser_net_cm", "radar_net_cm", "net_diff_cm"]]

    show = yt.copy()
    disp = show[["year", "revisit", "asc_window", "asc_net_cm", "asc_laser_cm", "des_window", "des_net_cm",
                 "des_laser_cm", "mean_net_cm", "error_typical_cm", "laser_autumn_rel_cm", "wt_P6_mean_cm",
                 "wt_P6_mayoct_cm", "wt_P6_flag", "rain_mm", "rain_pct", "rain_class", "rain_mayNov_class",
                 "cwb_mayNov_anom_mm", "era5_rain_mm", "open_water_ha"]]
    yrs_show = yr[yr.year >= 2017][["year", "rain_mm", "rain_pct", "rain_class", "rain_season_mm", "rain_season_pct",
                                    "rain_season_class", "et0_mm", "cwb_mm", "cwb_anom_mm", "cwb_season_anom_mm",
                                    "missing_days"]]
    readme = f"""# X-139 — each year by the national stations, the seasons against their water balance, the laser across years (generated)

Generated by `experiments/drying_x139.py` (`make run E=drying_x139`); do not edit. Pre-registration:
`../insar_full_test/preregistration/X-139.yaml`. Builds on X-138 (`../drying_x138/`).

## Answers (pre-registered sentences)

- **Water balance since 1991 (IMGW Piła + Poznań):** {long_sent}.
- **Season sinking against its water balance (2017–2025):** {sv_sent} (Spearman rho {rho:.2f}, p {p:.3f}, n {len(m)}).
- **The mat across years (laser, datum B):** {laser_sent}.
- **2026 so far ({td2026['window']}):** rain {td2026['rain_pct']:.0f} % of the 1991–2020 normal for the same window
  ({D.rain_class(td2026['rain_pct'])}); water balance {td2026['cwb_anom_mm']:+.0f} mm against the normal, the
  {td2026['rank_driest_cwb_since_1991']}{'st' if td2026['rank_driest_cwb_since_1991'] == 1 else 'th'} driest such window since 1991.

![summary](fig_x139.png)

## One row per year

Season windows are each track's first → last Sentinel-1 overpass inside May–November. Nets in cm (negative = the
surface went down between the two dates). Laser season nets are X-137's (its own window); laser autumn level is the
15 Oct – 15 Nov median on datum B relative to autumn 2022. Water table in cm (negative = below the floating surface).
Rain: IMGW two-station mean; % of 1991–2020; class by Kaczorowska. CWB = rain − ET0 (Hargreaves).

{md(disp)}

## Why there is no "how much since 2017" from the radar

Adding the seasons (each year's May → November net) is not a change since 2017: the mat **rises back every winter**,
and the radar cannot measure that rise (it loses the level over each winter). The laser's own winter rises (X-099,
plain chain epochs, cm):

{md(wsteps)}

The only change across years that is measured is the laser on one datum (from August 2022): see the answer above and
`laser_multiyear.csv`. Before August 2022 the laser is on another datum (re-levelled in the July 2022 gap) and before
May 2021 there is no laser.

## How the radar's error per year is made (X-138)

For each laser year and track, X-137 computed d = radar season net − laser season net at the same overpasses:

{md(x)}

The error band of a revisit class is the root-mean-square of its d (the typical error) and the largest |d| (the worst
seen): six-day class = 2021 and 2025 (4 values), twelve-day class = 2022–2024 (6 values, 2022 after the re-levelling
only). A year without a laser (2017–2020) gets the band of its class — the assumption is that the radar errs the same
way in years of the same revisit. Nothing else (water table, station, weather) enters the error. With 4 or 6 values the
root-mean-square is itself uncertain (95 % chi-square interval on a standard deviation):

{md(bands[['revisit', 'n_year_tracks', 'net_typical_cm', 'typical_95_lo_cm', 'typical_95_hi_cm', 'net_worst_cm']])}

## IMGW years against 1991–2020

Normal (1991–2020, {normal['years']} years): rain {normal['rain_mm']:.0f} mm/yr, May–November {normal['rain_season_mm']:.0f} mm;
ET0 {normal['et0_mm']:.0f} mm/yr; CWB {normal['cwb_mm']:.0f} mm/yr, May–November {normal['cwb_season_mm']:.0f} mm.

{md(yrs_show, "{:.0f}")}

## Trends and correlations

{md(pd.DataFrame(trends)[['series', 'years', 'n', 'slope', 'slope_lo', 'slope_hi', 'p', 'sentence']], "{:.3g}")}

{md(corr, "{:.3g}")}
"""
    (OUT / "README.md").write_text(readme)
    print(readme)


if __name__ == "__main__":
    main()

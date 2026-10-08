"""X-081 — the data cube 2017–2026 analysed: what it holds, how good it is, and what it can support (exploratory).

Reads only the cube (``06_data/cube``, C-053…C-057, C-060) and the C-059 acquisition list; computes everything:

1. **Coverage.** Monthly observations of every source, 2017 → today (radar, optical, thermal, weather, field truth).
2. **Radar quality through time.** Zone coherence of every consecutive and 12-day pair (the cube's per-pair zone
   medians), by year, track and revisit; the 12-day pairs exist in every year, so they compare all ten years on one
   footing; months; frozen / snow pairs; unwrapped share.
3. **The mat vs the grassland, per year** (H2's paired difference, re-done on every year the cube holds; exploratory —
   the registered H2 numbers stay those of T05).
4. **Backscatter.** Yearly VV / VH of the mat and grassland; satellite offsets (S1B vs S1A, S1C vs S1A, S1D vs S1C).
5. **Drivers.** ERA5 against the site station where both exist (2017–2019 rely on ERA5 alone); Sentinel-2 from the
   older archive against collection 1 over the grassland.
6. **Model-ready layer.** Rows per period and the share with each kind of truth.

Field values are read from the hub and only aggregates are written (to the hub). No coordinate or field value here.

    PYTHONPATH=src python scripts/cube_analysis_x081.py
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import xarray as xr  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from scipy.stats import pearsonr, wilcoxon  # noqa: E402

from insar_wetlands import field  # noqa: E402

HUB = field.field_root().parents[1]
CUBE = HUB / "06_data" / "cube"
SIL, GOLD, RAW = CUBE / "silver", CUBE / "gold", CUBE / "raw"
OUT = HUB / "08_deliverables" / "cube_analysis_x081"
ASF = HUB / "06_data" / "s1_archive" / "asf_acquisitions.csv"
TRACKS = ("ascending", "descending")
ZONES = {"mat": "mat (A)", "grassland": "grassland (C)", "lake": "lake (B)", "stable": "stable ground", "hard": "hard targets"}

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, MUTED, GRIDC, AXIS, SURFACE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"
SEQ = LinearSegmentedColormap.from_list("seq_blue", ["#f0efec", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])
ZCOL = {"mat": AQUA, "grassland": ORANGE, "stable": BLUE}
plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE, "axes.edgecolor": AXIS,
    "axes.labelcolor": INK2, "axes.titlecolor": INK, "axes.titlesize": 9.5, "axes.titleweight": "bold",
    "axes.titlelocation": "left", "axes.labelsize": 8.5, "xtick.color": MUTED, "ytick.color": MUTED,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "axes.grid": True, "grid.color": GRIDC, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False, "legend.fontsize": 7.5,
})


def md(df: pd.DataFrame, fmt="{:.2f}") -> str:
    f = lambda v: fmt.format(v) if isinstance(v, (float, np.floating)) and np.isfinite(v) else ("—" if isinstance(v, float) else str(v))  # noqa: E731
    return "\n".join(["| " + " | ".join(map(str, df.columns)) + " |", "|" + "---|" * len(df.columns)]
                     + ["| " + " | ".join(f(v) for v in r) + " |" for r in df.itertuples(index=False)])


def fp(p: float) -> str:
    return "—" if not np.isfinite(p) else ("< 0.001" if p < 0.001 else f"{p:.3f}")


# ------------------------------------------------------------------------------------------- loaders

def load_pairs() -> pd.DataFrame:
    """One row per (track, 40 m pair): the cube's per-pair zone medians, flags and drivers."""
    cols = ["period", "track", "grid", "pair", "t1", "t2", "revisit_days", "platforms"] + \
           [f"ref_{z}_coh" for z in ZONES] + [f"flag_{f}_{t}" for f in ("frozen_air", "snow_om", "snow_imgw", "wet_any")
                                              for t in ("t1", "t2")]
    p = pd.read_csv(GOLD / "plots_pair.csv.gz", usecols=cols, low_memory=False)
    p = p[p.grid == "40m"].drop_duplicates(["track", "pair"]).copy()
    for c in ("t1", "t2"):
        p[c] = pd.to_datetime(p[c], format="ISO8601", utc=True)
    for c in [c for c in p.columns if c.startswith("flag_")]:
        p[c] = p[c].map({True: True, False: False, "True": True, "False": False}).fillna(False).astype(bool)
    p["year"] = p.t1.dt.year
    p["month"] = p.t1.dt.month
    p["cold"] = p[[c for c in p.columns if c.startswith(("flag_frozen", "flag_snow"))]].any(axis=1)
    p["delta"] = p.ref_mat_coh - p.ref_grassland_coh
    # consecutive = no acquisition of the same track strictly between t1 and t2
    acq = pd.read_csv(ASF, parse_dates=["date"])
    p["consecutive"] = False
    for t in TRACKS:
        d = np.sort(acq[acq.track == t].date.values.astype("datetime64[D]"))
        sel = p.track == t
        a = p.loc[sel, "t1"].dt.tz_convert(None).values.astype("datetime64[D]")
        b = p.loc[sel, "t2"].dt.tz_convert(None).values.astype("datetime64[D]")
        p.loc[sel, "consecutive"] = np.searchsorted(d, b) - np.searchsorted(d, a, side="right") == 0
    return p


def load_site() -> pd.DataFrame:
    s = pd.read_csv(SIL / "site_hourly.csv.gz", index_col=0, low_memory=False)
    s.index = pd.to_datetime(s.index, format="ISO8601", utc=True)
    return s


# ------------------------------------------------------------------------------------------- 1 coverage

def coverage(p: pd.DataFrame, site: pd.DataFrame) -> pd.DataFrame:
    months = pd.period_range("2017-01", pd.Timestamp.now(tz="UTC").tz_localize(None).to_period("M"), freq="M")
    rows = {}
    acq = pd.read_csv(ASF, parse_dates=["date"])
    for t in TRACKS:
        rows[f"S1 acquisitions, {t}"] = acq[acq.track == t].date.dt.to_period("M").value_counts()
        c = p[(p.track == t) & p.consecutive]
        rows[f"consecutive pairs, {t}"] = c.t2.dt.tz_convert(None).dt.to_period("M").value_counts()
        r = xr.open_dataset(SIL / f"rtc_{t}.nc")
        rows[f"backscatter (RTC), {t}"] = pd.Series(pd.DatetimeIndex(r.time.values).to_period("M")).value_counts()

    def times(f):
        return pd.Series(pd.DatetimeIndex(xr.open_dataset(SIL / f).time.values).to_period("M")).value_counts()
    rows["Sentinel-2 scenes"] = times("s2_40m.nc")
    rows["Landsat scenes"] = times("lst_40m.nc")
    rows["ECOSTRESS scenes"] = times("eco_40m.nc")
    n = xr.open_dataset(SIL / "nisar_lband.nc")
    rows["NISAR L-band pairs"] = pd.Series(pd.DatetimeIndex(n.t2.values).to_period("M")).value_counts()
    hm = site.index.tz_convert(None).to_period("M")
    for lab, col in (("ERA5 weather (hours)", "om_era5_land_temperature_2m"), ("site station (hours)", "st_air_c"),
                     ("IMGW synop (hours)", "imgw_dew_h"), ("water table P6 (hours)", "wtd_P6_cm"),
                     ("laser, QC'd (hours)", "laser_ok")):
        v = site[col]
        ok = v.map({True: True, False: False, "True": True, "False": False}).fillna(False) if col == "laser_ok" else v.notna()
        rows[lab] = pd.Series(hm[ok.to_numpy()]).value_counts()
    cov = pd.DataFrame({k: v.reindex(months, fill_value=0) for k, v in rows.items()}).T
    return cov


def fig_coverage(cov: pd.DataFrame, out):
    norm = cov.div(cov.max(axis=1).replace(0, 1), axis=0)
    fig, ax = plt.subplots(figsize=(13, 5.6))
    ax.imshow(norm.to_numpy(), aspect="auto", cmap=SEQ, vmin=0, vmax=1, interpolation="nearest")
    ax.set_yticks(range(len(cov))), ax.set_yticklabels(cov.index, fontsize=7.5)
    yrs = [i for i, m in enumerate(cov.columns) if m.month == 1]
    ax.set_xticks(yrs), ax.set_xticklabels([str(cov.columns[i].year) for i in yrs])
    ax.grid(False)
    for i in yrs:
        ax.axvline(i - 0.5, color=SURFACE, lw=1.2)
    ax.set_title("What the cube holds, month by month (colour: share of the source's busiest month; white = none)")
    fig.tight_layout()
    fig.savefig(out / "fig1_coverage.png", dpi=170)
    plt.close(fig)


# ------------------------------------------------------------------------------------------- 2 radar

def yearly_coherence(p: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (t, rv, y), g in p[p.revisit_days.isin([6, 12])].groupby(["track", "revisit_days", "year"]):
        warm = g[~g.cold]
        rows.append({"track": t, "revisit_days": rv, "year": y, "n": len(g), "n_unfrozen": len(warm),
                     **{f"coh_{z}": g[f"ref_{z}_coh"].median() for z in ("mat", "grassland", "stable")},
                     **{f"coh_{z}_unfrozen": warm[f"ref_{z}_coh"].median() for z in ("mat", "grassland")}})
    return pd.DataFrame(rows)


def fig_coherence(yc: pd.DataFrame, mon: pd.DataFrame, out):
    fig, axs = plt.subplots(2, 3, figsize=(13, 7.2), sharey=True)
    for i, t in enumerate(TRACKS):
        for j, rv in enumerate((6, 12)):
            ax = axs[i, j]
            g = yc[(yc.track == t) & (yc.revisit_days == rv)].set_index("year").reindex(range(2017, 2027))
            for z in ("mat", "grassland", "stable"):   # reindexed: years without such pairs break the line
                ax.plot(g.index, g[f"coh_{z}"], "o-", color=ZCOL[z], ms=4.5, lw=1.4, label=ZONES[z])
            g = g.dropna(subset=["n"])
            for y, n in zip(g.index, g.n):
                ax.text(y, 0.02, str(int(n)), ha="center", fontsize=6.5, color=MUTED)
            ax.set_ylim(0, 1), ax.set_xlim(2016.5, 2026.5)
            ax.set_title(f"{t} · {rv}-day pairs · yearly median", fontsize=8.5)
        ax = axs[i, 2]
        for era, mk in (("2017–2021", "o"), ("2022–2024", "s"), ("2025–2026", "^")):
            g = mon[(mon.track == t) & (mon.era == era)]
            ax.plot(g.month, g.coh_mat, mk + "-", color=AQUA, ms=4, lw=1.2, label=f"mat {era}")
            ax.plot(g.month, g.coh_grassland, mk + "--", color=ORANGE, ms=4, lw=1.0, label=f"grassland {era}")
        ax.set_xticks(range(1, 13)), ax.set_xlim(0.5, 12.5)
        ax.set_title(f"{t} · 12-day pairs · by month", fontsize=8.5)
    axs[0, 0].legend(loc="lower left"), axs[0, 2].legend(loc="lower left", fontsize=6.3, ncol=2)
    for ax in axs[:, 0]:
        ax.set_ylabel("coherence")
    fig.suptitle("Radar coherence through ten years (grey numbers: pairs per year; gaps: no such pairs that year)", x=0.01, ha="left",
                 fontsize=11, fontweight="bold")
    fig.tight_layout()
    fig.savefig(out / "fig2_coherence_years.png", dpi=170)
    plt.close(fig)


def monthly(p: pd.DataFrame) -> pd.DataFrame:
    g = p[p.revisit_days == 12].copy()
    g["era"] = np.where(g.year <= 2021, "2017–2021", np.where(g.year <= 2024, "2022–2024", "2025–2026"))
    return (g.groupby(["track", "era", "month"]).agg(n=("pair", "size"), coh_mat=("ref_mat_coh", "median"),
                                                     coh_grassland=("ref_grassland_coh", "median")).reset_index())


def h2_by_year(p: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (t, rv, y), g in p[p.revisit_days.isin([6, 12])].groupby(["track", "revisit_days", "year"]):
        d = g.delta.dropna()
        pv = wilcoxon(d).pvalue if len(d) >= 8 and (d != 0).any() else np.nan
        rows.append({"track": t, "revisit_days": rv, "year": y, "n": len(d), "median_delta": d.median(),
                     "share_mat_lower": float((d < 0).mean()), "wilcoxon_p": pv})
    return pd.DataFrame(rows)


def fig_h2(h2: pd.DataFrame, p: pd.DataFrame, out):
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.3), sharey=True)
    for ax, t in zip(axs, TRACKS):
        for k, rv in enumerate((6, 12)):
            g = p[(p.track == t) & (p.revisit_days == rv)]
            yrs = sorted(g.year.unique())
            data = [g[g.year == y].delta.dropna().to_numpy() for y in yrs]
            pos = np.array(yrs) + (k - 0.5) * 0.36
            b = ax.boxplot(data, positions=pos, widths=0.3, patch_artist=True, showfliers=False, manage_ticks=False,
                           medianprops=dict(color=INK, lw=1.2), whiskerprops=dict(color=MUTED), capprops=dict(color=MUTED))
            for box in b["boxes"]:
                box.set(facecolor=BLUE if rv == 6 else ORANGE, alpha=0.55, edgecolor="none")
        ax.axhline(0, color=INK2, lw=0.8)
        ax.set_xticks(range(2017, 2027)), ax.set_xlim(2016.5, 2026.5)
        ax.set_title(f"{t} · per pair · blue 6-day, orange 12-day", fontsize=8.5)
    axs[0].set_ylabel("Δ coherence, mat − grassland")
    fig.suptitle("The mat as a distinct radar target, year by year (below 0 = the mat less coherent)", x=0.01,
                 ha="left", fontsize=11, fontweight="bold")
    fig.tight_layout()
    fig.savefig(out / "fig3_mat_vs_grassland.png", dpi=170)
    plt.close(fig)


# ------------------------------------------------------------------------------------------- 4 backscatter

def backscatter() -> tuple[pd.DataFrame, pd.DataFrame]:
    st = xr.open_dataset(SIL / "static_40m.nc")
    masks = {"mat": st.zone.values == 1, "grassland": st.zone.values == 3}
    yearly, plat = [], []
    for t in TRACKS:
        r = xr.open_dataset(SIL / f"rtc_{t}.nc")
        tm = pd.DatetimeIndex(r.time.values)
        per = pd.DataFrame({"time": tm, "platform": r.platform.values.astype(str)})
        for z, m in masks.items():
            per[f"vv_{z}"] = np.nanmedian(r.vv_db.values[:, m], axis=1)
            per[f"vh_{z}"] = np.nanmedian(r.vh_db.values[:, m], axis=1)
        per["year"] = per.time.dt.year
        y = per.groupby("year").agg(n=("time", "size"), **{c: (c, "median") for c in per.columns if c.startswith(("vv_", "vh_"))})
        yearly.append(y.reset_index().assign(track=t))
        for (yr), g in per.groupby("year"):
            ps = sorted(g.platform.unique())
            for a, b in (("S1A", "S1B"), ("S1A", "S1C"), ("S1C", "S1D")):
                if a in ps and b in ps:
                    ga, gb = g[g.platform == a], g[g.platform == b]
                    plat.append({"track": t, "year": yr, "pair": f"{b} − {a}", "n": f"{len(gb)}/{len(ga)}",
                                 "d_vv_grassland_db": gb.vv_grassland.median() - ga.vv_grassland.median(),
                                 "d_vh_grassland_db": gb.vh_grassland.median() - ga.vh_grassland.median()})
    return pd.concat(yearly)[["track", "year", "n", "vv_mat", "vv_grassland", "vh_mat", "vh_grassland"]], pd.DataFrame(plat)


def fig_backscatter(bs: pd.DataFrame, out):
    fig, axs = plt.subplots(1, 2, figsize=(13, 3.8), sharey=True)
    for ax, t in zip(axs, TRACKS):
        g = bs[bs.track == t]
        for z, col in (("mat", AQUA), ("grassland", ORANGE)):
            ax.plot(g.year, g[f"vv_{z}"], "o-", color=col, ms=4, lw=1.4, label=f"VV {ZONES[z]}")
            ax.plot(g.year, g[f"vh_{z}"], "s--", color=col, ms=4, lw=1.0, label=f"VH {ZONES[z]}")
        ax.set_xticks(range(2017, 2027)), ax.set_title(f"{t}: yearly median backscatter (dB)", fontsize=8.5)
    axs[0].legend(loc="center left", ncol=2), axs[0].set_ylabel("γ⁰ (dB)")
    fig.tight_layout()
    fig.savefig(out / "fig4_backscatter_years.png", dpi=170)
    plt.close(fig)


# ------------------------------------------------------------------------------------------- 5 drivers

def era5_vs_station(site: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for lab, a, b in (("air temperature (°C)", "st_air_c", "om_era5_land_temperature_2m"),
                      ("relative humidity (%)", "st_rh_pct", "om_era5_land_relative_humidity_2m")):
        d = site[[a, b]].apply(pd.to_numeric, errors="coerce").dropna()
        rows.append({"variable": lab, "resolution": "hourly", "n": len(d), "r": d[a].corr(d[b]),
                     "bias_era5_minus_station": (d[b] - d[a]).mean(), "rmse": float(np.sqrt(((d[b] - d[a]) ** 2).mean()))})
    d = site[["st_rain_mm", "om_era5_precipitation"]].apply(pd.to_numeric, errors="coerce")
    d = d[d.st_rain_mm.notna()].resample("D").sum(min_count=20).dropna()
    rows.append({"variable": "precipitation (mm/day)", "resolution": "daily", "n": len(d),
                 "r": d.st_rain_mm.corr(d.om_era5_precipitation),
                 "bias_era5_minus_station": (d.om_era5_precipitation - d.st_rain_mm).mean(),
                 "rmse": float(np.sqrt(((d.om_era5_precipitation - d.st_rain_mm) ** 2).mean()))})
    a = site[["st_air_c", "om_era5_land_temperature_2m"]].apply(pd.to_numeric, errors="coerce").dropna()
    agree = float(((a.st_air_c < 0) == (a.om_era5_land_temperature_2m < 0)).mean())
    rows.append({"variable": "frozen air (< 0 °C), same verdict", "resolution": "hourly", "n": len(a), "r": np.nan,
                 "bias_era5_minus_station": np.nan, "rmse": np.nan, "agreement": agree})
    return pd.DataFrame(rows)


def s2_archives() -> pd.DataFrame:
    sc = pd.read_csv(RAW / "s2" / "scenes.csv")
    sc["time"] = pd.to_datetime(sc.time, format="ISO8601", utc=True)
    s2 = xr.open_dataset(SIL / "s2_40m.nc")
    st = xr.open_dataset(SIL / "static_40m.nc")
    grass = st.zone.values == 3
    t = pd.DatetimeIndex(s2.time.values).tz_localize("UTC")
    ndvi = np.nanmedian(s2.ndvi.values[:, grass], axis=1)
    ndwi = np.nanmedian(s2.ndwi.values[:, grass], axis=1)
    valid = np.nanmedian(s2.valid_frac.values[:, grass], axis=1)
    d = pd.DataFrame({"time": t, "ndvi_grass": ndvi, "ndwi_grass": ndwi, "valid": valid})
    d["time"] = d.time.astype("datetime64[ns, UTC]")
    sc["time"] = sc.time.astype("datetime64[ns, UTC]")
    d = pd.merge_asof(d.sort_values("time"), sc[["time", "collection", "baseline"]].sort_values("time"), on="time",
                      tolerance=pd.Timedelta("2min"), direction="nearest")
    d = d[(d.time.dt.month.isin([6, 7, 8])) & (d.valid >= 0.8)]
    d["archive"] = np.where(d.collection == "sentinel-2-c1-l2a", "collection 1", "older L2A")
    return (d.groupby([d.time.dt.year.rename("year"), "archive"])
            .agg(n=("time", "size"), ndvi_grassland=("ndvi_grass", "median"), ndwi_grassland=("ndwi_grass", "median"))
            .reset_index())


def weather_by_year(site: pd.DataFrame, bs: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """ERA5 per year (the only weather source for every year) next to the yearly backscatter, and their correlation."""
    w = site[["om_era5_precipitation", "om_era5_land_soil_moisture_0_to_7cm", "om_era5_land_temperature_2m",
              "om_era5_land_snow_depth"]].apply(pd.to_numeric, errors="coerce")
    full = w.groupby(w.index.year).size() >= 8700                     # whole years only (the current one is partial)
    y = w.groupby(w.index.year).agg(precip_mm=("om_era5_precipitation", "sum"),
                                    soil_moisture=("om_era5_land_soil_moisture_0_to_7cm", "mean"),
                                    frozen_share=("om_era5_land_temperature_2m", lambda v: float((v < 0).mean())),
                                    snow_share=("om_era5_land_snow_depth", lambda v: float((v > 0.01).mean())))
    y = y[full.reindex(y.index).fillna(False).to_numpy()].rename_axis("year").reset_index()
    rows = []
    for t in TRACKS:
        b = bs[bs.track == t].merge(y, on="year")
        for z in ("vv_mat", "vv_grassland"):
            for x in ("soil_moisture", "precip_mm"):
                r, pv = pearsonr(b[x], b[z])
                rows.append({"track": t, "backscatter": z, "against": x, "n_years": len(b), "r": r, "p": pv})
    return y, pd.DataFrame(rows)


# ------------------------------------------------------------------------------------------- 6 gold

def gold_readiness() -> pd.DataFrame:
    g = pd.read_csv(GOLD / "plots_pair.csv.gz", usecols=["period", "grid", "plot", "unw_1x1", "d_wtd_plot_cm",
                                                         "d_laser_los_mm"], low_memory=False)
    g = g[g.grid == "40m"]
    return (g.groupby("period").agg(rows=("plot", "size"), phase=("unw_1x1", lambda s: s.notna().mean()),
                                    water_table=("d_wtd_plot_cm", lambda s: s.notna().mean()),
                                    laser_P6=("d_laser_los_mm", lambda s: s.notna().mean())).reset_index())


# ------------------------------------------------------------------------------------------- README

def readme(cov, p, yc, mon, h2, bs, plat, era, s2a, gr, wy, wr, out):
    yr = cov.T.groupby(cov.columns.year).sum().T
    yr.columns = [str(c) for c in yr.columns]
    yr = yr.reset_index().rename(columns={"index": "source"})
    c12 = yc[yc.revisit_days == 12]
    r17 = c12[c12.year <= 2019]
    r20 = c12[c12.year >= 2020]
    inside = bool(((r17.coh_mat >= r20.coh_mat.min()) & (r17.coh_mat <= r20.coh_mat.max())).all()
                  and ((r17.coh_grassland >= r20.coh_grassland.min()) & (r17.coh_grassland <= r20.coh_grassland.max())).all())
    h12 = h2[h2.revisit_days == 12]
    h6a = h2[(h2.revisit_days == 6) & (h2.track == "ascending")]
    h6d = h2[(h2.revisit_days == 6) & (h2.track == "descending")]
    wet = wy.sort_values("precip_mm", ascending=False).iloc[0]
    rsm = wr[(wr.against == "soil_moisture")]
    rg, rm = rsm[rsm.backscatter == "vv_grassland"], rsm[rsm.backscatter == "vv_mat"]
    brightest_wet = sum(int(bs[bs.track == t].sort_values("vv_grassland").iloc[-1].year == wet.year) for t in TRACKS)
    first = {t: p[p.track == t].t1.min().date() for t in TRACKS}
    cold = p[p.revisit_days == 12].groupby("cold").ref_mat_coh.median()
    q = p[p.revisit_days == 12].assign(season=lambda d: np.where(d.month.between(5, 9), "May–Sep", "Oct–Apr"))
    sea = q.groupby("season")[["ref_mat_coh", "ref_grassland_coh"]].median()
    e = era.set_index("variable")
    s2_17 = int(cov.loc["Sentinel-2 scenes", [m for m in cov.columns if m.year == 2017]].sum())
    eco_first = cov.columns[cov.loc["ECOSTRESS scenes"].to_numpy() > 0][0]
    las_first = cov.columns[cov.loc["laser, QC'd (hours)"].to_numpy() > 0][0]
    wtd = cov.loc["water table P6 (hours)"]
    wtd_span = (wtd.index[wtd.to_numpy() > 24][0], wtd.index[wtd.to_numpy() > 24][-1])
    L = ["# X-081 — the data cube 2017–2026, analysed (exploratory)", "",
         "Generated by `05_code/SAr-adapted/scripts/cube_analysis_x081.py` from `06_data/cube` (after C-060) and the "
         "C-059 acquisition list; every number is computed by it.", "",
         "## What this cube can and cannot support — summary", "",
         f"- **Radar is complete and homogeneous 2017 → 2026.** Interferograms from {first['ascending']} (ascending) and "
         f"{first['descending']} (descending); 12-day pairs exist in every year on both tracks, so radar behaviour can "
         f"be compared across all ten years on one footing (section 2).",
         f"- **2017–2019 {'fall within' if inside else 'fall partly outside'} the spread of 2020–2026** at 12 days (both "
         f"tracks): yearly median coherence of the mat {r17.coh_mat.min():.2f}–{r17.coh_mat.max():.2f} vs "
         f"{r20.coh_mat.min():.2f}–{r20.coh_mat.max():.2f}, of the grassland {r17.coh_grassland.min():.2f}–"
         f"{r17.coh_grassland.max():.2f} vs {r20.coh_grassland.min():.2f}–{r20.coh_grassland.max():.2f}. Year-to-year "
         "differences are as large inside each period as between periods.",
         f"- **At 12 days the mat is less coherent than the grassland in every year on both tracks** (median Δ < 0 in "
         f"{int((h12.median_delta < 0).sum())}/{len(h12)} track-years, Wilcoxon p < 0.05 in {int((h12.wilcoxon_p < 0.05).sum())}/"
         f"{len(h12)}; share of pairs with the mat lower {h12.share_mat_lower.min():.0%}–{h12.share_mat_lower.max():.0%}). "
         f"At 6 days the descending track shows the same ({int((h6d.wilcoxon_p < 0.05).sum())}/{len(h6d)} years), the "
         f"**ascending track does not** (median Δ {h6a.median_delta.min():+.2f} to {h6a.median_delta.max():+.2f}, p < 0.05 "
         f"in {int((h6a.wilcoxon_p < 0.05).sum())}/{len(h6a)} years). Exploratory replication of H2 (registered numbers: T05).",
         f"- **The growing season is when the radar loses the surface:** median 12-day coherence May–September "
         f"{sea.loc['May–Sep', 'ref_mat_coh']:.2f} (mat) / {sea.loc['May–Sep', 'ref_grassland_coh']:.2f} (grassland), "
         f"October–April {sea.loc['Oct–Apr', 'ref_mat_coh']:.2f} / {sea.loc['Oct–Apr', 'ref_grassland_coh']:.2f}; in summer "
         "the two zones nearly converge, so the mat–grassland contrast is mostly a cold-season contrast (figure 2, right).",
         f"- **Frozen or snowy pairs are not the incoherent ones:** 12-day pairs with frozen air or snow at either date "
         f"have a median mat coherence of {cold.get(True, np.nan):.2f}, the others {cold.get(False, np.nan):.2f} "
         "(plausibly because a frozen surface holds still between the two dates).",
         f"- **Yearly backscatter tracks the wetness of the year:** {int(wet.year)} is the wettest full year in ERA5 "
         f"({wet.precip_mm:.0f} mm) and is the brightest VV year on {brightest_wet}/2 tracks. Across "
         f"{int(rsm.n_years.iloc[0])} full years the grassland's yearly VV follows ERA5 topsoil moisture (r "
         f"{rg.r.min():.2f}–{rg.r.max():.2f}, p {fp(rg.p.max())} at most), the mat's more weakly (r {rm.r.min():.2f}–"
         f"{rm.r.max():.2f}, p {fp(rm.p.min())}–{fp(rm.p.max())}). Nothing points to a calibration jump in the older years.",
         f"- **Weather before 2020 rests on ERA5 alone** (no site station). Where both exist ERA5 follows the station "
         f"closely for air temperature (r {e.loc['air temperature (°C)', 'r']:.2f}, bias "
         f"{e.loc['air temperature (°C)', 'bias_era5_minus_station']:+.2f} °C) and the frozen/not-frozen verdict "
         f"({e.loc['frozen air (< 0 °C), same verdict', 'agreement']:.0%} of hours agree), less for humidity "
         f"(r {e.loc['relative humidity (%)', 'r']:.2f}) and daily rain (r {e.loc['precipitation (mm/day)', 'r']:.2f}): "
         "frozen/snow flags transfer to 2017–2019, dew and rain flags less reliably.",
         f"- **Field truth is the bottleneck:** the water table spans {wtd_span[0]} → {wtd_span[1]}, the QC'd laser starts "
         f"{las_first}. 2017–2019 radar, optical and weather are in the cube with no ground truth yet — the coming laser "
         "archive (2017–2026) is what makes them testable.",
         f"- **Thin or partial sources:** Sentinel-2 2017 has {s2_17} usable scenes (older processing); ECOSTRESS starts "
         f"{eco_first}; NISAR L-band covers only 2025–2026; WorldCover (2020, 2021) and LiDAR are single epochs; static "
         "layers deliberately keep their 2020–2025 definitions (C-060).", "",
         "## 1. Coverage", "", "![coverage](fig1_coverage.png)", "", "Observations per year:", "", md(yr, "{:.0f}"), "",
         "## 2. Radar coherence through time", "", "![coherence](fig2_coherence_years.png)", "",
         "Zone medians per pair as stored in the cube (`ref_<zone>_coh`), then the median over the year's pairs. "
         "`_unfrozen`: pairs with no frozen-air or snow flag at either date.", "", md(yc), "",
         "## 3. The mat against the grassland, per year", "", "![h2](fig3_mat_vs_grassland.png)", "",
         "Δ = median coherence of the mat − of the grassland, per pair; Wilcoxon signed-rank test of Δ = 0 within the year.",
         "", md(h2.assign(wilcoxon_p=h2.wilcoxon_p.map(fp))), "",
         "## 4. Backscatter", "", "![backscatter](fig4_backscatter_years.png)", "", md(bs), "",
         "Satellite offsets over the grassland, same year (dB; seasonal sampling balanced by the alternating orbits):", "",
         md(plat), "",
         "Yearly weather (ERA5, full years) and its correlation with the yearly backscatter:", "", md(wy, "{:.3f}"), "",
         md(wr.assign(p=wr.p.map(fp))), "",
         "## 5. Drivers", "", "ERA5 against the site station where both exist:", "", md(era), "",
         "Sentinel-2 summer (June–August) grassland indices by archive, scenes with ≥ 80 % valid grassland pixels "
         "(older L2A and collection 1 overlap in some years — a calibration jump would show here):", "", md(s2a), "",
         "## 6. The model-ready layer", "", "Per-plot pair rows per period (40 m) and the share with each value:", "",
         md(gr), "",
         "## Files", "", "| File | Content |", "|---|---|", "| `coverage_monthly.csv` | section 1 |",
         "| `pairs_zone_coherence.csv` | every pair: zone medians, flags, consecutive, Δ |",
         "| `coherence_by_year.csv`, `coherence_by_month.csv` | section 2 |", "| `mat_vs_grassland_by_year.csv` | section 3 |",
         "| `backscatter_by_year.csv`, `platform_offsets.csv` | section 4 |",
         "| `era5_vs_station.csv`, `s2_archives.csv` | section 5 |", "| `gold_readiness.csv` | section 6 |"]
    (out / "README.md").write_text("\n".join(L) + "\n")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    p = load_pairs()
    site = load_site()
    cov = coverage(p, site)
    yc, mon, h2 = yearly_coherence(p), monthly(p), h2_by_year(p)
    bs, plat = backscatter()
    era, s2a, gr = era5_vs_station(site), s2_archives(), gold_readiness()
    wy, wr = weather_by_year(site, bs)
    cov.to_csv(OUT / "coverage_monthly.csv")
    p.drop(columns=[c for c in p.columns if c.startswith("flag_")]).to_csv(OUT / "pairs_zone_coherence.csv", index=False)
    for name, df in (("coherence_by_year", yc), ("coherence_by_month", mon), ("mat_vs_grassland_by_year", h2),
                     ("backscatter_by_year", bs), ("platform_offsets", plat), ("era5_vs_station", era),
                     ("s2_archives", s2a), ("gold_readiness", gr), ("weather_by_year", wy), ("backscatter_vs_weather", wr)):
        df.to_csv(OUT / f"{name}.csv", index=False)
    fig_coverage(cov, OUT)
    fig_coherence(yc, mon, OUT)
    fig_h2(h2, p, OUT)
    fig_backscatter(bs, OUT)
    readme(cov, p, yc, mon, h2, bs, plat, era, s2a, gr, wy, wr, OUT)
    print(f"X-081 written to {OUT}")


if __name__ == "__main__":
    main()

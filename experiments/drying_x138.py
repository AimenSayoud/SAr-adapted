"""X-138 — is the Rzecin peatland drying, 2017–2026? Three lines, each read on its own, by rules fixed in advance.

(A) How far the mat at P6 sinks within each snow-free season by Sentinel-1 alone (the plain chain of X-099/X-137,
without the laser), 2017–2025, each year–track given the error band of its revisit class from X-137's laser years.
(B) The open-water area of the residual lake each summer (June–September) from the cube's Sentinel-2 scenes (SCL water
class, 20 m), 2017–2026. (C) ERA5 water balance and the water table of the nine plots, as context. Then the
synthesis sentence for the professor. The routines are ``insar_wetlands.drying`` (tested on synthetic answers).

Pre-registration: hub ``08_deliverables/insar_full_test/preregistration/X-138.yaml`` (checked before the run).
Reads the cube and the hub; writes only to the hub ``08_deliverables/drying_x138/``. No field value in this file.

    make run E=drying_x138
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
import xarray as xr  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import p6_validation_x082 as X  # noqa: E402  (the X-082 readers: pairs and laser, unchanged)

from insar_wetlands import drying as D  # noqa: E402
from insar_wetlands import figstyle as fs  # noqa: E402
from insar_wetlands.cube import cube_root  # noqa: E402
from insar_wetlands.cube.sources import SCL_VALID, indices_from_npz  # noqa: E402
from insar_wetlands.paths import hub_root  # noqa: E402
from insar_wetlands.scoring import Protocol, assert_preregistered, chain  # noqa: E402

HUB = hub_root()
CUBE = cube_root()
FT = HUB / "08_deliverables" / "insar_full_test"
OUT = HUB / "08_deliverables" / "drying_x138"
PREREG = FT / "preregistration" / "X-138.yaml"
PROTOCOL = FT / "protocol.yaml"
X137 = FT / "season_total_x137" / "season_totals.csv"
YEARS_A = list(range(2017, 2026))
YEARS_B = list(range(2017, 2027))
SCL_WATER, SCL_SNOW = 6, 11
PRIMARY_REF = "hard_targets"


# --------------------------------------------------------------------------------------------- line A
def radar_seasons(prot: Protocol) -> tuple[pd.DataFrame, dict, pd.DataFrame]:
    p = X.load_pairs()
    las_los, lev = X.laser()
    x137 = pd.read_csv(X137, parse_dates=["start", "end"])
    rows, chains = [], {}
    for track in X.TRACKS:
        cos = float((las_los[track] / lev).dropna().median())          # P6's line-of-sight up-component
        for ref in X.REFS:
            g = p[p.track == track].sort_values("t1")
            c = chain(g.assign(est_mm=g[f"s1_{ref}"]), "est_mm")
            c["time"] = pd.to_datetime(c.time, utc=True)
            chains[(track, ref)] = (c, cos)
            # stop rule: X-137's radar season totals are reproduced between its own dates
            for r in x137[(x137.track == track) & (x137.reference == ref)].itertuples():
                a = c[c.time.dt.date == r.start.date()]
                b = c[c.time.dt.date == r.end.date()]
                a, b = a[a.seg.isin(b.seg)], b[b.seg.isin(a.seg)]
                if len(a) != 1 or len(b) != 1:
                    raise SystemExit(f"stop: X-137 {r.year} {track} {ref} dates not on one segment here")
                net = (float(b.est.iloc[0]) - float(a.est.iloc[0])) / cos / 10
                if abs(net - r.radar_net_cm) > 5e-5:          # X-137 stores 4 decimals (amendment)
                    raise SystemExit(f"stop: X-137 {r.year} {track} {ref}: {net:.6f} here, {r.radar_net_cm:.6f} there")
            for y in YEARS_A + [2026]:
                run = D.season_run(c, y, prot.season_months)
                if len(run) < prot.min_epochs:
                    continue
                s = D.season_total(run, cos)
                rows.append({"year": y, "track": track, "reference": ref, **s,
                             "revisit": D.revisit_class(s["median_revisit_days"]),
                             "chain_gap": c[(c.time.dt.year == y) & c.time.dt.month.between(*prot.season_months)]
                             .seg.nunique() > 1})
    t = pd.DataFrame(rows)
    t.loc[t.year == 2026, "revisit"] = "2026 mixed (not banded)"
    las = x137[["year", "track", "reference", "laser_net_cm", "laser_drop_cm", "net_diff_cm", "drop_diff_cm",
                "revisit", "partial"]]
    t = t.merge(las.drop(columns="revisit"), on=["year", "track", "reference"], how="left")
    bands = []
    for ref in X.REFS:
        for rev in ("six-day", "twelve-day"):
            s = x137[(x137.reference == ref) & (x137.revisit == rev)]
            bn, bd = D.error_band(s.net_diff_cm), D.error_band(s.drop_diff_cm)
            bands.append({"reference": ref, "revisit": rev, "n_year_tracks": bn["n"],
                          "net_typical_cm": bn["typical_cm"], "net_worst_cm": bn["worst_cm"],
                          "drop_typical_cm": bd["typical_cm"], "drop_worst_cm": bd["worst_cm"],
                          "years": ",".join(str(v) for v in sorted(s.year.unique()))})
    bands = pd.DataFrame(bands)
    t = t.merge(bands[["reference", "revisit", "net_typical_cm", "net_worst_cm", "drop_typical_cm", "drop_worst_cm"]],
                on=["reference", "revisit"], how="left")
    return t.sort_values(["reference", "year", "track"]).reset_index(drop=True), chains, bands


def year_means(t: pd.DataFrame, ref: str) -> pd.DataFrame:
    s = t[t.reference == ref].pivot_table(index="year", columns="track", values=["net_cm", "drop_cm"])
    out = pd.DataFrame({"net_mean_cm": s["net_cm"].mean(axis=1), "drop_mean_cm": s["drop_cm"].mean(axis=1),
                        "net_asc_minus_desc_cm": s["net_cm"]["ascending"] - s["net_cm"]["descending"]})
    out["tracks_disagree"] = out.net_asc_minus_desc_cm.abs() > 2.0
    return out.reset_index()


# --------------------------------------------------------------------------------------------- line B
def region_20m() -> np.ndarray:
    z = xr.open_dataset(CUBE / "silver" / "static_40m.nc").zone.values
    return np.kron(np.isin(z, (1, 2)), np.ones((2, 2), bool))           # mat (A) + lake (B), 40 m → 2 × 2 at 20 m


def s2_summers(region: np.ndarray) -> tuple[pd.DataFrame, dict]:
    sc = pd.read_csv(CUBE / "raw" / "s2" / "scenes.csv")
    sc["time"] = pd.to_datetime(sc.time, format="ISO8601", utc=True)
    sc = sc[sc.time.dt.month.between(6, 9) & sc.file.notna() & (sc.file != "")]
    sc["day"] = sc.time.dt.date
    sc = sc.sort_values("valid_share", ascending=False).drop_duplicates("day")   # one scene per acquisition day
    rows, maps = [], {}
    for y in YEARS_B:
        flags = {"scl": [], "idx": [], "valid": []}
        for f in sc[sc.time.dt.year == y].sort_values("time").file:
            z = np.load(CUBE / "raw" / "s2" / f)
            scl = z["scl"]
            if scl.shape != region.shape:
                raise SystemExit(f"stop: {f} grid {scl.shape} is not the template's 20 m grid {region.shape}")
            valid = np.isin(scl, SCL_VALID) & (scl != SCL_SNOW)
            if valid[region].mean() < 0.8:
                continue
            idx = indices_from_npz(z)
            b03 = np.where(z["refl_B03"] == -32768, np.nan, z["refl_B03"] / 1e4)
            b11 = np.where(z["refl_B11"] == -32768, np.nan, z["refl_B11"] / 1e4)
            with np.errstate(invalid="ignore"):
                w_idx = valid & ((idx["mndwi"] > 0) | ((b03 > 0) & (b11 <= 0.02)))
            flags["scl"].append(valid & (scl == SCL_WATER))
            flags["idx"].append(w_idx)
            flags["valid"].append(valid)
        n = len(flags["valid"])
        r = {"year": y, "scenes_used": n}
        if n:
            v = np.stack(flags["valid"])
            for name, key, thr in (("area_ha", "scl", 0.5), ("area_f03_ha", "scl", 0.3), ("area_f07_ha", "scl", 0.7),
                                   ("area_index_ha", "idx", 0.5)):
                a, freq = D.water_area(np.stack(flags[key]), v, region, thr)
                r[name] = a
                if name == "area_ha":
                    maps[y] = np.where(region, freq, np.nan)
            vals = [r[k] for k in ("area_ha", "area_f03_ha", "area_f07_ha", "area_index_ha")]
            r.update(band_lo_ha=min(vals), band_hi_ha=max(vals), band_width_ha=max(vals) - min(vals))
        rows.append(r)
    return pd.DataFrame(rows), maps


def ortho_water() -> pd.DataFrame:
    o = xr.open_dataset(CUBE / "silver" / "ortho_cir_1m.nc")
    st = xr.open_dataset(CUBE / "silver" / "static_40m.nc")
    reg = st.zone.isin([1, 2]).astype("int8").reindex_like(o, method="nearest").values.astype(bool)
    px = float(abs(o.x.values[1] - o.x.values[0]) * abs(o.y.values[1] - o.y.values[0]))
    out = []
    for v in sorted(o.data_vars):
        a = o[v].values
        ok = reg & np.isfinite(a)
        out.append({"year": int(v.split("_")[-1]), "ortho_area_ha": float((ok & (a < 0)).sum() * px / 1e4),
                    "ortho_region_valid_share": float(ok.sum() / reg.sum())})
    return pd.DataFrame(out)


# --------------------------------------------------------------------------------------------- line C
def water_climate() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    plots = [f"P{k}" for k in range(1, 10)]
    cols = (["time_utc", "om_era5_precipitation", "om_era5_et0_fao_evapotranspiration", "st_rain_mm",
             "wtd_mat_mean_cm"] + [f"wtd_{p}_cm" for p in plots] + [f"censored_{p}" for p in plots])
    d = pd.read_csv(CUBE / "silver" / "site_hourly.csv.gz", usecols=cols, low_memory=False)
    d["time"] = pd.to_datetime(d.time_utc, format="ISO8601", utc=True)
    d["year"], d["month"] = d.time.dt.year, d.time.dt.month
    gs = d.month.between(5, 10)
    clim = []
    for y, g in d.groupby("year"):
        s = g[g.month.between(5, 10)]
        clim.append({"year": y, "hours": len(g), "era5_p_mm": g.om_era5_precipitation.sum(),
                     "era5_et0_mm": g.om_era5_et0_fao_evapotranspiration.sum(),
                     "era5_p_minus_et0_mm": g.om_era5_precipitation.sum() - g.om_era5_et0_fao_evapotranspiration.sum(),
                     "era5_p_mayoct_mm": s.om_era5_precipitation.sum(),
                     "era5_p_minus_et0_mayoct_mm": s.om_era5_precipitation.sum()
                     - s.om_era5_et0_fao_evapotranspiration.sum(),
                     "station_rain_mm": g.st_rain_mm.sum(min_count=1), "station_rain_hours": int(g.st_rain_mm.count())})
    clim = pd.DataFrame(clim)
    clim["complete_year"] = clim.hours >= 8760
    # 2026 to date against the same months of 2017–2025
    last = d.time.max()
    m = d[(d.month < last.month) | ((d.month == last.month) & (d.time.dt.day <= last.day))]
    todate = m.groupby("year").agg(era5_p_mm=("om_era5_precipitation", "sum"),
                                   era5_et0_mm=("om_era5_et0_fao_evapotranspiration", "sum")).reset_index()
    todate["era5_p_minus_et0_mm"] = todate.era5_p_mm - todate.era5_et0_mm
    todate["window"] = f"01-01 → {last:%m-%d %H:00}"
    wt = []
    for y in range(2020, 2025):
        g, s = d[d.year == y], d[(d.year == y) & gs]
        for p in plots + ["mat_mean"]:
            col = f"wtd_{p}_cm" if p != "mat_mean" else "wtd_mat_mean_cm"
            cen = int(pd.to_numeric(g.get(f"censored_{p}", pd.Series(dtype=float)), errors="coerce").fillna(0).sum())
            wt.append({"year": y, "plot": p, "hours": int(g[col].count()), "mean_cm": g[col].mean(),
                       "mayoct_mean_cm": s[col].mean(), "min_cm": g[col].min(), "max_cm": g[col].max(),
                       "censored_hours": cen, "flag": "true value lower" if cen else ""})
    return clim, todate, pd.DataFrame(wt)


# --------------------------------------------------------------------------------------------- figures
def figures(t, ym, bands, ow, maps, clim, wt) -> None:
    prim = t[(t.reference == PRIMARY_REF) & (t.year <= 2025)]
    with fs.style():
        fig, axes = plt.subplots(3, 1, figsize=fs.size("full", 200), gridspec_kw={"hspace": 0.55, "left": 0.09,
                                 "right": 0.90, "top": 0.92, "bottom": 0.06})
        ax = axes[0]
        ax.axhline(0, color=fs.INK, lw=0.6)
        for track, dx in (("ascending", -0.15), ("descending", 0.15)):
            s = prim[prim.track == track]
            ax.errorbar(s.year + dx, s.net_cm, yerr=s.net_typical_cm, fmt="o", ms=4, lw=0.8, capsize=2,
                        color=fs.TRACK[track], label=f"radar, {track} (bar: typical error of its revisit)")
            ax.errorbar(s.year + dx, s.net_cm, yerr=s.net_worst_cm, fmt="none", lw=0.4, color=fs.TRACK[track], alpha=0.5)
            ax.plot(s.year + dx, s.laser_net_cm, ls="none", marker="_", ms=9, mew=1.6, color=fs.CLASS["laser"],
                    label="laser at the same overpasses" if track == "ascending" else None)
        for y in (2022, 2023, 2024):
            ax.axvspan(y - 0.45, y + 0.45, color=fs.GRID, alpha=0.35, lw=0)
        ax.set_ylabel("May → Nov net change at P6 (cm)")
        ax.set_title("(a) How far the mat sank each snow-free season (Sentinel-1 alone; grey = twelve-day years)",
                     loc="left", pad=22)
        ax.set_xticks(YEARS_A)
        ax.legend(frameon=False, ncol=3, loc="lower center", bbox_to_anchor=(0.5, 1.12), fontsize=7)
        ax = axes[1]
        o = ow[ow.scenes_used >= 5]                       # a summer with fewer clear scenes is drawn hollow, unjoined
        thin = ow[(ow.scenes_used > 0) & (ow.scenes_used < 5)]
        ax.fill_between(o.year, o.band_lo_ha, o.band_hi_ha, color=fs.CLASS["open water"], alpha=0.3, lw=0,
                        label="band (thresholds, index rule)")
        ax.plot(o.year, o.area_ha, "o-", ms=4, color=fs.CLASS["open water"], label="Sentinel-2 water class, June–Sept")
        for r in thin.itertuples():
            ax.plot(r.year, r.area_ha, "o", ms=5, mfc="white", color=fs.CLASS["open water"])
            ax.annotate(f"{r.scenes_used} clear scene(s)", (r.year, r.area_ha), xytext=(4, 6),
                        textcoords="offset points", fontsize=7, color=fs.MUTED)
        ax.set_ylim(bottom=0)
        ax.set_ylabel("open water (ha)")
        ax.set_xticks(YEARS_B)
        ax.set_title("(b) Open water of the lake each summer (mat + lake zones)", loc="left")
        ax.legend(frameon=False, loc="center left", bbox_to_anchor=(0.06, 0.45))
        ax = axes[2]
        c = clim[clim.complete_year & (clim.year <= 2025)]
        ax.bar(c.year, c.era5_p_minus_et0_mm, color=[fs.CLASS["open water"] if v > 0 else fs.CLASS["grassland"]
                                                     for v in c.era5_p_minus_et0_mm], width=0.6)
        ax.axhline(0, color=fs.INK, lw=0.6)
        ax.set_ylabel("ERA5 P − ET0 (mm/yr)")
        ax2 = ax.twinx()
        w = wt[wt["plot"] == "P6"]
        ax2.plot(w.year, w.mean_cm, "o-", ms=4, color=fs.CLASS["mat"], label="P6 water table, yearly mean")
        ax2.set_ylabel("water table at P6 (cm, below surface < 0)")
        ax2.legend(frameon=False, loc="lower right")
        ax.set_xticks(YEARS_A)
        ax.set_title("(c) Water balance (bars) and water table (line)", loc="left")
        fs.export(fig, OUT / "fig_drying", provenance={"script": "experiments/drying_x138.py", "ticket": "X-138"})
        plt.close(fig)
        yrs = sorted(maps)
        fig, axs = plt.subplots(2, 5, figsize=fs.size("full", 85))
        rows, cols = np.where(np.isfinite(maps[yrs[0]]))
        win = (slice(rows.min(), rows.max() + 1), slice(cols.min(), cols.max() + 1))
        for ax, y in zip(axs.ravel(), yrs):
            im = ax.imshow(maps[y][win], vmin=0, vmax=1, cmap="Blues", interpolation="nearest")
            ax.set_title(str(y))
            ax.set_xticks([]), ax.set_yticks([])
        for ax in axs.ravel()[len(yrs):]:
            ax.axis("off")
        fig.colorbar(im, ax=axs, shrink=0.7, label="share of clear summer scenes classed water")
        fs.export(fig, OUT / "fig_water_frequency", provenance={"script": "experiments/drying_x138.py", "ticket": "X-138"})
        plt.close(fig)


# --------------------------------------------------------------------------------------------- README
def md(df: pd.DataFrame, fmt: str = "{:.1f}") -> str:
    cols = list(df.columns)
    out = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for r in df.itertuples(index=False):
        out.append("| " + " | ".join(fmt.format(v) if isinstance(v, (float, np.floating)) and np.isfinite(v)
                                     else ("" if isinstance(v, float) else str(v)) for v in r) + " |")
    return "\n".join(out)


def main() -> None:
    assert_preregistered(PREREG, [OUT / "radar_seasons.csv", OUT / "open_water.csv"])
    prot = Protocol.load(PROTOCOL)
    OUT.mkdir(parents=True, exist_ok=True)

    t, _, bands = radar_seasons(prot)
    ym = year_means(t[t.year <= 2025], PRIMARY_REF)
    six_worst_net = float(bands[(bands.reference == PRIMARY_REF) & (bands.revisit == "six-day")].net_worst_cm.iloc[0])
    six_worst_drop = float(bands[(bands.reference == PRIMARY_REF) & (bands.revisit == "six-day")].drop_worst_cm.iloc[0])
    trends = []
    sent = {}
    for what, col, thr in (("net change", "net_mean_cm", six_worst_net), ("largest drop", "drop_mean_cm", six_worst_drop)):
        tr = D.trend(ym.year, ym[col])
        s = D.trend_sentence(tr, thr, "the seasons' sinking grew over 2017–2025", "the seasons' sinking shrank",
                             "no trend in how far the seasons sank can be shown at this precision")
        sent[what] = s
        trends.append({"line": "A", "series": f"radar two-track mean, {what}", "years": "2017–2025", **tr,
                       "threshold": thr, "sentence": s, "role": "primary"})
    # A secondary
    six = ym[ym.year.isin([2017, 2018, 2019, 2020, 2021, 2025])]
    for what, col, thr in (("net change", "net_mean_cm", six_worst_net), ("largest drop", "drop_mean_cm", six_worst_drop)):
        tr = D.trend(six.year, six[col])
        trends.append({"line": "A", "series": f"radar two-track mean, {what}, six-day years only", "years":
                       "2017–2021, 2025", **tr, "threshold": thr, "role": "secondary",
                       "sentence": D.trend_sentence(tr, thr, "grew", "shrank", "no trend")})
    for ref in X.REFS:
        for track in X.TRACKS:
            s = t[(t.reference == ref) & (t.track == track) & (t.year <= 2025)]
            tr = D.trend(s.year, s.net_cm)
            trends.append({"line": "A", "series": f"radar {track}, {ref}, net change", "years": "2017–2025", **tr,
                           "threshold": six_worst_net, "role": "secondary",
                           "sentence": D.trend_sentence(tr, six_worst_net, "grew", "shrank", "no trend")})
    a_dry = True if any("grew" in v for v in sent.values()) else (False if any("shrank" in v for v in sent.values())
                                                                     else None)

    region = region_20m()
    ow, maps = s2_summers(region)
    ow = ow.merge(ortho_water(), on="year", how="left")
    o = ow[ow.scenes_used > 0]
    med_band = float(o.band_width_ha.median())
    tr = D.trend(o.year, o.area_ha)
    b_sent = D.trend_sentence(tr, med_band, "the open water shrank over 2017–2026", "the open water grew over 2017–2026",
                              "no trend in the open water can be shown at this precision")
    trends.append({"line": "B", "series": "Sentinel-2 open water, June–September", "years": "2017–2026", **tr,
                   "threshold": med_band, "sentence": b_sent, "role": "primary"})
    o2 = o[o.year != 2017]
    tr = D.trend(o2.year, o2.area_ha)
    trends.append({"line": "B", "series": "Sentinel-2 open water, without 2017", "years": "2018–2026", **tr,
                   "threshold": med_band, "role": "secondary",
                   "sentence": D.trend_sentence(tr, med_band, "shrank", "grew", "no trend")})
    b_dry = True if "shrank" in b_sent else (False if "grew" in b_sent else None)

    clim, todate, wt = water_climate()
    c = clim[clim.complete_year & clim.year.between(2017, 2025)]
    tr = D.trend(c.year, c.era5_p_minus_et0_mm)
    c_sent = D.trend_sentence(tr, 0.0, "the water balance became drier", "the water balance became wetter",
                              "no trend in the water balance at this length")
    trends.append({"line": "C", "series": "ERA5 annual P − ET0", "years": "2017–2025", **tr, "threshold": 0.0,
                   "sentence": c_sent, "role": "primary"})
    syn = D.synthesis(a_dry, b_dry)
    trends = pd.DataFrame(trends)

    t.to_csv(OUT / "radar_seasons.csv", index=False, float_format="%.4f")
    ym.to_csv(OUT / "radar_year_means.csv", index=False, float_format="%.4f")
    bands.to_csv(OUT / "error_bands.csv", index=False, float_format="%.4f")
    ow.to_csv(OUT / "open_water.csv", index=False, float_format="%.4f")
    clim.to_csv(OUT / "water_climate.csv", index=False, float_format="%.2f")
    todate.to_csv(OUT / "water_climate_to_date.csv", index=False, float_format="%.2f")
    wt.to_csv(OUT / "water_table.csv", index=False, float_format="%.3f")
    trends.to_csv(OUT / "trends.csv", index=False, float_format="%.5f")
    (OUT / "verdict.json").write_text(json.dumps({"synthesis": syn, "line_A": sent, "line_B": b_sent, "line_C": c_sent},
                                                 indent=2, ensure_ascii=False))
    figures(t, ym, bands, ow, maps, clim, wt)

    prim = t[(t.reference == PRIMARY_REF)].copy()
    prim["error (typical / worst)"] = [f"±{a:.1f} / ±{b:.1f}" if np.isfinite(a) else "—"
                                       for a, b in zip(prim.net_typical_cm, prim.net_worst_cm)]
    tab_a = prim[["year", "track", "revisit", "start", "end", "overpasses", "net_cm", "error (typical / worst)",
                  "laser_net_cm", "drop_cm", "drop_date", "laser_drop_cm"]].rename(
        columns={"laser_net_cm": "laser net (X-137 window)", "laser_drop_cm": "laser drop (X-137 window)"})
    tr_show = trends[["line", "series", "years", "n", "slope", "slope_lo", "slope_hi", "p", "threshold", "sentence"]]
    wt6 = wt[wt["plot"].isin(["P6", "mat_mean"])][["year", "plot", "mean_cm", "mayoct_mean_cm", "min_cm", "max_cm",
                                                 "censored_hours", "flag"]]
    readme = f"""# X-138 — is the Rzecin peatland drying, 2017–2026? (generated)

Generated by `experiments/drying_x138.py` (`make run E=drying_x138`); do not edit. Pre-registration:
`../insar_full_test/preregistration/X-138.yaml` (sentences fixed there before anything was computed). Stop rule
passed: the radar chain read here reproduces X-137's season totals 2021–2025 between X-137's dates.

## Answer for the professor (pre-registered synthesis)

**{syn}.**

- Line A, Sentinel-1 at P6 (how far the mat sank each season): {sent['net change']} (net change); {sent['largest drop']}
  (largest drop).
- Line B, open water of the lake (Sentinel-2): {b_sent}.
- Line C, context (ERA5 water balance 2017–2025): {c_sent}.

What this cannot say: the radar restarts every spring and cannot carry a level across a winter, so there is **no
2017→2026 cumulative level change** from InSAR — only how far each season sank. The 2017–2020 errors are borrowed
from the laser years (2021–2025) of the same revisit, not measured. Shrinking open water in an overgrowing lake can be
succession as well as drying. The water table is depth below a floating surface, which follows the water. One radar
point (P6).

![summary](fig_drying.png)

## Line A — each season by Sentinel-1 alone, with the error of its revisit (P6, hard-target reference, vertical cm)

Error bands from X-137 (radar − laser, same overpasses): typical = root-mean-square, worst = max |difference|.\nThe laser columns are X-137's, over X-137's window (first and last overpass with a laser reading, so a few days\ninside the radar window of this table in some years): they show the size of the truth, not a new comparison.

{md(bands[bands.reference == PRIMARY_REF].drop(columns='reference'))}

{md(tab_a)}

Two-track means per year (tracks disagree when |asc − desc| > 2 cm):

{md(ym)}

## Line B — open water of the lake each summer (Sentinel-2, June–September, mat + lake zones, 20 m)

Primary: pixels classed water (SCL 6) in at least half the clear scenes. Band: frequency 0.3 / 0.7 and the index rule
(MNDWI > 0 or dark SWIR). One scene per acquisition day (the most valid). CIR orthophoto areas beside them: flight
dates differ from the window, printed only.

{md(ow, "{:.2f}")}

![water frequency](fig_water_frequency.png)

## Line C — water balance and water table (context)

ERA5 per calendar year (2026 incomplete) and the site station where it ran:

{md(clim, "{:.0f}")}

2026 to date against the same window of every year:

{md(todate, "{:.0f}")}

Water table at P6 and the mat mean, 2020–2024 (cm, negative = below the floating surface; years ranked only, no
trend with five years). Censored hours: the well was dry at the logger, so the true value was lower.

{md(wt6)}

All nine plots: `water_table.csv`.

## Trends (Theil–Sen slope per year, 95 % interval, Kendall p; sentences by the pre-registered rules)

{md(tr_show, "{:.3g}")}
"""
    (OUT / "README.md").write_text(readme)
    print(readme)


if __name__ == "__main__":
    main()

"""C-071 — the supervisor's table: Sentinel-1 at the nine monitoring plots, one row per plot and acquisition, joined
with the field data (hourly water table and its antecedents, meteorology, the P6 laser, UAV/LAI campaigns).

The supervisor's request (`06_data/field/raw/2026-09_abdallah_delivery/mail prof.md` §2, §4): per plot and acquisition,
phase or displacement, coherence, VV, VH, RVI and the number of valid pixels, with the median or mean and the spread
(SD or IQR) inside a neighbourhood; the water table at the overpass with its 24-hour, 3-day and 7-day antecedents. The
first answer (`field_first`, 2026-09-25) predates the grid corrections (X-120 HyP3, X-122 RTC); this one reads the
corrected cube (silver for the window statistics, gold for everything already joined there).

Which Sentinel-1 data, and why:

* **Interferometric phase and coherence:** the HyP3 burst interferograms on the cube's 40 m grid, relocated under
  their labels (X-120), both tracks; a pixel is valid where its phase is finite — as in every validated read of the
  project (X-065, v4), the unwrapping component is not used as a mask but counted (pixels inside a connected
  component); only the shortest pairs, because only they keep the mat's motion (`SETTLED.md`
  §3, §4). Per acquisition: the consecutive pair ending there (6 days with two satellites, 12 days with one). Every
  pair up to ``MAX_PAIR_DAYS`` is kept in the pair table.
* **Displacement:** −(φ − φ_ref)·λ/4π, + toward the satellite, against the hard targets (the reference that matched
  the laser best, X-072 / X-082) and against the grassland (the first study's). Vertical = LOS / up-component of the
  line of sight, i.e. assuming the motion is vertical.
* **Backscatter:** the Planetary Computer RTC (γ⁰ VV, VH, VH/VV, RVI) re-sampled from its 10 m source at the corrected
  positions (X-122, ``*_placed``).
* **Window:** 3×3 cells (120 m) on the mat (zone A) — the plots are a few tens of metres apart, so a larger window
  repeats the same pixels; the plot's own cell (1×1) beside it.
* **Years:** 2020–2025, the years with field data (water table 2020–2024, laser 27 May 2021 – 2025). The field team's
  common period 2022–2024 is a column.

Not included, on purpose: per-pixel inversion series (SBAS, ISBAS, EVD), which do not recover the mat (`SETTLED.md`
§1); the first study's EVD temporal coherence is carried as a per-plot attribute only, read on the grid before the
X-120 correction (stated in its column).

Writes only to the hub (field values): ``08_deliverables/field_supervisor_workbook_c071/tables/``.

    make run E=supervisor_table_c071
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import xarray as xr

from insar_wetlands import field, plot_table as pt
from insar_wetlands.cube import core
from insar_wetlands.cube.build import plot_units
from insar_wetlands.field_link import season_residual
from insar_wetlands.paths import cube_root, hub_root
from insar_wetlands.scoring import Protocol, guard_lockbox

HUB = hub_root()
CUBE = cube_root()
OUT = HUB / "08_deliverables" / "field_supervisor_workbook_c071" / "tables"
PROTOCOL = HUB / "08_deliverables" / "insar_full_test" / "protocol.yaml"
PERIODS = ("2020_2021", "2022_2024", "2025")
YEARS = (2020, 2025)
COMMON = (2022, 2024)                      # the field team's common study period
MAX_PAIR_DAYS = 48
H = 1                                      # 3×3 window
UAV_MAX_LAG_DAYS = 12                      # nearest acquisition per track, at most one 12-day revisit away
REFS = {"hard": "hard targets", "grassland": "grassland (zone C)"}
SUPERVISOR_GROUP = {**{p: "hydrological reference (P1–P5)" for p in ("P1", "P2", "P3", "P4", "P5")},
                    **{p: "floating mat (P6–P9)" for p in ("P6", "P7", "P8", "P9")}}


def in_years(dates) -> np.ndarray:
    y = pd.to_datetime(pd.Series(dates)).dt.year
    return ((y >= YEARS[0]) & (y <= YEARS[1])).to_numpy()


def windows(px: pd.DataFrame, mat: np.ndarray, shape) -> dict:
    return {p: pt.window_index(r, c, shape, H, mat) for p, r, c in zip(px["plot"], px.row, px.col)}


def stats_cols(s: pd.DataFrame, prefix: str) -> pd.DataFrame:
    return s.rename(columns={k: f"{prefix}_{k}" for k in s.columns})


# ------------------------------------------------------------------------------ pairs

def pair_table(px, mat, gold_pairs: pd.DataFrame) -> pd.DataFrame:
    """Every pair up to MAX_PAIR_DAYS, per plot: window statistics of coherence and of displacement (two references),
    joined with what the cube already holds for the pair (drivers, flags, laser change)."""
    out = []
    for period in PERIODS:
        for track in ("ascending", "descending"):
            ds = xr.open_dataset(CUBE / "silver" / f"ifg_{period}_{track}.nc")
            keep = ds.revisit_days.values <= MAX_PAIR_DAYS
            t1 = pd.DatetimeIndex(ds.t1.values[keep], tz="UTC")
            t2 = pd.DatetimeIndex(ds.t2.values[keep], tz="UTC")
            cons = pt.consecutive_mask(ds.t1.values, ds.t2.values)[keep]
            U = ds.unw.values[keep].astype(float)
            C = ds.coh.values[keep].astype(float)
            K = ds.conncomp.values[keep]
            win = windows(px, mat, U.shape[1:])
            g = gold_pairs[(gold_pairs.period == period) & (gold_pairs.track == track)]
            for p in px["plot"]:
                rr, cc = win[p]
                base = pd.DataFrame({"period": period, "track": track, "pair": ds.pair.values[keep].astype(str),
                                     "plot": p, "t1": t1, "t2": t2, "revisit_days": ds.revisit_days.values[keep],
                                     "consecutive": cons})
                gp = g[g["plot"] == p].set_index("pair")
                j = gp.reindex(base.pair)
                coh = pt.window_stats_stack(C[:, rr, cc])
                unw = pt.window_stats_stack(U[:, rr, cc])
                cols = [base, stats_cols(coh, "coh_3x3")]
                for ref in REFS:
                    los = pt.phase_to_los_mm(unw, j[f"ref_{ref}_unw"].to_numpy())
                    cols.append(stats_cols(los.drop(columns="n_valid"), f"dlos_{ref}_mm_3x3"))
                cols.append(pd.DataFrame({"phase_n_valid_3x3": unw["n_valid"].to_numpy(),
                                          "phase_n_in_component_3x3": (K[:, rr, cc] > 0).sum(1)}))
                one = pd.DataFrame({
                    "coh_1x1": j["coh_1x1"].to_numpy(),
                    "dlos_hard_mm_1x1": -core.rad_to_los_mm(j["unw_1x1"] - j["ref_hard_unw"]),
                    "dlos_grassland_mm_1x1": -core.rad_to_los_mm(j["unw_1x1"] - j["ref_grassland_unw"]),
                    "consistency_3x3": j["consistency_3x3"].to_numpy(),
                    "ref_hard_coh": j["ref_hard_coh"].to_numpy(), "ref_grassland_coh": j["ref_grassland_coh"].to_numpy(),
                    "los_u": j["los_u"].to_numpy(), "platforms": j["platforms"].to_numpy()})
                drv = j[["d_wtd_plot_cm", "d_wtd_P6_cm", "rain_sum_mm", "et0_sum_mm", "d_vv_db_3x3", "d_vh_db_3x3",
                         "d_laser_los_mm", "d_laser_vertical_mm", "laser_cycle_risk",
                         "flag_frozen_air_t1", "flag_frozen_air_t2", "flag_snow_om_t1", "flag_snow_om_t2",
                         "flag_snow_imgw_t1", "flag_snow_imgw_t2", "flag_wet_any_t1", "flag_wet_any_t2",
                         "laser_snow_t1", "laser_snow_t2"]].reset_index(drop=True)
                out.append(pd.concat([c.reset_index(drop=True) for c in cols] + [one.reset_index(drop=True), drv],
                                     axis=1))
    d = pd.concat(out, ignore_index=True)
    d = d[in_years(d.t2.dt.tz_convert(None)) & in_years(d.t1.dt.tz_convert(None))]
    # a pair at a year-set boundary is in two sets; keep one
    d = d.drop_duplicates(["track", "pair", "plot"], keep="first").reset_index(drop=True)
    d["dvert_hard_mm_3x3"] = d["dlos_hard_mm_3x3_median"] / d["los_u"]
    truthy = lambda c: d[c].map({True: True, False: False, "True": True, "False": False}).fillna(False)  # noqa: E731
    d["frozen_or_snow"] = (truthy("flag_frozen_air_t1") | truthy("flag_frozen_air_t2") | truthy("flag_snow_om_t1")
                           | truthy("flag_snow_om_t2") | truthy("flag_snow_imgw_t1") | truthy("flag_snow_imgw_t2"))
    d["wet_surface"] = truthy("flag_wet_any_t1") | truthy("flag_wet_any_t2")
    return d


def laser_fill_crossed(pairs: pd.DataFrame, laser_hourly: pd.DataFrame) -> np.ndarray:
    """True where a filled (interpolated) laser stretch lies between the pair's dates: the change there is not
    measured, and the 2022 re-levelling hides in one (`field_p6/laser_qc/laser_qc.md` §5, §7)."""
    fill = pd.DatetimeIndex(pd.to_datetime(laser_hourly.loc[laser_hourly.laser_filled.astype(bool), "time_utc"],
                                           utc=True))
    f = np.sort(fill.as_unit("ns").asi8)
    a = pd.DatetimeIndex(pairs.t1).as_unit("ns").asi8
    b = pd.DatetimeIndex(pairs.t2).as_unit("ns").asi8
    return np.searchsorted(f, b) - np.searchsorted(f, a) > 0


# ------------------------------------------------------------------------------ acquisitions

def rtc_stats(px, mat) -> pd.DataFrame:
    out = []
    for track in ("ascending", "descending"):
        r = xr.open_dataset(CUBE / "silver" / f"rtc_{track}.nc")
        dates = pd.DatetimeIndex(r.time.values).normalize()
        win = windows(px, mat, r.vv_db.shape[1:])
        for p in px["plot"]:
            rr, cc = win[p]
            cols = [pd.DataFrame({"track": track, "date": dates.date, "plot": p})]
            for v in ("vv_db", "vh_db", "vh_vv_db", "rvi"):
                cols.append(stats_cols(pt.window_stats_stack(r[v].values[:, rr, cc]), f"{v}_3x3"))
            out.append(pd.concat(cols, axis=1))
    d = pd.concat(out, ignore_index=True)
    d["date"] = d["date"].astype(str)
    return d


def anomaly(d: pd.DataFrame, col: str, ok: pd.Series, dt_col: str | None = None) -> pd.Series:
    """The value minus its annual harmonic (and the pair length, for coherence), fitted per plot and track on the
    unflagged rows (`SETTLED.md` §5: two seasonal series always correlate); flagged rows stay empty."""
    res = pd.Series(np.nan, index=d.index)
    for _, g in d.groupby(["plot", "track"]):
        m = ok.loc[g.index] & g[col].notna() & (g[dt_col].notna() if dt_col else True)
        if m.sum() < 10:
            continue
        gg = g[m]
        res.loc[gg.index] = season_residual(gg[col].to_numpy(), gg["date"], gg[dt_col].to_numpy() if dt_col else None)
    return res


def by_date(acq: pd.DataFrame, rtc: pd.DataFrame, pairs: pd.DataFrame) -> pd.DataFrame:
    keep_acq = ["track", "date", "time_utc", "acq_time_exact_utc", "platform", "plot",
                "vv_db_1x1", "vh_db_1x1", "rvi_1x1",
                "wtd_at", "wtd_mean24h", "wtd_mean72h", "wtd_mean168h", "wtd_d72h", "wtd_d168h", "censored",
                "wtd_P6_cm", "laser_surface_cm", "laser_snow",
                "air_c", "st_rh_pct", "st_vpd_kpa", "dpd_c", "wind_ms", "rain_3h_mm", "rain_24h_mm", "rain_72h_mm",
                "hours_since_rain", "imgw_dew_h", "imgw_snow_depth_cm", "om_era5_land_snow_depth",
                "om_era5_land_soil_moisture_0_to_7cm", "om_era5_land_soil_temperature_0_to_7cm",
                "flag_frozen_air", "flag_frozen_soil_om", "flag_snow_om", "flag_snow_imgw", "flag_wet_rh95",
                "flag_dew_likely", "flag_wet_any",
                "s2_ndvi_3x3", "s2_ndmi_3x3", "s2_lag_days", "uav_ndvi", "uav_ndre", "uav_thermal_c", "uav_lai",
                "uav_lag_days"]
    d = acq[keep_acq].merge(rtc, on=["track", "date", "plot"], how="left")
    c = pairs[pairs.consecutive].copy()
    c["date"] = c.t2.dt.tz_convert(None).dt.date.astype(str)
    c["prev_date"] = c.t1.dt.tz_convert(None).dt.date.astype(str)
    c = c.sort_values("revisit_days").drop_duplicates(["track", "date", "plot"])
    pc = ["prev_date", "revisit_days", "platforms", "coh_3x3_n_valid", "coh_3x3_median", "coh_3x3_mean", "coh_3x3_sd",
          "coh_3x3_iqr", "coh_1x1", "phase_n_valid_3x3", "phase_n_in_component_3x3", "dlos_hard_mm_3x3_median",
          "dlos_hard_mm_3x3_mean", "dlos_hard_mm_3x3_sd", "dlos_hard_mm_3x3_iqr", "dlos_hard_mm_1x1", "dvert_hard_mm_3x3",
          "dlos_grassland_mm_3x3_median", "dlos_grassland_mm_1x1", "consistency_3x3", "d_wtd_plot_cm",
          "d_laser_los_mm", "laser_cycle_risk", "laser_crosses_fill", "frozen_or_snow"]
    d = d.merge(c[["track", "date", "plot"] + pc], on=["track", "date", "plot"], how="left")
    d.insert(1, "year", pd.to_datetime(d.date).dt.year)
    d.insert(2, "common_period_2022_2024", d.year.between(*COMMON))
    flag = d[["flag_frozen_air", "flag_snow_om", "flag_snow_imgw"]].fillna(False).astype(bool).any(axis=1)
    d["coh_3x3_anom"] = anomaly(d, "coh_3x3_median", ~flag & ~d.frozen_or_snow.fillna(True).astype(bool),
                                "revisit_days")
    for v in ("vv_db_3x3_median", "vh_db_3x3_median"):
        d[v.replace("_median", "_anom")] = anomaly(d, v, ~flag)
    d["wtd_at_anom"] = anomaly(d, "wtd_at", ~d.censored.fillna(True).astype(bool))
    d.insert(4, "plot_group", d["plot"].map(SUPERVISOR_GROUP))
    return d.sort_values(["track", "date", "plot"]).reset_index(drop=True)


# ------------------------------------------------------------------------------ plots, UAV, daily, coverage

def grid_cells(static: xr.Dataset, px: pd.DataFrame, pad: int = 4) -> pd.DataFrame:
    """The 40 m cells around the plot transect: centre coordinates and zone (for the map of plots on the grid)."""
    r0, r1 = px.row.min() - pad, px.row.max() + pad
    c0, c1 = px.col.min() - pad, px.col.max() + pad
    rows = [{"grid_row": r, "grid_col": c, "cell_E": float(static.x[c]), "cell_N": float(static.y[r]),
             "zone": int(static.zone.values[r, c]), "boardwalk_share": float(static.boardwalk_share.values[r, c])}
            for r in range(r0, r1 + 1) for c in range(c0, c1 + 1)]
    return pd.DataFrame(rows)


def plots_table(px, mat, acq) -> pd.DataFrame:
    note = pd.read_excel(field.field_root() / field.DELIVERY / "DataSet_All_RS_LAI_merged_with_WTD_Meteo.xlsx",
                         sheet_name="Note", usecols="F:H", header=1).dropna()
    note.columns = ["plot", "lon", "lat"]
    ll = note.assign(base=note["plot"].str[:2]).groupby("base").agg(
        lon=("lon", "mean"), lat=("lat", "mean"), coordinates=("plot", lambda s: ", ".join(s)))
    shape = mat.shape
    rows = []
    for p, r, c, e, n, share in zip(px["plot"], px.row, px.col, px.E, px.N, px.shares_3x3_with):
        rows.append({"plot": p, "plot_group": SUPERVISOR_GROUP[p], "lon": ll.lon.get(p), "lat": ll.lat.get(p),
                     "coordinates_from": ll.coordinates.get(p), "E_utm33": e, "N_utm33": n, "grid_row": r,
                     "grid_col": c, "mat_px_1x1": len(pt.window_index(r, c, shape, 0, mat)[0]),
                     "mat_px_3x3": len(pt.window_index(r, c, shape, 1, mat)[0]),
                     "mat_px_5x5": len(pt.window_index(r, c, shape, 2, mat)[0]), "shares_3x3_with": share})
    d = pd.DataFrame(rows)
    st = acq.drop_duplicates("plot").set_index("plot")
    keep = ["st_dist_to_mat_m", "st_depth_in_mat_m", "st_boardwalk_share", "st_wc_herb_wetland", "st_wc_grass",
            "st_wc_water", "st_lidar_dtm_m", "st_lidar_roughness_m", "st_lidar_microrelief_sd_m",
            "st_lidar_slope_deg", "st_lidar_vegh_mean_m", "st_lidar_vegh_p90_m", "st_incidence_asc_deg",
            "st_incidence_des_deg", "st_los_u_asc", "st_los_u_des", "st_coh6_asc_median", "st_coh12_asc_median",
            "st_coh6_des_median", "st_coh12_des_median"]
    d = d.join(st[keep], on="plot")
    # a terrain model held at one level with no relief is water-flattened there: roughness is not measured
    d["lidar_dtm_flattened"] = (d.st_lidar_roughness_m < 0.01) & (d.st_lidar_dtm_m.round(3) % 1 == 0)
    tc = pd.read_csv(HUB / "08_deliverables" / "field_first" / "s1_plot_temporal_coherence.csv")
    tc = tc[tc.window == "3x3"].set_index("plot")
    d["evd_temporal_coh_3x3_first_study_uncorrected_grid"] = d["plot"].map(tc["median"])
    return d


def uav_table(dates: pd.DataFrame) -> pd.DataFrame:
    u = field.load_uav_table()
    flagged = pd.Series(0, index=u.index)
    for r, _ in field.flagged_cells(field.field_root() / field.DELIVERY / "DataSet_All_RS_LAI_merged_with_WTD_Meteo.xlsx"):
        flagged.loc[u.index[r]] += 1
    nir, red, re = u["REMX_NIR_842_mean"], u["REMX_Red_668_mean"], u["REMX_Red_Edge_717_mean"]
    v = pd.DataFrame({"date": u.Date.dt.date.astype(str), "plot": u.base_plot, "subplot": u.Plot,
                      "ndvi": (nir - red) / (nir + red), "ndre": (nir - re) / (nir + re), "nir_842": nir,
                      "red_edge_717": re, "red_668": red, "thermal_c": u["Altum_Thermal_11um_mean"],
                      "ts_minus_air_c": u["Altum_Thermal_11um_mean"] - u["Air_2m_MEAN"], "lai": u["LAI"],
                      "fapar": u["LAI_fAPAR"], "flagged_cells": flagged})
    g = v.groupby(["plot", "date"])
    num = ["ndvi", "ndre", "nir_842", "red_edge_717", "red_668", "thermal_c", "ts_minus_air_c", "lai", "fapar"]
    d = pd.concat([g[num].mean(), g[["ndvi", "ndre", "thermal_c", "lai"]].std().add_suffix("_sd_between_subplots"),
                   g["subplot"].nunique().rename("n_subplots"), g[["ndvi", "thermal_c", "lai"]].count()
                   .add_prefix("n_with_"), g["flagged_cells"].sum().rename("flagged_cells_masked")], axis=1)
    day = u.assign(date=u.Date.dt.date.astype(str), plot=u.base_plot).groupby(["plot", "date"])[
        ["WTD_MEAN", "WTD_MIN", "WTD_MAX", "Air_2m_MEAN", "RH_2m_MEAN", "VPD_Kpa_MEAN", "Rain_mm_Tot_SUM"]].first()
    d = d.join(day.add_prefix("daily_")).reset_index()
    # the nearest acquisition of each track with a consecutive pair, and its radar values
    s = dates[dates.coh_3x3_median.notna() | dates.vv_db_3x3_median.notna()].copy()
    s["_t"] = pd.to_datetime(s.date)
    for track, tag in (("ascending", "asc"), ("descending", "des")):
        st = s[s.track == track]
        rows = []
        for _, r in d.iterrows():
            q = st[st["plot"] == r["plot"]]
            if q.empty:
                rows.append({})
                continue
            lag = (q._t - pd.Timestamp(r.date)).dt.days
            i = lag.abs().idxmin()
            if abs(lag[i]) > UAV_MAX_LAG_DAYS:
                rows.append({})
                continue
            rows.append({f"s1_{tag}_date": q.date[i], f"s1_{tag}_lag_days": int(lag[i]),
                         f"s1_{tag}_coh_3x3": q.coh_3x3_median[i], f"s1_{tag}_revisit_days": q.revisit_days[i],
                         f"s1_{tag}_vv_db_3x3": q.vv_db_3x3_median[i], f"s1_{tag}_vh_db_3x3": q.vh_db_3x3_median[i],
                         f"s1_{tag}_rvi_3x3": q.rvi_3x3_median[i], f"s1_{tag}_wtd_at": q.wtd_at[i]})
        d = pd.concat([d, pd.DataFrame(rows, index=d.index)], axis=1)
    d.insert(2, "plot_group", d["plot"].map(SUPERVISOR_GROUP))
    return d.sort_values(["date", "plot"]).reset_index(drop=True)


def daily_table(site: pd.DataFrame, laser_daily: pd.DataFrame) -> pd.DataFrame:
    s = site[(site.index.year >= YEARS[0]) & (site.index.year <= YEARS[1])]
    w = s[[f"wtd_P{i}_cm" for i in range(1, 10)]].resample("1D").mean()
    w.columns = [f"wtd_P{i}_daily_cm" for i in range(1, 10)]
    m = pd.DataFrame({"rain_mm": s["rain_mm"].resample("1D").sum(min_count=1),
                      "air_c_mean": s["air_c"].resample("1D").mean(),
                      "air_c_min": s["air_c"].resample("1D").min(),
                      "snow_depth_imgw_cm": s["imgw_snow_depth_cm"].resample("1D").max()})
    d = pd.concat([w, m], axis=1)
    d.index = d.index.tz_convert(None).date.astype(str)
    lz = laser_daily.set_index("date")[["surface_cm", "surface_sd_cm", "n_ok_hours"]].add_prefix("laser_")
    d = d.join(lz)
    d.index.name = "date"
    return d.reset_index()


def coverage(dates: pd.DataFrame, pairs: pd.DataFrame, uav: pd.DataFrame) -> pd.DataFrame:
    p6 = dates[dates["plot"] == "P6"]
    pp = pairs[(pairs["plot"] == "P6") & pairs.consecutive]
    rows = []
    for (y, track), g in p6.groupby(["year", "track"]):
        q = pp[(pp.t2.dt.year == y) & (pp.track == track)]
        ok = q.d_laser_los_mm.notna() & ~q.frozen_or_snow & ~q.laser_crosses_fill
        rows.append({"year": y, "track": track, "acquisitions": len(g),
                     "platforms": "/".join(sorted(set(g.platform.dropna()) - {"?"})),
                     "with_backscatter": int(g.vv_db_3x3_median.notna().sum()),
                     "consecutive_pairs": len(q), "median_revisit_days": q.revisit_days.median(),
                     "with_wtd": int(g.wtd_at.notna().sum()),
                     "with_laser_at_overpass": int(g.laser_surface_cm.notna().sum()),
                     "pairs_with_usable_laser_change": int(ok.sum()),
                     "of_which_beyond_quarter_wavelength": int((ok & (q.laser_cycle_risk == 1)).sum()),
                     "frozen_or_snow_pairs": int(q.frozen_or_snow.sum()),
                     "uav_campaigns": int(uav[pd.to_datetime(uav.date).dt.year == y].date.nunique())})
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    static = xr.open_dataset(CUBE / "silver" / "static_40m.nc")
    mat = static.zone.values == 1
    px = plot_units()
    acq = pd.read_csv(CUBE / "gold" / "plots_acq.csv.gz", low_memory=False)
    acq = acq[in_years(acq.date)].reset_index(drop=True)
    gp = pd.read_csv(CUBE / "gold" / "plots_pair.csv.gz", low_memory=False)
    gp = gp[gp.period.isin(PERIODS) & (gp.grid == "40m")]
    laser_h = pd.read_csv(CUBE / "gold" / "truth_laser_hourly.csv.gz")
    guard_lockbox(laser_h.time_utc, Protocol.load(PROTOCOL))

    pairs = pair_table(px, mat, gp)
    pairs["laser_crosses_fill"] = np.where(pairs["plot"] == "P6", laser_fill_crossed(pairs, laser_h), False)
    pairs.loc[pairs.laser_crosses_fill, ["d_laser_los_mm", "d_laser_vertical_mm"]] = np.nan
    rtc = rtc_stats(px, mat)
    dates = by_date(acq, rtc, pairs)
    uav = uav_table(dates)
    site = pd.read_csv(CUBE / "silver" / "site_hourly.csv.gz", index_col=0, parse_dates=True, low_memory=False)
    site.index = pd.DatetimeIndex(site.index, tz="UTC") if site.index.tz is None else site.index
    daily = daily_table(site, pd.read_csv(CUBE / "gold" / "truth_laser_daily.csv"))
    plots = plots_table(px, mat, acq)
    cov = coverage(dates, pairs, uav)
    cells = grid_cells(static, px)

    for name, d in (("s1_by_date", dates), ("s1_pairs", pairs)):
        d.to_csv(OUT / f"{name}.csv.gz", index=False)
    for name, d in (("plots", plots), ("uav_by_campaign", uav), ("daily_field", daily), ("coverage", cov),
                    ("grid_cells", cells)):
        d.to_csv(OUT / f"{name}.csv", index=False)
    print(f"s1_by_date {dates.shape}, s1_pairs {pairs.shape}, uav {uav.shape}, daily {daily.shape}, "
          f"plots {plots.shape}, coverage {cov.shape} → {OUT}")


if __name__ == "__main__":
    main()

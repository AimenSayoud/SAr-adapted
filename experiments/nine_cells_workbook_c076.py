"""C-076 — every value we hold for the nine 40 m cells around P6, 2017–2026, as tables for the professor's workbook.

Descriptive extraction, no claim: per cell (attributes, the C-074 summaries), per Sentinel-1 scene (VV, VH), per
consecutive pair (line-of-sight and vertical change against the hard-target reference, coherence, backscatter change),
per UAV flight (NDVI and thermal over the part of each cell inside the strip), and the field records of the three WTD
plots inside the nine cells (P5 south, P6 centre, P7 north) and the P6 laser — at each acquisition (the professor's
24 h / 3-day / 7-day means and 3-day / 7-day changes) and daily.

Reads the cube, the C-074 tables and the C-065 UAV stack; writes only to the hub
``08_deliverables/nine_cells_c076/tables/``. No field value in this file.

    make run E=nine_cells_workbook_c076
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from pyproj import Transformer

sys.path.insert(0, str(Path(__file__).resolve().parent))
import p6_validation_x082 as X  # noqa: E402  (the X-082 readers: consecutive pairs with the references)

from insar_wetlands.inversion.isbas import PHASE_TO_MM  # noqa: E402
from insar_wetlands.paths import cube_root, hub_root  # noqa: E402

HUB, CUBE = hub_root(), cube_root()
OUT = HUB / "08_deliverables" / "nine_cells_c076" / "tables"
C074 = HUB / "08_deliverables" / "supervisor_pack_c074" / "tables"
UAV = HUB / "06_data" / "field" / "processed" / "uav_0p5m.nc"
UAV_QC = HUB / "08_deliverables" / "field_uav_c065" / "qc_flags.csv"
NAMES = {(-1, -1): "NW", (-1, 0): "N", (-1, 1): "NE", (0, -1): "W", (0, 0): "P6", (0, 1): "E",
         (1, -1): "SW", (1, 0): "S", (1, 1): "SE"}
PERIODS = ["2017_2019", "2020_2021", "2022_2024", "2025", "2026"]
PLOTS_IN = {"P5": "S", "P6": "P6", "P7": "N"}
HALF = 20.0


def cells_table(st: xr.Dataset) -> pd.DataFrame:
    pu = pd.read_csv(CUBE / "silver" / "plot_units.csv").set_index("plot")
    r6, c6 = int(pu.loc["P6", "row"]), int(pu.loc["P6", "col"])
    to_ll = Transformer.from_crs(32633, 4326, always_xy=True)
    rows = []
    for (dr, dc), name in NAMES.items():
        r, c = r6 + dr, c6 + dc
        x, y = float(st.x[c]), float(st.y[r])
        lon, lat = to_ll.transform(x, y)
        plots = [p for p, rr in pu.iterrows() if int(rr.row) == r and int(rr["col"]) == c]
        rows.append({"cell": name, "row": r, "col": c, "centre_E_utm33": x, "centre_N_utm33": y, "lat": lat, "lon": lon,
                     "wtd_plots_in_cell": ", ".join(plots), "incidence_asc_deg": float(st.incidence_asc_deg.values[r, c]),
                     "incidence_des_deg": float(st.incidence_des_deg.values[r, c])})
    t = pd.DataFrame(rows)
    a = pd.read_csv(C074 / "nine_cells.csv").drop(columns=["row", "col", "x", "y"])
    return t.merge(a, on="cell", how="left")


def scene_times(p: pd.DataFrame) -> dict:
    """Acquisition time (UTC) of each scene date per track, from the pairs' t1/t2."""
    out = {}
    for track in X.TRACKS:
        g = p[p.track == track]
        t = pd.concat([g.t1, g.t2]).drop_duplicates()
        out[track] = {ts.date(): ts for ts in t}
    return out


def field_hourly() -> pd.DataFrame:
    cols = ["time_utc", "laser_surface_cm", "laser_ok", "wtd_P5_cm", "wtd_P6_cm", "wtd_P7_cm", "censored_P5",
            "censored_P6", "censored_P7", "rain_mm", "air_c", "flag_snow_imgw", "flag_frozen_soil_om", "st_rh_pct"]
    s = pd.read_csv(CUBE / "silver" / "site_hourly.csv.gz", usecols=cols, low_memory=False)
    s.index = pd.to_datetime(s.pop("time_utc"), format="ISO8601", utc=True)
    ok = s.laser_ok.map({True: True, False: False, "True": True, "False": False}).fillna(False)
    s["laser_cm"] = s.laser_surface_cm.where(ok)
    return s


def at_acquisition(s: pd.DataFrame, t: pd.Timestamp) -> dict:
    """The professor's WTD descriptors at time t for P5/P6/P7, the laser and the weather."""
    h = t.floor("h")
    r = {}
    for p in PLOTS_IN:
        w = s[f"wtd_{p}_cm"]
        at = w.get(h, np.nan)
        m24 = w[h - pd.Timedelta(hours=23):h].mean()
        m72 = w[h - pd.Timedelta(hours=71):h].mean()
        m168 = w[h - pd.Timedelta(hours=167):h].mean()
        r[f"wtd_{p}_at_cm"] = at
        r[f"wtd_{p}_mean24h_cm"], r[f"wtd_{p}_mean3d_cm"], r[f"wtd_{p}_mean7d_cm"] = m24, m72, m168
        r[f"wtd_{p}_change3d_cm"] = at - w.get(h - pd.Timedelta(hours=72), np.nan)
        r[f"wtd_{p}_change7d_cm"] = at - w.get(h - pd.Timedelta(hours=168), np.nan)
        r[f"wtd_{p}_censored"] = bool(pd.to_numeric(pd.Series([s[f"censored_{p}"].get(h, 0)]), errors="coerce").fillna(0).iloc[0])
    r["laser_P6_cm"] = s.laser_cm.get(h, np.nan)
    r["air_c"] = s.air_c.get(h, np.nan)
    r["rh_pct"] = s.st_rh_pct.get(h, np.nan)
    r["rain_24h_mm"] = s.rain_mm[h - pd.Timedelta(hours=23):h].sum(min_count=1)
    r["rain_3d_mm"] = s.rain_mm[h - pd.Timedelta(hours=71):h].sum(min_count=1)
    r["snow_imgw"] = str(s.flag_snow_imgw.get(h, "")) == "True"
    r["frozen_soil_era5"] = str(s.flag_frozen_soil_om.get(h, "")) == "True"
    return r


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    st = xr.open_dataset(CUBE / "silver" / "static_40m.nc")
    cells = cells_table(st)
    ry = xr.DataArray(cells.row.to_numpy(), dims="cell")
    cx = xr.DataArray(cells["col"].to_numpy(), dims="cell")
    names = list(cells.cell)
    cosi = {"ascending": np.cos(np.radians(cells.incidence_asc_deg.to_numpy())),
            "descending": np.cos(np.radians(cells.incidence_des_deg.to_numpy()))}
    p = X.load_pairs()
    times = scene_times(p)
    s = field_hourly()
    las_los, lev = X.laser()

    # per scene: backscatter
    scenes, acq = [], []
    for track in X.TRACKS:
        rtc = xr.open_dataset(CUBE / "silver" / f"rtc_{track}.nc")
        vv = rtc.vv_db.isel(y=ry, x=cx).values
        vh = rtc.vh_db.isel(y=ry, x=cx).values
        for k, d in enumerate(pd.DatetimeIndex(rtc.time.values)):
            t = times[track].get(d.date(), pd.Timestamp(d.date(), tz="UTC") + pd.Timedelta(hours=16 if track == "ascending" else 5))
            for j, c in enumerate(names):
                scenes.append({"track": track, "date": d.date(), "time_utc": t, "cell": c,
                               "vv_db": float(vv[k, j]), "vh_db": float(vh[k, j])})
            acq.append({"track": track, "date": d.date(), "time_utc": t, **at_acquisition(s, t)})
    scenes = pd.DataFrame(scenes)
    acq = pd.DataFrame(acq)

    # per consecutive pair: phase change, coherence, backscatter change
    rows, pf = [], []
    sc_idx = scenes.set_index(["track", "date", "cell"])
    for track in X.TRACKS:
        g = p[p.track == track].sort_values("t1").reset_index(drop=True)
        for per in PERIODS:
            ds = xr.open_dataset(CUBE / "silver" / f"ifg_{per}_{track}.nc")
            idx = pd.Index(ds.pair.values)
            sel = g[g.pair.isin(idx)]
            if sel.empty:
                continue
            pos = idx.get_indexer(sel.pair)
            unw = ds.unw.isel(y=ry, x=cx).values[pos].astype(float)
            coh = ds.coh.isel(y=ry, x=cx).values[pos].astype(float)
            los = (unw - sel.ref_hard_unw.to_numpy(float)[:, None]) * PHASE_TO_MM
            los_g = (unw - sel.ref_grassland_unw.to_numpy(float)[:, None]) * PHASE_TO_MM
            for k, q in enumerate(sel.itertuples()):
                lt1, lt2 = las_los[track].get(q.t1.floor("h"), np.nan), las_los[track].get(q.t2.floor("h"), np.nan)
                v1, v2 = lev.get(q.t1.floor("h"), np.nan), lev.get(q.t2.floor("h"), np.nan)
                rec = {"track": track, "period_file": per, "pair": q.pair, "t1": q.t1, "t2": q.t2,
                       "revisit_days": int(q.revisit_days), "laser_P6_los_mm": lt2 - lt1, "laser_P6_vertical_mm": v2 - v1}
                for pl in PLOTS_IN:
                    w = s[f"wtd_{pl}_cm"]
                    rec[f"wtd_{pl}_change_cm"] = w.get(q.t2.floor("h"), np.nan) - w.get(q.t1.floor("h"), np.nan)
                rec["rain_between_mm"] = s.rain_mm[q.t1.floor("h"):q.t2.floor("h")].sum(min_count=1)
                rec["snow_or_frozen_t1_or_t2"] = any(str(s[c].get(t.floor("h"), "")) == "True"
                                                      for c in ("flag_snow_imgw", "flag_frozen_soil_om") for t in (q.t1, q.t2))
                pf.append(rec)
                for j, c in enumerate(names):
                    def bs(var, t):
                        key = (track, t.date(), c)
                        return float(sc_idx.loc[key, var]) if key in sc_idx.index else np.nan
                    rows.append({"track": track, "pair": q.pair, "t1": q.t1, "t2": q.t2, "revisit_days": int(q.revisit_days),
                                 "cell": c, "los_mm": los[k, j], "vertical_mm": los[k, j] / cosi[track][j],
                                 "los_mm_grassland_ref": los_g[k, j], "coherence": coh[k, j],
                                 "d_vv_db": bs("vv_db", q.t2) - bs("vv_db", q.t1), "d_vh_db": bs("vh_db", q.t2) - bs("vh_db", q.t1)})
    pairs = pd.DataFrame(rows)
    pairs_field = pd.DataFrame(pf)

    # UAV: NDVI and thermal over the part of each cell inside the strip
    u = xr.open_dataset(UAV)
    qc = pd.read_csv(UAV_QC)
    ur = []
    for cam, red, nir, dim in (("Altum", "altum_Red_668", "altum_NIR_842", "date_altum"),
                               ("RedEdge-MX Dual", "remx_Red_668", "remx_NIR_842", "date_remx")):
        for d in u[dim].values:
            sub = u.sel({dim: d})
            for cr in cells.itertuples():
                box = sub.sel(x=slice(cr.centre_E_utm33 - HALF, cr.centre_E_utm33 + HALF),
                              y=slice(cr.centre_N_utm33 + HALF, cr.centre_N_utm33 - HALF))
                rr, nn = box[red].values, box[nir].values
                ok = np.isfinite(rr) & np.isfinite(nn) & ((rr + nn) > 0)
                share = float(ok.sum() / ((2 * HALF / 0.5) ** 2))
                ndvi = float(np.nanmean((nn[ok] - rr[ok]) / (nn[ok] + rr[ok]))) if ok.any() else np.nan
                th = float(np.nanmean(box["altum_Thermal_11um"].values)) if cam == "Altum" and ok.any() else np.nan
                q = qc[(qc.sensor.str.startswith(cam.split()[0][:5])) & (qc.date == str(pd.Timestamp(d).date()))]
                ur.append({"camera": cam, "date": pd.Timestamp(d).date(), "cell": cr.cell, "cell_share_covered": share,
                           "ndvi_mean": ndvi if share > 0 else np.nan, "thermal_mean_c": th if share > 0 else np.nan,
                           "use_ndvi": bool(q.use_ndvi.iloc[0]) if len(q) else None,
                           "use_thermal": bool(q.use_thermal.iloc[0]) if len(q) and cam == "Altum" else None})
    uav = pd.DataFrame(ur)
    uav = uav[uav.cell_share_covered > 0]

    # daily field
    day = pd.DataFrame({"laser_P6_cm": s.laser_cm.resample("D").mean(),
                        **{f"wtd_{p}_cm": s[f"wtd_{p}_cm"].resample("D").mean() for p in PLOTS_IN},
                        **{f"censored_{p}_hours": pd.to_numeric(s[f"censored_{p}"], errors="coerce").resample("D").sum() for p in PLOTS_IN},
                        "rain_mm": s.rain_mm.resample("D").sum(min_count=1), "air_c": s.air_c.resample("D").mean()})
    day = day[(day.index.year >= 2020) & (day.index.year <= 2025)]
    day.index = day.index.date

    # C-074 per-cell summaries (design years)
    by = pd.read_csv(C074 / "nine_cells_by_group.csv")
    bsum = pd.read_csv(C074 / "nine_cells_backscatter.csv")

    cells.to_csv(OUT / "cells.csv", index=False)
    scenes.to_csv(OUT / "s1_scenes.csv", index=False)
    acq.to_csv(OUT / "field_at_acquisition.csv", index=False)
    pairs.to_csv(OUT / "s1_pairs.csv", index=False)
    pairs_field.to_csv(OUT / "pairs_field.csv", index=False)
    uav.to_csv(OUT / "uav.csv", index=False)
    day.rename_axis("date").to_csv(OUT / "daily_field.csv")
    by.to_csv(OUT / "c074_by_group.csv", index=False)
    bsum.to_csv(OUT / "c074_backscatter.csv", index=False)
    print(f"cells {len(cells)}, scenes {scenes.date.nunique()} dates, pairs {pairs.pair.nunique()} "
          f"({pairs.t1.min():%Y-%m-%d} → {pairs.t2.max():%Y-%m-%d}), uav rows {len(uav)} "
          f"(cells {sorted(uav.cell.unique())}), daily {len(day)}")


if __name__ == "__main__":
    main()

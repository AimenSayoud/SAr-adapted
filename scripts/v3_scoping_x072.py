"""X-072: what the data cube says before fusion v3 is modelled (exploratory, Mac-safe, reads only the cube).

Eight questions the v3 plan depends on, each answered from ``06_data/cube/gold`` and ``silver``:

1. Truth size — how many P6 pairs have a usable laser change, by period, track, revisit, season.
2. Signal budget preview — at P6, short pairs: phase = m·laser_LOS + residual; what does the residual follow
   (dew, rain, backscatter change, soil moisture, temperature, water table)? Year-blocked fits, permutation p.
3. Detectability — share of pairs whose phase matches the laser within λ/8, by revisit, coherence, surface state,
   and laser change size (the cycle limit).
4. Reference — residual spread of the P6 phase vs laser with each reference (grassland, stable ground, hard targets).
5. Spatial modes — EOFs of the mat's departures from its own median (6-day pairs): variance share and loading
   correlation with the mat's structure (depth in the mat, LiDAR, vegetation), against size-matched stable patches.
6. Laser and vegetation — does the laser's residual (after the water table) follow canopy (NDRE, NDVI) in summer?
7. Laser day–night — does the surface move between dawn (05:09) and dusk (16:36) on the same day?
8. Connections — season-removed anomalies of every acquisition-level variable at P6, one Spearman matrix.

    PYTHONPATH=src python scripts/v3_scoping_x072.py
Outputs only to the hub (field values): ``08_deliverables/fusion_v3/scoping/``.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import xarray as xr  # noqa: E402
from scipy import stats  # noqa: E402

from insar_wetlands.cube import core  # noqa: E402

HUB = Path(__file__).resolve().parents[3]
CUBE = HUB / "06_data" / "cube"
G, S = CUBE / "gold", CUBE / "silver"
OUT = HUB / "08_deliverables" / "fusion_v3" / "scoping"
LAM = core.WAVELENGTH_MM
RNG = np.random.default_rng(72)
SEASON = {12: "DJF", 1: "DJF", 2: "DJF", 3: "MAM", 4: "MAM", 5: "MAM", 6: "JJA", 7: "JJA", 8: "JJA",
          9: "SON", 10: "SON", 11: "SON"}


def truthy(x: pd.Series) -> pd.Series:
    return x.map({True: True, False: False, "True": True, "False": False}).fillna(False).astype(bool)


def p6_pairs() -> pd.DataFrame:
    d = pd.read_csv(G / "plots_pair.csv.gz", low_memory=False, parse_dates=["t1", "t2"])
    p = d[(d["plot"] == "P6") & (d.grid == "40m")].copy()
    p["frozen"] = truthy(p.flag_frozen_air_t1) | truthy(p.flag_frozen_air_t2)
    p["snow"] = truthy(p.laser_snow_t1) | truthy(p.laser_snow_t2) | truthy(p.flag_snow_om_t1) | truthy(p.flag_snow_om_t2)
    p["wet1"], p["wet2"] = truthy(p.flag_wet_any_t1), truthy(p.flag_wet_any_t2)
    p["state"] = np.where(p.wet1 | p.wet2, "wet", "dry")
    p["usable"] = p.d_laser_los_mm.notna() & p.unw_1x1.notna() & ~p.frozen & ~p.snow
    mid = p.t1 + (p.t2 - p.t1) / 2
    p["season"] = mid.dt.month.map(SEASON)
    p["doy"] = mid.dt.dayofyear
    for ref in ("grassland", "stable", "hard"):
        c = f"ref_{ref}_unw"
        if c in p:
            p[f"phase_mm_{ref}"] = -core.rad_to_los_mm(p.unw_1x1 - p[c])
    return p


# ------------------------------------------------------------------------------ 1 truth size

def q1_truth(p: pd.DataFrame) -> pd.DataFrame:
    u = p[p.usable & (p.revisit_days <= 24)]
    t = u.groupby(["period", "track", "revisit_days", "season"]).size().unstack("season", fill_value=0)
    t["total"] = t.sum(1)
    return t.reset_index()


# ------------------------------------------------------------------------------ 2 signal budget

DRIVERS = {"dpd_c_t1": "dew-point depression at t1 (°C)", "dpd_c_t2": "dew-point depression at t2 (°C)",
           "wet_both": "number of wet ends (0–2)", "rain_sum_mm": "rain between the dates (mm)",
           "d_vv_db_3x3": "VV change (dB)", "d_vh_db_3x3": "VH change (dB)",
           "d_om_era5_land_soil_moisture_0_to_7cm": "ERA5-Land top-soil moisture change",
           "d_air_c": "air temperature change (°C)", "d_wtd_plot_cm": "P6 water-table change (cm)",
           "et0_sum_mm": "reference ET between the dates (mm)"}


def perm_r(x, y, n=4000):
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if len(x) < 10 or np.ptp(x) == 0 or np.ptp(y) == 0:      # a constant column carries no information
        return np.nan, np.nan, len(x)
    r = stats.spearmanr(x, y)[0]
    null = np.array([stats.spearmanr(x, RNG.permutation(y))[0] for _ in range(n)])
    return float(r), float((1 + np.sum(np.abs(null) >= abs(r))) / (1 + n)), int(len(x))


def q2_budget(p: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    u = p[p.usable & p.revisit_days.isin([6, 12])].copy()
    u["wet_both"] = u.wet1.astype(int) + u.wet2.astype(int)
    fits, rows = [], []
    for (tr, rv), g in u.groupby(["track", "revisit_days"]):
        if len(g) < 12:
            continue
        x, y = g.d_laser_los_mm.to_numpy(), g.phase_mm_stable.to_numpy()
        ok = np.isfinite(x) & np.isfinite(y)
        # year-blocked: m fitted without the pair's own year, residual from that fit
        res = np.full(len(g), np.nan)
        yrs = g.t1.dt.year.to_numpy()
        for yv in np.unique(yrs):
            tr_ = ok & (yrs != yv)
            te = ok & (yrs == yv)
            if tr_.sum() >= 6 and te.any():
                m = np.sum(x[tr_] * y[tr_]) / np.sum(x[tr_] ** 2)
                res[te] = y[te] - m * x[te]
        m_all = np.sum(x[ok] * y[ok]) / np.sum(x[ok] ** 2)
        var_y = np.nanvar(y[ok])
        fits.append({"track": tr, "revisit_days": rv, "n": int(ok.sum()), "m_all_years": m_all,
                     "r_phase_laser": float(np.corrcoef(x[ok], y[ok])[0, 1]),
                     "share_var_motion": float(1 - np.nanvar(y[ok] - m_all * x[ok]) / var_y),
                     "share_var_motion_year_blocked": float(1 - np.nanvar(res[np.isfinite(res)]) / var_y)})
        for k, lab in DRIVERS.items():
            if k not in g:
                continue
            r, pv, n = perm_r(g[k].to_numpy(float), res)
            rows.append({"track": tr, "revisit_days": rv, "driver": k, "label": lab, "n": n, "spearman_r": r, "perm_p": pv})
    return pd.DataFrame(fits), pd.DataFrame(rows)


# ------------------------------------------------------------------------------ 3 detectability

def q3_detect(p: pd.DataFrame) -> pd.DataFrame:
    u = p[p.usable & (p.revisit_days <= 24)].copy()
    u["match"] = (u.phase_mm_stable - u.d_laser_los_mm).abs() < LAM / 8
    u["coh_class"] = pd.cut(u.coh_3x3, [0, 0.3, 0.45, 0.6, 1.0])
    u["motion_class"] = pd.cut(u.d_laser_los_mm.abs(), [0, LAM / 8, LAM / 4, 1e3], labels=["< λ/8", "λ/8–λ/4", "> λ/4"])
    out = []
    for col in ("revisit_days", "coh_class", "state", "motion_class", "track", "season"):
        g = u.groupby(col, observed=True).match.agg(["mean", "size"]).reset_index().rename(columns={col: "level"})
        g.insert(0, "factor", col)
        g["level"] = g["level"].astype(str)
        out.append(g)
    return pd.concat(out, ignore_index=True).rename(columns={"mean": "share_matching", "size": "n"})


# ------------------------------------------------------------------------------ 4 reference

def q4_reference(p: pd.DataFrame) -> pd.DataFrame:
    u = p[p.usable & p.revisit_days.isin([6, 12])]
    rows = []
    for (tr, rv), g in u.groupby(["track", "revisit_days"]):
        for ref in ("grassland", "stable", "hard"):
            c = f"phase_mm_{ref}"
            if c not in g or g[c].notna().sum() < 10:
                continue
            ok = g[c].notna()
            x, y = g.d_laser_los_mm[ok], g[c][ok]
            m = np.sum(x * y) / np.sum(x ** 2)
            rows.append({"track": tr, "revisit_days": rv, "reference": ref, "n": int(ok.sum()), "m": m,
                         "r": float(np.corrcoef(x, y)[0, 1]), "residual_sd_mm": float(np.std(y - m * x))})
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------------ 5 spatial modes

def departures(ds: xr.Dataset, sel: np.ndarray, revisit: int, stable_ref: np.ndarray) -> np.ndarray:
    """(pairs × pixels) of a zone's unwrapped phase minus the zone's own median per pair (mm), frozen pairs kept."""
    keep = ds.revisit_days.values == revisit
    U = ds.unw.values[:, keep].T                              # pairs × pixels
    U = U - np.nanmedian(U[:, stable_ref], 1, keepdims=True)
    Z = U[:, sel]
    Z = Z - np.nanmedian(Z, 1, keepdims=True)
    good = np.isfinite(Z).mean(0) > 0.8
    Z = Z[:, good]
    Z = np.where(np.isfinite(Z), Z, 0.0)
    return -core.rad_to_los_mm(Z), good


def eof(Z: np.ndarray, k: int = 3):
    Zc = Z - Z.mean(0)
    u, s, vt = np.linalg.svd(Zc, full_matrices=False)
    var = s ** 2 / np.sum(s ** 2)
    return var[:k], vt[:k]


def q5_spatial(n_null: int = 200) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows, corr = [], []
    for period, rv in (("2020_2021", 6), ("2025", 6), ("2022_2024", 12)):
        for track in ("ascending", "descending"):
            f = G / f"mat_pairs_{period}_{track}.nc"
            if not f.exists():
                continue
            ds = xr.open_dataset(f)
            zone, stable = ds.st_zone.values, ds.st_stable.values == 1
            mat = zone == 1
            Z, good = departures(ds, mat, rv, stable)
            if Z.shape[0] < 10:
                continue
            var, vt = eof(Z)
            # size-matched null: compact stable patches with the same pixel count
            rr, cc = ds.row.values, ds.col.values
            st_idx = np.flatnonzero(stable)
            nulls = []
            n_px = int(good.sum())
            for _ in range(n_null):
                c0 = st_idx[RNG.integers(len(st_idx))]
                d2 = (rr[st_idx] - rr[c0]) ** 2 + (cc[st_idx] - cc[c0]) ** 2
                patch = np.zeros(len(zone), bool)
                patch[st_idx[np.argsort(d2)[:n_px]]] = True
                Zs, _ = departures(ds, patch, rv, stable & ~patch)
                if Zs.shape[1] >= n_px // 2:
                    nulls.append(eof(Zs, 1)[0][0])
            nulls = np.array(nulls)
            rows.append({"period": period, "track": track, "revisit_days": rv, "pairs": Z.shape[0], "mat_pixels": n_px,
                         "eof1_share": var[0], "eof2_share": var[1], "eof3_share": var[2],
                         "null_eof1_median": float(np.median(nulls)) if len(nulls) else np.nan,
                         "null_eof1_p95": float(np.percentile(nulls, 95)) if len(nulls) else np.nan,
                         "p_eof1": float((1 + np.sum(nulls >= var[0])) / (1 + len(nulls))) if len(nulls) else np.nan})
            pix = {v: ds[f"st_{v}"].values[mat][good] for v in ("depth_in_mat_m", "lidar_vegh_mean_m",
                                                                 "lidar_microrelief_sd_m", "lidar_dtm_m", "boardwalk_share")
                   if f"st_{v}" in ds}
            for k in range(3):
                for v, a in pix.items():
                    ok = np.isfinite(a)
                    if ok.sum() > 20:
                        corr.append({"period": period, "track": track, "eof": k + 1, "pixel_layer": v,
                                     "abs_spearman": abs(stats.spearmanr(vt[k][ok], a[ok])[0])})
    return pd.DataFrame(rows), pd.DataFrame(corr)


# ------------------------------------------------------------------------------ 6, 7 laser checks

def harmonic_resid(t: pd.Series, v: pd.Series) -> np.ndarray:
    """Residual of a two-harmonic annual fit (season removed)."""
    doy = t.dt.dayofyear.to_numpy() / 365.25 * 2 * np.pi
    X = np.column_stack([np.ones(len(t)), np.sin(doy), np.cos(doy), np.sin(2 * doy), np.cos(2 * doy)])
    y = v.to_numpy(float)
    ok = np.isfinite(y)
    out = np.full(len(y), np.nan)
    if ok.sum() > 10:
        b = np.linalg.lstsq(X[ok], y[ok], rcond=None)[0]
        out[ok] = y[ok] - X[ok] @ b
    return out


def q6_laser_vegetation() -> pd.DataFrame:
    a = pd.read_csv(G / "plots_acq.csv.gz", low_memory=False, parse_dates=["time_utc"])
    p = a[(a["plot"] == "P6") & a.laser_surface_cm.notna()].drop_duplicates("date").copy()
    rows = []
    # laser residual after the water table (linear, per year to absorb laser re-levelling)
    p["laser_resid_cm"] = np.nan
    for y, g in p.groupby(p.time_utc.dt.year):
        ok = g.wtd_at.notna()
        if ok.sum() > 8:
            b = np.polyfit(g.wtd_at[ok], g.laser_surface_cm[ok], 1)
            p.loc[g.index[ok], "laser_resid_cm"] = g.laser_surface_cm[ok] - np.polyval(b, g.wtd_at[ok])
    summer = p.time_utc.dt.month.isin([6, 7, 8, 9])
    # raw relations are confounded by the season (canopy and mat both have an annual cycle, SETTLED §5): the test that
    # counts is on season-removed anomalies of both series
    p["laser_resid_anom"] = harmonic_resid(p.time_utc, p.laser_resid_cm)
    for v in ("s2_ndre_3x3", "s2_ndvi_3x3", "s2_ndmi_3x3", "lst_c_3x3", "air_c", "uav_lai", "uav_ndre"):
        if v not in p:
            continue
        anom = harmonic_resid(p.time_utc, p[v]) if p[v].notna().sum() > 20 else np.full(len(p), np.nan)
        for lab, m, x_all, y_all in (("raw, all", np.ones(len(p), bool), p[v].to_numpy(float), p.laser_resid_cm.to_numpy(float)),
                                     ("raw, summer", summer.to_numpy(), p[v].to_numpy(float), p.laser_resid_cm.to_numpy(float)),
                                     ("season removed", np.ones(len(p), bool), anom, p.laser_resid_anom.to_numpy(float))):
            r, pv, n = perm_r(x_all[m], y_all[m], 2000)
            rows.append({"variable": v, "test": lab, "n": n, "spearman_r_with_laser_residual": r, "perm_p": pv})
    return pd.DataFrame(rows)


def q7_laser_daynight() -> pd.DataFrame:
    h = pd.read_csv(G / "truth_laser_hourly.csv.gz", index_col=0, parse_dates=True)
    h.index = pd.DatetimeIndex(h.index, tz="UTC") if h.index.tz is None else h.index
    ok = h.laser_surface_cm.where(truthy(h.laser_ok))
    days = pd.DataFrame({"dawn": ok.at_time("05:00"), })
    dawn = ok[ok.index.hour == 5]
    dusk = ok[ok.index.hour == 17]
    d = pd.DataFrame({"dawn": dawn.groupby(dawn.index.date).mean(), "dusk": dusk.groupby(dusk.index.date).mean()}).dropna()
    d["dusk_minus_dawn_mm"] = (d.dusk - d.dawn) * 10
    d["month"] = pd.to_datetime(d.index).month
    del days
    s = d.groupby("month").dusk_minus_dawn_mm.agg(["median", "mean", "std", "size"]).reset_index()
    s.loc[len(s)] = ["all", d.dusk_minus_dawn_mm.median(), d.dusk_minus_dawn_mm.mean(), d.dusk_minus_dawn_mm.std(), len(d)]
    return s


# ------------------------------------------------------------------------------ 8 connections

CONN = {"vv_db_3x3": "VV", "vh_db_3x3": "VH", "rvi_3x3": "RVI", "s2_ndvi_3x3": "NDVI", "s2_ndre_3x3": "NDRE",
        "s2_ndmi_3x3": "NDMI", "s2_ndwi_3x3": "NDWI", "lst_c_3x3": "LST", "air_c": "air T", "dpd_c": "dew-pt dep.",
        "om_era5_land_soil_moisture_0_to_7cm": "soil moist.", "om_era5_land_soil_temperature_0_to_7cm": "soil T",
        "rain_72h_mm": "rain 72 h", "wtd_at": "water table", "wtd_mean168h": "WT 7 d", "laser_surface_cm": "laser",
        "uav_ndre": "UAV NDRE", "uav_lai": "UAV LAI"}


def q8_connections() -> pd.DataFrame:
    a = pd.read_csv(G / "plots_acq.csv.gz", low_memory=False, parse_dates=["time_utc"])
    p = a[(a["plot"] == "P6") & (a.track == "ascending")].copy()
    an = pd.DataFrame({lab: harmonic_resid(p.time_utc, p[c]) for c, lab in CONN.items() if c in p})
    # 6-day coherence at P6 on the pair ending at each date (the radar's own view of change)
    pp = pd.read_csv(G / "plots_pair.csv.gz", low_memory=False, parse_dates=["t1", "t2"])
    c6 = pp[(pp["plot"] == "P6") & (pp.grid == "40m") & (pp.track == "ascending") & (pp.revisit_days.isin([6, 12]))]
    c6 = c6.sort_values("revisit_days").drop_duplicates("t2").set_index(c6.t2.dt.normalize().name if False else "t2")
    key = p.time_utc.dt.normalize()
    cmap = c6.coh_3x3.copy()
    cmap.index = pd.DatetimeIndex(cmap.index).normalize()
    an["coherence (6/12 d)"] = harmonic_resid(p.time_utc, pd.Series(key.map(cmap).to_numpy(), index=p.index))
    return an.corr(method="spearman", min_periods=25)


# ------------------------------------------------------------------------------ report

def md(df: pd.DataFrame, nd: int = 2) -> str:
    d = df.copy()
    for c in d.columns:
        if d[c].dtype.kind == "f":
            d[c] = d[c].round(nd)
    cols = [str(c) for c in d.columns]
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for r in d.itertuples(index=False):
        lines.append("| " + " | ".join("" if (isinstance(v, float) and np.isnan(v)) else str(v) for v in r) + " |")
    return "\n".join(lines)


def top_pairs(c: pd.DataFrame, thr: float = 0.3) -> pd.DataFrame:
    rows = []
    cols = list(c.columns)
    for i in range(len(cols)):
        for j in range(i + 1, len(cols)):
            v = c.iloc[i, j]
            if np.isfinite(v) and abs(v) >= thr:
                rows.append({"a": cols[i], "b": cols[j], "spearman_r": v})
    return pd.DataFrame(rows).sort_values("spearman_r", key=np.abs, ascending=False) if rows else pd.DataFrame()


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    p = p6_pairs()
    t1 = q1_truth(p)
    fits, budget = q2_budget(p)
    det = q3_detect(p)
    ref = q4_reference(p)
    sp, spc = q5_spatial()
    veg = q6_laser_vegetation()
    dn = q7_laser_daynight()
    conn = q8_connections()
    for name, df in (("q1_truth_size", t1), ("q2_motion_share", fits), ("q2_residual_drivers", budget),
                     ("q3_detectability", det), ("q4_reference", ref), ("q5_spatial_eof", sp),
                     ("q5_eof_vs_pixel_layers", spc), ("q6_laser_vs_vegetation", veg), ("q7_laser_dusk_minus_dawn", dn)):
        df.to_csv(OUT / f"{name}.csv", index=False)
    conn.to_csv(OUT / "q8_connections_spearman.csv")
    fig, ax = plt.subplots(figsize=(10, 8.5))
    im = ax.imshow(conn.values, cmap="RdBu_r", vmin=-0.8, vmax=0.8)
    ax.set_xticks(range(len(conn)), conn.columns, rotation=60, ha="right", fontsize=8)
    ax.set_yticks(range(len(conn)), conn.index, fontsize=8)
    for i in range(len(conn)):
        for j in range(len(conn)):
            v = conn.values[i, j]
            if np.isfinite(v) and i != j and abs(v) >= 0.3:
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=6)
    fig.colorbar(im, ax=ax, shrink=0.7, label="Spearman r of season-removed anomalies")
    ax.set_title("P6, ascending acquisitions 2020–2026: how the variables move together once the seasons are removed", fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / "fig_q8_connections.png", dpi=150)
    sig = budget[budget.perm_p < 0.05]
    det_r = det[det.factor == "revisit_days"]
    txt = ["# X-072 — what the data cube says before fusion v3 is modelled (exploratory)", "",
           "Generated by `05_code/SAr-adapted/scripts/v3_scoping_x072.py` from `06_data/cube` (every number below is "
           "computed). P6 phase in mm of line of sight, −φ·λ/4π, uplift positive; usable = laser change present, not "
           "frozen, no snow at either end.", "",
           "## 1. How much truth there is (usable P6 pairs ≤ 24 days)", "", md(t1), "",
           "## 2. Signal budget preview (P6, stable-ground reference)", "",
           "Share of the phase variance explained by the laser motion, all years and year-blocked (m fitted without the "
           "pair's year):", "", md(fits), "",
           "What the year-blocked residual follows (Spearman, permutation p, floor 1/4001):", "", md(budget), "",
           f"{len(sig)} of {len(budget)} driver tests at p < 0.05 (≈ {0.05 * len(budget):.1f} expected by chance).", "",
           "## 3. Detectability (phase within λ/8 of the laser change)", "", md(det), "",
           "## 4. Reference", "", md(ref), "",
           "## 5. Spatial modes of the mat (departures from its own median)", "",
           "EOF variance shares against compact stable-ground patches of the same pixel count (size-matched null):", "",
           md(sp), "", "Loading correlation with the mat's structure (|Spearman|):", "",
           md(spc.pivot_table(index=["period", "track", "eof"], columns="pixel_layer", values="abs_spearman").reset_index())
           if len(spc) else "(none)", "",
           "## 6. Laser and vegetation (laser residual after the water table, per year)", "",
           "Raw relations share the annual cycle; only the season-removed row tests whether the canopy enters the laser.", "",
           md(veg), "",
           "## 7. Does the surface move between dawn and dusk? (laser 17 h − 05 h UTC, same day, QC-passed)", "", md(dn), "",
           "## 8. Connections", "", "`fig_q8_connections.png`, `q8_connections_spearman.csv`: season-removed anomalies "
           "(two-harmonic annual fit removed) of every acquisition-level variable at P6, ascending.", "",
           "Strongest pairs (|r| ≥ 0.3, at least 25 dates):", "", md(top_pairs(conn)), "",
           "Note: IMGW Poznań records dew duration on almost no day of this record, so it is left out of the tests.", ""]
    _ = det_r
    (OUT / "README.md").write_text("\n".join(txt))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()

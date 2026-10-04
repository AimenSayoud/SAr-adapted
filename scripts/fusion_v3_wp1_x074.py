"""X-074 — fusion v3 WP1: is the laser a clean truth, and what lifts the mat? (exploratory, Mac-safe, reads the cube)

E1.1 Laser and canopy. Daily laser level at P6 (QC-passed hours) against the P6 water table, with and without the
     canopy (Sentinel-2 NDRE at P6) and the season; levels demeaned per year (the laser is re-levelled between years),
     slopes fitted on the other years (year-blocked). If the canopy adds beyond the water table and the season, the
     laser carries vegetation and the truth needs a canopy term.
E1.2 The driver of the mat's height. Six-day changes of the laser (every day, and on the radar's 6-day pairs) against
     the water-table change, with rising/falling asymmetry (hysteresis), soil-temperature change (gas in the peat),
     canopy change, a seasonal buoyancy; year-blocked.
E1.3 A water table where there is no logger (2025). A bucket model driven by ERA5 rain and temperature (Hamon PET),
     and a backscatter proxy (VV at P6), trained on 2020–2023, tested on 2024 (logger present), then used for 2025,
     where the laser — but no water table — exists: does the driver predict the laser's 6-day changes out of sample?

    PYTHONPATH=src python scripts/fusion_v3_wp1_x074.py
Outputs only to the hub (field values): ``08_deliverables/fusion_v3/wp1/``.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from insar_wetlands import fusion as fu  # noqa: E402

HUB = Path(__file__).resolve().parents[3]
G, S = HUB / "06_data" / "cube" / "gold", HUB / "06_data" / "cube" / "silver"
OUT = HUB / "08_deliverables" / "fusion_v3" / "wp1"
LAT = 52.76


def truthy(x: pd.Series) -> pd.Series:
    return x.map({True: True, False: False, "True": True, "False": False}).fillna(False).astype(bool)


def daily() -> pd.DataFrame:
    """One row per day at P6: laser (QC-passed mean), water table (uncensored), ERA5 / station weather, soil
    temperature, canopy (Sentinel-2 NDRE nearest within ±10 d), radar VV at P6 (both tracks, daily mean)."""
    h = pd.read_csv(S / "site_hourly.csv.gz", index_col=0, parse_dates=True, low_memory=False)
    h.index = pd.DatetimeIndex(h.index, tz="UTC") if h.index.tz is None else h.index
    las = h.laser_surface_cm.where(truthy(h.laser_ok))
    w = h.wtd_P6_cm.where(~truthy(h.censored_P6))
    d = pd.DataFrame({"laser_cm": las.resample("1D").mean(), "laser_n": las.resample("1D").count(),
                      "wtd_cm": w.resample("1D").mean(),
                      "rain_era5_mm": h.om_era5_precipitation.resample("1D").sum(min_count=20),
                      "air_era5_c": h.om_era5_land_temperature_2m.resample("1D").mean(),
                      "soil_t07_c": h.om_era5_land_soil_temperature_0_to_7cm.resample("1D").mean(),
                      "soil_t28_c": h.om_era5_land_soil_temperature_7_to_28cm.resample("1D").mean(),
                      "frozen": truthy(h.flag_frozen_soil_om).resample("1D").mean()})
    d.loc[d.laser_n < 12, "laser_cm"] = np.nan
    d.index = d.index.tz_convert(None).normalize()
    s2 = pd.read_csv(S / "s2_plots.csv.gz", parse_dates=["time_utc"])
    s2 = s2[(s2["plot"] == "P6") & (s2.valid_3x3 >= 0.5)].copy()
    s2["day"] = s2.time_utc.dt.tz_convert(None).dt.normalize()
    nd = s2.groupby("day").ndre_3x3.mean()
    d["ndre"] = nd.reindex(d.index).interpolate(limit=10, limit_direction="both")
    a = pd.read_csv(G / "plots_acq.csv.gz", low_memory=False, parse_dates=["time_utc"])
    a = a[a["plot"] == "P6"].copy()
    a["day"] = a.time_utc.dt.tz_convert(None).dt.normalize()
    vv = a.groupby(["day", "track"]).vv_db_3x3.mean().unstack()
    for t in vv:
        d[f"vv_{t[:3]}"] = vv[t].reindex(d.index)
    d["year"] = d.index.year
    d["doy"] = d.index.dayofyear
    return d


def harm(doy) -> np.ndarray:
    a = np.asarray(doy, float) / 365.25 * 2 * np.pi
    return np.column_stack([np.sin(a), np.cos(a), np.sin(2 * a), np.cos(2 * a)])


def yb_fit(X: np.ndarray, y: np.ndarray, years: np.ndarray, demean_by_year: bool) -> np.ndarray:
    """Year-blocked prediction: OLS on the other years; with ``demean_by_year`` both sides are centred per year."""
    X, y = X.copy(), y.copy()
    ok = np.isfinite(y) & np.all(np.isfinite(X), 1)
    if demean_by_year:
        for yv in np.unique(years):
            m = ok & (years == yv)
            if m.any():
                X[m] -= X[m].mean(0)
                y[m] -= y[m].mean()
    pred = np.full(len(y), np.nan)
    for yv in np.unique(years[ok]):
        tr, te = ok & (years != yv), ok & (years == yv)
        if tr.sum() > X.shape[1] + 5 and te.any():
            b = np.linalg.lstsq(X[tr], y[tr], rcond=None)[0]
            pred[te] = X[te] @ b
    return pred, y, ok


def score(pred, y) -> dict:
    ok = np.isfinite(pred) & np.isfinite(y)
    if ok.sum() < 10:
        return {"n": int(ok.sum()), "rmse_mm": np.nan, "r": np.nan}
    return {"n": int(ok.sum()), "rmse_mm": float(np.sqrt(np.mean((pred[ok] - y[ok]) ** 2)) * 10),
            "r": float(np.corrcoef(pred[ok], y[ok])[0, 1])}


# ------------------------------------------------------------------------------ E1.1

def e11(d: pd.DataFrame) -> pd.DataFrame:
    yrs = d.year.to_numpy()
    nd_anom = d.ndre.to_numpy() - harm(d.doy) @ np.linalg.lstsq(
        np.column_stack([np.ones(len(d)), harm(d.doy)])[d.ndre.notna()][:, 1:], d.ndre.dropna().to_numpy() - d.ndre.mean(),
        rcond=None)[0] - d.ndre.mean()
    models = {"water table": [d.wtd_cm.to_numpy()[:, None]],
              "water table + season": [d.wtd_cm.to_numpy()[:, None], harm(d.doy)],
              "water table + canopy (NDRE)": [d.wtd_cm.to_numpy()[:, None], d.ndre.to_numpy()[:, None]],
              "water table + season + canopy anomaly": [d.wtd_cm.to_numpy()[:, None], harm(d.doy), nd_anom[:, None]],
              "water table + soil temperature": [d.wtd_cm.to_numpy()[:, None], d.soil_t07_c.to_numpy()[:, None]],
              "water table + season + soil temperature": [d.wtd_cm.to_numpy()[:, None], harm(d.doy), d.soil_t07_c.to_numpy()[:, None]],
              "season only": [harm(d.doy)]}
    rows = []
    for name, parts in models.items():
        X = np.column_stack(parts)
        pred, y, ok = yb_fit(X, d.laser_cm.to_numpy(), yrs, demean_by_year=True)
        rows.append({"model (daily laser level, demeaned per year)": name, **score(pred, y)})
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------------ E1.2

def changes(d: pd.DataFrame, lag: int = 6) -> pd.DataFrame:
    c = pd.DataFrame(index=d.index)
    for v in ("laser_cm", "wtd_cm", "soil_t07_c", "soil_t28_c", "ndre", "air_era5_c"):
        c[f"d_{v}"] = d[v].shift(-lag) - d[v]
    c["frozen"] = d.frozen.rolling(lag + 1).max().shift(-lag)
    c["year"], c["doy"] = d.year, d.doy
    c["rain"] = d.rain_era5_mm.rolling(lag).sum().shift(-lag)
    return c[(c.frozen.fillna(1) == 0)]


def e12(d: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    c = changes(d)
    dw = c.d_wtd_cm.to_numpy()
    up, dn = np.clip(dw, 0, None), np.clip(dw, None, 0)
    hm = harm(c.doy)
    models = {"g·ΔW": [dw[:, None]], "g⁺ΔW⁺ + g⁻ΔW⁻ (hysteresis)": [up[:, None], dn[:, None]],
              "g·ΔW + k·ΔT_soil": [dw[:, None], c.d_soil_t07_c.to_numpy()[:, None]],
              "g·ΔW + c·ΔNDRE": [dw[:, None], c.d_ndre.to_numpy()[:, None]],
              "seasonal g: ΔW·(1, season)": [dw[:, None], dw[:, None] * hm[:, :2]],
              "hysteresis + ΔT_soil + ΔNDRE": [up[:, None], dn[:, None], c.d_soil_t07_c.to_numpy()[:, None],
                                               c.d_ndre.to_numpy()[:, None]],
              "ΔT_soil only (no water table)": [c.d_soil_t07_c.to_numpy()[:, None]]}
    rows, coefs = [], {}
    for name, parts in models.items():
        X = np.column_stack(parts)
        pred, y, ok = yb_fit(X, c.d_laser_cm.to_numpy(), c.year.to_numpy(), demean_by_year=False)
        okk = ok & np.isfinite(pred)
        rows.append({"model (6-day laser change)": name, **score(pred, y),
                     "rmse_zero_change_mm": float(np.sqrt(np.mean(y[okk] ** 2)) * 10) if okk.any() else np.nan})
        b = np.linalg.lstsq(X[ok], y[ok], rcond=None)[0]
        coefs[name] = b
    return pd.DataFrame(rows), coefs


# ------------------------------------------------------------------------------ E1.3

def e13(d: pd.DataFrame, g_mm_per_cm: float) -> tuple[pd.DataFrame, pd.DataFrame, pd.Series]:
    pet = fu.hamon_pet(d.air_era5_c.to_numpy(), d.doy.to_numpy(), LAT)
    rain = d.rain_era5_mm.fillna(0).to_numpy()
    w = d.wtd_cm.to_numpy()
    train = (d.year <= 2023).to_numpy()
    par = fu.fit_bucket(np.where(train, w, np.nan), rain, pet)
    floor = np.nanpercentile(w[train], 1)
    sim = pd.Series(np.nan, index=d.index)
    for start_year in (2024, 2025, 2026):
        m = (d.year >= start_year).to_numpy() & (d.year <= start_year).to_numpy()
        if not m.any():
            continue
        i0 = np.flatnonzero(m)[0]
        prev = w[:i0][np.isfinite(w[:i0])]
        w0 = prev[-1] if len(prev) else par["w0"]
        if start_year == 2026:                 # no logger since 2024: continue the 2025 simulation
            w0 = sim.iloc[i0 - 1] if np.isfinite(sim.iloc[i0 - 1]) else par["w0"]
        sim.iloc[np.flatnonzero(m)] = fu.simulate_bucket(par, w0, rain[m], pet[m], floor=floor)
    # in-sample simulation (diagnostic) for the training years
    sim_tr = fu.simulate_bucket(par, w[np.isfinite(w)][0], rain, pet, floor=floor)
    # backscatter proxy: W ~ VV (ascending and descending at P6) + season, trained ≤ 2023
    vv = d[["vv_asc", "vv_des"]].mean(1)
    X = np.column_stack([np.ones(len(d)), vv.to_numpy(), harm(d.doy)])
    okv = train & np.isfinite(w) & np.isfinite(vv.to_numpy())
    b = np.linalg.lstsq(X[okv], w[okv], rcond=None)[0]
    proxy_vv = pd.Series(X @ b, index=d.index).interpolate(limit=6, limit_direction="both")
    rows = []
    test = (d.year == 2024).to_numpy()
    for name, est in (("bucket (ERA5)", sim.to_numpy()), ("VV proxy", proxy_vv.to_numpy()),
                      ("bucket + VV (mean)", np.nanmean(np.column_stack([sim.to_numpy(), proxy_vv.to_numpy()]), 1))):
        ok = test & np.isfinite(est) & np.isfinite(w)
        rows.append({"proxy": name, "test year": 2024, "n_days": int(ok.sum()),
                     "rmse_cm": float(np.sqrt(np.mean((est[ok] - w[ok]) ** 2))) if ok.any() else np.nan,
                     "r": float(np.corrcoef(est[ok], w[ok])[0, 1]) if ok.sum() > 5 else np.nan,
                     "nse": fu.nse(w[ok], est[ok]) if ok.sum() > 5 else np.nan})
    rows.append({"proxy": "bucket (ERA5), in-sample 2020–2023", "test year": "train", "n_days": int((train & np.isfinite(w)).sum()),
                 "rmse_cm": float(np.sqrt(np.nanmean((sim_tr[train] - w[train]) ** 2))),
                 "r": float(pd.Series(sim_tr[train]).corr(pd.Series(w[train]))), "nse": fu.nse(w[train], sim_tr[train])})
    # 2025: the driver from each proxy against the laser's 6-day changes (out of sample: g from 2021–2024)
    out = []
    lag = 6
    dl = (d.laser_cm.shift(-lag) - d.laser_cm) * 10
    fr = d.frozen.rolling(lag + 1).max().shift(-lag).fillna(1) == 0
    for name, est in (("no driver (Δh = 0)", None), ("logger water table", d.wtd_cm), ("bucket (ERA5)", sim),
                      ("VV proxy", proxy_vv), ("bucket + VV (mean)", pd.concat([sim, proxy_vv], axis=1).mean(1))):
        for yr in (2024, 2025):
            m = (d.year == yr) & fr & dl.notna()
            if est is None:
                pred = pd.Series(0.0, index=d.index)
            else:
                pred = (est.shift(-lag) - est) * g_mm_per_cm
            ok = m & pred.notna()
            if ok.sum() < 10:
                continue
            p, y = pred[ok].to_numpy(), dl[ok].to_numpy()
            out.append({"driver": name, "year": yr, "n_days": int(ok.sum()), "rmse_mm": float(np.sqrt(np.mean((p - y) ** 2))),
                        "r": float(np.corrcoef(p, y)[0, 1]) if np.std(p) > 0 else np.nan,
                        "amplitude_ratio": float(np.std(p) / np.std(y))})
    return pd.DataFrame(rows), pd.DataFrame(out), sim


# ------------------------------------------------------------------------------ report

def md(df: pd.DataFrame, nd: int = 2) -> str:
    dd = df.copy()
    for c in dd.columns:
        if dd[c].dtype.kind == "f":
            dd[c] = dd[c].round(nd)
    lines = ["| " + " | ".join(map(str, dd.columns)) + " |", "|" + "---|" * len(dd.columns)]
    for r in dd.itertuples(index=False):
        lines.append("| " + " | ".join("" if (isinstance(v, float) and np.isnan(v)) else str(v) for v in r) + " |")
    return "\n".join(lines)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    d = daily()
    t11 = e11(d)
    t12, coefs = e12(d)
    g = float(coefs["g·ΔW"][0])          # cm of laser per cm of water table (pooled, all years) → mm per cm ×10
    t13a, t13b, sim = e13(d, g * 10)
    best = t12.sort_values("rmse_mm").iloc[0]
    for n, t in (("e11_laser_level_models", t11), ("e12_change_models", t12), ("e13_water_table_proxies", t13a),
                 ("e13_driver_vs_laser_2024_2025", t13b)):
        t.to_csv(OUT / f"{n}.csv", index=False)
    pd.DataFrame({k: pd.Series(v) for k, v in coefs.items()}).to_csv(OUT / "e12_coefficients.csv", index=False)
    fig, ax = plt.subplots(2, 1, figsize=(11, 6.5), sharex=True)
    ax[0].plot(d.index, d.wtd_cm, color="C0", lw=0.8, label="P6 water table (logger)")
    ax[0].plot(sim.index, sim, color="C1", lw=0.8, label="bucket model (ERA5), 2024 test, 2025–26 out of sample")
    ax[0].set_ylabel("water table (cm)")
    ax[0].legend(fontsize=8)
    ax[1].plot(d.index, d.laser_cm, color="0.3", lw=0.8, label="laser (QC-passed daily mean)")
    ax[1].set_ylabel("laser surface (cm)")
    ax[1].legend(fontsize=8)
    ax[0].set_title("P6: the water table, its proxy where the logger stops, and the laser", fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / "fig_driver.png", dpi=150)
    plt.close(fig)
    txt = ["# X-074 — fusion v3 WP1: the truth and the driver (exploratory)", "",
           "Generated by `05_code/SAr-adapted/scripts/fusion_v3_wp1_x074.py` from `06_data/cube`; every number below is "
           "computed. Year-blocked: each year predicted from a fit on the other years. Laser = QC-passed hours (snow, "
           "outliers, filled stretches out; 2025 snow from the calibrated Open-Meteo mask), ≥ 12 hours a day.", "",
           "## E1.1 Does the laser carry the canopy? (daily level, demeaned per year)", "", md(t11), "",
           "Read: a canopy term that lowers the error **beyond water table + season** means the laser sees vegetation; "
           "a canopy term that only replaces the season does not.", "",
           "## E1.2 What drives the mat's 6-day height change?", "", md(t12), "",
           f"Lowest year-blocked error: **{best.iloc[0]}**. Pooled buoyancy g = {g * 10:.2f} mm of surface per cm of water "
           "table (all years; `e12_coefficients.csv`).", "",
           "## E1.3 A water table without the logger", "", "Proxies, trained 2020–2023, tested on 2024:", "", md(t13a), "",
           "The driver (proxy × g) against the laser's 6-day changes — 2024 (logger exists, proxy tested) and **2025 "
           "(no logger; out of sample)**:", "", md(t13b), "", "![driver](fig_driver.png)", ""]
    (OUT / "README.md").write_text("\n".join(txt))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()

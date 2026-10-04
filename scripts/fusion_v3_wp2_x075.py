"""X-075 — fusion v3 WP2: how much of the phase is motion, when is motion detectable, can the rest be removed?
(exploratory, Mac-safe, reads the cube; the phase referenced to the hard targets, WP0)

E2.1 Signal budget at P6. phase = m·(laser change, LOS) + Σ e_k·driver_k + ε for 6- and 12-day pairs, per track.
     Drivers: backscatter change (VH, VV), dew-point depression at each end, wet ends, top-soil moisture change, rain,
     air-temperature change, surface-temperature change near the overpass (ECOSTRESS, where present). Bayesian ridge
     on standardised drivers; every share from year-blocked out-of-fold predictions: motion alone, motion + drivers,
     and what remains against the WP0 noise model.
E2.2 Detectability. Does a pair's phase match the laser within λ/8? Features known without the laser: coherence,
     revisit, the driver's predicted |Δh| (cycle risk; buoyancy × water-table change, or the bucket model where no
     logger), surface state, season, backscatter change. Logistic regression and gradient boosting, leave-one-year-out;
     AUC, Brier score, calibration, permutation importance.
E2.3 Moisture correction. The fitted driver part subtracted (year-blocked); the corrected 6-day chain against the raw
     chain at P6 (pairs and chained levels). Gate for v3: the corrected chain beats the raw one on both tracks.

    PYTHONPATH=src python scripts/fusion_v3_wp2_x075.py
Outputs only to the hub (field values): ``08_deliverables/fusion_v3/wp2/``.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from sklearn.ensemble import HistGradientBoostingClassifier  # noqa: E402
from sklearn.inspection import permutation_importance  # noqa: E402
from sklearn.linear_model import BayesianRidge, LogisticRegression  # noqa: E402
from sklearn.metrics import brier_score_loss, roc_auc_score  # noqa: E402
from sklearn.pipeline import make_pipeline  # noqa: E402
from sklearn.preprocessing import StandardScaler  # noqa: E402

from insar_wetlands import fusion3 as f3  # noqa: E402
from insar_wetlands.cube.core import WAVELENGTH_MM, rad_to_los_mm  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fusion_v3_wp1_x074 as W1  # noqa: E402

HUB = Path(__file__).resolve().parents[3]
G, S = HUB / "06_data" / "cube" / "gold", HUB / "06_data" / "cube" / "silver"
OUT = HUB / "08_deliverables" / "fusion_v3" / "wp2"
WP0 = HUB / "08_deliverables" / "fusion_v3" / "wp0"
SEASON = {12: "DJF", 1: "DJF", 2: "DJF", 3: "MAM", 4: "MAM", 5: "MAM", 6: "JJA", 7: "JJA", 8: "JJA", 9: "SON",
          10: "SON", 11: "SON"}
DRIVERS = {"d_vh_db_3x3": "VH change", "d_vv_db_3x3": "VV change", "dpd_c_t1": "dew-point depression t1",
           "dpd_c_t2": "dew-point depression t2", "wet_ends": "wet ends (0–2)",
           "d_om_era5_land_soil_moisture_0_to_7cm": "top-soil moisture change", "rain_sum_mm": "rain over the pair",
           "d_air_c": "air-temperature change", "d_eco_lst_c": "surface-temperature change (ECOSTRESS ±6 h)"}


def truthy(x: pd.Series) -> pd.Series:
    return x.map({True: True, False: False, "True": True, "False": False}).fillna(False).astype(bool)


def p6_pairs() -> pd.DataFrame:
    d = pd.read_csv(G / "plots_pair.csv.gz", low_memory=False, parse_dates=["t1", "t2"])
    p = d[(d["plot"] == "P6") & (d.grid == "40m")].copy()
    p["frozen"] = truthy(p.flag_frozen_air_t1) | truthy(p.flag_frozen_air_t2)
    p["snow"] = (truthy(p.laser_snow_t1) | truthy(p.laser_snow_t2) | truthy(p.flag_snow_om_t1) |
                 truthy(p.flag_snow_om_t2))
    p["wet_ends"] = truthy(p.flag_wet_any_t1).astype(int) + truthy(p.flag_wet_any_t2).astype(int)
    p["state"] = np.where(p.frozen | p.snow, "frozen/snow", np.where(p.wet_ends > 0, "wet", "dry"))
    p["usable"] = p.d_laser_los_mm.notna() & ~p.frozen & ~p.snow & ~p.platforms.str.contains("S1D")
    p["y"] = -rad_to_los_mm(p.unw_1x1 - p.ref_hard_unw)              # WP0's reference
    p["year"] = p.t1.dt.year
    mid = p.t1 + (p.t2 - p.t1) / 2
    p["season"] = mid.dt.month.map(SEASON)
    p["doy"] = mid.dt.dayofyear
    p["d_eco_lst_c"] = eco_change(p)
    return p


def eco_change(p: pd.DataFrame) -> np.ndarray:
    """Surface-temperature change at P6 between the pair's dates from ECOSTRESS within ±6 h of each overpass."""
    a = pd.read_csv(G / "plots_acq.csv.gz", low_memory=False, parse_dates=["time_utc"])
    if "eco_lst_c_3x3" not in a:
        return np.full(len(p), np.nan)
    a = a[a["plot"] == "P6"].set_index(["track", "time_utc"]).eco_lst_c_3x3
    v1 = [a.get((t, x), np.nan) for t, x in zip(p.track, p.t1)]
    v2 = [a.get((t, x), np.nan) for t, x in zip(p.track, p.t2)]
    return np.asarray(v2, float) - np.asarray(v1, float)


def driver_prediction(p: pd.DataFrame) -> pd.Series:
    """Predicted Δh (mm, LOS) known without the laser: buoyancy × water-table change, the WP1 bucket where no logger."""
    d = W1.daily()
    g = 10 * float(W1.e12(d)[1]["g·ΔW"][0])
    _, _, sim = W1.e13(d, g)
    w = d.wtd_cm.combine_first(sim)
    day1, day2 = p.t1.dt.tz_convert(None).dt.normalize(), p.t2.dt.tz_convert(None).dt.normalize()
    dw = w.reindex(day2).to_numpy() - w.reindex(day1).to_numpy()
    return pd.Series(g * dw * p.los_u.to_numpy(), index=p.index)


# ------------------------------------------------------------------------------ E2.1 and E2.3

def budget(p: pd.DataFrame, noise: pd.DataFrame):
    rows, coef_rows, corrected = [], [], []
    fit = {(r.rev_class, r.state): r for r in noise.itertuples()}
    for (track, rv), g in p[p.usable & p.revisit_days.isin([6, 12])].groupby(["track", "revisit_days"]):
        g = g.copy()
        cols = [c for c in DRIVERS if c in g and g[c].notna().mean() >= 0.6]
        X = g[["d_laser_los_mm"] + cols].to_numpy(float)
        X = np.where(np.isfinite(X), X, np.nanmedian(X, 0))         # rare gaps filled with the median (no information)
        y = g.y.to_numpy()
        yrs = g.year.to_numpy()
        pm, pf, pdrv = (np.full(len(g), np.nan) for _ in range(3))
        for yv in np.unique(yrs):
            tr, te = yrs != yv, yrs == yv
            if tr.sum() < 12:
                continue
            m = f3.slope_origin(X[tr, 0], y[tr])
            pm[te] = m * X[te, 0]
            mdl = make_pipeline(StandardScaler(), BayesianRidge()).fit(X[tr], y[tr])
            pf[te] = mdl.predict(X[te])
            Xd = X[te].copy()
            Xd[:, 0] = 0
            Xm = np.zeros_like(X[te])
            pdrv[te] = mdl.predict(Xd) - mdl.predict(Xm)        # the drivers' part alone
        ok = np.isfinite(pm) & np.isfinite(pf)
        vy = np.var(y[ok])
        rc = "≤ 7 d" if rv <= 7 else "12 d"
        sig2 = np.nanmean([f3.crb_sigma_mm(c, fit[(rc, s)].looks, fit[(rc, s)].floor_mm) ** 2 if (rc, s) in fit else np.nan
                           for c, s in zip(g.coh_3x3[ok], g.state[ok])])
        r2m = 1 - np.var(y[ok] - pm[ok]) / vy
        r2f = 1 - np.var(y[ok] - pf[ok]) / vy
        rows.append({"track": track, "revisit_days": rv, "n": int(ok.sum()), "drivers_used": len(cols),
                     "share_motion": r2m, "share_added_by_drivers": r2f - r2m,
                     "share_noise (WP0 model)": sig2 / vy, "share_unexplained_beyond_noise": max(1 - max(r2m, r2f) - sig2 / vy, 0),
                     "resid_sd_motion_only_mm": float(np.std(y[ok] - pm[ok])), "resid_sd_with_drivers_mm": float(np.std(y[ok] - pf[ok]))})
        full = make_pipeline(StandardScaler(), BayesianRidge()).fit(X, y)
        for name, c in zip(["laser change"] + [DRIVERS[c] for c in cols], full[-1].coef_):
            coef_rows.append({"track": track, "revisit_days": rv, "term": name, "standardised_coef_mm": float(c)})
        g["y_corrected"] = y - np.where(np.isfinite(pdrv), pdrv, 0)
        g["oof_ok"] = np.isfinite(pdrv)
        corrected.append(g)
    return pd.DataFrame(rows), pd.DataFrame(coef_rows), pd.concat(corrected)


def single_driver_gate(p: pd.DataFrame) -> pd.DataFrame:
    """Each driver alone: y − (m·x + e·z) with m, e fitted on the other years; does it beat motion alone?"""
    rows = []
    for (track, rv), g in \
            p[p.usable & p.revisit_days.isin([6])].groupby(["track", "revisit_days"]):
        x, y, yrs = g.d_laser_los_mm.to_numpy(), g.y.to_numpy(), g.year.to_numpy()
        base = np.nanstd(f3.year_blocked_residual(x, y, yrs))
        for c, lab in DRIVERS.items():
            if c not in g or g[c].notna().mean() < 0.6:
                continue
            z = g[c].to_numpy(float)
            z = np.where(np.isfinite(z), z, np.nanmedian(z))
            res = np.full(len(g), np.nan)
            for yv in np.unique(yrs):
                tr, te = yrs != yv, yrs == yv
                if tr.sum() < 12:
                    continue
                A = np.column_stack([x[tr], z[tr], np.ones(tr.sum())])
                b_ = np.linalg.lstsq(A, y[tr], rcond=None)[0]
                res[te] = y[te] - np.column_stack([x[te], z[te], np.ones(te.sum())]) @ b_
            rows.append({"track": track, "revisit_days": rv, "driver": lab, "resid_sd_motion_only_mm": float(base),
                         "resid_sd_motion_plus_driver_mm": float(np.nanstd(res)),
                         "helps": bool(np.nanstd(res) < base)})
    return pd.DataFrame(rows)


def correction_gate(c: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (track, rv), g in c[c.oof_ok].groupby(["track", "revisit_days"]):
        for name, col in (("raw", "y"), ("moisture-corrected", "y_corrected")):
            x, y = g.d_laser_los_mm.to_numpy(), g[col].to_numpy()
            res = f3.year_blocked_residual(x, y, g.year.to_numpy())
            run, (ly, lx) = f3.chain_runs(pd.DatetimeIndex(g.t1), pd.DatetimeIndex(g.t2), y, x, min_len=4)
            rows.append({"track": track, "revisit_days": rv, "series": name, "n": len(g), "r": float(np.corrcoef(x, y)[0, 1]),
                         "resid_sd_year_blocked_mm": float(np.nanstd(res)), **f3.level_metrics(ly, lx)})
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------------ E2.2

FEATURES = ["coh_3x3", "revisit_days", "abs_driver_dh_mm", "wet_ends", "is_summer", "abs_d_vh_db", "is_desc"]


def detectability(p: pd.DataFrame):
    q = p[p.usable & (p.revisit_days <= 24)].copy()
    q["match"] = ((q.y - q.d_laser_los_mm).abs() < WAVELENGTH_MM / 8).astype(int)
    q["abs_driver_dh_mm"] = q.driver_dh_mm.abs()
    q["abs_d_vh_db"] = q.d_vh_db_3x3.abs()
    q["is_summer"] = (q.season == "JJA").astype(int)
    q["is_desc"] = (q.track == "descending").astype(int)
    X = q[FEATURES].to_numpy(float)
    X = np.where(np.isfinite(X), X, np.nanmedian(X, 0))
    y = q.match.to_numpy()
    yrs = q.year.to_numpy()
    models = {"logistic": lambda: make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000, C=1.0)),
              "gradient boosting": lambda: HistGradientBoostingClassifier(max_depth=3, learning_rate=0.05, max_iter=200,
                                                                            l2_regularization=1.0, random_state=75)}
    rows, cal, imp = [], [], []
    q["p_logistic"] = np.nan
    for name, mk in models.items():
        prob = np.full(len(q), np.nan)
        for yv in np.unique(yrs):
            tr, te = yrs != yv, yrs == yv
            if len(np.unique(y[tr])) < 2:
                continue
            prob[te] = mk().fit(X[tr], y[tr]).predict_proba(X[te])[:, 1]
        ok = np.isfinite(prob)
        rows.append({"model": name, "n": int(ok.sum()), "base_rate": float(y[ok].mean()),
                     "auc_leave_one_year_out": float(roc_auc_score(y[ok], prob[ok])),
                     "brier": float(brier_score_loss(y[ok], prob[ok])),
                     "brier_climatology": float(brier_score_loss(y[ok], np.full(ok.sum(), y[ok].mean())))})
        b = pd.cut(prob[ok], [0, 0.2, 0.4, 0.6, 0.8, 1.0])
        cc = pd.DataFrame({"bin": b, "obs": y[ok]}).groupby("bin", observed=True).obs.agg(["mean", "size"]).reset_index()
        cc.insert(0, "model", name)
        cal.append(cc)
        full = mk().fit(X, y)
        pi = permutation_importance(full, X, y, n_repeats=30, random_state=75, scoring="roc_auc")
        for f, m_, s_ in zip(FEATURES, pi.importances_mean, pi.importances_std):
            imp.append({"model": name, "feature": f, "auc_drop_mean": float(m_), "auc_drop_sd": float(s_)})
        if name == "logistic":
            q["p_logistic"] = prob
    return pd.DataFrame(rows), pd.concat(cal, ignore_index=True), pd.DataFrame(imp), q


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
    p = p6_pairs()
    p["driver_dh_mm"] = driver_prediction(p)
    noise = pd.read_csv(WP0 / "noise_model_fits.csv")
    bud, coefs, corr = budget(p, noise)
    gate = correction_gate(corr)
    single = single_driver_gate(p)
    det, cal, imp, q = detectability(p)
    for n, t in (("e21_signal_budget", bud), ("e21_coefficients", coefs), ("e23_correction_gate", gate), ("e23_single_driver_gate", single),
                 ("e22_detectability_skill", det), ("e22_calibration", cal), ("e22_importance", imp)):
        t.to_csv(OUT / f"{n}.csv", index=False)
    q[["track", "pair", "t1", "t2", "revisit_days", "match", "p_logistic"] + FEATURES].to_csv(OUT / "e22_pairs.csv", index=False)
    six = gate[gate.revisit_days == 6].pivot_table(index="track", columns="series", values="resid_sd_year_blocked_mm")
    passed = bool(len(six) == 2 and (six["moisture-corrected"] < six["raw"]).all())
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    b6 = bud[bud.revisit_days == 6].set_index("track")
    parts = ["share_motion", "share_added_by_drivers", "share_noise (WP0 model)", "share_unexplained_beyond_noise"]
    lab = ["motion (laser)", "moisture / surface drivers", "noise (WP0 model)", "unexplained"]
    bottom = np.zeros(len(b6))
    for c, l_, col in zip(parts, lab, ["C3", "C0", "0.6", "0.85"]):
        v = np.clip(b6[c].to_numpy(), 0, None)
        ax[0].bar(b6.index, v, bottom=bottom, color=col, label=l_)
        bottom += v
    ax[0].set_ylabel("share of the 6-day phase variance at P6")
    ax[0].legend(fontsize=7)
    ax[0].set_title("Signal budget, 6-day pairs (year-blocked)", fontsize=9)
    cl = cal[cal.model == "logistic"]
    ax[1].plot([0, 1], [0, 1], "k-", lw=0.6)
    ax[1].plot([iv.mid for iv in cl.bin], cl["mean"], "o-", color="C2")
    ax[1].set_xlabel("predicted probability that the pair carries the motion")
    ax[1].set_ylabel("observed share matching the laser")
    ax[1].set_title("Detectability (logistic, leave-one-year-out)", fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / "fig_budget_detectability.png", dpi=150)
    plt.close(fig)
    eco_n = int(p.d_eco_lst_c.notna().sum())
    txt = ["# X-075 — fusion v3 WP2: signal budget, detectability, moisture correction (exploratory)", "",
           "Generated by `05_code/SAr-adapted/scripts/fusion_v3_wp2_x075.py` from `06_data/cube`; every number is computed. "
           "P6 phase referenced to the hard targets (WP0), −φ·λ/4π, uplift positive; usable = laser present, not frozen, "
           "no snow, no early-S1D pair. All skill year-blocked (leave one year out).", "",
           "## E2.1 Signal budget", "", md(bud), "",
           "Shares: `share_motion` = year-blocked R² of m·laser alone; `share_added_by_drivers` = what the drivers add "
           "out of fold; `share_noise` = the WP0 noise model's variance over the phase variance; the rest is unexplained.", "",
           "Standardised coefficients (mm per SD of each term, all data):", "", md(coefs), "",
           f"ECOSTRESS surface-temperature change available for {eco_n} P6 pairs (used only where ≥ 60 % of a group has it).", "",
           "## E2.2 Detectability", "", md(det), "", "Calibration:", "", md(cal), "", "Permutation importance (AUC drop):", "",
           md(imp), "",
           "## E2.3 Moisture correction — the gate", "", md(gate), "",
           f"**Gate (corrected 6-day chain lowers the year-blocked residual on both tracks): {'passed' if passed else 'not passed'}.**", "",
           "Each driver alone (6-day; slope and driver coefficient fitted on the other year):", "", md(single), "",
           f"Drivers that help on both tracks: {', '.join(sorted(set.intersection(*[set(g[g.helps].driver) for _, g in single.groupby('track')]))) or 'none'}.", "",
           "![budget and detectability](fig_budget_detectability.png)", ""]
    (OUT / "README.md").write_text("\n".join(txt))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()

"""X-062 — the hybrid fusion system: every dataset combined into a mat-motion map, judged by the field.

Branch fusion-x062 (DECISIONS 2026-09-28): the fusion block chains per-pixel consecutive 12-day
phase changes, which SETTLED §1 excludes on main. Blocks (insar_wetlands.fusion):

A  water model   — bucket model per plot from station rain and Hamon PET; leave-one-year-out skill,
                   against a gradient-boosting benchmark on the same inputs.
B  zones         — k-means on per-pixel behaviour (coherence, phase noise, backscatter, S2 NDWI,
                   DEM, distance to the mat edge, water sensitivity); dissimilarity to the plots.
C  transfer      — Gaussian process: each pixel's water-table response relative to the mat plots,
                   judged leave-one-plot-out.
D  phase physics — at P6, the 12-day LOS change = m·laser + e·ΔWTD + c; motion-only / moisture-only /
                   both compared by year-out cross-validation; buoyancy g (surface mm per cm WTD).
E  fusion        — per pixel and consecutive pair: forcing g·r·ΔWTD_ref (variance from P6) fused
                   with the radar increment (moisture removed, scaled by m; variance from coherence).
Validation: P6 laser held out year by year; P7–P9 leave-one-plot-out; stable ground.

Inputs: field data (hub 06_data/field), the first-deliverable tables, X-058 laser flags, X-050
wetness flags, the Drive snapshot (interferograms, RTC, S2, EVD). Writes only to the hub.

    INSAR_DRIVE_ROOT=../local/drive_pristine PYTHONPATH=src ../local/.venv/bin/python scripts/field_fusion_x062.py
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import xarray as xr  # noqa: E402
from scipy import ndimage  # noqa: E402

from insar_wetlands import field  # noqa: E402
from insar_wetlands import field_link as fl  # noqa: E402
from insar_wetlands import fusion as fu  # noqa: E402
from insar_wetlands.bootstrap import start  # noqa: E402
from insar_wetlands.inversion.isbas import PHASE_TO_MM  # noqa: E402
from insar_wetlands.stack import load_layer  # noqa: E402

HUB = Path(__file__).resolve().parents[3]
DLV = HUB / "08_deliverables"
LAT = 52.76
INC = {"ascending": 32.26, "descending": 39.17}
TRACKS = ("ascending", "descending")
PLOTS = [f"P{i}" for i in range(1, 10)]
MAT = ["P6", "P7", "P8", "P9"]                  # the supervisor: P6–P9 are the floating-mat plots
Q_STABLE_MM = 2.0                               # prior sd of a 12-day increment off the mat (method parameter)
COLOR = {"ascending": "tab:orange", "descending": "tab:blue"}


# ------------------------------------------------------------------------------ loading

def load(ctx):
    acq = pd.read_csv(DLV / "field_first" / "s1_acquisitions.csv")
    px = pd.read_csv(DLV / "field_first" / "plot_pixels.csv").set_index("plot")
    wtd_s1 = pd.read_csv(DLV / "field_first" / "wtd_at_s1.csv").set_index(["date", "track", "plot"])
    ov = pd.read_csv(DLV / "field_p6" / "laser_qc" / "laser_at_s1.csv").set_index(["date", "track"])
    runs = pd.read_csv(DLV / "field_p6" / "laser_qc" / "interpolated.csv", parse_dates=["start", "end"])
    runs = runs.assign(start=runs.start.dt.tz_convert(None), end=runs.end.dt.tz_convert(None))
    wet = pd.read_csv(DLV / "field_dew_x050" / "wetness_at_overpasses_2020_2024.csv").set_index(["date", "track"])
    crop = {"ascending": ctx.paths.cropped, "descending": ctx.paths.for_phase("field_fusion_x062", track="descending").cropped}
    zones = {k: ctx.zones[k].values.astype(bool) for k in "ABCD"}
    cube = {}
    for t in TRACKS:
        dates = sorted(pd.to_datetime(acq[acq.track == t].date).unique())
        pairs = [(a, b) for a, b in zip(dates[:-1], dates[1:])]
        ids = [f"{a:%Y%m%d}_{b:%Y%m%d}" for a, b in pairs]
        unw = load_layer(crop[t], "unw_phase", ids).values
        coh = load_layer(crop[t], "corr", ids).values
        refC = np.array([np.nanmedian(u[zones["C"] & np.isfinite(u)]) for u in unw])
        cube[t] = {"dates": dates, "pairs": pairs, "ids": ids, "dlos": (unw - refC[:, None, None]) * PHASE_TO_MM, "coh": coh}
    return acq, px, wtd_s1, ov, runs, wet, zones, cube


# ------------------------------------------------------------------------------ A. water model

def daily_frame(w: pd.DataFrame) -> pd.DataFrame:
    d = pd.DataFrame({"rain": w["Rain_mm_Tot"].resample("D").sum(), "air": w["Air_2m"].resample("D").mean()})
    d["pet"] = fu.hamon_pet(d.air, d.index.dayofyear, LAT)
    for p in PLOTS:
        cens = field.censored_flag(w[p])
        d[p] = w[p].resample("D").mean()
        d[f"{p}_ok"] = ~cens.resample("D").max().astype(bool)
    return d.dropna(subset=["rain", "air"])


def ml_features(d: pd.DataFrame, wcol: str) -> pd.DataFrame:
    f = pd.DataFrame(index=d.index)
    for k in range(4):
        f[f"rain_l{k}"] = d.rain.shift(k)
    f["rain7"], f["rain30"] = d.rain.rolling(7).sum(), d.rain.rolling(30).sum()
    f["pet"], f["pet7"] = d.pet, d.pet.rolling(7).sum()
    f["doy_s"], f["doy_c"] = np.sin(2 * np.pi * d.index.dayofyear / 365.25), np.cos(2 * np.pi * d.index.dayofyear / 365.25)
    f["w"] = d[wcol]
    return f


def block_a(d: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    from sklearn.ensemble import HistGradientBoostingRegressor
    rows, sims = [], {}
    years = sorted(set(d.index.year))
    for p in PLOTS:
        ok = d[f"{p}_ok"].to_numpy()
        floor = float(d[p][d[f"{p}_ok"]].min())
        feats = ml_features(d, p)
        target = d[p].shift(-1) - d[p]
        for y in years:
            tr, te = d.index.year != y, d.index.year == y
            par = fu.fit_bucket(d[p].where(tr), d.rain, d.pet, valid=ok & tr)
            dy = d[te]
            first = np.flatnonzero(dy[f"{p}_ok"].to_numpy())
            if len(first) < 60:
                continue
            i0 = first[0]
            sim_b = fu.simulate_bucket(par, dy[p].iloc[i0], dy.rain.iloc[i0:], dy.pet.iloc[i0:], floor=floor)
            m = (tr & ok & feats.notna().all(axis=1).to_numpy() & target.notna().to_numpy())
            gb = HistGradientBoostingRegressor(max_iter=300, learning_rate=0.05, random_state=0).fit(feats[m], target[m])
            fy = feats[te].iloc[i0:].copy()
            w = dy[p].iloc[i0]
            sim_m = np.empty(len(fy))
            for k in range(len(fy)):
                sim_m[k] = w
                row = fy.iloc[[k]].copy()
                row["w"] = w
                w = max(w + float(gb.predict(row.fillna(0))[0]), floor)
            obs = dy[p].iloc[i0:].where(dy[f"{p}_ok"].iloc[i0:])
            rows.append({"plot": p, "year": y, "n_days": int(obs.notna().sum()),
                         "nse_bucket": fu.nse(obs, sim_b), "nse_boosting": fu.nse(obs, sim_m),
                         "rmse_bucket_cm": float(np.sqrt(np.nanmean((obs - sim_b) ** 2))),
                         "rmse_boosting_cm": float(np.sqrt(np.nanmean((obs - sim_m) ** 2)))})
            sims[(p, y)] = (obs.index, obs.to_numpy(), sim_b, sim_m)
    cv = pd.DataFrame(rows)
    full = pd.DataFrame([{"plot": p, **fu.fit_bucket(d[p], d.rain, d.pet, valid=d[f"{p}_ok"].to_numpy())} for p in PLOTS])
    return cv, {"full": full, "sims": sims}


# ------------------------------------------------------------------------------ B. zones

def pixel_features(ctx, cube, zones) -> tuple[dict, np.ndarray]:
    D = ctx.paths.drive
    ny, nx = zones["A"].shape
    flat = lambda a: np.asarray(a, float).reshape(len(a), -1)  # noqa: E731
    F = {}
    for t in TRACKS:
        c = cube[t]
        mids = [a + (b - a) / 2 for a, b in c["pairs"]]
        ty = fl.years_since(mids)
        F[f"coh_mean_{t[:3]}"] = np.nanmean(c["coh"], axis=0)
        F[f"dlos_sd_{t[:3]}"] = np.nanstd(c["dlos"], axis=0)
        if t == "ascending":
            F["coh_amp_asc"] = fu.harmonic_amplitude_map(ty, flat(c["coh"])).reshape(ny, nx)
    F["tcoh_evd"] = ctx.to_grid(xr.open_dataset(D / "phaseE2_evd.nc")["temporal_coherence"]).values
    rtc = xr.open_dataset(D / "rtc_dualpol_stack.nc")
    vv, vh = ctx.to_grid(rtc["gamma0_vv_db"]), ctx.to_grid(rtc["gamma0_vh_db"])
    F["vv_mean"], F["vh_mean"] = vv.mean("time").values, vh.mean("time").values
    F["vv_amp"] = fu.harmonic_amplitude_map(fl.years_since(vv.time.values), flat(vv.values)).reshape(ny, nx)
    s2 = ctx.to_grid(xr.open_dataset(D / "s2_stack.nc")["ndwi"])
    F["ndwi_mean"] = s2.mean("time").values
    F["ndwi_amp"] = fu.harmonic_amplitude_map(fl.years_since(s2.time.values), flat(s2.values)).reshape(ny, nx)
    F["dem"] = ctx.dem.values if hasattr(ctx.dem, "values") else np.asarray(ctx.dem)
    inside = ndimage.distance_transform_edt(zones["A"]) * 40.0
    outside = ndimage.distance_transform_edt(~zones["A"]) * 40.0
    F["dist_mat_edge_m"] = np.where(zones["A"], inside, -outside)
    pm = xr.open_dataset(DLV / "field_pixel_x053" / "pixel_maps.nc")
    F["coh_wtd_r_asc"] = pm["coh_dwtd_r_ascending"].values
    F["wet_penalty_desc"] = pm["wet_penalty_descending"].values
    X = np.column_stack([F[k].reshape(-1) for k in F])
    return F, X


def block_b(F: dict, X: np.ndarray, zones, px) -> dict:
    ny, nx = zones["A"].shape
    ok = np.isfinite(X).all(axis=1)
    Z, _, _ = fu.standardize(X[ok])
    labels, k, scores = fu.choose_kmeans(Z)
    lab = np.full(ny * nx, -1)
    lab[ok] = labels
    idx_ok = np.flatnonzero(ok)
    plot_flat = sorted({int(r.row) * nx + int(r.col) for r in px.itertuples()})
    train = np.array([np.flatnonzero(idx_ok == f)[0] for f in plot_flat if f in set(idx_ok)])
    di, thr = fu.dissimilarity_index(Z, train)
    DI = np.full(ny * nx, np.nan)
    DI[ok] = di
    zone_of = np.full(ny * nx, "", dtype=object)
    for z, m in zones.items():
        zone_of[m.reshape(-1)] = z
    tab = pd.crosstab(pd.Series(lab[ok], name="cluster"), pd.Series(zone_of[ok], name="zone"))
    plot_cluster = {p: int(lab[int(r.row) * nx + int(r.col)]) for p, r in px.iterrows()}
    return {"labels": lab.reshape(ny, nx), "k": k, "silhouette": scores, "DI": DI.reshape(ny, nx), "thr": thr,
            "crosstab": tab, "plot_cluster": plot_cluster, "features": list(F)}


# ------------------------------------------------------------------------------ C. transfer

def wtd_pair_matrix(wtd_s1, cube, track) -> pd.DataFrame:
    """ΔWTD (cm) per consecutive pair (rows) and plot (columns); NaN where either end is censored."""
    out = {}
    for p in PLOTS:
        v = []
        for a, b in cube[track]["pairs"]:
            ra, rb = wtd_s1.loc[(f"{a:%Y-%m-%d}", track, p)], wtd_s1.loc[(f"{b:%Y-%m-%d}", track, p)]
            v.append(np.nan if (ra.wtd_censored or rb.wtd_censored) else rb.wtd_at - ra.wtd_at)
        out[p] = v
    return pd.DataFrame(out, index=cube[track]["ids"])


def block_c(F: dict, px, dW: pd.DataFrame, zones) -> dict:
    ny, nx = zones["A"].shape
    ref = dW[MAT].mean(axis=1)
    ratio = {p: float(dW[p].std() / ref[dW[p].notna()].std()) for p in PLOTS}
    units = px.groupby(["row", "col"]).apply(lambda g: list(g.index), include_groups=False)
    feats = ["dist_mat_edge_m", "coh_mean_asc", "ndwi_mean"]
    Xu = np.array([[F[f][r, c] for f in feats] for (r, c) in units.index])
    yu = np.array([np.mean([ratio[p] for p in ps]) for ps in units])
    Xall = np.column_stack([F[f].reshape(-1) for f in feats])
    okall = np.isfinite(Xall).all(axis=1)
    g = fu.gp_transfer(Xu, yu, np.where(okall[:, None], Xall, 0))
    mean = np.where(okall, g["mean"], np.nan).reshape(ny, nx)
    sd = np.where(okall, g["sd"], np.nan).reshape(ny, nx)
    unit_tab = pd.DataFrame({"plots": [",".join(ps) for ps in units], "row": [r for r, _ in units.index],
                             "col": [c for _, c in units.index], "ratio": yu, "loo_pred": g["loo_pred"]})
    return {"ratio": ratio, "units": unit_tab, "features": feats, "mean": mean, "sd": sd, "gp": g, "ref": ref,
            "gp_mean_map": mean, "gp_sd_map": sd}


# ------------------------------------------------------------------------------ D. phase physics

def p6_pairs(track, cube, px, ov, runs, dW, wet) -> pd.DataFrame:
    r, c = int(px.loc["P6", "row"]), int(px.loc["P6", "col"])
    rows = []
    for k, ((a, b), pid) in enumerate(zip(cube[track]["pairs"], cube[track]["ids"])):
        la, lb = ov.loc[(f"{a:%Y-%m-%d}", track)], ov.loc[(f"{b:%Y-%m-%d}", track)]
        across = bool(((runs.end > a) & (runs.start < b)).any())
        usable = bool(la.usable and lb.usable and not across)
        dh = (lb.laser_cm - la.laser_cm) * 10 if usable else np.nan
        rows.append({"k": k, "pair": pid, "t1": a, "t2": b, "year": b.year, "dh_mm": dh,
                     "laser_los": dh * np.cos(np.radians(INC[track])), "dwtd": dW.loc[pid, "P6"],
                     "s1": cube[track]["dlos"][k, r, c], "coh": cube[track]["coh"][k, r, c],
                     "dry": not (wet.loc[(f"{a:%Y-%m-%d}", track)].wet or wet.loc[(f"{b:%Y-%m-%d}", track)].wet
                                 or wet.loc[(f"{a:%Y-%m-%d}", track)].frozen or wet.loc[(f"{b:%Y-%m-%d}", track)].frozen)})
    return pd.DataFrame(rows)


MODELS = {"motion + moisture": ["laser_los", "dwtd"], "motion only": ["laser_los"], "moisture only": ["dwtd"]}


def fit_phase(d: pd.DataFrame, cols) -> dict:
    """Through the origin: no surface change and no water change give no phase change. (An
    intercept would be subtracted from every pair and chained into a drift the raw phase does
    not have — the first run did exactly that.)"""
    X = np.column_stack([d[c] for c in cols])
    return fu.bayes_linear(X, d.s1.to_numpy(), fu.phase_sigma_mm(d.coh.to_numpy()))


def block_d(P6: dict) -> dict:
    coef, cv, buoy = [], [], {}
    for t, d in P6.items():
        ok = d.dropna(subset=["laser_los", "dwtd", "s1", "coh"])
        for sub, dd in (("all", ok), ("dry", ok[ok.dry])):
            if len(dd) < 8:
                continue
            f = fit_phase(dd, MODELS["motion + moisture"])
            coef.append({"track": t, "subset": sub, "n": len(dd), "m": f["beta"][0], "m_sd": f["sd"][0],
                         "e_mm_per_cm": f["beta"][1], "e_sd": f["sd"][1], "noise_scale": f["noise_scale"],
                         "motion_detected": bool(f["beta"][0] - 2 * f["sd"][0] > 0)})
        for name, cols in MODELS.items():
            pred = np.full(len(ok), np.nan)
            for y in sorted(ok.year.unique()):
                tr, te = ok.year != y, ok.year == y
                if tr.sum() < 8 or te.sum() == 0:
                    continue
                f = fit_phase(ok[tr], cols)
                pred[te.to_numpy()] = np.column_stack([ok[te][c] for c in cols]) @ f["beta"]
            m = np.isfinite(pred)
            cv.append({"track": t, "model": name, "n": int(m.sum()),
                       "rmse_mm": float(np.sqrt(np.mean((pred[m] - ok.s1.to_numpy()[m]) ** 2))),
                       "r": float(np.corrcoef(pred[m], ok.s1.to_numpy()[m])[0, 1]),
                       "rmse_no_model_mm": float(np.sqrt(np.mean(ok.s1.to_numpy()[m] ** 2)))})
        w = d.dropna(subset=["dh_mm", "dwtd"])
        sl = np.polyfit(w.dwtd, w.dh_mm, 1)
        buoy[t] = {"g_mm_per_cm": float(sl[0]), "resid_sd_mm": float(np.std(w.dh_mm - np.polyval(sl, w.dwtd))), "n": len(w)}
    both = pd.concat([d.dropna(subset=["dh_mm", "dwtd"]) for d in P6.values()])
    sl = np.polyfit(both.dwtd, both.dh_mm, 1)
    buoy["both"] = {"g_mm_per_cm": float(sl[0]), "resid_sd_mm": float(np.std(both.dh_mm - np.polyval(sl, both.dwtd))), "n": len(both)}
    return {"coef": pd.DataFrame(coef), "cv": pd.DataFrame(cv), "buoyancy": buoy}


# ------------------------------------------------------------------------------ E. fusion

def run_fusion(track, cube, zones, dW_ref, r_map, r_sd, phys, g, q0, restrict=None):
    """Fused increments for every pixel (or the (row, col) in ``restrict``). Returns dict of
    (steps × pixels) arrays: forcing u, radar z, fused, variances, gate masks."""
    c = cube[track]
    cos = np.cos(np.radians(INC[track]))
    m, m_sd, e, e_sd = phys["m"], phys["m_sd"], phys["e"], phys["e_sd"]
    ny, nx = zones["A"].shape
    if restrict is None:
        sel = (slice(None), slice(None))
        mat, lake = zones["A"], zones["B"]
        r, rsd = r_map, r_sd
    else:
        sel = restrict
        mat, lake = zones["A"][sel], zones["B"][sel]
        r = r_map if np.ndim(r_map) == 0 else r_map[sel]        # a single pixel may pass its own scalar
        rsd = r_sd if np.ndim(r_sd) == 0 else r_sd[sel]
    dlos, coh = c["dlos"][(slice(None),) + sel], c["coh"][(slice(None),) + sel]
    ref = dW_ref.to_numpy()[:, None, None] if np.ndim(mat) else dW_ref.to_numpy()
    rr = np.where(np.isfinite(r), np.clip(r, 0.2, 3.0), 1.0)
    dW = rr * ref                                                   # cm, per pixel and pair
    u = np.where(mat, g * dW, 0.0)
    Q = np.where(mat, q0 ** 2 + (g * np.nan_to_num(rsd) * ref) ** 2, Q_STABLE_MM ** 2)
    Q = np.where(np.isfinite(ref), Q, np.where(mat, 10.0 ** 2, Q))  # no water record: wide forcing
    u = np.where(np.isfinite(u), u, 0.0)
    corr_moist = np.where(mat, e * np.nan_to_num(dW), 0.0)
    scale = np.where(mat, m * cos, cos)
    z = (dlos - corr_moist) / scale
    sig = fu.phase_sigma_mm(coh) * phys["noise_scale"]
    R = (sig / scale) ** 2 + np.where(mat, ((e_sd * np.nan_to_num(dW)) / scale) ** 2 + (z * m_sd / m) ** 2, 0.0)
    R = np.where(np.isfinite(z), R, np.nan)
    if not phys["motion_detected"]:                              # this track's phase does not carry the motion:
        R = np.where(mat, np.nan, R)                             # on the mat it is not a motion observation
    out = fu.fuse_increments(u, Q, z, R)
    out.update({"u": u, "Q": Q, "z": z, "R": R, "lake": lake})
    return out


def validate(track, cube, zones, px, dW, P6d, blockc, phys_by_year, buoy_by_year, phys_all, g_all, q_all):
    rows, series = [], {}
    r6, c6 = int(px.loc["P6", "row"]), int(px.loc["P6", "col"])
    # V1: P6, laser held out year by year (m, e, c, g, q fitted without that year)
    d = P6d[track]
    for y in sorted(d.year.unique()):
        ph, (g, q) = phys_by_year[(track, y)], buoy_by_year[(track, y)]
        if ph is None:
            continue
        o = run_fusion(track, cube, zones, dW["P6"], np.array(1.0), np.array(0.0), ph, g, q, restrict=(r6, c6))  # P6's own water
        te = (d.year == y) & d.dh_mm.notna()
        truth = d.dh_mm[te].to_numpy()
        for name, est in (("model only", o["u"]), ("radar only", o["z"]), ("fused", o["inc"])):
            e_ = np.asarray(est)[te.to_numpy()]
            m = np.isfinite(e_)
            rows.append({"track": track, "test": "V1 P6 laser, year held out", "fold": y, "estimate": name, "n": int(m.sum()),
                         "rmse_mm": float(np.sqrt(np.mean((e_[m] - truth[m]) ** 2))) if m.any() else np.nan,
                         "r": float(np.corrcoef(e_[m], truth[m])[0, 1]) if m.sum() > 2 else np.nan})
    # V2: P7–P9, leave that plot out of the forcing and the transfer
    units = blockc["units"]
    for p in ["P7", "P8", "P9"]:
        rp, cp = int(px.loc[p, "row"]), int(px.loc[p, "col"])
        ref = dW[[q for q in MAT if q != p]].mean(axis=1)
        u_row = units[units.plots.str.contains(p)].iloc[0]
        o = run_fusion(track, cube, zones, ref, np.array(u_row.loo_pred), np.array(0.3), phys_all, g_all, q_all, restrict=(rp, cp))
        truth = (g_all * dW[p]).to_numpy()
        ok = np.isfinite(truth)
        for name, est in (("model only", o["u"]), ("radar only", o["z"]), ("fused", o["inc"])):
            e_ = np.asarray(est)[ok]
            m = np.isfinite(e_)
            rows.append({"track": track, "test": "V2 plot left out (surface from its own WTD)", "fold": p, "estimate": name,
                         "n": int(m.sum()), "rmse_mm": float(np.sqrt(np.mean((e_[m] - truth[ok][m]) ** 2))),
                         "r": float(np.corrcoef(e_[m], truth[ok][m])[0, 1])})
    return pd.DataFrame(rows), series


# ------------------------------------------------------------------------------ figures

def extent(ctx):
    x, y = ctx.template.x.values, ctx.template.y.values
    return [x[0] - 20, x[-1] + 20, y[-1] - 20, y[0] + 20]


def mat_window(zones, pad=6):
    rr, cc = np.nonzero(zones["A"])
    return slice(max(rr.min() - pad, 0), rr.max() + pad + 1), slice(max(cc.min() - pad, 0), cc.max() + pad + 1)


def draw_plots(ax, px, ctx, sl=None):
    x, y = ctx.template.x.values, ctx.template.y.values
    for p, r in px.iterrows():
        ax.plot(x[int(r.col)], y[int(r.row)], "k^" if p in MAT else "kv", ms=5)
    ax.plot([], [], "k^", label="P6–P9 (mat)")
    ax.plot([], [], "kv", label="P1–P5 (reference)")


def fig_a(out, cv, sims):
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    g = cv.groupby("plot")[["nse_bucket", "nse_boosting"]].median().reindex(PLOTS)
    xx = np.arange(len(PLOTS))
    ax[0].bar(xx - 0.2, g.nse_bucket, 0.4, label="bucket (physics)")
    ax[0].bar(xx + 0.2, g.nse_boosting, 0.4, label="gradient boosting (ML)")
    ax[0].axhline(0, color="k", lw=0.6)
    ax[0].set_xticks(xx, PLOTS)
    ax[0].set_ylim(max(-2, np.nanmin(g.values) - 0.2), 1)
    ax[0].set_ylabel("NSE, year held out (median over years)")
    ax[0].legend(fontsize=8)
    ax[0].set_title("A · water model skill per plot", fontsize=9)
    for a, p in zip(ax[1:], ["P6", "P9"]):
        y = 2023 if (p, 2023) in sims else sorted(k[1] for k in sims if k[0] == p)[0]
        t, obs, sb, sm = sims[(p, y)]
        a.plot(t, obs, "k", lw=1, label="measured")
        a.plot(t, sb, lw=1, label="bucket")
        a.plot(t, sm, lw=1, label="boosting")
        a.set_title(f"{p}, {y} held out", fontsize=9)
        a.set_ylabel("WTD (cm)")
        a.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(out / "fig_A_water_model.png", dpi=150)
    plt.close(fig)


def fig_b(out, B, ctx, px, zones):
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch
    fig, ax = plt.subplots(1, 2, figsize=(14, 5.8), layout="constrained")
    ext = extent(ctx)
    x, y = ctx.template.x.values, ctx.template.y.values
    cols = ["#60a5fa", "#fb923c", "#22c55e", "#a855f7", "#f43f5e", "#eab308", "#14b8a6", "#78716c"][:B["k"]]
    lab = np.where(B["labels"] < 0, np.nan, B["labels"]).astype(float)
    ax[0].imshow(lab, extent=ext, cmap=ListedColormap(cols), vmin=-0.5, vmax=B["k"] - 0.5, interpolation="nearest")
    ax[0].contour(x, y, zones["A"].astype(float), [0.5], colors="k", linewidths=0.9)
    draw_plots(ax[0], px, ctx)
    mat_zone = B["plot_cluster"]["P6"]
    ax[0].legend(handles=[Patch(color=c, label=f"zone {i}" + (" — where the plots are" if i == mat_zone else "")) for i, c in enumerate(cols)]
                 + [Patch(color="white", ec="k", label="no complete features")], fontsize=8, loc="lower left")
    ax[0].set_title(f"B · {B['k']} behaviour zones (k-means; black: the mat, zone A)", fontsize=10)
    di = np.log10(np.clip(B["DI"], 1e-3, None))
    im = ax[1].imshow(di, extent=ext, cmap="viridis_r", interpolation="nearest")
    ax[1].contour(x, y, (B["DI"] > B["thr"]).astype(float), [0.5], colors="w", linewidths=0.7)
    ax[1].contour(x, y, zones["A"].astype(float), [0.5], colors="r", linewidths=0.8)
    draw_plots(ax[1], px, ctx)
    ax[1].legend(fontsize=8, loc="lower left")
    ax[1].set_title("How unlike the plot pixels each pixel is (white line: the threshold; red: the mat)", fontsize=10)
    fig.colorbar(im, ax=ax[1], shrink=0.8, label="log10 dissimilarity index")
    for a in ax:
        a.set_xticks([])
        a.set_yticks([])
    fig.savefig(out / "fig_B_zones.png", dpi=150)
    plt.close(fig)


def fig_c(out, C, ctx, px, zones):
    fig, ax = plt.subplots(1, 3, figsize=(16, 5))
    u = C["units"]
    ax[0].scatter(u.ratio, u.loo_pred, color="k")
    for _, r in u.iterrows():
        ax[0].annotate(r.plots, (r.ratio, r.loo_pred), fontsize=8, xytext=(3, 3), textcoords="offset points")
    lo, hi = min(u.ratio.min(), u.loo_pred.min()) - 0.1, max(u.ratio.max(), u.loo_pred.max()) + 0.1
    ax[0].plot([lo, hi], [lo, hi], "k--", lw=0.7)
    ax[0].set(xlabel="water-table response relative to the mat plots (measured)", ylabel="predicted with that plot left out")
    ax[0].set_title(f"C · leave-one-plot-out, skill {C['gp']['skill']:.2f} (1 = perfect, 0 = the mean)", fontsize=9)
    sl = mat_window(zones)
    ext = extent(ctx)
    for a, arr, t, cm in ((ax[1], C["gp_mean_map"], "GP-predicted response ratio", "magma"), (ax[2], C["gp_sd_map"], "its uncertainty (sd)", "Greys")):
        im = a.imshow(np.where(zones["A"] | zones["B"], arr, np.nan), extent=ext, cmap=cm, interpolation="nearest")
        a.contour(ctx.template.x.values, ctx.template.y.values, zones["A"].astype(float), [0.5], colors="c", linewidths=0.6)
        draw_plots(a, px, ctx)
        x, y = ctx.template.x.values, ctx.template.y.values
        a.set_xlim(x[sl[1].start], x[sl[1].stop - 1])
        a.set_ylim(y[sl[0].stop - 1], y[sl[0].start])
        a.set_title(t + " (mat and lake)", fontsize=9)
        a.set_xticks([])
        a.set_yticks([])
        plt.colorbar(im, ax=a, shrink=0.8)
    fig.tight_layout()
    fig.savefig(out / "fig_C_transfer.png", dpi=150)
    plt.close(fig)


def fig_d(out, D, P6d):
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
    cf = D["coef"]
    for i, (_, r) in enumerate(cf.iterrows()):
        ax[0].errorbar(r.m, i, xerr=2 * r.m_sd, fmt="o", color=COLOR[r.track])
        ax[1].errorbar(r.e_mm_per_cm, i, xerr=2 * r.e_sd, fmt="o", color=COLOR[r.track])
    labels = [f"{r.track[:4]} {r.subset} (n {r.n})" for _, r in cf.iterrows()]
    for a, t, ref in ((ax[0], "m — how much of the laser motion the phase carries", 1), (ax[1], "e — moisture term, mm LOS per cm ΔWTD", 0)):
        a.set_yticks(range(len(labels)), labels, fontsize=8)
        a.axvline(ref, color="k", ls="--", lw=0.8)
        a.set_title(t, fontsize=9)
    cv = D["cv"]
    for j, t in enumerate(TRACKS):
        c = cv[cv.track == t].set_index("model")
        xx = np.arange(len(c)) + j * 0.35
        ax[2].bar(xx, c.rmse_mm, 0.33, color=COLOR[t], label=t)
        ax[2].scatter(xx, c.rmse_no_model_mm, marker="_", s=200, color="k")
    ax[2].set_xticks(np.arange(len(MODELS)) + 0.17, list(MODELS), fontsize=8)
    ax[2].set_ylabel("RMSE predicting the 12-day phase, year held out (mm)")
    ax[2].set_title("D · which physics explains the phase? (black: no model)", fontsize=9)
    ax[2].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "fig_D_phase_physics.png", dpi=150)
    plt.close(fig)


def fig_e_maps(out, fused, dates, ctx, px, zones, B):
    track = "ascending"
    cum = fused[track]["cum"]
    sl = mat_window(zones)
    x, y = ctx.template.x.values, ctx.template.y.values
    ext = extent(ctx)
    picks = np.linspace(8, len(dates[track]) - 1, 6).astype(int)
    fig, ax = plt.subplots(2, 4, figsize=(19, 9), layout="constrained")
    vmax = np.nanpercentile(np.abs(cum[:, zones["A"]]), 97)
    for a, k in zip(ax.flat[:6], picks):
        im = a.imshow(np.where(zones["B"], np.nan, cum[k]), extent=ext, cmap="RdBu", vmin=-vmax, vmax=vmax, interpolation="nearest")
        a.set_title(f"{pd.Timestamp(dates[track][k]):%Y-%m-%d}", fontsize=10)
    fig.colorbar(im, ax=ax[:, :3], shrink=0.55, location="left", label="fused surface displacement since the first date (mm, + up)")
    im2 = ax.flat[6].imshow(np.where(zones["B"], np.nan, fused[track]["amp"]), extent=ext, cmap="viridis", interpolation="nearest")
    ax.flat[6].set_title("annual amplitude (mm)", fontsize=10)
    fig.colorbar(im2, ax=ax.flat[6], shrink=0.8)
    im3 = ax.flat[7].imshow(np.where(zones["B"], np.nan, np.sqrt(fused[track]["cum_var"][-1])), extent=ext, cmap="Greys", interpolation="nearest")
    ax.flat[7].contour(x, y, (B["DI"] > B["thr"]).astype(float), [0.5], colors="m", linewidths=0.6)
    ax.flat[7].set_title("uncertainty at the end (sd, mm)\nmagenta: unlike any plot", fontsize=10)
    fig.colorbar(im3, ax=ax.flat[7], shrink=0.8)
    for a in ax.flat:
        a.contour(x, y, zones["A"].astype(float), [0.5], colors="k", linewidths=0.5)
        draw_plots(a, px, ctx)
        a.set_xlim(x[sl[1].start], x[sl[1].stop - 1])
        a.set_ylim(y[sl[0].stop - 1], y[sl[0].start])
        a.set_xticks([])
        a.set_yticks([])
    fig.suptitle("E · fused surface displacement of the mat, ascending (lake masked; ▲ P6–P9, ▼ P1–P5)", fontsize=12)
    fig.savefig(out / "fig_E_maps.png", dpi=130)
    plt.close(fig)


def p6_levels(fused, dates, px, ov, runs) -> pd.DataFrame:
    """V1b — the three-year level at P6 against the laser (ascending), within each laser segment
    (the laser's level across a filled stretch is unknown): shape (r) and amplitude."""
    r6, c6 = int(px.loc["P6", "row"]), int(px.loc["P6", "col"])
    f = fused["ascending"]
    d = pd.to_datetime(dates["ascending"])
    s = pd.DataFrame({"fused": f["cum"][:, r6, c6],
                      "radar only": np.concatenate([[0], np.nancumsum(f["z"][:, r6, c6])]),
                      "model only": np.concatenate([[0], np.cumsum(f["u"][:, r6, c6])])}, index=d)
    o = ov.reset_index()
    o = o[(o.track == "ascending") & o.usable]
    o.index = pd.to_datetime(o.date)
    o["seg"] = [(runs.end < t).sum() for t in o.index]
    j = s.join(o[["laser_cm", "seg"]], how="inner")
    j["laser"] = j.laser_cm * 10
    rows = []
    b = j.groupby("seg").laser.transform(lambda v: v - v.mean())
    for col in ("model only", "radar only", "fused"):
        a = j.groupby("seg")[col].transform(lambda v: v - v.mean())
        rows.append({"estimate": col, "dates": len(j), "laser_segments": int(j.seg.nunique()), "r_with_laser_level": float(a.corr(b)),
                     "amplitude_vs_laser": float(a.std() / b.std()), "rmse_mm": float(np.sqrt(((a - b) ** 2).mean()))})
    return pd.DataFrame(rows)


def fig_e_validation(out, val, P6d, fused, px, dates, stable):
    fig, ax = plt.subplots(2, 2, figsize=(15, 9))
    for a, test in zip(ax[0], ["V1 P6 laser, year held out", "V2 plot left out (surface from its own WTD)"]):
        v = val[val.test == test].groupby(["track", "estimate"]).rmse_mm.mean().unstack()
        v = v[["model only", "radar only", "fused"]]
        v.plot.bar(ax=a, rot=0)
        a.set_ylabel("RMSE of the 12-day surface change (mm)")
        a.set_title(test, fontsize=9)
    t = "ascending"
    r6, c6 = int(px.loc["P6", "row"]), int(px.loc["P6", "col"])
    d = pd.to_datetime(dates[t])
    f = fused[t]
    ax[1, 0].plot(d, f["cum"][:, r6, c6], color="k", lw=1.5, label="fused")
    ax[1, 0].plot(d, np.concatenate([[0], np.nancumsum(f["z"][:, r6, c6])]), lw=0.9, label="radar only (chained)")
    ax[1, 0].plot(d, np.concatenate([[0], np.cumsum(f["u"][:, r6, c6])]), lw=0.9, label="model only (water forcing)")
    ax[1, 0].set_title("P6, ascending: fused vs its two inputs (laser changes in V1)", fontsize=9)
    ax[1, 0].set_ylabel("surface since the first date (mm)")
    ax[1, 0].legend(fontsize=8)
    s = stable
    ax[1, 1].bar(np.arange(len(s)), s.radar_only_sd_mm, 0.4, label="radar only")
    ax[1, 1].bar(np.arange(len(s)) + 0.4, s.fused_sd_mm, 0.4, label="fused")
    ax[1, 1].set_xticks(np.arange(len(s)) + 0.2, [f"{r.track[:4]} {r.zone}" for r in s.itertuples()], fontsize=8)
    ax[1, 1].set_ylabel("median over pixels of the series' SD (mm)")
    ax[1, 1].set_title("V3 · how much each zone 'moves' (stable ground should be small)", fontsize=9)
    ax[1, 1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "fig_E_validation.png", dpi=150)
    plt.close(fig)


# ------------------------------------------------------------------------------ README

def f2(v, n=2):
    return "—" if v is None or not np.isfinite(v) else f"{v:.{n}f}"


def readme(out, A, B, C, D, val, stable, phys, g, q, zones, lev):
    L = ["# X-062 — the hybrid fusion system (exploratory, branch fusion-x062)\n",
         "Generated by `05_code/SAr-adapted/scripts/field_fusion_x062.py` (outputs here only). Every number below is "
         "computed by the script. **Branch only**: the fusion step chains per-pixel consecutive 12-day phase changes, "
         "which `SETTLED.md` §1 excludes on main (DECISIONS 2026-09-28). No network inversion is run.\n",
         "## In short\n", *in_short(A, B, C, D, val, lev), "",
         "## What goes in\n",
         "| Source | Used as |", "|---|---|",
         "| Station rain and air temperature (hourly → daily), Hamon PET | A: drivers of the water model |",
         "| WTD at 9 plots (hourly), dry-well flags | A: target; C: each plot's response; E: forcing (P6–P9) |",
         "| P6 laser (X-058: measured, snow-free, no filled stretch) | D: the motion term and the buoyancy; V1: the judge |",
         "| Sentinel-1 consecutive 12-day pairs, both tracks: phase (P6 pixel minus grassland) and coherence | B: features; D; E: radar observation and its noise |",
         "| RTC VV/VH, S2 NDWI, EVD temporal coherence, DEM, distance to the mat edge, X-053 water sensitivity | B: zones; C: transfer features |",
         "| Wet/frozen flags at the overpasses (X-050) | D: dry subset |\n",
         "## A · Water model (bucket vs gradient boosting, year held out)\n",
         "| plot | NSE bucket (median over years) | NSE boosting | RMSE bucket (cm) | RMSE boosting (cm) |", "|---|---|---|---|---|"]
    ga = A.groupby("plot")[["nse_bucket", "nse_boosting", "rmse_bucket_cm", "rmse_boosting_cm"]].median().reindex(PLOTS)
    for p, r in ga.iterrows():
        L.append(f"| {p} | {f2(r.nse_bucket)} | {f2(r.nse_boosting)} | {f2(r.rmse_bucket_cm, 1)} | {f2(r.rmse_boosting_cm, 1)} |")
    L += ["\nNSE 1 = perfect, 0 = no better than the plot's mean over the held-out year. Weather is one station, "
          "so the water model explains *when* the water moves, not *where* it differs.\n",
          f"## B · Zones\n\n{B['k']} zones (k-means, silhouette {', '.join(f'k={k}: {v:.2f}' for k, v in B['silhouette'].items())}) "
          f"from {len(B['features'])} per-pixel features ({', '.join(B['features'])}). Plots' zones: "
          + ", ".join(f"{p} {z}" for p, z in B["plot_cluster"].items()) + ".\n",
          f"Pixels whose features are unlike any plot pixel (dissimilarity > {f2(B['thr'])}): "
          f"mat {pct(B, zones, 'A')}, lake {pct(B, zones, 'B')}, grassland {pct(B, zones, 'C')}, outside {pct(B, zones, 'D')} — "
          "results there are extrapolation.\n",
          "Pixels by zone (rows) and project zone (columns):\n", md(B["crosstab"]), "",
          "## C · Transfer of each plot's water-table response\n",
          f"Trait: the plot's 12-day ΔWTD standard deviation relative to the mean of P6–P9. Gaussian process on "
          f"{', '.join(C['features'])}; leave-one-plot-out RMSE {f2(C['gp']['rmse_loo'], 3)} vs {f2(C['gp']['rmse_mean_baseline'], 3)} "
          f"for the mean → skill **{f2(C['gp']['skill'])}** (1 perfect, ≤ 0 no better than the mean). "
          + ("Used in E." if C["used"] else "**Not used in E**: the transfer is no better than the mean, so every mat pixel takes "
             "the mat plots' mean response (ratio 1) with their spread as its uncertainty.") + "\n",
          md(C["units"].round(3), index=False), "",
          "## D · Phase physics at P6: motion, moisture, or both?\n",
          "Posterior of s1 = m·laser_LOS + e·ΔWTD, through the origin (coherence-weighted; ±1 sd):\n",
          "| track | subset | n | m | e (mm per cm) | noise scale | motion detected (m − 2 sd > 0) |", "|---|---|---|---|---|---|---|"]
    for r in D["coef"].itertuples():
        L.append(f"| {r.track} | {r.subset} | {r.n} | {f2(r.m)} ± {f2(r.m_sd)} | {f2(r.e_mm_per_cm)} ± {f2(r.e_sd)} | {f2(r.noise_scale)} | {'yes' if r.motion_detected else 'no'} |")
    L += ["\nPredicting the 12-day phase with a year held out:\n",
          "| track | model | n | RMSE (mm) | r | RMSE with no model = zero change (mm) |", "|---|---|---|---|---|---|"]
    for r in D["cv"].itertuples():
        L.append(f"| {r.track} | {r.model} | {r.n} | {f2(r.rmse_mm, 1)} | {f2(r.r)} | {f2(r.rmse_no_model_mm, 1)} |")
    L += ["\nBuoyancy (laser surface change per cm of WTD change, consecutive pairs): "
          + "; ".join(f"{k} {f2(v['g_mm_per_cm'])} mm/cm (residual sd {f2(v['resid_sd_mm'], 1)} mm, n {v['n']})" for k, v in D["buoyancy"].items()) + ".\n",
          "## E · Fusion and validation\n",
          "Used in the maps (all data): "
          + "; ".join(f"{t} m {f2(phys[t]['m'])}, e {f2(phys[t]['e'])} mm/cm — radar "
                      + ("used on the mat" if phys[t]["motion_detected"] else "**not** used on the mat (its phase does not carry the motion)")
                      for t in TRACKS)
          + f"; buoyancy g {f2(g)} mm/cm, forcing sd at P6 {f2(q, 1)} mm; off the mat the forcing is 0 ± {Q_STABLE_MM} mm.\n",
          "| track | test | estimate | mean RMSE (mm) | mean r | folds |", "|---|---|---|---|---|---|"]
    vs = val.groupby(["track", "test", "estimate"]).agg(rmse=("rmse_mm", "mean"), r=("r", "mean"), folds=("fold", "nunique")).reset_index()
    for r in vs.itertuples():
        L.append(f"| {r.track} | {r.test} | {r.estimate} | {f2(r.rmse, 1)} | {f2(r.r)} | {r.folds} |")
    L += ["\nV1b — the three-year **level** at P6 against the laser (ascending, within laser segments):\n",
          md(lev.round(2), index=False),
          "\nThe fused series follows the *shape* of the laser best, but at a fraction of its *amplitude*: its water "
          "forcing is a regression prediction, which is shrunk toward zero, and the fusion leans on it where the radar is "
          "noisy. For how much the mat moves, read the radar-only chain (with its noise); for when and where, the fused map.\n"]
    L += ["\nV3 — how much each zone moves (median over pixels of the series' SD). Off the mat the fused value is close "
          "to its prior (no motion) by construction; the test there is the radar-only column — the chained noise level "
          "that any mat signal has to exceed:\n", md(stable.round(2), index=False), "",
          "## How to read it\n",
          "- **The fused map is the radar's 12-day chain, corrected for moisture and scaled with what P6 taught, pulled toward "
          "the water forcing where the radar is noisy.** Its skill is what V1 says at P6 and V2 at P7–P9 — nowhere else is "
          "it tested.",
          "- V2's 'truth' is each plot's own water table × the P6 buoyancy: a consistency check, not a measurement.",
          "- The whole mat is assumed to float and to wet like P6 (one m, e, g). The dissimilarity map says where that is "
          "extrapolation; C says how well the plots' differences can be carried to other pixels.",
          "- Coherence, the lake, dew and snow set the radar's weight pixel by pixel; unwrapping outliers are rejected by the gate.\n",
          "## Files\n",
          "| File | Content |", "|---|---|",
          "| `water_model_cv.csv`, `water_model_params.csv` | A |", "| `zones_crosstab.csv`, `zones.json` | B |",
          "| `transfer_units.csv` | C |", "| `phase_physics_coef.csv`, `phase_physics_cv.csv`, `buoyancy.json` | D |",
          "| `validation.csv`, `stable_ground.csv` | E validation |",
          "| `fusion_ascending.nc`, `fusion_descending.nc` | fused / radar-only / model-only displacement, sd, zones, dissimilarity (not in git) |",
          "| `fig_A_…png` … `fig_E_…png` | one figure per block |"]
    (out / "README.md").write_text("\n".join(L) + "\n")


def in_short(A, B, C, D, val, lev) -> list[str]:
    """The answer in six lines, every number read from the tables."""
    ca = D["coef"].set_index(["track", "subset"])
    asc, des = ca.loc[("ascending", "all")], ca.loc[("descending", "all")]
    cv = D["cv"][D["cv"].track == "ascending"].set_index("model")
    v1 = val[(val.test.str.startswith("V1")) & (val.track == "ascending")].groupby("estimate").rmse_mm.mean()
    lv = lev.set_index("estimate")
    ga = A.groupby("plot")[["nse_bucket", "nse_boosting"]].median()
    z = B["plot_cluster"]["P6"]
    ct = B["crosstab"]
    beyond = int(ct.loc[z, [c for c in ct.columns if c in ("C", "D")]].sum()) if z in ct.index else 0
    return [
        f"- **D, the physics.** On 12-day ascending pairs at P6 the phase carries the laser motion (m {f2(asc.m)} ± {f2(asc.m_sd)}) "
        f"plus an opposite-sign water term (e {f2(asc.e_mm_per_cm)} ± {f2(asc.e_sd)} mm per cm). Year held out: motion + moisture "
        f"RMSE {f2(cv.loc['motion + moisture', 'rmse_mm'], 1)} mm, motion only {f2(cv.loc['motion only', 'rmse_mm'], 1)}, "
        f"moisture only {f2(cv.loc['moisture only', 'rmse_mm'], 1)}, no model {f2(cv.loc['motion only', 'rmse_no_model_mm'], 1)}. "
        f"Descending: no motion (m {f2(des.m)} ± {f2(des.m_sd)}) — dew.",
        f"- **E, the fusion.** At P6, years held out, the 12-day surface change: fused RMSE {f2(v1.get('fused'), 1)} mm, radar only "
        f"{f2(v1.get('radar only'), 1)}, water model only {f2(v1.get('model only'), 1)}. Over three years the fused level follows the "
        f"laser's shape best (r {f2(lv.loc['fused', 'r_with_laser_level'])}) but at {lv.loc['fused', 'amplitude_vs_laser']:.0%} of its "
        f"amplitude; the radar chain has {lv.loc['radar only', 'amplitude_vs_laser']:.0%} of it with a noisier shape.",
        f"- **C, the transfer, fails its test** (skill {f2(C['gp']['skill'])}): the plots' water-table differences cannot be "
        "carried to other pixels from what the satellites see; the maps use one mat-wide forcing.",
        f"- **B, the zones.** All plots fall in zone {z}; {beyond} pixels outside the mat and lake behave like it (a second "
        "mat-like area to ask the field team about). The plots sit on one line, so most of the mat is formally 'unlike any plot'.",
        f"- **A, the water.** A physical bucket model beats gradient boosting at {int((ga.nse_bucket > ga.nse_boosting).sum())} of "
        f"{len(ga)} plots with a year held out — one weather station is too little for machine learning to add anything.",
        "- **So:** a map of *when and where* the mat moves (ascending, 2022–2024), calibrated at P6 and tested there; its "
        "amplitude is conservative, and away from the plots it rests on the assumption that the mat floats and wets like P6.",
    ]


def md(df: pd.DataFrame, index: bool = True) -> str:
    """A markdown table without extra dependencies."""
    d = df.reset_index() if index else df
    fmt = lambda v: f"{v:.3f}" if isinstance(v, float) else str(v)  # noqa: E731
    head = "| " + " | ".join(map(str, d.columns)) + " |"
    return "\n".join([head, "|" + "---|" * len(d.columns)] + ["| " + " | ".join(fmt(v) for v in row) + " |" for row in d.itertuples(index=False)])


def pct(B, zones, z):
    m = zones[z] & np.isfinite(B["DI"])
    return f"{np.mean(B['DI'][m] > B['thr']):.0%}" if m.any() else "—"


# ------------------------------------------------------------------------------ main

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default=str(DLV / "field_fusion_x062"))
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    ctx = start("field_fusion_x062", mount=False, git=False)
    acq, px, wtd_s1, ov, runs, wet, zones, cube = load(ctx)
    ny, nx = zones["A"].shape
    print("loaded", {t: len(cube[t]["ids"]) for t in TRACKS}, "pairs")

    # A
    w = field.load_wtd_hourly()
    d = daily_frame(w)
    A, Ax = block_a(d)
    A.to_csv(out / "water_model_cv.csv", index=False)
    Ax["full"].to_csv(out / "water_model_params.csv", index=False)
    fig_a(out, A, Ax["sims"])
    print("A done")

    # B
    F, X = pixel_features(ctx, cube, zones)
    B = block_b(F, X, zones, px)
    B["crosstab"].to_csv(out / "zones_crosstab.csv")
    (out / "zones.json").write_text(json.dumps({"k": B["k"], "silhouette": B["silhouette"], "threshold": B["thr"],
                                                "plot_zone": B["plot_cluster"], "features": B["features"]}, indent=1))
    fig_b(out, B, ctx, px, zones)
    print("B done, k =", B["k"])

    # C (ascending pairs define the trait)
    dW = {t: wtd_pair_matrix(wtd_s1, cube, t) for t in TRACKS}
    C = block_c(F, px, dW["ascending"], zones)
    C["ref_by_track"] = {t: dW[t][MAT].mean(axis=1) for t in TRACKS}
    C["used"] = C["gp"]["skill"] > 0                             # a block that fails its own test does not feed the next
    if not C["used"]:
        mat_r = np.array([C["ratio"][p] for p in MAT])
        C["mean"] = np.where(np.isfinite(C["mean"]), 1.0, np.nan)
        C["sd"] = np.where(np.isfinite(C["sd"]), float(mat_r.std()), np.nan)
    C["units"].to_csv(out / "transfer_units.csv", index=False)
    fig_c(out, C, ctx, px, zones)
    print("C done, skill", round(C["gp"]["skill"], 3))

    # D
    P6d = {t: p6_pairs(t, cube, px, ov, runs, dW[t], wet) for t in TRACKS}
    D = block_d(P6d)
    D["coef"].to_csv(out / "phase_physics_coef.csv", index=False)
    D["cv"].to_csv(out / "phase_physics_cv.csv", index=False)
    (out / "buoyancy.json").write_text(json.dumps(D["buoyancy"], indent=1))
    fig_d(out, D, P6d)
    print("D done")

    # parameters for E: all data, and without each year (for V1)
    def phys_of(dd):
        dd = dd.dropna(subset=["laser_los", "dwtd", "s1", "coh"])
        if len(dd) < 8:
            return None
        f = fit_phase(dd, MODELS["motion + moisture"])
        return {"m": float(f["beta"][0]), "m_sd": float(f["sd"][0]), "e": float(f["beta"][1]), "e_sd": float(f["sd"][1]),
                "noise_scale": f["noise_scale"], "motion_detected": bool(f["beta"][0] - 2 * f["sd"][0] > 0)}

    def buoy_of(dd):
        ww = dd.dropna(subset=["dh_mm", "dwtd"])
        sl = np.polyfit(ww.dwtd, ww.dh_mm, 1)
        return float(sl[0]), float(np.std(ww.dh_mm - np.polyval(sl, ww.dwtd)))

    both = pd.concat(P6d.values())
    g_all, q_all = buoy_of(both)
    phys = {t: phys_of(P6d[t]) for t in TRACKS}
    phys_by_year, buoy_by_year = {}, {}
    for t in TRACKS:
        for y in sorted(P6d[t].year.unique()):
            phys_by_year[(t, y)] = phys_of(P6d[t][P6d[t].year != y])
            buoy_by_year[(t, y)] = buoy_of(both[both.year != y])

    # E: maps
    fused, dates = {}, {}
    for t in TRACKS:
        o = run_fusion(t, cube, zones, C["ref_by_track"][t], C["mean"], C["sd"], phys[t], g_all, q_all)
        dates[t] = cube[t]["dates"]
        ty = fl.years_since(dates[t])
        cum = o["cum"]
        cum[:, zones["B"]] = np.nan
        o["amp"] = fu.harmonic_amplitude_map(ty, cum.reshape(len(ty), -1)).reshape(ny, nx)
        fused[t] = o
        radar_cum = np.concatenate([np.zeros((1, ny, nx)), np.nancumsum(o["z"], axis=0)])
        model_cum = np.concatenate([np.zeros((1, ny, nx)), np.cumsum(o["u"], axis=0)])
        ds = xr.Dataset({"fused_mm": (("time", "y", "x"), cum.astype("float32")),
                         "fused_sd_mm": (("time", "y", "x"), np.sqrt(o["cum_var"]).astype("float32")),
                         "radar_only_mm": (("time", "y", "x"), radar_cum.astype("float32")),
                         "model_only_mm": (("time", "y", "x"), model_cum.astype("float32")),
                         "amplitude_mm": (("y", "x"), o["amp"].astype("float32")),
                         "zone": (("y", "x"), np.select([zones[k] for k in "ABCD"], list("ABCD"), "").astype("U1")),
                         "dissimilarity": (("y", "x"), B["DI"].astype("float32")),
                         "radar_step_used": (("step", "y", "x"), o["used"].astype("int8"))},
                        coords={"time": pd.to_datetime(dates[t]), "y": ctx.template.y.values, "x": ctx.template.x.values},
                        attrs={"track": t, "note": "X-062 branch fusion-x062; exploratory; lake masked; + up, mm"})
        ds.to_netcdf(out / f"fusion_{t}.nc")
    # validation
    val = pd.concat([validate(t, cube, zones, px, dW[t], P6d, C, phys_by_year, buoy_by_year, phys[t], g_all, q_all)[0] for t in TRACKS])
    val.to_csv(out / "validation.csv", index=False)
    st = []
    for t in TRACKS:
        o = fused[t]
        radar_cum = np.concatenate([np.zeros((1, ny, nx)), np.nancumsum(o["z"], axis=0)])
        for z in "ACD":
            m = zones[z]
            st.append({"track": t, "zone": z, "pixels": int(m.sum()),
                       "radar_only_sd_mm": float(np.nanmedian(np.nanstd(radar_cum[:, m], axis=0))),
                       "fused_sd_mm": float(np.nanmedian(np.nanstd(o["cum"][:, m], axis=0))),
                       "radar_steps_used": float(np.mean(o["used"][:, m])),
                       "radar_steps_rejected": float(np.mean(o["rejected"][:, m]))})
    stable = pd.DataFrame(st)
    stable.to_csv(out / "stable_ground.csv", index=False)
    lev = p6_levels(fused, dates, px, ov, runs)
    lev.to_csv(out / "p6_levels.csv", index=False)
    fig_e_maps(out, fused, dates, ctx, px, zones, B)
    fig_e_validation(out, val, P6d, fused, px, dates, stable)
    readme(out, A, B, C, D, val, stable, phys, g_all, q_all, zones, lev)
    print("wrote", out)


if __name__ == "__main__":
    main()

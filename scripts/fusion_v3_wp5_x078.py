"""X-078 — fusion v3 WP5: the estimator, assembled only from what passed WP0–WP4 (exploratory, branch fusion-x062).

Model (vertical surface height h at every acquisition instant of either track):
    h(t) = g·W(t) + r(t),   r ~ Ornstein–Uhlenbeck(τ, σ_r)                       (WP1: buoyancy is the driver)
    y_k / u_k = h(t_j) − h(t_i) + ε_k,   ε_k ~ N(0, σ_eff,k²)                      (6-day pairs, both tracks)
    y_k: P6 phase referenced to the hard targets, −φ·λ/4π (WP0); u_k: the look vector's up component;
    σ_eff = √(σ_noise² + (1 − p_detect)·(λ/2)²): the WP0 noise model inflated by the WP2 detectability;
    no moisture term (WP2 gate failed); SNAPHU's cycles kept for 6-day pairs (WP4); W = the P6 logger, the WP1 bucket
    model where no logger (2025–2026).
Solved exactly (Kalman/RTS smoother as one sparse system, ``fusion3.ou_increment_smoother``).
Space: the raft c(t) from the mat's median phase, and the WP3 edge gradient h(x, t) = c(t)·(1 + s·z(x)), z = the
standardised depth in the mat, s = WP3's depth slope (6-day, both tracks).

Validation at P6 against the laser, year-blocked: for a test year, g, τ, σ_r and the detectability model come from the
other years. Levels demeaned within laser segments; r, RMSE, amplitude ratio, ±1σ/±2σ coverage. Compared with the
plain 6-day chain (same reference) and the prior alone; fusion v2's published numbers read from its table.

    PYTHONPATH=src python scripts/fusion_v3_wp5_x078.py
Outputs only to the hub: ``08_deliverables/fusion_v3/wp5/`` (tables, series, mat height maps).
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import xarray as xr  # noqa: E402
from sklearn.linear_model import LogisticRegression  # noqa: E402
from sklearn.pipeline import make_pipeline  # noqa: E402
from sklearn.preprocessing import StandardScaler  # noqa: E402

from insar_wetlands import fusion3 as f3  # noqa: E402
from insar_wetlands.cube.core import WAVELENGTH_MM, rad_to_los_mm  # noqa: E402
from insar_wetlands.fusion2 import ou_from_series  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fusion_v3_wp1_x074 as W1  # noqa: E402
import fusion_v3_wp2_x075 as W2  # noqa: E402

HUB = Path(__file__).resolve().parents[3]
G, S = HUB / "06_data" / "cube" / "gold", HUB / "06_data" / "cube" / "silver"
D3 = HUB / "08_deliverables" / "fusion_v3"
OUT = D3 / "wp5"
RUNS = {"R1 2020–2021": (2020, 2021), "R3 2025": (2025, 2025), "R4 2026": (2026, 2026)}
LASER_YEARS = (2021, 2022, 2023, 2024, 2025)
FEATS = ["coh_3x3", "revisit_days", "abs_driver_dh_mm", "wet_ends", "is_summer", "abs_d_vh_db", "is_desc"]


# ------------------------------------------------------------------------------ inputs

def inputs():
    d = W1.daily()
    p = W2.p6_pairs()
    p["driver_dh_mm"] = W2.driver_prediction(p)
    p["abs_driver_dh_mm"] = p.driver_dh_mm.abs()
    p["abs_d_vh_db"] = p.d_vh_db_3x3.abs()
    p["is_summer"] = (p.season == "JJA").astype(int)
    p["is_desc"] = (p.track == "descending").astype(int)
    p["match"] = ((p.y - p.d_laser_los_mm).abs() < WAVELENGTH_MM / 8).astype(float)
    p.loc[p.d_laser_los_mm.isna(), "match"] = np.nan
    noise = pd.read_csv(D3 / "wp0" / "noise_model_fits.csv")
    return d, p, noise


def params(d: pd.DataFrame, exclude_year: int | None) -> dict:
    """g (mm/cm), τ and σ_r of the laser residual, from the laser years other than ``exclude_year``."""
    keep = (d.year != exclude_year) if exclude_year else np.ones(len(d), bool)
    dd = d[keep]
    c = W1.changes(dd)
    ok = c.d_laser_cm.notna() & c.d_wtd_cm.notna()
    g = float(np.sum(c.d_laser_cm[ok] * c.d_wtd_cm[ok]) / np.sum(c.d_wtd_cm[ok] ** 2)) * 10      # mm per cm
    res = []
    for _, gg in dd.groupby("year"):
        r = gg.laser_cm * 10 - g * gg.wtd_cm
        res.append(r - r.mean())
    rr = pd.concat(res)
    t = (rr.index - rr.index[0]).days.to_numpy(float)
    ou = ou_from_series(t, rr.to_numpy(), max_lag_days=60)
    return {"g_mm_per_cm": g, "tau_days": ou["tau_days"], "sigma_r_mm": ou["sigma"]}


def detect_model(p: pd.DataFrame, exclude_year: int | None):
    q = p[p.usable & (p.revisit_days <= 24) & p.match.notna()]
    if exclude_year:
        q = q[q.year != exclude_year]
    X = q[FEATS].to_numpy(float)
    med = np.nanmedian(X, 0)
    X = np.where(np.isfinite(X), X, med)
    m = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000)).fit(X, q.match.astype(int))
    return m, med


def driver_series(d: pd.DataFrame, g: float) -> pd.Series:
    _, _, sim = W1.e13(d, g)
    return d.wtd_cm.combine_first(sim)


# ------------------------------------------------------------------------------ one run

MISS_SCALES = {"cycle-scale misses (λ/2, first run)": WAVELENGTH_MM / 2,
               "miss scale λ/4 (6-day misses are not cycle slips, WP4)": WAVELENGTH_MM / 4}
PRIMARY = "miss scale λ/4 (6-day misses are not cycle slips, WP4)"


def observations(p: pd.DataFrame, noise: pd.DataFrame, years: tuple[int, int], dm, med,
                 miss_scale: float = WAVELENGTH_MM / 4) -> pd.DataFrame:
    fit = {(r.rev_class, r.state): r for r in noise.itertuples()}
    q = p[(p.t1.dt.year >= years[0]) & (p.t2.dt.year <= years[1]) & (p.revisit_days <= 7)].copy()
    q = q[~q.platforms.str.contains("S1D") & (q.state != "frozen/snow") & q.y.notna() & q.los_u.notna()]
    X = np.where(np.isfinite(q[FEATS].to_numpy(float)), q[FEATS].to_numpy(float), med)
    q["p_detect"] = dm.predict_proba(X)[:, 1]
    q["sigma_noise_mm"] = [f3.crb_sigma_mm(c, fit[("≤ 7 d", s)].looks, fit[("≤ 7 d", s)].floor_mm)
                           for c, s in zip(q.coh_3x3, q.state)]
    q["sigma_eff_mm"] = f3.detect_inflated_sigma(q.sigma_noise_mm, q.p_detect, miss_scale_mm=miss_scale)
    q["y_vert_mm"] = q.y / q.los_u
    q["sigma_vert_mm"] = q.sigma_eff_mm / q.los_u
    return q


def run_p6(q: pd.DataFrame, W: pd.Series, par: dict) -> pd.DataFrame:
    nodes = pd.DatetimeIndex(sorted(set(q.t1) | set(q.t2)))
    t_days = ((nodes - nodes[0]).total_seconds() / 86400).to_numpy()
    idx = {t: k for k, t in enumerate(nodes)}
    i, j = q.t1.map(idx).to_numpy(), q.t2.map(idx).to_numpy()
    w_nodes = W.reindex(nodes.tz_convert(None).normalize()).to_numpy()
    gW = par["g_mm_per_cm"] * w_nodes
    obs = q.y_vert_mm.to_numpy() - (gW[j] - gW[i])
    r, sd = f3.ou_increment_smoother(t_days, i, j, obs, q.sigma_vert_mm.to_numpy(), par["tau_days"], par["sigma_r_mm"])
    out = pd.DataFrame({"time_utc": nodes, "h_v3_mm": gW + r, "h_v3_sd_mm": sd, "prior_mm": gW})
    # baseline: the plain 6-day chain per track (cumulative vertical increments), as WP0 validated it
    for track, g in q.groupby("track"):
        g = g.sort_values("t1")
        run, (lv,) = f3.chain_runs(pd.DatetimeIndex(g.t1), pd.DatetimeIndex(g.t2), g.y_vert_mm.to_numpy(), min_len=2)
        ch = pd.Series(np.nan, index=nodes)
        cum = 0.0
        prev_t2 = None
        for t1_, t2_, y_ in zip(g.t1, g.t2, g.y_vert_mm):
            if prev_t2 is None or t1_ != prev_t2:
                cum = 0.0
                ch[t1_] = cum
            cum += y_
            ch[t2_] = cum
            prev_t2 = t2_
        out[f"chain_{track[:3]}_mm"] = ch.to_numpy()
    return out


def mat_raft(years: tuple[int, int], noise: pd.DataFrame) -> pd.DataFrame:
    """The mat's median phase per 6-day pair (hard-target reference), both tracks, as vertical increments."""
    fit = noise[(noise.rev_class == "≤ 7 d") & (noise.state == "dry")].iloc[0]
    rows = []
    for period in ("2020_2021", "2025", "2026"):
        for track in ("ascending", "descending"):
            f = G / f"mat_pairs_{period}_{track}.nc"
            if not f.exists():
                continue
            ds = xr.open_dataset(f)
            t1, t2 = pd.DatetimeIndex(ds.t1.values, tz="UTC"), pd.DatetimeIndex(ds.t2.values, tz="UTC")
            keep = ((ds.revisit_days.values <= 7) & (t1.year >= years[0]) & (t2.year <= years[1]) &
                    np.array(["S1D" not in x for x in ds.platforms.values.astype(str)]))
            if not keep.any():
                continue
            U = ds.unw.values[:, keep].T.astype(float)
            C = ds.coh.values[:, keep].T.astype(float)
            mat, hard = ds.st_zone.values == 1, ds.st_hard_target.values == 1
            ref = np.nanmedian(U[:, hard], 1)
            yl = np.nanmedian(-rad_to_los_mm(U[:, mat] - ref[:, None]), 1)
            losu = np.nanmedian(ds[f"st_los_u_{track[:3]}"].values[mat])
            coh = np.nanmedian(C[:, mat], 1)
            fz = [ds[c].values[keep] for c in ("drv_flag_frozen_air_t1", "drv_flag_frozen_air_t2")]
            frozen = (np.nan_to_num(fz[0]) > 0) | (np.nan_to_num(fz[1]) > 0)
            sig = f3.crb_sigma_mm(coh, fit.looks, fit.floor_mm)
            rows.append(pd.DataFrame({"track": track, "t1": t1[keep], "t2": t2[keep], "y_vert_mm": yl / losu,
                                      "sigma_vert_mm": np.where(frozen, np.nan, sig / losu), "coh": coh}))
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def run_raft(m: pd.DataFrame, W: pd.Series, par: dict) -> pd.DataFrame:
    nodes = pd.DatetimeIndex(sorted(set(m.t1) | set(m.t2)))
    t_days = ((nodes - nodes[0]).total_seconds() / 86400).to_numpy()
    idx = {t: k for k, t in enumerate(nodes)}
    i, j = m.t1.map(idx).to_numpy(), m.t2.map(idx).to_numpy()
    gW = par["g_mm_per_cm"] * W.reindex(nodes.tz_convert(None).normalize()).to_numpy()
    r, sd = f3.ou_increment_smoother(t_days, i, j, m.y_vert_mm.to_numpy() - (gW[j] - gW[i]), m.sigma_vert_mm.to_numpy(),
                                     par["tau_days"], par["sigma_r_mm"])
    return pd.DataFrame({"time_utc": nodes, "raft_mm": gW + r, "raft_sd_mm": sd})


# ------------------------------------------------------------------------------ validation

def laser_at(times: pd.DatetimeIndex) -> np.ndarray:
    h = pd.read_csv(G / "truth_laser_hourly.csv.gz", index_col=0, parse_dates=True)
    h.index = pd.DatetimeIndex(h.index, tz="UTC") if h.index.tz is None else h.index
    ok = h.laser_surface_cm.where(W2.truthy(h.laser_ok)) * 10
    from insar_wetlands.cube.core import value_at
    return value_at(ok, times)


def segments(times: pd.DatetimeIndex, valid: np.ndarray, max_gap_days: float = 15) -> np.ndarray:
    seg = np.full(len(times), -1)
    sid, prev = -1, None
    for k, (t, v) in enumerate(zip(times, valid)):
        if not v:
            continue
        if prev is None or (t - prev).days > max_gap_days:
            sid += 1
        seg[k] = sid
        prev = t
    return seg


def validate(series: pd.DataFrame, year: int, cols: list[str]) -> list[dict]:
    s = series[series.time_utc.dt.year == year].copy()
    s["laser_mm"] = laser_at(pd.DatetimeIndex(s.time_utc))
    s["seg"] = segments(pd.DatetimeIndex(s.time_utc), s.laser_mm.notna().to_numpy())
    rows = []
    for c in cols:
        if c not in s:
            continue
        g = s[(s.seg >= 0) & s[c].notna()].copy()
        g = g[g.groupby("seg").seg.transform("size") >= 4]
        if len(g) < 8:
            continue
        e = g[c] - g.groupby("seg")[c].transform("mean")
        t = g.laser_mm - g.groupby("seg").laser_mm.transform("mean")
        row = {"estimate": c, "year": year, "n_dates": len(g), "segments": int(g.seg.nunique()),
               "r": float(np.corrcoef(e, t)[0, 1]), "rmse_mm": float(np.sqrt(np.mean((e - t) ** 2))),
               "amplitude_ratio": float(np.std(e) / np.std(t))}
        sdc = c.replace("_mm", "_sd_mm")
        if sdc in g:
            # coverage of the level uncertainty: within-segment demeaning removes most of the shared offset
            cov = f3.coverage(((e - t) / g[sdc]).to_numpy())
            row.update({"within_1sd": cov["within_1sd"], "within_2sd": cov["within_2sd"]})
        rows.append(row)
    return rows


def depth_slope() -> float:
    t = pd.read_csv(D3 / "wp3" / "e31_shape_vs_raft.csv")
    t = t[(t["shape"] == "depth_in_mat") & (t.subset == "all mat pixels") & (t.grid == "40m") & (t.revisit_days == 6)]
    return float(t.slope.mean())


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
    d, p, noise = inputs()
    s_depth = depth_slope()
    st = xr.open_dataset(S / "static_40m.nc")
    mat = st.zone.values == 1
    depth = st.depth_in_mat_m.values
    z = (depth - np.nanmean(depth[mat])) / np.nanstd(depth[mat])
    rr, cc = np.nonzero(mat)
    z_p6 = float(z[68, 66])
    val, pars, all_series = [], [], []
    for name, years in RUNS.items():
        test_years = [y for y in range(years[0], years[1] + 1) if y in LASER_YEARS]
        # year-blocked for validation; for 2026 (no laser) parameters from every laser year
        blocks = test_years or [None]
        for ty in blocks:
          par = params(d, ty)
          dm, med = detect_model(p, ty)
          W = driver_series(d, par["g_mm_per_cm"])
          for vname, vscale in MISS_SCALES.items():
            q = observations(p, noise, years, dm, med, vscale)
            if q.empty:
                continue
            ser = run_p6(q, W, par)
            m = mat_raft(years, noise)
            if len(m):
                rf = run_raft(m, W, par)
                ser = ser.merge(rf, on="time_utc", how="left")
                ser["h_matmode_p6_mm"] = ser.raft_mm * (1 + s_depth * z_p6)
                ser["h_matmode_p6_sd_mm"] = ser.raft_sd_mm * abs(1 + s_depth * z_p6)
            ser.insert(0, "run", name)
            ser.insert(1, "test_year", ty if ty else "none")
            ser.insert(2, "variant", vname)
            all_series.append(ser)
            if vname == PRIMARY:
                pars.append({"run": name, "test_year": ty, **par, "pairs": len(q), "depth_slope": s_depth})
            if ty:
                cols = ["h_v3_mm", "h_matmode_p6_mm"] + (["chain_asc_mm", "chain_des_mm", "prior_mm"] if vname == PRIMARY else [])
                val += [{"run": name, "variant": vname if r["estimate"].startswith("h_") else "", **r}
                        for r in validate(ser, ty, cols)]
            # mat height map for the run (primary variant, parameters of the last block; 40 m)
            if len(m) and vname == PRIMARY and (ty is None or ty == blocks[-1]):
                rf_nodes = pd.DatetimeIndex(rf.time_utc)
                H = rf.raft_mm.to_numpy()[:, None] * (1 + s_depth * z[rr, cc])[None, :]
                Hs = rf.raft_sd_mm.to_numpy()[:, None] * np.abs(1 + s_depth * z[rr, cc])[None, :]
                grid = np.full((len(rf_nodes),) + mat.shape, np.nan, "float32")
                grid_sd = grid.copy()
                grid[:, rr, cc], grid_sd[:, rr, cc] = H, Hs
                xr.Dataset({"height_mm": (("time", "y", "x"), grid), "height_sd_mm": (("time", "y", "x"), grid_sd),
                            "raft_mm": ("time", rf.raft_mm.to_numpy()), "raft_sd_mm": ("time", rf.raft_sd_mm.to_numpy())},
                           coords={"time": rf_nodes.tz_convert(None), "x": st.x.values, "y": st.y.values},
                           attrs={"model": "h(x,t) = c(t)·(1 + s·z(x)); c = v3 raft (hard-target reference, 6-day pairs, "
                                           "noise × detectability, buoyancy prior); s = WP3 depth slope; z = standardised "
                                           "depth in the mat", "depth_slope": s_depth, "run": name,
                                  "units": "mm vertical, relative (levels arbitrary per run)"}) \
                    .to_netcdf(OUT / f"height_v3_{name.split()[0]}.nc")
    V = pd.DataFrame(val)
    P = pd.DataFrame(pars)
    A = pd.concat(all_series, ignore_index=True)
    V.to_csv(OUT / "validation_p6_laser.csv", index=False)
    P.to_csv(OUT / "parameters_by_block.csv", index=False)
    A.to_csv(OUT / "series_p6.csv", index=False)
    v2 = pd.read_csv(D3.parent / "fusion_v2" / "validation_summary.csv")
    v2 = v2[v2.estimate.isin(["v2 mat (P6 pixel)", "radar chain, ascending", "radar chain, descending", "prior only (g·W)"])]
    fig, ax = plt.subplots(2, 1, figsize=(11, 7))
    for k, (ty, a) in enumerate(zip((2021, 2025), ax)):
        s = A[(A.test_year == ty) & (A.variant == PRIMARY)]
        s = s[s.time_utc.dt.year == ty]
        lz = laser_at(pd.DatetimeIndex(s.time_utc))
        a.plot(s.time_utc, lz - np.nanmean(lz), "k.", ms=4, label="laser (QC-passed)")
        for c, col, lab in (("h_v3_mm", "C3", "v3 (P6)"), ("h_matmode_p6_mm", "C1", "v3 raft + edge gradient at P6"),
                            ("prior_mm", "0.6", "prior g·W")):
            v = s[c] - np.nanmean(s[c][np.isfinite(lz)])
            a.plot(s.time_utc, v, "-", color=col, lw=1, label=lab)
        if "h_v3_sd_mm" in s:
            v = s.h_v3_mm - np.nanmean(s.h_v3_mm[np.isfinite(lz)])
            a.fill_between(s.time_utc, v - s.h_v3_sd_mm, v + s.h_v3_sd_mm, color="C3", alpha=0.15)
        a.set_title(f"P6, {ty} held out (parameters from the other years)", fontsize=9)
        a.set_ylabel("vertical height (mm, relative)")
        a.legend(fontsize=7, ncol=4)
    fig.tight_layout()
    fig.savefig(OUT / "fig_p6_v3.png", dpi=150)
    plt.close(fig)
    txt = ["# X-078 — fusion v3 WP5: the estimator (exploratory, branch fusion-x062)", "",
           "Generated by `05_code/SAr-adapted/scripts/fusion_v3_wp5_x078.py` from `06_data/cube` and the WP0–WP4 outputs; "
           "every number is computed. Year-blocked: for each held-out laser year, g, τ, σ_r and the detectability model "
           "come from the other years. Levels demeaned within laser segments (≥ 4 dates, gaps ≤ 15 days).", "",
           "## Validation at P6 against the laser", "", md(V), "",
           "`h_v3_mm` = v3 at P6 (both tracks, one smoother); `h_matmode_p6_mm` = the mat's raft with the edge gradient, "
           "read at P6's depth (the map product); `chain_*` = the plain 6-day chain, same reference; `prior_mm` = g·W "
           "alone. Coverage = share of levels within ±1/±2 posterior SD.", "",
           "Two settings of the detectability inflation are reported, not one chosen silently: the first run treated a "
           "miss as a whole-cycle error (λ/2); WP4 showed 6-day misses are not cycle slips, so the primary setting uses "
           f"λ/4 (maps and figure: **{PRIMARY}**).", "",
           "Fusion v2 as published (its own reference and dates; `../../fusion_v2/validation_summary.csv`):", "", md(v2), "",
           "## Parameters by block", "", md(P), "", "![P6](fig_p6_v3.png)", "",
           "Mat height maps: `height_v3_R1.nc`, `height_v3_R3.nc`, `height_v3_R4.nc` (mm, vertical, relative; the edge "
           "gradient from WP3). 2026 has no laser: R4 is out of sample and unvalidated.", ""]
    (OUT / "README.md").write_text("\n".join(txt))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()

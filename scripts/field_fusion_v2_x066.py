"""X-066 — fusion v2: the floating mat's vertical height, every mat pixel and every acquisition of both tracks,
as one exact posterior (``insar_wetlands.fusion2``), judged by the P6 laser on blocks it was not fitted on.

Branch fusion-x062 (it chains per-pixel short phase changes; SETTLED §1 excludes that on main until Aymen
decides). No network inversion. Writes only to the hub.

Runs
  R1  2020–2021, consecutive 6-day pairs, ascending + descending (S1A + S1B; D-020 crops) — the primary map
  R2  2022–2024, consecutive 12-day pairs, ascending only (descending 12-day does not carry the motion, X-065)
Calibration at P6 (per validation fold, never on the held-out block)
  g   buoyancy: laser change vs P6 water-table change between consecutive acquisitions (field data only)
  τ, σ_c   OU time scale and size of the laser beyond buoyancy (daily laser − g·WTD, per laser segment)
  m, e, s  per track × revisit: s1 = m·cos θ·Δh_laser + e·ΔW, coherence-weighted (X-065), noise scale s
Mat model
  h(x,t) = g·W(t) + c(t) + d(x,t); W the mat plots' (P6–P9) water-table anomaly at the overpass; c shared
  (OU τ, σ_c); d departures (OU τ, σ_d, spatial coupling κ over zone A). σ_d and κ chosen by predicting
  the descending pairs from an ascending-only fit (no laser involved).
Validation
  V1  P6 level vs laser, held-out time blocks: v2 (mat run, P6 pixel), v2 single-pixel, prior only, radar chain
  V2  ascending-only vs descending-only: the same vertical motion?
  V3  stable ground (zone D block next to the mat, no buoyancy prior): how much "motion" is noise
Frozen dates (X-050) are dropped; wet dates get double noise.

    INSAR_DRIVE_ROOT=../local/drive_pristine PYTHONPATH=src python scripts/field_fusion_v2_x066.py
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

from insar_wetlands import field  # noqa: E402
from insar_wetlands import fusion as fu  # noqa: E402
from insar_wetlands import fusion2 as f2  # noqa: E402
from insar_wetlands.bootstrap import start  # noqa: E402
from insar_wetlands.inversion.isbas import PHASE_TO_MM  # noqa: E402
from insar_wetlands.paths import make_paths  # noqa: E402
from insar_wetlands.stack import list_pairs, load_layer  # noqa: E402

HUB = field.field_root().parents[1]
DLV = HUB / "08_deliverables"
INC = {"ascending": 32.26, "descending": 39.17}
HOUR = {"ascending": "16:36", "descending": "05:09"}
MAT = ["P6", "P7", "P8", "P9"]
WET_FACTOR = 2.0                       # noise multiplier on a pair with a wet date (method parameter)
MIN_LASER_PAIRS = 6                    # laser pairs needed to calibrate one track × revisit (two coefficients)
CACHE = HUB / "05_code" / "local" / "cache_fusion_v2"   # loaded pair cubes (local, not versioned): delete to reload


# ------------------------------------------------------------------------------ inputs

class Inputs:
    def __init__(self, ctx, extra: Path):
        self.ctx = ctx
        self.zones = {k: ctx.zones[k].values.astype(bool) for k in "ABCD"}
        px = pd.read_csv(DLV / "field_first" / "plot_pixels.csv").set_index("plot")
        self.p6 = (int(px.loc["P6", "row"]), int(px.loc["P6", "col"]))
        fl = pd.read_csv(DLV / "field_p6" / "laser_qc" / "laser_flags.csv", parse_dates=["time_utc"]).set_index("time_utc")
        fl["ok"] = fl.surface_cm.notna() & ~fl.snow_72h & ~fl.outlier & ~fl.filled
        self.laser = fl
        runs = pd.read_csv(DLV / "field_p6" / "laser_qc" / "interpolated.csv", parse_dates=["start", "end"])
        self.runs = runs
        self.wtd = field.load_wtd_hourly()
        self.cens = {p: field.censored_flag(self.wtd[p]) for p in field.PLOTS}
        w = pd.read_csv(DLV / "field_dew_x050" / "wetness_at_overpasses_2020_2024.csv", parse_dates=["date"])
        self.wet = w.set_index(["date", "track"])
        self.roots = {t: [make_paths(cfg=ctx.cfg, track=t).cropped, extra / f"hyp3_cropped_{t}"] for t in INC}
        # laser segments: a filled stretch may hide a re-levelling (X-058) — levels compare within a segment
        cuts = np.array(sorted(pd.to_datetime(runs.start, utc=True).dt.tz_convert(None)), dtype="datetime64[ns]")
        self.segment = lambda t: int(np.searchsorted(cuts, np.datetime64(pd.Timestamp(t).tz_convert(None) if pd.Timestamp(t).tzinfo else pd.Timestamp(t))))

    def when(self, d, track):
        return pd.Timestamp(f"{pd.Timestamp(d).date()} {HOUR[track]}", tz="UTC")

    def laser_at(self, t):
        i = self.laser.index.get_indexer([t], method="nearest")[0]
        if abs((self.laser.index[i] - t).total_seconds()) > 3600 or not self.laser.ok.iloc[i]:
            return np.nan
        return float(self.laser.surface_cm.iloc[i]) * 10            # mm

    def wtd_at(self, plot, t):
        return np.nan if bool(self.cens[plot].asof(t)) else field._interp_at(self.wtd[plot], t)

    def forcing(self, t):
        v = [self.wtd_at(p, t) for p in MAT]
        v = [x for x in v if np.isfinite(x)]
        return float(np.mean(v)) if len(v) >= 2 else np.nan

    def across_fill(self, ta, tb):
        return bool(((self.runs.end > ta) & (self.runs.start < tb)).any())

    def flags(self, d, track):
        try:
            r = self.wet.loc[(pd.Timestamp(d), track)]
            return bool(r.wet), bool(r.frozen)
        except KeyError:
            return False, True                                    # unknown → treat as frozen (dropped)


def consecutive_pairs(inp: Inputs, track: str, start: str, end: str) -> list[tuple[str, Path]]:
    have = {p: r for r in inp.roots[track] if r.is_dir() for p in list_pairs(r)}
    dates = sorted({pd.Timestamp(x) for p in have for x in p.split("_")})
    dates = [d for d in dates if pd.Timestamp(start) <= d <= pd.Timestamp(end)]
    out = []
    for a, b in zip(dates[:-1], dates[1:]):
        pid = f"{a:%Y%m%d}_{b:%Y%m%d}"
        if pid in have:
            out.append((pid, have[pid]))
    return out


def pair_cube(inp: Inputs, track: str, pairs) -> pd.DataFrame:
    """Per pair: LOS change cube (P − grassland median, mm), coherence cube, dates, Δt, flags, forcing, laser."""
    C = inp.zones["C"]
    rows, dl, co = [], [], []
    r6, c6 = inp.p6
    for pid, root in pairs:
        u = load_layer(root, "unw_phase", [pid]).isel(pair=0).values
        c = load_layer(root, "corr", [pid]).isel(pair=0).values
        d = (u - np.nanmedian(u[C & np.isfinite(u)])) * PHASE_TO_MM
        a, b = (pd.Timestamp(x) for x in pid.split("_"))
        ta, tb = inp.when(a, track), inp.when(b, track)
        wa, fa = inp.flags(a, track)
        wb, fb = inp.flags(b, track)
        la, lb = inp.laser_at(ta), inp.laser_at(tb)
        rows.append({"pair": pid, "track": track, "t1": a, "t2": b, "dt": (b - a).days, "wet": wa or wb, "frozen": fa or fb,
                     "W1": inp.forcing(ta), "W2": inp.forcing(tb),
                     "p6_w1": inp.wtd_at("P6", ta), "p6_w2": inp.wtd_at("P6", tb),
                     "laser1": la, "laser2": lb, "seg1": inp.segment(ta), "seg2": inp.segment(tb),
                     "across_fill": inp.across_fill(ta, tb), "s1_p6": d[r6, c6], "coh_p6": c[r6, c6]})
        dl.append(d.astype("float32"))
        co.append(c.astype("float32"))
    df = pd.DataFrame(rows)
    df["dW"] = df.W2 - df.W1
    df["dh_laser"] = np.where(df.across_fill | (df.seg1 != df.seg2), np.nan, df.laser2 - df.laser1)
    return df, np.stack(dl), np.stack(co)


# ------------------------------------------------------------------------------ calibration at P6

def calibrate(inp: Inputs, tables: dict, train) -> dict:
    """g, τ, σ_c from the laser and water table (field data only); m, e, s per track × revisit from the phase.
    ``train(t)`` says whether a time belongs to the calibration block."""
    allp = pd.concat(tables.values())
    lp = allp[allp.dh_laser.notna() & (allp.p6_w2 - allp.p6_w1).notna() & allp.t2.map(train)]
    g = float(np.polyfit(lp.p6_w2 - lp.p6_w1, lp.dh_laser, 1)[0])            # mm per cm
    # daily laser beyond buoyancy, per segment, in the calibration block
    lz = inp.laser[inp.laser.ok].surface_cm * 10
    daily = lz.resample("D").mean().dropna()
    p6 = inp.wtd["P6"].where(~inp.cens["P6"]).resample("D").mean()
    p6.index = p6.index.tz_localize("UTC") if p6.index.tz is None else p6.index
    res = (daily - g * p6.reindex(daily.index)).dropna()
    seg = np.array([inp.segment(t) for t in res.index])
    res = res - pd.Series(res.values, index=res.index).groupby(seg).transform("mean").values
    keep = np.array([train(t.tz_convert(None)) for t in res.index])
    ou = f2.ou_from_series((res.index[keep] - res.index[0]).days.to_numpy(float), res.values[keep], max_lag_days=60)
    phys = {}
    for (track, dt), d in allp.groupby(["track", "dt"]):
        d = d[d.dh_laser.notna() & d.dW.notna() & ~d.frozen & d.t2.map(train)]
        if len(d) < MIN_LASER_PAIRS:
            continue
        X = np.column_stack([d.dh_laser * np.cos(np.radians(INC[track])), d.dW])
        f = fu.bayes_linear(X, d.s1_p6.to_numpy(), fu.phase_sigma_mm(d.coh_p6.to_numpy()))
        phys[(track, int(dt))] = {"m": float(f["beta"][0]), "m_sd": float(f["sd"][0]), "e": float(f["beta"][1]),
                                  "s": f["noise_scale"], "n": len(d)}
    return {"g": g, "tau": ou["tau_days"], "sigma_c": ou["sigma"], "phys": phys, "n_buoyancy": len(lp)}


# ------------------------------------------------------------------------------ the model on a mask

def stable_noise(tables, cubes, zones) -> dict:
    """Per track × revisit × wet/dry: the median over stable-ground pixels near the mat (zone D, 0.2–1.5 km from it)
    of the SD of their pair-to-pair LOS changes. Stable ground does not move, so this is the radar's real noise
    (mostly atmosphere: about five times the coherence bound, and nearly independent from one pair to the next).
    Found 2026-10-01: the coherence bound × a noise scale fitted on 13 laser pairs gave about half of it."""
    from scipy import ndimage
    dist = ndimage.distance_transform_edt(~zones["A"]) * 40.0
    near = zones["D"] & (dist > 200) & (dist < 1500)
    out = {}
    for track, d in tables.items():
        dl = cubes[track][0]
        for (dt, wet), g in d[~d.frozen].groupby(["dt", "wet"]):
            if len(g) >= 8:
                out[(track, int(dt), bool(wet))] = float(np.nanmedian(np.nanstd(dl[g.index.to_numpy()][:, near], axis=0)))
    return out


def build_pairs(tables, cubes, mask, cal, T_index, prior_W=True, tracks=None):
    """fusion2.Pairs over the True pixels of ``mask`` for every usable pair (track × revisit calibrated)."""
    idx = np.flatnonzero(mask.ravel())
    out = []
    for track, df in tables.items():
        if tracks and track not in tracks:
            continue
        dl, co = cubes[track]
        for k, r in df.iterrows():
            ph = cal["phys"].get((track, int(r["dt"])))
            if ph is None or r.frozen or not np.isfinite(r.dW) or ph["m"] <= 0.1:
                continue
            y = dl[k].ravel()[idx] - (ph["e"] * r.dW if prior_W else 0.0)
            emp = cal["noise"].get((track, int(r["dt"]), bool(r.wet)), cal["noise"].get((track, int(r["dt"]), False)))
            crb = fu.phase_sigma_mm(co[k].ravel()[idx])
            sd = np.sqrt(emp ** 2 + crb ** 2) if emp is not None else crb * ph["s"] * (WET_FACTOR if r.wet else 1.0)
            ok = np.isfinite(y) & np.isfinite(sd)
            n = int(ok.sum())
            out.append(f2.Pairs(np.flatnonzero(ok), np.full(n, T_index[(track, r.t1)]), np.full(n, T_index[(track, r.t2)]),
                                np.full(n, ph["m"] * np.cos(np.radians(INC[track]))), y[ok], sd[ok]))
    return f2.Pairs(*(np.concatenate([getattr(p, f) for p in out]) for f in ("px", "i", "j", "A", "y", "sd"))), idx


def run_mask(tables, cubes, mask, cal, T, T_index, W_T, sigma_d_ratio, kappa, buoyant=True, tracks=None):
    p, idx = build_pairs(tables, cubes, mask, cal, T_index, prior_W=buoyant, tracks=tracks)
    n_px, n_t = len(idx), len(T)
    tdays = (pd.DatetimeIndex(T) - pd.Timestamp(T[0])).days.to_numpy(float)
    L, _ = f2.grid_laplacian(mask)
    Qc = f2.ou_precision(tdays, cal["tau"], cal["sigma_c"])
    Qd = f2.space_time_precision(f2.ou_precision(tdays, cal["tau"], sigma_d_ratio * cal["sigma_c"]), n_px, L, kappa)
    mu = np.tile((cal["g"] * np.nan_to_num(W_T - np.nanmean(W_T))) if buoyant else np.zeros(n_t), (n_px, 1))
    # ambiguities pixel by pixel (cheap direct solves), then the mat-wide system by conjugate gradients
    Qsingle = f2.ou_precision(tdays, cal["tau"], cal["sigma_c"] * np.sqrt(1 + sigma_d_ratio ** 2))
    k = f2.resolve_ambiguities_per_pixel(p, mu, Qsingle, n_t)
    p = f2.Pairs(p.px, p.i, p.j, p.A, p.y + k * f2.CYCLE_LOS_MM, p.sd)
    out = f2.smooth(p, mu, Qd, n_t, n_px, Qcommon=Qc, ambiguity_iters=0, solver="cg")
    out.update({"idx": idx, "pairs": p, "k": k})
    return out


def joint_axis(tables, tracks):
    keys = sorted({(t, d) for t, df in tables.items() if t in tracks for d in pd.concat([df.t1, df.t2])},
                  key=lambda k: pd.Timestamp(f"{k[1].date()} {HOUR[k[0]]}"))
    T = [pd.Timestamp(f"{d.date()} {HOUR[t]}") for t, d in keys]
    return T, {k: i for i, k in enumerate(keys)}


# ------------------------------------------------------------------------------ validation helpers

def laser_level_metrics(inp, T, h, name):
    """Level of an estimate at the dates with a laser value, against the laser, within laser segments."""
    rows = []
    for i, t in enumerate(T):
        tt = pd.Timestamp(t).tz_localize("UTC")
        lv = inp.laser_at(tt)
        if np.isfinite(lv) and np.isfinite(h[i]):
            rows.append({"t": t, "laser": lv, "est": h[i], "seg": inp.segment(tt)})
    d = pd.DataFrame(rows)
    if len(d) < 6:
        return {"estimate": name, "n": len(d)}, d
    for c in ("laser", "est"):
        d[c + "_c"] = d[c] - d.groupby("seg")[c].transform("mean")
    d = d[d.groupby("seg").t.transform("size") >= 3]
    r = float(np.corrcoef(d.laser_c, d.est_c)[0, 1])
    return {"estimate": name, "n": len(d), "r": r, "rmse_mm": float(np.sqrt(np.mean((d.est_c - d.laser_c) ** 2))),
            "amplitude_ratio": float(d.est_c.std() / d.laser_c.std())}, d


def radar_chain(table, cal, track, T, T_index):
    """The P6 pixel's increments divided by m·cos θ, chained — the radar alone, per track."""
    h = np.full(len(T), np.nan)
    df = table.sort_values("t1")
    level = 0.0
    for _, r in df.iterrows():
        ph = cal["phys"].get((track, int(r["dt"])))
        i, j = T_index[(track, r.t1)], T_index[(track, r.t2)]
        if np.isnan(level):
            level = 0.0                     # a new segment after a break (levels compare within laser segments)
        if np.isnan(h[i]):
            h[i] = level
        if ph is None or r.frozen or ph["m"] <= 0.1:
            level = np.nan
            continue
        level = h[i] + (r.s1_p6 - ph["e"] * np.nan_to_num(r.dW)) / (ph["m"] * np.cos(np.radians(INC[track])))
        h[j] = level
    return h


# ------------------------------------------------------------------------------ main

def md(df):
    f = lambda v: f"{v:.2f}" if isinstance(v, float) else str(v)  # noqa: E731
    return "\n".join(["| " + " | ".join(map(str, df.columns)) + " |", "|" + "---|" * len(df.columns)]
                     + ["| " + " | ".join(f(v) for v in r) + " |" for r in df.itertuples(index=False)])


def run(inp, name, period, tracks, folds, out):
    print(f"\n== {name}: {period} {tracks}")
    tables, cubes = {}, {}
    for t in tracks:
        cache = CACHE / f"{name}_{t}.pkl"
        if cache.exists():
            tables[t], dl, co = pd.read_pickle(cache)
        else:
            pairs = consecutive_pairs(inp, t, *period)
            tables[t], dl, co = pair_cube(inp, t, pairs)
            CACHE.mkdir(parents=True, exist_ok=True)
            pd.to_pickle((tables[t], dl, co), cache)
        cubes[t] = (dl, co)
        print(f"  {t}: {len(tables[t])} consecutive pairs, Δt {tables[t].dt.value_counts().to_dict()}", flush=True)
    T, T_index = joint_axis(tables, tracks)
    W_T = np.array([inp.forcing(pd.Timestamp(t).tz_localize("UTC")) for t in T])
    A = inp.zones["A"]
    r6, c6 = inp.p6
    p6_flat = r6 * A.shape[1] + c6
    if callable(folds):
        folds = folds(tables)
    noise = stable_noise(tables, cubes, inp.zones)
    print("  stable-ground noise (mm per pair):", {f"{k[0][:3]} {k[1]}d {'wet' if k[2] else 'dry'}": round(v, 2) for k, v in noise.items()}, flush=True)
    metrics, choice_rows, cal_rows = [], [], []
    series = {}
    for fold_name, train, test in folds:
        cal = calibrate(inp, tables, train)
        cal["noise"] = noise
        cal_rows.append({"fold": fold_name, "g_mm_per_cm": cal["g"], "tau_days": cal["tau"], "sigma_c_mm": cal["sigma_c"],
                         **{f"m_{t[:3]}_{dt}d": v["m"] for (t, dt), v in cal["phys"].items()},
                         **{f"e_{t[:3]}_{dt}d": v["e"] for (t, dt), v in cal["phys"].items()}})
        # σ_d, κ: predict descending pairs from an ascending-only fit (both tracks needed)
        best = (0.5, 16.0)
        if len(tracks) == 2 and any(k[0] == "descending" for k in cal["phys"]):
            pd_desc, idx = build_pairs(tables, cubes, A, cal, T_index, tracks=["descending"])
            scores = []
            for sdr in (0.25, 0.5):
                for kap in (4.0, 16.0):
                    o = run_mask(tables, cubes, A, cal, T, T_index, W_T, sdr, kap, tracks=["ascending"])
                    hx = o["h"]
                    pred = pd_desc.A * (hx[pd_desc.px, pd_desc.j] - hx[pd_desc.px, pd_desc.i])
                    sc = float(np.sqrt(np.mean(((pd_desc.y - pred) / pd_desc.sd) ** 2)))
                    scores.append((sc, sdr, kap))
                    choice_rows.append({"fold": fold_name, "sigma_d_ratio": sdr, "kappa": kap, "desc_pred_rms_z": sc})
            best = min(scores)[1:]
        o = run_mask(tables, cubes, A, cal, T, T_index, W_T, *best)
        k6 = int(np.flatnonzero(o["idx"] == p6_flat)[0])
        h_mat = o["h"][k6]
        single = run_mask(tables, cubes, np.eye(1, A.size, p6_flat, dtype=bool).reshape(A.shape), cal, T, T_index, W_T, 1.0, 0.0)
        prior = cal["g"] * np.nan_to_num(W_T - np.nanmean(W_T))
        ests = {"v2 mat (P6 pixel)": h_mat, "v2 single pixel": single["h"][0], "prior only (g·W)": prior}
        for t in tracks:
            ests[f"radar chain, {t}"] = radar_chain(tables[t], cal, t, T, T_index)
        held = np.array([test(pd.Timestamp(t)) for t in T])
        for en, h in ests.items():
            m, dd = laser_level_metrics(inp, [t for t, k in zip(T, held) if k], np.asarray(h)[held], en)
            metrics.append({"run": name, "fold": fold_name, **m})
        series[fold_name] = {"T": T, "ests": ests, "held": held, "cal": cal, "best": best, "out": o}
    # final model: calibrated on everything, for the maps and the track / stable-ground checks
    cal = calibrate(inp, tables, lambda t: True)
    cal["noise"] = noise
    best = series[folds[0][0]]["best"]
    full = run_mask(tables, cubes, A, cal, T, T_index, W_T, *best)
    checks = []
    if len(tracks) == 2:
        ha = run_mask(tables, cubes, A, cal, T, T_index, W_T, *best, tracks=["ascending"])["h"]
        hd = run_mask(tables, cubes, A, cal, T, T_index, W_T, *best, tracks=["descending"])["h"]
        da, dd = ha - ha.mean(axis=1, keepdims=True), hd - hd.mean(axis=1, keepdims=True)
        pix_r = [np.corrcoef(da[k], dd[k])[0, 1] for k in range(len(da))]
        checks.append({"run": name, "check": "V2 asc-only vs desc-only, per mat pixel: median r of the height series",
                       "value": float(np.nanmedian(pix_r))})
        # the same without the shared prior (g·W): only what each track's radar adds — the fair test
        prior = cal["g"] * np.nan_to_num(W_T - np.nanmean(W_T))
        ra, rd = da - (prior - prior.mean()), dd - (prior - prior.mean())
        pix_r2 = [np.corrcoef(ra[k], rd[k])[0, 1] for k in range(len(ra))]
        checks.append({"run": name, "check": "V2 asc-only vs desc-only, per mat pixel: median r of height minus the water prior",
                       "value": float(np.nanmedian(pix_r2))})
        checks.append({"run": name, "check": "V2 mat mean, height minus the water prior: r(asc-only, desc-only)",
                       "value": float(np.corrcoef(ra.mean(0), rd.mean(0))[0, 1])})
    # V3 stable ground: a block of zone D east of the mat, no buoyancy, same noise model
    rows_, cols_ = np.nonzero(A)
    block = np.zeros_like(A)
    r0, c0 = int(rows_.mean()), int(cols_.max()) + 8
    block[r0 - 10:r0 + 10, c0:c0 + 20] = True
    block &= inp.zones["D"]
    st = run_mask(tables, cubes, block, cal, T, T_index, W_T, *best, buoyant=False)
    prior_T = cal["g"] * np.nan_to_num(W_T - np.nanmean(W_T))
    sd_mat = float(np.median(np.std(full["h"] - prior_T[None, :], axis=1)))
    sd_st = float(np.median(np.std(st["h"], axis=1)))
    checks += [{"run": name, "check": "V3 median SD of the radar-driven height (minus the water prior), mat (mm)", "value": sd_mat},
               {"run": name, "check": "V3 median SD of the height, stable ground block, same model without buoyancy (mm)", "value": sd_st},
               {"run": name, "check": "V3 median SD of the full height incl. the water prior, mat (mm)",
                "value": float(np.median(np.std(full["h"], axis=1)))},
               {"run": name, "check": "phase ambiguities kept (cycles), mat run", "value": float(np.count_nonzero(full["k"]))},
               {"run": name, "check": "pair observations, mat run", "value": float(len(full["k"]))}]
    # maps
    grid = np.full((len(T),) + A.shape, np.nan, "float32")
    grid.reshape(len(T), -1)[:, full["idx"]] = full["h"].T
    ds = xr.Dataset({"height_mm": (("time", "y", "x"), grid)},
                    coords={"time": pd.DatetimeIndex(T), "y": inp.ctx.template.y.values, "x": inp.ctx.template.x.values})
    ds["shared_motion_mm"] = ("time", (full["c"] + cal["g"] * np.nan_to_num(W_T - np.nanmean(W_T))).astype("float32"))
    ds.attrs.update({"run": name, "g_mm_per_cm": cal["g"], "tau_days": cal["tau"], "sigma_c_mm": cal["sigma_c"],
                     "sigma_d_ratio": best[0], "kappa": best[1], "units": "mm, vertical, relative to the series mean"})
    ds.to_netcdf(out / f"height_{name}.nc")
    return {"metrics": metrics, "choice": choice_rows, "cal": cal_rows, "checks": checks, "series": series,
            "final_cal": cal, "best": best, "T": T, "full": full, "inp": inp}


def figure(out, res, name):
    s = res["series"]
    fig, axes = plt.subplots(len(s), 1, figsize=(12, 3.4 * len(s)), squeeze=False)
    for ax, (fold, v) in zip(axes[:, 0], s.items()):
        T = pd.DatetimeIndex(v["T"])
        held = v["held"]
        lz = [res["inp"].laser_at(pd.Timestamp(t).tz_localize("UTC")) for t in T]
        lz = np.array(lz)
        ax.plot(T[held], lz[held] - np.nanmean(lz[held]), "k.", ms=4, label="laser (held out)")
        for en, c in (("v2 mat (P6 pixel)", "C3"), ("prior only (g·W)", "C2"), ("v2 single pixel", "C1")):
            h = np.asarray(v["ests"][en])
            ax.plot(T[held], h[held] - np.nanmean(h[held]), "-", color=c, lw=1.2, label=en)
        ax.set_ylabel("vertical (mm)")
        ax.set_title(f"{name} — fold '{fold}': P6, held-out block", fontsize=9)
        ax.legend(fontsize=7, ncol=4)
    fig.tight_layout()
    fig.savefig(out / f"fig_p6_{name}.png", dpi=130)
    plt.close(fig)


def readme(out, summ, checks, cal, R1, R2):
    x62 = pd.read_csv(DLV / "field_fusion_x062" / "p6_levels.csv")
    L = ["# X-066 — fusion v2: vertical height of the floating mat (exploratory, branch fusion-x062)", "",
         "Generated by `05_code/SAr-adapted/scripts/field_fusion_v2_x066.py`; every number below is computed by it. "
         "Plan and data inventory: `PLAN.md`. Model: `insar_wetlands.fusion2` (tested on synthetic ground truth).", "",
         "## Model", "",
         "h(x,t) = g·W(t) + c(t) + d(x,t) — buoyancy × the mat plots' water-table anomaly, a motion the whole mat shares "
         "(Ornstein–Uhlenbeck, mean-reverting) and smooth per-pixel departures; observed through consecutive short pairs, "
         "s = m·cos θ·Δh + e·ΔW + noise; the noise per pair is the measured spread of stable ground near the mat for that track, revisit and wet/dry state (with the coherence bound added), frozen dates out. Ambiguities decided pixel by "
         "pixel by model comparison; the mat system solved exactly (preconditioned conjugate gradients). No network inversion.", "",
         "- **R1** 2020–2021, consecutive 6-day pairs, ascending + descending (S1A + S1B).",
         "- **R2** 2022–2024, consecutive 12-day pairs, ascending only.", "",
         "## Calibration (P6; per fold, never on the held-out block)", "", md(cal.round(3)), "",
         "g from laser vs water-table changes; τ and σ_c from the daily laser beyond buoyancy (field data only); m, e per "
         "track × revisit from the phase (X-065).", "",
         "## Validation at P6 against the laser, held-out blocks (levels within laser segments)", "", md(summ.round(3)), "",
         "`amplitude_ratio` = SD of the estimate / SD of the laser on the held-out dates (1 = full motion).", "",
         "For reference, X-062 (2022–2024, ascending, level over the whole record — not held out):", "",
         md(x62.round(3)), "",
         "## Checks", "", md(checks.round(3)), "",
         "- V2: ascending-only and descending-only fits see the same vertical motion if their height series agree; the "
         "fair version removes the water prior both share.",
         "- V3: the same model on a stable-ground block (no buoyancy) gives the noise level any mat motion must exceed.", "",
         "## Files", "", "| File | Content |", "|---|---|",
         "| `height_R1_2020_2021_6day_both.nc`, `height_R2_2022_2024_12day_asc.nc` | height per mat pixel and date (mm, vertical), shared motion |",
         "| `validation_p6_laser.csv`, `validation_summary.csv` | V1 per fold and summary |",
         "| `calibration_by_fold.csv`, `final_calibration.json`, `hyperparameter_choice.csv` | parameters |",
         "| `checks.csv` | V2, V3, ambiguities |", "| `fig_p6_*.png` | P6 against the laser on the held-out blocks |"]
    (out / "README.md").write_text("\n".join(L) + "\n")


def median_split_folds(tables):
    """Two folds split at the median date of the usable laser pairs (the laser covers only May–Dec 2021 here),
    so both blocks can calibrate and both can be tested."""
    d = pd.concat(tables.values())
    split = pd.Timestamp(d[d.dh_laser.notna() & ~d.frozen].t2.median()).normalize()
    tag = f"{split:%Y-%m-%d}"
    return [(f"calibrate before {tag}, test after", lambda t: t < split, lambda t: t >= split),
            (f"calibrate from {tag}, test before", lambda t: t >= split, lambda t: t < split)]


def year_fold(y: int):
    """(name, train, test) holding out calendar year ``y``."""
    return f"year {y} held out", lambda t: t.year != y, lambda t: t.year == y


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--extra", default=str(HUB / "05_code" / "local" / "s1_2020_2021"))
    ap.add_argument("--out", default=str(DLV / "fusion_v2"))
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    ctx = start("field_fusion_v2_x066", mount=False, git=False)
    inp = Inputs(ctx, Path(a.extra))
    R1 = run(inp, "R1_2020_2021_6day_both", ("2020-01-01", "2022-01-31"), ("ascending", "descending"), median_split_folds, out)
    yrs = (2022, 2023, 2024)
    R2 = run(inp, "R2_2022_2024_12day_asc", ("2022-01-01", "2024-12-31"), ("ascending",),
             [year_fold(y) for y in yrs], out)
    for R, n in ((R1, "R1_2020_2021_6day_both"), (R2, "R2_2022_2024_12day_asc")):
        figure(out, R, n)
    met = pd.DataFrame(R1["metrics"] + R2["metrics"])
    met.to_csv(out / "validation_p6_laser.csv", index=False)
    pd.DataFrame(R1["checks"] + R2["checks"]).to_csv(out / "checks.csv", index=False)
    pd.DataFrame(R1["cal"] + R2["cal"]).to_csv(out / "calibration_by_fold.csv", index=False)
    pd.DataFrame(R1["choice"]).to_csv(out / "hyperparameter_choice.csv", index=False)
    summ = met.groupby(["run", "estimate"]).agg(n=("n", "sum"), r=("r", "mean"), rmse_mm=("rmse_mm", "mean"),
                                                amplitude_ratio=("amplitude_ratio", "mean")).reset_index()
    summ.to_csv(out / "validation_summary.csv", index=False)
    (out / "final_calibration.json").write_text(json.dumps(
        {n: {"g_mm_per_cm": R["final_cal"]["g"], "tau_days": R["final_cal"]["tau"], "sigma_c_mm": R["final_cal"]["sigma_c"],
             "sigma_d_ratio": R["best"][0], "kappa": R["best"][1],
             "phys": {f"{t}|{dt}": v for (t, dt), v in R["final_cal"]["phys"].items()}} for n, R in (("R1", R1), ("R2", R2))}, indent=2))
    readme(out, summ, pd.DataFrame(R1["checks"] + R2["checks"]), pd.DataFrame(R1["cal"] + R2["cal"]), R1, R2)
    pd.set_option("display.width", 220)
    print(summ.round(3).to_string(index=False))
    print(pd.DataFrame(R1["checks"] + R2["checks"]).round(3).to_string(index=False))


if __name__ == "__main__":
    main()

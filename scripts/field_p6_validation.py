"""X-059 — P6/CR, the primary validation site: laser ↔ WTD ↔ Sentinel-1 phase ↔ coherence.

The supervisor's request of 2026-09-28: at P6 (the floating-mat plot with the laser), put the
laser surface, the water table, the Sentinel-1 apparent displacement and the coherence on one
graph for 2022–2024; quantify, **between consecutive Sentinel-1 acquisitions**, the laser
displacement, the WTD change, the apparent Sentinel-1 displacement and the coherence of the same
interferogram; and give the laser / Sentinel-1 ratio once both are in comparable units. He
adds a caution: the magnitude mismatch alone shows the phase is not the full motion, not that
the rest is dielectric — that needs the water/moisture relations with the lake and stable-ground
controls. Those controls are built here at the pixel scale of P6.

Conventions (X-058, ``08_deliverables/field_p6/laser_qc``): laser in cm, increase = up, as
delivered (no cos 10°), nearest **measured** hourly value within ±1 h of the overpass, snow-free
(72 h rule), outliers and filled stretches excluded. Sentinel-1: unwrapped phase of the P6
window minus the grassland (zone C) median, × −λ/4π → mm LOS, positive toward the satellite
(= up). Laser in LOS: Δh · cos(incidence), i.e. what a purely vertical motion gives.

Inputs: the first deliverable's plot extraction (``field_first``), the X-058 laser table, the
X-050 wetness flags, and — for the controls only — the cropped interferograms (Drive snapshot,
``INSAR_DRIVE_ROOT``; Mac-safe). Writes ONLY into the hub (``--out``).

    INSAR_DRIVE_ROOT=../local/drive_pristine PYTHONPATH=src python scripts/field_p6_validation.py
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import binomtest  # noqa: E402

from insar_wetlands import field  # noqa: E402
from insar_wetlands import field_link as fl  # noqa: E402
from insar_wetlands.bootstrap import start  # noqa: E402
from insar_wetlands.inversion.isbas import PHASE_TO_MM, WAVELENGTH_M  # noqa: E402
from insar_wetlands.stack import load_layer  # noqa: E402

HUB = Path(__file__).resolve().parents[3]
DLV = HUB / "08_deliverables"
REPO = Path(__file__).resolve().parents[1]
INCIDENCE = {"ascending": 32.26, "descending": 39.17}       # config.yaml, measured (X-042)
QUARTER_WAVE_MM = WAVELENGTH_M / 4 * 1000
WINDOWS = ("1x1", "3x3")
PRIMARY = "1x1"
PERIOD = ("2022-01-01", "2024-12-31")
COLOR = {"ascending": "tab:orange", "descending": "tab:blue"}
RNG = np.random.default_rng(59)


# ------------------------------------------------------------------------------ inputs

def load_inputs():
    acq = pd.read_csv(DLV / "field_first" / "s1_acquisitions.csv")
    acq["time_utc"] = pd.to_datetime(acq.date + " " + acq.overpass_utc).dt.tz_localize("UTC")
    by_pair = pd.read_csv(DLV / "field_first" / "s1_plot_by_pair.csv")
    p6 = by_pair[(by_pair["plot"] == "P6") & by_pair.window.isin(WINDOWS)]
    ov = pd.read_csv(DLV / "field_p6" / "laser_qc" / "laser_at_s1.csv")
    qc = pd.read_csv(DLV / "field_p6" / "laser_qc" / "summary.csv", index_col=0)["value"]
    interp = pd.read_csv(DLV / "field_p6" / "laser_qc" / "interpolated.csv", parse_dates=["start", "end"])
    wtd_s1 = pd.read_csv(DLV / "field_first" / "wtd_at_s1.csv")
    wtd_s1 = wtd_s1[wtd_s1["plot"] == "P6"].set_index(["date", "track"])
    wet = pd.read_csv(DLV / "field_dew_x050" / "wetness_at_overpasses_2020_2024.csv").set_index(["date", "track"])
    zc = pd.concat([pd.read_csv(REPO / "results" / "tables" / f).assign(track=t) for f, t in
                    (("phaseD_coh_by_zone.csv", "ascending"), ("phaseD_coh_by_zone_descending.csv", "descending"))])
    return acq, p6, ov, qc, interp, wtd_s1, wet, zc


def filled_runs(interp: pd.DataFrame) -> pd.DataFrame:
    """The filled stretches of the laser (X-058), naive UTC like the pair dates. A level may have
    been re-referenced inside any of them, so no laser change is taken across one."""
    return interp.assign(start=interp.start.dt.tz_convert(None), end=interp.end.dt.tz_convert(None))


def laser_segment(times, runs: pd.DataFrame) -> np.ndarray:
    """Segment id of each time: the number of filled stretches that end before it."""
    tt = pd.to_datetime(pd.Series(times))
    if tt.dt.tz is not None:
        tt = tt.dt.tz_convert(None)
    return np.array([(runs.end < x).sum() for x in tt])


# ------------------------------------------------------------------------------ the pair table

def pair_table(acq, p6, ov, wtd_s1, wet, wtd_hourly, runs) -> pd.DataFrame:
    """One row per consecutive pair × track × window: laser, WTD, CR level, Sentinel-1, flags."""
    ovi = ov.set_index(["date", "track"])
    rows = []
    for track, g in acq.groupby("track"):
        times = dict(zip(pd.to_datetime(g.date), g.time_utc))
        for t1, t2 in fl.consecutive_pairs(g.date):
            d1, d2 = f"{t1:%Y-%m-%d}", f"{t2:%Y-%m-%d}"
            pid = f"{t1:%Y%m%d}_{t2:%Y%m%d}"
            l1, l2 = ovi.loc[(d1, track)], ovi.loc[(d2, track)]
            w1, w2 = wtd_s1.loc[(d1, track)], wtd_s1.loc[(d2, track)]
            cr1, cr2 = (field._interp_at(wtd_hourly["CR_raw"], times[t]) for t in (t1, t2))
            usable = bool(l1.usable and l2.usable)
            dh = (l2.laser_cm - l1.laser_cm) * 10 if usable else np.nan
            rec = {"track": track, "pair": pid, "t1": d1, "t2": d2, "dt_days": (t2 - t1).days,
                   "mid": t1 + (t2 - t1) / 2, "laser_usable": usable, "dh_mm": dh,
                   "laser_dlos_mm": dh * np.cos(np.radians(INCIDENCE[track])),
                   "crosses_filled": bool(((runs.end > t1) & (runs.start < t2)).any()),
                   "dwtd_cm": w2.wtd_at - w1.wtd_at, "wtd_censored": bool(w1.wtd_censored or w2.wtd_censored),
                   "dcr_cm": cr2 - cr1,
                   "wet_either": bool(wet.loc[(d1, track)].wet or wet.loc[(d2, track)].wet),
                   "frozen_either": bool(wet.loc[(d1, track)].frozen or wet.loc[(d2, track)].frozen)}
            for w in WINDOWS:
                s = p6[(p6.pair == pid) & (p6.track == track) & (p6.window == w)]
                ok = len(s) == 1
                rows.append({**rec, "window": w, "in_network": ok,
                             "n_valid": int(s.n_valid.iloc[0]) if ok else 0,
                             "coh": float(s.coh_median.iloc[0]) if ok else np.nan,
                             "dphase_rad": float(s.dphase_vs_C_rad.iloc[0]) if ok else np.nan,
                             "s1_dlos_mm": float(s.dlos_vs_C_mm.iloc[0]) if ok else np.nan})
    df = pd.DataFrame(rows)
    df.loc[df.crosses_filled, ["dh_mm", "laser_dlos_mm"]] = np.nan
    df.loc[df.crosses_filled, "laser_usable"] = False
    df["laser_beyond_quarter_wave"] = df.laser_dlos_mm.abs() > QUARTER_WAVE_MM
    return df


# ------------------------------------------------------------------------------ statistics

def coh_residual(d: pd.DataFrame) -> np.ndarray:
    return fl.season_residual(d.coh.to_numpy(), d.mid, d.dt_days.to_numpy())


def pair_stats(d: pd.DataFrame, floor_mm: float) -> dict:
    """Every quantity the supervisor asked for, on one set of pairs (time-ordered)."""
    d = d.sort_values("mid")
    L = d[d.laser_usable & d.s1_dlos_mm.notna() & d.laser_dlos_mm.notna()]
    W = d[d.s1_dlos_mm.notna() & d.dwtd_cm.notna() & ~d.wtd_censored]
    out = {"n_pairs": len(d), "n_in_network": int(d.in_network.sum()), "n_laser": len(L), "n_wtd": len(W)}
    if len(L) >= 8:
        exp, obs = L.laser_dlos_mm.to_numpy(), L.s1_dlos_mm.to_numpy()
        sl = fl.slopes(exp, obs)
        r, p = fl.circular_shift_p(exp, obs)
        big = np.abs(exp) > floor_mm
        mid = big & (np.abs(exp) < QUARTER_WAVE_MM)
        agree = int(np.sum(np.sign(exp[big]) == np.sign(obs[big])))
        out.update({
            "median_abs_dh_mm": float(np.median(np.abs(L.dh_mm))), "max_abs_dh_mm": float(np.max(np.abs(L.dh_mm))),
            "median_abs_laser_dlos_mm": float(np.median(np.abs(exp))), "median_abs_s1_dlos_mm": float(np.median(np.abs(obs))),
            "share_laser_beyond_quarter_wave": float(np.mean(np.abs(exp) > QUARTER_WAVE_MM)),
            "r_s1_vs_laser": r, "p_s1_vs_laser": p, "slope_ols": sl["slope_ols"], "slope_tls": sl["slope_tls"],
            "rms_ratio_s1_over_laser": float(np.sqrt(np.mean(obs ** 2)) / np.sqrt(np.mean(exp ** 2))),
            "n_above_floor": int(big.sum()), "n_same_sign": agree,
            "share_same_sign": agree / big.sum() if big.sum() else np.nan,
            "p_sign_binomial": binomtest(agree, int(big.sum()), 0.5).pvalue if big.sum() else np.nan,
            "n_ratio": int(mid.sum()),
            "ratio_median": float(np.median(obs[mid] / exp[mid])) if mid.sum() else np.nan,
            "ratio_q25": float(np.percentile(obs[mid] / exp[mid], 25)) if mid.sum() else np.nan,
            "ratio_q75": float(np.percentile(obs[mid] / exp[mid], 75)) if mid.sum() else np.nan,
        })
        Lw = L[L.dwtd_cm.notna() & ~L.wtd_censored]
        rs, ps = fl.circular_shift_p(Lw.dwtd_cm.to_numpy(), Lw.dh_mm.to_numpy()) if len(Lw) >= 8 else (np.nan, np.nan)
        out.update({"r_dh_vs_dwtd": rs, "p_dh_vs_dwtd": ps})
        cr = L.coh.notna()
        if cr.sum() >= 8:
            Lc = L[cr]
            res = coh_residual(Lc)
            rch, pch = fl.circular_shift_p(np.abs(Lc.dh_mm.to_numpy()), res)
            ok = (Lc.dwtd_cm.notna() & ~Lc.wtd_censored).to_numpy()
            rcw, pcw = fl.circular_shift_p(np.abs(Lc.dwtd_cm.to_numpy()[ok]), res[ok])
            out.update({"r_coh_vs_abs_dh": rch, "p_coh_vs_abs_dh": pch, "r_coh_vs_abs_dwtd_laser_pairs": rcw,
                        "p_coh_vs_abs_dwtd_laser_pairs": pcw,
                        "partial_r_coh_dwtd_given_dh": fl.partial_r(res[ok], np.abs(Lc.dwtd_cm.to_numpy()[ok]),
                                                                    np.abs(Lc.dh_mm.to_numpy()[ok])),
                        "partial_r_coh_dh_given_dwtd": fl.partial_r(res[ok], np.abs(Lc.dh_mm.to_numpy()[ok]),
                                                                    np.abs(Lc.dwtd_cm.to_numpy()[ok]))})
    if len(W) >= 8:
        rw, pw = fl.circular_shift_p(W.dwtd_cm.to_numpy(), W.s1_dlos_mm.to_numpy())
        res = coh_residual(W[W.coh.notna()])
        rc, pc = fl.circular_shift_p(np.abs(W[W.coh.notna()].dwtd_cm.to_numpy()), res)
        out.update({"r_s1_vs_dwtd": rw, "p_s1_vs_dwtd": pw,
                    "slope_s1_mm_per_cm_wtd": fl.slopes(W.dwtd_cm, W.s1_dlos_mm)["slope_ols"],
                    "r_coh_vs_abs_dwtd": rc, "p_coh_vs_abs_dwtd": pc})
    return out


def seasonal(df: pd.DataFrame, acq, ov, runs) -> pd.DataFrame:
    """Annual semi-amplitude of the Sentinel-1 apparent displacement (consecutive pairs chained,
    a new segment at every missing link) and of the laser in LOS (at the usable overpasses, a new
    segment after the re-referencing), each with a linear trend and one intercept per segment."""
    rows = []
    for (track, w), d in df.groupby(["track", "window"]):
        d = d.sort_values("t1")
        dates = [d.t1.iloc[0]] + list(d.t2)
        s1, seg = fl.chain(d.s1_dlos_mm.to_numpy())
        a_s1 = fl.harmonic_amplitude(dates, s1, seg)
        o = ov[(ov.track == track) & ov.usable].sort_values("date")
        lz = o.laser_cm.to_numpy() * 10 * np.cos(np.radians(INCIDENCE[track]))
        lseg = laser_segment(o.time_utc, runs)
        a_l = fl.harmonic_amplitude(o.date, lz, lseg)
        rows.append({"track": track, "window": w, "s1_amplitude_mm": a_s1["amplitude"], "s1_n_dates": a_s1["n"],
                     "s1_segments": int(seg.max() + 1), "laser_amplitude_los_mm": a_l["amplitude"],
                     "laser_n_dates": a_l["n"], "ratio_s1_over_laser": a_s1["amplitude"] / a_l["amplitude"],
                     "s1_doy_max": a_s1["doy_max"], "laser_doy_max": a_l["doy_max"]})
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------------ controls

def controls(df: pd.DataFrame, ctx, p6_rc, n_sample: int = 300) -> tuple[pd.DataFrame, pd.DataFrame]:
    """The same per-pair statistics for single pixels elsewhere, each minus the grassland median
    of the same interferogram — the lake (B, all pixels), the grassland (C) and the ground outside
    (D, samples) and the rest of the mat (A) — and P6's 1×1 value among them. If stable ground
    follows the water table as P6 does, the link is regional or in the reference, not the mat's."""
    zones = {k: ctx.zones[k].values.astype(bool) for k in "ABCD"}
    r0, c0 = p6_rc
    near = np.zeros_like(zones["A"])
    near[r0 - 1:r0 + 2, c0 - 1:c0 + 2] = True
    picks = {}
    for z, m in zones.items():
        idx = np.argwhere(m & ~near)
        if len(idx) > n_sample:
            idx = idx[RNG.choice(len(idx), n_sample, replace=False)]
        picks[z] = idx
    crop = {"ascending": ctx.paths.cropped,
            "descending": ctx.paths.for_phase("field_p6_validation", track="descending").cropped}
    per_px, check = [], []
    for track, d in df[df.window == PRIMARY].groupby("track"):
        d = d[d.in_network].sort_values("mid").reset_index(drop=True)
        unw = load_layer(crop[track], "unw_phase", list(d.pair)).values
        corr = load_layer(crop[track], "corr", list(d.pair)).values
        refC = np.array([np.nanmedian(u[zones["C"] & np.isfinite(u)]) for u in unw])
        p6_dlos = (unw[:, r0, c0] - refC) * PHASE_TO_MM
        check.append({"track": track, "max_abs_diff_mm": float(np.nanmax(np.abs(p6_dlos - d.s1_dlos_mm)))})
        exp = d.laser_dlos_mm.to_numpy()
        dwtd = np.where(d.wtd_censored, np.nan, d.dwtd_cm.to_numpy())
        dh = d.dh_mm.to_numpy()
        for z, idx in picks.items():
            for r, c in idx:
                s1 = (unw[:, r, c] - refC) * PHASE_TO_MM
                coh = corr[:, r, c]
                per_px.append({"track": track, "zone": z, "row": r, "col": c,
                               **pixel_stats(s1, coh, exp, dwtd, dh, d)})
        per_px.append({"track": track, "zone": "P6", "row": r0, "col": c0,
                       **pixel_stats(p6_dlos, corr[:, r0, c0], exp, dwtd, dh, d)})
    return pd.DataFrame(per_px), pd.DataFrame(check)


STATS = {"r_s1_vs_laser": "S1 ΔLOS vs laser ΔLOS (r)", "slope_s1_vs_laser": "S1 ΔLOS vs laser ΔLOS (OLS slope)",
         "r_s1_vs_dwtd": "S1 ΔLOS vs ΔWTD (r)", "r_coh_vs_abs_dwtd": "coherence vs |ΔWTD| (r, season removed)",
         "r_coh_vs_abs_dh": "coherence vs |Δsurface| (r, season removed)"}


def pixel_stats(s1, coh, exp, dwtd, dh, d) -> dict:
    out = {}
    ok = np.isfinite(s1) & np.isfinite(exp)
    out["r_s1_vs_laser"] = float(np.corrcoef(s1[ok], exp[ok])[0, 1]) if ok.sum() >= 8 else np.nan
    out["slope_s1_vs_laser"] = fl.slopes(exp[ok], s1[ok])["slope_ols"] if ok.sum() >= 8 else np.nan
    ok = np.isfinite(s1) & np.isfinite(dwtd)
    out["r_s1_vs_dwtd"] = float(np.corrcoef(s1[ok], dwtd[ok])[0, 1]) if ok.sum() >= 8 else np.nan
    okc = np.isfinite(coh)
    if okc.sum() >= 8:
        res = np.full(len(coh), np.nan)
        res[okc] = fl.season_residual(coh[okc], d.mid[okc], d.dt_days.to_numpy()[okc])
        ok = okc & np.isfinite(dwtd)
        out["r_coh_vs_abs_dwtd"] = float(np.corrcoef(res[ok], np.abs(dwtd[ok]))[0, 1]) if ok.sum() >= 8 else np.nan
        ok = okc & np.isfinite(dh)
        out["r_coh_vs_abs_dh"] = float(np.corrcoef(res[ok], np.abs(dh[ok]))[0, 1]) if ok.sum() >= 8 else np.nan
    return out


def control_summary(px: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for track, g in px.groupby("track"):
        p6 = g[g.zone == "P6"].iloc[0]
        for stat in STATS:
            for z in "ABCD":
                v = g[(g.zone == z)][stat].dropna()
                rows.append({"track": track, "statistic": stat, "zone": z, "n_pixels": len(v),
                             "zone_median": float(v.median()) if len(v) else np.nan,
                             "zone_q05": float(v.quantile(0.05)) if len(v) else np.nan,
                             "zone_q95": float(v.quantile(0.95)) if len(v) else np.nan,
                             "p6": float(p6[stat]),
                             "p6_percentile": float((v < p6[stat]).mean() * 100) if len(v) else np.nan})
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------------ pair length and closure

def all_pairs(p6, ov, wtd_s1, runs, max_dt: int = 48) -> pd.DataFrame:
    """Every pair ≤ ``max_dt`` days (not only consecutive ones) with the X-058 laser rules — to see
    at which pair length the phase stops following the surface."""
    ovi = ov.set_index(["date", "track"])
    rows = []
    for r in p6[(p6.dt_days <= max_dt)].itertuples(index=False):
        a, b = ovi.loc[(r.ref_date, r.track)], ovi.loc[(r.sec_date, r.track)]
        if not (a.usable and b.usable) or not np.isfinite(r.dlos_vs_C_mm) \
                or ((runs.end > pd.Timestamp(r.ref_date)) & (runs.start < pd.Timestamp(r.sec_date))).any():
            continue
        wa, wb = wtd_s1.loc[(r.ref_date, r.track)], wtd_s1.loc[(r.sec_date, r.track)]
        rows.append({"track": r.track, "window": r.window, "pair": r.pair, "dt_days": r.dt_days,
                     "s1_dlos_mm": r.dlos_vs_C_mm, "coh": r.coh_median,
                     "laser_dlos_mm": (b.laser_cm - a.laser_cm) * 10 * np.cos(np.radians(INCIDENCE[r.track])),
                     "dwtd_cm": np.nan if (wa.wtd_censored or wb.wtd_censored) else wb.wtd_at - wa.wtd_at})
    return pd.DataFrame(rows)


def by_pair_length(ap: pd.DataFrame) -> pd.DataFrame:
    rows = []
    bins = {"12 d": (12, 12), "24 d": (24, 24), "24 d, laser < λ/4": (24, 24), "36–48 d": (36, 48)}
    for (t, w), g in ap.groupby(["track", "window"]):
        for lab, (lo, hi) in bins.items():
            d = g[(g.dt_days >= lo) & (g.dt_days <= hi)]
            if "λ/4" in lab:
                d = d[d.laser_dlos_mm.abs() < QUARTER_WAVE_MM]
            if len(d) < 8:
                continue
            sl = fl.slopes(d.laser_dlos_mm, d.s1_dlos_mm)
            dw = d.dropna(subset=["dwtd_cm"])
            rows.append({"track": t, "window": w, "pairs": lab, "n": len(d), "median_coh": float(d.coh.median()),
                         "r": sl["r"], "slope_ols": sl["slope_ols"], "slope_tls": sl["slope_tls"],
                         "share_same_sign": float(np.mean(np.sign(d.laser_dlos_mm) == np.sign(d.s1_dlos_mm))),
                         "partial_r_s1_laser_given_dwtd": fl.partial_r(dw.s1_dlos_mm, dw.laser_dlos_mm, dw.dwtd_cm)
                         if len(dw) >= 8 else np.nan,
                         "partial_r_s1_dwtd_given_laser": fl.partial_r(dw.s1_dlos_mm, dw.dwtd_cm, dw.laser_dlos_mm)
                         if len(dw) >= 8 else np.nan, "n_partial": len(dw)})
    return pd.DataFrame(rows)


def closure(p6, ov, runs) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Triplets t1, t1+12 d, t1+24 d: does the 24-day interferogram equal the sum of its two
    12-day ones? Motion closes; a phase bias that depends on the pair (fading) or a 24-day phase
    lost to decorrelation does not. Also: which of the two follows the laser's 24-day change."""
    ovi = ov.set_index(["date", "track"])
    rows = []
    for (t, w), g in p6.groupby(["track", "window"]):
        d = g.set_index(["ref_date", "sec_date"])
        for (a, c), r in d[d.dt_days == 24].iterrows():
            b = str((pd.Timestamp(a) + pd.Timedelta(days=12)).date())
            if (a, b) not in d.index or (b, c) not in d.index:
                continue
            s12 = d.loc[(a, b)].dlos_vs_C_mm + d.loc[(b, c)].dlos_vs_C_mm
            la, lc = ovi.loc[(a, t)], ovi.loc[(c, t)]
            across = ((runs.end > pd.Timestamp(a)) & (runs.start < pd.Timestamp(c))).any()
            lz = (lc.laser_cm - la.laser_cm) * 10 * np.cos(np.radians(INCIDENCE[t])) \
                if la.usable and lc.usable and not across else np.nan
            rows.append({"track": t, "window": w, "t1": a, "t3": c, "sum_12d_mm": s12, "s1_24d_mm": r.dlos_vs_C_mm,
                         "closure_mm": r.dlos_vs_C_mm - s12, "laser_24d_los_mm": lz,
                         "coh_24d": r.coh_median, "coh_12d_mean": (d.loc[(a, b)].coh_median + d.loc[(b, c)].coh_median) / 2})
    tri = pd.DataFrame(rows)
    summ = []
    for (t, w), g in tri.groupby(["track", "window"]):
        g = g.dropna(subset=["sum_12d_mm", "s1_24d_mm"])
        lz = g.dropna(subset=["laser_24d_los_mm"])
        r12, p12 = fl.circular_shift_p(lz.laser_24d_los_mm.to_numpy(), lz.sum_12d_mm.to_numpy())
        r24, p24 = fl.circular_shift_p(lz.laser_24d_los_mm.to_numpy(), lz.s1_24d_mm.to_numpy())
        summ.append({"track": t, "window": w, "n_triplets": len(g), "closure_median_abs_mm": float(g.closure_mm.abs().median()),
                     "closure_sd_mm": float(g.closure_mm.std()), "sd_sum_12d_mm": float(g.sum_12d_mm.std()),
                     "r_24d_vs_sum_12d": float(g.s1_24d_mm.corr(g.sum_12d_mm)),
                     "n_with_laser": len(lz), "r_sum12_vs_laser": r12, "p_sum12_vs_laser": p12,
                     "r_24d_vs_laser": r24, "p_24d_vs_laser": p24,
                     "r_closure_vs_laser": float(lz.closure_mm.corr(lz.laser_24d_los_mm)),
                     "median_coh_12d": float(g.coh_12d_mean.median()), "median_coh_24d": float(g.coh_24d.median())})
    return tri, pd.DataFrame(summ)


# ------------------------------------------------------------------------------ X-048 check

def x048_reproduction(acq, p6, ov, wtd_s1) -> pd.DataFrame:
    """X-048's slope with X-048's settings (all pairs ≤ 24 d, 3×3, laser including filled values,
    snow rule only) through this script's arithmetic — to show the difference in the new numbers
    comes from the new choices, not from the code."""
    ovi = ov.set_index(["date", "track"])
    old = pd.read_csv(DLV / "field_x047_x048" / "p6_laser_vs_radar.csv").set_index("track")
    rows = []
    for track in ("ascending", "descending"):
        d = p6[(p6.track == track) & (p6.window == "3x3") & (p6.dt_days <= 24)]
        e, o = [], []
        for r in d.itertuples(index=False):
            a, b = ovi.loc[(r.ref_date, track)], ovi.loc[(r.sec_date, track)]
            wa, wb = wtd_s1.loc[(r.ref_date, track)], wtd_s1.loc[(r.sec_date, track)]
            if a.snow_72h or b.snow_72h or not np.isfinite(a.laser_cm) or not np.isfinite(b.laser_cm) \
                    or wa.wtd_censored or wb.wtd_censored or not np.isfinite(r.dlos_vs_C_mm) or not np.isfinite(r.coh_median):
                continue
            e.append((b.laser_cm - a.laser_cm) * 10 * np.cos(np.radians(INCIDENCE[track])))
            o.append(r.dlos_vs_C_mm)
        s = fl.slopes(e, o)
        rows.append({"track": track, "n_pairs_here": s["n"], "n_pairs_x048": int(old.loc[track, "n_pairs"]),
                     "slope_here": s["slope_ols"], "slope_x048": float(old.loc[track, "slope_los_per_expected_los"])})
    return pd.DataFrame(rows)


def web_series(df: pd.DataFrame, ov: pd.DataFrame, runs: pd.DataFrame) -> pd.DataFrame:
    """The Sentinel-1 apparent displacement chained from consecutive pairs (mean 0) and the laser
    in LOS at the usable overpasses, its level matched to the chain within each laser segment
    (one offset per segment: the laser's level across a filled stretch is unknown). The same
    numbers as panel 3 of the time-series figure, for the web atlas."""
    rows = []
    for (t, w), d in df.groupby(["track", "window"]):
        d = d.sort_values("t1")
        dates = pd.to_datetime([d.t1.iloc[0]] + list(d.t2))
        s1, seg = fl.chain(d.s1_dlos_mm.to_numpy())
        s1 = pd.Series(s1, index=dates) - np.nanmean(s1)
        o = ov[(ov.track == t) & ov.usable].copy()
        o["day"] = pd.to_datetime(o.date)
        o["los"] = o.laser_cm * 10 * np.cos(np.radians(INCIDENCE[t]))
        o["seg"] = laser_segment(o.time_utc, runs)
        o["laser_matched"] = np.nan
        for k, g in o.groupby("seg"):
            s1g = s1.reindex(g.day.to_numpy())
            ok = s1g.notna().to_numpy()
            if ok.sum() >= 2:
                o.loc[g.index, "laser_matched"] = g.los - float(np.mean(g.los.to_numpy()[ok] - s1g.to_numpy()[ok]))
        lm = o.set_index("day")
        for dt, v in s1.items():
            rows.append({"date": f"{dt:%Y-%m-%d}", "track": t, "window": w, "s1_chain_mm": v,
                         "laser_los_matched_mm": lm.laser_matched.get(dt, np.nan),
                         "laser_segment": int(lm.seg.get(dt, -1)) if dt in lm.index else None})
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------------ figures

def fig_timeseries(out, df, ov, laser_flags, wtd, zc, runs):
    x0, x1 = pd.Timestamp(PERIOD[0]), pd.Timestamp(PERIOD[1] + " 23:59")
    fig, ax = plt.subplots(4, 1, figsize=(12, 13), sharex=True)
    laser_flags = laser_flags.assign(time_utc=laser_flags.time_utc.dt.tz_convert(None))
    ov = ov.assign(time_utc=pd.to_datetime(ov.time_utc).dt.tz_convert(None))
    wtd = wtd.set_axis(wtd.index.tz_convert(None))
    lf = laser_flags[(laser_flags.time_utc >= x0) & (laser_flags.time_utc <= x1)]
    meas = lf.surface_cm.where(~lf.snow_72h & ~lf.filled & ~lf.outlier)
    ax[0].plot(lf.time_utc, lf.surface_cm_raw.where(lf.snow_72h & ~lf.filled), lw=0.5, color="0.75", label="snow period")
    ax[0].plot(lf.time_utc, lf.surface_cm_raw.where(lf.filled), lw=1.2, ls="--", color="tab:purple", label="filled in the delivery (not data)")
    ax[0].plot(lf.time_utc, meas, lw=0.7, color="k", label="laser surface, measured, snow-free")
    for t, g in ov[ov.usable].groupby("track"):
        ax[0].scatter(pd.to_datetime(g.time_utc), g.laser_cm, s=14, color=COLOR[t], zorder=3, label=f"at {t} overpasses (usable)")
    for r in runs.itertuples(index=False):
        if x0 <= r.end <= x1:
            ax[0].axvline(r.end, color="tab:purple", lw=0.8, ls=":")
    ax[0].set_ylabel("laser surface (cm, + up)")
    ax[0].legend(fontsize=7, ncol=2, loc="lower left")
    ax[0].set_title("P6/CR — 1 laser surface · 2 water table · 3 Sentinel-1 apparent displacement · 4 coherence (2022–2024)", fontsize=10)
    w = wtd.loc[x0:x1]
    ax[1].plot(w.index, w["P6"], lw=0.7, color="tab:cyan", label="WTD P6 (cm, − below surface)")
    a1 = ax[1].twinx()
    a1.plot(w.index, w["CR_raw"], lw=0.6, color="tab:blue", alpha=0.6, label="raw CR water level (cm)")
    a1.set_ylabel("raw CR level (cm)", color="tab:blue")
    ax[1].set_ylabel("WTD P6 (cm)")
    ax[1].legend(fontsize=7, loc="lower left")
    for t in ("ascending", "descending"):
        d = df[(df.track == t) & (df.window == PRIMARY)].sort_values("t1")
        dates = pd.to_datetime([d.t1.iloc[0]] + list(d.t2))
        s1, seg = fl.chain(d.s1_dlos_mm.to_numpy())
        s1 = pd.Series(s1, index=dates)
        ax[2].plot(dates, s1 - s1.mean(), lw=1.0, color=COLOR[t], label=f"Sentinel-1 {t}, chained consecutive pairs (mean 0)")
        o = ov[(ov.track == t) & ov.usable].copy()
        o["los"] = o.laser_cm * 10 * np.cos(np.radians(INCIDENCE[t]))
        o["day"] = o.time_utc.dt.normalize()
        o["seg"] = laser_segment(o.time_utc, runs)
        first = True
        for k, g in o.groupby("seg"):
            s1g = s1.reindex(g.day.to_numpy())
            ok = s1g.notna().to_numpy()
            if ok.sum() < 2:
                continue
            shift = float(np.mean(g.los.to_numpy()[ok] - s1g.to_numpy()[ok])) + s1.mean()
            ax[2].plot(g.time_utc, g.los - shift, "o", ms=3, mfc="none", color=COLOR[t],
                       label=f"laser in {t} LOS, level matched per laser segment" if first else None)
            first = False
    ax[2].axhspan(-QUARTER_WAVE_MM, QUARTER_WAVE_MM, color="tab:cyan", alpha=0.08)
    ax[2].set_ylabel("LOS (mm, + toward satellite); S1 mean 0, laser level-matched")
    ax[2].legend(fontsize=7, ncol=2, loc="lower left")
    for t in ("ascending", "descending"):
        d = df[(df.track == t) & (df.window == PRIMARY)]
        ax[3].plot(pd.to_datetime(d["mid"]), d.coh, "o-", ms=3, lw=0.6, color=COLOR[t], label=f"P6 {t} (1×1)")
        c = zc[(zc.track == t) & (zc.zone == "C") & zc.pair.isin(d.pair)].merge(d[["pair", "mid"]], on="pair")
        ax[3].plot(pd.to_datetime(c["mid"]), c.mean_coh, lw=0.8, ls="--", color=COLOR[t], alpha=0.6, label=f"grassland C {t}")
    ax[3].set_ylabel("coherence, consecutive pair")
    ax[3].set_ylim(0, 1)
    ax[3].legend(fontsize=7, ncol=4, loc="lower left")
    ax[3].set_xlim(x0, x1)
    fig.tight_layout()
    fig.savefig(out / "fig_p6_timeseries.png", dpi=160)
    plt.close(fig)


def fig_pairs(out, df, floor_mm):
    d = df[(df.window == PRIMARY) & df.laser_usable & df.s1_dlos_mm.notna()]
    fig, ax = plt.subplots(2, 2, figsize=(11, 10))
    lim = 1.1 * np.nanmax(np.abs(np.r_[d.laser_dlos_mm, d.s1_dlos_mm, QUARTER_WAVE_MM]))
    for t, g in d.groupby("track"):
        ax[0, 0].scatter(g.laser_dlos_mm, g.s1_dlos_mm, s=16, color=COLOR[t], label=t)
        ax[0, 1].scatter(g.dwtd_cm, g.dh_mm, s=16, color=COLOR[t], label=t)
        ax[1, 0].scatter(g.dwtd_cm.where(~g.wtd_censored), g.s1_dlos_mm, s=16, color=COLOR[t], label=t)
    ax[0, 0].plot([-lim, lim], [-lim, lim], "k--", lw=0.8, label="1 : 1 (the phase follows the surface)")
    ax[0, 0].axvspan(-floor_mm, floor_mm, color="0.85", label="laser noise floor")
    for s in (-1, 1):
        ax[0, 0].axvline(s * QUARTER_WAVE_MM, color="tab:cyan", lw=0.8)
    ax[0, 0].set(xlim=(-lim, lim), ylim=(-lim, lim), xlabel="laser surface change in LOS (mm)",
                 ylabel="Sentinel-1 apparent displacement, P6 − grassland (mm)")
    ax[0, 0].set_title("Between consecutive acquisitions (cyan: ±λ/4)", fontsize=9)
    ax[0, 0].legend(fontsize=7)
    ax[0, 1].set(xlabel="ΔWTD at P6 (cm)", ylabel="laser surface change (mm, vertical)")
    ax[0, 1].set_title("The surface follows the water", fontsize=9)
    ax[1, 0].set(xlabel="ΔWTD at P6 (cm)", ylabel="Sentinel-1 apparent displacement (mm)")
    ax[1, 0].set_title("The phase against the water", fontsize=9)
    for t, g in df[(df.window == PRIMARY) & df.coh.notna()].groupby("track"):
        g = g.sort_values("mid")
        res = coh_residual(g)
        ax[1, 1].scatter(np.abs(g.dwtd_cm.where(~g.wtd_censored)), res, s=14, color=COLOR[t], label=f"|ΔWTD| (cm), {t}")
        ax[1, 1].scatter(np.abs(g.dh_mm) / 10, res, s=14, marker="x", color=COLOR[t], label=f"|Δsurface| (cm), {t}")
    ax[1, 1].set(xlabel="absolute change (cm)", ylabel="coherence, season removed")
    ax[1, 1].set_title("Coherence against water and surface change", fontsize=9)
    ax[1, 1].legend(fontsize=7)
    for a in ax.flat:
        a.axhline(0, color="0.6", lw=0.5)
        a.axvline(0, color="0.6", lw=0.5)
    fig.tight_layout()
    fig.savefig(out / "fig_p6_pairs.png", dpi=160)
    plt.close(fig)


def fig_controls(out, px):
    fig, ax = plt.subplots(2, len(STATS), figsize=(18, 7))
    zc = {"A": "tab:green", "B": "tab:blue", "C": "tab:olive", "D": "0.5"}
    zl = {"A": "mat (A)", "B": "lake (B)", "C": "grassland (C)", "D": "outside (D)"}
    for i, t in enumerate(("ascending", "descending")):
        g = px[px.track == t]
        for j, (stat, label) in enumerate(STATS.items()):
            a = ax[i, j]
            for z in "ABCD":
                v = g[g.zone == z][stat].dropna()
                if len(v):
                    a.hist(v, bins=25, histtype="step", color=zc[z], label=zl[z], density=True)
            a.axvline(g[g.zone == "P6"][stat].iloc[0], color="k", lw=2, label="P6")
            a.set_title(f"{t}: {label}", fontsize=8)
            if i == 0 and j == 0:
                a.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(out / "fig_p6_controls.png", dpi=150)
    plt.close(fig)


# ------------------------------------------------------------------------------ README

def f2(x, n=2):
    return "—" if x is None or not np.isfinite(x) else f"{x:.{n}f}"


def fp(p):
    return "—" if p is None or not np.isfinite(p) else ("≤ 0.001" if p <= 0.001 else f"{p:.3f}")


def readme(out, df, st, sea, ctl, check, x048, floor_mm, runs, bl, cs):
    s = st.set_index(["track", "window", "subset"])
    L = ["# X-059 — P6/CR, the primary validation site (exploratory)\n",
         "Generated by `05_code/SAr-adapted/scripts/field_p6_validation.py` (outputs here only; field data). The "
         "supervisor's request of 2026-09-28 (`01_admin/supervisor/2026-09-28_email_p6_priority.md`), design in "
         "`PLAN.md`, laser verified first in `laser_qc/` (X-058). Every number is computed by the script.\n",
         "## Settings\n",
         "- **Observations**: consecutive Sentinel-1 acquisitions of each track, 2022–2024 (pairs t_i → t_i+1; they do "
         "not overlap). Earlier work (X-048) used every pair ≤ 24 d.",
         "- **Extraction**: the P6 pixel (1×1, 40 m) as primary — the smallest area the grid allows; 3×3 (120 m, also "
         "covering P4, P5, P7–P9) for comparison. Finer needs reprocessing (D-021); boardwalk flags pending (X-060).",
         "- **Laser**: X-058 conventions — cm, increase = up, as delivered (no cos 10°), nearest measured hourly value "
         "within ±1 h of the overpass, snow-free (72 h), no outlier, no filled value; no pair across a filled stretch "
         f"(a level may have been re-referenced inside one: {', '.join(f'{r.start:%Y-%m-%d}–{r.end:%Y-%m-%d}' for r in runs.itertuples())}).",
         "- **Units made comparable**: laser Δh (mm, vertical) → LOS with cos(incidence) — what a purely vertical motion "
         "gives; Sentinel-1 unwrapped phase, P6 − grassland (C) median, × −λ/4π → mm LOS; both positive toward the "
         "satellite (up).",
         f"- **Ratio** only where the laser change is above its noise floor ({f2(floor_mm,1)} mm, X-058) and below λ/4 "
         f"({f2(QUARTER_WAVE_MM,1)} mm); beyond λ/4 the phase is ambiguous and a ratio means nothing.",
         "- p-values from circular shifts in time (floor 1/2001); sign agreement: binomial against 0.5. Exploratory: "
         "one site, tens of pairs, many statistics.\n",
         "## Counts\n",
         "| track | consecutive pairs | in the HyP3 network | with a usable laser change | with WTD (not censored) |",
         "|---|---|---|---|---|"]
    for t in ("ascending", "descending"):
        r = s.loc[(t, PRIMARY, "all")]
        L.append(f"| {t} | {int(r.n_pairs)} | {int(r.n_in_network)} | {int(r.n_laser)} | {int(r.n_wtd)} |")
    L += ["\n## 1. How much the surface moves between consecutive acquisitions\n",
          "| track | median \\|Δh\\| (mm) | max \\|Δh\\| (mm) | laser LOS change beyond λ/4 | surface vs ΔWTD r (p) |",
          "|---|---|---|---|---|"]
    for t in ("ascending", "descending"):
        r = s.loc[(t, PRIMARY, "all")]
        L.append(f"| {t} | {f2(r.median_abs_dh_mm,1)} | {f2(r.max_abs_dh_mm,1)} | {r.share_laser_beyond_quarter_wave:.0%} | "
                 f"{f2(r.r_dh_vs_dwtd)} ({fp(r.p_dh_vs_dwtd)}) |")
    L += ["\n## 2. Does the Sentinel-1 apparent displacement follow it? (same pairs, same units)\n",
          "| track | window | subset | n | median \\|laser LOS\\| (mm) | median \\|S1\\| (mm) | r (p) | slope OLS | slope TLS | RMS ratio S1/laser | same sign | per-pair ratio S1/laser, median [IQR] (n) |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for (t, w, sub), r in s.iterrows():
        if not np.isfinite(r.get("slope_ols", np.nan)):
            continue
        L.append(f"| {t} | {w} | {sub} | {int(r.n_laser)} | {f2(r.median_abs_laser_dlos_mm,1)} | {f2(r.median_abs_s1_dlos_mm,1)} | "
                 f"{f2(r.r_s1_vs_laser)} ({fp(r.p_s1_vs_laser)}) | {f2(r.slope_ols)} | {f2(r.slope_tls)} | "
                 f"{f2(r.rms_ratio_s1_over_laser)} | {int(r.n_same_sign)}/{int(r.n_above_floor)} (p {fp(r.p_sign_binomial)}) | "
                 f"{f2(r.ratio_median)} [{f2(r.ratio_q25)}, {f2(r.ratio_q75)}] ({int(r.n_ratio)}) |")
    L += ["\n`dry` = neither date wet (RH ≥ 95 % or rain in the previous 3 h) nor frozen (X-050 flags). The total-least-"
          "squares slope treats both variables as noisy; OLS is attenuated by noise in the laser change. "
          "A phase that recorded the motion would give slopes and ratios near 1 and agree in sign.\n",
          "**Seasonal amplitude** (annual harmonic + trend + one level per segment; Sentinel-1 chained from "
          "consecutive pairs, which accumulates noise — descriptive):\n",
          "| track | window | S1 amplitude (mm LOS) | laser amplitude (mm LOS) | ratio | S1 chain segments |",
          "|---|---|---|---|---|---|"]
    for r in sea.itertuples(index=False):
        L.append(f"| {r.track} | {r.window} | {f2(r.s1_amplitude_mm,1)} | {f2(r.laser_amplitude_los_mm,1)} | "
                 f"{f2(r.ratio_s1_over_laser)} | {r.s1_segments} |")
    L += ["\n## 3. Water table, surface and coherence on the same pairs\n",
          "| track | S1 vs ΔWTD r (p) | slope (mm per cm) | coherence vs \\|ΔWTD\\| r (p) | on laser pairs: coh vs \\|ΔWTD\\| r (p) | coh vs \\|Δsurface\\| r (p) | partial: WTD given surface | partial: surface given WTD |",
          "|---|---|---|---|---|---|---|---|"]
    for t in ("ascending", "descending"):
        r = s.loc[(t, PRIMARY, "all")]
        L.append(f"| {t} | {f2(r.r_s1_vs_dwtd)} ({fp(r.p_s1_vs_dwtd)}) | {f2(r.slope_s1_mm_per_cm_wtd)} | "
                 f"{f2(r.r_coh_vs_abs_dwtd)} ({fp(r.p_coh_vs_abs_dwtd)}) | {f2(r.r_coh_vs_abs_dwtd_laser_pairs)} "
                 f"({fp(r.p_coh_vs_abs_dwtd_laser_pairs)}) | {f2(r.r_coh_vs_abs_dh)} ({fp(r.p_coh_vs_abs_dh)}) | "
                 f"{f2(r.partial_r_coh_dwtd_given_dh)} | {f2(r.partial_r_coh_dh_given_dwtd)} |")
    L += ["\nCoherence with the annual cycle and the pair length regressed out.\n",
          "## 4. Controls: is any of this specific to P6?\n",
          "The same statistics for single pixels elsewhere (each minus the same grassland median): the lake (B, all "
          "pixels), the grassland (C), the ground outside the site (D) and the rest of the mat (A), ≤ 300 sampled per "
          "zone, P6's 3×3 excluded. P6's value, the zone median, and P6's percentile within the zone. A statistic "
          "that stable ground or the lake reproduces is not evidence about the mat.\n",
          "| track | statistic | P6 | mat A median (P6 pct) | lake B median (P6 pct) | grassland C median (P6 pct) | outside D median (P6 pct) |",
          "|---|---|---|---|---|---|---|"]
    for (t, stat), g in ctl.groupby(["track", "statistic"], sort=False):
        z = g.set_index("zone")
        L.append(f"| {t} | {STATS[stat]} | {f2(z.p6.iloc[0])} | " + " | ".join(
            f"{f2(z.loc[k, 'zone_median'])} ({f2(z.loc[k, 'p6_percentile'], 0)})" for k in "ABCD") + " |")
    L += ["\n## 5. At which pair length does the phase stop following the surface?\n",
          "All pairs ≤ 48 d (not only consecutive), same laser rules. Partial correlations remove the water-table "
          "change: if the phase followed the surface only because both follow the water, the first would vanish.\n",
          "| track | window | pairs | n | median coherence | r | slope OLS | slope TLS | same sign | partial r (S1, laser \\| ΔWTD) | partial r (S1, ΔWTD \\| laser) |",
          "|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in bl.itertuples(index=False):
        L.append(f"| {r.track} | {r.window} | {r.pairs} | {r.n} | {f2(r.median_coh)} | {f2(r.r)} | {f2(r.slope_ols)} | "
                 f"{f2(r.slope_tls)} | {r.share_same_sign:.0%} | {f2(r.partial_r_s1_laser_given_dwtd)} | "
                 f"{f2(r.partial_r_s1_dwtd_given_laser)} |")
    L += ["\n## 6. Closure: is the 24-day interferogram the sum of its two 12-day ones?\n",
          "Surface motion closes (φ₂₄ = φ₁₂ₐ + φ₁₂ᵦ). A pair-dependent bias (fading) or a 24-day phase lost to "
          "decorrelation does not. The laser tells which side is right.\n",
          "| track | window | triplets | closure median \\|·\\| (mm) | closure SD (mm) | SD of 12+12 (mm) | r(24 d, 12+12) | with laser | r(12+12, laser) (p) | r(24 d, laser) (p) | r(closure, laser) | median coherence 12 d / 24 d |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in cs.itertuples(index=False):
        L.append(f"| {r.track} | {r.window} | {r.n_triplets} | {f2(r.closure_median_abs_mm,1)} | {f2(r.closure_sd_mm,1)} | "
                 f"{f2(r.sd_sum_12d_mm,1)} | {f2(r.r_24d_vs_sum_12d)} | {r.n_with_laser} | {f2(r.r_sum12_vs_laser)} "
                 f"({fp(r.p_sum12_vs_laser)}) | {f2(r.r_24d_vs_laser)} ({fp(r.p_24d_vs_laser)}) | {f2(r.r_closure_vs_laser)} | "
                 f"{f2(r.median_coh_12d)} / {f2(r.median_coh_24d)} |")
    L += ["\n## Checks\n",
          f"- P6 1×1 phase read here from the interferograms equals the first-deliverable extraction: largest "
          f"difference {', '.join(f'{r.track} {r.max_abs_diff_mm:.1e} mm' for r in check.itertuples())}.",
          "- X-048's slope recomputed with X-048's settings (all pairs ≤ 24 d, 3×3, filled laser values kept): "
          + "; ".join(f"{r.track} {f2(r.slope_here)} here vs {f2(r.slope_x048)} in X-048 (n {r.n_pairs_here} vs "
                      f"{r.n_pairs_x048})" for r in x048.itertuples()) + ".\n",
          "## How to read this — and what it does not show\n",
          "- Section 2 answers *does the C-band phase follow the measured surface?* A slope and ratios far below 1, "
          "weak correlation and sign agreement near chance mean it does not follow the mechanical motion. That alone "
          "does **not** say what the remaining few millimetres are.",
          "- Section 3 and 4 are what the supervisor asks for before any dielectric reading: whether the phase and the "
          "coherence go with the water table, and whether the lake and stable ground do the same. If stable ground "
          "reproduces P6's link to ΔWTD, the link is regional (weather, the grassland reference), not a property of the "
          "mat.",
          "- Sections 5 and 6 test the two readings of a phase that follows the laser: real motion (closes over "
          "triplets; the 12+12 sum and the 24-day pair both follow the laser unless the 24-day pair has decorrelated) "
          "versus a moisture-driven short-pair bias (does not close; would also appear on stable ground, section 4, and "
          "would not survive removing the water-table change).",
          "- Not shown here: backscatter and UAV (next), P7–P9 (X-061), boardwalk flags (X-060), finer than 40 m (D-021).\n",
          "## Files\n",
          "| File | Content |", "|---|---|",
          "| `p6_pairs.csv` | every consecutive pair × track × window: laser, WTD, CR, Sentinel-1, coherence, flags |",
          "| `p6_stats.csv` | section 1–3 statistics per track × window × subset |",
          "| `p6_seasonal.csv` | seasonal amplitudes |",
          "| `p6_pair_length.csv`, `p6_closure_triplets.csv`, `p6_closure_summary.csv` | sections 5 and 6 |",
          "| `p6_controls_pixels.csv`, `p6_controls_summary.csv` | section 4 |",
          "| `p6_web_series.csv` | panel 3 of the time-series figure as numbers (for the web atlas) |",
          "| `p6_check_extraction.csv`, `p6_check_x048.csv` | the two checks |",
          "| `fig_p6_timeseries.png` | the four quantities on one time axis, 2022–2024 |",
          "| `fig_p6_pairs.png` | pair by pair: phase vs laser, surface vs water, phase vs water, coherence |",
          "| `fig_p6_controls.png` | P6 among the lake, grassland, outside and mat pixels |"]
    (out / "README.md").write_text("\n".join(L) + "\n")


# ------------------------------------------------------------------------------ main

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default=str(DLV / "field_p6"))
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    acq, p6, ov, qc, interp, wtd_s1, wet, zc = load_inputs()
    p6_all = p6
    runs = filled_runs(interp)
    floor_mm = float(qc["noise_floor_mm"])
    wtd = field.load_wtd_hourly()
    df = pair_table(acq, p6, ov, wtd_s1, wet, wtd, runs)
    cC = zc[zc.zone == "C"][["track", "pair", "mean_coh"]].rename(columns={"mean_coh": "coh_grassland_C"})
    df = df.merge(cC, on=["track", "pair"], how="left")
    df.to_csv(out / "p6_pairs.csv", index=False)
    web_series(df, ov, runs).round(4).to_csv(out / "p6_web_series.csv", index=False)

    rows = []
    for (t, w), d in df.groupby(["track", "window"]):
        rows.append({"track": t, "window": w, "subset": "all", **pair_stats(d, floor_mm)})
        rows.append({"track": t, "window": w, "subset": "dry", **pair_stats(d[~d.wet_either & ~d.frozen_either], floor_mm)})
    st = pd.DataFrame(rows)
    st.to_csv(out / "p6_stats.csv", index=False)
    sea = seasonal(df, acq, ov, runs)
    sea.to_csv(out / "p6_seasonal.csv", index=False)

    ctx = start("field_p6_validation", mount=False, git=False)
    pix = pd.read_csv(DLV / "field_first" / "plot_pixels.csv").set_index("plot")
    px, check = controls(df, ctx, (int(pix.loc["P6", "row"]), int(pix.loc["P6", "col"])))
    px.to_csv(out / "p6_controls_pixels.csv", index=False)
    ctl = control_summary(px)
    ctl.to_csv(out / "p6_controls_summary.csv", index=False)
    x048 = x048_reproduction(acq, p6, ov, wtd_s1)
    check.to_csv(out / "p6_check_extraction.csv", index=False)
    x048.to_csv(out / "p6_check_x048.csv", index=False)
    bl = by_pair_length(all_pairs(p6_all, ov, wtd_s1, runs))
    bl.to_csv(out / "p6_pair_length.csv", index=False)
    tri, cs = closure(p6_all, ov, runs)
    tri.to_csv(out / "p6_closure_triplets.csv", index=False)
    cs.to_csv(out / "p6_closure_summary.csv", index=False)

    lf = pd.read_csv(DLV / "field_p6" / "laser_qc" / "laser_flags.csv", parse_dates=["time_utc"])
    fig_timeseries(out, df, ov, lf, wtd, zc, runs)
    fig_pairs(out, df, floor_mm)
    fig_controls(out, px)
    readme(out, df, st, sea, ctl, check, x048, floor_mm, runs, bl, cs)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()

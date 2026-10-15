"""X-119 — the v4 architecture head-to-head at P6: plain chain, A (re-levelled by SPAMS), B (v3+), C (v3+ with the
SPAMS winter anchor), one scorer, the same epochs, the whole design span.

B is fusion v3's estimator (``fusion_v3_wp5_x078``: buoyancy prior + mean-reverting residual smoothed from the pairs,
noise by revisit class and state, detectability inflation) extended to 12-day pairs, with the prior's water table from
the ERA5 bucket in every year. C replaces B's prior step across each winter by SPAMS's. P6 is read at the half-cell-
corrected HyP3 position (T77). Folds: a season of year y from fits without y; a winter y/y+1 from fits without y and y+1.
Rule: ``08_deliverables/insar_full_test/preregistration/X-119.yaml``. Writes only to the hub:
``08_deliverables/insar_full_test/v4_head_to_head/``.

    make run E=v4_head_to_head_x119
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fusion_v3_wp1_x074 as W1  # noqa: E402
import fusion_v3_wp5_x078 as W5  # noqa: E402
import half_cell_p6_x118 as H  # noqa: E402
import p6_validation_x082 as X  # noqa: E402
import spams_baseline_x112 as SB  # noqa: E402

from insar_wetlands import fusion as fu  # noqa: E402
from insar_wetlands import fusion3 as f3  # noqa: E402
from insar_wetlands import spams  # noqa: E402
from insar_wetlands.cube import cube_root  # noqa: E402
from insar_wetlands.cube.core import WAVELENGTH_MM, rad_to_los_mm  # noqa: E402
from insar_wetlands.fusion2 import ou_from_series  # noqa: E402
from insar_wetlands.paths import hub_root  # noqa: E402
from insar_wetlands.scoring import (Protocol, assert_preregistered, chain, detrended_stats, guard_lockbox,  # noqa: E402
                                    level_views, winter_steps)

HUB = hub_root()
FT = HUB / "08_deliverables" / "insar_full_test"
OUT = FT / "v4_head_to_head"
PREREG = FT / "preregistration" / "X-119.yaml"
YEARS = (2021, 2022, 2023, 2024, 2025)
LAG = 6
ORDER = ["chain", "A: chain + SPAMS re-levelling", "B: v3+", "C: v3+ with SPAMS winter anchor"]


# ------------------------------------------------------------------------------ inputs

def pairs_corrected(p: pd.DataFrame) -> pd.DataFrame:
    """Consecutive P6 pairs (X-099's set), phase re-read at the corrected position (T77), v3's sign convention."""
    for c in ("t1", "t2"):
        p[c] = pd.to_datetime(p[c], format="ISO8601", utc=True)
    cons = X.load_pairs()
    keyset = set(zip(cons.track, cons.t1, cons.t2))
    p = p[[k in keyset for k in zip(p.track, p.t1, p.t2)]].copy()
    shift = pd.read_csv(FT / "hyp3_placement" / "decisive" / "verdict.csv").set_index("track")
    px = pd.read_csv(HUB / "08_deliverables" / "field_first" / "plot_pixels.csv").set_index("plot")
    xt, yt = float(px.loc["P6", "E"]), float(px.loc["P6", "N"])
    parts = []
    for track, g in p.groupby("track"):
        u = H.unw_at(track, xt - shift.loc[track, "median_dx_m"], yt - shift.loc[track, "median_dy_m"])
        key = g.t1.dt.strftime("%Y%m%d") + "_" + g.t2.dt.strftime("%Y%m%d")
        g = g.assign(y=-rad_to_los_mm(key.map(u).to_numpy() - g.ref_hard_unw.to_numpy()))
        parts.append(g)
    p = pd.concat(parts).sort_values(["track", "t1"])
    p["y_vert_mm"] = p.y / p.los_u
    return p[p.y.notna() & p.los_u.notna() & ~p.platforms.str.contains("S1D")]


def params(d: pd.DataFrame, excl: set) -> dict:
    dd = d[~d.year.isin(excl)]
    c = W1.changes(dd)
    ok = c.d_laser_cm.notna() & c.d_wtd_cm.notna()
    g = float(np.sum(c.d_laser_cm[ok] * c.d_wtd_cm[ok]) / np.sum(c.d_wtd_cm[ok] ** 2)) * 10
    res = []
    for _, gg in dd.groupby("year"):
        if gg.laser_cm.notna().sum() < 30:
            continue
        r = gg.laser_cm * 10 - g * gg.wtd_cm
        res.append(r - r.mean())
    rr = pd.concat(res).dropna()
    ou = ou_from_series((rr.index - rr.index[0]).days.to_numpy(float), rr.to_numpy(), max_lag_days=60)
    return {"g": g, "tau": ou["tau_days"], "sigma_r": ou["sigma"]}


def detect_model(p: pd.DataFrame, excl: set):
    q = p[p.usable & (p.revisit_days <= 24) & p.match.notna() & ~p.year.isin(excl)]
    X_ = q[W5.FEATS].to_numpy(float)
    med = np.nanmedian(X_, 0)
    m = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000)).fit(np.where(np.isfinite(X_), X_, med), q.match.astype(int))
    return m, med


def bucket_all(d: pd.DataFrame, excl: set) -> pd.Series:
    """The ERA5 bucket water table over the whole span, fitted on the logger of the training years."""
    rain, pet, w = d.rain_era5_mm.fillna(0).to_numpy(), d.pet_hamon.to_numpy(), d.wtd_cm.to_numpy()
    train = ~d.year.isin(excl).to_numpy()
    par = fu.fit_bucket(np.where(train, w, np.nan), rain, pet)
    floor = np.nanpercentile(w[train], 1)
    w0 = w[np.isfinite(w)][0]
    return pd.Series(fu.simulate_bucket(par, w0, rain, pet, floor=floor), index=d.index)


def spams_level(d: pd.DataFrame, excl: set) -> pd.Series:
    valid = (d.dh_mm.notna() & ~d.frozen_win & ~d.year.isin(excl)).to_numpy()
    P, E = d.rain_era5_mm.to_numpy(), d.et0_mm.to_numpy()
    par = spams.fit(P, E, d.dh_mm.to_numpy(), valid, LAG)
    return pd.Series(spams.reversible(P, E, par["x_p"], par["x_e"], par["tau"]), index=d.index)


# ------------------------------------------------------------------------------ one fold

def winter_bounds(nodes: pd.DatetimeIndex, prot: Protocol) -> list[tuple[int, int]]:
    m0, m1 = prot.season_months
    s0, s1 = prot.spring_window
    out = []
    for y in sorted(set(nodes.year)):
        a = np.flatnonzero((nodes.year == y) & (nodes.month >= m0) & (nodes.month <= m1))
        b = np.flatnonzero((nodes.year == y + 1) & (nodes.month >= s0) & (nodes.month <= s1))
        if len(a) and len(b):
            out.append((a[-1], b[0]))
    return out


def run_fold(d, p, noise, excl: set, prot: Protocol, variant: str, p_all: pd.DataFrame) -> pd.DataFrame:
    par = params(d, excl)
    dm, med = detect_model(p_all, excl)          # trained as v3: all P6 pairs up to 24 days
    W = bucket_all(d, excl)
    if variant == "B, logger water table":
        W = d.wtd_cm.combine_first(W)
    fit = {(r.rev_class, r.state): r for r in noise.itertuples()}
    q = p[p.state != "frozen/snow"].copy()
    Xf = q[W5.FEATS].to_numpy(float)
    q["p_detect"] = dm.predict_proba(np.where(np.isfinite(Xf), Xf, med))[:, 1]
    cls = np.where(q.revisit_days <= 7, "≤ 7 d", "12 d")
    q["sigma_noise_mm"] = [f3.crb_sigma_mm(c, fit[(k, s)].looks, fit[(k, s)].floor_mm) for c, k, s in zip(q.coh_3x3, cls, q.state)]
    q["sigma_eff_mm"] = (q.sigma_noise_mm if variant == "B, no detectability inflation"
                         else f3.detect_inflated_sigma(q.sigma_noise_mm, q.p_detect, miss_scale_mm=WAVELENGTH_MM / 4))
    q["sigma_vert_mm"] = q.sigma_eff_mm / q.los_u
    nodes = pd.DatetimeIndex(sorted(set(p.t1) | set(p.t2)))
    t_days = ((nodes - nodes[0]).total_seconds() / 86400).to_numpy()
    idx = {t: k for k, t in enumerate(nodes)}
    i, j = q.t1.map(idx).to_numpy(), q.t2.map(idx).to_numpy()
    prior = par["g"] * W.reindex(nodes.tz_convert(None).normalize()).to_numpy()
    if variant.startswith("C"):
        R = spams_level(d, excl).reindex(nodes.tz_convert(None).normalize()).to_numpy()
        offset = np.zeros(len(nodes))
        for a, b in winter_bounds(nodes, prot):
            offset[b:] += (R[b] - R[a]) - (prior[b] - prior[a])
        prior = prior + offset
    obs = q.y_vert_mm.to_numpy() - (prior[j] - prior[i])
    r, sd = f3.ou_increment_smoother(t_days, i, j, obs, q.sigma_vert_mm.to_numpy(), par["tau"], par["sigma_r"])
    return pd.DataFrame({"time": nodes, "h": prior + r, "sd": sd})


# ------------------------------------------------------------------------------ scoring

def score_series(c: pd.DataFrame, prot: Protocol, year: int | None = None, winter: str | None = None) -> tuple[list, pd.DataFrame]:
    """Season views of ``year`` (or all) and winter steps from a (time, seg, est, truth) table."""
    seasons = []
    for view, dd in level_views(c, prot):
        if view.startswith("season") and (year is None or view.startswith(f"season {year}")) and len(dd) >= 4:
            seasons.append({"view": view, **detrended_stats(dd.time, dd.truth, dd.est)})
    z = c.dropna(subset=["truth"]).copy()
    for _, cs in z.groupby("seg"):
        z.loc[cs.index, "truth"] -= cs.truth.iloc[0]
        z.loc[cs.index, "est"] -= cs.est.iloc[0]
    w = winter_steps(z, prot)
    if winter is not None and len(w):
        w = w[w.winter == winter]
    return seasons, w


def main() -> None:
    assert_preregistered(PREREG, [OUT / "verdict.csv"])
    OUT.mkdir(parents=True, exist_ok=True)
    prot = Protocol.load(FT / "protocol.yaml")
    d = SB.table()
    d = d[d.year.between(2020, 2025)]
    _, p_all, noise = W5.inputs()
    p = pairs_corrected(p_all.copy())
    p = p[p.t2.dt.year <= 2025]
    guard_lockbox(pd.concat([p.t1, p.t2]), prot)
    lev = X.laser()[1]
    truth_at = lambda t: lev.reindex(pd.DatetimeIndex(t).round("h")).to_numpy()
    winters = [f"{y}/{str(y + 1)[2:]}" for y in YEARS[:-1]]
    folds = {("season", y): {y} for y in YEARS} | {("winter", w): {int(w[:4]), int(w[:4]) + 1} for w in winters}
    variants = ["B: v3+", "C: v3+ with SPAMS winter anchor", "B, no detectability inflation", "B, logger water table"]
    runs = {(v, k): run_fold(d, p, noise, ex, prot, v.split(":")[0] if v.startswith("C") else v, p_all)
            for v in variants for k, ex in folds.items()}
    season_rows, winter_rows, cover = [], [], []
    for track in X.TRACKS:
        g = p[p.track == track].sort_values("t1")
        c0 = chain(g.assign(est_mm=g.y_vert_mm), "est_mm")          # vertical chain, its segments
        c0["truth"] = truth_at(c0.time)
        for y in YEARS:
            s, _ = score_series(c0, prot, year=y)
            season_rows += [{"track": track, "estimate": "chain", **r} for r in s]
        _, w = score_series(c0, prot)
        winter_rows += [{"track": track, "estimate": "chain", **r} for r in w.to_dict("records")]
        # A: the chain with each winter step set to SPAMS's (fit without the winter's two years)
        for wname in winters:
            ex = folds[("winter", wname)]
            R = spams_level(d, ex)
            cA = c0.copy()
            z = cA.dropna(subset=["truth"])
            yy = int(wname[:4])
            a_ = z[(z.time.dt.year == yy) & z.time.dt.month.between(*prot.season_months)]
            b_ = z[(z.time.dt.year == yy + 1) & z.time.dt.month.between(*prot.spring_window)]
            if len(a_) and len(b_):
                ia, ib = a_.index[-1], b_.index[0]
                step = R.reindex(pd.DatetimeIndex([cA.time[ib]]).tz_convert(None).normalize()).iloc[0] - \
                    R.reindex(pd.DatetimeIndex([cA.time[ia]]).tz_convert(None).normalize()).iloc[0]
                cA.loc[cA.index >= ib, "est"] += step - (cA.est[ib] - cA.est[ia])
                cA.loc[cA.index >= ib, "seg"] = cA.seg[ia]           # re-levelled: one segment across the winter
            _, w = score_series(cA, prot, winter=wname)
            winter_rows += [{"track": track, "estimate": ORDER[1], **r} for r in w.to_dict("records")]
        season_rows += [{**r, "estimate": ORDER[1]} for r in season_rows
                        if r["track"] == track and r["estimate"] == "chain"]  # A = the chain within seasons
        season_rows = [dict(r) for r in season_rows]
        for v in variants:
            for (kind, key), ex in folds.items():
                h = runs[(v, (kind, key))]
                hv = h.set_index("time").reindex(c0.time)
                c = pd.DataFrame({"time": c0.time, "seg": 0, "est": hv.h.to_numpy(), "truth": c0.truth.to_numpy()})
                if kind == "season":
                    s, _ = score_series(c, prot, year=key)
                    season_rows += [{"track": track, "estimate": v, **r} for r in s]
                    zz = c[c.time.dt.year == key].dropna()
                    if len(zz) > 3:
                        dev = (zz.est - zz.est.mean()) - (zz.truth - zz.truth.mean())
                        sdv = hv.sd.reindex(zz.time).to_numpy()
                        cover.append({"track": track, "estimate": v, "year": key,
                                      "within_1sd": float(np.mean(np.abs(dev) <= sdv)),
                                      "within_2sd": float(np.mean(np.abs(dev) <= 2 * sdv))})
                else:
                    _, w = score_series(c, prot, winter=key)
                    winter_rows += [{"track": track, "estimate": v, **r} for r in w.to_dict("records")]
    S, Wt, Cv = pd.DataFrame(season_rows).drop_duplicates(), pd.DataFrame(winter_rows), pd.DataFrame(cover)
    summ = (S.groupby(["estimate", "track"]).rmse_detrended_mm.median().rename("season_rmse_median_mm").to_frame()
            .join(Wt.assign(ae=Wt.difference_mm.abs()).groupby(["estimate", "track"]).ae.mean().rename("winter_mae_mm"))
            .join(Wt.groupby(["estimate", "track"]).size().rename("winter_steps")).reset_index())
    verdict = decide(summ)
    for name, df in (("seasons", S), ("winters", Wt), ("summary", summ), ("coverage", Cv), ("verdict", verdict)):
        df.to_csv(OUT / f"{name}.csv", index=False)
    (OUT / "README.md").write_text(readme(summ, verdict, Cv), encoding="utf-8")
    print(summ.to_string(index=False))
    print(verdict.to_string(index=False))


def decide(summ: pd.DataFrame) -> pd.DataFrame:
    s = summ.set_index(["estimate", "track"])
    rows = []
    for e in ORDER:
        if e not in s.index.get_level_values(0):
            continue
        elig = all(s.loc[(e, t), "season_rmse_median_mm"] <= 1.10 * s.loc[("chain", t), "season_rmse_median_mm"] for t in X.TRACKS)
        wm = float(np.mean([s.loc[(e, t), "winter_mae_mm"] for t in X.TRACKS]))
        rows.append({"estimate": e, "eligible": elig, "winter_mae_mean_mm": wm})
    v = pd.DataFrame(rows)
    chain_w = float(v.loc[v.estimate == "chain", "winter_mae_mean_mm"].iloc[0])
    el = v[v.eligible]
    best = el.winter_mae_mean_mm.min()
    tied = el[el.winter_mae_mean_mm <= 1.10 * best]
    pick = tied.iloc[0].estimate if len(tied) else "chain"          # ORDER is simplest first
    if best > 0.90 * chain_w:
        pick = "chain"
    v["selected"] = v.estimate == pick
    return v


def readme(summ, verdict, cov) -> str:
    md = lambda df: "\n".join(["| " + " | ".join(map(str, df.columns)) + " |", "|" + "---|" * len(df.columns)] +
                              ["| " + " | ".join(f"{x:.2f}" if isinstance(x, float) else str(x) for x in r) + " |"
                               for r in df.itertuples(index=False)])
    cv = cov.groupby(["estimate", "track"])[["within_1sd", "within_2sd"]].median().reset_index() if len(cov) else cov
    return f"""# X-119 — the v4 architecture head-to-head at P6 (generated)

Generated by `experiments/v4_head_to_head_x119.py` (`make run E=v4_head_to_head_x119`); do not edit. Pre-registration:
`../preregistration/X-119.yaml`. Vertical mm against the laser's vertical level; P6 at the half-cell-corrected HyP3
position (T77); consecutive 6- and 12-day pairs; hard-target reference. Season views from fits without that year;
winter steps from fits without the winter's two years. Selection: eligible if the within-season error is at most 10%
above the chain's on both tracks; the lowest mean winter error wins, ties (within 10%) go to the simpler; the chain
stays if nothing beats its winter error by 10%. Rows "B, no detectability inflation" and "B, logger water table" are
the ablation (secondary, not selectable).

## Verdict

{md(verdict)}

## Summary

{md(summ)}

## Posterior coverage of the smoother (median over held-out seasons)

{md(cv) if len(cv) else "—"}
"""


if __name__ == "__main__":
    main()

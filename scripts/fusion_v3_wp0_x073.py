"""X-073 — fusion v3 WP0: which phase reference, and how noisy is each pair? (exploratory, Mac-safe, reads the cube)

E0.1 Reference study. The P6 phase (−φ·λ/4π, uplift positive) referenced to: grassland, stable ground, the mat's own
     median (shown only to make the point that it removes the raft), the hard targets (built-up ≥ 0.5 and the track's
     median 12-day coherence ≥ 0.5), a sweep of those two thresholds, the near (village) and far clusters of targets,
     the k nearest built-up cells, every leave-one-target-out set and random subsets. Judged at P6 against the laser:
     correlation, slope, year-blocked residual spread (pairs), and chained levels within runs of consecutive pairs
     (r, RMSE, amplitude ratio — the quantities v2 was validated on). A laser-free check: the spread of stable ground
     after referencing.
E0.2 Noise model. Stable-ground pixels (not used as reference) after the chosen reference: robust spread binned by
     the pixel's coherence, per revisit and surface state, fitted with the Cramér–Rao form (effective looks + a floor);
     the reference's own floor from leave-one-target-out; platform pairs. Checked at P6: year-blocked residuals over the
     predicted σ — coverage of ±1σ and ±2σ, and the extra (non-motion) spread the model leaves.

    PYTHONPATH=src python scripts/fusion_v3_wp0_x073.py
Outputs only to the hub (field values): ``08_deliverables/fusion_v3/wp0/``.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import xarray as xr  # noqa: E402

from insar_wetlands import fusion3 as f3  # noqa: E402
from insar_wetlands.cube.core import rad_to_los_mm  # noqa: E402

HUB = Path(__file__).resolve().parents[3]
G, S = HUB / "06_data" / "cube" / "gold", HUB / "06_data" / "cube" / "silver"
OUT = HUB / "08_deliverables" / "fusion_v3" / "wp0"
P6 = (68, 66)
RNG = np.random.default_rng(73)
PERIODS = ("2020_2021", "2022_2024", "2025", "2026")
TRACKS = ("ascending", "descending")
NEAR_KM = 1.5
SEASON = {12: "DJF", 1: "DJF", 2: "DJF", 3: "MAM", 4: "MAM", 5: "MAM", 6: "JJA", 7: "JJA", 8: "JJA", 9: "SON",
          10: "SON", 11: "SON"}


def truthy(x: pd.Series) -> pd.Series:
    return x.map({True: True, False: False, "True": True, "False": False}).fillna(False).astype(bool)


def p6_truth() -> pd.DataFrame:
    d = pd.read_csv(G / "plots_pair.csv.gz", low_memory=False, parse_dates=["t1", "t2"])
    p = d[(d["plot"] == "P6") & (d.grid == "40m")].copy()
    p["frozen"] = truthy(p.flag_frozen_air_t1) | truthy(p.flag_frozen_air_t2)
    p["snow"] = (truthy(p.laser_snow_t1) | truthy(p.laser_snow_t2) | truthy(p.flag_snow_om_t1) |
                 truthy(p.flag_snow_om_t2))
    p["wet"] = truthy(p.flag_wet_any_t1) | truthy(p.flag_wet_any_t2)
    p["state"] = np.where(p.frozen | p.snow, "frozen/snow", np.where(p.wet, "wet", "dry"))
    p["usable"] = p.d_laser_los_mm.notna() & ~p.frozen & ~p.snow
    return p.set_index(["period", "track", "pair"])


def references(st: xr.Dataset, track: str) -> dict[str, np.ndarray]:
    z, built = st.zone.values, st.wc_built.values
    coh = st[f"coh12_{track[:3]}_median"].values
    rr, cc = np.indices(z.shape)
    dist_km = np.hypot(rr - P6[0], cc - P6[1]) * 0.04
    hard = (built >= 0.5) & (coh >= 0.5)
    refs = {"grassland": z == 3, "stable": st.stable.values == 1, "mat median (removes the raft)": z == 1,
            "hard targets": hard, "hard targets, near (≤ 1.5 km)": hard & (dist_km <= NEAR_KM),
            "hard targets, far (> 1.5 km)": hard & (dist_km > NEAR_KM)}
    for b in (0.25, 0.5, 0.75):
        for c in (0.4, 0.5, 0.6):
            refs[f"sweep built≥{b} coh≥{c}"] = (built >= b) & (coh >= c)
    for k in (1, 3, 5, 10):
        refs[f"nearest {k} built-up cells"] = f3.nearest_cells(built >= 0.5, *P6, k)
    return {k: v for k, v in refs.items() if v.sum() > 0}


def score(d: pd.DataFrame) -> list[dict]:
    """Metrics at P6 per revisit (6, 12) for a referenced series pooled over periods (columns y, d_laser_los_mm,
    usable, revisit_days, t1, t2, period). The slope is year-blocked over the pooled years; chained levels are built
    within runs of consecutive usable pairs. Rows for each period come from the same pooled residuals."""
    out = []
    for rv in (6, 12):
        sub = d[(d.revisit_days == rv) & d.usable.fillna(False).astype(bool) & np.isfinite(d.y)].copy()
        if len(sub) < 8:
            continue
        x, yy = sub.d_laser_los_mm.to_numpy(), sub.y.to_numpy()
        sub["res"] = f3.year_blocked_residual(x, yy, sub.t1.dt.year.to_numpy())
        run, (ly, lx) = f3.chain_runs(pd.DatetimeIndex(sub.t1), pd.DatetimeIndex(sub.t2), yy, x, min_len=4)
        sub["ly"], sub["lx"] = ly, lx
        for per, g in [("all", sub)] + list(sub.groupby("period")):
            if len(g) < 5:
                continue
            gx, gy = g.d_laser_los_mm.to_numpy(), g.y.to_numpy()
            out.append({"period": per, "revisit_days": rv, "n": len(g), "r": float(np.corrcoef(gx, gy)[0, 1]),
                        "m": f3.slope_origin(gx, gy), "resid_sd_year_blocked_mm": float(np.nanstd(g.res)),
                        **f3.level_metrics(g.ly.to_numpy(), g.lx.to_numpy())})
    return out


def stable_spread(U: np.ndarray, ref_vals: np.ndarray, stable: np.ndarray, rev: np.ndarray) -> dict:
    dev = rad_to_los_mm(U[:, stable] - ref_vals[:, None])
    sd = np.array([f3.robust_sd(-r) for r in dev])
    return {f"stable_spread_{rv}d_mm": float(np.nanmedian(sd[rev == rv])) for rv in (6, 12) if (rev == rv).any()}


def series(truth: pd.DataFrame, st: xr.Dataset, track: str, mask_fn) -> pd.DataFrame:
    """The P6 referenced phase over every period of a track, with the truth columns (mask_fn(refs) → mask)."""
    parts = []
    refs = references(st, track)
    for period in PERIODS:
        ds = stack(period, track)
        if ds is None:
            continue
        U = ds.unw
        pairs = ds.pair
        mask = mask_fn(refs, U)
        if mask is None or not mask.any():
            continue
        tr = truth.reindex(pd.MultiIndex.from_arrays([[period] * len(pairs), [track] * len(pairs), pairs])).reset_index()
        tr["period"], tr["track"], tr["pair"] = period, track, pairs
        tr["revisit_days"] = ds.revisit_days
        tr["y"] = -rad_to_los_mm(U[:, P6[0], P6[1]] - f3.reference_values(U, mask))
        parts.append(tr)
    return pd.concat(parts, ignore_index=True)


_STACKS: dict = {}


class Stack:
    """One interferogram stack held in memory (read once; the files are not reopened per reference)."""

    def __init__(self, f: Path):
        with xr.open_dataset(f) as ds:
            self.unw = ds.unw.values.astype("float32")
            self.coh = ds.coh.values.astype("float32")
            self.pair = ds.pair.values.astype(str)
            self.revisit_days = ds.revisit_days.values
            self.t1, self.t2 = pd.DatetimeIndex(ds.t1.values), pd.DatetimeIndex(ds.t2.values)
            self.platforms = ds.platforms.values.astype(str)


def stack(period: str, track: str):
    if (period, track) not in _STACKS:
        f = S / f"ifg_{period}_{track}.nc"
        _STACKS[(period, track)] = Stack(f) if f.exists() else None
    return _STACKS[(period, track)]


def e01(truth: pd.DataFrame, st: xr.Dataset) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    rows, jk, scat = [], [], {}
    stable = st.stable.values == 1
    for track in TRACKS:
        refs = references(st, track)
        for name, mask in refs.items():
            d = series(truth, st, track, lambda r, U, n=name: r[n])
            spread = {}
            for period in PERIODS:
                ds = stack(period, track)
                if ds is None:
                    continue
                U = ds.unw
                rv = f3.reference_values(U, mask)
                for k, v in stable_spread(U, rv, stable & ~mask, ds.revisit_days).items():
                    spread.setdefault(k, []).append(v)
            base = {"track": track, "reference": name, "cells": int(mask.sum()),
                    **{k: float(np.nanmedian(v)) for k, v in spread.items()}}
            for sc in score(d):
                rows.append({**base, **sc})
            if name in ("grassland", "stable", "hard targets"):
                m = (d.revisit_days == 6) & d.usable.fillna(False).astype(bool)
                scat[(track, name)] = [(d.d_laser_los_mm[m].to_numpy(), d.y[m].to_numpy())]
        hard = refs["hard targets"]
        cells = np.argwhere(hard)
        for i, (r, c) in enumerate(cells):
            m2 = hard.copy()
            m2[r, c] = False
            for sc in score(series(truth, st, track, lambda rf, U, m2=m2: m2)):
                if sc["period"] == "all":
                    jk.append({"track": track, "set": "leave one out", "draw": i, **sc})
        for k in (3, 5, 10):
            for draw in range(40):
                pick = cells[RNG.choice(len(cells), min(k, len(cells)), replace=False)]
                m2 = np.zeros_like(hard)
                m2[pick[:, 0], pick[:, 1]] = True
                for sc in score(series(truth, st, track, lambda rf, U, m2=m2: m2)):
                    if sc["period"] == "all":
                        jk.append({"track": track, "set": f"random {k}", "draw": draw, **sc})
    return pd.DataFrame(rows), pd.DataFrame(jk), scat


def e02(truth: pd.DataFrame, st: xr.Dataset, per_pair_sample: int = 600):
    """Noise of stable ground after the hard-target reference, binned by pixel coherence, per revisit and state."""
    samples, floors, plat = [], [], []
    for period in PERIODS:
        for track in TRACKS:
            ds = stack(period, track)
            if ds is None:
                continue
            U, C = ds.unw, ds.coh
            pairs, rev = ds.pair, ds.revisit_days
            hard = references(st, track)["hard targets"]
            stable = (st.stable.values == 1) & ~hard
            rv = f3.reference_values(U, hard)
            idx = pd.MultiIndex.from_arrays([[period] * len(pairs), [track] * len(pairs), pairs])
            state = truth.reindex(idx).state.fillna("dry").to_numpy()
            mid = ds.t1 + (ds.t2 - ds.t1) / 2
            st_idx = np.flatnonzero(stable.ravel())
            hc = np.argwhere(hard)
            for k in range(len(pairs)):
                pick = RNG.choice(st_idx, min(per_pair_sample, len(st_idx)), replace=False)
                e = -rad_to_los_mm(U[k].ravel()[pick] - rv[k])
                g = C[k].ravel()[pick]
                samples.append(pd.DataFrame({"period": period, "track": track, "revisit_days": rev[k],
                                             "state": state[k], "season": SEASON[mid[k].month],
                                             "platforms": ds.platforms[k], "coh": g, "dev_mm": e}))
                # the reference's own floor: each target against the median of the others
                vals = U[k][hard]
                if np.isfinite(vals).sum() > 3:
                    loo = [vals[i] - np.nanmedian(np.delete(vals, i)) for i in range(len(vals))]
                    floors.append({"period": period, "track": track, "revisit_days": rev[k],
                                   "loo_target_sd_mm": f3.robust_sd(-rad_to_los_mm(np.array(loo)))})
            del hc
    smp = pd.concat(samples, ignore_index=True)
    smp = smp[np.isfinite(smp.coh) & np.isfinite(smp.dev_mm)]
    smp["rev_class"] = pd.cut(smp.revisit_days, [0, 7, 13, 25, 1000], labels=["≤ 7 d", "12 d", "18–24 d", "≥ 36 d"])
    smp["coh_bin"] = (np.floor(smp.coh / 0.05) * 0.05 + 0.025).round(3)
    binned = (smp.groupby(["rev_class", "state", "coh_bin"], observed=True).dev_mm
              .agg(sd_mm=f3.robust_sd, n="size").reset_index())
    binned = binned[binned.n >= 200]
    fits = []
    for (rc, stt), g in binned.groupby(["rev_class", "state"], observed=True):
        if len(g) >= 5:
            fits.append({"rev_class": rc, "state": stt, "bins": len(g), **f3.fit_noise(g.coh_bin.to_numpy(),
                                                                                      g.sd_mm.to_numpy(), g.n.to_numpy())})
    plat = (smp[smp.rev_class == "≤ 7 d"].groupby(["track", "platforms"]).dev_mm
            .agg(sd_mm=f3.robust_sd, n="size").reset_index())
    season = smp.groupby(["rev_class", "season"], observed=True).dev_mm.agg(sd_mm=f3.robust_sd, n="size").reset_index()
    fl = pd.DataFrame(floors).groupby(["track", "revisit_days"]).loo_target_sd_mm.median().reset_index()
    fl = fl[fl.revisit_days <= 48]
    return binned, pd.DataFrame(fits), plat, season, fl


def validate(truth: pd.DataFrame, st: xr.Dataset, fits: pd.DataFrame) -> pd.DataFrame:
    """P6, hard-target reference, pooled over periods: year-blocked residual over the noise model's σ."""
    rows = []
    fit = {(r.rev_class, r.state): r for r in fits.itertuples()}
    for track in TRACKS:
        d = series(truth, st, track, lambda r, U: r["hard targets"])
        for rv, rc in ((6, "≤ 7 d"), (12, "12 d")):
            sub = d[(d.revisit_days == rv) & d.usable.fillna(False).astype(bool) & np.isfinite(d.y)].copy()
            if len(sub) < 8:
                continue
            sub["res"] = f3.year_blocked_residual(sub.d_laser_los_mm.to_numpy(), sub.y.to_numpy(), sub.t1.dt.year.to_numpy())
            sub["sig"] = [f3.crb_sigma_mm(c, fit[(rc, s)].looks, fit[(rc, s)].floor_mm) if (rc, s) in fit else np.nan
                          for c, s in zip(sub.coh_3x3, sub.state)]
            for per, g in [("all", sub)] + list(sub.groupby("period")):
                if len(g) < 5:
                    continue
                res, sig = g.res.to_numpy(), g.sig.to_numpy()
                extra = np.sqrt(max(np.nanvar(res) - np.nanmean(sig ** 2), 0)) if np.isfinite(res).sum() > 2 else np.nan
                rows.append({"track": track, "revisit_days": rv, "period": per, "median_sigma_mm": float(np.nanmedian(sig)),
                             "resid_sd_mm": float(np.nanstd(res)), "extra_non_noise_sd_mm": float(extra),
                             **f3.coverage(res / sig)})
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------------ report

def md(df: pd.DataFrame, nd: int = 2) -> str:
    d = df.copy()
    for c in d.columns:
        if d[c].dtype.kind == "f":
            d[c] = d[c].round(nd)
    lines = ["| " + " | ".join(map(str, d.columns)) + " |", "|" + "---|" * len(d.columns)]
    for r in d.itertuples(index=False):
        lines.append("| " + " | ".join("" if (isinstance(v, float) and np.isnan(v)) else str(v) for v in r) + " |")
    return "\n".join(lines)


def figures(ref: pd.DataFrame, scat: dict, binned: pd.DataFrame, fits: pd.DataFrame):
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.6), sharex=True, sharey=True)
    lim = 30
    for i, track in enumerate(TRACKS):
        for name, col in (("grassland", "C2"), ("stable", "0.55"), ("hard targets", "C3")):
            for x, y in scat.get((track, name), []):
                ax[i].plot(x, y, ".", color=col, ms=5, alpha=0.75, label=name)
        ax[i].plot([-lim, lim], [-lim, lim], "k-", lw=0.6)
        ax[i].set_xlim(-lim, lim)
        ax[i].set_ylim(-lim, lim)
        ax[i].set_title(f"P6, {track}, 6-day pairs (2020–2021, 2025)", fontsize=9)
        ax[i].set_xlabel("laser change, line of sight (mm)")
        h, lab = ax[i].get_legend_handles_labels()
        ax[i].legend(dict(zip(lab, h)).values(), dict(zip(lab, h)).keys(), fontsize=8)
    ax[0].set_ylabel("radar phase change at P6 (mm, uplift positive)")
    fig.tight_layout()
    fig.savefig(OUT / "fig_reference_scatter.png", dpi=150)
    plt.close(fig)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
    for i, (rc, g) in enumerate(binned[binned.state != "frozen/snow"].groupby("rev_class", observed=True)):
        if i > 1:
            break
        for stt, gg in g.groupby("state"):
            col = "C0" if stt == "dry" else "C1"
            ax[i].plot(gg.coh_bin, gg.sd_mm, "o", color=col, ms=4, label=f"{stt} (binned)")
            f = fits[(fits.rev_class == rc) & (fits.state == stt)]
            if len(f):
                c = np.linspace(0.1, 0.98, 100)
                ax[i].plot(c, f3.crb_sigma_mm(c, f.looks.iloc[0], f.floor_mm.iloc[0]), "-", color=col,
                           label=f"{stt}: looks {f.looks.iloc[0]:.0f}, floor {f.floor_mm.iloc[0]:.1f} mm")
        ax[i].set_title(f"stable ground after the hard-target reference, {rc}", fontsize=9)
        ax[i].set_xlabel("pixel coherence")
        ax[i].legend(fontsize=7)
    ax[0].set_ylabel("robust spread of the phase (mm, line of sight)")
    fig.tight_layout()
    fig.savefig(OUT / "fig_noise_model.png", dpi=150)
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    st = xr.open_dataset(S / "static_40m.nc")
    truth = p6_truth()
    ref, jk, scat = e01(truth, st)
    binned, fits, plat, season, floors = e02(truth, st)
    val = validate(truth, st, fits)
    for n, d in (("reference_study", ref), ("hard_target_robustness", jk), ("noise_binned", binned),
                 ("noise_model_fits", fits), ("noise_by_platform_pair", plat), ("noise_by_season", season),
                 ("reference_floor_leave_one_target_out", floors), ("noise_model_check_p6", val)):
        d.to_csv(OUT / f"{n}.csv", index=False)
    figures(ref, scat, binned, fits)
    key = ["grassland", "stable", "mat median (removes the raft)", "hard targets", "hard targets, near (≤ 1.5 km)",
           "hard targets, far (> 1.5 km)"]
    core = ref[ref.reference.isin(key) & ref.revisit_days.notna()]
    cols = ["track", "revisit_days", "period", "reference", "cells", "n", "r", "m", "resid_sd_year_blocked_mm",
            "n_levels", "r_level", "rmse_level_mm", "amplitude_ratio", "stable_spread_6d_mm", "stable_spread_12d_mm"]
    sweep = ref[~ref.reference.isin(key) & (ref.revisit_days == 6) & (ref.period == "all")]
    jks = (jk.groupby(["track", "revisit_days", "set"]).resid_sd_year_blocked_mm
           .describe(percentiles=[0.05, 0.5, 0.95])[["count", "5%", "50%", "95%"]].reset_index())
    # decision rule, fixed before looking: among the a-priori candidates (grassland, stable ground, hard targets at the
    # thresholds set in the cube), the lowest year-blocked residual for 6-day pairs pooled over the laser years, mean of
    # the two tracks. The mat median is excluded (it removes the motion it is meant to measure); near/far clusters,
    # the threshold sweep and the nearest cells are sensitivity tests, not candidates (choosing among them on the
    # laser would fit the reference to the one validation pixel).
    candidates = ["grassland", "stable", "hard targets"]
    six = (core[(core.revisit_days == 6) & (core.period == "all") & core.reference.isin(candidates)]
           .groupby("reference").resid_sd_year_blocked_mm.mean().sort_values())
    noisy = plat[(plat.n >= 300) & (plat.sd_mm > 2 * plat.sd_mm.median())]
    txt = ["# X-073 — fusion v3 WP0: phase reference and per-pair noise (exploratory)", "",
           "Generated by `05_code/SAr-adapted/scripts/fusion_v3_wp0_x073.py` from `06_data/cube` — every number below "
           "is computed. P6 phase −φ·λ/4π (uplift positive) against the laser's line-of-sight change; usable = laser "
           "present, not frozen, no snow. Year-blocked residual: slope fitted without the pair's year. Chained levels: "
           "within runs of consecutive usable pairs (≥ 4), centred per run — as v2 was validated.", "",
           "## E0.1 Reference study", "",
           md(core[[c for c in cols if c in core.columns]].sort_values(["track", "revisit_days", "period", "reference"])), "",
           "Year-blocked residual, 6-day pairs pooled over the laser years, mean of the two tracks (lower is better):", "",
           md(six.reset_index().rename(columns={"resid_sd_year_blocked_mm": "mean_resid_sd_mm"})), "",
           f"**Chosen reference: {six.index[0]}** (a-priori candidates only; see the rule in the script).", "",
           "Threshold sweep and nearest built-up cells (6-day):", "",
           md(sweep[["track", "reference", "cells", "n", "r", "resid_sd_year_blocked_mm", "r_level",
                     "rmse_level_mm", "amplitude_ratio"]]), "",
           "Robustness of the hard-target set (year-blocked residual SD, mm): leave one target out, random subsets:", "",
           md(jks), "", "![P6 against the laser by reference](fig_reference_scatter.png)", "",
           "## E0.2 Noise model (stable ground after the hard-target reference)", "",
           "Cramér–Rao form σ(γ) = √(floor² + (λ/4π)²·(1−γ²)/(2·L·γ²)) fitted to robust spreads binned by pixel coherence:", "",
           md(fits), "", "Platform pairs more than twice the median spread (≥ 300 samples) — exclude or down-weight:", "",
           md(noisy) if len(noisy) else "(none)", "", "The reference's own floor (each hard target against the median of the others, median over pairs):", "",
           md(floors), "", "By platform pair (≤ 7-day pairs):", "", md(plat), "", "By season:", "", md(season), "",
           "![noise model](fig_noise_model.png)", "",
           "## Check at P6: does the noise model explain the residual?", "",
           "Year-blocked residual divided by the predicted σ (P6 3×3 coherence, the pair's revisit and surface state). "
           "If σ were complete, ±1σ would hold 68 % and ±2σ 95 %; `extra_non_noise_sd_mm` is the spread the stable-ground "
           "noise does not account for — the mat's own non-motion signal, the part WP2 must model.", "",
           md(val), ""]
    (OUT / "README.md").write_text("\n".join(txt))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()

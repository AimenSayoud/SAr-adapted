"""X-146 — does dry-snow delay or freeze–thaw explain the chain's winter under-reading? (R5 M1, M2, m8)

Two parts. (1) At P6, each winter step's error (X-142's plain chain, hard targets) is split into the pairs the laser
judges (snow-free, valid laser change) and the stretches it does not, each stretch classed by snow and frost at its
dates. (2) Radar only, every winter 2017–2026: the grassland's referenced phase against the change of snow depth on
dry-snow pairs, set against the Guneriussen et al. (2001) dry-snow delay; and pairs spanning a freeze or thaw.

Pre-registration: hub ``08_deliverables/insar_full_test/preregistration/X-146.yaml`` (checked before the run).
Reads the hub and the cube; writes only to the hub: ``08_deliverables/insar_full_test/winter_snow_x146/``. No field
value in this file.

    make run E=winter_snow_x146
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import p6_validation_x082 as X  # noqa: E402  (the X-082 laser reader, unchanged)

from insar_wetlands.cube.core import rad_to_los_mm  # noqa: E402
from insar_wetlands.paths import cube_root, hub_root  # noqa: E402
from insar_wetlands.scoring import assert_preregistered  # noqa: E402

HUB, CUBE = hub_root(), cube_root()
FT = HUB / "08_deliverables" / "insar_full_test"
OUT = FT / "winter_snow_x146"
PREREG = FT / "preregistration" / "X-146.yaml"
X142 = FT / "winter_cycles_x142"
PERIODS = ("2017_2019", "2020_2021", "2022_2024", "2025", "2026")
TRACKS = ("ascending", "descending")
WINTER_MONTHS = (11, 12, 1, 2, 3, 4)
QUARTER_WAVE_MM = 55.465763 / 4
N_BOOT, SEED = 2000, 146
RHO = (0.15, 0.35)                                    # dry-snow density range, g cm⁻³


def delay_slope(rho: float, theta_deg: float) -> float:
    """Uplift read per mm of dry snow depth (negative: the delay reads as sinking), Guneriussen et al. (2001)."""
    eps = 1 + 1.6 * rho + 1.86 * rho ** 3
    th = np.radians(theta_deg)
    return -float(np.sqrt(eps - np.sin(th) ** 2) - np.cos(th))


# ------------------------------------------------------------------------------------------- inputs

def acquisitions(track: str) -> pd.DataFrame:
    """Site conditions at each acquisition of a track, by date (snow depth in mm)."""
    a = xr.open_dataset(CUBE / "gold" / f"mat_acq_{track}.nc")
    d = pd.DataFrame({"time": pd.to_datetime(a.time.values, utc=True),
                      "platform": a.platform.values.astype(str),
                      "h_mm": a.site_om_era5_land_snow_depth.values * 1000.0,
                      "h_imgw_mm": a.site_imgw_snow_depth_cm.values * 10.0,
                      "air_c": a.site_air_c.values,
                      "frozen_soil": a.site_flag_frozen_soil_om.values == 1,
                      "frozen_air": a.site_flag_frozen_air.values == 1,
                      "snow": (a.site_flag_snow_om.values == 1) | (a.site_flag_snow_imgw.values == 1)})
    d["date"] = d.time.dt.normalize()
    return d.drop_duplicates("date").set_index("date")


def grassland_incidence(track: str) -> float:
    d = xr.open_dataset(CUBE / "gold" / "mat_pairs_2022_2024_ascending.nc")
    col = "st_incidence_asc_deg" if track == "ascending" else "st_incidence_des_deg"
    return float(np.nanmedian(d[col].values[d.st_zone.values == 3]))


def pairs() -> pd.DataFrame:
    """Consecutive six- and twelve-day winter pairs at P6 with the scene references and the site conditions."""
    cols = ["plot", "grid", "period", "track", "t1", "t2", "revisit_days", "platforms",
            "ref_grassland_unw", "ref_mat_unw", "ref_hard_unw"]
    q = pd.read_csv(CUBE / "gold" / "plots_pair.csv.gz", usecols=cols, low_memory=False)
    q = q[(q["plot"] == "P6") & (q.grid == "40m") & q.period.isin(PERIODS) & q.revisit_days.isin([6, 12])
          & ~q.platforms.str.contains("S1D")].drop(columns=["plot", "grid"])
    for c in ("t1", "t2"):
        q[c] = pd.to_datetime(q[c], format="ISO8601", utc=True)
    q = q.drop_duplicates(["track", "t1", "t2"])
    out = []
    for track in TRACKS:
        a = acquisitions(track)
        days = np.sort(a.index[~a.platform.str.contains("S1D")].to_numpy())
        g = q[q.track == track].copy()
        d1, d2 = g.t1.dt.normalize().to_numpy(), g.t2.dt.normalize().to_numpy()
        between = np.searchsorted(days, d2, side="left") - np.searchsorted(days, d1, side="right")
        g = g[between == 0]
        for end in ("t1", "t2"):
            s = a.reindex(g[end].dt.normalize()).reset_index(drop=True)
            for c in ("h_mm", "h_imgw_mm", "air_c", "frozen_soil", "frozen_air", "snow"):
                g[f"{c}_{end}"] = s[c].to_numpy()
        out.append(g)
    p = pd.concat(out, ignore_index=True)
    p = p[p.t1.dt.month.isin(WINTER_MONTHS)].copy()
    p["month"] = p.t1.dt.strftime("%Y-%m")
    p["grass_mm"] = -rad_to_los_mm(p.ref_grassland_unw - p.ref_hard_unw)
    p["mat_mm"] = -rad_to_los_mm(p.ref_mat_unw - p.ref_hard_unw)
    p["dh_mm"] = p.h_mm_t2 - p.h_mm_t1
    p["dh_imgw_mm"] = p.h_imgw_mm_t2 - p.h_imgw_mm_t1
    any_snow = (p.h_mm_t1 > 0) | (p.h_mm_t2 > 0)
    p["snow_pair"] = any_snow
    p["dry_snow"] = any_snow & (p.air_c_t1 <= 0) & (p.air_c_t2 <= 0)
    p["transition"] = p.frozen_soil_t1.astype(bool) != p.frozen_soil_t2.astype(bool)
    return p


# ------------------------------------------------------------------------------------------- part 2: radar only

def on_track(p: pd.DataFrame, track: str) -> pd.Series:
    return pd.Series(True, index=p.index) if track == "both" else p.track == track


def ols(x, y) -> float:
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    return float(np.polyfit(x[ok], y[ok], 1)[0]) if ok.sum() > 3 and np.std(x[ok]) > 0 else np.nan


def month_boot(g: pd.DataFrame, stat, rng) -> tuple[float, float]:
    pos = {m: np.flatnonzero(g.month.to_numpy() == m) for m in np.unique(g.month)}
    keys = list(pos)
    draws = []
    for _ in range(N_BOOT):
        idx = np.concatenate([pos[keys[j]] for j in rng.integers(0, len(keys), len(keys))])
        draws.append(stat(g.iloc[idx]))
    lo, hi = np.nanpercentile(draws, [2.5, 97.5])
    return float(lo), float(hi)


def snow_slopes(p: pd.DataFrame, theta: dict[str, float], rng) -> pd.DataFrame:
    rows = []
    sets = {"dry snow, ERA5-Land depth (primary)": ("dry_snow", "dh_mm"),
            "dry snow, IMGW depth": ("dry_snow", "dh_imgw_mm"),
            "all snow pairs, ERA5-Land depth": ("snow_pair", "dh_mm")}
    for label, (sel, xcol) in sets.items():
        for zone in ("grass_mm", "mat_mm"):
            for track in ("both",) + TRACKS:
                g = p[p[sel] & on_track(p, track)].dropna(subset=[xcol, zone])
                th = np.mean(list(theta.values())) if track == "both" else theta[track]
                lo_m, hi_m = sorted((delay_slope(RHO[1], th), delay_slope(RHO[0], th)))
                b = ols(g[xcol], g[zone])
                lo, hi = month_boot(g, lambda d: ols(d[xcol], d[zone]), rng) if len(g) > 5 else (np.nan, np.nan)
                reading = ("visible, consistent with the model" if hi < 0 and lo <= hi_m and hi >= lo_m else
                           "visible" if hi < 0 else "not detected")
                rows.append({"set": label, "zone": zone.replace("_mm", ""), "track": track, "pairs": len(g),
                             "months": g.month.nunique(), "slope": b, "lo95": lo, "hi95": hi,
                             "model_lo": lo_m, "model_hi": hi_m, "incidence_deg": th, "reading": reading})
    return pd.DataFrame(rows)


def freeze_thaw(p: pd.DataFrame, rng) -> pd.DataFrame:
    rows = []
    for zone in ("grass_mm", "mat_mm"):
        for track in ("both",) + TRACKS:
            g = p[on_track(p, track)].dropna(subset=[zone])

            def diff(d):
                a, b = d.loc[d.transition, zone], d.loc[~d.transition, zone]
                return float(a.mean() - b.mean()) if len(a) and len(b) else np.nan
            lo, hi = month_boot(g, diff, rng)
            rows.append({"zone": zone.replace("_mm", ""), "track": track, "transition_pairs": int(g.transition.sum()),
                         "other_winter_pairs": int((~g.transition).sum()), "mean_difference_mm": diff(g),
                         "lo95": lo, "hi95": hi,
                         "reading": "a step at freeze or thaw" if lo > 0 or hi < 0 else "no step detected"})
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------------------------- part 1: P6

def attribution(slopes: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, float]:
    pk = pd.read_csv(X142 / "pairs_k.csv")
    for c in ("t1", "t2"):
        pk[c] = pd.to_datetime(pk[c], format="ISO8601", utc=True)
    w = pd.read_csv(X142 / "winters.csv")
    w = w[(w.estimate == "plain chain") & (w.reference == "hard_targets")].rename(columns={"from": "t_from"})
    las = X.laser()[0]
    rows, runs, worst = [], [], 0.0
    for r in w.itertuples():
        a = acquisitions(r.track)
        lt = las[r.track].dropna()
        t0, t1 = pd.Timestamp(r.t_from, tz="UTC"), pd.Timestamp(r.to, tz="UTC") + pd.Timedelta(days=1)
        g = pk[(pk.track == r.track) & (pk.t1 >= t0) & (pk.t2 <= t1)].sort_values("t1")
        worst = max(worst, abs(g.s1_hard_targets.sum() - r.insar_change_mm))
        e = r.insar_change_mm - r.laser_change_mm
        valid = g.d_laser_los_mm.notna()
        e_free = float((g.s1_hard_targets - g.d_laser_los_mm)[valid].sum())
        amb = valid & (g.d_laser_los_mm.abs() > QUARTER_WAVE_MM)
        e_amb = float((g.s1_hard_targets - g.d_laser_los_mm)[amb].sum())
        # stretches of pairs without a valid laser change, each closed by valid laser readings at its ends
        cls_sum = {"snow": 0.0, "frost only": 0.0, "neither": 0.0, "unresolved": 0.0}
        g = g.assign(run=(valid != valid.shift()).cumsum())
        for _, s in g[~valid].groupby("run"):
            ends = []
            for t in (s.t1.iloc[0], s.t2.iloc[-1]):
                k = lt.index.get_indexer([t], method="nearest", tolerance=pd.Timedelta(hours=1))[0]
                ends.append(lt.iloc[k] if k >= 0 else np.nan)
            d = a.reindex(pd.concat([s.t1, s.t2]).dt.normalize())
            cls = "snow" if d.snow.any() else "frost only" if (d.frozen_soil | d.frozen_air).any() else "neither"
            err = float(s.s1_hard_targets.sum() - (ends[1] - ends[0])) if np.all(np.isfinite(ends)) else np.nan
            key = cls if np.isfinite(err) else "unresolved"
            cls_sum[key] += err if np.isfinite(err) else 0.0
            runs.append({"track": r.track, "winter": r.winter, "from": s.t1.iloc[0].date(), "to": s.t2.iloc[-1].date(),
                         "pairs": len(s), "class": cls, "radar_mm": float(s.s1_hard_targets.sum()),
                         "laser_mm": float(ends[1] - ends[0]) if np.isfinite(err) else np.nan, "error_mm": err,
                         "max_snow_depth_mm": float(d.h_mm.max())})
        e_rest = e - e_free
        sl = slopes[(slopes.set.str.contains("primary")) & (slopes.zone == "grass") & (slopes.track == r.track)].iloc[0]
        steep = min(sl.model_lo, sl.lo95) if np.isfinite(sl.lo95) else sl.model_lo
        hmax = float(a.loc[(a.index >= t0.normalize()) & (a.index < t1.normalize()), "h_mm"].max())
        rows.append({"track": r.track, "winter": r.winter, "E_mm": e, "E_free_mm": e_free, "E_rest_mm": e_rest,
                     "E_free_beyond_quarter_wave_mm": e_amb, "pairs": len(g), "pairs_valid": int(valid.sum()),
                     "rest_snow_mm": cls_sum["snow"], "rest_frost_only_mm": cls_sum["frost only"],
                     "rest_neither_mm": cls_sum["neither"],
                     "rest_unresolved_mm": e_rest - cls_sum["snow"] - cls_sum["frost only"] - cls_sum["neither"],
                     "max_snow_depth_mm": hmax, "snow_bound_mm": steep * hmax})
    return pd.DataFrame(rows), pd.DataFrame(runs), worst


def verdict(att: pd.DataFrame, slopes: pd.DataFrame) -> dict:
    se, sf, sr = att.E_mm.sum(), att.E_free_mm.sum(), att.E_rest_mm.sum()
    if sr < 0 and sr / se >= 2 / 3:
        a = "the under-reading is carried by the pairs that touch snow or frost"
    elif sf / se >= 2 / 3:
        a = "the under-reading is carried by the snow-free pairs the laser judges"
    else:
        a = "shared"
    prim = slopes[(slopes.set.str.contains("primary")) & (slopes.zone == "grass") & (slopes.track == "both")].iloc[0]
    return {"attribution": a, "sum_E_mm": se, "sum_E_free_mm": sf, "sum_E_rest_mm": sr,
            "share_rest": sr / se, "snow_delay": prim.reading, "snow_slope": prim.slope,
            "snow_slope_ci": [prim.lo95, prim.hi95], "model_range": [prim.model_lo, prim.model_hi],
            "dry_snow_pairs": int(prim.pairs)}


def md(df: pd.DataFrame) -> str:
    d = df.copy()
    for c in d.columns:
        if pd.api.types.is_float_dtype(d[c]):
            d[c] = d[c].round(3)
    head = "| " + " | ".join(map(str, d.columns)) + " |\n|" + "---|" * len(d.columns) + "\n"
    return head + "\n".join("| " + " | ".join(map(str, row)) + " |" for row in d.itertuples(index=False))


def main() -> None:
    assert_preregistered(PREREG, [OUT / "verdict.json"])
    rng = np.random.default_rng(SEED)
    p = pairs()
    pk = pd.read_csv(X142 / "pairs_k.csv")
    for c in ("t1", "t2"):
        pk[c] = pd.to_datetime(pk[c], format="ISO8601", utc=True)
    m = pk.merge(p[["track", "t1", "t2", "ref_hard_unw"]], on=["track", "t1", "t2"], suffixes=("", "_cube"))
    ref_diff = float((m.ref_hard_unw - m.ref_hard_unw_cube).abs().max()) if len(m) else np.nan
    theta = {t: grassland_incidence(t) for t in TRACKS}
    slopes = snow_slopes(p, theta, rng)
    ft = freeze_thaw(p, rng)
    att, runs, worst = attribution(slopes)
    if not worst < 1e-6:
        raise SystemExit(f"stop: X-142's winter steps not reproduced (max {worst:.2e} mm)")
    if len(m) and not ref_diff < 1e-9:
        raise SystemExit(f"stop: the cube's hard-target reference differs from X-142's (max {ref_diff:.2e} rad)")
    v = verdict(att, slopes)
    v.update({"x142_reproduced_to_mm": worst, "reference_matched_pairs": len(m), "reference_max_diff_rad": ref_diff,
              "grassland_incidence_deg": theta})
    OUT.mkdir(parents=True, exist_ok=True)
    p.drop(columns=[c for c in p.columns if c.startswith("ref_")]).to_csv(OUT / "pairs.csv", index=False)
    slopes.to_csv(OUT / "snow_slope.csv", index=False)
    ft.to_csv(OUT / "freeze_thaw.csv", index=False)
    att.to_csv(OUT / "attribution.csv", index=False)
    runs.to_csv(OUT / "stretches.csv", index=False)
    att[["track", "winter", "E_mm", "max_snow_depth_mm", "snow_bound_mm"]].to_csv(OUT / "bound.csv", index=False)
    (OUT / "verdict.json").write_text(json.dumps(v, indent=1, default=float), encoding="utf-8")
    (OUT / "README.md").write_text(
        "# X-146 — dry-snow delay and freeze–thaw in the winter chain (R5 M1, M2, m8)\n\n"
        "Generated by `05_code/SAr-adapted/experiments/winter_snow_x146.py`; pre-registered in "
        f"`preregistration/X-146.yaml`. Stop checks passed: X-142's winter steps reproduced to {worst:.1e} mm; the "
        f"cube's hard-target reference equals X-142's on {len(m)} pairs (max {ref_diff:.1e} rad).\n\n"
        f"## Verdict\n\n- Attribution at P6: **{v['attribution']}** (ΣE {v['sum_E_mm']:.1f} mm = snow-free pairs "
        f"{v['sum_E_free_mm']:.1f} + the rest {v['sum_E_rest_mm']:.1f}).\n- Dry-snow delay on the grassland: "
        f"**{v['snow_delay']}** (slope {v['snow_slope']:.3f} [{v['snow_slope_ci'][0]:.3f}, {v['snow_slope_ci'][1]:.3f}] "
        f"mm per mm of snow depth, {v['dry_snow_pairs']} pairs; model {v['model_range'][0]:.2f} to "
        f"{v['model_range'][1]:.2f}).\n\n"
        "## Attribution at P6 (estimate − laser, mm)\n\n" + md(att) +
        "\n\n## Stretches without a valid laser change\n\n" + md(runs) +
        "\n\n## Snow slopes\n\n" + md(slopes) + "\n\n## Freeze or thaw transitions\n\n" + md(ft) + "\n",
        encoding="utf-8")
    print(json.dumps(v, indent=1, default=float))


if __name__ == "__main__":
    main()

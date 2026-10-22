"""X-129 — the P6 laser's 72-h frost rule against a snow-evidence rule: do the laser changes the second rule recovers
agree with Sentinel-1 as well as the kept ones, and does the by-revisit verdict change?

X-058 drops every laser hour with sub-zero station air in the previous 72 h (rule A), a proxy for snow: the SDMS40
reads snow as surface. Rule B drops an hour only where there is evidence of snow in the previous 72 h (IMGW snow depth
at Piła or Poznań, ERA5-Land snow depth, ERA5 snowfall). The laser at an overpass and its change over a pair are made
exactly as the cube makes them (``gold.site_at``, ``core.value_at``); rule A must reproduce the cube's
``d_laser_los_mm`` or the run stops. Pairs as X-126: P6, 40 m, the labelled cell, hard targets (grassland beside),
consecutive, Sentinel-1D out. The laser test is made on radar-clean pairs (no frost or snow at either overpass).

Pre-registration: hub ``08_deliverables/insar_full_test/preregistration/X-129.yaml`` (checked before the run). Design
years only. Reads the cube; writes only to the hub (``08_deliverables/field_p6/laser_frost_rule_x129/``). No field
value in this file.

    make run E=laser_frost_rule_x129
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import p6_validation_x082 as X  # noqa: E402  (cube root, tracks, references)
import x065_hard_targets_x126 as H  # noqa: E402  (acquisitions between a pair's dates, the bootstrap of r)
from insar_wetlands import field  # noqa: E402
from insar_wetlands.cube import core  # noqa: E402
from insar_wetlands.inversion.isbas import PHASE_TO_MM  # noqa: E402
from insar_wetlands.paths import hub_root  # noqa: E402
from insar_wetlands.scoring import (Protocol, assert_preregistered, guard_lockbox, increment_stats,  # noqa: E402
                                    normalised_sq_dev)

HUB = hub_root()
FT = HUB / "08_deliverables" / "insar_full_test"
OUT = HUB / "08_deliverables" / "field_p6" / "laser_frost_rule_x129"
PREREG = FT / "preregistration" / "X-129.yaml"
PROTOCOL = FT / "protocol.yaml"
X126 = FT / "x065_hard_targets" / "by_revisit.csv"
REGULAR = ("2020_2021", "2022_2024", "2025")
LASER_START = pd.Timestamp("2021-05-27", tz="UTC")
EVIDENCE = {"imgw230_snow_depth_cm": 0.0, "imgw330_snow_depth_cm": 0.0, "om_era5_land_snow_depth": 0.001,
            "om_era5_snowfall": 0.0}
RADAR = ["flag_frozen_air", "flag_snow_om", "flag_snow_imgw"]
VERDICT_GROUPS = {("2021", 6), ("2022–2024", 12), ("2025", 6)}
N_BOOT, SEED, MIN_RECOVERED = 2000, 0, 5


def boolean(s: pd.Series) -> pd.Series:
    return s.map({True: True, False: False, "True": True, "False": False})


def load_site() -> pd.DataFrame:
    cols = ["time_utc", "laser_surface_cm", "laser_ok", "laser_snow_72h", "laser_snow_72h_source", "laser_outlier",
            "laser_filled", "air_c", "cr_raw"] + list(EVIDENCE)
    s = pd.read_csv(X.CUBE / "silver" / "site_hourly.csv.gz", index_col=0, low_memory=False, usecols=cols)
    s.index = pd.to_datetime(s.index, format="ISO8601", utc=True)
    for c in ("laser_ok", "laser_snow_72h", "laser_outlier", "laser_filled"):
        s[c] = boolean(s[c])
    return s


def rules(s: pd.DataFrame) -> pd.DataFrame:
    """Usable laser hours under A (the cube's laser_ok) and B (snow evidence instead of frost)."""
    snow_b = field.snow_evidence_mask(s[list(EVIDENCE)], s.index, hours=72, thresholds=EVIDENCE)
    clean = ~s[["laser_outlier", "laser_filled"]].fillna(True).any(axis=1)
    return pd.DataFrame({"A": s.laser_ok.fillna(False).astype(bool),
                         "B": s.laser_surface_cm.notna() & clean & ~snow_b, "snow_b": snow_b}, index=s.index)


def load_pairs() -> pd.DataFrame:
    cols = ["period", "track", "grid", "plot", "pair", "t1", "t2", "revisit_days", "platforms", "unw_1x1", "los_u",
            "d_laser_los_mm"] + list(X.REFS.values()) + [f"{f}_{t}" for f in RADAR for t in ("t1", "t2")]
    p = pd.read_csv(X.CUBE / "gold" / "plots_pair.csv.gz", usecols=cols, low_memory=False)
    p = p[(p["plot"] == "P6") & (p.grid == "40m") & p.period.astype(str).isin(REGULAR)].copy()
    p = p[~p.platforms.astype(str).str.contains("S1D")]
    for c in ("t1", "t2"):
        p[c] = pd.to_datetime(p[c], format="ISO8601", utc=True)
    for c in [f"{f}_{t}" for f in RADAR for t in ("t1", "t2")]:
        p[c] = boolean(p[c])
    return p.reset_index(drop=True)


def laser_change(s: pd.DataFrame, ok: pd.Series, p: pd.DataFrame) -> np.ndarray:
    """LOS laser change over each pair, as the cube makes it (gold.site_at + pair_drivers + los_u)."""
    lev = s.laser_surface_cm.where(ok)
    a, b = core.value_at(lev, p.t1), core.value_at(lev, p.t2)
    los_u = p.los_u.to_numpy(np.float32).astype(float)       # the cube multiplies by the float32 layer; its CSV prints 8 digits
    return (b - a) * 10 * los_u


def reproduction(p: pd.DataFrame) -> pd.DataFrame:
    g, m = p.d_laser_los_mm.to_numpy(float), p.laser_A.to_numpy(float)
    same_nan = np.isnan(g) == np.isnan(m)
    both = ~np.isnan(g) & ~np.isnan(m)
    return pd.DataFrame([{"pairs": len(p), "with_value_cube": int((~np.isnan(g)).sum()),
                          "with_value_rule_A": int((~np.isnan(m)).sum()), "nan_pattern_equal": bool(same_nan.all()),
                          "max_abs_diff_mm": float(np.max(np.abs(g[both] - m[both]))) if both.any() else 0.0}])


def classify(p: pd.DataFrame) -> pd.DataFrame:
    p = p.copy()
    p["between"] = H.acquisitions_between(p)
    p = p[(p.between == 0) & (p.t1 >= LASER_START)].copy()
    p["year"] = p.t1.dt.year
    p["group"] = np.select([p.year == 2021, p.year.between(2022, 2024), p.year == 2025],
                           ["2021", "2022–2024", "2025"], "other")
    for name, col in X.REFS.items():
        p[f"s1_{name}"] = (p.unw_1x1 - p[col]) * PHASE_TO_MM
    p["radar_clean"] = ~p[[f"{f}_{t}" for f in RADAR for t in ("t1", "t2")]].fillna(True).any(axis=1)
    a, b = p.laser_A.notna(), p.laser_B.notna()
    p["status"] = np.select([a & b, ~a & b, a & ~b], ["kept by both", "recovered by B", "dropped by B only"], "neither")
    p["verdict_group"] = [(g, int(r)) in VERDICT_GROUPS for g, r in zip(p.group, p.revisit_days)]
    return p


def nsd_table(p: pd.DataFrame, ref: str, clean: bool) -> pd.DataFrame:
    """The primary statistic per track: recovered pairs against the kept pairs of their group (laser under B)."""
    rng = np.random.default_rng(SEED)
    rows = []
    q = p[p.verdict_group & p.radar_clean] if clean else p[p.verdict_group & ~p.radar_clean]
    kept_all = p[p.verdict_group & p.radar_clean & (p.status == "kept by both")]
    for track in X.TRACKS:
        rec = q[(q.track == track) & (q.status == "recovered by B")].dropna(subset=["laser_B", f"s1_{ref}"])
        kept = kept_all[kept_all.track == track].dropna(subset=["laser_B", f"s1_{ref}"])
        rm = {g: increment_stats(k.laser_B, k[f"s1_{ref}"])["rmse_increments_mm"] for g, k in kept.groupby("group")}
        dev = (rec[f"s1_{ref}"] - rec.laser_B).to_numpy()
        den = rec.group.map(rm).to_numpy(float)
        est = normalised_sq_dev(dev, den)
        draws = []
        kg = {g: k for g, k in kept.groupby("group")}
        for _ in range(N_BOOT if len(rec) else 0):
            i = rng.integers(0, len(rec), len(rec))
            rmb = {}
            for g, k in kg.items():
                j = rng.integers(0, len(k), len(k))
                rmb[g] = increment_stats(k.laser_B.to_numpy()[j], k[f"s1_{ref}"].to_numpy()[j])["rmse_increments_mm"]
            draws.append(normalised_sq_dev(dev[i], rec.group.map(rmb).to_numpy(float)[i]))
        lo, hi = (np.nanpercentile(draws, [2.5, 97.5]) if draws else (np.nan, np.nan))
        rows.append({"track": track, "reference": ref, "radar": "clean" if clean else "frost or snow at an overpass",
                     "recovered_pairs": len(rec), "kept_pairs": len(kept),
                     "recovered_by_group": ", ".join(f"{g}: {n}" for g, n in rec.group.value_counts().sort_index().items()),
                     "nsd": est, "nsd_lo95": float(lo), "nsd_hi95": float(hi)})
    return pd.DataFrame(rows)


def verdict(n: pd.DataFrame) -> pd.DataFrame:
    h = n[(n.reference == "hard_targets") & (n.radar == "clean")].copy()
    h["track_result"] = np.select(
        [h.recovered_pairs < MIN_RECOVERED, h.nsd_lo95 > 1], ["inconclusive (too few pairs)", "worse"], "not worse")
    overall = ("B REJECTED" if (h.track_result == "worse").any() else
               "B REPLACES A" if (h.track_result == "not worse").all() else "INCONCLUSIVE: A stays")
    h["verdict"] = overall
    return h[["track", "recovered_pairs", "nsd", "nsd_lo95", "nsd_hi95", "track_result", "verdict"]]


def by_revisit(p: pd.DataFrame) -> pd.DataFrame:
    rng = np.random.default_rng(SEED)
    rows = []
    for (track, group, days), g in p.groupby(["track", "group", "revisit_days"]):
        if group == "other":
            continue
        for rule in ("A", "B"):
            m = g[[f"laser_{rule}", "s1_hard_targets"]].dropna()
            if len(m) < 5:
                continue
            a, b = m[f"laser_{rule}"].to_numpy(), m.s1_hard_targets.to_numpy()
            lo, hi = H.boot_r(a, b, rng)
            rows.append({"track": track, "period": group, "revisit_days": int(days), "rule": rule,
                         **increment_stats(a, b), "r_lo95": lo, "r_hi95": hi})
    return pd.DataFrame(rows)


def x126_check(t: pd.DataFrame) -> pd.DataFrame:
    """Rule A's consecutive rows beside X-126's (same pairs, hard targets): a consistency check, not a stop."""
    x = pd.read_csv(X126)
    x = x[(x.reference == "hard_targets") & (x.pairs_kind == "consecutive")]
    a = t[t.rule == "A"].merge(x, on=["track", "period", "revisit_days"], suffixes=("", "_x126"))
    return a[["track", "period", "revisit_days", "pairs", "pairs_x126", "pearson_increments",
              "pearson_increments_x126"]]


def hours_table(s: pd.DataFrame, r: pd.DataFrame) -> pd.DataFrame:
    m = s.laser_surface_cm.notna() & ~s[["laser_outlier", "laser_filled"]].fillna(True).any(axis=1)
    m &= s.index >= LASER_START
    d = pd.DataFrame({"year": s.index.year, "measured": m, "A": r.A & m, "B": r.B & m})[m]
    out = d.groupby("year")[["measured", "A", "B"]].sum()
    out["share_A"], out["share_B"] = out.A / out.measured, out.B / out.measured
    return out.reset_index()


def residual_classes(s: pd.DataFrame, r: pd.DataFrame) -> pd.DataFrame:
    """The exploratory look, repeated on record: laser minus b·CR raw level, minus its 15-day running median of warm
    hours, by class (warm / frost without snow evidence / snow evidence), ±8 days of a filled stretch's ends left out."""
    ok = s.laser_surface_cm.notna() & ~s[["laser_outlier", "laser_filled"]].fillna(True).any(axis=1)
    cold = (s.air_c < 0).astype(float).rolling("72h", min_periods=1).max() > 0
    cls = pd.Series("warm", index=s.index).mask(cold, "frost, no snow evidence").mask(r.snow_b, "snow evidence")
    L, W = s.laser_surface_cm, s.cr_raw
    w = ok & (cls == "warm")
    dL, dW = L[w].resample("D").mean().diff(), W[w].resample("D").mean().diff()
    k = dL.notna() & dW.notna()
    b = np.polyfit(dW[k], dL[k], 1)[0]
    res = L - b * W
    anom = (res - res.where(w).rolling("15D", center=True, min_periods=48).median()).where(ok)
    f = s.laser_filled.fillna(False).astype(bool)
    near = pd.Series(False, index=s.index)
    for e in f.index[f.ne(f.shift())][1:]:
        near |= (s.index > e - pd.Timedelta("8D")) & (s.index < e + pd.Timedelta("8D"))
    a = anom.where(~near & (s.index >= LASER_START))
    g = a.groupby(cls)
    out = g.describe(percentiles=[0.05, 0.5, 0.95])[["count", "5%", "50%", "95%"]]
    out["share_above_0.5cm"] = a.gt(0.5).groupby(cls).sum() / g.count()
    out["share_above_1cm"] = a.gt(1.0).groupby(cls).sum() / g.count()
    out["laser_per_cm_of_cr"] = b
    return out.reset_index(names="class")


def md(df: pd.DataFrame, nd: int = 2) -> str:
    return H.md(df, nd)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    names = ("reproduction.csv", "pairs.csv", "nsd.csv", "verdict.csv", "by_revisit.csv", "x126_check.csv",
             "hours.csv", "residual_classes.csv", "README.md")
    assert_preregistered(PREREG, [OUT / f for f in names])
    s = load_site()
    if (s.get("laser_snow_72h_source", pd.Series(dtype=str)) == "snow_evidence").any():
        sys.exit("X-129 compares rule A with rule B on a rule-A cube; the cube is on rule B since X-130 "
                 "(its result stands: T87; the rule-A cube files are in 06_data/cube_snapshots/laser_rule_A_2026-10-22)")
    r = rules(s)
    p = load_pairs()
    guard_lockbox(pd.concat([p.t1, p.t2]), Protocol.load(PROTOCOL))
    p["laser_A"], p["laser_B"] = laser_change(s, r.A, p), laser_change(s, r.B, p)
    rep = reproduction(p)
    rep.to_csv(OUT / "reproduction.csv", index=False)
    if not (rep.nan_pattern_equal.all() and (rep.max_abs_diff_mm < 1e-9).all()):
        print(rep.to_string())
        sys.exit("stop (pre-registered): rule A does not reproduce the cube's d_laser_los_mm; see reproduction.csv")
    q = classify(p)
    q.drop(columns=["unw_1x1"] + list(X.REFS.values())).to_csv(OUT / "pairs.csv", index=False)
    n = pd.concat([nsd_table(q, ref, clean) for ref in ("hard_targets", "grassland") for clean in (True, False)])
    n.to_csv(OUT / "nsd.csv", index=False)
    v = verdict(n)
    v.to_csv(OUT / "verdict.csv", index=False)
    t = by_revisit(q)
    t.to_csv(OUT / "by_revisit.csv", index=False)
    xc = x126_check(t)
    xc.to_csv(OUT / "x126_check.csv", index=False)
    hrs = hours_table(s, r)
    hrs.to_csv(OUT / "hours.csv", index=False)
    rc = residual_classes(s, r)
    rc.to_csv(OUT / "residual_classes.csv", index=False)
    status = q.groupby(["track", "radar_clean", "status"]).size().rename("pairs").reset_index()
    side = t.pivot_table(index=["track", "period", "revisit_days"], columns="rule",
                         values=["pairs", "pearson_increments", "r_lo95", "r_hi95"]).reset_index()
    side.columns = [" ".join(str(x) for x in c if x != "").strip() for c in side.columns]
    (OUT / "README.md").write_text(f"""# X-129 — the P6 laser's frost rule against a snow-evidence rule (generated)

Generated by SAr-adapted `experiments/laser_frost_rule_x129.py` (`make run E=laser_frost_rule_x129`); do not edit.
Pre-registration: `../../insar_full_test/preregistration/X-129.yaml`. Rule A: no sub-zero station air in the previous
72 h (X-058, the cube's `laser_ok`). Rule B: no snow evidence in the previous 72 h (IMGW snow depth at Piła or Poznań,
ERA5-Land snow depth > 1 mm, ERA5 snowfall). Pairs: P6, 40 m, labelled cell, consecutive, Sentinel-1D out, from
27 May 2021; radar clean = no frost or snow at either overpass. Bootstrap {N_BOOT} resamples, seed {SEED}.

## Verdict (pre-registered rule: hard targets, radar-clean pairs)

{md(v)}

NSD = mean of (S1 − laser)² / RMSE² of the kept pairs of the same group: 1 = the recovered pairs agree with Sentinel-1
as well as the kept ones.

## The statistic, every reading

{md(n)}

## Pairs by status

{md(status)}

## Pair phase against the laser by revisit, rule A and rule B (hard targets, all consecutive pairs)

{md(side)}

## Rule A against X-126 (consistency)

{md(xc)}

## Laser hours kept, per year

{md(hrs)}

## Hourly laser residual after the water level, by class (the exploratory look, on record)

{md(rc, 3)}

## Reproduction of the cube's laser change under rule A (pre-registered stop)

{md(rep, 12)}
""", encoding="utf-8")
    print(v.to_string())


if __name__ == "__main__":
    main()

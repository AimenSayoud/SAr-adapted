"""X-136: the revisit headline read without truth — do the two tracks agree on consecutive six-day pairs, and much less
on pairs that skip an acquisition, in the radar years the analysis never read for motion (2017–2019)?

No laser exists before 2021-05-27 (hub DECISIONS 2026-10-26), so X-131 cannot run. The two tracks see the same mat from
different orbits a day and a half apart: what they share over the same days is the surface signal. Each ascending pair
is matched with the descending pair 1.5 ± 0.6 days later at both ends (same revisit, same kind), both read as X-126
(P6, 40 m, labelled cell of the corrected cube, hard-target reference, Sentinel-1D out), each turned to vertical by its
own line-of-sight up-component; a matched pair is kept when none of its four dates has snow evidence (the laser's own
rule, made from IMGW and ERA5 weather, so the same in every year). The laser years 2021 (from 2021-05-27) and 2025 are
read first, as a gate: there the truth says whether cross-track agreement marks the headline.

Pre-registration: hub ``insar_full_test/preregistration/X-136.yaml``. Writes only to the hub
(``insar_full_test/cross_track_x136/``). No field value in this file.

    make run E=cross_track_x136
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import p6_validation_x082 as X  # noqa: E402
import x065_hard_targets_x126 as X126  # noqa: E402

from insar_wetlands.inversion.isbas import PHASE_TO_MM  # noqa: E402
from insar_wetlands.paths import hub_root  # noqa: E402
from insar_wetlands.scoring import (  # noqa: E402
    Protocol,
    assert_preregistered,
    guard_lockbox,
    increment_stats,
)

HUB = hub_root()
FT = HUB / "08_deliverables" / "insar_full_test"
OUT = FT / "cross_track_x136"
PREREG = FT / "preregistration" / "X-136.yaml"
LASER_START = pd.Timestamp("2021-05-27", tz="UTC")
OFFSET_D, TOL_D = 1.5, 0.6                 # descending dates after ascending, at both ends
N_BOOT, SEED, MIN_PAIRS = 2000, 0, 15
COLS = ["period", "track", "grid", "plot", "t1", "t2", "revisit_days", "platforms", "unw_1x1", "los_u", "coh_3x3",
        "d_laser_los_mm", "d_laser_vertical_mm", "d_wtd_P6_cm", "d_om_era5_land_soil_moisture_0_to_7cm",
        "laser_snow_t1", "laser_snow_t2"] + list(X.REFS.values())


def load_pairs() -> pd.DataFrame:
    """X-126's reader without its laser-span window: every regular set, consecutive or skipping one acquisition."""
    p = pd.read_csv(X.CUBE / "gold" / "plots_pair.csv.gz", usecols=COLS, low_memory=False)
    p = p[(p["plot"] == "P6") & (p.grid == "40m") & p.period.astype(str).isin(X126.REGULAR)].copy()
    p = p[~p.platforms.astype(str).str.contains("S1D")]
    for c in ("t1", "t2"):
        p[c] = pd.to_datetime(p[c], format="ISO8601", utc=True)
    p["between"] = X126.acquisitions_between(p)
    p = p[p.between.isin([0, 1])].copy()
    for name, col in X.REFS.items():
        p[f"s1_{name}"] = (p.unw_1x1 - p[col]) * PHASE_TO_MM
        p[f"v_{name}"] = p[f"s1_{name}"] / p.los_u
    p["snow"] = p[["laser_snow_t1", "laser_snow_t2"]].astype(str).isin(["True", "1", "1.0"]).any(axis=1)
    return p.drop_duplicates(["track", "t1", "t2"]).sort_values(["track", "t1"]).reset_index(drop=True)


def reproduce_x126(p: pd.DataFrame) -> pd.DataFrame:
    """Stop rule (1): per track, the consecutive six-day pairs of 2021 and 2025 inside the laser's span (X-126's
    window) give X-126's promoted rows (hard targets), pairs and r to 1e-9."""
    ref = pd.read_csv(FT / "x065_hard_targets" / "by_revisit.csv")
    ref = ref[(ref.reference == "hard_targets") & (ref.revisit_days == 6) & (ref.pairs_kind == "consecutive")
              & ref.period.astype(str).isin(["2021", "2025"])]
    las = X.laser()[0]
    rows = []
    for track in X.TRACKS:
        lt = las[track].dropna()
        g = p[(p.track == track) & (p.t1 >= lt.index.min() - pd.Timedelta(days=13)) &
              (p.t2 <= lt.index.max() + pd.Timedelta(days=13)) & (p.between == 0) & (p.revisit_days == 6)]
        for y in (2021, 2025):
            m = g[g.t1.dt.year == y][["d_laser_los_mm", "s1_hard_targets"]].dropna()
            mine = increment_stats(m.d_laser_los_mm, m.s1_hard_targets)
            th = ref[(ref.track == track) & (ref.period.astype(str) == str(y))].iloc[0]
            rows.append({"track": track, "year": y, "pairs": mine["pairs"], "pairs_x126": int(th.pairs),
                         "r": mine["pearson_increments"], "r_x126": float(th.pearson_increments)})
    t = pd.DataFrame(rows)
    t["identical"] = (t.pairs == t.pairs_x126) & ((t.r - t.r_x126).abs() < 1e-9)
    return t


def match(p: pd.DataFrame) -> pd.DataFrame:
    """Each ascending pair with the descending pair 1.5 ± 0.6 days later at both ends, same revisit and kind."""
    a = p[p.track == "ascending"].add_prefix("a_")
    d = p[p.track == "descending"].add_prefix("d_")
    m = a.merge(d, left_on=["a_revisit_days", "a_between"], right_on=["d_revisit_days", "d_between"])
    dt1 = (m.d_t1 - m.a_t1).dt.total_seconds() / 86400
    dt2 = (m.d_t2 - m.a_t2).dt.total_seconds() / 86400
    m = m[(dt1 - OFFSET_D).abs().le(TOL_D) & (dt2 - OFFSET_D).abs().le(TOL_D)].copy()
    m = m.sort_values("a_t1").drop_duplicates(["a_t1", "a_t2"]).drop_duplicates(["d_t1", "d_t2"])
    m["kind"] = np.where(m.a_between == 0, "consecutive", "skips one")
    m["revisit_days"] = m.a_revisit_days.astype(int)
    m["snow"] = m.a_snow | m.d_snow
    m["month"] = m.a_t1.dt.month
    t1, t2 = m.a_t1, m.d_t2
    m["period"] = np.select(
        [t1.dt.year.between(2017, 2019), (t1.dt.year >= 2020) & (t2 < LASER_START),
         (t1 >= LASER_START) & (t1.dt.year == 2021), t1.dt.year.between(2022, 2024), t1.dt.year == 2025,
         t1.dt.year == 2026],
        ["2017–2019", "2020 – May 2021", "2021 (laser)", "2022–2024 (laser)", "2025 (laser)", "2026"], "other")
    return m[~m.snow].reset_index(drop=True)


# ------------------------------------------------------------------------------------------- statistics

def r(a, b) -> float:
    a, b = np.asarray(a, float), np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    return float(np.corrcoef(a[ok], b[ok])[0, 1]) if ok.sum() > 2 and np.std(a[ok]) > 0 and np.std(b[ok]) > 0 else np.nan


def boot_r(a, b, rng) -> tuple[float, float]:
    a, b = np.asarray(a, float), np.asarray(b, float)
    draws = [r(a[i], b[i]) for i in (rng.integers(0, len(a), len(a)) for _ in range(N_BOOT))]
    lo, hi = np.nanpercentile(draws, [2.5, 97.5])
    return float(lo), float(hi)


def boot_diff_by_month(g6: pd.DataFrame, g12: pd.DataFrame, ref: str, rng) -> tuple[float, float]:
    """Interval of r(consecutive) − r(skip one), both kinds resampled within each calendar month."""
    draws = []
    m6, m12 = dict(tuple(g6.groupby("month"))), dict(tuple(g12.groupby("month")))
    for _ in range(N_BOOT):
        s6 = pd.concat([x.iloc[rng.integers(0, len(x), len(x))] for x in m6.values()])
        s12 = pd.concat([x.iloc[rng.integers(0, len(x), len(x))] for x in m12.values()])
        draws.append(r(s6[f"a_v_{ref}"], s6[f"d_v_{ref}"]) - r(s12[f"a_v_{ref}"], s12[f"d_v_{ref}"]))
    lo, hi = np.nanpercentile(draws, [2.5, 97.5])
    return float(lo), float(hi)


def residual_r(g: pd.DataFrame, ref: str, cov_a: str, cov_d: str) -> float:
    """Cross-track r after each track's increment is regressed (with intercept) on its own covariate."""
    ok = g[[f"a_v_{ref}", f"d_v_{ref}", cov_a, cov_d]].dropna()
    if len(ok) < 5:
        return np.nan
    res = []
    for v, c in ((f"a_v_{ref}", cov_a), (f"d_v_{ref}", cov_d)):
        A = np.column_stack([np.ones(len(ok)), ok[c]])
        beta, *_ = np.linalg.lstsq(A, ok[v], rcond=None)
        res.append(ok[v] - A @ beta)
    return r(*res)


def signature(m: pd.DataFrame, period: str, ref: str, short: int, long_: int) -> dict:
    """(a), (b), (c) for one period: consecutive pairs of ``short`` days against pairs of ``long_`` days skipping one."""
    rng = np.random.default_rng(SEED)
    g = m[m.period == period]
    g6 = g[(g.kind == "consecutive") & (g.revisit_days == short)]
    g12 = g[(g.kind == "skips one") & (g.revisit_days == long_)]
    out = {"period": period, "reference": ref, "consecutive_days": short, "skip_days": long_,
           "n_consecutive": len(g6), "n_skip": len(g12)}
    if len(g6) < 3 or len(g12) < 3:
        return out
    out["r_consecutive"] = r(g6[f"a_v_{ref}"], g6[f"d_v_{ref}"])
    out["r_consecutive_lo95"], out["r_consecutive_hi95"] = boot_r(g6[f"a_v_{ref}"], g6[f"d_v_{ref}"], rng)
    out["r_skip"] = r(g12[f"a_v_{ref}"], g12[f"d_v_{ref}"])
    out["r_skip_lo95"], out["r_skip_hi95"] = boot_r(g12[f"a_v_{ref}"], g12[f"d_v_{ref}"], rng)
    out["difference"] = out["r_consecutive"] - out["r_skip"]
    out["difference_lo95"], out["difference_hi95"] = boot_diff_by_month(g6, g12, ref, rng)
    for kind, gg in (("consecutive", g6), ("skip", g12)):
        out[f"r_{kind}_moisture_removed"] = residual_r(gg, ref, "a_d_om_era5_land_soil_moisture_0_to_7cm",
                                                       "d_d_om_era5_land_soil_moisture_0_to_7cm")
        out[f"median_coh_{kind}"] = float(np.nanmedian(np.r_[gg.a_coh_3x3, gg.d_coh_3x3]))
    if "laser" in period:
        for kind, gg in (("consecutive", g6), ("skip", g12)):
            out[f"r_asc_laser_{kind}"] = r(gg[f"a_v_{ref}"], gg.a_d_laser_vertical_mm)
            out[f"r_desc_laser_{kind}"] = r(gg[f"d_v_{ref}"], gg.d_d_laser_vertical_mm)
            out[f"r_{kind}_laser_removed"] = residual_r(gg, ref, "a_d_laser_vertical_mm", "d_d_laser_vertical_mm")
    return out


def reading(row: dict | pd.Series) -> str:
    if row.get("n_consecutive", 0) < MIN_PAIRS:
        return "too few to judge"
    if not (row["r_consecutive_lo95"] > 0):
        return "not replicated"
    return "replicates" if row["difference_lo95"] > 0 else "replicates in part"


def water_table(m: pd.DataFrame) -> pd.DataFrame:
    """Secondary: r of the two tracks' mean vertical increment (consecutive pairs) with the P6 logger's change."""
    rows = []
    g = m[m.kind == "consecutive"].copy()
    g["mean_v"] = (g.a_v_hard_targets + g.d_v_hard_targets) / 2
    for label, sel in (("2020 – May 2021 (no laser)", g.period == "2020 – May 2021"),
                       ("2021-05-27 → 2024 (laser years with a logger)", g.period.isin(["2021 (laser)", "2022–2024 (laser)"]))):
        gg = g[sel][["mean_v", "a_d_wtd_P6_cm", "a_d_laser_vertical_mm"]]
        rows.append({"span": label, "pairs": int(gg[["mean_v", "a_d_wtd_P6_cm"]].dropna().shape[0]),
                     "r_mean_track_vs_water_table": r(gg.mean_v, gg.a_d_wtd_P6_cm),
                     "r_laser_vs_water_table": r(gg.a_d_laser_vertical_mm, gg.a_d_wtd_P6_cm)})
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    assert_preregistered(PREREG, [OUT / "verdict.csv"])
    p = load_pairs()
    guard_lockbox(pd.concat([p.t1, p.t2]), Protocol.load(FT / "protocol.yaml"))
    rep = reproduce_x126(p)
    rep.to_csv(OUT / "reproduction_x126.csv", index=False)
    if not rep.identical.all():
        print(rep.to_string())
        sys.exit("STOP (X-136 stop rule 1): the reader does not reproduce X-126; see reproduction_x126.csv")
    m = match(p)
    m["calib"] = m.period.isin(["2021 (laser)", "2025 (laser)"])
    cal = m.assign(period=np.where(m.calib, "2021 + 2025 (laser)", m.period))
    n_cal = int(((cal.period == "2021 + 2025 (laser)") & (cal.kind == "consecutive") & (cal.revisit_days == 6)).sum())
    if n_cal < MIN_PAIRS:
        sys.exit(f"STOP (X-136 stop rule 2): {n_cal} consecutive matched pairs in the calibration years")
    # the gate first, then everything else
    rows = []
    gate = signature(cal, "2021 + 2025 (laser)", "hard_targets", 6, 12)
    gate["role"] = "calibration gate"
    rows.append(gate)
    gate_ok = gate.get("r_consecutive_lo95", -1) > 0 and gate.get("difference_lo95", -1) > 0
    for period, role in (("2017–2019", "unseen (primary)"), ("2020 – May 2021", "apart: read against the water table before"),
                         ("2026", "apart: Sentinel-1A + 1C")):
        rows.append({**signature(cal, period, "hard_targets", 6, 12), "role": role})
    for y in (2017, 2018, 2019):
        sub = cal[(cal.period == "2017–2019") & (cal.a_t1.dt.year == y)].assign(period=str(y))
        rows.append({**signature(sub, str(y), "hard_targets", 6, 12), "role": "unseen, per year"})
    for period, role in (("2021 + 2025 (laser)", "calibration, grassland"), ("2017–2019", "unseen, grassland")):
        rows.append({**signature(cal, period, "grassland", 6, 12), "role": role})
    rows.append({**signature(cal, "2022–2024 (laser)", "hard_targets", 12, 24), "role": "second calibration (12 vs 24 days)"})
    S = pd.DataFrame(rows)
    S["reading"] = [reading(x) for x in S.to_dict("records")]
    prim = S[S.role == "unseen (primary)"].iloc[0]
    verdict = pd.DataFrame([{"gate_passes": gate_ok, "gate_reading": reading(gate),
                             "unseen_reading": reading(prim) if gate_ok else "uninformative (gate failed)",
                             "n_consecutive_unseen": int(prim.n_consecutive), "n_skip_unseen": int(prim.n_skip)}])
    W = water_table(m)
    keep = [c for c in m.columns if not c.endswith(("unw_1x1", "ref_hard_unw", "ref_grassland_unw"))]
    m[keep].to_csv(OUT / "matched_pairs.csv.gz", index=False)
    for name, df in (("signature", S), ("water_table", W), ("verdict", verdict)):
        df.to_csv(OUT / f"{name}.csv", index=False)
    (OUT / "README.md").write_text(readme(S, verdict, W, rep, m), encoding="utf-8")
    print(verdict.T.to_string())
    print(S[["role", "period", "n_consecutive", "n_skip", "r_consecutive", "r_consecutive_lo95", "r_skip", "difference",
             "difference_lo95", "reading"]].round(2).to_string())


def fmt(x, nd=2):
    return "" if x is None or (isinstance(x, float) and np.isnan(x)) else (f"{x:.{nd}f}" if isinstance(x, float) else str(x))


def readme(S, v, W, rep, m) -> str:
    v0 = v.iloc[0]
    cols = ["role", "period", "reference", "consecutive_days", "skip_days", "n_consecutive", "n_skip", "r_consecutive",
            "r_consecutive_lo95", "r_consecutive_hi95", "r_skip", "r_skip_lo95", "r_skip_hi95", "difference",
            "difference_lo95", "difference_hi95", "reading"]
    extra = ["role", "period", "reference", "r_consecutive_moisture_removed", "r_skip_moisture_removed",
             "median_coh_consecutive", "median_coh_skip"]
    las = ["role", "period", "r_asc_laser_consecutive", "r_desc_laser_consecutive", "r_consecutive_laser_removed",
           "r_asc_laser_skip", "r_desc_laser_skip", "r_skip_laser_removed"]

    def tab(df, cs):
        cs = [c for c in cs if c in df.columns]
        return "\n".join(["| " + " | ".join(cs) + " |", "|" + "---|" * len(cs)] +
                         ["| " + " | ".join(fmt(x) for x in row) + " |" for row in df[cs].itertuples(index=False)])
    L = ["# X-136 — the revisit headline across the two tracks, without truth (generated)", "",
         "Generated by `experiments/cross_track_x136.py` (`make run E=cross_track_x136`); do not edit. Pre-registration:",
         "`../preregistration/X-136.yaml`. P6, 40 m, labelled cell of the corrected cube, hard-target reference,",
         "Sentinel-1D out; each ascending pair matched with the descending pair 1.5 ± 0.6 days later at both ends; every",
         "increment in vertical; matched pairs with snow evidence on any of their four dates dropped (the laser's rule,",
         f"from IMGW and ERA5 weather). Bootstrap {N_BOOT} resamples, seed {SEED}; the difference resampled within months.", "",
         "## Verdict", "",
         f"- Calibration gate (2021 from 2021-05-27, and 2025, where the laser confirms the headline): "
         f"**{'passes' if v0.gate_passes else 'fails'}** ({v0.gate_reading}).",
         f"- Unseen radar years 2017–2019: **{v0.unseen_reading}** ({v0.n_consecutive_unseen} consecutive six-day and"
         f" {v0.n_skip_unseen} skip-one matched pairs).", "",
         "## The signature, period by period", "", tab(S, cols), "",
         "## What the agreement is made of", "",
         "Moisture control (each track's increment regressed on the ERA5-Land soil-moisture change over its window) and"
         " coherence:", "", tab(S, extra), "",
         "Laser years only — each track against the laser, and the cross-track r once each track's own laser change is"
         " regressed out (near zero: the agreement is the mat's motion):", "",
         tab(S[S.period.str.contains("laser")], las), "",
         "## The water table where the laser is missing (secondary)", "", tab(W, list(W.columns)), "",
         "## Stop rule", "", tab(rep, list(rep.columns)), "",
         f"Matched pairs kept (snow-free): {len(m)}; by period and kind in `matched_pairs.csv.gz`.", ""]
    return "\n".join(L)


if __name__ == "__main__":
    main()

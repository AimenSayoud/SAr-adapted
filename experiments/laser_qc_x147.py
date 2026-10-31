"""X-147 — laser QC amendment (R6 M3): Table 2 without the pairs whose laser change crosses a filled stretch, and the
article's across-year laser numbers with the April 2024 step that is not the mat removed.

(A) Table 2's twelve rows, grouped as X-126 groups them, with and without the pairs whose laser change crosses a filled
stretch of the delivered record (laser QC ``interpolated.csv``); r, and the wider of a pair bootstrap and a
calendar-month bootstrap, as Table 2 quotes it. (B) The April 2024 step: the change of the residual laser − g·water
(g the article's buoyancy coefficient) between the week before 2024-04-08 and the week after 2024-04-16, from the measured
hours of the cube's hourly laser; X-140's 2023/24 laser winter, 2023 loop and 2023 → 2025 check with that step removed.
Pre-registered in the hub's ``insar_full_test/preregistration/X-147.yaml``. Reads hub tables and the cube only; writes
only to the hub ``insar_full_test/laser_qc_x147/``. No field value in this file.

    make run E=laser_qc_x147
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from insar_wetlands.paths import cube_root, hub_root
from insar_wetlands.scoring import assert_preregistered

HUB, CUBE = hub_root(), cube_root()
FT = HUB / "08_deliverables" / "insar_full_test"
OUT = FT / "laser_qc_x147"
PREREG = FT / "preregistration" / "X-147.yaml"
QC = HUB / "08_deliverables" / "field_p6" / "laser_qc"
X140 = HUB / "08_deliverables" / "drying_x140"
TABLE2 = HUB / "08_deliverables" / "article_plan" / "tables" / "revisit_correlations.csv"
G = HUB / "08_deliverables" / "fusion_v3" / "wp1" / "e12_coefficients.csv"
ROWS = (("2021", 6, 0), ("2021", 12, 1), ("2022–2024", 12, 0), ("2022–2024", 24, 1), ("2025", 6, 0), ("2025", 12, 1))
N_BOOT, SEED = 2000, 147
BEFORE, AFTER = ("2024-04-01 00:00", "2024-04-08 08:00"), ("2024-04-16 00:00", "2024-04-23 00:00")
STEP_1, STEP_2 = ("2024-04-08 08:00", "2024-04-08 11:00"), ("2024-04-15 11:00", "2024-04-15 23:00")


def r(a, b) -> float:
    return float(np.corrcoef(a, b)[0, 1]) if len(a) > 2 and np.std(a) > 0 and np.std(b) > 0 else np.nan


def intervals(g: pd.DataFrame, rng: np.random.Generator) -> dict:
    """Pair bootstrap and calendar-month bootstrap of r (every pair whose first date falls in a drawn month)."""
    a, b = g.d_laser_los_mm.to_numpy(), g.s1_hard_targets.to_numpy()
    pi = [r(a[i], b[i]) for i in (rng.integers(0, len(a), len(a)) for _ in range(N_BOOT))]
    pos = {m: np.flatnonzero(g.month.to_numpy() == m) for m in np.unique(g.month)}
    keys = list(pos)
    mi = []
    for _ in range(N_BOOT):
        i = np.concatenate([pos[keys[k]] for k in rng.integers(0, len(keys), len(keys))])
        mi.append(r(a[i], b[i]))
    p_lo, p_hi = np.nanpercentile(pi, [2.5, 97.5])
    m_lo, m_hi = np.nanpercentile(mi, [2.5, 97.5])
    wider = (m_hi - m_lo) > (p_hi - p_lo)
    return {"lo95": float(m_lo if wider else p_lo), "hi95": float(m_hi if wider else p_hi),
            "resampled": "months" if wider else "pairs"}


def table2(p: pd.DataFrame, fills: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    hit = np.zeros(len(p), bool)
    for _, f in fills.iterrows():
        hit |= ((p.t1 < f.end) & (p.t2 > f.start)).to_numpy()
    p = p.assign(crosses_fill=hit)
    rng = np.random.default_rng(SEED)
    rows = []
    for track in ("ascending", "descending"):
        for period, days, between in ROWS:
            g = p[(p.track == track) & (p.group == period) & (p.revisit_days == days) & (p.between == between)]
            k = g[~g.crosses_fill]
            iw, ik = intervals(g, rng), intervals(k, rng)
            rows.append({"track": track, "period": period, "revisit_days": days, "between": between,
                         "n_all": len(g), "r_all": r(g.d_laser_los_mm, g.s1_hard_targets),
                         "lo95_all": iw["lo95"], "hi95_all": iw["hi95"], "resampled_all": iw["resampled"],
                         "n_kept": len(k), "r_kept": r(k.d_laser_los_mm, k.s1_hard_targets),
                         "lo95_kept": ik["lo95"], "hi95_kept": ik["hi95"], "resampled_kept": ik["resampled"]})
    t = pd.DataFrame(rows)
    t["delta_r"] = t.r_kept - t.r_all
    return t, p[p.crosses_fill]


def residual(h: pd.DataFrame, g: float, a: str, b: str) -> pd.Series:
    w = h.loc[a:b]
    return (w.laser_surface_raw_cm - g * w.cr_raw).dropna()


def main() -> None:
    assert_preregistered(PREREG, [OUT / "verdict.json"])
    OUT.mkdir(parents=True, exist_ok=True)
    fills = pd.read_csv(QC / "interpolated.csv", parse_dates=["start", "end"])
    p = pd.read_csv(FT / "x065_hard_targets" / "pairs.csv", dtype={"group": str})
    p = p.dropna(subset=["d_laser_los_mm", "s1_hard_targets"])
    for c in ("t1", "t2"):
        p[c] = pd.to_datetime(p[c], utc=True)
    p["month"] = p.t1.dt.strftime("%Y-%m")
    t2, crossing = table2(p, fills)
    pub = pd.read_csv(TABLE2, dtype={"period": str})
    chk = t2.merge(pub, on=["track", "period", "revisit_days"])
    worst_r = float((chk.r_all - chk.r).abs().max())
    if len(chk) != len(t2) or worst_r > 1e-9:
        raise SystemExit(f"X-147 stop: Table 2's r not reproduced (worst {worst_r:.2e})")
    t2.to_csv(OUT / "table2_fill.csv", index=False)
    t1 = pd.Timestamp(STEP_1[0], tz="UTC")
    over8 = p[(p.t1 < t1) & (p.t2 > pd.Timestamp(STEP_1[1], tz="UTC"))]
    pd.concat([crossing.assign(why="crosses a filled stretch"),
               over8.assign(why="spans the 2024-04-08 step")])[
        ["track", "pair", "revisit_days", "between", "group", "d_laser_los_mm", "s1_hard_targets", "why"]
    ].to_csv(OUT / "pairs_crossing.csv", index=False)

    # (B) the April 2024 step from the measured hours
    g = float(pd.read_csv(G).iloc[0]["g·ΔW"]) / 10                       # mm per cm → fraction
    h = pd.read_csv(CUBE / "gold" / "truth_laser_hourly.csv.gz",
                    usecols=["time_utc", "laser_surface_raw_cm", "laser_filled", "cr_raw", "laser_snow_evidence_72h"])
    h["time_utc"] = pd.to_datetime(h.time_utc, utc=True)
    h = h[~h.laser_filled.astype(str).str.lower().isin(["true", "1", "1.0"])].set_index("time_utc").sort_index()
    steps = []
    for days in (3, 7, 14):
        b0 = pd.Timestamp(BEFORE[1], tz="UTC")
        a0 = pd.Timestamp(AFTER[0], tz="UTC")
        before = residual(h, g, str(b0 - pd.Timedelta(days=days)), str(b0)).median()
        after = residual(h, g, str(a0), str(a0 + pd.Timedelta(days=days))).median()
        steps.append({"part": "8–16 April 2024", "window_days": days, "step_cm": float(after - before)})
    for name, (a, b) in (("8 April 2024", STEP_1), ("15 April 2024", STEP_2)):
        before = residual(h, g, str(pd.Timestamp(a, tz="UTC") - pd.Timedelta(hours=12)), a).median()
        after = residual(h, g, b, str(pd.Timestamp(b, tz="UTC") + pd.Timedelta(hours=12))).median()
        steps.append({"part": name, "window_days": 0.5, "step_cm": float(after - before)})
    st = pd.DataFrame(steps)
    st.to_csv(OUT / "step.csv", index=False)
    step = float(st[(st.part == "8–16 April 2024") & (st.window_days == 7)].step_cm.iloc[0])
    rng_lo = float(st[st.part == "8–16 April 2024"].step_cm.min())
    rng_hi = float(st[st.part == "8–16 April 2024"].step_cm.max())

    loops = pd.read_csv(X140 / "loops.csv").set_index("y")
    v140 = json.loads((X140 / "verdict.json").read_text())
    lmean = loops[["asc_laser_loop_cm", "des_laser_loop_cm"]].mean(axis=1)
    if abs(float(lmean[2023] + lmean[2024]) - float(v140["check_2023_2025_laser_cm"])) > 1e-6:
        raise SystemExit("X-147 stop: X-140's 2023 → 2025 laser check is not the sum of its mean loops")
    w = pd.read_csv(X140 / "winters.csv")
    w = w[(w.reference == "hard_targets") & w.laser_cm.notna()]
    lz = w.groupby("year").laser_cm.mean()
    lz_c = lz.copy()
    lz_c[2023] = lz[2023] - step
    radar = float(v140["check_2023_2025_radar_cm"])
    laser_c = float(v140["check_2023_2025_laser_cm"]) - step
    loop_err = float(loops.loop_err_cm.dropna().iloc[0])
    ay = pd.DataFrame([{
        "step_cm": step, "step_lo_cm": rng_lo, "step_hi_cm": rng_hi,
        "laser_winter_2023_cm": float(lz[2023]), "laser_winter_2023_corrected_cm": float(lz_c[2023]),
        "laser_loop_2023_cm": float(lmean[2023]), "laser_loop_2023_corrected_cm": float(lmean[2023] - step),
        "laser_loop_2024_cm": float(lmean[2024]),
        "check_radar_cm": radar, "check_laser_cm": float(v140["check_2023_2025_laser_cm"]),
        "check_laser_corrected_cm": laser_c, "check_difference_corrected_cm": radar - laser_c,
        "loop_err_cm": loop_err,
        "winter_2023_rank_corrected": int(lz_c.rank(ascending=False)[2023]), "laser_winters": int(len(lz_c)),
    }])
    ay.to_csv(OUT / "across_years.csv", index=False)

    # (secondary) December 2021, descriptive: snow-free days either side of the fill
    sf = h[h.laser_snow_evidence_72h == 0]
    d = (sf.laser_surface_raw_cm - g * sf.cr_raw).resample("1D").median().dropna()
    pre = d.loc["2021-11-17":"2021-11-30"]
    post = d.loc["2021-12-24":].iloc[:14]
    fill = fills[fills.start.dt.strftime("%Y-%m") == "2021-12"].iloc[0]
    pd.DataFrame([{"before_days": len(pre), "before_median_cm": float(pre.median()),
                   "after_days": len(post), "after_first": str(post.index.min().date()),
                   "after_last": str(post.index.max().date()), "after_median_cm": float(post.median()),
                   "difference_cm": float(post.median() - pre.median()),
                   "laser_change_across_fill_cm": float(fill.laser_change_cm)}]).to_csv(OUT / "december_2021.csv",
                                                                                         index=False)

    a_ok = bool((t2.delta_r.abs() < 0.10).all())
    for period in ("2021", "2025"):
        for track in ("ascending", "descending"):
            q = t2[(t2.track == track) & (t2.period == period)]
            a_ok &= bool(q[q.between == 0].r_kept.iloc[0] > q[q.between == 1].r_kept.iloc[0])
    b_ok = (laser_c - radar) > loop_err
    verdict = {
        "table2": "Table 2 unchanged in its reading" if a_ok else "Table 2 to be rebuilt without the pairs",
        "pairs_crossing_fill": int(len(crossing)), "max_abs_delta_r": float(t2.delta_r.abs().max()),
        "step_cm": step, "step_range_cm": [rng_lo, rng_hi],
        "losses_do_not_cancel": "kept" if b_ok else "withdrawn",
        "check_difference_corrected_cm": radar - laser_c, "loop_err_cm": loop_err,
        "winter_2023_rank_corrected": int(lz_c.rank(ascending=False)[2023]),
        "stop_checks": {"table2_r_reproduced_to": worst_r},
    }
    (OUT / "verdict.json").write_text(json.dumps(verdict, indent=1, ensure_ascii=False))
    md = lambda df: df.to_markdown(index=False, floatfmt=".3f")  # noqa: E731
    (OUT / "README.md").write_text(
        "# X-147 — laser QC amendment: pairs across filled stretches, the April 2024 step (R6 M3)\n\n"
        "Generated by `05_code/SAr-adapted/experiments/laser_qc_x147.py`; pre-registered in "
        f"`preregistration/X-147.yaml`. Stop checks passed: Table 2's r reproduced to {worst_r:.1e}; X-140's 2023 → 2025 "
        "laser check equals its mean loops.\n\n"
        f"## Verdicts\n\n```json\n{json.dumps(verdict, indent=1, ensure_ascii=False)}\n```\n\n"
        "## (A) Table 2 with and without the pairs whose laser change crosses a filled stretch\n\n" + md(t2) +
        "\n\n## Pairs crossing a filled stretch or the 8 April 2024 step\n\n" + md(pd.read_csv(OUT / "pairs_crossing.csv")) +
        "\n\n## (B) The April 2024 step (laser − g·water, measured hours)\n\n" + md(st) +
        "\n\n## Across-year laser numbers with the step removed\n\n" + md(ay) +
        "\n\n## December 2021 (descriptive)\n\n" + md(pd.read_csv(OUT / "december_2021.csv")) + "\n",
        encoding="utf-8")
    print(json.dumps(verdict, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()

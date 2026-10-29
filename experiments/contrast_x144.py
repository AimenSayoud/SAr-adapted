"""X-144 — consecutive minus skip-one as a paired test, Table 2 with month-block intervals, the cross-track slope (R4
M4, M6, m6).

Pairs that share an acquisition are not independent, and one pair's error enters two adjacent increments. So the
bootstrap here resamples calendar months (every pair whose first date falls in a drawn month), and a contrast draws the
same months for both of its groups. Inputs are X-126's pairs (Table 2) and X-136's matched pairs (cross-track).

Pre-registration: hub ``08_deliverables/insar_full_test/preregistration/X-144.yaml`` (checked before the run).
Reads the hub; writes only to the hub: ``08_deliverables/insar_full_test/contrast_x144/``. No field value in this file.

    make run E=contrast_x144
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from insar_wetlands.paths import hub_root
from insar_wetlands.scoring import assert_preregistered, increment_stats

HUB = hub_root()
FT = HUB / "08_deliverables" / "insar_full_test"
OUT = FT / "contrast_x144"
PREREG = FT / "preregistration" / "X-144.yaml"
PAIRS = FT / "x065_hard_targets" / "pairs.csv"
BY_REVISIT = FT / "x065_hard_targets" / "by_revisit.csv"
MATCHED = FT / "cross_track_x136" / "matched_pairs.csv.gz"
N_BOOT, SEED = 2000, 144
TRACKS = ("ascending", "descending")
CONTRASTS = (("2021", 6, 12), ("2025", 6, 12), ("2022–2024", 12, 24))   # period, consecutive days, skip-one days


def r(a, b) -> float:
    a, b = np.asarray(a, float), np.asarray(b, float)
    return float(np.corrcoef(a, b)[0, 1]) if len(a) > 2 and np.std(a) > 0 and np.std(b) > 0 else np.nan


def slope0(truth, est) -> float:
    return increment_stats(truth, est)["slope_through_origin"]


def rma(x, y) -> float:
    """Reduced-major-axis slope of y on x (both noisy)."""
    rr = r(x, y)
    return float(np.sign(rr) * np.std(y) / np.std(x)) if np.isfinite(rr) else np.nan


def month_draws(months: np.ndarray, rng: np.random.Generator, span: int = 1):
    """Yield, per draw, the drawn block labels (with replacement); ``span`` months per block for the sensitivity."""
    u = np.unique(months)
    blocks = [u[i:i + span] for i in range(0, len(u), span)]
    for _ in range(N_BOOT):
        yield [blocks[j] for j in rng.integers(0, len(blocks), len(blocks))]


class Months:
    """Row positions of a table by month, built once, so each draw only concatenates index arrays."""

    def __init__(self, g: pd.DataFrame, cols: tuple[str, str]):
        self.x, self.y = g[cols[0]].to_numpy(float), g[cols[1]].to_numpy(float)
        self.pos = {m: np.flatnonzero(g.month.to_numpy() == m) for m in np.unique(g.month)}

    def take(self, drawn) -> tuple[np.ndarray, np.ndarray]:
        parts = [self.pos[m] for blk in drawn for m in blk if m in self.pos]
        i = np.concatenate(parts) if parts else np.array([], int)
        return self.x[i], self.y[i]


def load() -> pd.DataFrame:
    p = pd.read_csv(PAIRS)
    p["t1"] = pd.to_datetime(p.t1, format="ISO8601", utc=True)
    p["month"] = p.t1.dt.strftime("%Y-%m")
    p["group"] = p.group.astype(str)
    return p


def contrasts(p: pd.DataFrame, ref: str, span: int, rng: np.random.Generator) -> pd.DataFrame:
    rows = []
    for period, dc, ds in CONTRASTS:
        for track in TRACKS:
            g = p[(p.group == period) & (p.track == track)].dropna(subset=["d_laser_los_mm", f"s1_{ref}"])
            gc = g[(g.between == 0) & (g.revisit_days == dc)]
            gs = g[(g.between == 1) & (g.revisit_days == ds)]
            d0 = r(gc.d_laser_los_mm, gc[f"s1_{ref}"]) - r(gs.d_laser_los_mm, gs[f"s1_{ref}"])
            draws = []
            mc, ms = Months(gc, ("d_laser_los_mm", f"s1_{ref}")), Months(gs, ("d_laser_los_mm", f"s1_{ref}"))
            for drawn in month_draws(g.month.to_numpy(), rng, span):
                draws.append(r(*mc.take(drawn)) - r(*ms.take(drawn)))
            lo, hi = np.nanpercentile(draws, [2.5, 97.5])
            rows.append({"period": period, "track": track, "reference": ref, "block_months": span,
                         "consecutive_days": dc, "skip_days": ds, "n_consecutive": len(gc), "n_skip": len(gs),
                         "months": g.month.nunique(), "delta_r": d0, "lo95": float(lo), "hi95": float(hi),
                         "reading": "consecutive pairs follow the laser better than skip-one pairs" if lo > 0
                         else "the difference cannot be told from zero"})
    return pd.DataFrame(rows)


def table2_blocks(p: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    rows = []
    for (track, group, between, days), g in p.groupby(["track", "group", "between", "revisit_days"]):
        g = g.dropna(subset=["d_laser_los_mm", "s1_hard_targets"])
        if group == "other" or len(g) < 5:
            continue
        rs, ss = [], []
        mg = Months(g, ("d_laser_los_mm", "s1_hard_targets"))
        for drawn in month_draws(g.month.to_numpy(), rng):
            a, b = mg.take(drawn)
            rs.append(r(a, b))
            ss.append(slope0(a, b))
        rows.append({"track": track, "period": group, "revisit_days": int(days),
                     "pairs_kind": "consecutive" if between == 0 else "skips one", "n": len(g),
                     "months": g.month.nunique(), "r": r(g.d_laser_los_mm, g.s1_hard_targets),
                     "r_lo95_month": float(np.nanpercentile(rs, 2.5)), "r_hi95_month": float(np.nanpercentile(rs, 97.5)),
                     "slope0": slope0(g.d_laser_los_mm, g.s1_hard_targets),
                     "slope0_lo95_month": float(np.nanpercentile(ss, 2.5)),
                     "slope0_hi95_month": float(np.nanpercentile(ss, 97.5))})
    t = pd.DataFrame(rows)
    b = pd.read_csv(BY_REVISIT)
    b = b[b.reference == "hard_targets"][["track", "period", "revisit_days", "pearson_increments", "r_lo95", "r_hi95"]]
    b["period"] = b.period.astype(str)
    return t.merge(b.rename(columns={"r_lo95": "r_lo95_iid", "r_hi95": "r_hi95_iid", "pearson_increments": "r_x126"}),
                   on=["track", "period", "revisit_days"], how="left")


def cross_track(rng: np.random.Generator) -> pd.DataFrame:
    m = pd.read_csv(MATCHED)
    m["a_t1"] = pd.to_datetime(m.a_t1, format="ISO8601", utc=True)
    m = m[(m.kind == "consecutive") & (m.revisit_days == 6) & m.a_t1.dt.year.isin([2021, 2025]) & ~m.snow.astype(bool)]
    m = m.dropna(subset=["a_v_hard_targets", "d_v_hard_targets"])
    m["month"] = m.a_t1.dt.strftime("%Y-%m")
    x, y = m.a_v_hard_targets.to_numpy(), m.d_v_hard_targets.to_numpy()
    stats = {"rma_desc_on_asc": rma, "ols_desc_on_asc": lambda a, b: float(np.polyfit(a, b, 1)[0]),
             "ols_asc_on_desc": lambda a, b: float(np.polyfit(b, a, 1)[0]), "r": r}
    draws = {k: [] for k in stats}
    mm = Months(m, ("a_v_hard_targets", "d_v_hard_targets"))
    for drawn in month_draws(m.month.to_numpy(), rng):
        a, b = mm.take(drawn)
        for k, f in stats.items():
            draws[k].append(f(a, b))
    return pd.DataFrame([{"statistic": k, "value": f(x, y), "lo95": float(np.nanpercentile(draws[k], 2.5)),
                          "hi95": float(np.nanpercentile(draws[k], 97.5)), "pairs": len(m), "months": m.month.nunique()}
                         for k, f in stats.items()])


def md(df: pd.DataFrame) -> str:
    d = df.copy()
    for c in d.columns:
        if pd.api.types.is_float_dtype(d[c]):
            d[c] = d[c].round(2)
    head = "| " + " | ".join(map(str, d.columns)) + " |\n|" + "---|" * len(d.columns) + "\n"
    return head + "\n".join("| " + " | ".join(map(str, row)) + " |" for row in d.itertuples(index=False))


def main() -> None:
    assert_preregistered(PREREG, [OUT / "verdict.json"])
    rng = np.random.default_rng(SEED)
    p = load()
    t2 = table2_blocks(p, rng)
    diff = float((t2.r - t2.r_x126).abs().max())
    if not diff < 1e-9:
        raise SystemExit(f"stop: Table 2's r not reproduced (max difference {diff:.2e})")
    c = pd.concat([contrasts(p, "hard_targets", 1, rng), contrasts(p, "hard_targets", 2, rng),
                   contrasts(p, "grassland", 1, rng)], ignore_index=True)
    xt = cross_track(rng)
    OUT.mkdir(parents=True, exist_ok=True)
    c.to_csv(OUT / "contrasts.csv", index=False)
    t2.to_csv(OUT / "table2_blocks.csv", index=False)
    xt.to_csv(OUT / "cross_track_slope.csv", index=False)
    prim = c[(c.reference == "hard_targets") & (c.block_months == 1)]
    v = {"primary": prim[["period", "track", "delta_r", "lo95", "hi95", "reading"]].to_dict("records"),
         "table2_reproduced_to": diff}
    (OUT / "verdict.json").write_text(json.dumps(v, indent=1), encoding="utf-8")
    (OUT / "README.md").write_text(
        "# X-144 — consecutive minus skip-one, paired; Table 2 with month blocks; the cross-track slope (R4 M4, M6, m6)\n\n"
        "Generated by `05_code/SAr-adapted/experiments/contrast_x144.py`; pre-registered in `preregistration/X-144.yaml`. "
        f"Month-block bootstrap ({N_BOOT} draws, seed {SEED}); Table 2's r reproduced to {diff:.1e}.\n\n"
        "## Primary: Δr = r(consecutive) − r(skip-one), hard targets, one-month blocks\n\n" + md(prim.drop(columns=["reference", "block_months"])) +
        "\n\n## All contrasts (two-month blocks and grassland as sensitivity)\n\n" + md(c) +
        "\n\n## Table 2 with month-block intervals (r and slope through the origin) beside the iid ones\n\n" + md(t2) +
        "\n\n## Cross-track slope, consecutive six-day pairs 2021 + 2025, snow-free (vertical, hard targets)\n\n" + md(xt) + "\n",
        encoding="utf-8")
    print(json.dumps(v, indent=1))


if __name__ == "__main__":
    main()

"""X-127 — why pairs that skip one acquisition lose the mat's motion: ambiguity alone, or more noise than two short pairs
carry?

For each track and group of skip-one pairs at P6 (2021 twelve-day, 2022–2024 twenty-four-day, 2025 twelve-day), the
observed correlation of the pair phase with the laser (X-126's pairs: hard targets, corrected cube) is set against what
ambiguity alone would leave: each pair's own laser change plus the noise two consecutive pairs would carry (√2 times the
robust sd of the consecutive pairs' residual in the same years), wrapped to one cycle (``insar_wetlands.phase_sim``).
An observed correlation below the simulated range means the skip-one pair carries more noise than two short pairs.

Pre-registration: hub ``08_deliverables/insar_full_test/preregistration/X-127.yaml`` (checked before the run). Design
years only. Reads the hub; writes only to the hub (``08_deliverables/insar_full_test/skip_pairs/``). No field value in
this file.

    make run E=skip_pairs_x127
"""
from __future__ import annotations

import sys

import numpy as np
import pandas as pd

from insar_wetlands import phase_sim as ps
from insar_wetlands.paths import hub_root
from insar_wetlands.scoring import Protocol, assert_preregistered, guard_lockbox

HUB = hub_root()
FT = HUB / "08_deliverables" / "insar_full_test"
OUT = FT / "skip_pairs"
PREREG = FT / "preregistration" / "X-127.yaml"
PROTOCOL = FT / "protocol.yaml"
PAIRS = FT / "x065_hard_targets" / "pairs.csv"
EST = "s1_hard_targets"
GROUPS = (("2021", 12, 6), ("2022–2024", 24, 12), ("2025", 12, 6))   # (period, skip-one length, consecutive length)
N_SIM, SEED, MIN_PAIRS = 2000, 0, 10
QW = ps.WAVELENGTH_MM / 4


def robust_sd(x) -> float:
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return float(1.4826 * np.median(np.abs(x - np.median(x)))) if len(x) else float("nan")


def r(a, b) -> float:
    return float(np.corrcoef(a, b)[0, 1]) if len(a) > 2 and np.std(a) > 0 and np.std(b) > 0 else float("nan")


def simulated_r(truth: np.ndarray, sigma: float, rng: np.random.Generator) -> np.ndarray:
    """r of the wrapped (truth + noise) against the truth, over N_SIM realisations."""
    out = np.empty(N_SIM)
    for k in range(N_SIM):
        wrapped, _ = ps.simulate(truth, sigma, rng=rng)
        out[k] = r(truth, wrapped)
    return out


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    outs = [OUT / f for f in ("by_group.csv", "README.md")]
    assert_preregistered(PREREG, outs)
    p = pd.read_csv(PAIRS).dropna(subset=["d_laser_los_mm", EST])
    for c in ("t1", "t2"):
        p[c] = pd.to_datetime(p[c], format="ISO8601", utc=True)
    guard_lockbox(pd.concat([p.t1, p.t2]), Protocol.load(PROTOCOL))
    p["group"] = p.group.astype(str)
    rng = np.random.default_rng(SEED)
    rows = []
    for track in ("ascending", "descending"):
        for period, skip_days, cons_days in GROUPS:
            g = p[(p.track == track) & (p.group == period)]
            cons = g[(g.between == 0) & (g.revisit_days == cons_days)]
            skip = g[(g.between == 1) & (g.revisit_days == skip_days)]
            sigma = robust_sd(cons[EST] - cons.d_laser_los_mm)
            row = {"track": track, "period": period, "skip_days": skip_days, "consecutive_days": cons_days,
                   "consecutive_pairs": len(cons), "sigma_short_mm": sigma, "skip_pairs": len(skip),
                   "beyond_quarter_wave": float((skip.d_laser_los_mm.abs() > QW).mean()) if len(skip) else np.nan}
            if len(skip) < MIN_PAIRS or not np.isfinite(sigma):
                rows.append({**row, "verdict": "not judged (too few pairs)"})
                continue
            truth, obs = skip.d_laser_los_mm.to_numpy(), skip[EST].to_numpy()
            r_obs = r(truth, obs)
            sim = simulated_r(truth, np.sqrt(2) * sigma, rng)
            strict = simulated_r(truth, sigma, rng)
            within = np.abs(truth) <= QW
            lo, mid, hi = np.nanpercentile(sim, [2.5, 50, 97.5])
            rows.append({**row, "r_observed": r_obs, "r_sim_p2_5": lo, "r_sim_median": mid, "r_sim_p97_5": hi,
                         "share_sim_below_observed": float(np.mean(sim <= r_obs)),
                         "r_strict_p2_5": float(np.nanpercentile(strict, 2.5)),
                         "residual_sd_within_qw_mm": robust_sd(obs[within] - truth[within]),
                         "residual_ratio_to_sigma_short": robust_sd(obs[within] - truth[within]) / sigma,
                         "verdict": "beyond ambiguity" if r_obs < lo else "ambiguity enough"})
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "by_group.csv", index=False)
    cols = ["track", "period", "skip_days", "skip_pairs", "beyond_quarter_wave", "sigma_short_mm", "r_observed",
            "r_sim_p2_5", "r_sim_median", "r_sim_p97_5", "residual_ratio_to_sigma_short", "verdict"]
    d = t[[c for c in cols if c in t]].copy()
    for c in d.columns:
        if pd.api.types.is_float_dtype(d[c]):
            d[c] = d[c].map(lambda x: f"{x:.2f}")
    table = "\n".join(["| " + " | ".join(d.columns) + " |", "|" + "---|" * len(d.columns)] +
                      ["| " + " | ".join(str(x) for x in row) + " |" for row in d.itertuples(index=False)])
    (OUT / "README.md").write_text(f"""# X-127 — skip-one pairs at P6: ambiguity alone, or more noise than two short pairs? (generated)

Generated by `experiments/skip_pairs_x127.py` (`make run E=skip_pairs_x127`); do not edit. Pre-registration:
`../preregistration/X-127.yaml`. Pairs: `../x065_hard_targets/pairs.csv` (X-126). Null H_amb: each pair's laser change +
N(0, √2·σ_short), wrapped to one cycle; {N_SIM} realisations, seed {SEED}. σ_short = robust sd of the consecutive pairs'
residual (observed − laser) in the same track and years. "beyond ambiguity" = observed r below the simulated 2.5th
percentile.

{table}

Secondary columns (strict null with σ_short alone, share of simulations below the observed r) are in `by_group.csv`.
""", encoding="utf-8")
    print(t[["track", "period", "r_observed", "r_sim_p2_5", "verdict"]].to_string())


if __name__ == "__main__":
    sys.exit(main())

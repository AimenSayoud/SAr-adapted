"""X-126 — X-065's by-revisit comparison at P6 re-read on the hard-target reference, from the HyP3-corrected cube, by
the common scorer.

X-065 read the HyP3 crops directly, referenced the P6 phase to the grassland median and computed its own statistics.
The article's §3.1 references to hard targets throughout and §3.2 scores every laser comparison with one scorer, so
Table 2 and Fig. 4 are re-read here the way X-099 reads its pairs: plot P6, 40 m, the labelled cell of the corrected
cube (gold ``plots_pair``), both references, the cube's QC'd laser change along the line of sight. Groups as Table 2
(revisit × consecutive or skipping one acquisition × period), plus the consecutive six-day pairs of 2025.

Pre-registration: hub ``08_deliverables/insar_full_test/preregistration/X-126.yaml`` (checked before the run). The
2021 consecutive pairs must reproduce X-099's 2021 increments exactly, or the run stops. Design years only.
Reads the cube and the hub; writes only to the hub (``08_deliverables/insar_full_test/x065_hard_targets/``). No field
value in this file.

    make run E=x065_hard_targets_x126
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import p6_validation_x082 as X  # noqa: E402  (the X-082 readers: refs, acquisitions, laser)
from insar_wetlands.inversion.isbas import PHASE_TO_MM  # noqa: E402
from insar_wetlands.paths import hub_root  # noqa: E402
from insar_wetlands.scoring import Protocol, assert_preregistered, guard_lockbox, increment_stats  # noqa: E402

HUB = hub_root()
FT = HUB / "08_deliverables" / "insar_full_test"
OUT = FT / "x065_hard_targets"
PREREG = FT / "preregistration" / "X-126.yaml"
PROTOCOL = FT / "protocol.yaml"
X099_YEARS = FT / "reread_corrected" / "baseline_x099" / "years.csv"   # X-099 on the corrected cube (X-121)
REGULAR = ("2017_2019", "2020_2021", "2022_2024", "2025", "2026")       # not the annual or matched-water-table sets
N_BOOT, SEED = 2000, 0


def acquisitions_between(p: pd.DataFrame) -> pd.Series:
    """How many acquisitions of the pair's track lie strictly between its dates (0 = consecutive)."""
    acq = pd.read_csv(X.ASF, parse_dates=["date"])
    out = pd.Series(-1, index=p.index)
    for t in X.TRACKS:
        d = np.sort(acq[acq.track == t].date.values.astype("datetime64[D]"))
        s = p.track == t
        a = p.loc[s, "t1"].dt.tz_convert(None).values.astype("datetime64[D]")
        b = p.loc[s, "t2"].dt.tz_convert(None).values.astype("datetime64[D]")
        out[s] = np.searchsorted(d, b) - np.searchsorted(d, a, side="right")
    return out


def load_pairs() -> pd.DataFrame:
    """P6 pairs on the 40 m grid, consecutive or skipping one acquisition, inside the laser's span, both references."""
    cols = ["period", "track", "grid", "plot", "pair", "t1", "t2", "revisit_days", "platforms", "unw_1x1",
            "d_laser_los_mm"] + list(X.REFS.values())
    p = pd.read_csv(X.CUBE / "gold" / "plots_pair.csv.gz", usecols=cols, low_memory=False)
    p = p[(p["plot"] == "P6") & (p.grid == "40m") & p.period.astype(str).isin(REGULAR)].copy()
    p = p[~p.platforms.astype(str).str.contains("S1D")]                        # as v4: Sentinel-1D left out
    for c in ("t1", "t2"):
        p[c] = pd.to_datetime(p[c], format="ISO8601", utc=True)
    p["between"] = acquisitions_between(p)
    p = p[p.between.isin([0, 1])]
    for name, col in X.REFS.items():
        p[f"s1_{name}"] = (p.unw_1x1 - p[col]) * PHASE_TO_MM
    las = X.laser()[0]
    keep = []
    for track in X.TRACKS:
        lt = las[track].dropna()
        g = p[(p.track == track) & (p.t1 >= lt.index.min() - pd.Timedelta(days=13)) &
              (p.t2 <= lt.index.max() + pd.Timedelta(days=13))]
        keep.append(g)
    p = pd.concat(keep).sort_values(["track", "t1", "t2"]).reset_index(drop=True)
    p["year"] = p.t1.dt.year
    p["group"] = np.select(
        [p.year == 2021, p.year.between(2022, 2024), p.year == 2025], ["2021", "2022–2024", "2025"], "other")
    return p


def boot_r(a: np.ndarray, b: np.ndarray, rng: np.random.Generator) -> tuple[float, float]:
    """Bootstrap 95 % interval of Pearson's r over the pairs (as Table 2 was made: 2000 resamples)."""
    draws = []
    for _ in range(N_BOOT):
        i = rng.integers(0, len(a), len(a))
        if np.std(a[i]) > 0 and np.std(b[i]) > 0:
            draws.append(np.corrcoef(a[i], b[i])[0, 1])
    lo, hi = np.nanpercentile(draws, [2.5, 97.5])
    return float(lo), float(hi)


def by_revisit(p: pd.DataFrame) -> pd.DataFrame:
    rng = np.random.default_rng(SEED)
    rows = []
    for (track, group, between, days), g in p.groupby(["track", "group", "between", "revisit_days"]):
        if group == "other":
            continue
        for ref in X.REFS:
            m = g[["d_laser_los_mm", f"s1_{ref}"]].dropna()
            if len(m) < 5:
                continue
            a, b = m.d_laser_los_mm.to_numpy(), m[f"s1_{ref}"].to_numpy()
            lo, hi = boot_r(a, b, rng)
            rows.append({"track": track, "period": group, "revisit_days": int(days),
                         "pairs_kind": "consecutive" if between == 0 else "skips one", "reference": ref,
                         **increment_stats(a, b), "r_lo95": lo, "r_hi95": hi})
    return pd.DataFrame(rows)


def reproduce_x099(p: pd.DataFrame) -> pd.DataFrame:
    """The 2021 consecutive pairs against X-099's 2021 plain-chain increments (corrected cube): must be identical."""
    x = pd.read_csv(X099_YEARS)
    x = x[(x.year == 2021) & (x.estimate == "plain chain")]
    rows = []
    for track in X.TRACKS:
        g = p[(p.track == track) & (p.year == 2021) & (p.between == 0)]
        for ref in X.REFS:
            mine = increment_stats(g.d_laser_los_mm, g[f"s1_{ref}"])
            theirs = x[(x.track == track) & (x.reference == ref)].iloc[0]
            rows.append({"track": track, "reference": ref, "pairs": mine["pairs"], "pairs_x099": int(theirs.pairs),
                         "r": mine["pearson_increments"], "r_x099": float(theirs.pearson_increments),
                         "revisits": ",".join(str(int(v)) for v in sorted(g.revisit_days.unique()))})
    t = pd.DataFrame(rows)
    t["identical"] = (t.pairs == t.pairs_x099) & (np.abs(t.r - t.r_x099) < 1e-9)
    return t


def verdict(t: pd.DataFrame) -> pd.DataFrame:
    """The pre-registered verdict, per track, on the hard-target reading."""
    h = t[t.reference == "hard_targets"]
    rows = []
    for track in X.TRACKS:
        six = h[(h.track == track) & (h.period == "2021") & (h.revisit_days == 6) & (h.pairs_kind == "consecutive")]
        long_ = h[(h.track == track) & (h.period == "2022–2024") & (h.revisit_days == 24)]
        six_ok = len(six) == 1 and six.r_lo95.iloc[0] > 0
        long_ok = len(long_) == 1 and long_.r_lo95.iloc[0] <= 0 <= long_.r_hi95.iloc[0]
        rows.append({"track": track, "six_day_2021_above_zero": six_ok, "twenty_four_day_contains_zero": long_ok})
    v = pd.DataFrame(rows)
    v["verdict"] = np.where(v.six_day_2021_above_zero & v.twenty_four_day_contains_zero, "holds", "changes")
    return v


def md(df: pd.DataFrame, nd: int = 2) -> str:
    d = df.copy()
    for c in d.columns:
        if pd.api.types.is_float_dtype(d[c]):
            d[c] = d[c].map(lambda x: f"{x:.{nd}f}")
    return "\n".join(["| " + " | ".join(d.columns) + " |", "|" + "---|" * len(d.columns)] +
                     ["| " + " | ".join(str(x) for x in r) + " |" for r in d.itertuples(index=False)])


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    outs = [OUT / f for f in ("pairs.csv", "by_revisit.csv", "reproduction_x099.csv", "verdict.csv", "README.md")]
    assert_preregistered(PREREG, outs)
    p = load_pairs()
    guard_lockbox(pd.concat([p.t1, p.t2]), Protocol.load(PROTOCOL))
    rep = reproduce_x099(p)
    rep.to_csv(OUT / "reproduction_x099.csv", index=False)
    if not rep.identical.all():
        print(rep.to_string())
        sys.exit("stop (pre-registered): the 2021 consecutive pairs do not reproduce X-099; see reproduction_x099.csv")
    p.drop(columns=["unw_1x1"] + list(X.REFS.values())).to_csv(OUT / "pairs.csv", index=False)
    t = by_revisit(p)
    t.to_csv(OUT / "by_revisit.csv", index=False)
    v = verdict(t)
    v.to_csv(OUT / "verdict.csv", index=False)
    cols = ["track", "period", "revisit_days", "pairs_kind", "reference", "pairs", "pearson_increments", "r_lo95",
            "r_hi95", "rmse_increments_mm", "slope_through_origin", "within_quarter_wave"]
    side = t.pivot_table(index=["track", "period", "revisit_days", "pairs_kind"], columns="reference",
                         values="pearson_increments").reset_index()
    (OUT / "README.md").write_text(f"""# X-126 — X-065's by-revisit comparison on the hard-target reference (generated)

Generated by `experiments/x065_hard_targets_x126.py` (`make run E=x065_hard_targets_x126`); do not edit.
Pre-registration: `../preregistration/X-126.yaml`. Input: the corrected cube (X-120, X-122), gold `plots_pair`, P6,
40 m, the labelled cell; statistics by `insar_wetlands.scoring.increment_stats`; bootstrap {N_BOOT} resamples, seed
{SEED}. Design years only.

## Verdict (pre-registered rule, hard-target reference)

{md(v)}

## Reproduction check: 2021 consecutive pairs against X-099 (corrected cube)

{md(rep)}

## Pair phase against the laser at P6, by revisit and reference

{md(t[cols])}

## The two references side by side (Pearson r of the increments)

{md(side)}
""", encoding="utf-8")
    print(v.to_string())


if __name__ == "__main__":
    main()

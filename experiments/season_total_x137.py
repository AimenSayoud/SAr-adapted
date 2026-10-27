"""X-137 — how far P6 moves from spring to autumn each year: the season total by the laser and by the plain chain alone.

Descriptive. The plain chain of consecutive pairs and the laser are read exactly as X-099 reads them (the X-082
readers, ``insar_wetlands.scoring``'s chain and season views); each season view becomes one row: first and last
overpass of the snow-free window on one chain segment and one laser datum, the net change and the largest drop below
the start, both turned to vertical (cm). The stop rule checks every season view against X-099's levels.csv first.

Pre-registration: hub ``08_deliverables/insar_full_test/preregistration/X-137.yaml`` (checked before the run).
Reads the cube and the hub; writes only to the hub ``08_deliverables/insar_full_test/season_total_x137/``. No field
value in this file.

    make run E=season_total_x137
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import p6_validation_x082 as X  # noqa: E402  (the X-082 readers: pairs and laser, unchanged)

from insar_wetlands.paths import hub_root  # noqa: E402
from insar_wetlands.scoring import Protocol, assert_preregistered, chain, detrended_stats, level_views  # noqa: E402

HUB = hub_root()
FT = HUB / "08_deliverables" / "insar_full_test"
OUT = FT / "season_total_x137"
PREREG = FT / "preregistration" / "X-137.yaml"
PROTOCOL = FT / "protocol.yaml"
X099 = FT / "reread_corrected" / "baseline_x099" / "levels.csv"
SKIP = {"season 2022 A"}                      # before the 2022 re-levelling: the pre-registration reads 2022 after it


def verdict(diff_mm: float) -> str:
    d = abs(diff_mm)
    return "within 1 cm" if d <= 10 else "within 2 cm" if d <= 20 else "off by more than 2 cm"


def season_rows(track: str, ref: str, p: pd.DataFrame, las_los: dict, lev: pd.Series, prot: Protocol,
                x099: pd.DataFrame) -> list[dict]:
    lt = las_los[track].dropna()
    cos = float((las_los[track] / lev).dropna().median())          # P6's line-of-sight up-component, this track
    g = p[(p.track == track) & (p.t1 >= lt.index.min() - pd.Timedelta(days=13)) &
          (p.t2 <= lt.index.max() + pd.Timedelta(days=13))].sort_values("t1")
    c = chain(g.assign(est_mm=g[f"s1_{ref}"]), "est_mm")
    c["truth"] = lt.reindex(pd.DatetimeIndex(c.time).round("h")).to_numpy()
    rows = []
    for view, d in level_views(c, prot):
        if not view.startswith("season") or len(d) < prot.min_epochs:
            continue
        ref_row = x099[(x099.view == view) & (x099.track == track) & (x099.reference == ref)
                       & (x099.estimate == "plain chain")]
        st = detrended_stats(d.time, d.truth, d.est)
        if len(ref_row) != 1 or int(ref_row.epochs.iloc[0]) != len(d) or \
                abs(float(ref_row.pearson_detrended.iloc[0]) - st["pearson_detrended"]) > 1e-9:
            raise SystemExit(f"stop: {view} {track} {ref} does not reproduce X-099 ({len(d)} epochs, "
                             f"r {st['pearson_detrended']:.12f} against {ref_row.to_dict('records')})")
        if view in SKIP:
            continue
        run = d[d.seg == d.seg.value_counts().idxmax()].copy()           # one unbroken chain segment
        run["truth"] -= run.truth.iloc[0]
        run["est"] -= run.est.iloc[0]
        z_l, z_r = run.truth / cos / 10, run.est / cos / 10              # vertical, cm
        i_l, i_r = int(np.argmin(z_l.to_numpy())), int(np.argmin(z_r.to_numpy()))
        net_l, net_r = float(z_l.iloc[-1]), float(z_r.iloc[-1])
        drop_l, drop_r = float(z_l.min()), float(z_r.min())
        rows.append({
            "year": int(view.split()[1]), "track": track, "reference": ref, "view": view,
            "start": run.time.iloc[0].date(), "end": run.time.iloc[-1].date(), "overpasses": len(run),
            "partial": "partial season (after the 2022 re-levelling)" if view == "season 2022 B" else
                       ("chain gap: longest run" if d.seg.nunique() > 1 else ""),
            "laser_net_cm": net_l, "radar_net_cm": net_r, "net_diff_cm": net_r - net_l,
            "net_verdict": verdict(10 * (net_r - net_l)),
            "laser_drop_cm": drop_l, "laser_drop_date": run.time.iloc[i_l].date(),
            "radar_drop_cm": drop_r, "radar_drop_date": run.time.iloc[i_r].date(),
            "drop_diff_cm": drop_r - drop_l, "drop_verdict": verdict(10 * (drop_r - drop_l)),
            "revisit": "six-day" if int(view.split()[1]) in (2021, 2025) else "twelve-day"})
    return rows


def laser_daily(lev: pd.Series, prot: Protocol) -> pd.DataFrame:
    """Secondary: the laser alone over the snow-free window at daily resolution (first and last valid day)."""
    day = lev.dropna().resample("D").mean().dropna()
    datum = np.asarray(prot.datum(pd.Series(day.index)))
    m0, m1 = prot.season_months
    out = []
    for y in sorted(set(day.index.year)):
        for k in sorted(set(datum)):
            s = day[(day.index.year == y) & day.index.month.isin(range(m0, m1 + 1)) & (datum == k)]
            if len(s) < 30 or f"season {y} {k}" in SKIP:
                continue
            z = (s - s.iloc[0]) / 10
            out.append({"year": y, "datum": k, "first_day": s.index[0].date(), "last_day": s.index[-1].date(),
                        "days": len(s), "net_cm": float(z.iloc[-1]), "drop_cm": float(z.min()),
                        "drop_date": z.idxmin().date(), "range_cm": float(z.max() - z.min())})
    return pd.DataFrame(out)


def readme(t: pd.DataFrame, daily: pd.DataFrame) -> str:
    prim = t[t.reference == "hard_targets"]
    def f(v): return f"{v:+.1f}"
    lines = ["| year | track | dates | overpasses | laser net | radar net | difference | reading | laser largest drop | "
             "radar largest drop | reading |", "|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in prim.itertuples():
        lines.append(f"| {r.year}{' *' if r.partial else ''} | {r.track} | {r.start} → {r.end} | {r.overpasses} | "
                     f"{f(r.laser_net_cm)} | {f(r.radar_net_cm)} | {f(r.net_diff_cm)} | {r.net_verdict} | "
                     f"{f(r.laser_drop_cm)} ({r.laser_drop_date}) | {f(r.radar_drop_cm)} ({r.radar_drop_date}) | "
                     f"{r.drop_verdict} |")
    head = []
    for rev in ("six-day", "twelve-day"):
        s = prim[prim.revisit == rev]
        for col, what in (("net_verdict", "net change"), ("drop_verdict", "largest drop")):
            head.append(f"- {rev} years, {what}: within 1 cm in {int((s[col] == 'within 1 cm').sum())} of {len(s)}, "
                        f"within 2 cm in {int(s[col].isin(['within 1 cm', 'within 2 cm']).sum())} of {len(s)}")
    dl = ["| year | datum | first → last day | days | net | largest drop (date) | range |", "|---|---|---|---|---|---|---|"]
    for r in daily.itertuples():
        dl.append(f"| {r.year} | {r.datum} | {r.first_day} → {r.last_day} | {r.days} | {f(r.net_cm)} | "
                  f"{f(r.drop_cm)} ({r.drop_date}) | {r.range_cm:.1f} |")
    return f"""# X-137 — how far P6 moves from spring to autumn (generated)

Generated by `experiments/season_total_x137.py` (`make run E=season_total_x137`); do not edit. Pre-registration:
`../preregistration/X-137.yaml`. The plain chain and the laser of X-099 (stop rule passed: every season view
reproduces X-099's levels.csv). Vertical, cm, positive = up; each season starts at zero at its first overpass.
Hard-target reference (primary); grassland in `season_totals.csv`. `*` = partial season (2022: after the laser's
re-levelling of 5 July – 8 August).

## Season totals at the overpasses

{chr(10).join(lines)}

Reading (pre-registered: |radar − laser| ≤ 1 cm, ≤ 2 cm, more):

{chr(10).join(head)}

## The laser alone, daily, over the snow-free window (May–November)

What the overpass sampling leaves out: first and last valid laser day of the window on one datum.

{chr(10).join(dl)}
"""


def main() -> None:
    assert_preregistered(PREREG, [OUT / "season_totals.csv"])
    prot = Protocol.load(PROTOCOL)
    p = X.load_pairs()
    las_los, lev = X.laser()
    x099 = pd.read_csv(X099)
    rows = [r for track in X.TRACKS for ref in X.REFS for r in season_rows(track, ref, p, las_los, lev, prot, x099)]
    t = pd.DataFrame(rows).sort_values(["reference", "year", "track"])
    daily = laser_daily(lev, prot)
    OUT.mkdir(parents=True, exist_ok=True)
    t.to_csv(OUT / "season_totals.csv", index=False, float_format="%.4f")
    daily.to_csv(OUT / "laser_daily_season.csv", index=False, float_format="%.4f")
    (OUT / "README.md").write_text(readme(t, daily))
    print(readme(t, daily))


if __name__ == "__main__":
    main()

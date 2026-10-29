"""X-142 — are whole-cycle unwrapping offsets between P6 and the hard targets behind the winter level loss? (R4 M1)

For each consecutive pair at P6 (the X-082 reader, unchanged): u = unwrapped P6 − unwrapped hard targets, w = the same
for the wrapped phase, wrapped into (−π, π], k = round((u − w) / 2π). The plain chain cumulates u (X-099); the
cycle-free chain cumulates u − 2πk (same pairs, same fractional phase, every whole cycle removed). Both are scored on
the laser winters by the common scorer; the laser says, pair by pair, whether a cycle can be motion (|Δ| ≥ λ/4).

Pre-registration: hub ``08_deliverables/insar_full_test/preregistration/X-142.yaml`` (checked before the run).
Reads the cube and the hub; writes only to the hub: ``08_deliverables/insar_full_test/winter_cycles_x142/``. No field
value in this file.

    make run E=winter_cycles_x142
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import p6_validation_x082 as X  # noqa: E402  (the X-082 readers: pairs and laser, unchanged)

from insar_wetlands.paths import cube_root, hub_root  # noqa: E402
from insar_wetlands.scoring import Protocol, assert_preregistered, chain, level_views, detrended_stats, score  # noqa: E402

HUB, CUBE = hub_root(), cube_root()
FT = HUB / "08_deliverables" / "insar_full_test"
OUT = FT / "winter_cycles_x142"
PREREG = FT / "preregistration" / "X-142.yaml"
PROTOCOL = FT / "protocol.yaml"
BASE = FT / "reread_corrected" / "baseline_x099" / "winters.csv"
QUARTER_WAVE_MM = 55.465763 / 4          # Sentinel-1 C-band, line of sight
CYCLE_MM = abs(X.PHASE_TO_MM) * 2 * np.pi
REF_COLS = {"hard_targets": ("ref_hard_unw", "ref_hard_wrapped"),
            "grassland": ("ref_grassland_unw", "ref_grassland_wrapped")}


def pairs() -> pd.DataFrame:
    """X-082's consecutive P6 pairs, with the wrapped phase and the component mask beside them."""
    p = X.load_pairs()
    extra = ["track", "t1", "t2", "wrapped_1x1", "conncomp_1x1", "coh_1x1"] + [w for _, w in REF_COLS.values()]
    q = pd.read_csv(CUBE / "gold" / "plots_pair.csv.gz", low_memory=False,
                    usecols=extra + ["plot", "grid"])
    q = q[(q["plot"] == "P6") & (q.grid == "40m")].drop(columns=["plot", "grid"])
    for c in ("t1", "t2"):
        q[c] = pd.to_datetime(q[c], format="ISO8601", utc=True)
    p = p.merge(q.drop_duplicates(["track", "t1", "t2"]), on=["track", "t1", "t2"], how="left", validate="1:1")
    for ref, (ucol, wcol) in REF_COLS.items():
        u = p.unw_1x1 - p[ucol]
        w = np.angle(np.exp(1j * (p.wrapped_1x1 - p[wcol])))
        p[f"k_{ref}"] = np.round((u - w) / (2 * np.pi))
        p[f"free_{ref}"] = (u - 2 * np.pi * p[f"k_{ref}"]) * X.PHASE_TO_MM
    return p


def winter_spans(sc_winters: pd.DataFrame) -> list[tuple[str, pd.Timestamp, pd.Timestamp]]:
    """(winter, first day, day after the last) of each scored winter step, so pairs dated within it are caught."""
    return [(r["winter"], pd.Timestamp(r["from"], tz="UTC"), pd.Timestamp(r["to"], tz="UTC") + pd.Timedelta(days=1))
            for _, r in sc_winters.iterrows()]


def adjudicate(g: pd.DataFrame, ref: str) -> pd.Series:
    """Per pair with k ≠ 0: against (|laser| < λ/4: the cycle cannot be motion), with, or not adjudicable."""
    k, las = g[f"k_{ref}"], g.d_laser_los_mm
    lab = pd.Series("k = 0", index=g.index)
    nz = k != 0
    lab[nz & las.isna()] = "not adjudicable (no laser)"
    lab[nz & las.notna() & (las.abs() < QUARTER_WAVE_MM)] = "cycle against the laser"
    lab[nz & las.notna() & (las.abs() >= QUARTER_WAVE_MM)] = "cycle with the laser (|laser| ≥ λ/4)"
    return lab


def run() -> dict[str, pd.DataFrame]:
    prot = Protocol.load(PROTOCOL)
    p, las = pairs(), X.laser()[0]
    win, adj, levs, pk = [], [], [], []
    for track in X.TRACKS:
        lt = las[track].dropna()
        g = p[(p.track == track) & (p.t1 >= lt.index.min() - pd.Timedelta(days=13)) &
              (p.t2 <= lt.index.max() + pd.Timedelta(days=13))].sort_values("t1").copy()
        for ref in REF_COLS:
            g[f"adj_{ref}"] = adjudicate(g, ref)
            against = g[f"adj_{ref}"] == "cycle against the laser"
            g[f"adjud_{ref}"] = np.where(against, g[f"free_{ref}"], g[f"s1_{ref}"])
            chains = {"plain chain": f"s1_{ref}", "cycle-free chain": f"free_{ref}",
                      "laser-adjudicated chain (diagnostic)": f"adjud_{ref}"}
            spans = None
            for name, col in chains.items():
                r = score(g.assign(est_mm=g[col], truth_inc_mm=g.d_laser_los_mm), lt, prot)
                w = r["winters"].assign(estimate=name, track=track, reference=ref)
                win.append(w)
                levs.append(r["levels"].assign(estimate=name, track=track, reference=ref))
                if name == "plain chain":
                    spans = winter_spans(r["winters"])
            for wname, a, b in spans:
                s = g[(g.t1 >= a) & (g.t2 <= b)]
                for lab, n in s[f"adj_{ref}"].value_counts().items():
                    adj.append({"track": track, "reference": ref, "span": wname, "label": lab, "pairs": int(n),
                                "net_cycles": int(s.loc[s[f"adj_{ref}"] == lab, f"k_{ref}"].sum())})
            for lab, n in g[f"adj_{ref}"].value_counts().items():
                adj.append({"track": track, "reference": ref, "span": "laser years, all consecutive pairs",
                            "label": lab, "pairs": int(n),
                            "net_cycles": int(g.loc[g[f"adj_{ref}"] == lab, f"k_{ref}"].sum())})
        pk.append(g)
    return {"winters": pd.concat(win, ignore_index=True), "adjudication": pd.DataFrame(adj),
            "levels": pd.concat(levs, ignore_index=True), "pairs_k": pd.concat(pk, ignore_index=True)}


def stop_check(w: pd.DataFrame) -> float:
    b = pd.read_csv(BASE)
    b = b[(b.reference == "hard_targets") & (b.estimate == "plain chain")]
    m = w[(w.reference == "hard_targets") & (w.estimate == "plain chain")]
    j = b.merge(m, on=["track", "winter"], suffixes=("_x099", "_x142"), validate="1:1")
    if len(j) != len(b):
        raise SystemExit(f"stop: winters differ from X-099 ({len(j)} matched of {len(b)})")
    d = float((j.insar_change_mm_x099 - j.insar_change_mm_x142).abs().max())
    if d > 1e-3:
        raise SystemExit(f"stop: plain chain winter steps differ from X-099 by up to {d:.4f} mm")
    return d


def verdict(w: pd.DataFrame) -> dict:
    h = w[w.reference == "hard_targets"].assign(ae=lambda d: d.difference_mm.abs())
    mae = h.groupby(["estimate", "track"]).ae.mean().unstack("track")
    red = {t: float(1 - mae.loc["cycle-free chain", t] / mae.loc["plain chain", t]) for t in X.TRACKS}
    if all(r >= 0.5 for r in red.values()):
        v = "(a) whole-cycle unwrapping offsets explain most of the winter level loss"
    elif any(r >= 0.25 for r in red.values()):
        v = "(b) whole-cycle offsets explain part of the winter level loss; accumulated noise the rest"
    else:
        v = "(c) whole-cycle offsets do not explain the winter level loss; accumulated noise remains the account"
    return {"mae_mm": {e: {t: float(mae.loc[e, t]) for t in X.TRACKS} for e in mae.index},
            "reduction": red, "verdict": v,
            "winters_per_track": {t: int((h.track == t).sum() / h.estimate.nunique()) for t in X.TRACKS}}


def md(df: pd.DataFrame, nd: int = 1) -> str:
    d = df.copy()
    for c in d.columns:
        if pd.api.types.is_float_dtype(d[c]):
            d[c] = d[c].round(nd)
    head = "| " + " | ".join(map(str, d.columns)) + " |\n|" + "---|" * len(d.columns) + "\n"
    return head + "\n".join("| " + " | ".join(map(str, r)) + " |" for r in d.itertuples(index=False))


def readme(t: dict[str, pd.DataFrame], v: dict, d: float) -> str:
    w = t["winters"]
    hw = w[w.reference == "hard_targets"][["track", "winter", "estimate", "laser_change_mm", "insar_change_mm",
                                           "difference_mm"]]
    a = t["adjudication"]
    ha = a[a.reference == "hard_targets"].pivot_table(index=["track", "span"], columns="label", values="pairs",
                                                     aggfunc="sum", fill_value=0).reset_index()
    return f"""# X-142 — whole cycles at P6 and the winter level loss (R4 M1)

Generated by `05_code/SAr-adapted/experiments/winter_cycles_x142.py`; pre-registered in
`preregistration/X-142.yaml`. Consecutive pairs at P6 (40 m), hard-target reference unless stated; a cycle is
{CYCLE_MM:.1f} mm in the line of sight, a quarter wavelength {QUARTER_WAVE_MM:.2f} mm.

Stop check passed: the plain chain reproduces X-099's winter steps to {d:.1e} mm.

## Verdict

**{v["verdict"]}.** Mean absolute winter-step error (mm): plain chain {v["mae_mm"]["plain chain"]["ascending"]:.1f}
ascending, {v["mae_mm"]["plain chain"]["descending"]:.1f} descending; cycle-free chain
{v["mae_mm"]["cycle-free chain"]["ascending"]:.1f} and {v["mae_mm"]["cycle-free chain"]["descending"]:.1f}
(reduction {v["reduction"]["ascending"]:.0%} and {v["reduction"]["descending"]:.0%}); laser-adjudicated chain
(diagnostic, uses the truth) {v["mae_mm"]["laser-adjudicated chain (diagnostic)"]["ascending"]:.1f} and
{v["mae_mm"]["laser-adjudicated chain (diagnostic)"]["descending"]:.1f}.

## Winter steps (estimate − laser)

{md(hw)}

## Pairs with a whole cycle, read against the laser

{md(ha, 0)}

## Files

| file | content |
|---|---|
| `winters.csv` | winter steps of the three chains, both references |
| `adjudication.csv` | pairs with k ≠ 0 by winter span and label, with their net cycles |
| `levels.csv` | the scorer's level views of the three chains |
| `pairs_k.csv` | every consecutive pair: u, w, k, the laser increment, the label |
| `verdict.json` | the registered verdict |
"""


def main() -> None:
    assert_preregistered(PREREG, [OUT / "verdict.json"])
    t = run()
    d = stop_check(t["winters"])
    v = verdict(t["winters"])
    OUT.mkdir(parents=True, exist_ok=True)
    for name, df in t.items():
        df.to_csv(OUT / f"{name}.csv", index=False)
    (OUT / "verdict.json").write_text(json.dumps({**v, "stop_check_max_diff_mm": d}, indent=1), encoding="utf-8")
    (OUT / "README.md").write_text(readme(t, v, d), encoding="utf-8")
    print(json.dumps(v, indent=1))


if __name__ == "__main__":
    main()

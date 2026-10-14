"""X-117 — does the mat follow its water table with memory? Buoyancy on a smoothed water table plus a laser-drift rate.

Fitted once on X-116's within-season intervals of every length (equal weight per length), year-blocked; checked for
interval invariance and on the X-099 winter steps. Rule: ``08_deliverables/insar_full_test/preregistration/X-117.yaml``.
Writes only to the hub: ``08_deliverables/insar_full_test/memory_response/``.

    make run E=memory_response_x117
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import spams_baseline_x112 as B  # noqa: E402
import winter_forward_x111 as W  # noqa: E402
import winter_physics_x116 as P  # noqa: E402

from insar_wetlands import forward  # noqa: E402
from insar_wetlands.paths import hub_root  # noqa: E402
from insar_wetlands.scoring import Protocol, assert_preregistered, guard_lockbox  # noqa: E402

HUB = hub_root()
FT = HUB / "08_deliverables" / "insar_full_test"
OUT = FT / "memory_response"
PREREG = FT / "preregistration" / "X-117.yaml"
VARIANTS = {"memory + drift": dict(drift=True), "memory, no drift": dict(drift=False),
            "no memory, drift": dict(drift=True, taus=(1,)), "plain g (no memory, no drift)": dict(drift=False, taus=(1,))}


def pooled_intervals(d: pd.DataFrame, prot: Protocol) -> pd.DataFrame:
    parts = []
    for L in P.LENGTHS:
        x = P.intervals(d, L, prot)
        x["L"] = L
        parts.append(x)
    x = pd.concat(parts, ignore_index=True)
    x["i0"] = d.index.get_indexer(x.t0)
    x["i1"] = d.index.get_indexer(x.t1)
    x["w"] = 1.0 / x.groupby("L").L.transform("size")
    return x[(x.i0 >= 0) & (x.i1 >= 0)]


def main() -> None:
    assert_preregistered(PREREG, [OUT / "verdict.csv"])
    OUT.mkdir(parents=True, exist_ok=True)
    prot = Protocol.load(FT / "protocol.yaml")
    d = B.table()
    d = d[d.year.between(2020, 2025)]
    guard_lockbox(d.index[d.year >= 2021], prot)
    w = d.wtd_cm.to_numpy()
    x = pooled_intervals(d, prot)
    preds = {k: np.full(len(x), np.nan) for k in VARIANTS}
    folds = []
    for y in P.YEARS:
        tr, te = (x.year != y).to_numpy(), (x.year == y).to_numpy()
        for name, kw in VARIANTS.items():
            par = forward.fit(w, x.i0[tr], x.i1[tr], x.dh_mm[tr], weights=x.w[tr], **kw)
            preds[name][te] = forward.predict(w, par, x.i0[te], x.i1[te])
            folds.append({"held_out_year": y, "variant": name, **{k: par[k] for k in ("tau", "g", "c", "n")}})
    inv = []
    for name, p in preds.items():
        for L in P.LENGTHS:
            m = (x.L == L).to_numpy() & np.isfinite(p)
            t, e = x.dh_mm.to_numpy()[m], p[m]
            inv.append({"variant": name, "length_days": L, "intervals": int(m.sum()),
                        "slope_laser_on_model": float(np.sum(e * t) / np.sum(e * e)) if m.any() else np.nan,
                        "rmse_mm": float(np.sqrt(np.mean((e - t) ** 2))), "pearson": float(np.corrcoef(e, t)[0, 1])})
    inv = pd.DataFrame(inv)
    base = pd.read_csv(FT / "baseline" / "winters.csv", parse_dates=["from", "to"])
    base = base[(base.estimate == "plain chain") & (base.reference == "hard_targets")]
    st = []
    los = W.los_factor()
    for _, r in base.iterrows():
        yy = int(r.winter[:4])
        if not (np.isfinite(W.at(d.wtd_cm, r["from"])) and np.isfinite(W.at(d.wtd_cm, r["to"]))):
            continue
        i0, i1 = d.index.get_indexer([r["from"]])[0], d.index.get_indexer([r["to"]])[0]
        tr = ~x.year.isin([yy, yy + 1]).to_numpy()
        rec = {"winter": r.winter, "track": r.track, "laser_step_mm": r.laser_change_mm, "chain_step_mm": r.insar_change_mm}
        for name, kw in VARIANTS.items():
            par = forward.fit(w, x.i0[tr], x.i1[tr], x.dh_mm[tr], weights=x.w[tr], **kw)
            rec[f"{name}"] = float(forward.predict(w, par, [i0], [i1])[0]) * los[r.track]   # vertical → LOS
        st.append(rec)
    steps = pd.DataFrame(st)
    g6 = pd.read_csv(FT / "winter_physics" / "winter_steps.csv")
    steps = steps.merge(g6[["winter", "track", "pred_g6_mm", "pred_g24_mm"]], on=["winter", "track"], how="left")
    mae = lambda col: float(np.mean(np.abs(steps[col] - steps.laser_step_mm)))
    wm = pd.DataFrame([{"estimate": c, "mae_mm": mae(c)} for c in list(VARIANTS) + ["pred_g6_mm", "pred_g24_mm", "chain_step_mm"]])
    full = inv[inv.variant == "memory + drift"].slope_laser_on_model
    verdict = pd.DataFrame([{"steps": len(steps), "mae_memory_drift_mm": mae("memory + drift"), "mae_chain_mm": mae("chain_step_mm"),
                             "mae_g6_mm": mae("pred_g6_mm"),
                             "carries_the_winter": bool(mae("memory + drift") <= 0.5 * mae("chain_step_mm")
                                                        and mae("memory + drift") < mae("pred_g6_mm")),
                             "slope_min": float(full.min()), "slope_max": float(full.max()),
                             "removes_length_dependence": bool(full.between(0.75, 1.33).all())}])
    for name, df in (("invariance", inv), ("folds", pd.DataFrame(folds)), ("winter_steps", steps), ("winter_mae", wm),
                     ("verdict", verdict)):
        df.to_csv(OUT / f"{name}.csv", index=False)
    (OUT / "README.md").write_text(readme(inv, pd.DataFrame(folds), steps, wm, verdict), encoding="utf-8")
    print(verdict.to_string(index=False))
    print(f"wrote {OUT}")


def _md(df: pd.DataFrame, nd: int = 2) -> str:
    d = df.copy()
    for col in d.columns:
        if pd.api.types.is_float_dtype(d[col]):
            d[col] = d[col].map(lambda x: "—" if pd.isna(x) else f"{x:.{nd}f}")
    return "\n".join(["| " + " | ".join(map(str, d.columns)) + " |", "|" + "---|" * len(d.columns)] +
                     ["| " + " | ".join(map(str, r)) + " |" for r in d.itertuples(index=False)])


def readme(inv, folds, steps, wm, verdict) -> str:
    piv = inv.pivot(index="length_days", columns="variant", values="slope_laser_on_model").reset_index()
    fm = folds.groupby("variant")[["tau", "g", "c"]].median().reset_index()
    return f"""# X-117 — a memory-response forward model of the mat (generated)

Generated by `experiments/memory_response_x117.py` (`make run E=memory_response_x117`); do not edit. Pre-registration:
`../preregistration/X-117.yaml`. Δh = g·ΔW_f + c·days, W_f the logger water table with exponential memory τ; fitted on
X-116's within-season intervals (6–96 days, equal weight per length), year-blocked. Interval invariance: slope of the
laser change on the model's held-out prediction, by length (1 = right at that length). Winter: X-099's logger-covered
steps, fits without the winter's two years; predictions projected to each track's line of sight.

## Verdict

{_md(verdict)}

## Slope of the laser change on the model, by interval length (held-out years)

{_md(piv)}

## Fitted parameters (median over held-out years)

{_md(fm, 3)}

## Winter steps — mean absolute error (pred_g6 / pred_g24: X-116's coefficients)

{_md(wm, 1)}

## Each winter step (mm)

{_md(steps, 1)}
"""


if __name__ == "__main__":
    main()

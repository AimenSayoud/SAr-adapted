"""X-113 — spring re-levelling: the 6-day chain at P6 with each winter's level step replaced by an anchor's.

Within seasons the chain keeps its increments; across each winter (X-099's definition) the step is set to a declared
anchor's (projected to line of sight): A1 buoyancy with the 6-day g, A2 SPAMS (both X-111), A3 buoyancy with the 24-day
g (X-116, post-hoc). Scored by the common scorer's level views. Rule:
``08_deliverables/insar_full_test/preregistration/X-113.yaml``. Writes only to the hub:
``08_deliverables/insar_full_test/relevel/``.

    make run E=relevel_x113
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import p6_validation_x082 as X  # noqa: E402  (pairs and the LOS laser, as X-099)

from insar_wetlands.paths import hub_root  # noqa: E402
from insar_wetlands.scoring import Protocol, assert_preregistered, chain, detrended_stats, level_views  # noqa: E402

HUB = hub_root()
FT = HUB / "08_deliverables" / "insar_full_test"
OUT = FT / "relevel"
PREREG = FT / "preregistration" / "X-113.yaml"


def anchors() -> pd.DataFrame:
    """Anchor steps (mm, LOS) per winter and track, from the corrected X-111 and X-116 tables."""
    f = pd.read_csv(FT / "winter_forward" / "steps.csv")
    a1 = f[f.candidate == "logger water table × g"][["winter", "track", "predicted_step_mm"]].rename(columns={"predicted_step_mm": "A1 buoyancy, 6-day g"})
    a2 = f[f.candidate == "SPAMS"][["winter", "track", "predicted_step_mm"]].rename(columns={"predicted_step_mm": "A2 SPAMS"})
    p = pd.read_csv(FT / "winter_physics" / "winter_steps.csv")[["winter", "track", "pred_g24_mm"]]
    p = p.rename(columns={"pred_g24_mm": "A3 buoyancy, 24-day g (post-hoc)"})
    return a1.merge(a2, on=["winter", "track"], how="outer").merge(p, on=["winter", "track"], how="left")


def winters(c: pd.DataFrame, prot: Protocol) -> list[tuple[str, int, int]]:
    """(winter, index of the last in-season epoch of y, index of the first spring-window epoch of y+1), same datum."""
    z = c.dropna(subset=["truth"])
    z = z.assign(datum=prot.datum(z.time))
    m0, m1 = prot.season_months
    s0, s1 = prot.spring_window
    out = []
    for y in sorted(z.time.dt.year.unique()):
        a = z[(z.time.dt.year == y) & z.time.dt.month.between(m0, m1)]
        b = z[(z.time.dt.year == y + 1) & z.time.dt.month.between(s0, s1)]
        if len(a) and len(b) and a.datum.iloc[-1] == b.datum.iloc[0]:
            out.append((f"{y}/{str(y + 1)[2:]}", a.index[-1], b.index[0]))
    return out


def relevel(c: pd.DataFrame, steps: dict[str, float], wins) -> pd.DataFrame:
    c = c.copy()
    for w, ia, ib in wins:
        v = steps.get(w)
        if v is None or not np.isfinite(v):
            continue                                    # no anchor for this winter: the chain's own step is kept
        delta = v - (c.est.loc[ib] - c.est.loc[ia])
        c.loc[c.index >= ib, "est"] += delta
    return c


def main() -> None:
    assert_preregistered(PREREG, [OUT / "verdict.csv"])
    OUT.mkdir(parents=True, exist_ok=True)
    prot = Protocol.load(FT / "protocol.yaml")
    p, las = X.load_pairs(), X.laser()[0]
    anc = anchors()
    names = [c for c in anc.columns if c not in ("winter", "track")]
    rows = []
    for track in X.TRACKS:
        lt = las[track].dropna()
        g = p[(p.track == track) & (p.t1 >= lt.index.min() - pd.Timedelta(days=13)) &
              (p.t2 <= lt.index.max() + pd.Timedelta(days=13))].sort_values("t1")
        for ref in X.REFS:
            gg = g.assign(est_mm=g[f"s1_{ref}"])
            for col in ("t1", "t2"):
                gg[col] = pd.to_datetime(gg[col], utc=True)
            c = chain(gg, "est_mm")
            c["truth"] = lt.reindex(pd.DatetimeIndex(c.time).round("h")).to_numpy()
            wins = winters(c, prot)
            a = anc[anc.track == track].set_index("winter")
            variants = {"plain chain": c}
            for n in names:
                variants[n] = relevel(c, a[n].to_dict(), wins)
            for vname, cv in variants.items():
                for view, d in level_views(cv, prot):
                    if view in ("full span", "datum A", "datum B") and len(d) >= prot.min_epochs:
                        rows.append({"track": track, "reference": ref, "estimate": vname, "view": view,
                                     **detrended_stats(d.time, d.truth, d.est)})
    lev = pd.DataFrame(rows)
    verdict = decide(lev, names)
    for name, df in (("levels", lev), ("verdict", verdict), ("anchors", anc)):
        df.to_csv(OUT / f"{name}.csv", index=False)
    (OUT / "README.md").write_text(readme(lev, verdict, anc), encoding="utf-8")
    print(verdict.to_string(index=False))
    print(f"wrote {OUT}")


def decide(lev: pd.DataFrame, names) -> pd.DataFrame:
    b = lev[(lev.view == "datum B") & (lev.reference == "hard_targets")].set_index(["estimate", "track"]).rmse_detrended_mm
    rows = []
    for n in names:
        gain = {t: 1 - b.loc[(n, t)] / b.loc[("plain chain", t)] for t in X.TRACKS}
        rows.append({"anchor": n, **{f"rmse_{t[:4]}_mm": b.loc[(n, t)] for t in X.TRACKS},
                     **{f"plain_{t[:4]}_mm": b.loc[("plain chain", t)] for t in X.TRACKS},
                     **{f"gain_{t[:4]}": gain[t] for t in X.TRACKS},
                     "re_levels": bool(all(v >= 0.25 for v in gain.values())),
                     "post_hoc": "post-hoc" in n})
    return pd.DataFrame(rows)


def _md(df: pd.DataFrame, nd: int = 2) -> str:
    d = df.copy()
    for col in d.columns:
        if pd.api.types.is_float_dtype(d[col]):
            d[col] = d[col].map(lambda x: "—" if pd.isna(x) else f"{x:.{nd}f}")
    return "\n".join(["| " + " | ".join(map(str, d.columns)) + " |", "|" + "---|" * len(d.columns)] +
                     ["| " + " | ".join(map(str, r)) + " |" for r in d.itertuples(index=False)])


def readme(lev, verdict, anc) -> str:
    h = lev[lev.reference == "hard_targets"][["track", "estimate", "view", "epochs", "pearson_detrended", "rmse_detrended_mm"]]
    return f"""# X-113 — spring re-levelling of the 6-day chain (generated)

Generated by `experiments/relevel_x113.py` (`make run E=relevel_x113`); do not edit. Pre-registration:
`../preregistration/X-113.yaml`. P6, hard-target and grassland references; within seasons the chain's own increments,
across each winter the step of a declared anchor (line of sight). A3's interval was chosen after looking (X-116): its
result is design-year evidence only. An anchor re-levels if datum B's detrended RMSE falls by at least 25% on both
tracks.

## Verdict

{_md(verdict)}

## Levels, hard-target reference

{_md(h)}

## Anchor steps (mm, line of sight)

{_md(anc, 1)}
"""


if __name__ == "__main__":
    main()

"""X-135: does the procedure that produced v4 — choose the architecture by X-119's rule, then fit it — carry the mat's
level across a winter it has never seen?

The laser archive that X-107 was to score on will not come (hub DECISIONS 2026-10-26). X-119 scored each winter with
fits made without that winter's two years, but chose the architecture while seeing every winter. Here each of the four
laser winters is left out in turn: X-119's ``decide`` is applied to scores that never use its two years, and the
architecture it picks is scored on it with fits made without them. Pre-registration: the hub's
``insar_full_test/preregistration/X-135.yaml``.

X-119's code is imported, not copied — its inputs, folds (``run_fold``, ``spams_level``), scoring (``score_series``) and
rule (``decide``). The one part X-119 keeps inside ``main`` (the chain re-levelled by SPAMS across one winter) is
``a_winter`` below, and the stop rule checks that, with no extra exclusion, X-119's summary and verdict come back to
1e-9. Run on the cube X-119 ran on (the as-delivered snapshot of 2026-10-15, RZECIN_CUBE_ROOT).

    RZECIN_CUBE_ROOT=<hub>/06_data/cube_snapshots/as_delivered_2026-10-15 PYTHONPATH=src python experiments/nested_winter_x135.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import p6_validation_x082 as X  # noqa: E402
import spams_baseline_x112 as SB  # noqa: E402
import v4_head_to_head_x119 as H  # noqa: E402

from insar_wetlands import spams  # noqa: E402
from insar_wetlands.paths import cube_root, hub_root  # noqa: E402
from insar_wetlands.scoring import (  # noqa: E402
    Protocol,
    assert_preregistered,
    chain,
    guard_lockbox,
)

HUB = hub_root()
FT = HUB / "08_deliverables" / "insar_full_test"
OUT = FT / "nested_winter_x135"
PREREG = FT / "preregistration" / "X-135.yaml"
X119 = FT / "v4_head_to_head"
LOYO = FT / "v4" / "spams_loyo_params.yaml"
YEARS = H.YEARS
WINTERS = [f"{y}/{str(y + 1)[2:]}" for y in YEARS[:-1]]
ARCH = H.ORDER                                    # chain, A, B, C — simplest first, as X-119
FITTED = {"B: v3+": "B: v3+", "C: v3+ with SPAMS winter anchor": "C"}   # X-119's variant names for run_fold


def years_of(w: str) -> set[int]:
    return {int(w[:4]), int(w[:4]) + 1}


def pin_cube() -> None:
    """fusion v3's WP0, WP1, WP2, WP4 and WP5 scripts read ``<hub>/06_data/cube`` by a fixed path, not ``cube_root()``:
    point them at the cube this run reads (the snapshot X-119 ran on), as X-119's own run read it on 2026-10-15."""
    import fusion_v3_wp0_x073 as W0
    import fusion_v3_wp1_x074 as W1
    import fusion_v3_wp2_x075 as W2
    import fusion_v3_wp4_x077 as W4
    import fusion_v3_wp5_x078 as W5
    for m in (W0, W1, W2, W4, W5):
        m.G, m.S = cube_root() / "gold", cube_root() / "silver"


class Inputs:
    """Everything X-119's main builds before it scores, and a cache of its fits keyed by what they leave out."""

    def __init__(self):
        pin_cube()
        self.prot = Protocol.load(FT / "protocol.yaml")
        d = SB.table()
        self.d = d[d.year.between(2020, 2025)]
        import fusion_v3_wp5_x078 as W5
        _, self.p_all, self.noise = W5.inputs()
        p = H.pairs_corrected(self.p_all.copy())
        self.p = p[p.t2.dt.year <= 2025]
        guard_lockbox(pd.concat([self.p.t1, self.p.t2]), self.prot)
        lev = X.laser()[1]
        self.c0 = {}
        for track in X.TRACKS:
            g = self.p[self.p.track == track].sort_values("t1")
            c = chain(g.assign(est_mm=g.y_vert_mm), "est_mm")
            c["truth"] = lev.reindex(pd.DatetimeIndex(c.time).round("h")).to_numpy()
            self.c0[track] = c
        self._fold, self._spams = {}, {}

    def fold(self, variant: str, excl: frozenset):
        k = (variant, excl)
        if k not in self._fold:
            self._fold[k] = H.run_fold(self.d, self.p, self.noise, set(excl), self.prot, variant, self.p_all)
        return self._fold[k]

    def spams(self, excl: frozenset) -> pd.Series:
        if excl not in self._spams:
            self._spams[excl] = H.spams_level(self.d, set(excl))
        return self._spams[excl]


def a_winter(c0: pd.DataFrame, R: pd.Series, wname: str, prot: Protocol) -> pd.DataFrame:
    """X-119's A across one winter: the chain with that winter's step set to SPAMS's (``R``); as in X-119's main."""
    cA = c0.copy()
    z = cA.dropna(subset=["truth"])
    yy = int(wname[:4])
    a_ = z[(z.time.dt.year == yy) & z.time.dt.month.between(*prot.season_months)]
    b_ = z[(z.time.dt.year == yy + 1) & z.time.dt.month.between(*prot.spring_window)]
    if len(a_) and len(b_):
        ia, ib = a_.index[-1], b_.index[0]
        step = R.reindex(pd.DatetimeIndex([cA.time[ib]]).tz_convert(None).normalize()).iloc[0] - \
            R.reindex(pd.DatetimeIndex([cA.time[ia]]).tz_convert(None).normalize()).iloc[0]
        cA.loc[cA.index >= ib, "est"] += step - (cA.est[ib] - cA.est[ia])
        cA.loc[cA.index >= ib, "seg"] = cA.seg[ia]
    _, w = H.score_series(cA, prot, winter=wname)
    return w


def season_rows(I: Inputs, arch: str, track: str, year: int, excl: frozenset) -> list[dict]:
    c0 = I.c0[track]
    if arch in ARCH[:2]:                                  # chain and A are the chain within seasons: no fit
        s, _ = H.score_series(c0, I.prot, year=year)
    else:
        h = I.fold(FITTED[arch], excl).set_index("time").reindex(c0.time)
        c = pd.DataFrame({"time": c0.time, "seg": 0, "est": h.h.to_numpy(), "truth": c0.truth.to_numpy()})
        s, _ = H.score_series(c, I.prot, year=year)
    return [{"track": track, "estimate": arch, **r} for r in s]


def winter_rows(I: Inputs, arch: str, track: str, wname: str, excl: frozenset) -> list[dict]:
    c0 = I.c0[track]
    if arch == ARCH[0]:
        _, w = H.score_series(c0, I.prot, winter=wname)
    elif arch == ARCH[1]:
        w = a_winter(c0, I.spams(excl), wname, I.prot)
    else:
        h = I.fold(FITTED[arch], excl).set_index("time").reindex(c0.time)
        c = pd.DataFrame({"time": c0.time, "seg": 0, "est": h.h.to_numpy(), "truth": c0.truth.to_numpy()})
        _, w = H.score_series(c, I.prot, winter=wname)
    return [{"track": track, "estimate": arch, **r} for r in w.to_dict("records")]


def summarise(S: pd.DataFrame, W: pd.DataFrame) -> pd.DataFrame:
    """X-119's summary table (same groupings and statistics)."""
    return (S.groupby(["estimate", "track"]).rmse_detrended_mm.median().rename("season_rmse_median_mm").to_frame()
            .join(W.assign(ae=W.difference_mm.abs()).groupby(["estimate", "track"]).ae.mean().rename("winter_mae_mm"))
            .join(W.groupby(["estimate", "track"]).size().rename("winter_steps")).reset_index())


def scores(I: Inputs, outer: frozenset, skip: str | None, drops: list, tag: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Every architecture's rows on the seasons of years outside ``outer`` and on every winter but ``skip``, each fit
    leaving out the unit's years and ``outer`` (X-135 amendment: the neighbouring winters are kept)."""
    S, W = [], []
    for arch in ARCH:
        for track in X.TRACKS:
            for y in YEARS:
                if y in outer:
                    continue
                try:
                    S += season_rows(I, arch, track, y, frozenset({y}) | outer)
                except Exception as e:  # noqa: BLE001  (a fit with too few years left: dropped and reported)
                    drops.append({"fold": tag, "unit": f"season {y}", "estimate": arch, "track": track, "why": repr(e)[:200]})
            for w in WINTERS:
                if w == skip:
                    continue
                try:
                    W += winter_rows(I, arch, track, w, frozenset(years_of(w)) | outer)
                except Exception as e:  # noqa: BLE001
                    drops.append({"fold": tag, "unit": f"winter {w}", "estimate": arch, "track": track, "why": repr(e)[:200]})
    return pd.DataFrame(S).drop_duplicates(), pd.DataFrame(W)


def reproduce(I: Inputs) -> pd.DataFrame:
    """The stop rule: with nothing extra left out, X-119's summary and verdict come back."""
    S, W = scores(I, frozenset(), None, [], "reproduction")
    mine = summarise(S, W).set_index(["estimate", "track"]).sort_index()
    ref = pd.read_csv(X119 / "summary.csv").set_index(["estimate", "track"]).loc[mine.index]
    diff = (mine - ref[mine.columns]).abs().max()
    v_mine, v_ref = H.decide(mine.reset_index()), pd.read_csv(X119 / "verdict.csv")
    same_pick = v_mine.loc[v_mine.selected, "estimate"].iloc[0] == v_ref.loc[v_ref.selected, "estimate"].iloc[0]
    out = pd.DataFrame({"statistic": diff.index, "max_abs_diff": diff.to_numpy()})
    out.loc[len(out)] = {"statistic": "same architecture selected", "max_abs_diff": 0.0 if same_pick else 1.0}
    if (out.max_abs_diff > 1e-9).any():
        out.to_csv(OUT / "reproduction.csv", index=False)
        raise SystemExit(f"STOP: X-119 not reproduced (X-135 stop rule):\n{out.to_string(index=False)}")
    return out


def anchor_spread(I: Inputs, track: str, wname: str) -> tuple[float, float]:
    """A's step across ``wname`` with each of the five leave-one-year-out SPAMS fits (X-107 amendment_3): min, max."""
    fits = yaml.safe_load(LOYO.read_text())["fits"]
    P, E = I.d.rain_era5_mm.to_numpy(), I.d.et0_mm.to_numpy()
    steps = []
    for f in fits.values():
        R = pd.Series(spams.reversible(P, E, f["x_p"], f["x_e"], int(f["tau_days"])), index=I.d.index)
        w = a_winter(I.c0[track], R, wname, I.prot)
        if len(w):
            steps.append(float(w.insar_change_mm.iloc[0]))
    return (min(steps), max(steps)) if steps else (np.nan, np.nan)


def main() -> None:
    assert_preregistered(PREREG, [OUT / "verdict.csv"])
    OUT.mkdir(parents=True, exist_ok=True)
    I = Inputs()
    rep = reproduce(I)
    rep.to_csv(OUT / "reproduction.csv", index=False)
    drops, choices, nested = [], [], []
    for w in WINTERS:
        outer = frozenset(years_of(w))
        S, W = scores(I, outer, w, drops, w)
        summ = summarise(S, W)
        if W.empty:                                   # no inner winter left: X-119's simplicity order (the chain)
            pick, verdict = ARCH[0], pd.DataFrame()
        else:
            verdict = H.decide(summ)
            pick = verdict.loc[verdict.selected, "estimate"].iloc[0]
        choices.append(summ.assign(outer_winter=w, chosen=pick).merge(
            verdict[["estimate", "eligible", "winter_mae_mean_mm"]] if len(verdict) else pd.DataFrame(columns=["estimate"]),
            on="estimate", how="left"))
        for track in X.TRACKS:
            for arch, role in ((pick, "nested choice"), (ARCH[0], "chain"), (ARCH[1], "X-119's A")):
                rows = winter_rows(I, arch, track, w, outer)
                for r in rows:
                    lo, hi = anchor_spread(I, track, w) if arch == ARCH[1] else (np.nan, np.nan)
                    nested.append({"outer_winter": w, "role": role, "anchor_step_min_mm": lo, "anchor_step_max_mm": hi, **r})
    C, N, D = pd.concat(choices, ignore_index=True), pd.DataFrame(nested), pd.DataFrame(drops)
    N["abs_error_mm"] = N.difference_mm.abs()
    per_track = N.groupby(["role", "track"]).agg(winter_mae_mm=("abs_error_mm", "mean"), winters=("abs_error_mm", "size")).reset_index()
    both = per_track.groupby("role").winter_mae_mm.mean()
    picks = C.drop_duplicates("outer_winter").set_index("outer_winter").chosen
    nested_mae, chain_mae = float(both["nested choice"]), float(both["chain"])
    if nested_mae <= 0.90 * chain_mae:
        reading = "carries across unseen winters"
    elif nested_mae < chain_mae:
        reading = "carries in part"
    else:
        reading = "does not carry"
    verdict = pd.DataFrame([{"nested_winter_mae_mm": nested_mae, "chain_winter_mae_mm": chain_mae,
                             "x119_A_winter_mae_mm": float(both["X-119's A"]), "ratio_to_chain": nested_mae / chain_mae,
                             "choice_stable": picks.nunique() == 1, "choices": "; ".join(f"{k}: {v}" for k, v in picks.items()),
                             "reading": reading, "laser_rule": "frost rule of 2026-10-15 (as the frozen v4 and X-107)"}])
    for name, df in (("choices", C), ("nested_steps", N), ("per_track", per_track), ("drops", D), ("verdict", verdict)):
        df.to_csv(OUT / f"{name}.csv", index=False)
    (OUT / "README.md").write_text(readme(verdict, per_track, N, C, D, rep), encoding="utf-8")
    print(verdict.T.to_string())


def readme(v, per_track, N, C, D, rep) -> str:
    r = v.iloc[0]
    L = ["# X-135 — the v4 procedure scored by nested leave-one-winter-out (generated)", "",
         "Generated by `experiments/nested_winter_x135.py` (`make run E=nested_winter_x135`); do not edit. Pre-registration:",
         "`../preregistration/X-135.yaml`. Laser years only; the laser on the frost rule of 2026-10-15, as the frozen v4 and",
         "X-107; the cube X-119 ran on (as-delivered snapshot 2026-10-15). For each laser winter: X-119's choice made",
         "without its two years, every fit without them, then the chosen architecture scored on it.", "",
         "## Verdict", "",
         f"**{r.reading}** — nested winter-step error {r.nested_winter_mae_mm:.1f} mm against the chain's"
         f" {r.chain_winter_mae_mm:.1f} mm (ratio {r.ratio_to_chain:.2f}; the registered margin is 0.90), averaged over the"
         f" two tracks and the four winters. X-119's A, whose choice saw every winter: {r.x119_A_winter_mae_mm:.1f} mm.",
         f"Choice stable across the folds: {'yes' if r.choice_stable else 'no'} ({r.choices}).", "",
         "## Per track", "", "| role | track | winter_mae_mm | winters |", "|---|---|---|---|",
         *[f"| {x.role} | {x.track} | {x.winter_mae_mm:.1f} | {x.winters} |" for x in per_track.itertuples()], "",
         "## Winter by winter (vertical mm)", "",
         "| outer winter | track | role | estimate | laser step | estimated step | error | anchor spread (A) |",
         "|---|---|---|---|---|---|---|---|",
         *[f"| {x.outer_winter} | {x.track} | {x.role} | {x.estimate} | {x.laser_change_mm:.0f} | {x.insar_change_mm:.0f} |"
           f" {x.difference_mm:+.0f} | {'' if np.isnan(x.anchor_step_min_mm) else f'{x.anchor_step_min_mm:.0f} to {x.anchor_step_max_mm:.0f}'} |"
           for x in N.sort_values(["outer_winter", "track", "role"]).itertuples()], "",
         "## Stop rule", "",
         f"With nothing extra left out, X-119's summary is reproduced: largest difference {rep.max_abs_diff.max():.1g};"
         " same architecture selected.", "",
         f"Inner units dropped for lack of data: {len(D)}" + (" (`drops.csv`)." if len(D) else "."), "",
         "## Files", "", "| File | Content |", "|---|---|",
         "| `verdict.csv` | the registered reading |", "| `per_track.csv` | nested, chain and X-119's A, per track |",
         "| `nested_steps.csv` | every scored winter step |", "| `choices.csv` | each fold's inner summary and choice |",
         "| `reproduction.csv` | the stop rule's check against X-119 |", "| `drops.csv` | inner units dropped |"]
    return "\n".join(L) + "\n"


if __name__ == "__main__":
    main()

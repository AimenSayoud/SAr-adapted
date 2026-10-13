"""X-095 — phase simulated from the P6 laser: what the cycle choice and segment-and-tie can recover in principle.

conroy2022's benchmark design at Rzecin: the laser increment of every real consecutive 40 m pair at the P6 cell (both
tracks, design years, D-028) is given realistic noise, a moisture phase and a short-pair bias, then wrapped
(``insar_wetlands.phase_sim``). Processing choices then try to recover it:

- cycle: C0 the wrapped value; CP the cycle closest to the water-table prior g·ΔWTD (buoyancy, g fitted on the laser
  leave-one-year-out, T47); CO the true cycle (oracle);
- integration: I0 through every pair; I1 cut at pairs flagged frozen or snow and tie each segment to the prior level
  (X-093's approach); IO tie to the truth (oracle).

Scored by the common scorer's rules (C-061, ``protocol.yaml``): the primary — the across-winter level error — over 500
realisations per scenario; level views and increments by the full scorer on 20. Pre-registration:
``08_deliverables/insar_full_test/preregistration/X-095.yaml``. Reads the cube and the hub; writes only to the hub:
``08_deliverables/insar_full_test/phase_sim/``.

    make run E=phase_sim_x095
"""
from __future__ import annotations

import sys
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import p6_validation_x082 as X  # noqa: E402  (pairs filter, laser in LOS: unchanged)

from insar_wetlands import phase_sim as ps  # noqa: E402
from insar_wetlands.paths import hub_root  # noqa: E402
from insar_wetlands.scoring import SEASON_OF_MONTH, Protocol, assert_preregistered, chain, score, winter_steps  # noqa: E402

HUB = hub_root()
FT = HUB / "08_deliverables" / "insar_full_test"
OUT = FT / "phase_sim"
PREREG = FT / "preregistration" / "X-095.yaml"
R, R_FULL = 500, 20
BIASES = (0.1, 0.3, 1.0)
GOLD = ["track", "grid", "plot", "pair", "t1", "t2", "revisit_days", "d_laser_los_mm", "d_laser_vertical_mm",
        "d_wtd_P6_cm", "flag_frozen_air_t1", "flag_frozen_air_t2", "flag_snow_imgw_t1", "flag_snow_imgw_t2",
        "unw_1x1", "ref_hard_unw"]
C = ps.cycle_mm()


def truthy(s: pd.Series) -> np.ndarray:
    return s.map({True: True, "True": True}).fillna(False).to_numpy(bool)


def robust_sd(v) -> float:
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    return float(1.4826 * np.median(np.abs(v - np.median(v)))) if len(v) else np.nan


def load() -> tuple[pd.DataFrame, dict[str, pd.Series], pd.Series, pd.Series]:
    p = pd.read_csv(X.CUBE / "gold" / "plots_pair.csv.gz", usecols=GOLD, low_memory=False)
    p = p[(p["plot"] == "P6") & (p.grid == "40m")].copy()
    for c in ("t1", "t2"):
        p[c] = pd.to_datetime(p[c], format="ISO8601", utc=True)
    p = p[X.consecutive(p)].sort_values(["track", "t1"]).reset_index(drop=True)
    p["obs_mm"] = (p.unw_1x1 - p.ref_hard_unw) * X.PHASE_TO_MM
    p["cut"] = truthy(p.flag_frozen_air_t1) | truthy(p.flag_frozen_air_t2) | truthy(p.flag_snow_imgw_t1) | \
        truthy(p.flag_snow_imgw_t2)
    las, lev = X.laser()
    site = pd.read_csv(X.CUBE / "silver" / "site_hourly.csv.gz", index_col=0, usecols=["time_utc", "wtd_P6_cm"])
    site.index = pd.to_datetime(site.index, format="ISO8601", utc=True)
    return p, las, lev, site.wtd_P6_cm


def buoyancy_loyo(p: pd.DataFrame) -> dict[int, float]:
    """g (mm of vertical surface per cm of water table), through the origin, fitted on every year but the one used."""
    d = p.dropna(subset=["d_laser_vertical_mm", "d_wtd_P6_cm"]).drop_duplicates("pair")
    yrs = d.t1.dt.year.to_numpy()
    out = {}
    for y in range(2020, 2027):
        m = yrs != y
        x, z = d.d_wtd_P6_cm.to_numpy()[m], d.d_laser_vertical_mm.to_numpy()[m]
        out[y] = float(np.sum(x * z) / np.sum(x * x)) if len(x) else np.nan
    return out


def track_inputs(p: pd.DataFrame, track: str, las: pd.Series, lev: pd.Series, wtd: pd.Series,
                 g: dict[int, float]) -> dict:
    lt = las.dropna()
    q = p[(p.track == track) & (p.t1 >= lt.index.min()) & (p.t2 <= lt.index.max()) &
          ~((p.t1 < X.RELEVEL[1]) & (p.t2 > X.RELEVEL[0]))].sort_values("t1").reset_index(drop=True)
    # truth increment from the laser level interpolated in time (pre-registration amendment): winter pairs have no
    # valid laser hour; the winter steps are still scored on real laser values at both ends
    tx = lt.index.asi8.astype(float)
    lvl = lambda t: np.interp(pd.DatetimeIndex(t).asi8.astype(float), tx, lt.to_numpy())  # noqa: E731
    q["truth_inc"] = lvl(q.t2.dt.round("h")) - lvl(q.t1.dt.round("h"))     # the hour the scorer reads
    q["truth_interpolated"] = q.d_laser_los_mm.isna()
    cosi = float((las / lev).median())                         # LOS / vertical at P6 for this track
    yr = q.t1.dt.year.to_numpy()
    gy = np.array([g.get(int(y), np.nan) for y in yr])
    q["prior_inc"] = gy * q.d_wtd_P6_cm.to_numpy() * cosi
    # residual model of the observed phase: moisture term e·ΔWTD, then noise by season (pre-registration S1, S2)
    res = (q.obs_mm - q.d_laser_los_mm).to_numpy()
    dw = q.d_wtd_P6_cm.to_numpy()
    m = np.isfinite(res) & np.isfinite(dw) & (np.abs(res) < C)
    e = float(np.sum(dw[m] * res[m]) / np.sum(dw[m] ** 2)) if m.sum() > 5 else 0.0
    res_clean = res - e * np.nan_to_num(dw)
    q["season"] = q.t1.dt.month.map(SEASON_OF_MONTH)
    sig = {s: robust_sd(res_clean[(q.season == s).to_numpy()]) for s in q.season.unique()}
    q["sigma"] = q.season.map(sig).fillna(robust_sd(res_clean))
    q["moist"] = e * np.nan_to_num(dw)
    epochs = pd.DatetimeIndex(sorted(set(q.t1) | set(q.t2)))
    yr_e = epochs.year.to_numpy()
    w_e = wtd.reindex(epochs.round("h")).to_numpy(float)
    prior_level = pd.Series(np.array([g.get(int(y), np.nan) for y in yr_e]) * w_e * cosi, index=epochs)
    have = q.d_laser_los_mm.notna()
    check = float(np.max(np.abs(q.truth_inc[have] - q.d_laser_los_mm[have]))) if have.any() else np.nan
    return {"q": q, "truth_level": lt, "prior_level": prior_level, "e": e, "sigma": sig, "cos": cosi,
            "max_diff_interp_vs_pair_laser_mm": check}


def blocks(q: pd.DataFrame) -> list[np.ndarray]:
    """Index runs of contiguous pairs (t1 equals the previous t2)."""
    out, cur = [], [0]
    for i in range(1, len(q)):
        if q.t1.iloc[i] == q.t2.iloc[i - 1]:
            cur.append(i)
        else:
            out.append(np.array(cur))
            cur = [i]
    out.append(np.array(cur))
    return out


def integrate(q: pd.DataFrame, inc: np.ndarray, how: str, prior_level: pd.Series, truth_level: pd.Series) -> np.ndarray:
    if how == "I0":
        return inc
    ref = prior_level if how == "I1" else pd.Series(truth_level.reindex(
        pd.DatetimeIndex(sorted(set(q.t1) | set(q.t2))).round("h")).to_numpy(),
        index=pd.DatetimeIndex(sorted(set(q.t1) | set(q.t2))))
    out = inc.copy()
    keep = ~q.cut.to_numpy()
    for b in blocks(q):
        out[b] = ps.tie_segments(q.t1.iloc[b], q.t2.iloc[b], inc[b], keep[b], ref)
    return out


def winter_design(q: pd.DataFrame, truth_level: pd.Series, prot: Protocol):
    """The scorer's chain and winter steps for this pair set: a matrix A with chain levels = A @ increments (the
    scorer's ``chain``: a gap starts a new segment at zero), and the (winter, row from, row to, truth change) of each
    winter step that lies in one segment."""
    c = chain(q.assign(est_mm=0.0))
    c["truth"] = truth_level.reindex(pd.DatetimeIndex(c.time).round("h")).to_numpy()
    A = np.zeros((len(c), len(q)))
    row, members, prev = -1, [], None
    for i in range(len(q)):                      # same walk as scoring.chain (q is sorted by t1)
        if prev is None or q.t1.iloc[i] != prev:
            row += 1
            members = []
        members.append(i)
        row += 1
        A[row, members] = 1.0
        prev = q.t2.iloc[i]
    assert row == len(c) - 1
    z = c.dropna(subset=["truth"]).copy()
    w = winter_steps(z.assign(est=0.0), prot)
    rows = []
    for _, s in w.iterrows():
        if not s.same_segment:
            continue
        d = pd.DatetimeIndex(z.time).date
        rows.append((s.winter, int(z.index[d == s["from"]][-1]), int(z.index[d == s["to"]][0]), s.laser_change_mm))
    return A, rows


def winter_errors(A: np.ndarray, rows, E: np.ndarray) -> np.ndarray:
    """(R, n_winters) across-winter level error, estimate minus truth."""
    L = E @ A.T
    return np.column_stack([(L[:, jb] - L[:, ja]) - dl for _, ja, jb, dl in rows]) if rows else np.empty((len(E), 0))


def main(argv: list[str] | None = None) -> None:
    global R, R_FULL, OUT
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--realisations", type=int, default=R, help="smoke runs only; the pre-registered run uses 500")
    ap.add_argument("--out", default=None, help="smoke runs only: write elsewhere than the deliverable")
    a = ap.parse_args(argv)
    if a.realisations != R:
        R, R_FULL = a.realisations, min(R_FULL, a.realisations)
    if a.out:
        OUT = Path(a.out)
    else:
        assert_preregistered(PREREG, [OUT / "primary.csv"])
    OUT.mkdir(parents=True, exist_ok=True)
    prot = Protocol.load(FT / "protocol.yaml")
    p, las, lev, wtd = load()
    g = buoyancy_loyo(p)
    scen = {"S0": (0, 0, 0.0), "S1": (1, 0, 0.0), "S2": (1, 1, 0.0), **{f"S3 b={b}": (1, 1, b) for b in BIASES}}
    prim, cyc, lev_rows, params = [], [], [], []
    for track in X.TRACKS:
        ti = track_inputs(p, track, las[track], lev, wtd, g)
        q = ti["q"]
        A, wrows = winter_design(q, ti["truth_level"], prot)
        probe = np.random.default_rng(7).normal(0, 5, len(q))      # A reproduces the scorer's chain exactly
        assert np.allclose(A @ probe, chain(q.assign(est_mm=probe)).est.to_numpy())
        params.append({"track": track, "pairs": len(q), "pairs_truth_interpolated": int(q.truth_interpolated.sum()),
                       "max_diff_interp_vs_pair_laser_mm": ti["max_diff_interp_vs_pair_laser_mm"],
                       "pairs_cut": int(q.cut.sum()), "winters_scored": len(wrows),
                       "moisture_e_mm_per_cm": ti["e"], "los_over_vertical": ti["cos"],
                       **{f"sigma_{k}_mm": v for k, v in ti["sigma"].items()},
                       **{f"g_without_{y}": v for y, v in g.items() if 2021 <= y <= 2025}})
        truth = q.truth_inc.to_numpy()
        for sname, (use_noise, use_moist, bias) in scen.items():
            rng = np.random.default_rng(zlib.crc32(f"{track}|{sname}".encode()))
            n_r = 1 if sname == "S0" else R
            extra = (q.moist.to_numpy() if use_moist else 0.0) + bias
            W, F = [], []
            for _ in range(n_r):
                w, f = ps.simulate(truth, q.sigma.to_numpy() if use_noise else 0.0, extra, rng if use_noise else None)
                W.append(w)
                F.append(f)
            W, F = np.array(W), np.array(F)
            refs = {"C0": None, "CP": q.prior_inc.to_numpy(), "CO": F}
            for cname, ref in refs.items():
                U = np.array([ps.choose_cycle(W[i], ref[i] if cname == "CO" else ref) for i in range(n_r)])
                right = np.isclose(U, F, atol=1e-6)
                for (L, s_), m in q.groupby(["revisit_days", "season"]).groups.items():
                    cyc.append({"track": track, "scenario": sname, "cycle": cname, "revisit_days": L, "season": s_,
                                "pairs": len(m), "right_cycle_share": float(right[:, m].mean())})
                for iname in ("I0", "I1", "IO"):
                    E = np.array([integrate(q, U[i], iname, ti["prior_level"], ti["truth_level"]) for i in range(n_r)])
                    we = winter_errors(A, wrows, E)
                    med = np.median(np.abs(we), axis=1) if we.size else np.full(n_r, np.nan)
                    prim.append({"track": track, "scenario": sname, "cycle": cname, "integration": iname,
                                 "realisations": n_r, "winters": we.shape[1],
                                 "median_abs_winter_error_mm": float(np.median(med)),
                                 "p05": float(np.quantile(med, 0.05)), "p95": float(np.quantile(med, 0.95))})
                    for i in range(min(n_r, R_FULL)):
                        sc = score(q.assign(est_mm=E[i], truth_inc_mm=truth), ti["truth_level"], prot)
                        for _, r in sc["levels"].iterrows():
                            if r.view in ("full span", "datum A", "datum B"):
                                lev_rows.append({"track": track, "scenario": sname, "cycle": cname,
                                                 "integration": iname, "view": r.view,
                                                 "rmse_detrended_mm": r.rmse_detrended_mm,
                                                 "pearson_detrended": r.pearson_detrended})
    prim, cyc, lev, params = map(pd.DataFrame, (prim, cyc, lev_rows, params))
    lev = lev.groupby(["track", "scenario", "cycle", "integration", "view"]).median().reset_index()
    gate = gate_table(prim)
    for name, df in (("primary", prim), ("cycles", cyc), ("levels", lev), ("parameters", params), ("gate", gate)):
        df.to_csv(OUT / f"{name}.csv", index=False)
    (OUT / "README.md").write_text(readme(prim, cyc, lev, params, gate), encoding="utf-8")
    print(f"wrote {OUT}")


def gate_table(prim: pd.DataFrame) -> pd.DataFrame:
    """The pre-registered gate: under S2 with the prior cycle, does I1 cut the across-winter error by >= 30 %?"""
    rows = []
    for track, g in prim[(prim.scenario == "S2") & (prim.cycle == "CP")].groupby("track"):
        i0 = g[g.integration == "I0"].median_abs_winter_error_mm.iloc[0]
        i1 = g[g.integration == "I1"].median_abs_winter_error_mm.iloc[0]
        rows.append({"track": track, "I0_mm": i0, "I1_mm": i1, "reduction": 1 - i1 / i0 if i0 else np.nan,
                     "passes": bool(i0 and (1 - i1 / i0) >= 0.30)})
    return pd.DataFrame(rows)


def _md(df: pd.DataFrame, nd: int = 2) -> str:
    d = df.copy()
    for c in d.columns:
        if pd.api.types.is_float_dtype(d[c]):
            d[c] = d[c].map(lambda x: "—" if pd.isna(x) else f"{x:.{nd}f}")
    return "\n".join(["| " + " | ".join(map(str, d.columns)) + " |", "|" + "---|" * len(d.columns)] +
                     ["| " + " | ".join(map(str, r)) + " |" for r in d.itertuples(index=False)])


def readme(prim, cyc, lev, params, gate) -> str:
    pv = prim.pivot_table(index=["track", "scenario", "cycle"], columns="integration",
                          values="median_abs_winter_error_mm").reset_index()
    cy = cyc[cyc.scenario.isin(["S1", "S2"])].groupby(["track", "scenario", "cycle", "revisit_days"]).apply(
        lambda d: np.average(d.right_cycle_share, weights=d.pairs), include_groups=False).rename(
        "right_cycle_share").reset_index()
    lv = lev[(lev.view == "datum B") & lev.scenario.isin(["S1", "S2"])].pivot_table(
        index=["track", "scenario", "cycle"], columns="integration", values="rmse_detrended_mm").reset_index()
    pt = params.T.reset_index().rename(columns={"index": "parameter"})
    return f"""# X-095 — phase simulated from the P6 laser (generated)

Generated by `experiments/phase_sim_x095.py` (`make run E=phase_sim_x095`); do not edit. Pre-registration:
`../preregistration/X-095.yaml`. Simulator: `insar_wetlands.phase_sim` (tested: no noise returns the laser exactly).
Truth: the P6 laser in LOS on the real consecutive 40 m pairs, design years 2021–2025 (D-028). **What is recoverable
in principle on this sampling and this prior — not a result on the real phase.** Tried on the laser in hand only.

Scenarios: S0 no noise; S1 empirical noise by season; S2 S1 + moisture term e·ΔWTD; S3 S2 + bias b mm per pair.
Cycle: C0 wrapped value; CP closest to the water-table prior; CO oracle. Integration: I0 through every pair; I1 cut at
frozen/snow pairs and tie to the prior level; IO tie to the truth.

## Gate for X-093(a) (pre-registered: I1 cuts the across-winter error by >= 30 % under S2 with CP, both tracks)

{_md(gate)}

## Across-winter level error: median |estimate − laser| over winters, median over realisations (mm)

{_md(pv)}

## Share of pairs given the right cycle (S1, S2; pair-weighted over seasons)

{_md(cy)}

## Level RMSE, datum B, detrended (mm; median over {R_FULL} realisations)

{_md(lv)}

## Parameters taken from the observed data

{_md(pt, 3)}

Files: `primary.csv` (every scenario × cycle × integration, with the 5–95 % range), `cycles.csv` (by pair length and
season), `levels.csv` (scorer level views), `parameters.csv`, `gate.csv`.
"""


if __name__ == "__main__":
    main()

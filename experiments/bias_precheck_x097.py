"""X-097 — pre-check for X-085: does the mat meet what the C-band short-pair bias corrections need?

The corrections (maghsoudi2022, zheng2022, ma2024, ansari2020) each state the conditions they need: coherent
connections of given lengths, a closure-phase amplitude they can work with, closures that unwrap consistently, coherent
long pairs, unfiltered input. This measures those conditions on the HyP3 40 m pairs of the cube, on the mat, at the P6
cell, on grassland and on hard targets, both tracks, 2017–2026, and gives a verdict per correction by the rule fixed in
the pre-registration (``08_deliverables/insar_full_test/preregistration/X-097.yaml``).

Measures per pixel: M1 mean coherence of the pairs of one connection level; M2 the closure amplitude ``|tau|`` and
phase ``arg(tau)``, ``tau`` the time mean of ``exp(i Phi)`` over the sequential closures of level n (n consecutive base
pairs minus the n-span pair, from wrapped phase); M3 the share of level-2 triplets whose unwrapped closure differs from
the wrapped one by a whole cycle (pairs referenced to hard targets); M4 mean coherence of pairs of 48 days or more.

Open water (WorldCover water share >= 0.9) is carried as a reference class, added after the first run and not part
of the pre-registered rule: it holds no phase, so where its HyP3 coherence or closure amplitude equals the mat's, a
condition measured on the mat cannot be told from the filtered product's floor.

No laser is read. Reads the cube's silver interferograms and static layers; writes only to the hub:
``08_deliverables/insar_full_test/bias_precheck/``.

    make run E=bias_precheck_x097
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import xarray as xr

from insar_wetlands.cube import cube_root
from insar_wetlands.paths import hub_root
from insar_wetlands.scoring import assert_preregistered

HUB = hub_root()
CUBE = cube_root()
FT = HUB / "08_deliverables" / "insar_full_test"
OUT = FT / "bias_precheck"
PREREG = FT / "preregistration" / "X-097.yaml"
SETS = {"2017_2019": 6, "2020_2021": 6, "2022_2024": 12, "2025": 6, "2026": 6}
TRACKS = ("ascending", "descending")
WINTER = {12, 1, 2, 3, 4}
COH_MIN, TAU_AMP_MIN, TAU_ARG_MIN = 0.3, 0.3, 0.3
UNW_ZHENG_MAX, UNW_MA_MAX = 0.10, 0.5
EPOCH_SHARE_MIN, CLASS_SHARE_MIN = 0.8, 0.5
LONG_DAYS = 48
WATER_SHARE = 0.9          # open-water reference class (post hoc, descriptive)
MIN_CLOSURES = 5


def season(ts) -> np.ndarray:
    return np.where(pd.DatetimeIndex(ts).month.isin(list(WINTER)), "winter", "snow-free")


def load(name: str, track: str) -> xr.Dataset | None:
    f = CUBE / "silver" / f"ifg_{name}_{track}.nc"
    return xr.open_dataset(f) if f.exists() else None


def index(d: xr.Dataset) -> dict[tuple[pd.Timestamp, pd.Timestamp], int]:
    t1, t2 = pd.to_datetime(d.t1.values), pd.to_datetime(d.t2.values)
    return {(a, b): i for i, (a, b) in enumerate(zip(t1, t2))}


def connection_mean(d: xr.Dataset, days: int, var: str = "coh") -> dict[str, np.ndarray]:
    """Mean of ``var`` over the pairs of one length, per season and over all: {season: (y, x)}."""
    dd = (pd.to_datetime(d.t2.values) - pd.to_datetime(d.t1.values)).days
    sel = np.flatnonzero(dd == days)
    if not len(sel):
        return {}
    v = d[var].isel(pair=sel).values.astype("float32")
    s = season(d.t2.values[sel])
    out = {"all": np.nanmean(v, axis=0)}
    for k in ("winter", "snow-free"):
        if (s == k).any():
            out[k] = np.nanmean(v[s == k], axis=0)
    return out


def epochs_with_levels(d: xr.Dataset, step: int, levels=(1, 2, 3)) -> float:
    """Share of epochs from which pairs of every listed level exist (maghsoudi2022's network)."""
    ix = index(d)
    starts = sorted({a for a, _ in ix})
    if not starts:
        return float("nan")
    ok = [all((a, a + pd.Timedelta(days=step * n)) in ix for n in levels) for a in starts]
    return float(np.mean(ok))


def closures(d: xr.Dataset, step: int, n: int) -> tuple[np.ndarray, np.ndarray, list[tuple[int, ...]]]:
    """Sequential closures of level n from wrapped phase: (k, y, x) closure phase, last dates, the pair indices."""
    ix = index(d)
    w = d.wrapped.values.astype("float32")
    phis, ends, used = [], [], []
    for a in sorted({a for a, _ in ix}):
        chain = [(a + pd.Timedelta(days=step * j), a + pd.Timedelta(days=step * (j + 1))) for j in range(n)]
        span = (a, a + pd.Timedelta(days=step * n))
        if all(c in ix for c in chain) and span in ix:
            ids = [ix[c] for c in chain]
            z = np.exp(1j * (w[ids].sum(axis=0) - w[ix[span]]))
            phis.append(np.angle(z).astype("float32"))
            ends.append(span[1])
            used.append((*ids, ix[span]))
    if not phis:
        return np.empty((0,) + w.shape[1:], "float32"), np.array([], "datetime64[ns]"), []
    return np.stack(phis), np.array(ends, "datetime64[ns]"), used


def tau(phi: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    z = np.nanmean(np.exp(1j * phi), axis=0)
    return np.abs(z), np.angle(z)


def unwrap_inconsistency(d: xr.Dataset, used: list[tuple[int, ...]], phi_w: np.ndarray, hard: np.ndarray,
                         ) -> np.ndarray:
    """(k, y, x) boolean: the level-2 triplet's unwrapped closure departs from its wrapped closure by >= 1 cycle.

    Each pair's unwrapped phase is referenced to its median over hard targets; the remaining per-triplet constant is
    removed as the hard-target median of the difference, so a cycle is counted relative to the hard targets."""
    unw = d.unw.values.astype("float32")
    cc = d.conncomp.values
    unw = np.where(cc > 0, unw, np.nan)
    out = np.zeros(phi_w.shape, bool)
    valid = np.zeros(phi_w.shape, bool)
    for k, (i12, i23, i13) in enumerate(used):
        p = [unw[i] - np.nanmedian(unw[i][hard]) for i in (i12, i23, i13)]
        diff = p[0] + p[1] - p[2] - phi_w[k]
        diff = diff - np.nanmedian(diff[hard])
        cyc = np.rint(diff / (2 * np.pi))
        valid[k] = np.isfinite(cyc)
        out[k] = np.abs(np.nan_to_num(cyc)) >= 1
    return np.where(valid, out, np.nan)


def class_value(m: np.ndarray, masks: dict[str, np.ndarray], p6: tuple[int, int], holds) -> dict[str, tuple[float, float]]:
    """Per class: (median of the measure, share of pixels where it holds); at P6 the cell's value and 0/1."""
    out = {}
    for k, msk in masks.items():
        v = m[msk]
        v = v[np.isfinite(v)]
        out[k] = (float(np.median(v)), float(np.mean(holds(v)))) if len(v) else (np.nan, np.nan)
    v = m[p6]
    out["P6"] = (float(v), float(holds(np.array([v]))[0])) if np.isfinite(v) else (np.nan, np.nan)
    return out


def main() -> None:
    assert_preregistered(PREREG, [OUT / "verdicts.csv"])
    OUT.mkdir(parents=True, exist_ok=True)
    st = xr.open_dataset(CUBE / "silver" / "static_40m.nc")
    zone, hard = st.zone.values, st.hard_target.values.astype(bool)
    masks = {"mat": zone == 1, "grassland": zone == 3, "hard_targets": hard,
             "open_water": st.wc_water.values >= WATER_SHARE}
    px = pd.read_csv(HUB / "08_deliverables" / "field_first" / "plot_pixels.csv").set_index("plot")
    p6 = (int(px.loc["P6", "row"]), int(px.loc["P6", "col"]))
    cond, conn, maps, unw_n = [], [], {}, []

    def add(corr, c, track, setname, seas, m, thr, holds, measure):
        for cls, (val, share) in class_value(m, masks, p6, holds).items():
            ok = (share >= 1.0) if cls == "P6" else (share >= CLASS_SHARE_MIN)
            cond.append({"correction": corr, "condition": c, "measure": measure, "track": track, "set": setname,
                         "season": seas, "class": cls, "median": val, "share_holds": share, "threshold": thr,
                         "holds": bool(ok) if np.isfinite(share) else None})

    for track in TRACKS:
        long_parts = []
        for name, step in SETS.items():
            d = load(name, track)
            if d is None:
                continue
            dd = (pd.to_datetime(d.t2.values) - pd.to_datetime(d.t1.values)).days
            # M1 coherence by connection level (all levels present, for the table)
            for n in range(1, 5):
                for seas, m in connection_mean(d, step * n).items():
                    for cls, (val, _) in class_value(m, masks, p6, lambda v: v > COH_MIN).items():
                        conn.append({"track": track, "set": name, "days": step * n, "level": n, "season": seas,
                                     "class": cls, "median_coherence": val})
            # maghsoudi2022 (a) network, (b) level-3 coherence
            share_ep = epochs_with_levels(d, step)
            for cls in (*masks, "P6"):
                cond.append({"correction": "maghsoudi2022", "condition": "a: connections 1-3 from >= 80 % of epochs",
                             "measure": "share of epochs", "track": track, "set": name, "season": "all", "class": cls,
                             "median": share_ep, "share_holds": share_ep, "threshold": EPOCH_SHARE_MIN,
                             "holds": bool(share_ep >= EPOCH_SHARE_MIN)})
            for seas, m in connection_mean(d, step * 3).items():
                add("maghsoudi2022", f"b: coherence at level 3 ({step * 3} d) > 0.3", track, name, seas, m, COH_MIN,
                    lambda v: v > COH_MIN, "M1")
            # ma2024 (a) level-1 coherence not low
            for seas, m in connection_mean(d, step).items():
                add("ma2024", f"a: coherence at level 1 ({step} d) > 0.3", track, name, seas, m, COH_MIN,
                    lambda v: v > COH_MIN, "M1")
            # closures: level 2, 3, 4 wrapped (zheng a, b at the highest level with enough closures); level-2 unwrapping
            top = None
            for n in (2, 3, 4):
                phi, ends, used = closures(d, step, n)
                if len(phi) < MIN_CLOSURES:
                    continue
                top = (n, phi, ends)
                amp, arg = tau(phi)
                maps[f"tau_amp_{track}_{name}_n{n}"] = amp
                maps[f"tau_arg_{track}_{name}_n{n}"] = arg
                if n == 2:
                    bad = unwrap_inconsistency(d, used, phi, hard)
                    unw_n.append({"track": track, "set": name, "triplets": len(used),
                                  "with_hard_target_reference": int(np.isfinite(bad[:, hard]).any(axis=1).sum())})
                    s = season(ends)
                    for seas, sel in (("all", np.ones(len(s), bool)), ("winter", s == "winter"),
                                      ("snow-free", s == "snow-free")):
                        if sel.sum() < MIN_CLOSURES:
                            continue
                        share = np.nanmean(bad[sel], axis=0)
                        maps.setdefault(f"unw_inconsistent_{track}_{name}", np.nanmean(bad, axis=0))
                        add("zheng2022", "c: unwrapped level-2 closures consistent (share <= 0.10)", track, name, seas,
                            share, UNW_ZHENG_MAX, lambda v: v <= UNW_ZHENG_MAX, "M3")
                        add("ma2024", "b: unwrapping errors not in nearly every interferogram (share < 0.5)", track,
                            name, seas, share, UNW_MA_MAX, lambda v: v < UNW_MA_MAX, "M3")
            if top is not None:
                n, phi, ends = top
                s = season(ends)
                for seas, sel in (("all", np.ones(len(s), bool)), ("winter", s == "winter"),
                                  ("snow-free", s == "snow-free")):
                    if sel.sum() < MIN_CLOSURES:
                        continue
                    amp, arg = tau(phi[sel])
                    add("zheng2022", f"a: |tau| >= 0.3 (level {n})", track, name, seas, amp, TAU_AMP_MIN,
                        lambda v: v >= TAU_AMP_MIN, "M2")
                    add("zheng2022", f"b: |arg tau| >= 0.3 rad (level {n})", track, name, seas, np.abs(arg),
                        TAU_ARG_MIN, lambda v: v >= TAU_ARG_MIN, "M2")
            if name == "2022_2024":
                sel = np.flatnonzero((dd >= LONG_DAYS) & (dd < 300))
                long_parts.append((d.coh.isel(pair=sel).values.astype("float32"), d.t2.values[sel]))
        a = load("annual", track)
        if a is not None:
            long_parts.append((a.coh.values.astype("float32"), a.t2.values))
        if long_parts:
            v = np.concatenate([p[0] for p in long_parts])
            s = season(np.concatenate([p[1] for p in long_parts]))
            for seas, sel in (("all", np.ones(len(s), bool)), ("winter", s == "winter"), ("snow-free", s == "snow-free")):
                if sel.any():
                    add("ansari2020", "a: coherence of pairs >= 48 d > 0.3", track, "long pairs", seas,
                        np.nanmean(v[sel], axis=0), COH_MIN, lambda v: v > COH_MIN, "M4")
    cond = pd.DataFrame(cond).astype({"median": float, "share_holds": float, "threshold": float})
    conn = pd.DataFrame(conn).astype({"median_coherence": float})
    pd.DataFrame(unw_n).to_csv(OUT / "unwrap_triplets.csv", index=False)
    verd = verdicts(cond)
    cond.to_csv(OUT / "conditions.csv", index=False)
    conn.to_csv(OUT / "coherence_by_connection.csv", index=False)
    verd.to_csv(OUT / "verdicts.csv", index=False)
    xr.Dataset({k: (("y", "x"), v.astype("float32")) for k, v in maps.items()},
               coords={"y": st.y, "x": st.x}).to_netcdf(OUT / "closure_maps.nc")
    (OUT / "README.md").write_text(readme(cond, conn, verd, pd.DataFrame(unw_n)), encoding="utf-8")
    print(f"wrote {OUT}")


def verdicts(cond: pd.DataFrame) -> pd.DataFrame:
    """The pre-registered rule: applicable if every measured condition holds in every set and season; partly if every
    one holds in some (set, season); otherwise not. maghsoudi2022 is capped (HyP3 input is filtered)."""
    rows = []
    for (corr, track, cls), g in cond.groupby(["correction", "track", "class"]):
        g = g.dropna(subset=["holds"])
        need = sorted({c.split(":")[0] for c in g.condition})
        everywhere = all(g[g.condition.str.startswith(c + ":")].holds.all() for c in need)
        # a (set, season) cell passes when every condition defined for it holds; set-level conditions apply to all seasons
        cells = []
        for (sname, seas), h in g[g.season != "all"].groupby(["set", "season"]):
            setwide = g[(g.set == sname) & (g.season == "all") & ~g.condition.isin(h.condition)]
            hh = pd.concat([h, setwide[setwide.condition.str.startswith("a: connections")]])
            if {c.split(":")[0] for c in hh.condition} == set(need) and hh.holds.all():
                cells.append(f"{sname} {seas}")
        if everywhere:
            v = "applicable"
        elif cells:
            v = "partly applicable"
        else:
            v = "not applicable"
        if corr == "maghsoudi2022" and v != "not applicable":
            v = "partly applicable (filtered input; needs LiCSAR or own SLC pairs)"
        failing = sorted({c for c in g[~g.holds.astype(bool)].condition})
        rows.append({"correction": corr, "track": track, "class": cls, "verdict": v,
                     "where_all_hold": "; ".join(cells) if cells else "", "conditions_failing": " | ".join(failing)})
    return pd.DataFrame(rows)


def _md(df: pd.DataFrame, nd: int = 2) -> str:
    d = df.copy()
    for c in d.columns:
        if pd.api.types.is_float_dtype(d[c]):
            d[c] = d[c].map(lambda x: "—" if pd.isna(x) else f"{x:.{nd}f}")
    return "\n".join(["| " + " | ".join(map(str, d.columns)) + " |", "|" + "---|" * len(d.columns)] +
                     ["| " + " | ".join(map(str, r)) + " |" for r in d.itertuples(index=False)])


def readme(cond: pd.DataFrame, conn: pd.DataFrame, verd: pd.DataFrame, unw_n: pd.DataFrame) -> str:
    key = verd[verd["class"].isin(["mat", "P6", "grassland", "hard_targets"])]
    floor = cond[(cond.season == "all") & cond["class"].isin(["mat", "open_water"])].pivot_table(
        index=["correction", "condition", "track", "set"], columns="class", values="median").reset_index()
    c_all = cond[cond.season == "all"].pivot_table(index=["correction", "condition", "track", "set"], columns="class",
                                                   values="share_holds").reset_index()
    cc = conn[(conn.season == "all")].pivot_table(index=["track", "set", "days"], columns="class",
                                                  values="median_coherence").reset_index()
    return f"""# X-097 — bias-correction pre-check (generated)

Generated by `experiments/bias_precheck_x097.py` (`make run E=bias_precheck_x097`); do not edit. Pre-registration:
`../preregistration/X-097.yaml` (conditions and verdict rule fixed there). No laser is read. Inputs: HyP3 40 m pairs
of the cube (Goldstein-filtered; coherence biased high at low coherence, so a pass on coherence is an upper bound).

## Verdict per correction

{_md(key)}

## Share of pixels where each condition holds (all seasons; P6 = 1 or 0)

{_md(c_all)}

## Reference: open water (added after the first run; not part of the pre-registered rule)

Median of each measure on the mat and on open water (WorldCover water share >= {WATER_SHARE}). Open water holds no
phase; where its value matches the mat's, the condition is not measurable on this product.

{_md(floor)}

## Median coherence by connection length (all seasons)

{_md(cc)}

## Unwrapped triplets with a hard-target reference

A triplet whose pairs have no unwrapped value on the hard targets cannot be referenced and is left out of M3.

{_md(unw_n)}

Files: `verdicts.csv`, `unwrap_triplets.csv`, `conditions.csv` (every condition by set, season, class: median, share of pixels where it holds,
holds), `coherence_by_connection.csv` (by season too), `closure_maps.nc` (|tau|, arg tau per set and level; share of
inconsistent unwrapped triplets per pixel).
"""


if __name__ == "__main__":
    main()

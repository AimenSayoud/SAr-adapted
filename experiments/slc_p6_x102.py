"""X-102 — at P6, does the 6-day phase from our own coregistered SLCs follow the laser better than HyP3's?

Four variants of each consecutive 2021 pair at P6, both tracks, as wrapped phase increments (cycle 0):
(a) the single CSLC sample nearest P6's surveyed position; (b) its 2 × 1 neighbours summed (10 m × 10 m, unfiltered);
(c) the P6 40 m cell, 32 looks, Goldstein-filtered like HyP3 (X-101 QC); (d) HyP3's wrapped phase at the P6 cell.
Each is referenced to its own interferogram's mean over the hard-target cells and scored against the laser by the
common scorer (C-061). Pre-registration: ``08_deliverables/insar_full_test/preregistration/X-102.yaml``. Caveat from
X-101 (DECISIONS 2026-10-13): the CSLCs are placed right, but their hard-target phase differs from HyP3's by processing.

Reads the CSLCs from the Drive root (``INSAR_DRIVE_ROOT``: on the Mac the mirror holding ``slc_pilot/2021/cslc``) and
the cube; writes only to the hub: ``08_deliverables/insar_full_test/slc_p6/``.

    INSAR_DRIVE_ROOT=<hub>/05_code/local/drive_mirror make run E=slc_p6_x102
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import p6_validation_x082 as X  # noqa: E402  (laser in LOS, PHASE_TO_MM)
import slc_qc_x101 as Q  # noqa: E402  (CSLC reading, cell looks, Goldstein, HyP3 pairs)

from insar_wetlands.cube import cube_root  # noqa: E402
from insar_wetlands.paths import hub_root, resolve_drive  # noqa: E402
from insar_wetlands.scoring import Protocol, assert_preregistered, score  # noqa: E402

HUB = hub_root()
CUBE = cube_root()
FT = HUB / "08_deliverables" / "insar_full_test"
OUT = FT / "slc_p6"
PREREG = FT / "preregistration" / "X-102.yaml"
VARIANTS = {"a": "single look", "b": "2 x 1 looks, unfiltered", "c": "40 m cell + Goldstein (ours)", "d": "HyP3 (wrapped)"}
MARGIN = 0.20


def circ_mean(z: np.ndarray) -> complex:
    z = z[np.isfinite(z)]
    return np.mean(z / np.maximum(np.abs(z), 1e-12)) if len(z) else np.nan


def variants(s1, s2, x, y, gx, gy, shape, hard, p6_xy, cells):
    """Reference-corrected phase (rad) of each variant at P6 and, for (c), at the nine cells."""
    rows, cols = Q.axis_index(x, y, gx, gy)
    z40, _ = Q.look(s1, s2, rows, cols, shape)
    zf = Q.goldstein(z40)
    ref_u, ref_f = circ_mean(z40[hard]), circ_mean(zf[hard])
    ix, iy = int(np.argmin(np.abs(x - p6_xy[0]))), int(np.argmin(np.abs(y - p6_xy[1])))
    ifg = s1 * np.conj(s2)
    a = ifg[iy, ix]
    b = ifg[iy, max(ix - 1, 0):ix + 1].sum()
    out = {"a": np.angle(a * np.conj(ref_u)), "b": np.angle(b * np.conj(ref_u))}
    r6, c6 = cells[4]
    out["c"] = np.angle(zf[r6, c6] * np.conj(ref_f))
    out["c_cells"] = [np.angle(zf[r, c] * np.conj(ref_f)) for r, c in cells]
    out["c_unfiltered"] = np.angle(z40[r6, c6] * np.conj(ref_u))
    return out, z40, zf


def main() -> None:
    assert_preregistered(PREREG, [OUT / "verdict.csv"])
    OUT.mkdir(parents=True, exist_ok=True)
    prot = Protocol.load(FT / "protocol.yaml")
    st = xr.open_dataset(CUBE / "silver" / "static_40m.nc")
    gx, gy = st.x.values, st.y.values
    shape = (len(gy), len(gx))
    hard = st.hard_target.values.astype(bool)
    px = pd.read_csv(HUB / "08_deliverables" / "field_first" / "plot_pixels.csv").set_index("plot")
    p6_xy = (float(px.loc["P6", "E"]), float(px.loc["P6", "N"]))
    r6, c6 = int(px.loc["P6", "row"]), int(px.loc["P6", "col"])
    cells = [(r6 + dr, c6 + dc) for dr in (-1, 0, 1) for dc in (-1, 0, 1)]
    gold = pd.read_csv(CUBE / "gold" / "plots_pair.csv.gz", low_memory=False,
                       usecols=["track", "grid", "plot", "pair", "d_laser_los_mm", "unw_1x1", "ref_hard_unw"])
    gold = gold[(gold["plot"] == "P6") & (gold.grid == "40m")].set_index(["track", "pair"])
    las = X.laser()[0]
    pair_rows, cell_rows, closure_rows, scores = [], [], [], []
    for track in X.TRACKS:
        folder = resolve_drive() / "slc_pilot" / "2021" / "cslc" / track
        files = {Q.DATE_RE.search(f.name).group(1): f for f in sorted(folder.glob("*.h5"))
                 if "static" not in f.name.lower() and Q.DATE_RE.search(f.name)}
        dates = sorted(files)
        hy = Q.hyp3_pairs(track)
        cache: dict = {}

        def get(d):
            if d not in cache:
                cache[d] = Q.cslc(files[d])
            return cache[d]

        z_store = {}
        for a_, b_ in zip(dates[:-1], dates[1:]):
            if (a_, b_) not in hy:
                continue
            (s1, x, y), (s2, _, _) = get(a_), get(b_)
            v, z40, zf = variants(s1, s2, x, y, gx, gy, shape, hard, p6_xy, cells)
            z_store[(a_, b_)] = (z40[r6, c6], zf[r6, c6])
            w_h = hy[(a_, b_)][0]
            v["d"] = np.angle(np.exp(1j * w_h[r6, c6]) * np.conj(circ_mean(np.exp(1j * w_h[hard]))))
            key = f"{a_}_{b_}"
            g = gold.loc[(track, key)] if (track, key) in gold.index else None
            rec = {"track": track, "t1": pd.Timestamp(a_, tz="UTC"), "t2": pd.Timestamp(b_, tz="UTC"),
                   "truth_inc_mm": float(g.d_laser_los_mm) if g is not None else np.nan,
                   "hyp3_unw_mm": float((g.unw_1x1 - g.ref_hard_unw) * X.PHASE_TO_MM) if g is not None else np.nan}
            for k in VARIANTS:
                rec[f"est_{k}_mm"] = float(v[k]) * X.PHASE_TO_MM
            rec["est_c_unfiltered_mm"] = float(v["c_unfiltered"]) * X.PHASE_TO_MM
            pair_rows.append(rec)
            for (r, c), ph in zip(cells, v["c_cells"]):
                cell_rows.append({**{k: rec[k] for k in ("track", "t1", "t2", "truth_inc_mm")},
                                  "row": r, "col": c, "est_c_mm": float(ph) * X.PHASE_TO_MM})
        # filter bias: closure of 6+6-12 triplets at the P6 cell, unfiltered vs filtered
        for a_, b_, c_ in zip(dates[:-2], dates[1:-1], dates[2:]):
            if (a_, b_) in z_store and (b_, c_) in z_store:
                (s1, x, y), (s3, _, _) = get(a_), get(c_)
                rows, cols = Q.axis_index(x, y, gx, gy)
                z13, _ = Q.look(s1, s3, rows, cols, shape)
                z13f = Q.goldstein(z13)
                cu = np.angle(z_store[(a_, b_)][0] * z_store[(b_, c_)][0] * np.conj(z13[r6, c6]))
                cf = np.angle(z_store[(a_, b_)][1] * z_store[(b_, c_)][1] * np.conj(z13f[r6, c6]))
                closure_rows.append({"track": track, "t1": a_, "t3": c_, "closure_unfiltered_rad": float(cu),
                                     "closure_filtered_rad": float(cf)})
        p = pd.DataFrame([r for r in pair_rows if r["track"] == track])
        for k, name in [*VARIANTS.items(), ("c_unfiltered", "40 m cell, unfiltered (ours)"), ("hyp3_unw", "HyP3 unwrapped (X-099)")]:
            col = f"est_{k}_mm" if k != "hyp3_unw" else "hyp3_unw_mm"
            q = p.dropna(subset=[col])
            sc = score(q.assign(est_mm=q[col]), las[track], prot)
            inc = sc["increments"].iloc[0].to_dict() if len(sc.get("increments", [])) else {}
            seas = sc["seasons"] if "seasons" in sc else pd.DataFrame()
            within = seas[seas.season.isin(["spring", "summer", "autumn"])] if len(seas) else seas
            scores.append({"track": track, "variant": k, "name": name, **inc,
                           "pearson_in_season_median": float(within.pearson_increments.median()) if len(within) else np.nan})
    pairs, cells_t, clos, sc = map(pd.DataFrame, (pair_rows, cell_rows, closure_rows, scores))
    v = verdict(sc)
    for name, df in (("pairs", pairs), ("cells", cells_t), ("closure", clos), ("scores", sc), ("verdict", v)):
        df.to_csv(OUT / f"{name}.csv", index=False)
    (OUT / "README.md").write_text(readme(sc, v, cells_t, clos), encoding="utf-8")
    print(v.to_string(index=False))
    print(f"wrote {OUT}")


def verdict(sc: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for k in ("a", "b"):
        rel, share_ok = [], []
        for track in X.TRACKS:
            ours = sc[(sc.track == track) & (sc.variant == k)].iloc[0]
            hy = sc[(sc.track == track) & (sc.variant == "d")].iloc[0]
            rel.append(ours.rmse_increments_mm / hy.rmse_increments_mm - 1)
            share_ok.append(ours.within_quarter_wave >= hy.within_quarter_wave)
        if all(r <= -MARGIN for r in rel) and all(share_ok):
            res = "beats HyP3"
        elif all(abs(r) <= MARGIN for r in rel):
            res = "matches HyP3"
        else:
            res = "loses to HyP3" if all(r > MARGIN for r in rel) else "mixed"
        rows.append({"variant": k, "name": VARIANTS[k], "rmse_change_asc": rel[0], "rmse_change_desc": rel[1],
                     "share_not_lower_both": all(share_ok), "result": res})
    return pd.DataFrame(rows)


def _md(df: pd.DataFrame, nd: int = 2) -> str:
    d = df.copy()
    for c in d.columns:
        if pd.api.types.is_float_dtype(d[c]):
            d[c] = d[c].map(lambda x: "—" if pd.isna(x) else f"{x:.{nd}f}")
    return "\n".join(["| " + " | ".join(map(str, d.columns)) + " |", "|" + "---|" * len(d.columns)] +
                     ["| " + " | ".join(map(str, r)) + " |" for r in d.itertuples(index=False)])


def readme(sc, v, cells_t, clos) -> str:
    cols = ["track", "name", "pairs", "pearson_increments", "rmse_increments_mm", "slope_through_origin",
            "within_quarter_wave", "pearson_in_season_median"]
    cm = cells_t.dropna(subset=["truth_inc_mm"]).groupby(["track", "row", "col"]).apply(
        lambda g: np.sqrt(np.mean((g.est_c_mm - g.truth_inc_mm) ** 2)), include_groups=False).rename(
        "rmse_mm").reset_index()
    cl = clos.groupby("track")[["closure_unfiltered_rad", "closure_filtered_rad"]].agg(
        lambda s: float(np.median(np.abs(s)))).reset_index()
    return f"""# X-102 — the 6-day phase at P6 from our own SLCs vs HyP3 (generated)

Generated by `experiments/slc_p6_x102.py` (`make run E=slc_p6_x102`); do not edit. Pre-registration:
`../preregistration/X-102.yaml`. 2021 consecutive pairs, both tracks, wrapped increments (cycle 0), each variant
referenced to its own hard-target mean, scored against the P6 laser by the common scorer. **Caveat (X-101, DECISIONS
2026-10-13):** the CSLCs are placed right but their hard-target phase differs from HyP3's by processing. Tried on P6,
2021 only (D-028).

## Verdict (a variant beats HyP3 if its increment RMSE is lower by >= {MARGIN:.0%} on both tracks, share within
lambda/4 not lower)

{_md(v)}

## Scores by variant

{_md(sc[[c for c in cols if c in sc.columns]])}

## The nine cells around P6, variant (c): increment RMSE (mm)

{_md(cm)}

## Closure of 6+6−12 triplets at the P6 cell (median |closure|, rad): the filter's effect

{_md(cl)}

Files: `verdict.csv`, `scores.csv`, `pairs.csv` (every pair, every variant), `cells.csv`, `closure.csv`.
"""


if __name__ == "__main__":
    main()

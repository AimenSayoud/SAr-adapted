"""X-134 — around P6, shrubs on the floating mat: does shrub cover (0.5–2 m, in cells without trees) change coherence,
backscatter or how the six-day phase follows the P6 laser, once depth into the mat is held fixed?

Cells: X-132's frame (40 m cells within 1 km of P6's cell, boardwalk share < 0.2, no hard targets) on the mat, with tree
share < 1 % (1 m canopy ≥ 2 m in the 2024 laser scan). Shrub share = the share of the cell's 1 m pixels at 0.5–2 m.
Outcomes O1–O5 per cell from X-132's code, unchanged. (i) Per track and set, least squares over the cells of each
outcome on shrub share (per 10 percentage points), depth into the mat, distance to P6 and boardwalk share; (ii) shrubby
(≥ 15 %) against open (< 3 %) cells matched on depth. 95 % intervals from 2000 resamples of the 200 m spatial blocks.

Pre-registration: hub ``08_deliverables/insar_full_test/preregistration/X-134.yaml`` (checked before the run). Stops
unless P6's own cell reproduces the cube's gold coherence and X-126's consecutive six-day r of 2021 and 2025 (1e-9).
Design years only. Reads the cube and the hub; writes only to the hub (``insar_full_test/shrub_x134/``).

    make run E=shrub_coherence_x134
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import landcover_coherence_x132 as X132  # noqa: E402  (X-132's frame, pairs and per-cell outcomes)
from insar_wetlands import contrasts as C  # noqa: E402
from insar_wetlands.paths import cube_root, hub_root  # noqa: E402
from insar_wetlands.scoring import Protocol, assert_preregistered, guard_lockbox  # noqa: E402

HUB, CUBE = hub_root(), cube_root()
FT = HUB / "08_deliverables" / "insar_full_test"
OUT = FT / "shrub_x134"
PREREG = FT / "preregistration" / "X-134.yaml"
PROTOCOL = FT / "protocol.yaml"
TRACKS, SETS = X132.TRACKS, X132.SETS
TREE_M, SHRUB_M, TALL_M = 2.0, 0.5, 1.0
TREE_FREE = 0.01
OPEN_LT, SHRUBBY_GE = 0.03, 0.15
NINE_RANGE = 1.5                     # 0 → 15 % shrub, in units of 10 percentage points
MARGIN = 0.03
CALIPER_M, N_BOOT, SEED = 40.0, 2000, 0
COVARIATES = ("depth_in_mat_m", "dist_p6_m", "boardwalk_share")
NINE = {(-1, -1): "NW", (-1, 0): "N", (-1, 1): "NE", (0, -1): "W", (0, 0): "P6", (0, 1): "E",
        (1, -1): "SW", (1, 0): "S", (1, 1): "SE"}


# ------------------------------------------------------------------------------ land cover

def cover_shares(cells: pd.DataFrame) -> pd.DataFrame:
    """Per 40 m cell, the shares of its 1 m canopy pixels in each height band, 2024 (leaf-on) and 2011 (leaf-off)."""
    li = xr.open_dataset(CUBE / "silver" / "lidar_dated_1m.nc")
    out = {}
    for yr in (2024, 2011):
        chm = li[f"chm_{yr}_m"]
        rows = []
        for x, y in zip(cells.x, cells.y):
            v = chm.sel(x=slice(x - 20, x + 20), y=slice(y + 20, y - 20)).values
            v = v[np.isfinite(v)]
            if not v.size:
                rows.append((np.nan, np.nan, np.nan, 0))
                continue
            rows.append(((v >= TREE_M).mean(), ((v >= SHRUB_M) & (v < TREE_M)).mean(),
                         ((v >= TALL_M) & (v < TREE_M)).mean(), v.size))
        a = np.array(rows, float)
        out[f"tree_{yr}"], out[f"shrub_{yr}"], out[f"shrub_tall_{yr}"], out[f"n1m_{yr}"] = a.T
    return pd.DataFrame(out, index=cells.index)


def classify(s: pd.Series, open_lt: float = OPEN_LT, shrubby_ge: float = SHRUBBY_GE) -> np.ndarray:
    return np.where(s < open_lt, "open", np.where(s >= shrubby_ge, "shrubby", np.where(s.notna(), "mid", "none")))


# ------------------------------------------------------------------------------ (i) dose–response

def dose(cells: pd.DataFrame, v: pd.Series, shrub_col: str = "shrub_2024", adjust: bool = True) -> dict:
    """Least squares over the cells of ``v`` on shrub share (per 10 points) and, when ``adjust``, the standardised
    covariates; the shrub coefficient, its 0 → 15 % effect and the shrub share of R², with block-bootstrap intervals."""
    d = cells.assign(y=v.reindex(cells.index)).dropna(subset=["y", shrub_col, *COVARIATES])
    X = [np.ones(len(d)), d[shrub_col].to_numpy(float) * 10]
    if adjust:
        for c in COVARIATES:
            x = d[c].to_numpy(float)
            X.append((x - x.mean()) / x.std() if x.std() > 0 else np.zeros(len(d)))
    X = np.column_stack(X)
    S = C.block_sums(X, d.y.to_numpy(float), d.spatial_block.to_numpy())
    res = C.ols_block_bootstrap(S, {"shrub": [1]}, N_BOOT, SEED)
    lo, hi = C.boot_interval(res["beta_draws"][:, 1])
    dlo, dhi = C.boot_interval(res["drop_draws"]["shrub"])
    return {"cells": len(d), "blocks": len(S["blocks"]), "per_10pts": res["beta"][1], "lo95": lo, "hi95": hi,
            "effect_0_15": NINE_RANGE * res["beta"][1], "effect_lo95": NINE_RANGE * lo,
            "effect_hi95": NINE_RANGE * hi, "r2": res["r2"], "shrub_r2_drop": res["drop"]["shrub"],
            "shrub_r2_drop_lo95": dlo, "shrub_r2_drop_hi95": dhi, "beta": res["beta"], "index": d.index, "X": X}


def dose_table(cells: pd.DataFrame, oc: pd.DataFrame, variant: str, shrub_col: str = "shrub_2024",
               adjust: bool = True, outcomes=None) -> tuple[pd.DataFrame, dict]:
    wide = oc.pivot_table(index="cell", columns=["track", "set", "outcome"], values="value", dropna=False)
    wide = wide.dropna(axis=1, how="all")
    rows, fits = [], {}
    for key in wide.columns:
        track, st, outcome = key
        if outcomes is not None and outcome not in outcomes:
            continue
        f = dose(cells, wide[key], shrub_col, adjust)
        fits[key] = f
        rows.append({"variant": variant, "track": track, "set": st, "outcome": outcome,
                     **{k: f[k] for k in f if k not in ("beta", "index", "X")}})
    return pd.DataFrame(rows), fits


# ------------------------------------------------------------------------------ (ii) matched contrast

def contrast(cells: pd.DataFrame, oc: pd.DataFrame, variant: str, cls_col: str = "cls",
             outcomes=None) -> pd.DataFrame:
    s, o = cells[cells[cls_col] == "shrubby"], cells[cells[cls_col] == "open"]
    M = C.match_nearest(s, o, on="depth_in_mat_m", caliper=CALIPER_M, tie="dist_p6_m", order="dist_p6_m")
    M["block"] = cells.loc[M.treated, "spatial_block"].to_numpy()
    ub, counts = C.block_counts(M.block.to_numpy(), N_BOOT, SEED)
    W = C.unit_weights(M.block.to_numpy(), ub, counts)
    wide = oc.pivot_table(index="cell", columns=["track", "set", "outcome"], values="value", dropna=False)
    wide = wide.dropna(axis=1, how="all")
    rows = []
    for key in wide.columns:
        track, st, outcome = key
        if outcomes is not None and outcome not in outcomes:
            continue
        v = wide[key]
        d = v.reindex(M.treated).to_numpy() - v.reindex(M.control).to_numpy()
        est = C.weighted_median(d, np.ones(len(d)))
        lo, hi = C.boot_interval([C.weighted_median(d, w) for w in W])
        rows.append({"variant": variant, "track": track, "set": st, "outcome": outcome, "shrubby": len(s),
                     "open": len(o), "pairs": int(np.isfinite(d).sum()), "shrubby_minus_open": est,
                     "lo95": lo, "hi95": hi})
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------------ reading

def reading(dt: pd.DataFrame, ct: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """The pre-registered reading rules on the primary variant, and (ii) against (i) set by set."""
    p = dt[(dt.variant == "primary") & (dt.outcome == "O1_coh")]

    def count(test):
        k = 0
        for st in SETS:
            g = p[p.set == st]
            if len(g) == 2 and all(test(r) for r in g.itertuples()):
                k += 1
        return k
    below = lambda r: r.effect_hi95 < 0  # noqa: E731
    above = lambda r: r.effect_lo95 > 0  # noqa: E731
    inside = lambda r: r.effect_lo95 > -MARGIN and r.effect_hi95 < MARGIN  # noqa: E731
    o3 = dt[(dt.variant == "primary") & (dt.outcome == "O3_r_laser") & (dt.set == "2021+2025")]
    nb, na, ni = count(below), count(above), count(inside)
    verdict = ("shrubs lower coherence on the mat" if nb >= 2 else "shrubs raise coherence" if na >= 2 else
               "no material difference over the nine cells' range" if ni >= 2 else "not determined at this precision")
    rows = [
        {"rule": "shrubs lower coherence on the mat", "test": "0→15 % interval below zero, both tracks, ≥ 2 of 3 sets",
         "sets_meeting": nb, "holds": nb >= 2},
        {"rule": "shrubs raise coherence", "test": "0→15 % interval above zero, both tracks, ≥ 2 of 3 sets",
         "sets_meeting": na, "holds": na >= 2},
        {"rule": "no material difference over the nine cells' range",
         "test": f"0→15 % interval inside ±{MARGIN}, both tracks, ≥ 2 of 3 sets", "sets_meeting": ni, "holds": ni >= 2},
        {"rule": "shrubs weaken the six-day motion read", "test": "0→15 % O3 interval below zero, both tracks",
         "sets_meeting": int(len(o3) == 2 and all(o3.effect_hi95 < 0)),
         "holds": bool(len(o3) == 2 and all(o3.effect_hi95 < 0))},
        {"rule": "verdict on O1", "test": "the first of the three O1 rules that holds, else not determined",
         "sets_meeting": "", "holds": verdict}]
    q = ct[(ct.variant == "primary") & (ct.outcome == "O1_coh")].set_index(["track", "set"])
    agree = []
    for r in p.itertuples():
        c = q.loc[(r.track, r.set)]
        excl_i = r.effect_lo95 > 0 or r.effect_hi95 < 0
        excl_ii = c.lo95 > 0 or c.hi95 < 0
        agree.append({"track": r.track, "set": r.set, "dose_0_15": r.effect_0_15, "matched": c.shrubby_minus_open,
                      "same_sign": bool(np.sign(r.effect_0_15) == np.sign(c.shrubby_minus_open)),
                      "dose_excludes_zero": excl_i, "matched_excludes_zero": excl_ii,
                      "agree": bool(np.sign(r.effect_0_15) == np.sign(c.shrubby_minus_open) and excl_i == excl_ii)})
    return pd.DataFrame(rows), pd.DataFrame(agree)


def nine_cells(cells: pd.DataFrame, oc: pd.DataFrame, fits: dict, p6: tuple[int, int]) -> pd.DataFrame:
    """P6 and its eight neighbours on the primary O1 fit: observed against predicted (cells in the frame only)."""
    rows = []
    for (dr, dc), name in NINE.items():
        cell = f"r{p6[0] + dr}c{p6[1] + dc}"
        rec = {"cell": name, "id": cell, "in_frame": cell in cells.index}
        if cell in cells.index:
            rec |= {"shrub_2024": cells.loc[cell, "shrub_2024"], "tree_2024": cells.loc[cell, "tree_2024"],
                    "depth_in_mat_m": cells.loc[cell, "depth_in_mat_m"]}
            for track in TRACKS:
                for st in SETS:
                    f = fits.get((track, st, "O1_coh"))
                    if f is None or cell not in f["index"]:
                        continue
                    i = f["index"].get_loc(cell)
                    obs = oc[(oc.cell == cell) & (oc.track == track) & (oc.set == st) & (oc.outcome == "O1_coh")]
                    rec[f"{track[:3]} {st} obs"] = float(obs.value.iloc[0])
                    rec[f"{track[:3]} {st} pred"] = float(f["X"][i] @ f["beta"])
        rows.append(rec)
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------------ figures

def fig_dose(cells: pd.DataFrame, oc: pd.DataFrame, fits: dict, p6: tuple[int, int], out: Path) -> None:
    """O1 against shrub share, adjusted for the covariates (partial residuals), with the fitted line; nine cells marked."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    nine = [f"r{p6[0] + dr}c{p6[1] + dc}" for dr, dc in NINE]
    fig, axs = plt.subplots(2, 3, figsize=(15, 9), sharex=True)
    for i, track in enumerate(TRACKS):
        for j, st in enumerate(SETS):
            ax = axs[i, j]
            f = fits.get((track, st, "O1_coh"))
            if f is None:
                ax.set_axis_off()
                continue
            v = oc[(oc.track == track) & (oc.set == st) & (oc.outcome == "O1_coh")].set_index("cell").value
            y = v.reindex(f["index"]).to_numpy(float)
            x = f["X"][:, 1]
            part = y - f["X"][:, 2:] @ f["beta"][2:]
            ax.scatter(x * 10, part, s=10, c="#7f7f7f", alpha=0.6, label="tree-free mat cells")
            m = np.isin(f["index"], nine)
            ax.scatter(x[m] * 10, part[m], s=40, facecolors="none", edgecolors="#d95f02", label="the nine cells")
            xs = np.linspace(0, x.max(), 50)
            ax.plot(xs * 10, f["beta"][0] + f["beta"][1] * xs, c="#1b9e77", lw=2, label="fit")
            ax.axvspan(0, 15, color="#d95f02", alpha=0.07)
            r = f["beta"][1] * NINE_RANGE
            ax.set_title(f"{track} {st}: 0→15 % shrub {r:+.3f}", fontsize=9)
            if i == 1:
                ax.set_xlabel("shrub cover 0.5–2 m (% of the 40 m cell, 2024 scan)")
            if j == 0:
                ax.set_ylabel("median coherence, adjusted for depth,\ndistance to P6 and boardwalk")
    axs[0, 0].legend(fontsize=7, loc="lower right")
    fig.suptitle("X-134: coherence of consecutive pairs against shrub cover, tree-free cells on the mat within 1 km of "
                 "P6 (shaded: the nine cells' range)", fontsize=10)
    fig.savefig(out, dpi=120, bbox_inches="tight")
    plt.close(fig)


def fig_map(all_cells: pd.DataFrame, cells: pd.DataFrame, p6: tuple[int, int], out: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8, 8))
    on = all_cells[all_cells.on_mat == 1]
    ax.scatter(all_cells.x, all_cells.y, s=9, marker="s", c="#e0e0e0", label="other cells in X-132's frame")
    ax.scatter(on.x, on.y, s=9, marker="s", c="#bdbdbd", label="mat cells with trees ≥ 1 %")
    sc = ax.scatter(cells.x, cells.y, s=12, marker="s", c=cells.shrub_2024 * 100, cmap="YlGn", vmin=0, vmax=60,
                    label="tree-free mat cells (colour: shrub %)")
    me = all_cells[(all_cells.row == p6[0]) & (all_cells["col"] == p6[1])]
    ax.scatter(me.x, me.y, marker="*", s=120, c="w", edgecolors="k", label="P6", zorder=5)
    fig.colorbar(sc, ax=ax, shrink=0.7, label="shrub cover 0.5–2 m (%), 2024 scan")
    ax.set_aspect("equal")
    ax.legend(fontsize=7, loc="lower left")
    ax.set_title("X-134 cells within 1 km of P6 (UTM 33N, m)", fontsize=9)
    fig.savefig(out, dpi=120, bbox_inches="tight")
    plt.close(fig)


# ------------------------------------------------------------------------------ main

def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    names = ["cells.csv", "outcomes.csv.gz", "dose.csv", "matched.csv", "reading.csv", "agreement.csv",
             "nine_cells.csv", "reproduction.csv", "README.md", "fig_dose_response.png", "fig_shrub_map.png"]
    assert_preregistered(PREREG, [OUT / f for f in names])
    allc, p6 = X132.load_cells()
    allc = allc.join(cover_shares(allc))
    off = allc[(allc.on_mat == 0) & (allc.tree_2024 < TREE_FREE)]
    off_note = f"{int((off.shrub_2024 >= 0.10).sum())} of {len(off)}"
    cells = allc[(allc.on_mat == 1) & (allc.tree_2024 < TREE_FREE)].copy()
    me = f"r{p6[0]}c{p6[1]}"
    if me not in cells.index:
        sys.exit("stop (pre-registered): P6's cell is outside the frame")
    cells["cls"] = classify(cells.shrub_2024)
    cells["cls_2011"] = classify(cells.shrub_2011)
    cells["cls_loose"] = classify(cells.shrub_2024, 0.02, 0.10)
    cells["cls_tall"] = classify(cells.shrub_tall_2024)
    X132.P6, X132.NDVI_MONTH = p6, X132.ndvi_month(cells)
    g = X132.gold_p6()
    guard_lockbox(pd.concat([g.t1, g.t2]), Protocol.load(PROTOCOL))
    oc, _, checks = X132.outcomes(cells, g, X132.laser_pairs_x126())
    rep = pd.DataFrame(checks)
    rep.to_csv(OUT / "reproduction.csv", index=False)
    if not len(rep) or not rep.ok.all():
        print(rep.to_string())
        sys.exit("stop (pre-registered): P6's own cell does not reproduce the cube / X-126; see reproduction.csv")
    cells.to_csv(OUT / "cells.csv")
    oc.to_csv(OUT / "outcomes.csv.gz", index=False)

    # (i) and (ii): primary, unadjusted, and the pre-registered sensitivities (O1 and O3 there)
    d_prim, fits = dose_table(cells, oc, "primary")
    d_unadj, _ = dose_table(cells, oc, "unadjusted", adjust=False, outcomes={"O1_coh", "O3_r_laser"})
    c_prim = contrast(cells, oc, "primary")
    key = {"O1_coh", "O3_r_laser"}
    stable = cells[(cells.cls == cells.cls_2011) & (cells.tree_2011 < TREE_FREE)]
    sens = [("open < 2 % / shrubby ≥ 10 %", cells, "shrub_2024", "cls_loose"),
            ("tall shrub 1–2 m", cells, "shrub_tall_2024", "cls_tall"),
            ("tree share < 0.1 %", cells[cells.tree_2024 < 0.001], "shrub_2024", "cls"),
            ("leaf-off 2011 agrees", stable, "shrub_2024", "cls"),
            ("CIR-stable 2021–2023", cells[cells.cir_stable], "shrub_2024", "cls"),
            ("radius 500 m", cells[cells.dist_p6_m <= 500], "shrub_2024", "cls")]
    dts, cts = [d_prim, d_unadj], [c_prim]
    for name, cc, scol, ccol in sens:
        o = oc[oc.cell.isin(cc.index)]
        dts.append(dose_table(cc, o, name, scol, outcomes=key)[0])
        cts.append(contrast(cc, o, name, ccol, outcomes=key))
        print(f"sensitivity {name}: {len(cc)} cells", flush=True)
    dt, ct = pd.concat(dts, ignore_index=True), pd.concat(cts, ignore_index=True)
    dt.to_csv(OUT / "dose.csv", index=False)
    ct.to_csv(OUT / "matched.csv", index=False)
    rd, ag = reading(dt, ct)
    rd.to_csv(OUT / "reading.csv", index=False)
    ag.to_csv(OUT / "agreement.csv", index=False)
    nine = nine_cells(cells, oc, fits, p6)
    nine.to_csv(OUT / "nine_cells.csv", index=False)
    fig_dose(cells, oc, fits, p6, OUT / "fig_dose_response.png")
    fig_map(allc, cells, p6, OUT / "fig_shrub_map.png")

    md = X132.md
    cols = ["track", "set", "outcome", "cells", "per_10pts", "lo95", "hi95", "effect_0_15", "effect_lo95",
            "effect_hi95", "shrub_r2_drop"]
    o1 = d_prim[d_prim.outcome.str.startswith("O1")][cols]
    o3 = d_prim[d_prim.outcome == "O3_r_laser"][cols]
    rest = d_prim[~d_prim.outcome.str.startswith(("O1", "O3"))][cols]
    s1 = dt[(dt.outcome == "O1_coh")].pivot_table(index=["track", "set"], columns="variant",
                                                    values="effect_0_15").reset_index()
    s2 = ct[(ct.outcome == "O1_coh")].pivot_table(index=["track", "set"], columns="variant",
                                                    values="shrubby_minus_open").reset_index()
    n_s = ct.drop_duplicates("variant")[["variant", "shrubby", "open", "pairs"]]
    sh = cells.shrub_2024
    (OUT / "README.md").write_text(f"""# X-134 — shrubs on the mat around P6 (generated)

Generated by `experiments/shrub_coherence_x134.py` (`make run E=shrub_coherence_x134`); do not edit.
Pre-registration: `../preregistration/X-134.yaml`. Cells: X-132's frame (40 m, within 1000 m of P6, boardwalk share
< 0.2, no hard targets) on the mat, tree share < {TREE_FREE:.0%} (1 m canopy ≥ 2 m, 2024 scan): {len(cells)} cells, shrub
share (0.5–2 m) median {sh.median():.2f}, range {sh.min():.2f}–{sh.max():.2f}, correlation with depth
{sh.corr(cells.depth_in_mat_m):.2f}. Off the mat, tree-free cells with shrub share ≥ 0.10: {off_note} (not tested).
(i) least squares of each cell's outcome on shrub share (per 10 percentage points), depth into the mat, distance to P6
and boardwalk share; the 0 → 15 % effect is the nine cells' range. (ii) shrubby (≥ {SHRUBBY_GE:.0%}) − open
(< {OPEN_LT:.0%}), matched on depth (caliper {CALIPER_M:.0f} m). {N_BOOT} resamples of the 200 m spatial blocks, seed {SEED}.
Design years only.

## Reading rules (pre-registered)

{md(rd)}

## (ii) against (i), O1, set by set

{md(ag)}

## Reproduction check (P6's own cell)

{md(rep, 12)}

## (i) O1: coherence of consecutive pairs (all, summer, winter)

{md(o1)}

## (i) O3: how the cell's six-day phase follows the P6 laser

{md(o3)}

## (i) O2, O4, O5

{md(rest)}

## (ii) matched contrast, primary

{md(c_prim.drop(columns=["variant"]))}

## Sensitivities: O1, the 0 → 15 % effect of (i) per variant

{md(s1)}

## Sensitivities: O1, shrubby − open of (ii) per variant

{md(s2)}

{md(n_s)}

## The nine cells around P6 on the primary fit (O1)

{md(nine)}

## Figures

- `fig_dose_response.png` — O1 against shrub share, adjusted for the covariates, per track and set; the nine cells
  circled; the shaded band is their range (0–15 %).
- `fig_shrub_map.png` — the tested cells coloured by shrub share within 1 km of P6.
""", encoding="utf-8")
    print(rd.to_string())
    print(o1[o1.outcome == "O1_coh"].to_string())


if __name__ == "__main__":
    main()

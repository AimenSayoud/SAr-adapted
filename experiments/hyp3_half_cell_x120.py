"""X-120 — is the HyP3 half-cell shift a labelling error of the products? Checked against references independent of HyP3
and of our SLCs (GUGiK LiDAR 1 m, ESA WorldCover 10 m), against the shift predicted from the processing code.

(1) Per year set and track, HyP3's median coherence is modelled by the land cover of the 40 m (or 20 m) cell centred at
L + (dx, dy) — open water, trees (WorldCover), canopy height (LiDAR) — and the shift with the best R² is found, with a
bootstrap over pairs. (2) The products' own DEM layer is matched to the LiDAR DSM the same way; the code predicts it is
placed right. Rule: ``preregistration/X-120.yaml``. Writes to the hub, ``insar_full_test/hyp3_half_cell_confirm/``.

    make run E=hyp3_half_cell_x120
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from insar_wetlands import placement as P
from insar_wetlands.cube import cube_root
from insar_wetlands.paths import hub_root
from insar_wetlands.scoring import assert_preregistered

HUB = hub_root()
FT = HUB / "08_deliverables" / "insar_full_test"
OUT = FT / "hyp3_half_cell_confirm"
PREREG = FT / "preregistration" / "X-120.yaml"
SETS = ("2017_2019", "2020_2021", "2022_2024", "2025", "2026", "2026_20m", "annual", "matched_wt")
TRACKS = ("ascending", "descending")
STEPS = np.arange(-40.0, 40.01, 2.5)
TOL = 7.5
N_BOOT = 1000
LAT_SOUTH, LAT_SITE = 52.6, 52.763      # the burst's southern edge (approx.; the prediction moves < 0.1 m over 52.5–52.7)
SEED = 120


def references():
    """Fine references on one 2.5 m grid over the LiDAR footprint: water and tree fractions (WorldCover 2021, nearest
    from 10 m), canopy height and DSM (LiDAR, 2.5 m block means of the 1 m rasters)."""
    c = cube_root()
    li = xr.open_dataset(c / "silver" / "lidar_1m.nc")
    dsm, dtm = li.dsm_m.values.astype(float), li.dtm_m.values.astype(float)
    lx, ly = li.x.values, li.y.values
    st = xr.open_dataset(c / "silver" / "static_40m.nc")
    tx, ty = st.x.values, st.y.values
    wc = np.load(sorted((c / "raw" / "worldcover").glob("*2021*.npz"))[0])["classes"]
    wx0, wy0 = tx[0] - 20 + 5, ty[0] + 20 - 5                       # WorldCover 10 m cell centres (cube build)
    gx = np.arange(lx[0] - 0.5 + 1.25, lx[-1] + 0.5, 2.5)            # 2.5 m grid inside the LiDAR footprint
    gy = np.arange(ly[0] + 0.5 - 1.25, ly[-1] - 0.5, -2.5)
    iw = np.clip(np.round((gx - wx0) / 10).astype(int), 0, wc.shape[1] - 1)
    jw = np.clip(np.round((wy0 - gy) / 10).astype(int), 0, wc.shape[0] - 1)
    cl = wc[np.ix_(jw, iw)]
    blk = lambda a: P.FineField(a, lx, ly).box_mean(*np.meshgrid(gx, gy), 1.25, min_cover=0.5)  # noqa: E731
    canopy = np.clip(dsm - dtm, 0, 40)
    return gx, gy, {"water": (cl == 80).astype(float), "trees": (cl == 10).astype(float),
                    "canopy": blk(canopy), "dsm": blk(dsm)}


def ifg(set_: str, track: str) -> xr.Dataset | None:
    name = f"ifg_2026_{track}_20m.nc" if set_ == "2026_20m" else f"ifg_{set_}_{track}.nc"
    f = cube_root() / "silver" / name
    return xr.open_dataset(f) if f.exists() else None


def coherence_case(ds: xr.Dataset, fields: dict, rng) -> dict:
    res = float(ds.x[1] - ds.x[0])
    cx, cy = np.meshgrid(ds.x.values, ds.y.values)
    cx, cy = cx.ravel(), cy.ravel()
    feats = {s: np.stack([fields[k].box_mean(cx + s[0], cy + s[1], res / 2) for k in ("water", "trees", "canopy")], 1)
             for s in [(dx, dy) for dy in STEPS for dx in STEPS]}
    keys = list(feats)
    ok = np.all([np.all(np.isfinite(feats[k]), 1) for k in keys], 0)
    coh = ds.coh.values.reshape(ds.sizes["pair"], -1)[:, ok].astype(float)
    ok_c = np.isfinite(np.nanmedian(coh, 0))
    coh = coh[:, ok_c]
    # one orthonormal basis of the centred features per shift: R² = |Qᵀ ỹ|² / |ỹ|²
    Q = []
    for k in keys:
        X = feats[k][ok][ok_c]
        X = X - X.mean(0)
        q, _ = np.linalg.qr(X)
        Q.append(q.T)
    Q = np.stack(Q)                                                  # shifts × 3 × cells

    def r2(Y):                                                       # Y: cells × B
        Yc = Y - Y.mean(0)
        proj = np.einsum("skc,cb->skb", Q, Yc)
        return (proj ** 2).sum(1) / (Yc ** 2).sum(0)                 # shifts × B

    y0 = np.nanmedian(coh, 0)[:, None]
    s0 = r2(y0)[:, 0]
    kbest = int(np.argmax(s0))
    kzero = keys.index((0.0, 0.0))
    npair = coh.shape[0]
    boots = []
    for _ in range(N_BOOT // 50):
        Y = np.stack([np.nanmedian(coh[rng.integers(0, npair, npair)], 0) for _ in range(50)], 1)
        boots.append(np.argmax(r2(Y), 0))
    bk = np.concatenate(boots)
    bdx = np.array([keys[k][0] for k in bk])
    bdy = np.array([keys[k][1] for k in bk])
    return {"pairs": npair, "cells": int(coh.shape[1]), "dx_m": keys[kbest][0], "dy_m": keys[kbest][1],
            "dx_lo": float(np.percentile(bdx, 2.5)), "dx_hi": float(np.percentile(bdx, 97.5)),
            "dy_lo": float(np.percentile(bdy, 2.5)), "dy_hi": float(np.percentile(bdy, 97.5)),
            "r2_best": float(s0[kbest]), "r2_zero": float(s0[kzero])}


def dem_case(ds: xr.Dataset, dsm: P.FineField) -> dict:
    res = float(ds.x[1] - ds.x[0])
    cx, cy = np.meshgrid(ds.x.values, ds.y.values)
    dem = ds.dem.values.ravel().astype(float)
    dem[dem == 0] = np.nan
    dx, dy, r, _ = P.best_shift(lambda a, b: P.corr(dem, dsm.box_mean(cx.ravel() + a, cy.ravel() + b, res / 2)), STEPS)
    r0 = P.corr(dem, dsm.box_mean(cx.ravel(), cy.ravel(), res / 2))
    return {"dx_m": dx, "dy_m": dy, "r_best": r, "r_zero": r0}


def main() -> None:
    assert_preregistered(PREREG, [OUT / "verdict.csv"])
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    gx, gy, ref = references()
    F = {k: P.FineField(v, gx, gy) for k, v in ref.items()}
    pred40 = P.hyp3_isce2_predicted_shift(30.0, LAT_SOUTH, LAT_SITE)
    pred20 = P.hyp3_isce2_predicted_shift(20.0, LAT_SOUTH, LAT_SITE)
    rows, dems = [], []
    for s in SETS:
        for t in TRACKS:
            ds = ifg(s, t)
            if ds is None or ds.sizes["pair"] < 10:
                rows.append({"set": s, "track": t, "pairs": 0 if ds is None else ds.sizes["pair"], "counted": False})
                continue
            pdx, pdy = pred20 if s.startswith("2026") else pred40   # the 40 m 2026 set is averaged from 20 m
            r = coherence_case(ds, F, rng)
            r.update({"set": s, "track": t, "pred_dx_m": round(pdx, 2), "pred_dy_m": round(pdy, 2), "counted": True,
                      "within_tol": abs(r["dx_m"] - pdx) <= TOL and abs(r["dy_m"] - pdy) <= TOL,
                      "ci_contains_zero": r["dx_lo"] <= 0 <= r["dx_hi"] and r["dy_lo"] <= 0 <= r["dy_hi"]})
            rows.append(r)
            if s in ("2022_2024", "2026_20m"):
                dems.append({"set": s, "track": t, **dem_case(ds, F["dsm"])})
            print(s, t, {k: r[k] for k in ("pairs", "dx_m", "dy_m", "r2_best", "r2_zero", "within_tol")}, flush=True)
    C, D = pd.DataFrame(rows), pd.DataFrame(dems)
    cnt = C[C.counted]
    dem_ok = bool(((D.dx_m.abs() <= TOL) & (D.dy_m.abs() <= TOL)).all())
    share = float(cnt.within_tol.mean())
    any_zero = bool(cnt.ci_contains_zero.any())
    confirmed = dem_ok and share >= 0.8 and not any_zero
    V = pd.DataFrame([{"dem_placed_right": dem_ok, "share_within_tol": share, "any_interval_contains_zero": any_zero,
                       "cases": len(cnt), "verdict": "confirmed" if confirmed else "not confirmed",
                       "pred_40m_dx_m": round(pred40[0], 2), "pred_40m_dy_m": round(pred40[1], 2),
                       "pred_20m_dx_m": round(pred20[0], 2), "pred_20m_dy_m": round(pred20[1], 2)}])
    t75 = pd.read_csv(FT / "hyp3_placement" / "decisive" / "verdict.csv")
    for name, df in (("cases", C), ("dem", D), ("verdict", V)):
        df.to_csv(OUT / f"{name}.csv", index=False)
    md = lambda df: "\n".join(["| " + " | ".join(map(str, df.columns)) + " |", "|" + "---|" * len(df.columns)] +  # noqa: E731
                              ["| " + " | ".join(f"{x:.3f}" if isinstance(x, float) else str(x) for x in r) + " |"
                               for r in df.itertuples(index=False)])
    (OUT / "README.md").write_text(f"""# X-120 — the HyP3 half-cell shift against independent references (generated)

Generated by `experiments/hyp3_half_cell_x120.py`; do not edit. Rule: `../preregistration/X-120.yaml`. The content of a
cell labelled L lies at L + (dx, dy) (east, north, metres). Prediction from the processing code (hyp3_isce2 3.0.2 + ISCE2
geocode: values computed at the corner of their geocoding cell, labelled as the cell; the cell is the DEM grid, square
in degrees): `placement.hyp3_isce2_predicted_shift`. References: GUGiK LiDAR 1 m (DSM, canopy height), ESA WorldCover
2021 10 m (water, trees), both independent of HyP3 and of our SLCs. No laser data read.

## Verdict

{md(V)}

## Coherence against land cover, per year set and track

Best shift by R² of median coherence on water fraction, tree fraction and canopy height of the shifted cell; 95% interval
from {N_BOOT} bootstrap resamples of the pairs; 2.5 m steps within ±40 m.

{md(C[[c for c in ["set", "track", "pairs", "cells", "dx_m", "dx_lo", "dx_hi", "dy_m", "dy_lo", "dy_hi", "pred_dx_m", "pred_dy_m", "r2_best", "r2_zero", "within_tol", "ci_contains_zero"] if c in C]].fillna("")) }

## The products' own DEM layer against the LiDAR DSM (predicted placed right)

{md(D)}

## T75 (our 2021 SLC coherence as the reference), for comparison

{md(t75)}
""", encoding="utf-8")
    print(V.to_string(index=False))


if __name__ == "__main__":
    sys.exit(main())

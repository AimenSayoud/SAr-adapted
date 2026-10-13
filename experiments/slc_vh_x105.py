"""X-105 — cross-polarised (VH) 6-day interferograms beside VV from our 2021 CSLCs: by scattering regime, and at P6.

ledger2024 used cross-polarised interferograms to separate multiple-bounce from double-bounce returns on a peat swamp;
HyP3 is VV only. From X-101's VV and VH CSLCs (same grid, placed right): at P6, the VH and VV consecutive-pair increments
against the laser and against the water table; on the mat and at the plot cells, VV and VH coherence and the VV − VH
phase difference by X-098's regime class. Rules: ``08_deliverables/insar_full_test/preregistration/X-105.yaml``.

Reads the CSLCs from the Drive root (on the Mac the mirror: ``slc_pilot/2021/cslc`` and ``cslc_vh``) and the cube;
writes only to the hub: ``08_deliverables/insar_full_test/slc_vh/``.

    INSAR_DRIVE_ROOT=<hub>/05_code/local/drive_mirror make run E=slc_vh_x105
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import p6_validation_x082 as X  # noqa: E402
import scattering_regime_x098 as R  # noqa: E402  (regime classes per date)
import slc_qc_x101 as Q  # noqa: E402

from insar_wetlands import slc_stats as S  # noqa: E402
from insar_wetlands.cube import cube_root  # noqa: E402
from insar_wetlands.paths import hub_root, resolve_drive  # noqa: E402
from insar_wetlands.scoring import assert_preregistered  # noqa: E402

HUB = hub_root()
CUBE = cube_root()
FT = HUB / "08_deliverables" / "insar_full_test"
OUT = FT / "slc_vh"
PREREG = FT / "preregistration" / "X-105.yaml"
R_MIN, N_MIN = 0.3, 15
K = 1.5


def files_of(folder: Path) -> dict[str, Path]:
    return {Q.DATE_RE.search(f.name).group(1): f for f in sorted(folder.glob("*.h5"))
            if "static" not in f.name.lower() and Q.DATE_RE.search(f.name)}


def slope_r(x: np.ndarray, y: np.ndarray) -> tuple[float, float, int]:
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if len(x) < 3 or not np.any(x):
        return np.nan, np.nan, int(len(x))
    return float(np.sum(x * y) / np.sum(x * x)), float(np.corrcoef(x, y)[0, 1]), int(len(x))


def main() -> None:
    assert_preregistered(PREREG, [OUT / "p6.csv"])
    OUT.mkdir(parents=True, exist_ok=True)
    st = xr.open_dataset(CUBE / "silver" / "static_40m.nc")
    gx, gy = st.x.values, st.y.values
    shape = (len(gy), len(gx))
    hard = st.hard_target.values.astype(bool)
    mat = st.zone.values == 1
    px = pd.read_csv(HUB / "08_deliverables" / "field_first" / "plot_pixels.csv").set_index("plot")
    r6, c6 = int(px.loc["P6", "row"]), int(px.loc["P6", "col"])
    plot_mask = np.zeros(shape, bool)
    for r, c in zip(px.row, px.col):
        plot_mask[int(r), int(c)] = True
    cos = {"ascending": float(np.cos(np.radians(np.nanmedian(st.incidence_asc_deg.values)))),
           "descending": float(np.cos(np.radians(np.nanmedian(st.incidence_des_deg.values))))}
    gold = pd.read_csv(CUBE / "gold" / "plots_pair.csv.gz", low_memory=False,
                       usecols=["track", "grid", "plot", "pair", "d_laser_los_mm", "d_wtd_P6_cm"])
    gold = gold[(gold["plot"] == "P6") & (gold.grid == "40m")].set_index(["track", "pair"])
    floor32 = S.noise_floor(32)
    p6_rows, pair_rows, cls_rows = [], [], []
    for track in X.TRACKS:
        fv = files_of(resolve_drive() / "slc_pilot" / "2021" / "cslc" / track)
        fh = files_of(resolve_drive() / "slc_pilot" / "2021" / "cslc_vh" / track)
        dates = sorted(set(fv) & set(fh))
        cls = R.regime(track, K)
        cdates = {d.strftime("%Y%m%d"): i for i, d in enumerate(pd.DatetimeIndex(cls.date.values))}
        carr = cls.values
        prev = None
        for a, b in zip(dates[:-1], dates[1:]):
            if prev is None:
                prev = (Q.cslc(fv[a], "VV"), Q.cslc(fh[a], "VH"))
            cur = (Q.cslc(fv[b], "VV"), Q.cslc(fh[b], "VH"))
            out = {}
            for k, pol in ((0, "VV"), (1, "VH")):
                (s1, x, y), (s2, _, _) = prev[k], cur[k]
                rr, cc = Q.axis_index(x, y, gx, gy)
                z, coh = Q.look(s1, s2, rr, cc, shape)
                ref = np.nanmean(z[hard] / np.maximum(np.abs(z[hard]), 1e-12))
                out[pol] = (z * np.conj(ref), coh)
            key = f"{a}_{b}"
            g = gold.loc[(track, key)] if (track, key) in gold.index else None
            pair_rows.append({"track": track, "t1": a, "t2": b,
                              "vv_mm": float(np.angle(out["VV"][0][r6, c6])) * X.PHASE_TO_MM,
                              "vh_mm": float(np.angle(out["VH"][0][r6, c6])) * X.PHASE_TO_MM,
                              "laser_mm": float(g.d_laser_los_mm) if g is not None else np.nan,
                              "wtd_cm": float(g.d_wtd_P6_cm) if g is not None else np.nan,
                              "vv_coh_mat": float(np.nanmedian(out["VV"][1][mat])),
                              "vh_coh_mat": float(np.nanmedian(out["VH"][1][mat])),
                              "vv_above_floor_mat": float(np.nanmean(out["VV"][1][mat] > floor32)),
                              "vh_above_floor_mat": float(np.nanmean(out["VH"][1][mat] > floor32))})
            ia, ib = cdates.get(a), cdates.get(b)
            if ia is not None and ib is not None:
                ca, cb = carr[ia], carr[ib]
                same = np.where(ca == cb, ca, np.nan)
                dphi = np.angle(out["VV"][0] * np.conj(out["VH"][0]))
                for code, name in R.CODE.items():
                    for where, m in (("mat", mat), ("plot cells", plot_mask)):
                        sel = m & (same == code)
                        if sel.any():
                            cls_rows.append({"track": track, "where": where, "regime": name, "cells": int(sel.sum()),
                                             "vv_coh": float(np.nanmedian(out["VV"][1][sel])),
                                             "vh_coh": float(np.nanmedian(out["VH"][1][sel])),
                                             "abs_vv_minus_vh_rad": float(np.nanmedian(np.abs(dphi[sel])))})
            prev = cur
        p = pd.DataFrame([r for r in pair_rows if r["track"] == track])
        for pol in ("vv", "vh"):
            sl_l, r_l, n_l = slope_r(p.laser_mm.to_numpy(), p[f"{pol}_mm"].to_numpy())
            sl_w, r_w, n_w = slope_r(p.wtd_cm.to_numpy(), p[f"{pol}_mm"].to_numpy())
            p6_rows.append({"track": track, "pol": pol.upper(), "pairs_laser": n_l, "r_laser": r_l, "slope_laser": sl_l,
                            "pairs_wtd": n_w, "r_wtd": r_w, "slope_wtd_mm_per_cm": sl_w,
                            "water_level_expectation_mm_per_cm": 10 * cos[track],
                            "reading": reading(r_l, r_w, sl_w, min(n_l, n_w), cos[track])})
    p6, pairs = pd.DataFrame(p6_rows), pd.DataFrame(pair_rows)
    by_cls = pd.DataFrame(cls_rows).groupby(["track", "where", "regime"]).agg(
        pair_cells=("cells", "sum"), vv_coh=("vv_coh", "median"), vh_coh=("vh_coh", "median"),
        abs_vv_minus_vh_rad=("abs_vv_minus_vh_rad", "median")).reset_index() if cls_rows else pd.DataFrame()
    for name, df in (("p6", p6), ("pairs", pairs), ("by_regime", by_cls)):
        df.to_csv(OUT / f"{name}.csv", index=False)
    (OUT / "README.md").write_text(readme(p6, pairs, by_cls), encoding="utf-8")
    print(p6.to_string(index=False))
    print(f"wrote {OUT}")


def reading(r_l, r_w, sl_w, n, cosi) -> str:
    if n < N_MIN:
        return "too few"
    if np.isfinite(r_l) and r_l >= R_MIN and (not np.isfinite(r_w) or r_l > r_w):
        return "follows the laser"
    if np.isfinite(r_w) and r_w >= R_MIN and r_w > (r_l if np.isfinite(r_l) else -1):
        return "follows the water table" if abs(sl_w - 10 * cosi) < abs(sl_w - 1.25 * cosi) else "neither"
    return "neither"


def _md(df: pd.DataFrame, nd: int = 2) -> str:
    d = df.copy()
    for c in d.columns:
        if pd.api.types.is_float_dtype(d[c]):
            d[c] = d[c].map(lambda x: "—" if pd.isna(x) else f"{x:.{nd}f}")
    return "\n".join(["| " + " | ".join(map(str, d.columns)) + " |", "|" + "---|" * len(d.columns)] +
                     ["| " + " | ".join(map(str, r)) + " |" for r in d.itertuples(index=False)])


def readme(p6, pairs, by_cls) -> str:
    mat = pairs.groupby("track")[["vv_coh_mat", "vh_coh_mat", "vv_above_floor_mat", "vh_above_floor_mat"]].median().reset_index()
    return f"""# X-105 — VH beside VV from our 2021 CSLCs (generated)

Generated by `experiments/slc_vh_x105.py` (`make run E=slc_vh_x105`); do not edit. Pre-registration:
`../preregistration/X-105.yaml`. 40 m cells, 32 looks, unfiltered, each polarisation referenced to its own hard-target
mean, cycle 0. Caveats: X-101 QC v2 phase residual; ascending regime classes come from the HyP3 RTC whose ascending grid is
offset by about a cell (X-109). Tried on P6, 2021.

## At P6: increments against the laser and the water table

{_md(p6)}

## Over the mat, median over consecutive pairs (coherence, share above the 32-look noise floor)

{_md(mat)}

## By regime class (pairs whose two dates share the class): coherence and |VV − VH| phase

{_md(by_cls) if len(by_cls) else "—"}

Files: `p6.csv`, `pairs.csv` (every pair), `by_regime.csv`.
"""


if __name__ == "__main__":
    main()

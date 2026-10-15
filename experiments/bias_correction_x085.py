"""X-085 — the Maghsoudi short-pair bias correction on unfiltered 2021 SLC pairs at P6, scored against the laser.

1-, 2- and 3-connection pairs of each track's 2021 CSLC stack (X-101), the complex 40 m cell sum (unfiltered) at P6 and
the eight cells around, referenced to the hard-target cells; 12–6 and 18–6 closures; the 1-connection biases by least
squares (``insar_wetlands.phase_bias``); the chain of corrected vs uncorrected 1-connection pairs scored by the common
scorer. Rule: ``08_deliverables/insar_full_test/preregistration/X-085.yaml``. Reads the CSLCs from the Drive root
(``INSAR_DRIVE_ROOT``; on the Mac the mirror) and the cube; writes only to the hub:
``08_deliverables/insar_full_test/bias_correction/``.

    INSAR_DRIVE_ROOT=<hub>/05_code/local/drive_mirror make run E=bias_correction_x085
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import p6_validation_x082 as X  # noqa: E402  (laser in LOS, PHASE_TO_MM)
import slc_qc_x101 as Q  # noqa: E402  (CSLC reading, cell looks)

from insar_wetlands import phase_bias  # noqa: E402
from insar_wetlands.cube import cube_root  # noqa: E402
from insar_wetlands.paths import hub_root, resolve_drive  # noqa: E402
from insar_wetlands.scoring import Protocol, assert_preregistered, score  # noqa: E402

HUB = hub_root()
FT = HUB / "08_deliverables" / "insar_full_test"
OUT = FT / "bias_correction"
PREREG = FT / "preregistration" / "X-085.yaml"
A = {"paper (0.47, 0.31)": (0.47, 0.31), "low (0.35, 0.20)": (0.35, 0.20), "high (0.60, 0.45)": (0.60, 0.45)}
OVERPASS = {"ascending": "16:36", "descending": "05:09"}


def unit_mean(z: np.ndarray) -> complex:
    z = z[np.isfinite(z)]
    return np.mean(z / np.maximum(np.abs(z), 1e-12)) if len(z) else np.nan


def main() -> None:
    assert_preregistered(PREREG, [OUT / "verdict.csv"])
    OUT.mkdir(parents=True, exist_ok=True)
    prot = Protocol.load(FT / "protocol.yaml")
    st = xr.open_dataset(cube_root() / "silver" / "static_40m.nc")
    gx, gy = st.x.values, st.y.values
    shape = (len(gy), len(gx))
    hard, mat = st.hard_target.values.astype(bool), st.zone.values == 1
    px = pd.read_csv(HUB / "08_deliverables" / "field_first" / "plot_pixels.csv").set_index("plot")
    r6, c6 = int(px.loc["P6", "row"]), int(px.loc["P6", "col"])
    cells = [(r6 + dr, c6 + dc) for dr in (-1, 0, 1) for dc in (-1, 0, 1)]
    las = X.laser()[0]
    lev_rows, clos_rows, gate_rows, pair_rows = [], [], [], []
    for track in X.TRACKS:
        folder = resolve_drive() / "slc_pilot" / "2021" / "cslc" / track
        files = {Q.DATE_RE.search(f.name).group(1): f for f in sorted(folder.glob("*.h5"))
                 if "static" not in f.name.lower() and Q.DATE_RE.search(f.name)}
        dates = sorted(files)
        cache: dict = {}

        def get(d):
            if d not in cache:
                cache[d] = Q.cslc(files[d])
            return cache[d]

        z = {1: [], 2: [], 3: []}
        coh18 = []
        for k in (1, 2, 3):
            for i in range(len(dates) - k):
                (s1, x, y), (s2, _, _) = get(dates[i]), get(dates[i + k])
                rows, cols = Q.axis_index(x, y, gx, gy)
                zz, cc = Q.look(s1, s2, rows, cols, shape)
                ref = unit_mean(zz[hard])
                z[k].append([zz[r, c] * np.conj(ref) for r, c in cells])
                if k == 3:
                    coh18.append(cc)
        z = {k: np.array(v) for k, v in z.items()}                     # (pairs, 9 cells)
        c18 = np.nanmean(np.stack(coh18), 0)
        gate_rows.append({"track": track, "mat_share_coh18_above_0.3": float(np.nanmean(c18[mat] > 0.3)),
                          "p6_coh18": float(c18[r6, c6]), "pairs_18d": len(coh18)})
        t = pd.DatetimeIndex([pd.Timestamp(f"{d[:4]}-{d[4:6]}-{d[6:]} {OVERPASS[track]}", tz="UTC") for d in dates])
        base = pd.DataFrame({"t1": t[:-1], "t2": t[1:]})
        lt = las[track]
        base["truth_inc_mm"] = (lt.reindex(base.t2.dt.round("h")).to_numpy() - lt.reindex(base.t1.dt.round("h")).to_numpy())
        for ci, (r, c) in enumerate(cells):
            z1, z2, z3 = z[1][:, ci], z[2][:, ci], z[3][:, ci]
            p1 = np.angle(z1)
            c2, c3 = phase_bias.closures(z1, z2, z3)
            variants = {"uncorrected": p1}
            for name, (a1, a2) in A.items():
                d = phase_bias.solve(c2, c3, a1, a2)
                variants[f"corrected, {name}"] = p1 - d
                if ci == 4 and name.startswith("paper"):
                    res2 = c2 - (a1 - 1) * (d[:-1] + d[1:])
                    res3 = c3 - (a2 - 1) * (d[:-2] + d[1:-1] + d[2:])
                    clos_rows.append({"track": track, "closure": "12–6", "cumulative_before_rad": float(np.nansum(c2)),
                                      "cumulative_after_rad": float(np.nansum(res2)), "sd_before_rad": float(np.nanstd(c2)),
                                      "sd_after_rad": float(np.nanstd(res2))})
                    clos_rows.append({"track": track, "closure": "18–6", "cumulative_before_rad": float(np.nansum(c3)),
                                      "cumulative_after_rad": float(np.nansum(res3)), "sd_before_rad": float(np.nanstd(c3)),
                                      "sd_after_rad": float(np.nanstd(res3))})
                    pair_rows += [{"track": track, "t1": base.t1[i], "t2": base.t2[i], "phase_rad": p1[i],
                                   "bias_rad": d[i], "truth_inc_mm": base.truth_inc_mm[i]} for i in range(len(p1))]
            for vname, ph in variants.items():
                q = base.assign(est_mm=ph * X.PHASE_TO_MM)
                sc = score(q, lt, prot)
                for _, lv in sc["levels"].iterrows():
                    if lv.view in ("datum A", "full span"):
                        lev_rows.append({"track": track, "cell": ci, "row": r, "col": c, "variant": vname,
                                         "view": lv.view, "epochs": lv.epochs, "pearson_detrended": lv.pearson_detrended,
                                         "rmse_detrended_mm": lv.rmse_detrended_mm})
    lev, clos, gate, pairs = map(pd.DataFrame, (lev_rows, clos_rows, gate_rows, pair_rows))
    v = verdict(lev)
    for name, df in (("levels", lev), ("closures", clos), ("gate", gate), ("pairs_p6", pairs), ("verdict", v)):
        df.to_csv(OUT / f"{name}.csv", index=False)
    (OUT / "README.md").write_text(readme(lev, clos, gate, v), encoding="utf-8")
    print(gate.to_string(index=False))
    print(v.to_string(index=False))
    print(f"wrote {OUT}")


def verdict(lev: pd.DataFrame) -> pd.DataFrame:
    p6 = lev[(lev.cell == 4) & (lev.view == "datum A")].set_index(["variant", "track"]).rmse_detrended_mm
    rows = []
    for vname in sorted(set(lev.variant) - {"uncorrected"}):
        ch = {t: p6.loc[(vname, t)] / p6.loc[("uncorrected", t)] - 1 for t in X.TRACKS}
        rows.append({"variant": vname, **{f"rmse_{t[:4]}_mm": p6.loc[(vname, t)] for t in X.TRACKS},
                     **{f"uncorrected_{t[:4]}_mm": p6.loc[("uncorrected", t)] for t in X.TRACKS},
                     **{f"change_{t[:4]}": ch[t] for t in X.TRACKS},
                     "helps": bool(all(c <= -0.20 for c in ch.values())),
                     "harms": bool(any(c >= 0.20 for c in ch.values()))})
    return pd.DataFrame(rows)


def _md(df: pd.DataFrame, nd: int = 2) -> str:
    d = df.copy()
    for col in d.columns:
        if pd.api.types.is_float_dtype(d[col]):
            d[col] = d[col].map(lambda x: "—" if pd.isna(x) else f"{x:.{nd}f}")
    return "\n".join(["| " + " | ".join(map(str, d.columns)) + " |", "|" + "---|" * len(d.columns)] +
                     ["| " + " | ".join(map(str, r)) + " |" for r in d.itertuples(index=False)])


def readme(lev, clos, gate, v) -> str:
    p6 = lev[(lev.cell == 4)][["track", "variant", "view", "epochs", "pearson_detrended", "rmse_detrended_mm"]]
    nine = (lev[lev.view == "datum A"].groupby(["track", "variant"]).rmse_detrended_mm.median()
            .reset_index(name="median_rmse_nine_cells_mm"))
    return f"""# X-085 — the Maghsoudi short-pair bias correction on unfiltered SLC pairs (generated)

Generated by `experiments/bias_correction_x085.py` (`make run E=bias_correction_x085`); do not edit. Pre-registration:
`../preregistration/X-085.yaml`. 2021, both tracks; 1-, 2-, 3-connection pairs of the X-101 CSLC stack, complex 40 m
cell sum (unfiltered), hard-target reference; biases by least squares from the 12–6 and 18–6 closures
(maghsoudi2022, a1 = 0.47, a2 = 0.31 as published, two sensitivity pairs). Epochs at each track's overpass hour.
Caveat of X-101 (processing differs from HyP3). One year, one datum: the winter loss is not tested here.

## Gate the paper states (recorded, not decisive)

{_md(gate)}

## Verdict (datum A, P6 cell)

{_md(v)}

## Closures at P6 before and after (paper constants)

{_md(clos)}

## Levels at P6

{_md(p6)}

## Median over the nine cells (datum A)

{_md(nine)}
"""


if __name__ == "__main__":
    main()

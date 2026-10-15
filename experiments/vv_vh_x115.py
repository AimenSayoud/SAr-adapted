"""X-115 — VV and VH 6-day phase combined into one estimate of the mat's motion at P6 (2021 CSLCs).

Per consecutive pair, the 40 m cell sum at P6 and the eight cells around (unfiltered), each polarisation referenced to
its own hard-target mean; weights 1/σ² from each polarisation's cell coherence (Cramér–Rao, 32 looks); the weighted
complex sum. Chains scored by the common scorer. Rule: ``08_deliverables/insar_full_test/preregistration/X-115.yaml``.
Reads the CSLCs (``INSAR_DRIVE_ROOT``: ``slc_pilot/2021/cslc`` and ``cslc_vh``) and the cube; writes only to the hub:
``08_deliverables/insar_full_test/vv_vh/``.

    INSAR_DRIVE_ROOT=<hub>/05_code/local/drive_mirror make run E=vv_vh_x115
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import bias_correction_x085 as BC  # noqa: E402  (overpass hours, unit_mean)
import p6_validation_x082 as X  # noqa: E402
import slc_qc_x101 as Q  # noqa: E402

from insar_wetlands.cube import cube_root  # noqa: E402
from insar_wetlands.paths import hub_root, resolve_drive  # noqa: E402
from insar_wetlands.scoring import Protocol, assert_preregistered, score  # noqa: E402

HUB = hub_root()
FT = HUB / "08_deliverables" / "insar_full_test"
OUT = FT / "vv_vh"
PREREG = FT / "preregistration" / "X-115.yaml"
LOOKS = 32


def files(folder: Path) -> dict[str, Path]:
    return {Q.DATE_RE.search(f.name).group(1): f for f in sorted(folder.glob("*.h5"))
            if "static" not in f.name.lower() and Q.DATE_RE.search(f.name)}


def main() -> None:
    assert_preregistered(PREREG, [OUT / "verdict.csv"])
    OUT.mkdir(parents=True, exist_ok=True)
    prot = Protocol.load(FT / "protocol.yaml")
    st = xr.open_dataset(cube_root() / "silver" / "static_40m.nc")
    gx, gy = st.x.values, st.y.values
    shape = (len(gy), len(gx))
    hard = st.hard_target.values.astype(bool)
    px = pd.read_csv(HUB / "08_deliverables" / "field_first" / "plot_pixels.csv").set_index("plot")
    r6, c6 = int(px.loc["P6", "row"]), int(px.loc["P6", "col"])
    cells = [(r6 + dr, c6 + dc) for dr in (-1, 0, 1) for dc in (-1, 0, 1)]
    las = X.laser()[0]
    lev, pair_rows = [], []
    for track in X.TRACKS:
        fv = files(resolve_drive() / "slc_pilot" / "2021" / "cslc" / track)
        fh = files(resolve_drive() / "slc_pilot" / "2021" / "cslc_vh" / track)
        dates = sorted(set(fv) & set(fh))
        ph = {"VV": [], "VH": []}
        co = {"VV": [], "VH": []}
        for a, b in zip(dates[:-1], dates[1:]):
            for pol, ff in (("VV", fv), ("VH", fh)):
                s1, x, y = Q.cslc(ff[a], pol)
                s2, _, _ = Q.cslc(ff[b], pol)
                rr, cc = Q.axis_index(x, y, gx, gy)
                z, coh = Q.look(s1, s2, rr, cc, shape)
                ref = BC.unit_mean(z[hard])
                ph[pol].append([z[r, c] * np.conj(ref) for r, c in cells])
                co[pol].append([coh[r, c] for r, c in cells])
        zv, zh = np.array(ph["VV"]), np.array(ph["VH"])
        wv = 1 / Q.crb(np.array(co["VV"], float), LOOKS) ** 2
        wh = 1 / Q.crb(np.array(co["VH"], float), LOOKS) ** 2
        est = {"VV": np.angle(zv), "VH": np.angle(zh),
               "VV+VH": np.angle(wv * zv / np.abs(zv) + wh * zh / np.abs(zh))}
        t = pd.DatetimeIndex([pd.Timestamp(f"{d[:4]}-{d[4:6]}-{d[6:]} {BC.OVERPASS[track]}", tz="UTC") for d in dates])
        base = pd.DataFrame({"t1": t[:-1], "t2": t[1:]})
        lt = las[track]
        base["truth_inc_mm"] = lt.reindex(base.t2.dt.round("h")).to_numpy() - lt.reindex(base.t1.dt.round("h")).to_numpy()
        pair_rows += [{"track": track, "t1": base.t1[i], "w_vv": float(wv[i, 4]), "w_vh": float(wh[i, 4])}
                      for i in range(len(base))]
        for ci, (r, c) in enumerate(cells):
            for name, e in est.items():
                sc = score(base.assign(est_mm=e[:, ci] * X.PHASE_TO_MM), lt, prot)
                inc = sc["increments"].iloc[0].to_dict() if len(sc.get("increments", [])) else {}
                for _, lv in sc["levels"].iterrows():
                    if lv.view == "datum A":
                        lev.append({"track": track, "cell": ci, "estimate": name, "epochs": lv.epochs,
                                    "pearson_detrended": lv.pearson_detrended, "rmse_detrended_mm": lv.rmse_detrended_mm,
                                    **{f"inc_{k}": v for k, v in inc.items() if k in ("pearson_increments", "rmse_increments_mm")}})
    lev, pairs = pd.DataFrame(lev), pd.DataFrame(pair_rows)
    p6 = lev[lev.cell == 4].set_index(["estimate", "track"]).rmse_detrended_mm
    rows = []
    for tr in X.TRACKS:
        rows.append({"track": tr, "rmse_vv_mm": p6.loc[("VV", tr)], "rmse_vh_mm": p6.loc[("VH", tr)],
                     "rmse_combined_mm": p6.loc[("VV+VH", tr)],
                     "change_vs_vv": p6.loc[("VV+VH", tr)] / p6.loc[("VV", tr)] - 1,
                     "share_pairs_vh_weight_above_vv": float((pairs[pairs.track == tr].w_vh > pairs[pairs.track == tr].w_vv).mean())})
    v = pd.DataFrame(rows)
    v["helps_both_tracks"] = bool((v.change_vs_vv <= -0.15).all())
    nine = lev.groupby(["track", "estimate"]).rmse_detrended_mm.median().reset_index(name="median_rmse_nine_cells_mm")
    for name, df in (("levels", lev), ("pairs", pairs), ("verdict", v), ("nine_cells", nine)):
        df.to_csv(OUT / f"{name}.csv", index=False)
    (OUT / "README.md").write_text(f"""# X-115 — VV and VH combined at P6 (generated)

Generated by `experiments/vv_vh_x115.py`; do not edit. Pre-registration: `../preregistration/X-115.yaml`. 2021 CSLCs,
40 m cell, unfiltered, each polarisation referenced to its own hard-target mean, inverse-variance (Cramér–Rao) weights;
datum A. Caveat of X-101 QC v2.

## Verdict (P6)

{BC._md(v)}

## P6 cell, all estimates

{BC._md(lev[lev.cell == 4].drop(columns='cell'))}

## Median over the nine cells

{BC._md(nine)}
""", encoding="utf-8")
    print(v.to_string(index=False))


if __name__ == "__main__":
    main()

"""X-124 — is the HyP3 half-cell gain placement or interpolation smoothing?

(a) X-120 refine's shift search on HyP3 coherence pre-smoothed by a 3×3 cell mean, and without (must reproduce refine);
(b) the plain P6 chain (hard-target reference, X-118's scorer) with HyP3's unwrapped phase interpolated at the corrected
position and at its three mirrors, which carry the same interpolation weights. Both on the cube as delivered (the
2026-10-15 snapshot, which X-118 and X-120 measured), set here as the cube root. Rule: ``preregistration/X-124.yaml``.

    INSAR_DRIVE_ROOT=<hub>/05_code/local/drive_mirror make run E=hyp3_smoothing_check_x124
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

from insar_wetlands.paths import hub_root

HUB = hub_root()
SNAPSHOT = HUB / "06_data" / "cube_snapshots" / "as_delivered_2026-10-15"
os.environ["RZECIN_CUBE_ROOT"] = str(SNAPSHOT)          # before the modules below resolve the cube

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import xarray as xr  # noqa: E402
from scipy.ndimage import uniform_filter  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import half_cell_p6_x118 as H  # noqa: E402
import hyp3_half_cell_refine_x120 as R  # noqa: E402
import p6_validation_x082 as X  # noqa: E402
import slc_qc_x101 as Q  # noqa: E402

from insar_wetlands.cube import cube_root  # noqa: E402
from insar_wetlands.paths import resolve_drive  # noqa: E402
from insar_wetlands.scoring import Protocol, assert_preregistered, score  # noqa: E402

FT = HUB / "08_deliverables" / "insar_full_test"
OUT = FT / "hyp3_smoothing_check"
PREREG = FT / "preregistration" / "X-124.yaml"


def smooth3(a: np.ndarray) -> np.ndarray:
    """3×3 mean, NaN-aware (a cell is NaN only where all nine are)."""
    ok = np.isfinite(a)
    num = uniform_filter(np.where(ok, a, 0.0), 3, mode="nearest")
    den = uniform_filter(ok.astype(float), 3, mode="nearest")
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > 0, num / den, np.nan)


def part_a() -> pd.DataFrame:
    st = xr.open_dataset(cube_root() / "silver" / "static_40m.nc")
    gx, gy = st.x.values, st.y.values
    rows = []
    for track in ("ascending", "descending"):
        folder = resolve_drive() / "slc_pilot" / "2021" / "cslc" / track
        files = {Q.DATE_RE.search(f.name).group(1): f for f in sorted(folder.glob("*.h5"))
                 if "static" not in f.name.lower() and Q.DATE_RE.search(f.name)}
        dates = sorted(files)
        hy = Q.hyp3_pairs(track)
        for a, b in zip(dates[:-1], dates[1:]):
            if (a, b) not in hy:
                continue
            s1, x, y = Q.cslc(files[a])
            s2, _, _ = Q.cslc(files[b])
            rr, cc = Q.axis_index(x, y, gx, gy)
            _, coh = Q.look(s1, s2, rr, cc, (len(gy), len(gx)))
            h = np.asarray(hy[(a, b)][1], float)
            for name, img in (("as delivered", h), ("pre-smoothed 3x3", smooth3(h))):
                dx, dy, r, _ = R.best(img, coh, gx, gy)
                rows.append({"track": track, "t1": a, "t2": b, "hyp3": name, "dx_m": dx, "dy_m": dy, "r_best": r})
            print(track, a, b, rows[-2]["dx_m"], rows[-2]["dy_m"], rows[-1]["dx_m"], rows[-1]["dy_m"], flush=True)
    return pd.DataFrame(rows)


def part_b() -> pd.DataFrame:
    prot = Protocol.load(FT / "protocol.yaml")
    shift = pd.read_csv(FT / "hyp3_placement" / "decisive" / "verdict.csv").set_index("track")
    px = pd.read_csv(HUB / "08_deliverables" / "field_first" / "plot_pixels.csv").set_index("plot")
    xt, yt = float(px.loc["P6", "E"]), float(px.loc["P6", "N"])
    p, las = X.load_pairs(), X.laser()[0]
    g = pd.read_csv(cube_root() / "gold" / "plots_pair.csv.gz", low_memory=False,
                    usecols=["track", "grid", "plot", "pair", X.REFS["hard_targets"]])
    g = g[(g["plot"] == "P6") & (g.grid == "40m")].set_index(["track", "pair"])
    rows = []
    for track in X.TRACKS:
        dx, dy = float(shift.loc[track, "median_dx_m"]), float(shift.loc[track, "median_dy_m"])
        lt = las[track].dropna()
        q0 = p[(p.track == track) & (p.t1 >= lt.index.min() - pd.Timedelta(days=13)) &
               (p.t2 <= lt.index.max() + pd.Timedelta(days=13))].sort_values("t1").copy()
        keys = q0.t1.dt.strftime("%Y%m%d") + "_" + q0.t2.dt.strftime("%Y%m%d")
        refv = np.asarray([g.loc[(track, k), X.REFS["hard_targets"]] if (track, k) in g.index else np.nan
                           for k in keys], float)
        reads = {"corrected (laser − shift)": (xt - dx, yt - dy), "mirror (laser + shift)": (xt + dx, yt + dy),
                 "mirror east–west only": (xt + dx, yt - dy), "mirror north–south only": (xt - dx, yt + dy)}
        for name, (xq, yq) in reads.items():
            u = keys.map(H.unw_at(track, xq, yq))
            q = q0.assign(est_mm=(u.to_numpy(float) - refv) * X.PHASE_TO_MM).dropna(subset=["est_mm"])
            sc = score(q.assign(truth_inc_mm=q.d_laser_los_mm), lt, prot)
            for _, lv in sc["levels"].iterrows():
                if lv.view in ("datum A", "datum B"):
                    rows.append({"track": track, "read": name, "read_x": xq, "read_y": yq, "view": lv.view,
                                 "epochs": lv.epochs, "rmse_detrended_mm": lv.rmse_detrended_mm})
        lab = pd.read_csv(FT / "half_cell_p6" / "levels.csv")
        lab = lab[(lab.track == track) & (lab.reference == "hard_targets") & (lab["read"] == "labelled cell")
                  & lab.view.isin(["datum A", "datum B"])]
        for _, lv in lab.iterrows():
            rows.append({"track": track, "read": "labelled cell (X-118)", "read_x": np.nan, "read_y": np.nan,
                         "view": lv.view, "epochs": lv.epochs, "rmse_detrended_mm": lv.rmse_detrended_mm})
    return pd.DataFrame(rows)


def main() -> None:
    assert_preregistered(PREREG, [OUT / "verdict.csv"])
    OUT.mkdir(parents=True, exist_ok=True)
    A = part_a()
    A.to_csv(OUT / "pairs.csv", index=False)
    ref = pd.read_csv(FT / "hyp3_half_cell_confirm" / "refine" / "verdict.csv").set_index("track")
    va = (A.groupby(["track", "hyp3"]).agg(pairs=("t1", "size"), median_dx_m=("dx_m", "median"),
                                          median_dy_m=("dy_m", "median"), r_best_median=("r_best", "median"))
          .reset_index())
    va["refine_dx_m"] = va.track.map(ref.median_dx_m)
    va["refine_dy_m"] = va.track.map(ref.median_dy_m)
    sm = va[va.hyp3 == "pre-smoothed 3x3"]
    a_ok = bool(len(sm) == 2 and (sm.pairs >= 30).all() and ((sm.median_dx_m - sm.refine_dx_m).abs() <= 5).all()
                and ((sm.median_dy_m - sm.refine_dy_m).abs() <= 5).all())
    B = part_b()
    B.to_csv(OUT / "reads.csv", index=False)
    four = B[~B["read"].str.startswith("labelled")]
    best = four.loc[four.groupby(["track", "view"]).rmse_detrended_mm.idxmin()]
    n_best = int((best["read"] == "corrected (laser − shift)").sum())
    b_out = ("placement supported" if n_best >= 3 else "smoothing explains the gain" if n_best <= 1 else "mixed")
    mir = four[four["read"].str.startswith("mirror")].groupby(["track", "view"]).rmse_detrended_mm.mean()
    lab = B[B["read"].str.startswith("labelled")].set_index(["track", "view"]).rmse_detrended_mm
    cor = four[four["read"].str.startswith("corrected")].set_index(["track", "view"]).rmse_detrended_mm
    tb = pd.DataFrame({"labelled_cell": lab, "mirrors_mean": mir, "corrected": cor}).reset_index()
    V = pd.DataFrame([{"part": "a", "statistic": "pre-smoothed median within 5 m of refine on both axes, both tracks",
                       "value": a_ok, "outcome": "shift is placement" if a_ok else "size in doubt"},
                      {"part": "b", "statistic": "views where the corrected read is the lowest of four",
                       "value": n_best, "outcome": b_out}])
    va.to_csv(OUT / "shift_summary.csv", index=False)
    tb.to_csv(OUT / "smoothing_gain.csv", index=False)
    V.to_csv(OUT / "verdict.csv", index=False)
    md = lambda df: "\n".join(["| " + " | ".join(map(str, df.columns)) + " |", "|" + "---|" * len(df.columns)] +  # noqa: E731
                              ["| " + " | ".join(f"{x:.2f}" if isinstance(x, float) else str(x) for x in r) + " |"
                               for r in df.itertuples(index=False)])
    (OUT / "README.md").write_text(f"""# X-124 — the HyP3 half-cell gain: placement or interpolation smoothing? (generated)

Generated by `experiments/hyp3_smoothing_check_x124.py`; do not edit. Rule: `../preregistration/X-124.yaml`. On the
cube as delivered (`06_data/cube_snapshots/as_delivered_2026-10-15`). (a) X-120 refine's shift search with HyP3's
coherence as delivered and pre-smoothed by a 3×3 cell mean. (b) the plain P6 chain, hard-target reference, read at the
corrected position and at its three mirrors (same interpolation weights); detrended level RMSE, datum A and B.

## Verdict

{md(V)}

## (a) Shift per track

{md(va)}

## (b) Reads at P6 (mm)

{md(four[["track", "view", "read", "rmse_detrended_mm"]])}

## (b) The smoothing gain alone: the labelled cell, the mean of the mirrors, the corrected read (mm)

{md(tb)}
""", encoding="utf-8")
    print(V.to_string(index=False))
    print(va.to_string(index=False))
    print(tb.to_string(index=False))


if __name__ == "__main__":
    main()

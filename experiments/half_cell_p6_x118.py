"""X-118 — the P6 chain read at the half-cell-corrected HyP3 position (T75) vs the labelled cell.

Per consecutive pair of X-099, HyP3's unwrapped phase bilinearly interpolated at P6 moved by the track's T75 shift,
referenced as X-099 (hard targets, grassland); both chains scored by the common scorer. Rule:
``08_deliverables/insar_full_test/preregistration/X-118.yaml``. Reads the cube; writes only to the hub:
``08_deliverables/insar_full_test/half_cell_p6/``.

    make run E=half_cell_p6_x118
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import p6_validation_x082 as X  # noqa: E402

from insar_wetlands.cube import cube_root  # noqa: E402
from insar_wetlands.paths import hub_root  # noqa: E402
from insar_wetlands.scoring import Protocol, assert_preregistered, score  # noqa: E402

HUB = hub_root()
FT = HUB / "08_deliverables" / "insar_full_test"
OUT = FT / "half_cell_p6"
PREREG = FT / "preregistration" / "X-118.yaml"
SETS = ("2020_2021", "2022_2024", "2025")


def unw_at(track: str, xq: float, yq: float) -> dict[str, float]:
    """HyP3 unwrapped phase (rad) bilinearly interpolated at (xq, yq), per pair key."""
    out = {}
    for s in SETS:
        f = cube_root() / "silver" / f"ifg_{s}_{track}.nc"
        if not f.exists():
            continue
        with xr.open_dataset(f) as d:
            v = d.unw.interp(x=xq, y=yq, method="linear").values
            t1 = pd.to_datetime(d.t1.values).strftime("%Y%m%d")
            t2 = pd.to_datetime(d.t2.values).strftime("%Y%m%d")
            for k, (a, b) in enumerate(zip(t1, t2)):
                out[f"{a}_{b}"] = float(v[k])
    return out


def main() -> None:
    assert_preregistered(PREREG, [OUT / "verdict.csv"])
    OUT.mkdir(parents=True, exist_ok=True)
    prot = Protocol.load(FT / "protocol.yaml")
    shift = pd.read_csv(FT / "hyp3_placement" / "decisive" / "verdict.csv").set_index("track")
    px = pd.read_csv(HUB / "08_deliverables" / "field_first" / "plot_pixels.csv").set_index("plot")
    st = xr.open_dataset(cube_root() / "silver" / "static_40m.nc")
    r6, c6 = int(px.loc["P6", "row"]), int(px.loc["P6", "col"])
    x6, y6 = float(st.x.values[c6]), float(st.y.values[r6])          # the labelled cell's centre, as X-099 reads it
    xt, yt = float(px.loc["P6", "E"]), float(px.loc["P6", "N"])       # the laser's surveyed position
    p, las = X.load_pairs(), X.laser()[0]
    g = pd.read_csv(cube_root() / "gold" / "plots_pair.csv.gz", low_memory=False,
                    usecols=["track", "grid", "plot", "pair", *X.REFS.values()])
    g = g[(g["plot"] == "P6") & (g.grid == "40m")].set_index(["track", "pair"])
    rows, steps, pos = [], [], []
    for track in X.TRACKS:
        dx, dy = float(shift.loc[track, "median_dx_m"]), float(shift.loc[track, "median_dy_m"])
        xq, yq = xt - dx, yt - dy
        pos.append({"track": track, "label_x": x6, "label_y": y6, "laser_x": xt, "laser_y": yt,
                    "read_x": xq, "read_y": yq, "t75_dx_m": dx, "t75_dy_m": dy})
        corr = unw_at(track, xq, yq)
        lt = las[track].dropna()
        q = p[(p.track == track) & (p.t1 >= lt.index.min() - pd.Timedelta(days=13)) &
              (p.t2 <= lt.index.max() + pd.Timedelta(days=13))].sort_values("t1").copy()
        keys = q.t1.dt.strftime("%Y%m%d") + "_" + q.t2.dt.strftime("%Y%m%d")
        q["unw_corr"] = keys.map(corr)
        for ref, col in X.REFS.items():
            refv = [g.loc[(track, k), col] if (track, k) in g.index else np.nan for k in keys]
            q[f"c_{ref}"] = (q.unw_corr - np.asarray(refv, float)) * X.PHASE_TO_MM
            for name, est in (("labelled cell", f"s1_{ref}"), ("corrected position", f"c_{ref}")):
                qq = q.dropna(subset=[est])
                sc = score(qq.assign(est_mm=qq[est], truth_inc_mm=qq.d_laser_los_mm), lt, prot)
                for _, lv in sc["levels"].iterrows():
                    if lv.view in ("full span", "datum A", "datum B"):
                        rows.append({"track": track, "reference": ref, "read": name, "view": lv.view, "epochs": lv.epochs,
                                     "pearson_detrended": lv.pearson_detrended, "rmse_detrended_mm": lv.rmse_detrended_mm})
                w = sc["winters"].assign(track=track, reference=ref, read=name)
                steps.append(w)
    lev = pd.DataFrame(rows)
    win = pd.concat(steps, ignore_index=True)
    h = lev[(lev.reference == "hard_targets") & lev.view.isin(["datum A", "datum B"])]
    piv = h.pivot_table(index=["track", "view"], columns="read", values="rmse_detrended_mm").reset_index()
    piv["change"] = piv["corrected position"] / piv["labelled cell"] - 1
    adopted = bool((piv.change <= -0.10).sum() >= 3 and (piv.change <= 0.10).all())
    v = piv.assign(adopt_corrected=adopted)
    for name, df in (("levels", lev), ("winters", win), ("positions", pd.DataFrame(pos)), ("verdict", v)):
        df.to_csv(OUT / f"{name}.csv", index=False)
    md = lambda df: "\n".join(["| " + " | ".join(map(str, df.columns)) + " |", "|" + "---|" * len(df.columns)] +
                              ["| " + " | ".join(f"{x:.2f}" if isinstance(x, float) else str(x) for x in r) + " |"
                               for r in df.itertuples(index=False)])
    (OUT / "README.md").write_text(f"""# X-118 — the P6 chain at the half-cell-corrected HyP3 position (generated)

Generated by `experiments/half_cell_p6_x118.py`; do not edit. Pre-registration: `../preregistration/X-118.yaml`. HyP3's
unwrapped phase interpolated at the laser's position moved by each track's T75 shift, against the labelled cell; the
plain chain, both references; common scorer.

## Verdict (hard-target reference)

{md(v)}

## Read positions

{md(pd.DataFrame(pos))}

## Levels

{md(lev)}
""", encoding="utf-8")
    print(v.to_string(index=False))


if __name__ == "__main__":
    main()

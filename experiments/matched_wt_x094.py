"""X-094 — do long HyP3 pairs whose dates share the water-table state keep phase on the mat?

SETTLED §3's reopen condition "same-water-table pairs", on the T08 product. HyP3 pairs (D-029) 60–400 days long between
dry, snow-free dates of 2020–2024, matched (|dWTD| <= 2 cm) or unmatched (>= 10 cm), the same number per baseline bin
(``05_code/local/scripts/s1_matched_wt_run.py``; cube set ``matched_wt``). Per pair, the noise floor is HyP3's own: the
95th percentile of coherence over open water in that pair (T58: open water scores the mat's HyP3 coherence). Rule:
``08_deliverables/insar_full_test/preregistration/X-094.yaml``.

Reads the cube and the candidate list (hub); writes only to the hub: ``08_deliverables/insar_full_test/matched_wt/``.

    make run E=matched_wt_x094
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import p6_validation_x082 as X  # noqa: E402  (laser in LOS, PHASE_TO_MM)

from insar_wetlands.cube import cube_root  # noqa: E402
from insar_wetlands.paths import hub_root  # noqa: E402
from insar_wetlands.scoring import QUARTER_WAVE_MM, assert_preregistered  # noqa: E402

HUB = hub_root()
CUBE = cube_root()
FT = HUB / "08_deliverables" / "insar_full_test"
OUT = FT / "matched_wt"
PREREG = FT / "preregistration" / "X-094.yaml"
CAND = HUB / "06_data" / "local_small" / "s1_matched_wt"
WATER, KEEP, MARGIN = 0.9, 0.5, 0.2


def main() -> None:
    assert_preregistered(PREREG, [OUT / "verdict.csv"])
    OUT.mkdir(parents=True, exist_ok=True)
    st = xr.open_dataset(CUBE / "silver" / "static_40m.nc")
    classes = {"mat": st.zone.values == 1, "grassland": st.zone.values == 3,
               "hard_targets": st.hard_target.values.astype(bool)}
    water = st.wc_water.values >= WATER
    hard = classes["hard_targets"]
    px = pd.read_csv(HUB / "08_deliverables" / "field_first" / "plot_pixels.csv").set_index("plot")
    r6, c6 = int(px.loc["P6", "row"]), int(px.loc["P6", "col"])
    las = X.laser()[0]
    rows, p6 = [], []
    for track in X.TRACKS:
        cand = pd.read_csv(CAND / f"candidate_pairs_{track}.csv").set_index("pair")
        f = CUBE / "silver" / f"ifg_matched_wt_{track}.nc"
        if not f.exists():
            raise SystemExit(f"{f.name} missing: run the cube build (harmonise) after the X-094 download")
        d = xr.open_dataset(f)
        lt = las[track].dropna()
        for i, pair in enumerate(map(str, d.pair.values)):
            if pair not in cand.index:
                continue
            c = cand.loc[pair]
            coh = d.coh.isel(pair=i).values
            floor = float(np.nanquantile(coh[water], 0.95))
            rec = {"track": track, "pair": pair, "kind": c.kind, "bin": c["bin"], "days": int(c.dt_days),
                   "dwtd_cm": float(c.dwtd_cm), "floor": floor}
            for k, m in classes.items():
                v = coh[m]
                rec[f"share_above_{k}"] = float(np.nanmean(v > floor))
                rec[f"median_coh_{k}"] = float(np.nanmedian(v))
            rows.append(rec)
            if coh[r6, c6] > floor:
                w = d.wrapped.isel(pair=i).values
                ref = np.angle(np.nanmean(np.exp(1j * w[hard])))
                est = float(np.angle(np.exp(1j * (w[r6, c6] - ref)))) * X.PHASE_TO_MM
                t1, t2 = pd.Timestamp(d.t1.values[i], tz="UTC"), pd.Timestamp(d.t2.values[i], tz="UTC")
                a, b = lt.reindex([t1.round("h")]).iloc[0], lt.reindex([t2.round("h")]).iloc[0]
                p6.append({"track": track, "pair": pair, "kind": c.kind, "est_mm": est,
                           "laser_mm": float(b - a) if np.isfinite(a) and np.isfinite(b) else np.nan})
    t = pd.DataFrame(rows)
    p6 = pd.DataFrame(p6)
    summ = t.groupby(["track", "kind"]).agg(pairs=("pair", "size"), median_floor=("floor", "median"),
                                           **{f"median_share_above_{k}": (f"share_above_{k}", "median") for k in classes}
                                           ).reset_index()
    by_bin = t.groupby(["track", "kind", "bin"]).agg(pairs=("pair", "size"),
                                                    median_share_above_mat=("share_above_mat", "median")).reset_index()
    v = verdict(summ)
    p6s = p6_stats(p6)
    for name, df in (("pairs", t), ("summary", summ), ("by_bin", by_bin), ("p6", p6), ("p6_stats", p6s), ("verdict", v)):
        df.to_csv(OUT / f"{name}.csv", index=False)
    (OUT / "README.md").write_text(readme(summ, by_bin, p6s, v), encoding="utf-8")
    print(v.to_string(index=False))
    print(f"wrote {OUT}")


def verdict(summ: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for track, g in summ.groupby("track"):
        g = g.set_index("kind")
        m = g.loc["matched", "median_share_above_mat"] if "matched" in g.index else np.nan
        u = g.loc["unmatched", "median_share_above_mat"] if "unmatched" in g.index else np.nan
        rows.append({"track": track, "matched_share_mat": m, "unmatched_share_mat": u,
                     "keeps_phase": bool(np.isfinite(m) and np.isfinite(u) and m >= KEEP and m - u >= MARGIN)})
    v = pd.DataFrame(rows)
    v["reopens_settled_3_on_hyp3"] = bool(v.keeps_phase.all()) if len(v) else False
    return v


def p6_stats(p6: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (track, kind), g in p6.dropna(subset=["laser_mm"]).groupby(["track", "kind"]):
        r = float(np.corrcoef(g.est_mm, g.laser_mm)[0, 1]) if len(g) >= 3 else np.nan
        rows.append({"track": track, "kind": kind, "pairs_p6_above_floor": len(g), "pearson": r,
                     "within_quarter_wave": float(np.mean(np.abs(g.est_mm - g.laser_mm) <= QUARTER_WAVE_MM))})
    return pd.DataFrame(rows)


def _md(df: pd.DataFrame, nd: int = 2) -> str:
    d = df.copy()
    for c in d.columns:
        if pd.api.types.is_float_dtype(d[c]):
            d[c] = d[c].map(lambda x: "—" if pd.isna(x) else f"{x:.{nd}f}")
    return "\n".join(["| " + " | ".join(map(str, d.columns)) + " |", "|" + "---|" * len(d.columns)] +
                     ["| " + " | ".join(map(str, r)) + " |" for r in d.itertuples(index=False)])


def readme(summ, by_bin, p6s, v) -> str:
    return f"""# X-094 — long HyP3 pairs matched on the water-table state (generated)

Generated by `experiments/matched_wt_x094.py` (`make run E=matched_wt_x094`); do not edit. Pre-registration:
`../preregistration/X-094.yaml`. HyP3 INSAR_ISCE_BURST 40 m (D-029), 60–400-day pairs between dry snow-free dates of
2020–2024, matched |dWTD| <= 2 cm vs unmatched >= 10 cm, same count per baseline bin. Noise floor per pair = the 95th
percentile of HyP3 coherence over open water in that pair. Matched pairs keep phase if the mat's median share above the
floor is >= {KEEP} and >= {MARGIN} above the unmatched, both tracks.

## Verdict

{_md(v)}

## Median over pairs by track and kind

{_md(summ)}

## By baseline bin (mat)

{_md(by_bin)}

## P6, pairs where its cell clears the floor: phase vs laser change (2021–2024, design data)

{_md(p6s) if len(p6s) else "no pair clears the floor at P6"}
"""


if __name__ == "__main__":
    main()

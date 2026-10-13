"""X-104 — do long pairs keep phase on the mat when coherence is estimated at full resolution with a small window?

SETTLED §3: on the 40 m HyP3 products long pairs sat at the decorrelated level (T08); its reopen conditions include a
full-resolution, small-window coherence (tampuu2020 kept coherence on an open bog that way) and same-water-table pairs.
On the 2021 CSLCs of X-101: long pairs (>= 48 days, snow-free) from the highest-water-table date, and between dates
of matched water table; coherence in 3 × 3 and 5 × 3 boxcars and in X-103's SHP windows; the share of mat pixels above
the estimator's noise floor (``slc_stats.noise_floor``), against grassland and hard targets. Rules:
``08_deliverables/insar_full_test/preregistration/X-104.yaml``.

Reads the CSLCs from the Drive root (on the Mac the mirror) and the cube; writes only to the hub:
``08_deliverables/insar_full_test/slc_long_pairs/``.

    INSAR_DRIVE_ROOT=<hub>/05_code/local/drive_mirror make run E=slc_longpairs_x104
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import p6_validation_x082 as X  # noqa: E402
import slc_ps_shp_x103 as P  # noqa: E402  (stack loader, sample classes, date sets, SHP window)

from insar_wetlands import slc_stats as S  # noqa: E402
from insar_wetlands.cube import cube_root  # noqa: E402
from insar_wetlands.paths import hub_root  # noqa: E402
from insar_wetlands.scoring import Protocol, assert_preregistered  # noqa: E402

HUB = hub_root()
CUBE = cube_root()
FT = HUB / "08_deliverables" / "insar_full_test"
OUT = FT / "slc_long_pairs"
PREREG = FT / "preregistration" / "X-104.yaml"
MIN_DAYS, WTD_MATCH_CM = 48, 2.0
WINDOWS = {"3 x 3 (9 looks)": (3, 3), "5 x 3 (15 looks)": (3, 5)}
BINS = [(48, 96), (96, 192), (192, 10000)]
KEEP_SHARE = 0.5


def pair_sets(dates: list[str], w: np.ndarray, prot: Protocol) -> dict[str, list[tuple[int, int]]]:
    t = pd.DatetimeIndex(dates)
    m0, m1 = prot.season_months
    free = [i for i in range(len(dates)) if m0 <= t[i].month <= m1 and np.isfinite(w[i])]
    free21 = [i for i in free if t[i].year == 2021]
    ref = max(free21, key=lambda i: w[i])
    single = [(min(ref, j), max(ref, j)) for j in free if abs((t[j] - t[ref]).days) >= MIN_DAYS]
    matched = [(i, j) for k, i in enumerate(free) for j in free[k + 1:]
               if (t[j] - t[i]).days >= MIN_DAYS and abs(w[i] - w[j]) <= WTD_MATCH_CM]
    return {"single reference (highest water table)": single, "matched water table": matched}


def main() -> None:
    assert_preregistered(PREREG, [OUT / "verdict.csv"])
    OUT.mkdir(parents=True, exist_ok=True)
    prot = Protocol.load(FT / "protocol.yaml")
    st = xr.open_dataset(CUBE / "silver" / "static_40m.nc")
    px = pd.read_csv(HUB / "08_deliverables" / "field_first" / "plot_pixels.csv").set_index("plot")
    plots = {(int(r), int(c)) for r, c in zip(px.row, px.col)}
    floors = {lk: S.noise_floor(lk) for lk in (9, 15)}
    rows_out, pair_info = [], []
    for track in X.TRACKS:
        dates, stack, x, y = P.load_stack(track)
        rr, cc = P.sample_cells(x, y, st)
        tg = P.targets(rr, cc, st, plots)
        rows, cols = tg.row.to_numpy(), tg.col.to_numpy()
        sets = P.date_sets(dates, track, prot)
        mask, offs = S.shp_mask(np.abs(stack), rows, cols, *P.SHP_HALF)
        nshp = mask.sum(axis=1)
        rel = nshp >= P.SHP_MIN
        for lk in np.unique(nshp[rel] + 1):
            floors.setdefault(int(lk), S.noise_floor(int(lk)))
        fl_shp = np.array([floors[int(v + 1)] if r else np.nan for v, r in zip(nshp, rel)])
        t = pd.DatetimeIndex(dates)
        for pset, pairs in pair_sets(dates, sets["_wtd"], prot).items():
            for a, b in pairs:
                days = int((t[b] - t[a]).days)
                pair_info.append({"track": track, "pair_set": pset, "t1": dates[a], "t2": dates[b], "days": days,
                                  "wtd1_cm": float(sets["_wtd"][a]), "wtd2_cm": float(sets["_wtd"][b])})
                cohs = {name: (S.boxcar_coherence(stack[a], stack[b], size)[rows, cols], floors[size[0] * size[1]])
                        for name, size in WINDOWS.items()}
                g_shp = S.shp_coherence(stack[a], stack[b], rows, cols, mask, offs)
                cohs["SHP window"] = (np.where(rel, g_shp, np.nan), fl_shp)
                for wname, (g, fl) in cohs.items():
                    for cls, idx in tg.groupby("cls").groups.items():
                        gi, fi = g[idx], np.broadcast_to(fl, g.shape)[idx]
                        ok = np.isfinite(gi)
                        rows_out.append({"track": track, "pair_set": pset, "window": wname, "class": cls,
                                         "t1": dates[a], "t2": dates[b], "days": days, "pixels": int(ok.sum()),
                                         "share_above_floor": float(np.mean(gi[ok] > fi[ok])) if ok.any() else np.nan,
                                         "mean_coherence": float(np.nanmean(gi)) if ok.any() else np.nan})
        del stack
    per_pair = pd.DataFrame(rows_out)
    pairs = pd.DataFrame(pair_info)
    summ = per_pair.groupby(["track", "pair_set", "window", "class"]).agg(
        pairs=("t1", "size"), median_share_above_floor=("share_above_floor", "median"),
        median_mean_coherence=("mean_coherence", "median")).reset_index()
    per_pair["baseline_bin"] = pd.cut(per_pair.days, [b[0] for b in BINS] + [BINS[-1][1]], right=False,
                                      labels=[f"{a}-{b if b < 10000 else ''}" for a, b in BINS])
    by_bin = per_pair.groupby(["track", "pair_set", "window", "class", "baseline_bin"], observed=True).agg(
        pairs=("t1", "size"), median_share_above_floor=("share_above_floor", "median")).reset_index()
    v = verdict(summ)
    for name, df in (("pairs", pairs), ("per_pair", per_pair), ("summary", summ), ("by_baseline", by_bin),
                     ("verdict", v)):
        df.to_csv(OUT / f"{name}.csv", index=False)
    pd.DataFrame([{"looks": k, "noise_floor_95": f} for k, f in sorted(floors.items())]).to_csv(OUT / "noise_floor.csv",
                                                                                               index=False)
    (OUT / "README.md").write_text(readme(summ, by_bin, v, pairs), encoding="utf-8")
    print(v.to_string(index=False))
    print(f"wrote {OUT}")


def verdict(summ: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (pset, window), g in summ.groupby(["pair_set", "window"]):
        keeps = []
        for track in X.TRACKS:
            gg = g[g.track == track].set_index("class")
            if "mat" not in gg.index or "grassland" not in gg.index:
                keeps.append(False)
                continue
            m, gr = gg.loc["mat", "median_share_above_floor"], gg.loc["grassland", "median_share_above_floor"]
            keeps.append(bool(m >= KEEP_SHARE and m >= 0.5 * gr))
        rows.append({"pair_set": pset, "window": window, "mat_keeps_phase_both_tracks": all(keeps)})
    return pd.DataFrame(rows)


def _md(df: pd.DataFrame, nd: int = 2) -> str:
    d = df.copy()
    for c in d.columns:
        if pd.api.types.is_float_dtype(d[c]):
            d[c] = d[c].map(lambda x: "—" if pd.isna(x) else f"{x:.{nd}f}")
    return "\n".join(["| " + " | ".join(map(str, d.columns)) + " |", "|" + "---|" * len(d.columns)] +
                     ["| " + " | ".join(map(str, r)) + " |" for r in d.itertuples(index=False)])


def readme(summ, by_bin, v, pairs) -> str:
    n = pairs.groupby(["track", "pair_set"]).size().rename("pairs").reset_index()
    return f"""# X-104 — long pairs at full resolution on the mat (generated)

Generated by `experiments/slc_longpairs_x104.py` (`make run E=slc_longpairs_x104`); do not edit. Pre-registration:
`../preregistration/X-104.yaml`. CSLCs of X-101 (5 m × 10 m), caveat of its QC v2 phase residual. Long pairs >=
{MIN_DAYS} days, snow-free; single reference = the 2021 date with the highest mat water table; matched = water tables
within {WTD_MATCH_CM} cm. "Above the floor" = coherence above the 95th percentile of an incoherent estimate with the
same number of looks (`noise_floor.csv`), so the estimator's bias at few looks is not counted as phase.

## Verdict (the mat keeps phase if its median share above the floor is >= {KEEP_SHARE} and >= half the grassland's,
both tracks — SETTLED §3's reopen condition for that window)

{_md(v)}

## Pairs

{_md(n)}

## Median over pairs, by class

{_md(summ)}

## By baseline

{_md(by_bin)}

Files: `verdict.csv`, `summary.csv`, `by_baseline.csv`, `per_pair.csv`, `pairs.csv` (with the water tables),
`noise_floor.csv`.
"""


if __name__ == "__main__":
    main()

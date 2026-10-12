"""X-096 — the loss-of-lock atlas: where and when no coherent interferogram spans an acquisition gap.

Definition (Conroy et al. 2023): for one track, consecutive acquisitions e_i < e_{i+1}; the gap is *bridged* at a
pixel if any held pair with t1 <= e_i and t2 >= e_{i+1} has coherence >= THRESHOLD there. Where a gap is not bridged
the chain of phase is broken (loss of lock): integrating across it carries an unknown offset. Every held 40 m pair set
of the cube is used (consecutive, 12-day, longer, annual pairs). Coherence is HyP3's, as delivered.

Pre-registration: ``08_deliverables/insar_full_test/preregistration/X-096.yaml``. Reads the cube's silver interferogram
layers and static layers (no laser); writes only to the hub: ``08_deliverables/insar_full_test/loss_of_lock/``. The P6
segment list it writes is the input of X-093.

    make run E=loss_of_lock_x096
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from insar_wetlands.cube import cube_root
from insar_wetlands.paths import hub_root
from insar_wetlands.scoring import SEASON_OF_MONTH, assert_preregistered

HUB = hub_root()
CUBE = cube_root()
FT = HUB / "08_deliverables" / "insar_full_test"
OUT = FT / "loss_of_lock"
PREREG = FT / "preregistration" / "X-096.yaml"
THRESHOLD = 0.30
SENSITIVITY = (0.25, 0.40)
SETS = ("2017_2019", "2020_2021", "2022_2024", "2025", "2026", "annual")
TRACKS = ("ascending", "descending")
ZONE_A, ZONE_C = 1, 3


def pairs_and_coherence(track: str) -> tuple[pd.DataFrame, np.ndarray]:
    """Every held 40 m pair of one track (deduplicated) and its coherence (pair, y, x)."""
    metas, cohs, seen = [], [], set()
    for s in SETS:
        f = CUBE / "silver" / f"ifg_{s}_{track}.nc"
        if not f.exists():
            continue
        with xr.open_dataset(f) as d:
            keep = [i for i, p in enumerate(d.pair.values) if p not in seen]
            seen.update(d.pair.values[keep])
            metas.append(pd.DataFrame({"pair": d.pair.values[keep], "t1": pd.to_datetime(d.t1.values[keep]),
                                       "t2": pd.to_datetime(d.t2.values[keep])}))
            cohs.append(d.coh.isel(pair=keep).values.astype("float32"))
    meta = pd.concat(metas, ignore_index=True)
    return meta, np.concatenate(cohs, axis=0)


def bridged(meta: pd.DataFrame, coh: np.ndarray, threshold: float) -> tuple[pd.DatetimeIndex, np.ndarray]:
    """Epochs, and for each gap (e_i, e_{i+1}) a (y, x) boolean map: some pair spanning it is coherent."""
    epochs = pd.DatetimeIndex(sorted(set(meta.t1) | set(meta.t2)))
    t1, t2 = meta.t1.to_numpy(), meta.t2.to_numpy()
    ok = np.nan_to_num(coh, nan=0.0) >= threshold
    out = np.zeros((len(epochs) - 1,) + coh.shape[1:], dtype=bool)
    for i in range(len(epochs) - 1):
        span = (t1 <= epochs[i].to_datetime64()) & (t2 >= epochs[i + 1].to_datetime64())
        if span.any():
            out[i] = ok[span].any(axis=0)
    return epochs, out


def segments(epochs: pd.DatetimeIndex, br: np.ndarray) -> pd.DataFrame:
    """Runs of epochs joined by bridged gaps, at one pixel (``br`` 1-D over gaps)."""
    rows, start = [], 0
    for i, b in enumerate(br):
        if not b:
            rows.append({"start": epochs[start], "end": epochs[i], "epochs": i - start + 1})
            start = i + 1
    rows.append({"start": epochs[start], "end": epochs[-1], "epochs": len(epochs) - start})
    return pd.DataFrame(rows)


def main() -> None:
    assert_preregistered(PREREG, [OUT / "gaps.csv"])
    OUT.mkdir(parents=True, exist_ok=True)
    st = xr.open_dataset(CUBE / "silver" / "static_40m.nc")
    zone = st.zone.values
    hard = st.hard_target.values.astype(bool) if "hard_target" in st else np.zeros_like(zone, bool)
    px = pd.read_csv(HUB / "08_deliverables" / "field_first" / "plot_pixels.csv").set_index("plot")
    r6, c6 = int(px.loc["P6", "row"]), int(px.loc["P6", "col"])        # the P6 cell, as X-082 reads it
    masks = {"mat": zone == ZONE_A, "grassland": zone == ZONE_C, "hard_targets": hard}
    gap_rows, seg_rows, rel_rows, sens_rows = [], [], [], []
    for track in TRACKS:
        meta, coh = pairs_and_coherence(track)
        for thr in (THRESHOLD, *SENSITIVITY):
            epochs, br = bridged(meta, coh, thr)
            share = {k: br[:, m].mean(axis=1) for k, m in masks.items() if m.any()}
            if thr != THRESHOLD:
                sens_rows.append({"track": track, "threshold": thr,
                                  **{f"bridged_{k}": float(v.mean()) for k, v in share.items()},
                                  "p6_breaks": int((~br[:, r6, c6]).sum()) if r6 is not None else None})
                continue
            for i in range(len(epochs) - 1):
                gap_rows.append({"track": track, "from": epochs[i], "to": epochs[i + 1],
                                 "days": (epochs[i + 1] - epochs[i]).days,
                                 "season": SEASON_OF_MONTH[epochs[i + 1].month],
                                 **{f"bridged_{k}": float(v[i]) for k, v in share.items()},
                                 "bridged_p6": bool(br[i, r6, c6])})
            if r6 is not None:
                seg_rows.append(segments(epochs, br[:, r6, c6]).assign(track=track))
            seasons = np.array([SEASON_OF_MONTH[e.month] for e in epochs[1:]])
            years = np.array([e.year for e in epochs[1:]])
            rel = xr.Dataset({f"reliability_{s}": (("y", "x"), br[seasons == s].mean(axis=0))
                              for s in ("winter", "spring", "summer", "autumn")},
                             coords={"y": st.y, "x": st.x})
            rel["breaks_per_year"] = (("y", "x"), (~br).sum(axis=0) / len(set(years)))
            rel.attrs.update(track=track, threshold=thr, definition="share of acquisition gaps bridged by a pair with "
                             f"coherence >= {thr} (X-096)")
            rel.to_netcdf(OUT / f"reliability_{track}.nc")
            for s in ("winter", "spring", "summer", "autumn"):
                for k, m in masks.items():
                    if m.any():
                        rel_rows.append({"track": track, "season": s, "class": k,
                                         "bridged_share": float(br[seasons == s][:, m].mean())})
                if r6 is not None:
                    rel_rows.append({"track": track, "season": s, "class": "P6",
                                     "bridged_share": float(br[seasons == s][:, r6, c6].mean())})
    gaps = pd.DataFrame(gap_rows)
    segs = pd.concat(seg_rows, ignore_index=True) if seg_rows else pd.DataFrame()
    rel = pd.DataFrame(rel_rows)
    sens = pd.DataFrame(sens_rows)
    gaps.to_csv(OUT / "gaps.csv", index=False)
    segs.to_csv(OUT / "segments_p6.csv", index=False)
    rel.to_csv(OUT / "reliability_by_season.csv", index=False)
    sens.to_csv(OUT / "sensitivity.csv", index=False)
    (OUT / "README.md").write_text(readme(gaps, segs, rel, sens), encoding="utf-8")
    print(f"wrote {OUT}")


def _md(df: pd.DataFrame, nd: int = 2) -> str:
    d = df.copy()
    for c in d.columns:
        if pd.api.types.is_float_dtype(d[c]):
            d[c] = d[c].map(lambda x: f"{x:.{nd}f}")
    return "\n".join(["| " + " | ".join(map(str, d.columns)) + " |", "|" + "---|" * len(d.columns)] +
                     ["| " + " | ".join(map(str, r)) + " |" for r in d.itertuples(index=False)])


def readme(gaps: pd.DataFrame, segs: pd.DataFrame, rel: pd.DataFrame, sens: pd.DataFrame) -> str:
    piv = rel.pivot_table(index=["track", "class"], columns="season", values="bridged_share").reset_index()
    per_year = gaps.assign(year=pd.to_datetime(gaps["to"]).dt.year).groupby(["track", "year"]).agg(
        gaps=("bridged_p6", "size"), p6_breaks=("bridged_p6", lambda s: int((~s.astype(bool)).sum()))).reset_index()
    breaks = gaps[~gaps.bridged_p6.astype(bool)][["track", "from", "to", "days", "season"]]
    return f"""# X-096 — loss-of-lock atlas (generated)

Generated by `experiments/loss_of_lock_x096.py` (`make run E=loss_of_lock_x096`); do not edit. Pre-registration:
`../preregistration/X-096.yaml`. A gap between consecutive acquisitions is **bridged** at a pixel when a held pair
spanning it has coherence ≥ {THRESHOLD} there; an unbridged gap is a loss of lock (Conroy et al. 2023).

## Share of gaps bridged, by season

{_md(piv)}

## P6: breaks per year

{_md(per_year)}

## P6: every unbridged gap

{_md(breaks) if len(breaks) else "none"}

## P6 segments (input of X-093)

{_md(segs[["track", "start", "end", "epochs"]]) if len(segs) else "—"}

## Sensitivity to the threshold

{_md(sens) if len(sens) else "—"}

Files: `gaps.csv` (every gap: bridged shares by class, P6), `segments_p6.csv`, `reliability_by_season.csv`,
`reliability_<track>.nc` (per pixel: share of gaps bridged per season, breaks per year).
"""


if __name__ == "__main__":
    main()

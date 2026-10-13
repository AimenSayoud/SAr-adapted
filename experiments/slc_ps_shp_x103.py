"""X-103 — persistent scatterers and statistically homogeneous neighbours on the mat, from our full-resolution CSLCs.

PS were found on intact Sphagnum bog (fiaschi2019) and not on vegetation round a lake (osmanolu2015); SHP tests
assume a stationary amplitude (parizzi2010), while the mat's switches with flooding. On the 2021 CSLCs of X-101
(5 m × 10 m, placed right): PS candidates by amplitude dispersion and phase stability; SHP neighbourhoods by the KS test
on amplitude (``insar_wetlands.slc_stats``); per date set — all dates, snow-free wet, snow-free dry (median split of
the mat water table at overpass); mat vs matched grassland and hard targets. Rules and thresholds:
``08_deliverables/insar_full_test/preregistration/X-103.yaml``.

Reads the CSLCs from the Drive root (on the Mac the mirror) and the cube; writes only to the hub:
``08_deliverables/insar_full_test/slc_ps_shp/``. Also used by X-104 (``load_stack``, ``date_sets``, ``targets``).

    INSAR_DRIVE_ROOT=<hub>/05_code/local/drive_mirror make run E=slc_ps_shp_x103
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from scipy.ndimage import uniform_filter

sys.path.insert(0, str(Path(__file__).resolve().parent))
import p6_validation_x082 as X  # noqa: E402
import slc_qc_x101 as Q  # noqa: E402

from insar_wetlands import slc_stats as S  # noqa: E402
from insar_wetlands.cube import cube_root  # noqa: E402
from insar_wetlands.paths import hub_root, resolve_drive  # noqa: E402
from insar_wetlands.scoring import Protocol, assert_preregistered, score  # noqa: E402

HUB = hub_root()
CUBE = cube_root()
FT = HUB / "08_deliverables" / "insar_full_test"
OUT = FT / "slc_ps_shp"
PREREG = FT / "preregistration" / "X-103.yaml"
HOUR = {"ascending": "16:36", "descending": "05:09"}
DA_MAX, GAMMA_MIN = 0.25, 0.7
LOWPASS = (11, 21)                 # rows (10 m) x cols (5 m): about 110 m x 105 m
SHP_HALF = (2, 5)                  # 5 x 11 samples: 50 m x 55 m
SHP_MIN = 20
BOX = (3, 5)                       # 5 x 3 boxcar for pixels that are not SHP-reliable (rows x cols = 30 m x 25 m)
PS_PRESENT_MIN = 10
CLASSES = {1: "mat", 3: "grassland"}


def load_stack(track: str) -> tuple[list[str], np.ndarray, np.ndarray, np.ndarray]:
    folder = resolve_drive() / "slc_pilot" / "2021" / "cslc" / track
    files = {Q.DATE_RE.search(f.name).group(1): f for f in sorted(folder.glob("*.h5"))
             if "static" not in f.name.lower() and Q.DATE_RE.search(f.name)}
    dates = sorted(files)
    first = Q.cslc(files[dates[0]])
    stack = np.empty((len(dates),) + first[0].shape, "complex64")
    for i, d in enumerate(dates):
        stack[i] = Q.cslc(files[d])[0]
    return dates, stack, first[1], first[2]


def sample_cells(x: np.ndarray, y: np.ndarray, st: xr.Dataset) -> tuple[np.ndarray, np.ndarray]:
    """Cube (row, col) of every CSLC sample, as 2-D arrays (-1 outside)."""
    r, c = Q.axis_index(x, y, st.x.values, st.y.values)
    return np.meshgrid(r, c, indexing="ij")


def targets(rr: np.ndarray, cc: np.ndarray, st: xr.Dataset, plots: set[tuple[int, int]]) -> pd.DataFrame:
    """Every CSLC sample in a mat, grassland or hard-target cell, with its class and flags."""
    zone, hard = st.zone.values, st.hard_target.values.astype(bool)
    bw = st.boardwalk_share.values > 0
    inside = (rr >= 0) & (cc >= 0)
    rs, cs = np.nonzero(inside)
    cr, ccol = rr[rs, cs], cc[rs, cs]
    cls = np.full(len(rs), "", object)
    z = zone[cr, ccol]
    for k, name in CLASSES.items():
        cls[z == k] = name
    cls[hard[cr, ccol]] = "hard_targets"
    keep = cls != ""
    return pd.DataFrame({"row": rs[keep], "col": cs[keep], "cell_row": cr[keep], "cell_col": ccol[keep],
                         "cls": cls[keep], "boardwalk": bw[cr[keep], ccol[keep]],
                         "plot_cell": [(a, b) in plots for a, b in zip(cr[keep], ccol[keep])]})


def date_sets(dates: list[str], track: str, prot: Protocol) -> dict[str, list[int]]:
    site = pd.read_csv(CUBE / "silver" / "site_hourly.csv.gz", index_col=0, usecols=["time_utc", "wtd_mat_mean_cm"])
    site.index = pd.to_datetime(site.index, format="ISO8601", utc=True)
    t = pd.DatetimeIndex([pd.Timestamp(f"{d} {HOUR[track]}", tz="UTC") for d in dates])
    w = site.wtd_mat_mean_cm.reindex(t.round("h")).to_numpy()
    m0, m1 = prot.season_months
    free = np.flatnonzero((t.month >= m0) & (t.month <= m1) & np.isfinite(w))
    med = np.median(w[free])
    return {"all": list(range(len(dates))), "wet": [i for i in free if w[i] > med],
            "dry": [i for i in free if w[i] <= med], "_wtd": w}


def phase_stability(stack: np.ndarray, idx: list[int], rows: np.ndarray, cols: np.ndarray) -> np.ndarray:
    acc = np.zeros(len(rows), complex)
    n = 0
    for a, b in zip(idx[:-1], idx[1:]):
        ifg = np.nan_to_num(stack[a] * np.conj(stack[b]))
        lp = uniform_filter(ifg.real, LOWPASS) + 1j * uniform_filter(ifg.imag, LOWPASS)
        z = ifg[rows, cols] * np.conj(lp[rows, cols])
        acc += np.where(np.abs(z) > 0, z / np.maximum(np.abs(z), 1e-12), 0)
        n += 1
    return np.abs(acc) / max(n, 1)


def main() -> None:
    assert_preregistered(PREREG, [OUT / "counts.csv"])
    OUT.mkdir(parents=True, exist_ok=True)
    prot = Protocol.load(FT / "protocol.yaml")
    st = xr.open_dataset(CUBE / "silver" / "static_40m.nc")
    px = pd.read_csv(HUB / "08_deliverables" / "field_first" / "plot_pixels.csv").set_index("plot")
    plots = {(int(r), int(c)) for r, c in zip(px.row, px.col)}
    r6, c6 = int(px.loc["P6", "row"]), int(px.loc["P6", "col"])
    p6_cells = {(r6 + a, c6 + b) for a in (-1, 0, 1) for b in (-1, 0, 1)}
    las = X.laser()[0]
    gold = pd.read_csv(CUBE / "gold" / "plots_pair.csv.gz", low_memory=False,
                       usecols=["track", "grid", "plot", "pair", "d_laser_los_mm"])
    gold = gold[(gold["plot"] == "P6") & (gold.grid == "40m")].set_index(["track", "pair"]).d_laser_los_mm
    floors: dict[int, float] = {}
    counts, rel_rows, p6_rows, maps = [], [], [], {}
    for track in X.TRACKS:
        dates, stack, x, y = load_stack(track)
        rr, cc = sample_cells(x, y, st)
        tg = targets(rr, cc, st, plots)
        rows, cols = tg.row.to_numpy(), tg.col.to_numpy()
        hard_rows, hard_cols = rows[tg.cls == "hard_targets"], cols[tg.cls == "hard_targets"]
        sets = date_sets(dates, track, prot)
        amp_all = np.abs(stack)
        for name in ("all", "wet", "dry"):
            idx = sets[name]
            da = S.amplitude_dispersion(stack[idx])[rows, cols]
            gam = phase_stability(stack, idx, rows, cols)
            mask, offs = S.shp_mask(amp_all[idx], rows, cols, *SHP_HALF)
            nshp = mask.sum(axis=1)
            ps = (da <= DA_MAX) & (gam >= GAMMA_MIN)
            for cls, g in tg.assign(da=da, gam=gam, nshp=nshp, ps=ps).groupby("cls"):
                flagged = g.boardwalk | g.plot_cell
                counts.append({"track": track, "date_set": name, "dates": len(idx), "class": cls, "samples": len(g),
                               "ps_candidates": int(g.ps.sum()), "ps_unflagged": int((g.ps & ~flagged).sum()),
                               "median_dispersion": float(np.nanmedian(g.da)),
                               "shp_reliable_share": float(np.mean(g.nshp >= SHP_MIN)),
                               "median_shp": float(np.median(g.nshp))})
            if name == "all":
                full = np.full(stack.shape[1:], np.nan, "float32")
                for k, v in (("dispersion", da), ("shp_count", nshp), ("ps", ps.astype(float)), ("stability", gam)):
                    m = full.copy()
                    m[rows, cols] = v
                    maps[f"{k}_{track}"] = m
                all_mask, all_offs, all_nshp, all_ps = mask, offs, nshp, ps
        # reliability per date set: consecutive pairs whose coherence clears the floor for its looks
        rel_ok = all_nshp >= SHP_MIN
        for name in ("all", "wet", "dry"):
            idx = sets[name]
            above = []
            for a, b in zip(idx[:-1], idx[1:]):
                g_shp = S.shp_coherence(stack[a], stack[b], rows, cols, all_mask, all_offs)
                g_box = S.boxcar_coherence(stack[a], stack[b], BOX)[rows, cols]
                looks = np.where(rel_ok, all_nshp + 1, BOX[0] * BOX[1])
                for lk in np.unique(looks):
                    floors.setdefault(int(lk), S.noise_floor(int(lk)))
                fl = np.vectorize(lambda v: floors[int(v)])(looks)
                above.append(np.where(rel_ok, g_shp, g_box) > fl)
            share = np.mean(above, axis=0) if above else np.full(len(rows), np.nan)
            for cls, g in tg.assign(share=share).groupby("cls"):
                rel_rows.append({"track": track, "date_set": name, "class": cls, "pairs": len(above),
                                 "median_share_above_floor": float(np.nanmedian(g.share))})
        # PS candidates in P6's cell and its neighbours: increments against the laser
        sel = tg.assign(ps=all_ps)[lambda d: d.ps & [(a, b) in p6_cells for a, b in zip(d.cell_row, d.cell_col)]]
        for _, g in sel.iterrows():
            recs = []
            for a, b in zip(range(len(dates) - 1), range(1, len(dates))):
                ifg = stack[a] * np.conj(stack[b])
                ref = np.nanmean(ifg[hard_rows, hard_cols] / np.maximum(np.abs(ifg[hard_rows, hard_cols]), 1e-12))
                key = f"{dates[a]}_{dates[b]}"
                recs.append({"t1": pd.Timestamp(dates[a], tz="UTC"), "t2": pd.Timestamp(dates[b], tz="UTC"),
                             "est_mm": float(np.angle(ifg[g.row, g.col] * np.conj(ref))) * X.PHASE_TO_MM,
                             "truth_inc_mm": float(gold.get((track, key), np.nan))})
            q = pd.DataFrame(recs).dropna(subset=["truth_inc_mm"])
            inc = score(q, las[track], prot)["increments"].iloc[0].to_dict() if len(q) >= 5 else {}
            p6_rows.append({"track": track, "row": g.row, "col": g.col, "cell": (g.cell_row, g.cell_col), **inc})
        del stack, amp_all
    counts, rel, p6 = map(pd.DataFrame, (counts, rel_rows, p6_rows))
    read = readings(counts)
    for name, df in (("counts", counts), ("reliability", rel), ("p6_ps", p6), ("readings", read)):
        df.to_csv(OUT / f"{name}.csv", index=False)
    xr.Dataset({k: (("sample_row", "sample_col"), v) for k, v in maps.items()}).to_netcdf(OUT / "maps.nc")
    (OUT / "README.md").write_text(readme(counts, rel, p6, read), encoding="utf-8")
    print(read.to_string(index=False))
    print(f"wrote {OUT}")


def readings(c: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for track, g in c.groupby("track"):
        m = g[(g["class"] == "mat")].set_index("date_set")
        gr = g[(g["class"] == "grassland")].set_index("date_set")
        rows.append({"track": track,
                     "ps_present_on_mat": bool(m.loc["all", "ps_unflagged"] >= PS_PRESENT_MIN),
                     "shp_at_grassland_level": bool(m.loc["all", "shp_reliable_share"] >= 0.5 * gr.loc["all", "shp_reliable_share"]),
                     "subseason_changes_neighbourhoods": bool(abs(m.loc["wet", "median_shp"] - m.loc["dry", "median_shp"])
                                                             >= 0.25 * max(m.loc["dry", "median_shp"], 1))})
    return pd.DataFrame(rows)


def _md(df: pd.DataFrame, nd: int = 2) -> str:
    d = df.copy()
    for c in d.columns:
        if pd.api.types.is_float_dtype(d[c]):
            d[c] = d[c].map(lambda x: "—" if pd.isna(x) else f"{x:.{nd}f}")
    return "\n".join(["| " + " | ".join(map(str, d.columns)) + " |", "|" + "---|" * len(d.columns)] +
                     ["| " + " | ".join(map(str, r)) + " |" for r in d.itertuples(index=False)])


def readme(c, rel, p6, read) -> str:
    return f"""# X-103 — PS and SHP on the mat from our full-resolution CSLCs (generated)

Generated by `experiments/slc_ps_shp_x103.py` (`make run E=slc_ps_shp_x103`); do not edit. Pre-registration:
`../preregistration/X-103.yaml`. CSLCs of X-101 (5 m × 10 m, 2021–2022-02), caveat of its QC v2 phase residual.
PS candidate: dispersion <= {DA_MAX}, phase stability >= {GAMMA_MIN}. SHP: KS test, 5 × 11-sample window, reliable at
>= {SHP_MIN} neighbours. Wet / dry: snow-free dates split at the median mat water table. Flagged = boardwalk or plot cells.

## Readings

{_md(read)}

## Counts by track, date set and class

{_md(c)}

## Reliability: median share of consecutive pairs above the noise floor (SHP window where reliable, else 5 × 3)

{_md(rel)}

## PS candidates in P6's cell and its neighbours, against the laser

{_md(p6) if len(p6) else "none"}

Files: `readings.csv`, `counts.csv`, `reliability.csv`, `p6_ps.csv`, `maps.nc` (per sample, all dates: dispersion,
SHP count, PS candidates, phase stability).
"""


if __name__ == "__main__":
    main()

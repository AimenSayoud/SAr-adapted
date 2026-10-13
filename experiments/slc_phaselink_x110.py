"""X-110 — phase linking (Dolphin, EMI / EVD) on the 2021 CSLCs, scored at P6 against the laser.

X-103 (T65) found the mat a homogeneous distributed scatterer at full resolution — the neighbourhoods phase linking on
SLCs needs. Two modes, because phase linking runs only on Colab and the laser lives only in the hub:

- ``link`` (Colab): per track and estimator, Dolphin on the 66 VV CSLCs of X-101 (KS neighbourhoods, half-window 2 × 5
  samples, no unwrapping); the linked phase of every date is then summed per 40 m cube cell (complex) and kept on Drive
  ``slc_pilot/2021/phase_linking/<track>_<estimator>.nc`` with the cell-mean temporal coherence.
- ``score`` (Mac): at P6, the linked series and the plain chain of consecutive pairs from the same CSLCs (same 40 m
  cell, unfiltered), each referenced to the hard-target cells, chained by cumulated wrapped differences and scored
  against the laser by the common scorer. Rule: ``08_deliverables/insar_full_test/preregistration/X-110.yaml``.

    (Colab) python -m insar_wetlands.register run slc_phaselink_x110 -- link
    INSAR_DRIVE_ROOT=<hub>/05_code/local/drive_mirror make run E=slc_phaselink_x110 ARGS="score"
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import slc_qc_x101 as Q  # noqa: E402

from insar_wetlands.cube import cube_root  # noqa: E402
from insar_wetlands.paths import resolve_drive  # noqa: E402

TRACKS = ("ascending", "descending")
ESTIMATORS = {"emi": False, "evd": True}           # use_evd
HALF_WINDOW = {"y": 2, "x": 5}                     # 5 x 11 samples, as X-103
WORK = Path("/content/phase_linking")
MARGIN = 0.20


# ------------------------------------------------------------------------------------------------ link (Colab)

def dolphin_config(files: list[Path], work: Path, use_evd: bool) -> Path:
    import yaml  # noqa: PLC0415
    work.mkdir(parents=True, exist_ok=True)
    cfg_path = work / "dolphin_config.yaml"
    subprocess.run(["dolphin", "config", "--outfile", str(cfg_path), "--work-directory", str(work),
                    "--sds", "/data/VV", "--slc-files", *map(str, files)], check=True)
    cfg = yaml.safe_load(cfg_path.read_text())
    pl = cfg["phase_linking"]
    pl["half_window"] = dict(HALF_WINDOW)
    pl["shp_method"] = "ks"
    pl["use_evd"] = bool(use_evd)
    pl["output_reference_idx"] = 0
    cfg["unwrap_options"]["run_unwrap"] = False
    cfg.setdefault("timeseries_options", {})["run_inversion"] = False
    cfg["timeseries_options"]["run_velocity"] = False
    cfg["worker_settings"]["threads_per_worker"] = 2
    cfg_path.write_text(yaml.safe_dump(cfg, sort_keys=False))
    return cfg_path


def cell_sums(raster: Path, st: xr.Dataset) -> np.ndarray:
    """Complex (or real) raster on the CSLC grid summed per 40 m cube cell; NaN where a cell has few samples."""
    import rasterio  # noqa: PLC0415
    with rasterio.open(raster) as r:
        a = r.read(1)
        t = r.transform
        x = t.c + t.a * (np.arange(r.width) + 0.5)
        y = t.f + t.e * (np.arange(r.height) + 0.5)
    gx, gy = st.x.values, st.y.values
    rows, cols = Q.axis_index(x, y, gx, gy)
    rr, cc = np.meshgrid(rows, cols, indexing="ij")
    ok = (rr >= 0) & (cc >= 0) & np.isfinite(a) & (a != 0)
    idx = (rr * len(gx) + cc)[ok]
    n = len(gy) * len(gx)
    cnt = np.bincount(idx, minlength=n).reshape(len(gy), len(gx))
    if np.iscomplexobj(a):
        z = a / np.maximum(np.abs(a), 1e-12)
        s = (np.bincount(idx, z[ok].real, n) + 1j * np.bincount(idx, z[ok].imag, n)).reshape(len(gy), len(gx))
    else:
        s = (np.bincount(idx, a[ok], n) / np.maximum(np.bincount(idx, minlength=n), 1)).reshape(len(gy), len(gx))
    return np.where(cnt >= 16, s, np.nan)


def link() -> None:
    st = xr.open_dataset(cube_root() / "silver" / "static_40m.nc")
    dest = resolve_drive() / "slc_pilot" / "2021" / "phase_linking"
    dest.mkdir(parents=True, exist_ok=True)
    for track in TRACKS:
        src = resolve_drive() / "slc_pilot" / "2021" / "cslc" / track
        local = WORK / "cslc" / track
        local.mkdir(parents=True, exist_ok=True)
        files = []
        for f in sorted(src.glob("*.h5")):
            if "static" in f.name.lower() or not Q.DATE_RE.search(f.name):
                continue
            out = local / f.name
            if not out.exists():
                shutil.copy2(f, out)
            files.append(out)
        for est, use_evd in ESTIMATORS.items():
            target = dest / f"{track}_{est}.nc"
            if target.exists():
                print(f"{target.name} on Drive; skipped", flush=True)
                continue
            work = WORK / f"{track}_{est}"
            cfg = dolphin_config(files, work, use_evd)
            subprocess.run(["dolphin", "run", str(cfg)], check=True)
            linked = sorted((work / "linked_phase").glob("*.tif"))
            linked = [f for f in linked if Q.DATE_RE.search(f.name)]
            dates = [Q.DATE_RE.findall(f.name)[-1] for f in linked]
            z = np.stack([cell_sums(f, st) for f in linked])
            tc = sorted(work.rglob("temporal_coherence*.tif"))
            tcoh = cell_sums(tc[0], st) if tc else np.full(z.shape[1:], np.nan)
            ds = xr.Dataset({"linked_re": (("date", "y", "x"), z.real.astype("float32")),
                             "linked_im": (("date", "y", "x"), z.imag.astype("float32")),
                             "temporal_coherence": (("y", "x"), tcoh.astype("float32"))},
                            coords={"date": dates, "y": st.y.values, "x": st.x.values},
                            attrs={"track": track, "estimator": est, "half_window": str(HALF_WINDOW),
                                   "shp_method": "ks", "dolphin_config": cfg.read_text()})
            ds.to_netcdf(target)
            print(f"{target.name}: {len(dates)} dates", flush=True)


# ------------------------------------------------------------------------------------------------ score (Mac)

def chain_from_phase(dates: list[str], phase: np.ndarray) -> pd.DataFrame:
    """Consecutive increments (mm) of a referenced phase series: wrapped differences (cycle 0)."""
    import p6_validation_x082 as X  # noqa: PLC0415
    d = np.angle(np.exp(1j * (phase[1:] - phase[:-1])))
    return pd.DataFrame({"t1": [pd.Timestamp(a, tz="UTC") for a in dates[:-1]],
                         "t2": [pd.Timestamp(b, tz="UTC") for b in dates[1:]], "est_mm": d * X.PHASE_TO_MM})


def score_mode() -> None:
    import p6_validation_x082 as X  # noqa: PLC0415

    from insar_wetlands.paths import hub_root  # noqa: PLC0415
    from insar_wetlands.scoring import Protocol, assert_preregistered, score  # noqa: PLC0415
    hub = hub_root()
    ft = hub / "08_deliverables" / "insar_full_test"
    out = ft / "slc_phase_linking"
    assert_preregistered(ft / "preregistration" / "X-110.yaml", [out / "verdict.csv"])
    out.mkdir(parents=True, exist_ok=True)
    prot = Protocol.load(ft / "protocol.yaml")
    st = xr.open_dataset(cube_root() / "silver" / "static_40m.nc")
    hard = st.hard_target.values.astype(bool)
    zone = st.zone.values
    px = pd.read_csv(hub / "08_deliverables" / "field_first" / "plot_pixels.csv").set_index("plot")
    r6, c6 = int(px.loc["P6", "row"]), int(px.loc["P6", "col"])
    cells = [(r6 + a, c6 + b) for a in (-1, 0, 1) for b in (-1, 0, 1)]
    gold = pd.read_csv(cube_root() / "gold" / "plots_pair.csv.gz", low_memory=False,
                       usecols=["track", "grid", "plot", "pair", "d_laser_los_mm"])
    gold = gold[(gold["plot"] == "P6") & (gold.grid == "40m")].set_index(["track", "pair"]).d_laser_los_mm
    las = X.laser()[0]
    rows, tcoh_rows = [], []

    def scored(track, name, q, cell):
        q = q.assign(truth_inc_mm=[float(gold.get((track, f"{a:%Y%m%d}_{b:%Y%m%d}"), np.nan))
                                   for a, b in zip(q.t1, q.t2)])
        sc = score(q, las[track], prot)
        lev = sc["levels"]
        season = lev[lev.view.str.startswith("season")]
        full = lev[lev.view == "full span"]
        inc = sc["increments"].iloc[0].to_dict() if len(sc.get("increments", [])) else {}
        return {"track": track, "estimate": name, "cell": cell,
                "season_rmse_median_mm": float(season.rmse_detrended_mm.median()) if len(season) else np.nan,
                "season_views": len(season),
                "full_span_rmse_mm": float(full.rmse_detrended_mm.iloc[0]) if len(full) else np.nan,
                "full_span_pearson": float(full.pearson_detrended.iloc[0]) if len(full) else np.nan,
                "increment_rmse_mm": inc.get("rmse_increments_mm", np.nan),
                "increment_pearson": inc.get("pearson_increments", np.nan)}

    for track in TRACKS:
        folder = resolve_drive() / "slc_pilot" / "2021" / "cslc" / track
        files = {Q.DATE_RE.search(f.name).group(1): f for f in sorted(folder.glob("*.h5"))
                 if "static" not in f.name.lower() and Q.DATE_RE.search(f.name)}
        dates = sorted(files)
        # the plain chain of consecutive pairs from the same CSLCs (40 m cells, unfiltered), hard-target referenced
        gx, gy = st.x.values, st.y.values
        shape = (len(gy), len(gx))
        prev, recs = None, {c: [] for c in cells}
        for a, b in zip(dates[:-1], dates[1:]):
            s1, x, y = Q.cslc(files[a]) if prev is None else prev
            s2 = Q.cslc(files[b])
            r_, c_ = Q.axis_index(x, y, gx, gy)
            z, _ = Q.look(s1, s2[0], r_, c_, shape)
            ref = np.nanmean(z[hard] / np.maximum(np.abs(z[hard]), 1e-12))
            for cell in cells:
                recs[cell].append({"t1": pd.Timestamp(a, tz="UTC"), "t2": pd.Timestamp(b, tz="UTC"),
                                   "est_mm": float(np.angle(z[cell] * np.conj(ref))) * X.PHASE_TO_MM})
            prev = s2
        for cell in cells:
            rows.append(scored(track, "chain of consecutive pairs (same CSLCs)", pd.DataFrame(recs[cell]), cell))
        for est in ESTIMATORS:
            f = resolve_drive() / "slc_pilot" / "2021" / "phase_linking" / f"{track}_{est}.nc"
            if not f.exists():
                continue
            d = xr.open_dataset(f)
            z = d.linked_re.values + 1j * d.linked_im.values
            ld = list(d.date.values.astype(str))
            ref = np.array([np.nanmean(z[i][hard] / np.maximum(np.abs(z[i][hard]), 1e-12)) for i in range(len(ld))])
            ph = np.angle(z * np.conj(ref)[:, None, None])
            for cell in cells:
                rows.append(scored(track, f"phase linking ({est.upper()})", chain_from_phase(ld, ph[:, cell[0], cell[1]]),
                                   cell))
            tc = d.temporal_coherence.values
            for k, m in (("mat", zone == 1), ("grassland", zone == 3), ("hard_targets", hard)):
                tcoh_rows.append({"track": track, "estimator": est, "class": k,
                                  "median_temporal_coherence": float(np.nanmedian(tc[m]))})
    sc = pd.DataFrame(rows)
    tco = pd.DataFrame(tcoh_rows)
    v = verdict(sc, (r6, c6))
    sc.to_csv(out / "scores.csv", index=False)
    tco.to_csv(out / "temporal_coherence.csv", index=False)
    v.to_csv(out / "verdict.csv", index=False)
    (out / "README.md").write_text(readme(sc, tco, v, (r6, c6)), encoding="utf-8")
    print(v.to_string(index=False))
    print(f"wrote {out}")


def verdict(sc: pd.DataFrame, p6) -> pd.DataFrame:
    p = sc[sc.cell == p6]
    rows = []
    for est in sorted(set(p.estimate) - {"chain of consecutive pairs (same CSLCs)"}):
        rel, full_ok = [], []
        for track in TRACKS:
            a = p[(p.track == track) & (p.estimate == est)]
            b = p[(p.track == track) & (p.estimate == "chain of consecutive pairs (same CSLCs)")]
            if not len(a) or not len(b):
                rel.append(np.nan)
                full_ok.append(False)
                continue
            rel.append(a.season_rmse_median_mm.iloc[0] / b.season_rmse_median_mm.iloc[0] - 1)
            full_ok.append(a.full_span_rmse_mm.iloc[0] <= b.full_span_rmse_mm.iloc[0])
        if all(np.isfinite(rel)) and all(r <= -MARGIN for r in rel) and all(full_ok):
            res = "beats the chain"
        elif all(np.isfinite(rel)) and all(abs(r) <= MARGIN for r in rel):
            res = "matches the chain"
        else:
            res = "loses to the chain" if all(np.isfinite(rel)) and all(r > MARGIN for r in rel) else "mixed"
        rows.append({"estimate": est, "season_rmse_change_asc": rel[0], "season_rmse_change_desc": rel[1],
                     "full_span_not_higher_both": all(full_ok), "result": res})
    return pd.DataFrame(rows)


def _md(df: pd.DataFrame, nd: int = 2) -> str:
    d = df.copy()
    for c in d.columns:
        if pd.api.types.is_float_dtype(d[c]):
            d[c] = d[c].map(lambda x: "—" if pd.isna(x) else f"{x:.{nd}f}")
    return "\n".join(["| " + " | ".join(map(str, d.columns)) + " |", "|" + "---|" * len(d.columns)] +
                     ["| " + " | ".join(map(str, r)) + " |" for r in d.itertuples(index=False)])


def readme(sc, tco, v, p6) -> str:
    at = sc[sc.cell == p6].drop(columns=["cell"])
    nine = sc.groupby(["track", "estimate"])[["season_rmse_median_mm", "full_span_rmse_mm"]].median().reset_index()
    return f"""# X-110 — phase linking on the CSLCs, scored at P6 (generated)

Generated by `experiments/slc_phaselink_x110.py` (`link` on Colab, `score` on the Mac); do not edit. Pre-registration:
`../preregistration/X-110.yaml`. Dolphin (KS neighbourhoods, half-window {HALF_WINDOW}, no unwrapping) on the 66 VV CSLCs
per track of X-101 (caveat of its QC v2 phase residual); linked phase summed per 40 m cell, referenced to the hard-target
cells, chained by wrapped consecutive differences, scored against the P6 laser like the plain chain of the same CSLCs.
Tried on P6, 2021 (D-028).

## Verdict (beats the chain if the within-season level RMSE is lower by >= {MARGIN:.0%} on both tracks and the full-span
RMSE is not higher)

{_md(v)}

## At the P6 cell

{_md(at)}

## Median over the nine cells around P6

{_md(nine)}

## Temporal coherence of the linking (median by class)

{_md(tco)}
"""


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("mode", choices=["link", "score"])
    a = ap.parse_args(argv)
    link() if a.mode == "link" else score_mode()


if __name__ == "__main__":
    main()

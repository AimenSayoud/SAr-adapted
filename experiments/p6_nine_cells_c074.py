"""C-074 — the nine cells around P6 on the current cube, with dated tree cover: backscatter, coherence and how each
cell's Sentinel-1 phase follows the P6 laser (the professor's request; X-079 read them on earlier inputs).

Per cell × track × group (X-126's consecutive pairs: six-day 2021 from 27 May, twelve-day 2022–2024, six-day 2025;
hard-target reference; the cube's laser on the snow-evidence rule): r and slope through the origin against the laser,
each with a pair-bootstrap interval; median coherence; per cell × track × period, VV / VH medians and the temporal SD
of VV; per cell, tree cover from the 2024 laser scan (1 m canopy shares) and 2011, boardwalk, depth into the mat.
P6's eight neighbours are also read as one mean. Descriptive: nine cells cannot separate trees from depth (X-079);
X-132 answers that with matched cells.

Pre-registration: hub ``08_deliverables/insar_full_test/preregistration/C-074.yaml``. Stops unless P6's own cell
reproduces X-126. Writes only to the hub (``08_deliverables/supervisor_pack_c074/tables/``).

    make run E=p6_nine_cells_c074
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import x065_hard_targets_x126 as X126  # noqa: E402  (X-126's pairs, groups and laser)
from insar_wetlands.inversion.isbas import PHASE_TO_MM  # noqa: E402
from insar_wetlands.paths import cube_root, hub_root  # noqa: E402
from insar_wetlands.scoring import Protocol, assert_preregistered, guard_lockbox, increment_stats  # noqa: E402

HUB, CUBE = hub_root(), cube_root()
FT = HUB / "08_deliverables" / "insar_full_test"
OUT = HUB / "08_deliverables" / "supervisor_pack_c074" / "tables"
PREREG = FT / "preregistration" / "C-074.yaml"
X126_TABLE = FT / "x065_hard_targets" / "by_revisit.csv"
TRACKS = ("ascending", "descending")
GROUPS = {"2021": 6, "2022–2024": 12, "2025": 6}
GROUP_PERIOD = {"2021": "2020_2021", "2022–2024": "2022_2024", "2025": "2025"}
PERIOD_YEARS = {"2021": (2021, 2021), "2022–2024": (2022, 2024), "2025": (2025, 2025)}
NAMES = {(-1, -1): "NW", (-1, 0): "N", (-1, 1): "NE", (0, -1): "W", (0, 0): "P6", (0, 1): "E",
         (1, -1): "SW", (1, 0): "S", (1, 1): "SE"}
N_BOOT, SEED = 2000, 0
TREE_M, SHRUB_M = 2.0, 0.5


def nine(st) -> pd.DataFrame:
    pu = pd.read_csv(CUBE / "silver" / "plot_units.csv").set_index("plot")
    r6, c6 = int(pu.loc["P6", "row"]), int(pu.loc["P6", "col"])
    rows = []
    for (dr, dc), name in NAMES.items():
        r, c = r6 + dr, c6 + dc
        rows.append({"cell": name, "row": r, "col": c, "x": float(st.x[c]), "y": float(st.y[r])})
    return pd.DataFrame(rows).set_index("cell")


def descriptors(cells: pd.DataFrame, st: xr.Dataset) -> pd.DataFrame:
    """Per cell: the 2024 laser scan's 1 m canopy shares (≥ 2 m tree, 0.5–2 m shrub, < 0.5 m low), p90 2024 and 2011,
    boardwalk share, depth into the mat, zone."""
    li = xr.open_dataset(CUBE / "silver" / "lidar_dated_1m.nc")
    out = cells.copy()
    for k in ("zone", "boardwalk_share", "depth_in_mat_m", "canopy2024_p90_m", "canopy2011_p90_m", "canopy2024_cover"):
        out[k] = [float(st[k].values[r, c]) for r, c in zip(cells.row, cells["col"])]
    for name, cr in cells.iterrows():
        sub = li.chm_2024_m.sel(x=slice(cr.x - 20, cr.x + 20), y=slice(cr.y + 20, cr.y - 20)).values
        v = sub[np.isfinite(sub)]
        out.loc[name, "tree_share_2024"] = float((v >= TREE_M).mean())
        out.loc[name, "shrub_share_2024"] = float(((v >= SHRUB_M) & (v < TREE_M)).mean())
        out.loc[name, "low_share_2024"] = float((v < SHRUB_M).mean())
        out.loc[name, "chm_cells_1m"] = int(v.size)
        s11 = li.chm_2011_m.sel(x=slice(cr.x - 20, cr.x + 20), y=slice(cr.y + 20, cr.y - 20)).values
        s11 = s11[np.isfinite(s11)]
        out.loc[name, "tree_share_2011"] = float((s11 >= TREE_M).mean())
    out.attrs["dates"] = {k: li.attrs[k] for k in li.attrs if k.startswith("date_")}
    return out


def boot(a: np.ndarray, b: np.ndarray, rng) -> dict:
    m = np.isfinite(a) & np.isfinite(b)
    a, b = a[m], b[m]
    r, s = [], []
    for _ in range(N_BOOT):
        i = rng.integers(0, len(a), len(a))
        if np.std(a[i]) > 0 and np.std(b[i]) > 0:
            r.append(np.corrcoef(a[i], b[i])[0, 1])
        if np.sum(a[i] ** 2) > 0:
            s.append(np.sum(a[i] * b[i]) / np.sum(a[i] ** 2))
    return {"r_lo95": float(np.percentile(r, 2.5)), "r_hi95": float(np.percentile(r, 97.5)),
            "slope_lo95": float(np.percentile(s, 2.5)), "slope_hi95": float(np.percentile(s, 97.5))}


def nearest_rtc(rtc: xr.Dataset, times: pd.Series, var: str, ry, cx, tol_h: float = 12.0) -> np.ndarray:
    """The RTC value of each cell at the scene nearest each time (within ``tol_h``), else NaN; (times × cells)."""
    rt = pd.DatetimeIndex(rtc.time.values).tz_localize("UTC")
    t = pd.DatetimeIndex(times)
    pos = rt.get_indexer(t, method="nearest")
    ok = np.abs((rt[pos] - t).total_seconds()) <= tol_h * 3600
    a = rtc[var].isel(y=ry, x=cx).values                       # (time, cell)
    out = a[pos].astype(float)
    out[~ok] = np.nan
    return out


def backscatter(rtc, ry, cx, names, track: str, y0: int, y1: int) -> list[dict]:
    """Per cell: median VV and VH, the temporal SD of VV and the median VH − VV over the scenes of years y0–y1."""
    sel = (rtc.time.dt.year >= y0) & (rtc.time.dt.year <= y1)
    vv = rtc.vv_db.isel(y=ry, x=cx).sel(time=sel).values
    vh = rtc.vh_db.isel(y=ry, x=cx).sel(time=sel).values
    return [{"cell": lab, "track": track, "years": f"{y0}–{y1}" if y1 > y0 else str(y0),
             "scenes": int(np.isfinite(vv[:, j]).sum()), "vv_median_db": float(np.nanmedian(vv[:, j])),
             "vh_median_db": float(np.nanmedian(vh[:, j])), "vv_sd_db": float(np.nanstd(vv[:, j])),
             "vh_minus_vv_db": float(np.nanmedian(vh[:, j] - vv[:, j]))} for j, lab in enumerate(names)]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    names = ["nine_cells.csv", "nine_cells_by_group.csv", "nine_cells_backscatter.csv", "nine_cells_pairs.csv",
             "reproduction.csv"]
    assert_preregistered(PREREG, [OUT / f for f in names])
    st = xr.open_dataset(CUBE / "silver" / "static_40m.nc")
    cells = nine(st)
    desc = descriptors(cells, st)
    p = X126.load_pairs()
    p = p[(p.between == 0) & p.group.isin(list(GROUPS))]
    p = p[[GROUPS[g] == d for g, d in zip(p.group, p.revisit_days)]].reset_index(drop=True)
    guard_lockbox(pd.concat([p.t1, p.t2]), Protocol.load(FT / "protocol.yaml"))
    ry = xr.DataArray(cells.row.to_numpy(), dims="cell")
    cx = xr.DataArray(cells["col"].to_numpy(), dims="cell")
    ref = pd.read_csv(X126_TABLE)
    rng = np.random.default_rng(SEED)
    by, pairs_rows, checks, bs = [], [], [], []
    for track in TRACKS:
        rtc = xr.open_dataset(CUBE / "silver" / f"rtc_{track}.nc")
        for grp, rev in GROUPS.items():
            g = p[(p.track == track) & (p.group == grp)].sort_values("t1").reset_index(drop=True)
            ds = xr.open_dataset(CUBE / "silver" / f"ifg_{GROUP_PERIOD[grp]}_{track}.nc")
            pos = pd.Index(ds.pair.values).get_indexer(g.pair)
            unw = ds.unw.isel(y=ry, x=cx).values[pos].astype(float)
            coh = ds.coh.isel(y=ry, x=cx).values[pos].astype(float)
            inc = (unw - g.ref_hard_unw.to_numpy(float)[:, None]) * PHASE_TO_MM
            nb = [i for i, c in enumerate(cells.index) if c != "P6"]
            inc_all = np.column_stack([inc, np.nanmean(inc[:, nb], 1)])
            coh_all = np.column_stack([coh, np.nanmean(coh[:, nb], 1)])
            labels = list(cells.index) + ["mean of the 8 neighbours"]
            a = g.d_laser_los_mm.to_numpy(float)
            vv1, vv2 = (nearest_rtc(rtc, g[c].dt.tz_convert("UTC"), "vv_db", ry, cx) for c in ("t1", "t2"))
            vh1, vh2 = (nearest_rtc(rtc, g[c].dt.tz_convert("UTC"), "vh_db", ry, cx) for c in ("t1", "t2"))
            for j, lab in enumerate(labels):
                stt = increment_stats(a, inc_all[:, j])
                by.append({"cell": lab, "track": track, "group": grp, "revisit_days": rev, **stt, **boot(a, inc_all[:, j], rng),
                           "median_coh": float(np.nanmedian(coh_all[:, j])),
                           "same_sign_share": float(np.mean(np.sign(a[np.isfinite(a) & np.isfinite(inc_all[:, j])]) ==
                                                            np.sign(inc_all[:, j][np.isfinite(a) & np.isfinite(inc_all[:, j])])))})
                if lab == "P6":
                    th = ref[(ref.track == track) & (ref.period.astype(str) == grp) & (ref.revisit_days == rev) &
                             (ref.pairs_kind == "consecutive") & (ref.reference == "hard_targets")].iloc[0]
                    checks.append({"check": f"X-126 {track} {grp}", "pairs": stt["pairs"], "pairs_x126": int(th.pairs),
                                   "abs_diff_r": abs(stt["pearson_increments"] - th.pearson_increments),
                                   "ok": bool(stt["pairs"] == th.pairs and
                                              abs(stt["pearson_increments"] - th.pearson_increments) < 1e-9)})
            for k in range(len(g)):
                for j, lab in enumerate(cells.index):
                    pairs_rows.append({"track": track, "group": grp, "pair": g.pair[k], "t1": g.t1[k], "t2": g.t2[k],
                                       "revisit_days": rev, "cell": lab, "s1_los_mm": inc[k, j], "coherence": coh[k, j],
                                       "laser_los_mm": a[k], "d_vv_db": vv2[k, j] - vv1[k, j],
                                       "d_vh_db": vh2[k, j] - vh1[k, j]})
            bs += backscatter(rtc, ry, cx, cells.index, track, *PERIOD_YEARS[grp])
        bs += backscatter(rtc, ry, cx, cells.index, track, 2021, 2025)      # all the laser years together
    rep = pd.DataFrame(checks)
    rep.to_csv(OUT / "reproduction.csv", index=False)
    if not rep.ok.all():
        print(rep.to_string())
        sys.exit("stop (pre-registered): P6's own cell does not reproduce X-126; see reproduction.csv")
    t = pd.DataFrame(by)
    p6s = t[t.cell == "P6"].set_index(["track", "group"]).slope_through_origin
    t["r_above_zero"] = t.r_lo95 > 0
    t["slope_interval_contains_p6"] = [(r.slope_lo95 <= p6s[(r.track, r.group)] <= r.slope_hi95) for r in t.itertuples()]
    desc.reset_index().to_csv(OUT / "nine_cells.csv", index=False)
    t.to_csv(OUT / "nine_cells_by_group.csv", index=False)
    pd.DataFrame(bs).to_csv(OUT / "nine_cells_backscatter.csv", index=False)
    pd.DataFrame(pairs_rows).to_csv(OUT / "nine_cells_pairs.csv", index=False)
    print(rep.to_string())
    print(desc[["tree_share_2024", "shrub_share_2024", "boardwalk_share", "depth_in_mat_m"]].round(2).to_string())
    print(t[["cell", "track", "group", "pairs", "pearson_increments", "r_lo95", "slope_through_origin", "slope_lo95",
             "slope_hi95", "median_coh"]].round(2).to_string())


if __name__ == "__main__":
    main()

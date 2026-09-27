"""X-053: where on the mat does the radar follow the water? Per-pixel maps, 2020–2024.

The plots are one edge-to-interior gradient (X-049): 8 radar units cannot separate vegetation from
hydrology. The radar itself has hundreds of mat pixels. For every pixel, on both tracks, with the
2022–2024 stacks and the 2020–2021 extension (D-020):

1. **Coherence vs water-table change** — pairs ≤ 24 d, not frozen: coherence with its annual cycle
   and baseline trend regressed out (the pair analogue of anomalies, as `field_link.season_residual`)
   against |ΔWTD| (median over the plots, dry wells out). r, slope, circular-shift p (pairs ordered
   by mid-date), Benjamini–Hochberg q over the mat.
2. **VV vs water-table level** (ascending RTC, 2022–2024) — anomaly correlation, circular-shift p.
3. **Wet penalty** — season-cleaned coherence of pairs whose two dates were dry minus pairs with a
   wet date (RH ≥ 95 % or rain); p from the per-date wet flags circularly shifted (as X-050).
4. **Transect** — the pixel column through the plots, edge to edge: every map above plus mean
   coherence, VH, temporal coherence, NDWI and distance to the edge.
5. **Monthly atlas** — mean short-pair coherence per calendar month, both tracks.
6. **Arc table** — per pair: zone A and C mean coherence, ΔWTD, wet / frozen flags (for the arc view).

The regional water table stands in for every pixel: WTD is only measured at the plots, but its
changes are spatially coherent over the mat. Outputs only to the private hub (derived from field
data).

    INSAR_DRIVE_ROOT=<snapshot> PYTHONPATH=src python scripts/field_pixel_x053.py
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import xarray as xr  # noqa: E402

from insar_wetlands import field  # noqa: E402
from insar_wetlands import field_link as fl  # noqa: E402
from insar_wetlands.bootstrap import start  # noqa: E402
from insar_wetlands.stack import list_pairs, load_layer  # noqa: E402
from insar_wetlands.stratify import signed_distance_to_aoi  # noqa: E402

OVERPASS = {"ascending": "16:36", "descending": "05:09"}
N_SHIFT = 1000
MIN_VALID = 0.8           # a pixel needs coherence in ≥ 80 % of the pairs used
rng = np.random.default_rng(0)


def iso(d: str) -> str:
    return f"{d[:4]}-{d[4:6]}-{d[6:]}"


def fill_columns(y: np.ndarray, min_valid: float = MIN_VALID) -> tuple[np.ndarray, np.ndarray]:
    """NaN → column mean where ≥ min_valid of a column is valid; mask of usable columns."""
    ok = np.isfinite(y).mean(axis=0) >= min_valid
    mean = np.nanmean(np.where(np.isfinite(y), y, np.nan), axis=0)
    out = np.where(np.isfinite(y), y, mean)
    out[:, ~ok] = np.nan
    return out, ok


def residualize(y: np.ndarray, m: np.ndarray) -> np.ndarray:
    """y minus its least-squares fit on design m (all columns at once)."""
    b, *_ = np.linalg.lstsq(m, np.nan_to_num(y), rcond=None)
    return y - m @ b


def shift_corr(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return fl.circular_shift_p_many(x, y, n_shift=N_SHIFT, rng=rng)


def bh_q(p: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """Benjamini–Hochberg q-values over the pixels in mask (the mat)."""
    q = np.full(p.shape, np.nan)
    q[mask] = fl.bh_qvalues(p[mask])
    return q


def wet_penalty(wet_date: np.ndarray, members: np.ndarray, res: np.ndarray) -> np.ndarray:
    """Per pixel: mean season-cleaned coherence of pairs with two dry dates minus pairs with a wet one."""
    hit = wet_date[members].any(axis=1)
    w = np.where(hit, -1 / max(hit.sum(), 1), 1 / max((~hit).sum(), 1))
    return w @ np.nan_to_num(res)


def season_design(mid: pd.Series, dt: np.ndarray | None = None) -> np.ndarray:
    doy = mid.dt.dayofyear.to_numpy() / 365.25 * 2 * np.pi
    cols = [np.ones(len(mid)), np.cos(doy), np.sin(doy)]
    if dt is not None:
        cols.append(dt.astype(float))
    return np.column_stack(cols)


def load_track(ctx, track: str, extra: Path) -> tuple[pd.DataFrame, np.ndarray]:
    """All pairs of a track (2022–2024 + 2020–2021) and their coherence cube (pairs, y, x)."""
    metas, cubes = [], []
    for stack, d in (("2022–2024", ctx.paths.for_phase("x", track=track).cropped), ("2020–2021", extra)):
        ps = list_pairs(d)
        cubes.append(load_layer(d, "corr", ps).values.astype("float32"))
        metas.append(pd.DataFrame({"pair": ps, "stack": stack}))
    meta = pd.concat(metas, ignore_index=True)
    meta["date1"] = [iso(p.split("_")[0]) for p in meta.pair]
    meta["date2"] = [iso(p.split("_")[1]) for p in meta.pair]
    d1, d2 = pd.to_datetime(meta.date1), pd.to_datetime(meta.date2)
    meta["dt_days"] = (d2 - d1).dt.days
    meta["mid"] = d1 + (d2 - d1) / 2
    return meta, np.concatenate(cubes)


def regional_dwtd(meta: pd.DataFrame, track: str, wtd: pd.DataFrame, cens: dict) -> np.ndarray:
    out = []
    for a, b in zip(meta.date1, meta.date2):
        t1, t2 = (pd.Timestamp(f"{x} {OVERPASS[track]}", tz="UTC") for x in (a, b))
        dw = [field._interp_at(wtd[p], t2) - field._interp_at(wtd[p], t1) for p in field.PLOTS
              if not (bool(cens[p].asof(t1)) or bool(cens[p].asof(t2)))]
        out.append(float(np.median(dw)) if len(dw) >= 5 else np.nan)
    return np.array(out)


def main(argv=None):
    hub = field.field_root().parents[1]
    local = hub / "05_code" / "local"
    ap = argparse.ArgumentParser()
    ap.add_argument("--extra", default=str(local / "s1_2020_2021"))
    ap.add_argument("--wetness", default=str(hub / "08_deliverables" / "field_dew_x050" / "wetness_at_overpasses_2020_2024.csv"))
    ap.add_argument("--out", default=str(hub / "08_deliverables" / "field_pixel_x053"))
    ap.add_argument("--summary-only", action="store_true", help="re-derive summary_by_zone.csv from pixel_maps.nc")
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    ctx = start("field_pixel_x053", mount=False, git=False)
    if a.summary_only:
        with xr.open_dataset(out / "pixel_maps.nc") as ds:
            zs = zone_summary(ds.load(), ctx.zones)
        zs.to_csv(out / "summary_by_zone.csv", index=False)
        print(zs.round(3).to_string(index=False))
        return
    D, cfg, tpl, zones = ctx.paths.drive, ctx.cfg, ctx.template, ctx.zones
    shape = tpl.shape
    zone_a, zone_c = zones["A"].values, zones["C"].values
    wtd = field.load_wtd_hourly()
    cens = {p: field.censored_flag(wtd[p]) for p in field.PLOTS}
    wt = pd.read_csv(a.wetness)

    maps: dict[str, np.ndarray] = {}
    summary, arcs, monthly = [], [], {}
    for track in OVERPASS:
        meta, cube = load_track(ctx, track, Path(a.extra) / f"hyp3_cropped_{track}")
        flags = wt[wt.track == track].set_index("date")
        meta["wet_any"] = [bool(flags.wet.get(x, False) or flags.wet.get(y, False)) for x, y in zip(meta.date1, meta.date2)]
        meta["frozen_any"] = [bool(flags.frozen.get(x, True) or flags.frozen.get(y, True)) for x, y in zip(meta.date1, meta.date2)]
        meta["dwtd_median_cm"] = regional_dwtd(meta, track, wtd, cens)
        flat = cube.reshape(len(meta), -1)
        meta["coh_A"] = np.nanmean(flat[:, zone_a.ravel()], axis=1)
        meta["coh_C"] = np.nanmean(flat[:, zone_c.ravel()], axis=1)
        arcs.append(meta.assign(track=track).drop(columns=["mid"]))

        # 1. coherence vs |ΔWTD| (pairs ≤ 24 d, not frozen, with a regional ΔWTD)
        sel = (meta.dt_days <= 24) & ~meta.frozen_any & meta.dwtd_median_cm.notna()
        order = np.argsort(meta.mid[sel].to_numpy())
        m = meta[sel].iloc[order]
        y, ok = fill_columns(flat[sel.to_numpy()][order])
        res = residualize(y, season_design(m.mid, m.dt_days.to_numpy()))
        x = np.abs(m.dwtd_median_cm.to_numpy())
        r, p = shift_corr(x, res)
        slope = r * np.nanstd(res, axis=0) / x.std()
        for k, v in (("r", r), ("p", p), ("slope", slope)):
            maps[f"coh_dwtd_{k}_{track}"] = v.reshape(shape)
        maps[f"coh_dwtd_q_{track}"] = bh_q(maps[f"coh_dwtd_p_{track}"], zone_a)

        # 3. wet penalty (pairs ≤ 24 d, not frozen; ΔWTD not needed)
        sel3 = (meta.dt_days <= 24) & ~meta.frozen_any
        o3 = np.argsort(meta.mid[sel3].to_numpy())
        m3 = meta[sel3].iloc[o3]
        y3, _ = fill_columns(flat[sel3.to_numpy()][o3])
        res3 = residualize(y3, season_design(m3.mid, m3.dt_days.to_numpy()))
        dates = sorted(set(m3.date1) | set(m3.date2))
        idx = {d: i for i, d in enumerate(dates)}
        mem = np.array([[idx[p], idx[q]] for p, q in zip(m3.date1, m3.date2)])
        wet_date = np.array([bool(flags.wet.get(d, False)) for d in dates])
        pen = wet_penalty(wet_date, mem, res3)
        exceed = np.zeros_like(pen)
        for k in rng.integers(3, len(dates) - 3, N_SHIFT):
            exceed += np.abs(wet_penalty(np.roll(wet_date, k), mem, res3)) >= np.abs(pen) - 1e-12
        ppen = (exceed + 1) / (N_SHIFT + 1)
        pen[~np.isfinite(res3).all(axis=0)] = np.nan
        maps[f"wet_penalty_{track}"] = pen.reshape(shape)
        maps[f"wet_penalty_p_{track}"] = ppen.reshape(shape)
        maps[f"coh_short_mean_{track}"] = np.nanmean(flat[sel3.to_numpy()], axis=0).reshape(shape)

        # 5. monthly atlas
        mo = meta[sel3].mid.dt.month.to_numpy()
        monthly[track] = np.stack([np.nanmean(flat[sel3.to_numpy()][mo == k], axis=0).reshape(shape)
                                   if (mo == k).any() else np.full(shape, np.nan) for k in range(1, 13)])

        a_mask = zone_a.ravel() & np.isfinite(r)
        for key, rr, pp in (("coherence vs |ΔWTD|", r, p), ("wet penalty (dry − wet)", pen, ppen)):
            sig = pp[a_mask] < 0.05
            summary.append({"track": track, "test": key, "mat_pixels": int(a_mask.sum()),
                            "median_value": float(np.nanmedian(rr[a_mask])),
                            "share_p_lt_0.05": float(sig.mean()), "expected_by_chance": 0.05,
                            "share_negative" if "ΔWTD" in key else "share_positive":
                                float((rr[a_mask] < 0).mean() if "ΔWTD" in key else (rr[a_mask] > 0).mean()),
                            "n_pairs": int(sel.sum() if "ΔWTD" in key else sel3.sum())})
        q = maps[f"coh_dwtd_q_{track}"].ravel()[a_mask]
        summary[-2]["share_q_lt_0.05"] = float((q < 0.05).mean())
        del cube, flat

    # 2. VV vs WTD level (ascending RTC, 2022–2024), anomalies
    rtc = xr.open_dataset(D / "rtc_dualpol_stack.nc")
    vv = ctx.to_grid(rtc["gamma0_vv_db"])
    vh = ctx.to_grid(rtc["gamma0_vh_db"])
    times = pd.to_datetime(vv.time.values)
    lvl = []
    for t in times:
        tt = pd.Timestamp(f"{t.date()} {OVERPASS['ascending']}", tz="UTC")
        v = [field._interp_at(wtd[p], tt) for p in field.PLOTS if not bool(cens[p].asof(tt))]
        lvl.append(float(np.median(v)) if len(v) >= 5 else np.nan)
    lvl = np.array(lvl)
    keep = np.isfinite(lvl)
    yv, _ = fill_columns(vv.values.reshape(len(times), -1)[keep])
    tyr = (times[keep] - times[keep][0]).days.to_numpy() / 365.25
    dm = np.column_stack([np.ones(keep.sum()), tyr, np.cos(2 * np.pi * tyr), np.sin(2 * np.pi * tyr)])
    xv = residualize(lvl[keep][:, None], dm)[:, 0]
    rv, pv = shift_corr(xv, residualize(yv, dm))
    maps["vv_wtd_r_ascending"], maps["vv_wtd_p_ascending"] = rv.reshape(shape), pv.reshape(shape)
    am = zone_a.ravel() & np.isfinite(rv)
    summary.append({"track": "ascending", "test": "VV vs WTD level (anomalies)", "mat_pixels": int(am.sum()),
                    "median_value": float(np.nanmedian(rv[am])), "share_p_lt_0.05": float((pv[am] < 0.05).mean()),
                    "expected_by_chance": 0.05, "share_positive": float((rv[am] > 0).mean()), "n_dates": int(keep.sum())})

    # context maps for the transect and the spatial comparison
    sd = signed_distance_to_aoi(tpl, cfg).values
    ndwi = ctx.to_grid(xr.open_dataset(D / "s2_stack.nc")["ndwi"]).mean("time").values
    tcoh = ctx.to_grid(xr.open_dataset(D / "phaseE2_evd.nc")["temporal_coherence"]).values
    vh_mean = vh.mean("time").values

    # 4. transect: the plots' pixel column, edge to edge
    px = pd.read_csv(hub / "08_deliverables" / "field_first" / "plot_pixels.csv")
    col = int(px["col"].mode()[0])
    rows_a = np.flatnonzero(zone_a[:, col])
    rows = np.arange(max(rows_a.min() - 4, 0), min(rows_a.max() + 5, shape[0]))
    zlab = np.full(shape, "D", dtype=object)
    for z in "ABC":
        zlab[zones[z].values] = z
    tr = pd.DataFrame({"row": rows, "zone": zlab[rows, col], "signed_distance_m": sd[rows, col],
                       "plots": [", ".join(px.loc[px["row"] == r_, "plot"]) for r_ in rows],
                       "coh_short_mean_ascending": maps["coh_short_mean_ascending"][rows, col],
                       "coh_short_mean_descending": maps["coh_short_mean_descending"][rows, col],
                       "coh_dwtd_r_ascending": maps["coh_dwtd_r_ascending"][rows, col],
                       "coh_dwtd_r_descending": maps["coh_dwtd_r_descending"][rows, col],
                       "wet_penalty_ascending": maps["wet_penalty_ascending"][rows, col],
                       "wet_penalty_descending": maps["wet_penalty_descending"][rows, col],
                       "vv_wtd_r_ascending": maps["vv_wtd_r_ascending"][rows, col],
                       "vh_mean_db": vh_mean[rows, col], "temporal_coherence": tcoh[rows, col], "ndwi_mean": ndwi[rows, col]})
    tr.to_csv(out / "transect.csv", index=False)

    # spatial comparison over the mat (Spearman; pixels are spatially autocorrelated — descriptive)
    from scipy.stats import spearmanr
    ctx_rows = []
    for mname in ("coh_dwtd_r_ascending", "coh_dwtd_r_descending", "wet_penalty_descending", "vv_wtd_r_ascending"):
        for cname, cv in (("distance inside the mat (m)", -sd), ("S2 NDWI mean", ndwi), ("temporal coherence", tcoh), ("VH mean (dB)", vh_mean)):
            mm = zone_a & np.isfinite(maps[mname]) & np.isfinite(cv)
            rho = spearmanr(maps[mname][mm], cv[mm]).statistic
            ctx_rows.append({"map": mname, "against": cname, "rho_over_mat_pixels": float(rho), "n_pixels": int(mm.sum())})
    pd.DataFrame(ctx_rows).to_csv(out / "spatial_context.csv", index=False)
    pd.DataFrame(summary).to_csv(out / "summary.csv", index=False)
    pd.concat(arcs, ignore_index=True).to_csv(out / "arc_pairs.csv", index=False)

    ds = xr.Dataset({k: (("y", "x"), v.astype("float32")) for k, v in maps.items()}, coords={"y": tpl.y, "x": tpl.x})
    for track, cube_m in monthly.items():
        ds[f"coh_month_{track}"] = (("month", "y", "x"), cube_m.astype("float32"))
    ds = ds.assign_coords(month=np.arange(1, 13))
    ds.to_netcdf(out / "pixel_maps.nc")
    zone_summary(ds, zones).to_csv(out / "summary_by_zone.csv", index=False)
    figures(out, ds, zones, px, tr)
    print(pd.DataFrame(summary).round(3).to_string(index=False))
    print(pd.DataFrame(ctx_rows).round(3).to_string(index=False))


ZONE_TESTS = (("coherence vs |ΔWTD|", "coh_dwtd_r", "coh_dwtd_p", ("ascending", "descending")),
              ("wet penalty (dry − wet)", "wet_penalty", "wet_penalty_p", ("ascending", "descending")),
              ("VV vs WTD level (anomalies)", "vv_wtd_r", "vv_wtd_p", ("ascending",)))


def zone_summary(ds: xr.Dataset, zones) -> pd.DataFrame:
    """Every map summarised per zone (A mat, B lake, C grassland, D outside, whole scene): is the
    sensitivity specific to the mat, or regional? Pixels are spatially autocorrelated, so the
    shares are descriptive, not counts of independent tests."""
    masks = {z: zones[z].values for z in "ABCD"}
    masks["scene"] = np.ones(ds["coh_dwtd_r_ascending"].shape, bool)
    rows = []
    for test, rk, pk, tracks in ZONE_TESTS:
        for track in tracks:
            r, p = ds[f"{rk}_{track}"].values, ds[f"{pk}_{track}"].values
            for z, m in masks.items():
                mm = m & np.isfinite(r)
                rows.append({"track": track, "test": test, "zone": z, "pixels": int(mm.sum()),
                             "median_value": float(np.median(r[mm])), "share_p_lt_0.05": float((p[mm] < 0.05).mean()),
                             "share_same_sign_as_mat": float((np.sign(r[mm]) == np.sign(np.median(r[masks["A"] & np.isfinite(r)]))).mean())})
    return pd.DataFrame(rows)


def figures(out: Path, ds: xr.Dataset, zones, px: pd.DataFrame, tr: pd.DataFrame):
    outline = zones["A"].values.astype(float)
    def panel(ax, arr, title, cmap="RdBu", vlim=None):
        v = vlim or np.nanpercentile(np.abs(arr), 98)
        im = ax.imshow(arr, cmap=cmap, vmin=-v, vmax=v)
        ax.contour(outline, levels=[0.5], colors="k", linewidths=0.6)
        ax.plot(px["col"], px["row"], "k.", ms=3)
        ax.set_title(title, fontsize=9); ax.set_xticks([]); ax.set_yticks([])
        plt.colorbar(im, ax=ax, fraction=0.04)
    fig, ax = plt.subplots(2, 3, figsize=(15, 9))
    panel(ax[0, 0], ds.coh_dwtd_r_ascending.values, "coherence vs |ΔWTD|, ascending (r)")
    panel(ax[0, 1], ds.coh_dwtd_r_descending.values, "coherence vs |ΔWTD|, descending (r)")
    panel(ax[0, 2], ds.vv_wtd_r_ascending.values, "VV vs WTD level, ascending (r, anomalies)")
    panel(ax[1, 0], ds.wet_penalty_ascending.values, "wet penalty, ascending (dry − wet coherence)")
    panel(ax[1, 1], ds.wet_penalty_descending.values, "wet penalty, descending (dry − wet coherence)")
    sig = np.where(ds.coh_dwtd_q_ascending.values < 0.05, 1.0, np.where(np.isfinite(ds.coh_dwtd_q_ascending.values), 0.0, np.nan))
    panel(ax[1, 2], sig, "coherence vs |ΔWTD| asc: FDR q < 0.05 (1) / not (0)", cmap="Greys", vlim=1)
    fig.suptitle("X-053 — per-pixel sensitivity, 2020–2024 (mat outline; dots = plots)", fontsize=10)
    fig.tight_layout(); fig.savefig(out / "fig_sensitivity_maps.png", dpi=140); plt.close(fig)

    months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
    fig, ax = plt.subplots(2, 12, figsize=(22, 4.4))
    for i, track in enumerate(("ascending", "descending")):
        cube = ds[f"coh_month_{track}"].values
        for k in range(12):
            ax[i, k].imshow(cube[k], cmap="viridis", vmin=0.15, vmax=0.8)
            ax[i, k].contour(outline, levels=[0.5], colors="w", linewidths=0.4)
            ax[i, k].set_xticks([]); ax[i, k].set_yticks([])
            if i == 0:
                ax[i, k].set_title(months[k], fontsize=8)
        ax[i, 0].set_ylabel(track, fontsize=8)
    fig.suptitle("Mean short-pair coherence by month (2020–2024; mat outlined)", fontsize=10)
    fig.tight_layout(); fig.savefig(out / "fig_monthly_coherence.png", dpi=120); plt.close(fig)

    fig, ax = plt.subplots(4, 1, figsize=(10, 9), sharex=True)
    xr_ = -tr.signed_distance_m
    ax[0].plot(tr.row, tr.coh_short_mean_ascending, "o-", label="asc"); ax[0].plot(tr.row, tr.coh_short_mean_descending, "o-", label="desc")
    ax[0].set_ylabel("coherence"); ax[0].legend(fontsize=7)
    ax[1].plot(tr.row, tr.coh_dwtd_r_ascending, "o-", label="coh vs |ΔWTD| asc"); ax[1].plot(tr.row, tr.coh_dwtd_r_descending, "o-", label="desc")
    ax[1].plot(tr.row, tr.vv_wtd_r_ascending, "s--", label="VV vs WTD asc"); ax[1].axhline(0, c="k", lw=0.5); ax[1].set_ylabel("r"); ax[1].legend(fontsize=7)
    ax[2].plot(tr.row, tr.vh_mean_db, "o-", c="C3"); ax[2].set_ylabel("VH (dB)")
    ax[3].plot(tr.row, xr_, "o-", c="0.3"); ax[3].set_ylabel("distance inside\nthe mat (m)"); ax[3].set_xlabel("pixel row (north ← → south)")
    for a_ in ax:
        for r_, lab in zip(tr.row, tr.plots):
            if lab:
                a_.axvline(r_, c="0.85", lw=0.8, zorder=0)
    for r_, lab in zip(tr.row, tr.plots):
        if lab:
            ax[0].text(r_, ax[0].get_ylim()[1], lab, fontsize=7, ha="center", va="bottom")
    fig.suptitle("Transect through the plots' pixel column", fontsize=10)
    fig.tight_layout(); fig.savefig(out / "fig_transect.png", dpi=140); plt.close(fig)


if __name__ == "__main__":
    main()

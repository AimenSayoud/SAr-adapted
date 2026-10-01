"""Build the web-atlas data package from the Drive products and the repo.

    INSAR_DRIVE_ROOT=<drive or local mirror> python scripts/build_web_atlas.py --out <dir>

Reads, never writes, the data root. Everything it emits goes under `--out`:
quantised rasters on the native 40 m grid, GeoJSON outlines, chart JSON, a
figure gallery and one `manifest.json` — the format `web_export.AtlasWriter`
documents. It also runs the consistency checks that the atlas depends on and
writes their verdicts to `checks.json`, so a layer is never shown on a claim
nobody verified.

Deliberately NOT done here (CLAUDE.md): no MintPy, no phase linking, no null
realisations. Every layer below is read from a product that already exists, or
is a closed-form transformation of one (per-pixel harmonic fit, D_A, seasonal
means) with a ground-truth test in `tests/test_web_export.py`.
"""
from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import yaml

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from insar_wetlands.bootstrap import start  # noqa: E402
from insar_wetlands.web_export import (  # noqa: E402
    AtlasWriter,
    Provenance,
    amplitude_dispersion_from_db,
    labels_to_geojson,
    pair_midpoint_season,
    seasonal_fit_stack,
    seasonal_mean_maps,
    write_json,
)
from insar_wetlands.zone_viz import ZONE_COLORS  # noqa: E402

HUB = REPO.parents[1]
CHECKS: list[dict] = []


def check(name: str, ok: bool | None, detail: str, **data) -> None:
    """Record a verification. ok=None means 'informational / needs a human'."""
    verdict = "pass" if ok else ("info" if ok is None else "FAIL")
    CHECKS.append({"check": name, "verdict": verdict, "detail": detail, **data})
    print(f"  [{verdict:4s}] {name}: {detail}")


def git_sha(repo: Path) -> str:
    return subprocess.run(["git", "-C", str(repo), "rev-parse", "--short", "HEAD"],
                          capture_output=True, text=True).stdout.strip()


def zone_medians(stack: np.ndarray, zone_masks: dict) -> dict:
    """Per-date median over each zone's pixels (descriptive; NOT the project's
    double-difference aggregation, which lives in `aggregate`)."""
    out = {}
    for z, m in zone_masks.items():
        v = stack[:, m]
        with np.errstate(all="ignore"):
            med = np.nanmedian(v, axis=1)
            n = np.isfinite(v).sum(1)
        out[z] = {"median": med, "n_finite": n}
    return out


def to_dates(ds: xr.Dataset, dim: str = "time") -> list[str]:
    return [pd.Timestamp(t).strftime("%Y-%m-%d") for t in ds[dim].values]


def gis_vector_to_lonlat(path: Path) -> dict:
    """Reproject an EPSG:32633 GeoJSON (the QGIS package) to lon/lat."""
    from pyproj import Transformer
    from shapely.geometry import mapping, shape
    from shapely.ops import transform
    tr = Transformer.from_crs("EPSG:32633", "EPSG:4326", always_xy=True)
    fc = json.loads(path.read_text())
    feats = [{"type": "Feature", "properties": f["properties"],
              "geometry": mapping(transform(tr.transform, shape(f["geometry"])))}
             for f in fc["features"]]
    return {"type": "FeatureCollection", "features": feats}


def spans(flag: pd.Series) -> list[list[str]]:
    """[start, end] (inclusive, ISO dates) of the consecutive True stretches of a daily boolean series."""
    f = flag.fillna(False).astype(bool)
    grp = (f != f.shift()).cumsum()
    return [[str(g.index[0].date()), str(g.index[-1].date())] for _, g in f.groupby(grp) if g.iloc[0]]


MONTHS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]


def add_field_pixel(W, F7: Path, P, hub: Path) -> None:
    """X-053: per-pixel sensitivity maps (atlas layers) and the transect / arc / summary tables."""
    if not (F7 / "pixel_maps.nc").exists():
        return
    pm = xr.open_dataset(F7 / "pixel_maps.nc")
    prov = P(F7 / "pixel_maps.nc", root=hub)
    hours = {"ascending": "16:36", "descending": "05:09"}
    for track in ("ascending", "descending"):
        r = pm[f"coh_dwtd_r_{track}"].values.astype(float)
        p = pm[f"coh_dwtd_p_{track}"].values.astype(float)
        W.raster(f"field_coh_dwtd_r_{track}", r, title=f"Coherence vs |Δ water table|, per pixel ({track})",
                 group="Field data", status="exploratory", units="r", colormap="RdBu", symmetric=True,
                 display=(-0.4, 0.4), prov=prov,
                 description=f"Anomaly correlation, pair by pair (≤ 24 d, not frozen, {hours[track]} UTC), between each "
                             "pixel's season-cleaned coherence and the plots' median |ΔWTD|. Red = coherence falls "
                             "when the water table changes (almost everywhere: a regional effect, strongest on the mat at dusk). Circular-shift p per pixel; exploratory.")
        W.raster(f"field_coh_dwtd_r_sig_{track}", np.where(p < 0.05, r, np.nan),
                 title=f"Coherence vs |Δ water table|, p < 0.05 only ({track})", group="Field data",
                 status="exploratory", units="r", colormap="RdBu", symmetric=True, display=(-0.4, 0.4), prov=prov,
                 description="The same map keeping only pixels with a circular-shift p below 0.05 (uncorrected: about "
                             "5 % of pixels pass by chance alone).")
        W.raster(f"field_wet_penalty_{track}", pm[f"wet_penalty_{track}"].values.astype(float),
                 title=f"Wet penalty: coherence dry − wet dates ({track})", group="Field data", status="exploratory",
                 units="Δγ", colormap="RdBu_r", symmetric=True, prov=prov,
                 description="Per pixel, season-cleaned coherence of pairs with two dry dates minus pairs with a wet "
                             "date (RH ≥ 95 % or rain in the previous 3 h), 2020–2024, ≤ 24 d. Red = wetness costs "
                             "coherence here.")
        W.raster(f"field_coh_month_{track}", pm[f"coh_month_{track}"].values.astype(float),
                 title=f"Coherence by month ({track}, ≤ 24 d, 2020–2024)", group="Field data", status="derived",
                 units="γ", colormap="viridis", display=(0.2, 0.8), times=MONTHS, time_label="month", prov=prov,
                 description="Mean coherence of short pairs by the month of their midpoint, 2020–2024, frozen dates out.")
    for track in ("ascending", "descending"):
        if f"vv_wtd_r_{track}" not in pm:
            continue
        W.raster(f"field_vv_wtd_r_{track}", pm[f"vv_wtd_r_{track}"].values.astype(float),
                 title=f"VV backscatter vs water-table level, per pixel ({track})", group="Field data",
                 status="exploratory", units="r", colormap="RdBu", symmetric=True, display=(-0.5, 0.5), prov=prov,
                 description=f"Anomaly correlation over the RTC dates 2020–2024 (C-048) between each pixel's VV and the "
                             f"plots' median water-table level at {hours[track]} UTC, Sentinel-1B dates on their own level. "
                             "Blue = brighter when wetter.")
    W.chart("field_pixel", {
        "summary": pd.read_csv(F7 / "summary_by_zone.csv").round(4).to_dict("records"),
        "spatial_context": pd.read_csv(F7 / "spatial_context.csv").round(4).to_dict("records"),
        "transect": pd.read_csv(F7 / "transect.csv").round(4).astype(object).where(lambda d: d.notna(), None).to_dict("records")},
        title="Per-pixel sensitivity: summary, spatial context, transect (X-053)", group="Field data",
        status="exploratory", prov=P(F7 / "summary_by_zone.csv", F7 / "transect.csv", root=hub))
    arcs = pd.read_csv(F7 / "arc_pairs.csv")
    W.chart("field_arcs", arcs[["track", "pair", "stack", "date1", "date2", "dt_days", "wet_any", "frozen_any",
                                "dwtd_median_cm", "coh_A", "coh_C"]].round(4),
            title="Every interferogram 2020–2024 with its coherence and water-table change", group="Field data",
            status="derived", description="Zone A (mat) and C (grassland) mean coherence per pair, the plots' median "
                                          "ΔWTD, and whether a date was wet or frozen.",
            prov=P(F7 / "arc_pairs.csv", root=hub))


FUSION_ZONE_COLORS = ["#60a5fa", "#fb923c", "#22c55e", "#a855f7", "#f43f5e", "#eab308", "#14b8a6", "#78716c"]   # = fig_B_zones


def add_fusion(W, F11: Path, hub: Path, P, gallery: list, gal: Path) -> None:
    """X-062 (branch fusion-x062): the fused mat-displacement stacks and the block maps as atlas
    layers, and the block tables as one chart for the /fusion/ page. Reads the hub deliverable."""
    nc = F11 / "fusion_ascending.nc"
    if not nc.exists():
        print("  fusion: no field_fusion_x062 deliverable — skipped")
        return
    print("\n== fusion system (X-062, branch)")
    ds = xr.open_dataset(nc)
    times = [pd.Timestamp(t).strftime("%Y-%m-%d") for t in ds.time.values]
    prov = P(nc, F11 / "README.md", root=hub)
    lim = float(np.nanpercentile(np.abs(ds.fused_mm.where(ds.zone == "A").values), 98))
    W.raster("fusion_fused_ascending", ds.fused_mm.values, title="Fused mat displacement (ascending, 2022–2024)",
             group="Fusion system (X-062)", status="exploratory", units="mm", colormap="RdBu", display=(-lim, lim),
             times=times, prov=prov, x=ds.x.values, y=ds.y.values,
             description="Surface displacement since the first date, + up: the 12-day radar chain (moisture removed, "
                         "scaled with what P6 taught) fused with the water forcing; lake masked. Branch fusion-x062.")
    W.raster("fusion_radar_only_ascending", ds.radar_only_mm.where(ds.zone != "B").values,
             title="Radar-only chained displacement (ascending)", group="Fusion system (X-062)", status="exploratory",
             units="mm", colormap="RdBu", display=(-lim * 1.5, lim * 1.5), times=times, prov=prov, x=ds.x.values, y=ds.y.values,
             description="The calibrated 12-day chain alone (moisture removed, scaled by m on the mat), before fusion.")
    W.raster("fusion_amplitude_ascending", ds.amplitude_mm.values, title="Annual amplitude of the fused displacement",
             group="Fusion system (X-062)", status="exploratory", units="mm", colormap="viridis", prov=prov,
             x=ds.x.values, y=ds.y.values, description="Semi-amplitude of an annual harmonic fitted to each pixel's fused series.")
    W.raster("fusion_sd_ascending", ds.fused_sd_mm.isel(time=-1).where(ds.zone != "B").values,
             title="Fused displacement uncertainty (end of 2024)", group="Fusion system (X-062)", status="exploratory",
             units="mm", colormap="Greys", prov=prov, x=ds.x.values, y=ds.y.values,
             description="1-σ of the chained fused displacement at the last date.")
    zj = json.loads((F11 / "zones.json").read_text())
    W.raster("fusion_dissimilarity", np.log10(np.clip(ds.dissimilarity.values, 1e-3, None)), title="How unlike the plots each pixel is",
             group="Fusion system (X-062)", status="exploratory", units="log10 index", colormap="viridis_r", prov=prov,
             x=ds.x.values, y=ds.y.values,
             description=f"Dissimilarity of each pixel's 15 behaviour features to the plot pixels; above {zj['threshold']:.2f} "
                         "(log10 below) nothing measured is like it.")
    zf = F11 / "zones_map.npy"
    if zf.exists():
        zl = np.load(zf).astype(float)
        W.raster("fusion_zones", np.where(zl < 0, np.nan, zl), title="Behaviour zones (k-means)", group="Fusion system (X-062)",
                 status="exploratory", categories={i: {"label": f"zone {i}" + (" (the plots)" if i == zj["plot_zone"]["P6"] else ""), "color": FUSION_ZONE_COLORS[i % len(FUSION_ZONE_COLORS)]}
                             for i in range(zj["k"])}, colormap="tab10", prov=prov,
                 x=ds.x.values, y=ds.y.values, description="k-means zones of per-pixel behaviour; the plots all fall in one.")
    rd = lambda f: pd.read_csv(F11 / f)  # noqa: E731
    recs = lambda df: df.astype(object).where(pd.notna(df), None).to_dict("records")  # noqa: E731
    W.chart("fusion_x062", {
        "water": recs(rd("water_model_cv.csv")), "zones": zj, "zones_crosstab": recs(rd("zones_crosstab.csv")),
        "transfer": recs(rd("transfer_units.csv")), "phase_coef": recs(rd("phase_physics_coef.csv")),
        "phase_cv": recs(rd("phase_physics_cv.csv")), "buoyancy": json.loads((F11 / "buoyancy.json").read_text()),
        "validation": recs(rd("validation.csv")), "levels": recs(rd("p6_levels.csv")), "stable": recs(rd("stable_ground.csv")),
    }, title="Fusion system (X-062): block tables", group="Fusion system (X-062)", status="exploratory", prov=prov)
    for f in sorted(F11.glob("fig_*.png")):
        shutil.copy2(f, gal / f"fusion_{f.name}")
        gallery.append({"file": f"figures/fusion_{f.name}", "source": "08_deliverables/field_fusion_x062", "status": "exploratory"})


def add_field_p6(W, F10: Path, root: Path, hub: Path, P) -> None:
    """P6/CR, the primary validation site (X-058 laser verification, X-059 comparison): one chart
    file with everything the /p6/ page draws. Reads the hub's deliverables; holds nothing."""
    from insar_wetlands import field
    Q = F10 / "laser_qc"
    if not (F10 / "p6_stats.csv").exists() or not (Q / "summary.csv").exists():
        print("  field P6: no field_p6 deliverable — skipped")
        return
    rd = lambda f: pd.read_csv(f)  # noqa: E731
    recs = lambda df: df.astype(object).where(pd.notna(df), None).to_dict("records")  # noqa: E731
    fl = pd.read_csv(Q / "laser_flags.csv", parse_dates=["time_utc"]).set_index("time_utc")
    fl = fl.loc["2022-01-01":"2024-12-31 23:59"]
    meas = fl.surface_cm.where(~fl.snow_72h & ~fl.filled & ~fl.outlier)
    wtd = field.load_wtd_hourly().loc["2022-01-01":"2024-12-31 23:59"]
    day = lambda s: s.resample("D").mean()  # noqa: E731
    series = pd.DataFrame({"laser_measured_cm": day(meas), "laser_filled_cm": day(fl.surface_cm_raw.where(fl.filled)),
                           "laser_snow_cm": day(fl.surface_cm_raw.where(fl.snow_72h & ~fl.filled)),
                           "wtd_p6_cm": day(wtd["P6"]), "cr_raw_cm": day(wtd["CR_raw"])}).round(3)
    at = rd(Q / "laser_at_s1.csv")
    cov = at.groupby("track").agg(acquisitions=("date", "size"), within_1h=("laser_cm", lambda x: int(x.notna().sum())),
                                  filled=("filled", "sum"), usable=("usable", "sum"),
                                  snowfree_24h=("snow_24h", lambda x: int((~x).sum()))).reset_index()
    summ = pd.read_csv(Q / "summary.csv", index_col=0)["value"]
    data = {
        "qc": {"summary": {k: float(v) for k, v in summ.items()},
               "filled": recs(rd(Q / "interpolated.csv")), "grid": recs(rd(Q / "grid_offsets.csv")),
               "mount": recs(rd(Q / "mount_sensitivity.csv")), "years": recs(rd(Q / "coverage_by_year.csv")),
               "at_s1": recs(cov), "steps": recs(rd(Q / "steps.csv"))},
        "series": {"dates": [d.strftime("%Y-%m-%d") for d in series.index],
                   **{c: [None if pd.isna(v) else float(v) for v in series[c]] for c in series.columns}},
        "overpass": recs(at[["date", "track", "time_utc", "laser_cm", "usable", "filled", "snow_72h", "outlier"]]),
        "chain": recs(rd(F10 / "p6_web_series.csv")),
        "pairs": recs(rd(F10 / "p6_pairs.csv").drop(columns=["in_network"]).round(4)),
        "stats": recs(rd(F10 / "p6_stats.csv").round(5)),
        "seasonal": recs(rd(F10 / "p6_seasonal.csv").round(4)),
        "pair_length": recs(rd(F10 / "p6_pair_length.csv").round(4)),
        "closure": recs(rd(F10 / "p6_closure_summary.csv").round(4)),
        "triplets": recs(rd(F10 / "p6_closure_triplets.csv").round(3)),
        "controls": recs(rd(F10 / "p6_controls_summary.csv").round(4)),
        "control_pixels": recs(rd(F10 / "p6_controls_pixels.csv").round(4)),
        "checks": {"extraction": recs(rd(F10 / "p6_check_extraction.csv")), "x048": recs(rd(F10 / "p6_check_x048.csv"))},
        "params": {"quarter_wave_mm": 55.5465763 / 4, "noise_floor_mm": float(summ["noise_floor_mm"]),
                   "incidence_deg": {"ascending": 32.26, "descending": 39.17}, "window_primary": "1x1"},
    }
    W.chart("field_p6", data, title="P6/CR — the primary validation site (X-058, X-059)", group="Field data",
            status="exploratory", description="Laser verified (units, sign, 10° mount, gaps, filled stretches, "
                                                "re-levellings), then laser ↔ WTD ↔ Sentinel-1 phase ↔ coherence between "
                                                "consecutive acquisitions, pair length, closure and pixel controls.",
            prov=P(F10 / "p6_stats.csv", F10 / "p6_pairs.csv", Q / "summary.csv", root=hub))


def add_backscatter_x063(W, hub: Path, P) -> None:
    """Backscatter on both tracks, 2020–2024 (C-048 stacks, X-063 tables): per-date VV / VH, summary maps,
    and one chart file for the dusk-vs-dawn page. The stacks live in the hub's local mirror until copied to Drive."""
    from insar_wetlands.stratify import dual_pol_rvi
    mirror, F = hub / "05_code" / "local" / "drive_mirror", hub / "08_deliverables" / "backscatter_x063"
    names = {"ascending": "rtc_dualpol_2020_2024.nc", "descending": "rtc_dualpol_2020_2024_descending.nc"}
    if not all((mirror / n).exists() for n in names.values()) or not (F / "zone_backscatter_series.csv").exists():
        print("  backscatter 2020–2024: stacks or X-063 outputs missing — skipped")
        return
    print("\n== backscatter, both tracks 2020–2024 (C-048, X-063)")
    for track, name in names.items():
        ds = xr.load_dataset(mirror / name)
        prov, dd = P(mirror / name, root=hub), to_dates(ds)
        for v, t in (("gamma0_vv_db", "σ⁰ VV"), ("gamma0_vh_db", "σ⁰ VH")):
            W.raster(f"rtc5_{track}_{v}", ds[v].values.astype(float), title=f"{t} — per date ({track}, 2020–2024)",
                     group="Backscatter", status="pipeline", units="dB", colormap="gray", times=dd, prov=prov)
        for v, t in (("gamma0_vv_db", "σ⁰ VV"), ("gamma0_vh_db", "σ⁰ VH"), ("ratio_vh_vv_db", "VH − VV")):
            W.raster(f"rtc5_{track}_{v}_mean", np.nanmean(ds[v].values, 0).astype(float),
                     title=f"{t} — mean ({track}, 2020–2024)", group="Backscatter", status="supporting",
                     units="dB", colormap="gray", prov=prov)
        W.raster(f"rtc5_{track}_rvi", dual_pol_rvi(ds).values.astype(float),
                 title=f"Radar vegetation index RVI — median ({track}, 2020–2024)", group="Backscatter",
                 status="derived", colormap="YlGn", display=(0.4, 1.4), prov=prov,
                 description="RVI = 4·VH / (VV + VH) in power, per date, median over 2020–2024. Higher = more volume "
                             "(vegetation) scattering.")
        W.raster(f"rtc5_{track}_amplitude_dispersion", amplitude_dispersion_from_db(ds["gamma0_vv_db"].values),
                 title=f"Amplitude dispersion D_A, VV ({track}, 2020–2024)", group="Backscatter", status="derived",
                 colormap="magma", display=(0, 1), prov=prov,
                 description=f"σ_A/μ_A over {ds.time.size} RTC dates. < 0.25 = persistent-scatterer candidate.")
    ser = pd.read_csv(F / "zone_backscatter_series.csv")
    tab = lambda n: pd.read_csv(F / n).round(4).astype(object).where(lambda d: d.notna(), None).to_dict("records")  # noqa: E731
    W.chart("backscatter_x063", {
        "series": {track: {z: {"date": g.date.tolist(), **{v: g[v].round(3).astype(object).where(g[v].notna(), None).tolist()
                                                           for v in ("vv_db", "vh_db", "ratio_vh_vv_db", "rvi")}}
                           for z, g in st.groupby("zone")} for track, st in ser.groupby("track")},
        "summary": tab("zone_backscatter_summary.csv"), "wet_dry": tab("wet_dry_contrast.csv"),
        "satellite_offset": tab("satellite_offset_2020_2021.csv"), "lake_backscatter": tab("clean_lake_backscatter.csv"),
        "lake_coherence": tab("clean_lake_coherence.csv"), "checks": tab("reproduction_check.csv")},
        title="Backscatter by zone, both tracks, 2020–2024 (X-063)", group="Charts", status="exploratory",
        prov=P(F / "zone_backscatter_series.csv", F / "zone_backscatter_summary.csv", root=hub))


def add_field(W, gallery: list, gal: Path, template, P, drive: Path | None = None) -> None:
    """Field data (unpublished) → atlas layers, when the research hub's 06_data/field and the
    field deliverables are present next to this repository. Reads them; holds none of them."""
    from pyproj import Transformer
    from shapely.geometry import box, mapping
    from shapely.ops import transform

    from insar_wetlands import field
    root = field.field_root()
    hub = root.parents[1]
    F1, F2 = hub / "08_deliverables" / "field_first", hub / "08_deliverables" / "field_x047_x048"
    if not (root / field.DELIVERY).exists() or not (F1 / "plot_s1_wtd.csv").exists():
        print("  field: no field data or deliverables next to this repository — skipped")
        return
    print("\n== field data (hub, private)")
    to_ll = Transformer.from_crs("EPSG:32633", "EPSG:4326", always_xy=True).transform
    plots = field.load_plots()
    pts = [{"type": "Feature", "properties": {"name": r.plot, "base_plot": r.base_plot, "treatment": r.treatment,
                                              "color": "#22d3ee"},
            "geometry": mapping(transform(to_ll, r.geometry.centroid))} for r in plots.itertuples()]
    W.vector("field_plots", {"type": "FeatureCollection", "features": pts},
             title="Field monitoring plots (WTD, laser at P6)", group="Field data", status="field",
             description="9 water-table plots (P5/P6 with 3 UAV replicates; P8 moved in 2023), from the field "
                         "team's shapefile. All lie in one column of 40 m radar pixels.",
             prov=P(root / field.DELIVERY / "plots" / "Plots_PlanetScope.shp", root=hub))
    px = pd.read_csv(F1 / "plot_pixels.csv")
    xs, ys = template.x.values, template.y.values
    cells = []
    for r in px.itertuples():
        x0, y0 = xs[r.col], ys[r.row]
        cells.append({"type": "Feature", "properties": {"name": r.plot, "color": "#22d3ee",
                                                         "row": int(r.row), "col": int(r.col),
                                                         "mat_px_3x3": int(r.mat_px_3x3)},
                      "geometry": mapping(transform(to_ll, box(x0 - 20, y0 - 20, x0 + 20, y0 + 20)))})
    W.vector("field_plot_pixels", {"type": "FeatureCollection", "features": cells},
             title="Radar pixel of each field plot", group="Field data", status="derived",
             description="The 40 m Sentinel-1 pixel containing each plot (extraction windows are 3×3 around it).",
             prov=P(F1 / "plot_pixels.csv", root=hub))

    wtd = field.load_wtd_hourly()
    daily = wtd[field.PLOTS].resample("D").mean()
    cens = pd.DataFrame({p: field.censored_flag(wtd[p]).resample("D").mean() > 0.5 for p in field.PLOTS})
    W.chart("field_wtd_daily", {"dates": [d.strftime("%Y-%m-%d") for d in daily.index],
                                "plots": {p: {"wtd_cm": daily[p].round(2).tolist(), "censored": cens[p].tolist()}
                                          for p in field.PLOTS},
                                "meteo": {"rain_mm": wtd["Rain_mm_Tot"].resample("D").sum().round(2).tolist(),
                                          "air_c": wtd["Air_2m"].resample("D").mean().round(2).tolist(),
                                          "rh_pct": wtd["RH_2m"].resample("D").mean().round(1).tolist()}},
            title="Water-table depth, daily, 9 plots", group="Field data", status="field", units="cm",
            description="Daily means of the hourly series (cm, negative below the surface). `censored` = the "
                        "plot sat at its sensor floor (field.censored_flag): a bound, not a measurement.",
            prov=P(root / field.DELIVERY / "WTD_hourly_2020-2024_9plots_filled_meteo.csv", root=hub))
    laser = field.load_laser()
    snow = field.snow_mask(wtd["Air_2m"], laser.index)
    ld = pd.DataFrame({"surface_cm": laser.surface_cm.where(~snow), "raw_level_cm": wtd["CR_raw"].reindex(laser.index),
                       "wtd_p6_cm": wtd["P6"].reindex(laser.index)}).resample("D").mean()
    W.chart("field_laser_daily", {"dates": [d.strftime("%Y-%m-%d") for d in ld.index],
                                  **{c: ld[c].round(2).where(ld[c].notna(), None).tolist() for c in ld}},
            title="Laser surface at P6/CR (snow-free) with the water level", group="Field data", status="field",
            units="cm", description="SDMS40 surface position (cm, snow-masked: frost within 72 h), the raw CR "
                                    "water level and WTD_P6 corrected to the surface; daily means.",
            prov=P(root / field.DELIVERY / "Laser_Sensor.xlsx", root=hub))
    # --- what exists when (coverage timeline): acquisitions with their surface state, UAV days,
    #     laser snow-free periods, dry-well periods per plot, Sentinel-2 dates
    F3w = hub / "08_deliverables" / "field_dew_x050" / "wetness_at_overpasses_2020_2024.csv"
    cov = {"window": [str(wtd.index[0].date()), str(laser.index[-1].date())]}
    if F3w.exists():
        wt = pd.read_csv(F3w)
        wt["state"] = np.where(wt.frozen, "frozen", np.where(wt.wet, "wet", "dry"))
        cov["s1"] = wt[["date", "track", "state"]].to_dict("records")
    uav_days = field.load_uav_table()
    cov["uav"] = sorted(uav_days.Date.dt.strftime("%Y-%m-%d").unique().tolist())
    ok_hours = (laser.surface_cm.notna() & ~snow).resample("D").mean() >= 0.5
    cov["laser_snowfree"] = spans(ok_hours)
    cov["dry_wells"] = {p: spans(field.censored_flag(wtd[p]).resample("D").mean() >= 0.5) for p in field.PLOTS}
    if drive is not None and (drive / "s2_stack.nc").exists():
        cov["s2"] = [str(t)[:10] for t in xr.open_dataset(drive / "s2_stack.nc")["time"].values]
    W.chart("field_coverage", cov, title="What data exist when (2020–2025)", group="Field data", status="field",
            description="Sentinel-1 acquisitions (dry / wet / frozen at the overpass), UAV days, laser snow-free "
                        "periods, dry wells per plot, Sentinel-2 dates.", prov=P(root / field.DELIVERY / "SOURCE.md", root=hub))

    # --- the day in 24 hours (2022–2024): medians and quartiles per season × hour (UTC) of RH, air
    #     temperature, and the departure of WTD and of the snow-free laser surface from their daily mean
    h = wtd["2022":"2024"].copy()
    h["wtd_mean"] = h[field.PLOTS].where(~pd.DataFrame({q: field.censored_flag(wtd[q]) for q in field.PLOTS})).mean(axis=1)
    h["wtd_dev"] = h.wtd_mean - h.wtd_mean.groupby(h.index.floor("D")).transform("mean")
    s_ok = laser.surface_cm.where(~snow)["2022":"2024"]
    h["surface_dev"] = (s_ok - s_ok.groupby(s_ok.index.floor("D")).transform("mean")).reindex(h.index)
    season = {12: "DJF", 1: "DJF", 2: "DJF", 3: "MAM", 4: "MAM", 5: "MAM", 6: "JJA", 7: "JJA", 8: "JJA",
              9: "SON", 10: "SON", 11: "SON"}
    h["season"] = h.index.month.map(season)
    h["hour"] = h.index.hour + h.index.minute / 60   # UTC centre of the hour
    diurnal = {}
    for sname, g in h.groupby("season"):
        q = g.groupby("hour")[["RH_2m", "Air_2m", "wtd_dev", "surface_dev"]].quantile([0.25, 0.5, 0.75]).unstack()
        diurnal[sname] = {"hour": q.index.tolist(),
                          **{v: {k: q[(v, qq)].round(3).tolist() for k, qq in (("q1", 0.25), ("med", 0.5), ("q3", 0.75))}
                             for v in ("RH_2m", "Air_2m", "wtd_dev", "surface_dev")}}
    W.chart("field_diurnal", {"overpass_utc": {"ascending": 16 + 36 / 60, "descending": 5 + 9 / 60}, "seasons": diurnal},
            title="The day in 24 hours (station, WTD, laser; 2022–2024)", group="Field data", status="field",
            units="% / °C / cm", description="Per season and hour (UTC): RH, air temperature, and WTD (plot mean, dry "
                                             "wells out) and laser surface as departures from their daily mean.",
            prov=P(root / field.DELIVERY / "WTD_hourly_2020-2024_9plots_filled_meteo.csv", root=hub))

    j = pd.read_csv(F1 / "plot_s1_wtd.csv")
    keep = ["date", "track", "plot", "n_valid", "coh_pairs_le24d", "temporal_coherence", "vv_db", "vh_db", "rvi",
            "wtd_at", "mean_24h", "mean_prev3d", "mean_prev7d", "change_3d", "change_7d", "wtd_censored"]
    W.chart("field_plot_s1", j[keep].sort_values(["plot", "track", "date"]).round(4),
            title="Sentinel-1 at the plots with the WTD at each overpass", group="Field data", status="exploratory",
            description="3×3 mat-pixel medians per acquisition (field_first deliverable).",
            prov=P(F1 / "plot_s1_wtd.csv", root=hub))
    acq = pd.read_csv(F1 / "s1_acquisitions.csv")
    uav = field.load_uav_table()
    inv = {"wtd": {"start": str(wtd.index[0].ceil("D").date()), "end": str(wtd.index[-1].date()), "hours": len(wtd),
                   "plots": field.PLOTS, "censored_share": {p: round(float(field.censored_flag(wtd[p]).mean()), 4)
                                                            for p in field.PLOTS}},
           "laser": {"start": str(laser.index[0].date()), "end": str(laser.index[-1].date()),
                     "hours": len(laser), "hours_with_value": int(laser.surface_cm.notna().sum()),
                     "hours_snow_free": int((laser.surface_cm.notna() & ~snow).sum())},
           "uav": {"campaigns": sorted(uav.Date.dt.strftime("%Y-%m-%d").unique().tolist()), "rows": len(uav),
                   "subplots": int(uav.Plot.nunique()),
                   "days_with": {k: int(uav.dropna(subset=[c]).Date.nunique()) for k, c in
                                 (("multispectral", "REMX_NIR_842_mean"), ("thermal", "Altum_Thermal_11um_mean"),
                                  ("lai", "LAI"), ("wtd", "WTD_MEAN"))}},
           "acquisitions": acq.groupby("track").agg(n=("date", "size"), laser=("laser_within_1h", "sum"),
                                                    laser_snowfree=("laser_snowfree_within_1h", "sum")).reset_index()
                              .astype(object).to_dict("records"),
           "plot_pixels": px.astype(object).where(px.notna(), None).to_dict("records")}
    W.chart("field_inventory", inv, title="Field delivery inventory", group="Field data", status="field",
            prov=P(root / field.DELIVERY / "SOURCE.md", root=hub))
    if (F2 / "per_pair_correlations.csv").exists():
        pdc = pd.read_csv(F2 / "per_date_correlations.csv")
        res = {"per_pair": pd.read_csv(F2 / "per_pair_correlations.csv").to_dict("records"),
               "p6_laser": pd.read_csv(F2 / "p6_laser_vs_radar.csv").to_dict("records"),
               "p6_surface_vs_level": pd.read_csv(F2 / "p6_surface_vs_water_level.csv", index_col=0)["value"].to_dict(),
               "per_date_wtd_at": pdc[pdc.wtd == "wtd_at"].to_dict("records")}
        W.chart("field_results", res, title="Sentinel-1 vs WTD and laser (X-047, X-048)", group="Field data",
                status="exploratory", description="Anomaly correlations (annual cycle removed, circular-shift p) "
                                                  "and the P6 laser test.",
                prov=P(F2 / "per_pair_correlations.csv", F2 / "p6_laser_vs_radar.csv", root=hub))
    if (F2 / "per_pair_table.csv").exists():
        pairs = pd.read_csv(F2 / "per_pair_table.csv")
        pairs["mid"] = pairs["mid"].str[:10]
        W.chart("field_pairs", pairs[["plot", "track", "pair", "mid", "dt_days", "coh", "dlos_mm", "dwtd_cm",
                                      "dsurface_cm", "wtd_censored"]].round(3),
                title="Interferometric pairs at the plots with the WTD / surface change", group="Field data",
                status="exploratory", units="mm / cm",
                description="Coherence and LOS change per pair at each plot, with ΔWTD and, at P6, the "
                            "snow-free laser Δsurface.", prov=P(F2 / "per_pair_table.csv", root=hub))
    F3, F4 = hub / "08_deliverables" / "field_dew_x050", hub / "08_deliverables" / "field_t10_x051"
    if (F3 / "plot_phase_split_test.csv").exists():
        wt = pd.read_csv(F3 / "wetness_at_overpasses.csv")
        rh95 = lambda d: d[d.rh_threshold == 95.0].to_dict("records")  # noqa: E731 — the main threshold
        dew = {"wetness": [{"track": t, "n": len(g), "wet": int(g.wet.sum()), "frozen": int(g.frozen.sum()),
                            "wet_by_rain": int((g.rain_prev_mm > 0).sum()), "rh_median": float(g.rh.median())}
                           for t, g in wt.groupby("track")],
               "zone_coherence": rh95(pd.read_csv(F3 / "zone_coherence_by_wetness.csv")),
               "plot_phase": rh95(pd.read_csv(F3 / "plot_phase_by_wetness.csv")),
               "plot_phase_split": pd.read_csv(F3 / "plot_phase_split_test.csv").to_dict("records"),
               "series": rh95(pd.read_csv(F3 / "series_residual_by_wetness.csv"))}
        if (F3 / "zone_phase_by_wetness.csv").exists():
            zp = pd.read_csv(F3 / "zone_phase_by_wetness.csv")
            dew["zone_phase"] = zp[zp.max_dt == 24].to_dict("records")
        if (F3 / "zone_phase_dry_pairs_2020_2024.csv").exists():   # D-020: with 2020–2021
            dew["dry_pairs_2020_2024"] = pd.read_csv(F3 / "zone_phase_dry_pairs_2020_2024.csv").to_dict("records")
        W.chart("field_dew", dew,
                title="Dusk vs dawn: surface wetness at the overpass (X-050)", group="Field data",
                status="exploratory", description="Station RH / rain / frost at each overpass; coherence and "
                                                  "phase on pairs with two dry dates vs a wet one.",
                prov=P(F3 / "wetness_at_overpasses.csv", F3 / "plot_phase_split_test.csv", root=hub))
    if (F4 / "t10_proxy_vs_measured.csv").exists():
        W.chart("field_t10", {"table": pd.read_csv(F4 / "t10_proxy_vs_measured.csv").to_dict("records"),
                              "report": json.loads((F4 / "report.json").read_text())},
                title="T10 with the measured water table (X-051)", group="Field data", status="exploratory",
                description="Coherence vs |Δ water table| per zone: ERA5 precipitation proxy (committed T10) vs "
                            "the measured WTD; raw and season-cleaned.",
                prov=P(F4 / "t10_proxy_vs_measured.csv", F4 / "report.json", root=hub))
    F5, F6 = hub / "08_deliverables" / "field_methods_laser", hub / "08_deliverables" / "field_uav_x049"
    if (F5 / "methods_vs_laser.csv").exists():
        W.chart("field_methods", pd.read_csv(F5 / "methods_vs_laser.csv").round(5),
                title="The inversion outputs against the P6 laser (X-052)", group="Field data", status="exploratory",
                description="SBAS, EVD, ISBAS, hybrid and zone aggregation vs the laser surface at P6: seasonal "
                            "amplitude, series on anomalies, and pair changes.", prov=P(F5 / "methods_vs_laser.csv", root=hub))
    if (F6 / "plot_summary.csv").exists():
        W.chart("field_uav_analysis", {
            "summary": pd.read_csv(F6 / "plot_summary.csv").round(5).to_dict("records"),
            "spearman": pd.read_csv(F6 / "between_plot_spearman.csv").round(5).to_dict("records"),
            "intercorrelation": pd.read_csv(F6 / "explanatory_intercorrelation.csv").round(4).to_dict("records"),
            "two_way": pd.read_csv(F6 / "campaign_two_way.csv").round(5).to_dict("records")},
            title="UAV / LAI against the plots' radar behaviour (X-049)", group="Field data", status="exploratory",
            description="Between plots (exact Spearman, 9 plots / 8 radar units) and date-by-date with plot and "
                        "campaign means removed.", prov=P(F6 / "plot_summary.csv", F6 / "between_plot_spearman.csv", root=hub))
    F7, F8, F9 = (hub / "08_deliverables" / n for n in ("field_pixel_x053", "field_series_x054", "field_events_x055"))
    F10 = hub / "08_deliverables" / "field_p6"
    add_field_p6(W, F10, root, hub, P)
    add_field_pixel(W, F7, P, hub)
    zt_path, wt_path = hub / "08_deliverables" / "field_dew_x050" / "zone_pair_table_2020_2024.csv", F3w
    if (F7 / "arc_pairs.csv").exists() and zt_path.exists() and wt_path.exists():
        # one row per pair (both tracks, 2020–2024) for the cross-filter explorer: mat − grassland phase,
        # zone coherence, water-table / surface change, the wetter of the two overpasses
        zt = pd.read_csv(zt_path)
        zt = zt[zt.zones == "A-C"].drop_duplicates(["track", "pair"])
        arcs = pd.read_csv(F7 / "arc_pairs.csv")[["track", "pair", "coh_A", "coh_C"]]
        rh = pd.read_csv(wt_path).set_index(["track", "date"])["rh"]
        ex = zt.merge(arcs, on=["track", "pair"], how="left")
        ex["rh_max"] = [max(rh.get((t, a), np.nan), rh.get((t, b), np.nan)) for t, a, b in zip(ex.track, ex.date1, ex.date2)]
        ex["wet_any"] = ex.wet_any_rh95
        W.chart("field_explore", ex[["track", "stack", "date1", "date2", "dt_days", "dwtd_median_cm", "dsurface_p6_cm",
                                     "rh_max", "wet_any", "frozen_any", "ddisp_mm", "coh_A", "coh_C"]].round(4),
                title="Every pair 2020–2024 for the cross-filter explorer", group="Field data", status="derived",
                description="Per pair and track: mat − grassland (A − C) LOS change, zone A and C coherence, the plots' "
                            "median ΔWTD, the P6 laser Δsurface, the higher RH of the two overpasses, wet / frozen flags.",
                prov=P(zt_path, F7 / "arc_pairs.csv", wt_path, root=hub))
    if (F8 / "series_2020_2024.csv").exists():
        ser = pd.read_csv(F8 / "series_2020_2024.csv")
        W.chart("field_series5", {
            "series": ser.round(3).to_dict("records"),
            "fits": pd.read_csv(F8 / "seasonal_fits.csv").round(4).to_dict("records"),
            "vs_wtd": pd.read_csv(F8 / "series_vs_wtd.csv").round(4).to_dict("records"),
            "reproduction": pd.read_csv(F8 / "reproduction_check.csv").to_dict("records")},
            title="Aggregated zone series over five years, 2020–2024 (X-054)", group="Field data", status="exploratory",
            units="mm", description="Zone phase of every interferogram (2022–2024 + the 2020–2021 extension) inverted "
                                    "as one super-pixel per zone pair; seasonal fits over five years and year by year.",
            prov=P(F8 / "series_2020_2024.csv", F8 / "seasonal_fits.csv", root=hub))
    if (F9 / "epoch_response.csv").exists():
        ev = {"events": pd.read_csv(F9 / "events.csv").round(2).to_dict("records"),
              "response": pd.read_csv(F9 / "epoch_response.csv").round(3).to_dict("records"),
              "peaks": pd.read_csv(F9 / "peak_response.csv").round(3).to_dict("records")}
        if (F9 / "radar_event_test.csv").exists():
            ev["radar"] = pd.read_csv(F9 / "radar_event_test.csv").round(5).to_dict("records")
        W.chart("field_events", ev, title="What a rain event does (X-055)", group="Field data", status="exploratory",
                units="cm", description="Superposed epochs of rain events (≥ 10 mm / 24 h, April–October): water table "
                                        "and laser surface from 2 days before to 10 days after; radar pairs spanning "
                                        "an event against shifted event dates.",
                prov=P(F9 / "events.csv", F9 / "epoch_response.csv", root=hub))
    # The supervisor's first deliverable, whole, and every field table as a download (local site).
    files = []
    for folder in (F1, F2, F3, F4, F5, F6, F7, F8, F9, F10, F10 / "laser_qc"):
        if not folder.exists():
            continue
        for f in sorted(folder.iterdir()):
            if f.suffix.lower() in (".csv", ".md", ".json") and f.name != "laser_flags.csv":   # hourly flags: 2 MB, not a table to read
                rel = folder.relative_to(hub / "08_deliverables").as_posix()      # field_p6/laser_qc keeps its parent
                dst = W.out / "downloads" / "field" / rel / f.name
                dst.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(f, dst)
                rows = sum(1 for _ in f.open()) - 1 if f.suffix == ".csv" else None
                files.append({"folder": rel, "file": f.name, "path": f"downloads/field/{rel}/{f.name}",
                              "bytes": f.stat().st_size, "rows": rows})
    W.chart("field_deliverables", {
        "acquisitions": acq.astype(object).to_dict("records"),
        "laser_coverage": pd.read_csv(F1 / "laser_coverage.csv").to_dict("records"),
        "temporal_coherence": pd.read_csv(F1 / "s1_plot_temporal_coherence.csv").to_dict("records"),
        "files": files}, title="First deliverable for the supervisor, and field downloads", group="Field data",
        status="field", prov=P(F1 / "s1_acquisitions.csv", F1 / "laser_coverage.csv", root=hub))
    ext = pd.read_csv(F1 / "s1_plot_by_date.csv")
    W.chart("field_extraction", ext.round(5), title="Plot-level Sentinel-1 extraction (all windows)", group="Field data",
            status="derived", description="Per acquisition × plot × variable × window (1×1, 3×3, 5×5 mat pixels): n valid, "
                                          "median, mean, SD, IQR.", prov=P(F1 / "s1_plot_by_date.csv", root=hub))
    if (F2 / "per_date_correlations.csv").exists():
        W.chart("field_correlations", pd.read_csv(F2 / "per_date_correlations.csv").round(5),
                title="Sentinel-1 variables vs the six WTD variables, per plot (raw and anomalies)", group="Field data",
                status="exploratory", prov=P(F2 / "per_date_correlations.csv", root=hub))
    ndvi = (uav["REMX_NIR_842_mean"] - uav["REMX_Red_668_mean"]) / (uav["REMX_NIR_842_mean"] + uav["REMX_Red_668_mean"])
    ndre = (uav["REMX_NIR_842_mean"] - uav["REMX_Red_Edge_717_mean"]) / (uav["REMX_NIR_842_mean"] + uav["REMX_Red_Edge_717_mean"])
    uv = uav.assign(ndvi=ndvi, ndre=ndre).groupby([uav.Date.dt.strftime("%Y-%m-%d"), "base_plot"]).agg(
        ndvi=("ndvi", "median"), ndre=("ndre", "median"), thermal_c=("Altum_Thermal_11um_mean", "median"),
        lai=("LAI", "median"), wtd_mean_cm=("WTD_MEAN", "median"), subplots=("Plot", "nunique")).reset_index()
    uv.columns = ["date", "plot", *uv.columns[2:]]
    W.chart("field_uav", uv.round(4), title="UAV / LAI preview per plot and campaign (not yet analysed)", group="Field data",
            status="field", description="Medians over each plot's sub-plots of NDVI (NIR 842 / red 668), NDRE (NIR 842 / red "
                                        "edge 717), Altum surface temperature (°C) and SunScan LAI. Cells the field team "
                                        "marked as doubtful are not yet excluded.",
            prov=P(root / field.DELIVERY / "DataSet_All_RS_LAI_merged_with_WTD_Meteo.xlsx", root=hub))
    for folder, tag in ((F1, "field_first"), (F2, "field_mechanism"), (F3, "field_dew"), (F4, "field_t10"),
                        (F5, "field_methods"), (F6, "field_uav"), (F7, "field_pixel"), (F8, "field_series5"),
                        (F9, "field_events"), (F10, "field_p6"), (F10 / "laser_qc", "field_p6")):
        if not folder.exists():
            continue
        for f in sorted(folder.glob("*.png")):
            shutil.copy2(f, gal / f"{tag}_{f.name}")
            gallery.append({"file": f"figures/{tag}_{f.name}", "source": f"08_deliverables/{folder.name}",
                            "status": "exploratory"})


S2_SEASONS = {"spring": ("2023-04-01", "2023-05-31"), "summer": ("2023-06-01", "2023-08-31"),
              "autumn": ("2022-09-01", "2022-10-31"), "winter": ("2023-12-01", "2024-02-29")}


def s2_basemaps(W: AtlasWriter, x: np.ndarray, y: np.ndarray) -> None:
    """Least-cloudy Sentinel-2 L2A true-colour (TCI) crop per season, read as a
    COG window from Earth Search. Tile 33UWU is in the grid's own UTM zone and
    10 m-aligned, so the crop is an exact window — no resampling."""
    import rasterio
    from pystac_client import Client
    from rasterio.windows import from_bounds

    from insar_wetlands.web_export import grid_edges
    w_, s_, e_, n_ = grid_edges(x, y)
    from pyproj import Transformer
    tr = Transformer.from_crs("EPSG:32633", "EPSG:4326", always_xy=True)
    lo0, la0 = tr.transform(w_, s_)
    lo1, la1 = tr.transform(e_, n_)
    try:
        cat = Client.open("https://earth-search.aws.element84.com/v1")
    except Exception as e:                                          # noqa: BLE001
        check("S2 basemap", False, f"STAC unreachable: {e}")
        return
    for season, (d0, d1) in S2_SEASONS.items():
        items = list(cat.search(collections=["sentinel-2-l2a"], bbox=[lo0, la0, lo1, la1],
                                datetime=f"{d0}/{d1}", query={"eo:cloud_cover": {"lt": 10}},
                                max_items=80).items())
        items = [i for i in items if "_33UWU_" in i.id]     # the tile in the grid's UTM zone
        if not items:
            check(f"S2 basemap {season}", False, "no scene under 10 % cloud")
            continue
        it = min(items, key=lambda i: i.properties["eo:cloud_cover"])
        href = it.assets["visual"].href
        with rasterio.open(href) as src:
            if src.crs.to_epsg() != 32633:
                check(f"S2 basemap {season}", False, f"CRS {src.crs}")
                continue
            win = from_bounds(w_, s_, e_, n_, src.transform)
            rgb = src.read([1, 2, 3], window=win, boundless=True).transpose(1, 2, 0)
        W.image(f"s2_truecolour_{season}", rgb.astype("uint8"),
                title=f"Sentinel-2 true colour — {season} ({it.datetime:%Y-%m-%d})",
                group="Basemap", status="context",
                description=f"{it.id}, cloud {it.properties['eo:cloud_cover']:.2f} %. "
                            "10 m TCI, cropped to the radar grid.",
                prov=Provenance([{"file": href, "sha256": None, "stac_id": it.id}]),
                extra={"date": f"{it.datetime:%Y-%m-%d}", "cloud_cover": it.properties["eo:cloud_cover"]})
        print(f"  {season}: {it.id} ({it.properties['eo:cloud_cover']:.2f} % cloud) {rgb.shape}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--skip-pairs", action="store_true", help="skip the per-pair browser stacks")
    ap.add_argument("--attach", action="append", default=[], metavar="ID=JSON",
                    help="attach an existing JSON (e.g. a reproduction report) as a chart")
    a = ap.parse_args()
    out = Path(a.out)
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)

    ctx = start("web_atlas", mount=False, git=False)
    D = ctx.paths.drive
    cfg = ctx.cfg
    print("data root:", D)
    if "CloudStorage" in str(D) or str(D).startswith("/content"):
        print("  note: reading the LIVE Drive (read-only use)")

    tpl = ctx.template
    x, y = tpl.x.values, tpl.y.values
    W = AtlasWriter(out, x, y)
    FIG = REPO / "docs" / "paper" / "figures"      # figure images
    TAB = REPO / "results" / "tables"              # results tables (T*.csv) and series

    def P(*files, root=D):
        pv = Provenance()
        for f in files:
            pv.add(f, root=root)
        return pv

    def on_grid(da: xr.DataArray, name: str) -> np.ndarray:
        if not (np.allclose(da.x.values, x) and np.allclose(da.y.values, y)):
            raise ValueError(f"{name} is not on the template grid")
        return da.values

    # ------------------------------------------------------------------ zones
    print("\n== zones")
    zones = ctx.zones
    t01 = pd.read_csv(TAB / "T01_zones.csv").set_index("zone")
    lab = np.zeros(tpl.shape, "uint8")
    for code, z in enumerate("ABCD", start=1):
        lab[zones[z].values] = code
    counts = {z: int(zones[z].sum()) for z in "ABCD"}
    ok = all(counts[z] == int(t01.loc[z, "n_px"]) for z in "ABCD")
    check("zones reproduce T01", ok,
          f"recomputed {counts} vs T01 {t01['n_px'].to_dict()}", recomputed=counts)
    ZNAME = {"A": "Floating mat", "B": "Residual lake", "C": "Matched grassland (control)",
             "D": "Other cover"}
    zone_cats = {i: {"label": f"{z} — {ZNAME[z]}", "color": ZONE_COLORS[z]}
                 for i, z in enumerate("ABCD", start=1)}
    W.raster("zones", lab.astype(float) * np.where(lab > 0, 1, np.nan), title="Zones A–D",
             group="Zones & sampling", status="core", categories=zone_cats,
             description="Zone stratification recomputed with the pipeline's own "
                         "`define_zones` (AOI ∩ non-water = A; AOI ∩ water = B; matched "
                         "WorldCover + S2 features outside = C; rest = D).",
             prov=P(D / "water_mask.nc", D / "s2_stack.nc"))
    zfc = labels_to_geojson(lab, x, y, {i: {"zone": z, "name": ZNAME[z], "color": ZONE_COLORS[z]}
                                         for i, z in enumerate("ABCD", start=1)})
    W.vector("zones_outline", zfc, title="Zone outlines", group="Zones & sampling",
             status="core", style={"by": "zone", "colors": ZONE_COLORS})
    aoi_path = REPO / cfg["site"]["aoi_geojson"]
    aoi = json.loads(aoi_path.read_text())
    for f in aoi["features"]:        # drop the KML Z coordinate
        f["geometry"]["coordinates"] = [[c[:2] for c in ring] for ring in f["geometry"]["coordinates"]]
        f["properties"] = {"name": "Rzecin peatland (AOI)"}
    aoi.pop("bbox", None)
    W.vector("aoi", aoi, title="Rzecin peatland outline (AOI)", group="Zones & sampling",
             status="core", prov=P(aoi_path, root=REPO))
    zmask = {z: zones[z].values for z in "ABCD"}

    GIS = REPO / "data" / "gis"
    W.vector("lake_erosion_rings", gis_vector_to_lonlat(GIS / "lake_erosion_rings.geojson"),
             title="Lake erosion rings (X-010)", group="Zones & sampling",
             status="supporting", prov=P(GIS / "lake_erosion_rings.geojson", root=REPO))
    W.vector("marginal_27_pixels", gis_vector_to_lonlat(GIS / "marginal_27_pixels.geojson"),
             title="27 surviving edge pixels (X-011 / L7)", group="Zones & sampling",
             status="supporting", prov=P(GIS / "marginal_27_pixels.geojson", root=REPO))
    W.vector("field_stations", gis_vector_to_lonlat(GIS / "field_monitoring_stations.geojson"),
             title="Field stations (flux tower, piezometers, core)", group="Zones & sampling",
             status="unverified",
             description="No recorded source (ticket D-018). The flux tower agrees with the "
                         "published PolWET position (52.76 N, 16.31 E) only to that source's "
                         "2-decimal precision; P1, P2 and Core-1 sit on round 3-decimal "
                         "coordinates and look hand-placed. Core-1 is attributed to Lamentowicz "
                         "2008 here but to Milecka 2017 in the GIS README.",
             prov=P(GIS / "field_monitoring_stations.geojson", root=REPO))

    # ------------------------------------------------------- optical basemap
    print("\n== Sentinel-2 true-colour basemaps")
    s2_basemaps(W, x, y)

    # ---------------------------------------------------------------- context
    print("\n== context")
    dem = on_grid(ctx.dem, "dem").astype(float)
    W.raster("dem", dem, title="Elevation (HyP3 DEM)", group="Context", status="context",
             units="m (ellipsoidal)", colormap="terrain",
             description="Copernicus GLO-30 as delivered inside the HyP3 products, 40 m. "
                         "Heights are WGS84-ellipsoidal (HyP3 convention), so ~35–40 m "
                         "above sea level values.")
    gy, gx = np.gradient(dem)
    slope = np.degrees(np.arctan(np.hypot(gx, gy) / 40.0))
    W.raster("slope", slope, title="Slope", group="Context", status="context", units="deg",
             colormap="magma", description="From the HyP3 DEM; zones C/D exclude slopes > 5°.")
    try:
        from insar_wetlands.stratify import load_worldcover
        wc = load_worldcover(tpl, cfg, cache_dir=ctx.paths.cache)
        wcv = on_grid(wc, "worldcover").astype(float)
        WC = {10: ("Tree cover", "#006400"), 20: ("Shrubland", "#ffbb22"),
              30: ("Grassland", "#ffff4c"), 40: ("Cropland", "#f096ff"),
              50: ("Built-up", "#fa0000"), 60: ("Bare", "#b4b4b4"),
              70: ("Snow/ice", "#f0f0f0"), 80: ("Open water", "#0064c8"),
              90: ("Herbaceous wetland", "#0096a0"), 95: ("Mangroves", "#00cf75"),
              100: ("Moss/lichen", "#fae6a0")}
        present = sorted(int(v) for v in np.unique(wcv[np.isfinite(wcv)]) if int(v) in WC)
        W.raster("worldcover", np.where(np.isin(wcv, present), wcv, np.nan),
                 title="ESA WorldCover 2021", group="Context", status="context",
                 categories={k: {"label": WC[k][0], "color": WC[k][1]} for k in present},
                 description="10 m WorldCover v200, mode-resampled to the radar grid by the "
                             "pipeline's `load_worldcover`.")
    except Exception as e:                                      # noqa: BLE001
        check("WorldCover layer", False, f"not exported: {e}")
    from insar_wetlands.stratify import signed_distance_to_aoi
    sd = signed_distance_to_aoi(tpl, cfg)
    W.raster("signed_distance", on_grid(sd, "sd").astype(float),
             title="Signed distance to the peatland edge", group="Zones & sampling",
             status="supporting", units="m", colormap="RdBu", symmetric=True,
             description="Negative inside the AOI. The x-axis of the radial profiles (S8).")
    for track in ("ascending", "descending"):
        pth = ctx.paths.for_phase("web_atlas", track=track)
        from insar_wetlands.stack import load_static_layer
        th = load_static_layer(pth.cropped, "lv_theta")
        inc = np.degrees(np.pi / 2 - on_grid(th, "lv_theta").astype(float))
        inc[inc > 89] = np.nan
        W.raster(f"incidence_{track}", inc, title=f"Incidence angle ({track})",
                 group="Context", status="context", units="deg", colormap="cividis",
                 description="90° − HyP3 lv_theta (geometry.incidence_angle).")
        check(f"incidence {track}: grid median vs config", None,
              f"median {np.nanmedian(inc):.2f}° vs config "
              f"{cfg['sentinel1']['tracks'][track]['incidence_angle_deg']}°")

    # ------------------------------------------------------- radar quality
    print("\n== coherence and quality")
    cm = {}
    for track, fn in (("ascending", "phaseD_coh_mean.nc"), ("descending", "phaseD_coh_mean_descending.nc")):
        ds = xr.open_dataset(D / fn)
        cm[track] = on_grid(ds["coh_mean"], fn).astype(float)
        W.raster(f"coh_mean_{track}", cm[track], title=f"Mean coherence ({track})",
                 group="Radar quality", units="γ", colormap="viridis", display=(0.2, 0.8),
                 status="core" if track == "ascending" else "exploratory",
                 description=("Mean interferometric coherence over every pair (phaseD)."
                              if track == "ascending" else
                              "Same product on the descending (dawn) network."),
                 prov=P(D / fn))
    W.raster("coh_mean_asc_minus_desc", cm["ascending"] - cm["descending"],
             title="Coherence: dusk (asc) − dawn (desc)", group="Diurnal: dusk vs dawn",
             status="derived", units="Δγ", colormap="RdBu", symmetric=True,
             description="Difference of the two phaseD mean-coherence maps. Not pair-matched: "
                         "each track averages its own network. Descriptive — see X-037.",
             prov=P(D / "phaseD_coh_mean.nc", D / "phaseD_coh_mean_descending.nc"))
    gis_coh = REPO / "data/gis/mean_coherence_spatial.tif"
    import rasterio
    g = rasterio.open(gis_coh).read(1).astype(float)
    m = np.isfinite(g) & np.isfinite(cm["ascending"])
    check("QGIS mean_coherence_spatial.tif == phaseD_coh_mean.nc", bool(np.nanmax(np.abs(g[m] - cm["ascending"][m])) < 1e-5),
          f"max |Δ| = {np.nanmax(np.abs(g[m]-cm['ascending'][m])):.2e}")

    e2 = xr.open_dataset(D / "phaseE2_evd.nc")
    tcoh = on_grid(e2["temporal_coherence"], "tcoh").astype(float)
    W.raster("tcoh_evd", tcoh, title="Temporal coherence (EVD phase linking)",
             group="Radar quality", status="core", units="γ_t", colormap="viridis",
             display=(0.3, 1.0), description="phaseE2. 0.55 noise floor, 0.7 reliability "
             "threshold.", prov=P(D / "phaseE2_evd.nc"))
    W.raster("tcoh_usable", np.where(np.isfinite(tcoh), (tcoh >= 0.7).astype(float), np.nan),
             title="Usable pixels (γ_t ≥ 0.7)", group="Radar quality", status="core",
             categories={0: {"label": "γ_t < 0.7", "color": "#bdbdbd"},
                         1: {"label": "γ_t ≥ 0.7 (usable)", "color": "#1a9850"}},
             prov=P(D / "phaseE2_evd.nc"))
    qi = xr.open_dataset(D / "quality_index.nc")
    for v, t, u, cmap, st in (("W", "Quality index W", "", "viridis", "supporting"),
                              ("coh_all_pairs", "Coherence, all pairs (phase07)", "γ", "viridis", "pipeline"),
                              ("coh_conditional_dry", "Coherence, dry-condition pairs", "γ", "viridis", "pipeline"),
                              ("n_dry_pairs", "Number of dry pairs", "pairs", "cividis", "pipeline")):
        W.raster(f"qi_{v}", on_grid(qi[v], v).astype(float), title=t, group="Radar quality",
                 status=st, units=u, colormap=cmap, prov=P(D / "quality_index.nc"))
    clo = xr.open_dataset(D / "phase_closure.nc")
    W.raster("closure_error", on_grid(clo["closure_error_fraction"], "closure").astype(float),
             title="Closure-error fraction (518 triplets)", group="Radar quality",
             status="supporting", colormap="inferno",
             description="Fraction of closed triplets with |closure| above threshold (phase09).",
             prov=P(D / "phase_closure.nc"))

    # ------------------------------------------------------- estimators (H1)
    print("\n== per-pixel estimators (H1)")
    W.raster("evd_velocity", on_grid(e2["velocity_mm_yr"], "v").astype(float),
             title="EVD velocity", group="Per-pixel estimators (H1)", status="failed-estimator",
             units="mm/yr", colormap="RdBu_r", symmetric=True,
             description="A velocity fitted to a periodic signal is ~0 by construction; shown "
                         "because it is what a standard map would publish.",
             prov=P(D / "phaseE2_evd.nc"))
    stacks = {}
    e2d = to_dates(e2)
    stacks["evd"] = (e2["displacement_mm"].values.astype(float), e2d, "phaseE2_evd.nc", "LOS mm")
    isb = xr.open_dataset(D / "ts_isbas_ref_only.nc")
    stacks["isbas"] = (on_grid(isb["los_displacement_mm"], "isbas").astype(float), to_dates(isb),
                       "ts_isbas_ref_only.nc", "LOS mm")
    W.raster("isbas_rms", on_grid(isb["rms_residual_rad"], "rms").astype(float),
             title="ISBAS RMS residual", group="Per-pixel estimators (H1)",
             status="failed-estimator", units="rad", colormap="magma", prov=P(D / "ts_isbas_ref_only.nc"))
    W.raster("isbas_nvalid", on_grid(isb["n_valid_pairs"], "nv").astype(float),
             title="ISBAS valid pairs per pixel", group="Per-pixel estimators (H1)",
             status="failed-estimator", units="pairs", colormap="cividis", prov=P(D / "ts_isbas_ref_only.nc"))

    # hybrid: two files that look alike — decide which is which, don't guess
    hv = xr.open_dataset(D / "hybrid_vertical_ts.nc")
    hr = xr.open_dataset(D / "insar_hybrid_results.nc")
    dv = hv["vertical_mm"].values - hr["vertical_mm"].values
    same_nan = np.array_equal(np.isnan(hv["vertical_mm"].values), np.isnan(hr["vertical_mm"].values))
    check("hybrid_vertical_ts.nc vs insar_hybrid_results.nc", None,
          f"same NaN pattern: {same_nan}; max |Δ| = {np.nanmax(np.abs(dv)):.3f} mm; "
          f"median |Δ| = {np.nanmedian(np.abs(dv)):.3f} mm; "
          f"n_valid_pairs equal: {np.array_equal(hv['n_valid_pairs'].values, hr['n_valid_pairs'].values)}",
          max_abs_diff_mm=float(np.nanmax(np.abs(dv))))
    stacks["hybrid"] = (on_grid(hv["vertical_mm"], "hyb").astype(float), to_dates(hv),
                        "hybrid_vertical_ts.nc", "vertical mm")

    # SBAS (MintPy): coordinates are pixel CORNERS; shift by half a pixel and verify
    sb = xr.open_dataset(D / "ts_sbas_ref_only.nc")
    sx, sy = sb.x.values + 20.0, sb.y.values - 20.0
    aligned = np.allclose(sx, x) and np.allclose(sy, y)
    check("SBAS grid is corner-referenced (half-pixel offset)", aligned,
          f"raw x0={sb.x.values[0]}, template x0={x[0]}; after +20/−20 m shift aligned={aligned}")
    sbv = sb["los_displacement_mm"].values.astype(float)
    check("SBAS time series coverage", None,
          f"{np.isfinite(sbv).mean()*100:.1f}% finite — MintPy keeps only the reference-only mask")
    stacks["sbas"] = (sbv, to_dates(sb), "ts_sbas_ref_only.nc", "LOS mm")

    import h5py
    mp = D / "mintpy_ref_only"
    with h5py.File(mp / "velocity.h5") as h:
        at = dict(h.attrs)
        xf, yf = float(at["X_FIRST"]), float(at["Y_FIRST"])
        check("MintPy X_FIRST/Y_FIRST are corners of the template grid", bool(abs(xf + 20 - x[0]) < 1e-6 and abs(yf - 20 - y[0]) < 1e-6),
              f"X_FIRST={xf}, Y_FIRST={yf}, template centre x0={x[0]}, y0={y[0]}")
    for fn, ds, t, u, cmap in (("velocity.h5", "velocity", "MintPy SBAS velocity", "m/yr", "RdBu_r"),
                               ("temporalCoherence.h5", "temporalCoherence", "MintPy temporal coherence", "", "viridis"),
                               ("avgSpatialCoh.h5", "coherence", "MintPy average spatial coherence", "", "viridis"),
                               ("numInvIfgram.h5", "mask", "MintPy number of inverted ifgrams", "", "cividis"),
                               ("maskTempCoh.h5", "mask", "MintPy temporal-coherence mask", "", "gray")):
        try:
            with h5py.File(mp / fn) as h:
                key = ds if ds in h else list(h.keys())[0]
                arr = h[key][()].astype(float)
            if arr.ndim == 3:
                arr = arr[0]
            if fn == "velocity.h5":
                arr = arr * 1000.0
                u = "mm/yr"
            W.raster(f"mintpy_{Path(fn).stem}", arr, title=t, group="Per-pixel estimators (H1)",
                     status="failed-estimator", units=u, colormap=cmap,
                     symmetric=(cmap == "RdBu_r"), prov=P(mp / fn),
                     description="MintPy SBAS (phase08), reference-only network. "
                                 "Grid shifted by half a pixel to pixel centres.")
        except Exception as e:                                      # noqa: BLE001
            check(f"MintPy layer {fn}", False, f"not exported: {e}")

    # phase04b ensemble rates live on the outputs/phase04b branch
    for name, t in (("ensemble_vertical_rate", "Ensemble vertical rate (phase04b)"),
                    ("wls_vertical_rate", "WLS vertical rate (phase04b)")):
        blob = subprocess.run(["git", "-C", str(REPO), "show",
                               f"refs/remotes/origin/outputs/phase04b:outputs/phase04b/{name}.nc"],
                              capture_output=True)
        if blob.returncode == 0:
            tmp = out / f"_{name}.nc"
            tmp.write_bytes(blob.stdout)
            ds = xr.open_dataset(tmp)
            v = list(ds.data_vars)[0]
            arr = ds[v]
            if "x" in arr.dims and np.allclose(arr.x.values, x) and np.allclose(arr.y.values, y):
                W.raster(f"phase04b_{name}", arr.values.astype(float).squeeze(), title=t,
                         group="Per-pixel estimators (H1)", status="failed-estimator",
                         units=str(ds[v].attrs.get("units", "mm/yr")), colormap="RdBu_r", symmetric=True,
                         description="From branch outputs/phase04b (July 2026 vintage).")
            else:
                check(f"phase04b {name}", False, f"grid mismatch or dims {arr.dims}")
            ds.close()
            tmp.unlink()

    zone_series = {"zones": list("ABCD")}
    for key, (arr, dates, fn, units) in stacks.items():
        title = {"evd": "EVD displacement", "isbas": "ISBAS displacement",
                 "hybrid": "Hybrid vertical displacement (phase12)",
                 "hybrid_results": "Hybrid inversion results (phaseA/04b)",
                 "sbas": "MintPy SBAS displacement"}[key]
        src = D / fn
        W.raster(f"ts_{key}", arr, title=title + " — time series", group="Per-pixel estimators (H1)",
                 status="failed-estimator", units=units, colormap="RdBu_r", symmetric=True,
                 times=dates, prov=P(src),
                 description="Per-pixel inversion. Scrub the time slider: the fields are "
                             "dominated by per-pixel noise (H1).")
        if key in ("evd", "isbas", "hybrid"):
            fit = seasonal_fit_stack(arr, dates, min_obs=20)
            for f, u, cmap, sym in (("amplitude_mm", "mm", "magma", False),
                                    ("phase_doy", "day of year", "twilight", False),
                                    ("r2_seasonal", "R²", "viridis", False),
                                    ("trend_mm_yr", "mm/yr", "RdBu_r", True)):
                W.raster(f"fit_{key}_{f}", fit[f], title=f"{title}: per-pixel seasonal {f.split('_')[0]}",
                         group="Per-pixel seasonal fit", status="derived", units=u,
                         colormap=cmap, symmetric=sym, prov=P(src),
                         display=(0, 365) if f == "phase_doy" else ((0, 1) if f == "r2_seasonal" else None),
                         description="y = c + d·t + a·cos2πt + b·sin2πt at every pixel "
                                     "(web_export.seasonal_fit_stack ≡ aggregate.seasonal_amplitude). "
                                     "Compare with the zone aggregate: the per-pixel map is noise.")
        zm = zone_medians(arr, zmask)
        zone_series[key] = {"dates": dates, "units": units,
                            **{z: {"median": zm[z]["median"], "n": zm[z]["n_finite"]} for z in "ABCD"}}

    # ------------------------------------------------------- hydrology & optical
    print("\n== hydrology, optical, backscatter")
    wm = xr.open_dataset(D / "water_mask.nc")
    wor = wm["water_or_hidden"].values.astype(float)
    flood = ctx.flooded_fraction.values.astype(float)     # the pipeline's own definition
    W.raster("flooded_fraction", flood, title="Inundated-time fraction", group="Hydrology & optical",
             status="supporting", units="fraction", colormap="Blues", display=(0, 1),
             description="Mean of the dynamic water mask (open + hidden water), S3.",
             prov=P(D / "water_mask.nc"))
    gis_ff = rasterio.open(REPO / "data/gis/water_flooded_fraction.tif").read(1).astype(float)
    check("QGIS water_flooded_fraction.tif == pipeline flooded_fraction", bool(np.nanmax(np.abs(gis_ff - flood)) < 1e-6),
          f"max |Δ| = {np.nanmax(np.abs(gis_ff - flood)):.2e}")
    W.raster("water_stack", wm["water_or_hidden"].values.astype(float),
             title="Water mask (open + hidden) — per date", group="Hydrology & optical",
             status="pipeline", categories={0: {"label": "dry", "color": "#f0e6d2"},
                                            1: {"label": "water", "color": "#2166ac"}},
             times=to_dates(wm), prov=P(D / "water_mask.nc"))
    zone_series["water"] = {"dates": to_dates(wm), "units": "fraction",
                            **{z: {"median": np.nanmean(wor[:, zmask[z]], axis=1),
                                   "n": np.isfinite(wor[:, zmask[z]]).sum(1)} for z in "ABCD"}}
    s2 = xr.open_dataset(D / "s2_stack.nc")
    for v in ("ndwi", "mndwi"):
        arr = s2[v].values.astype(float)
        W.raster(f"{v}_stack", arr, title=f"Sentinel-2 {v.upper()} — per date",
                 group="Hydrology & optical", status="pipeline", colormap="BrBG",
                 display=(-0.6, 0.6), times=to_dates(s2), prov=P(D / "s2_stack.nc"))
        W.raster(f"{v}_mean", np.nanmean(arr, 0), title=f"Sentinel-2 {v.upper()} — mean",
                 group="Hydrology & optical", status="pipeline", colormap="BrBG",
                 display=(-0.6, 0.6), prov=P(D / "s2_stack.nc"))
        zm = zone_medians(arr, zmask)
        zone_series[v] = {"dates": to_dates(s2), "units": "index",
                          **{z: {"median": zm[z]["median"], "n": zm[z]["n_finite"]} for z in "ABCD"}}
    bc = xr.open_dataset(D / "behavior_classes.nc")["behavior_class"].values.astype(float)
    BC = {1: ("Stable ground", "#8c6d31"), 2: ("Stable vegetation", "#9ecae1"),
          3: ("Peat core", "#d62728"), 4: ("Transition / flooded", "#6baed6"),
          5: ("Permanent water", "#08306b")}
    W.raster("behaviour_classes", bc, title="Behaviour classes (phase06)", group="Hydrology & optical",
             status="pipeline", categories={k: {"label": v[0], "color": v[1]} for k, v in BC.items()},
             prov=P(D / "behavior_classes.nc"))

    rtc = xr.open_dataset(D / "rtc_dualpol_stack.nc")
    rd = to_dates(rtc)
    for v, t in (("gamma0_vv_db", "σ⁰ VV"), ("gamma0_vh_db", "σ⁰ VH"), ("ratio_vh_vv_db", "VH/VV ratio")):
        arr = rtc[v].values.astype(float)
        W.raster(f"rtc_{v}", arr, title=f"{t} — per date (ascending, 2022–2024)", group="Backscatter", status="pipeline",
                 units="dB", colormap="gray", times=rd, prov=P(D / "rtc_dualpol_stack.nc"))
        W.raster(f"rtc_{v}_mean", np.nanmean(arr, 0), title=f"{t} — mean (ascending, 2022–2024)", group="Backscatter",
                 status="supporting", units="dB", colormap="gray", prov=P(D / "rtc_dualpol_stack.nc"))
        zm = zone_medians(arr, zmask)
        zone_series[v] = {"dates": rd, "units": "dB",
                          **{z: {"median": zm[z]["median"], "n": zm[z]["n_finite"]} for z in "ABCD"}}
    W.raster("amplitude_dispersion", amplitude_dispersion_from_db(rtc["gamma0_vv_db"].values),
             title="Amplitude dispersion D_A (VV)", group="Backscatter", status="derived",
             colormap="magma", display=(0, 1),
             description="σ_A/μ_A over 85 RTC dates. < 0.25 = persistent-scatterer candidate.",
             prov=P(D / "rtc_dualpol_stack.nc"))
    W.chart("zone_series", zone_series, title="Per-zone medians through time", group="Charts",
            status="derived", description="Plain per-date zone medians of each stack. Descriptive "
            "only — the project's aggregate is the double-difference in `aggregate`, not this.")

    # ---------------------------------------------- per-pair browser & seasons
    pair_meta = {}
    for track in ("ascending", "descending"):
        print(f"\n== pairs ({track})")
        pth = ctx.paths.for_phase("web_atlas", track=track)
        from insar_wetlands.stack import list_pairs, load_layer
        pairs = list_pairs(pth.cropped)
        corr = load_layer(pth.cropped, "corr", pairs)
        on_grid(corr.isel(pair=0), "corr")
        cv = corr.values.astype(float)
        sm = seasonal_mean_maps(cv, pairs, max_dt_days=48)
        for s in ("DJF", "MAM", "JJA", "SON"):
            if s in sm:
                W.raster(f"coh_{s}_{track}", sm[s], title=f"Coherence {s} ({track}, Δt ≤ 48 d)",
                         group="Diurnal: dusk vs dawn", status="derived", units="γ",
                         colormap="viridis", display=(0.2, 0.8),
                         extra={"n_pairs": sm[f"{s}_n_pairs"]},
                         description="Mean coherence of pairs whose midpoint falls in the season.")
        pair_meta[track] = {"pairs": pairs, "seasons": sm}
        sel = pth.artifacts / "pair_selection.csv"
        ps = pd.read_csv(sel) if sel.exists() else pd.DataFrame()
        if not a.skip_pairs:
            wrapped = load_layer(pth.cropped, "wrapped_phase", pairs).values.astype(float)
            unw = load_layer(pth.cropped, "unw_phase", pairs).values.astype(float)
            W.raster(f"pairs_corr_{track}", cv, title=f"Per-pair coherence ({track})",
                     group="Pair browser", status="pipeline", units="γ", colormap="viridis",
                     display=(0, 1), times=pairs, time_label="pair", dtype="uint8",
                     prov=Provenance([{"file": f"{pth.cropped.name}/*/*_corr.tif", "sha256": None,
                                       "n_files": len(pairs)}]))
            W.raster(f"pairs_wrapped_{track}", wrapped, title=f"Per-pair wrapped phase ({track})",
                     group="Pair browser", status="pipeline", units="rad", colormap="twilight",
                     display=(-np.pi, np.pi), times=pairs, time_label="pair", dtype="uint8",
                     prov=Provenance([{"file": f"{pth.cropped.name}/*/*_wrapped_phase.tif",
                                       "sha256": None, "n_files": len(pairs)}]))
            W.raster(f"pairs_unw_{track}", unw, title=f"Per-pair unwrapped phase ({track})",
                     group="Pair browser", status="pipeline", units="rad", colormap="RdBu_r",
                     symmetric=True, times=pairs, time_label="pair",
                     prov=Provenance([{"file": f"{pth.cropped.name}/*/*_unw_phase.tif",
                                       "sha256": None, "n_files": len(pairs)}]))
            del wrapped, unw
        # pair table (network graph)
        seas = pair_midpoint_season(pairs)
        rows = []
        psi = ps.set_index("pair") if len(ps) else None
        for p, s in zip(pairs, seas):
            r0, r1 = p.split("_")
            rec = {"pair": p, "ref": r0, "sec": r1,
                   "dt_days": (pd.Timestamp(r1) - pd.Timestamp(r0)).days, "season_mid": s,
                   "coh_aoi_mean": None, "keep": None}
            if psi is not None and p in psi.index:
                rec["coh_aoi_mean"] = float(psi.loc[p, "coh_aoi_mean"])
                rec["keep"] = bool(psi.loc[p, "keep"])
                rec["season_ref"] = psi.loc[p, "season"]
            rec["coh_mean_A"] = float(np.nanmean(cv[pairs.index(p)][zmask["A"]]))
            rec["coh_mean_C"] = float(np.nanmean(cv[pairs.index(p)][zmask["C"]]))
            rows.append(rec)
        bridges = pth.artifacts / "bridge_pairs_needed.csv"
        W.chart(f"network_{track}", {"pairs": rows,
                                     "bridges_requested": (pd.read_csv(bridges).to_dict("records")
                                                           if bridges.exists() else [])},
                title=f"Interferogram network ({track})", group="Charts",
                status="supporting" if track == "ascending" else "exploratory",
                prov=P(sel) if sel.exists() else None)
        check(f"{track} network size", None,
              f"{len(pairs)} cropped pairs; {int(ps.keep.sum()) if len(ps) else 'n/a'} kept by phase03 QC")
        del cv, corr

    for s in ("DJF", "MAM", "JJA", "SON"):
        if s in pair_meta["ascending"]["seasons"] and s in pair_meta["descending"]["seasons"]:
            W.raster(f"coh_{s}_asc_minus_desc", pair_meta["ascending"]["seasons"][s] - pair_meta["descending"]["seasons"][s],
                     title=f"Coherence {s}: dusk − dawn", group="Diurnal: dusk vs dawn", status="derived",
                     units="Δγ", colormap="RdBu", symmetric=True,
                     description="Seasonal mean coherence, ascending minus descending (Δt ≤ 48 d).")

    # ------------------------------------------------------------------- charts
    print("\n== charts & consistency checks")
    committed = pd.read_csv(TAB / "phaseG_aggregate_series.csv")
    check("phaseG series: Drive == committed results copy",
          (D / "phaseG_aggregate_series.csv").read_bytes() == (TAB / "phaseG_aggregate_series.csv").read_bytes(),
          "byte comparison")
    fc = D / "figures_cache"
    cache_ac = pd.read_csv(fc / "series_AC.csv")
    same = len(cache_ac) == len(committed) and np.allclose(cache_ac.disp_mm, committed.disp_mm)
    check("figures_cache/series_AC.csv (July 25) == committed A−C series", same,
          "vintage check for series_BC/AB/NULL, which only exist in figures_cache")
    check("phaseD_coh_by_zone.csv (Drive) == figures_cache/coh_perpair.csv",
          (D / "phaseD_coh_by_zone.csv").read_bytes() == (fc / "coh_perpair.csv").read_bytes(),
          "closes X-039's 'no version of record' concern if equal")
    from insar_wetlands.aggregate import seasonal_amplitude
    t07 = pd.read_csv(TAB / "T07_seasonal_amplitudes.csv", keep_default_na=False, na_values=[""]).set_index("series")
    for s, f in (("A−C", committed), ("B−C", pd.read_csv(fc / "series_BC.csv")),
                 ("A−B", pd.read_csv(fc / "series_AB.csv")), ("NULL", pd.read_csv(fc / "series_NULL.csv"))):
        amp = seasonal_amplitude(f)["amplitude_mm"]
        check(f"T07 amplitude {s} reproduced from series", abs(amp - t07.loc[s, "amplitude_mm"]) < 0.0015,
              f"{amp:.3f} vs T07 {t07.loc[s, 'amplitude_mm']}")
    desc_series = pd.read_csv(D / "phaseG_aggregate_series_descending.csv")
    W.chart("aggregate_series", {
        "A−C (ascending, dusk)": committed.to_dict("list"),
        "B−C (ascending)": pd.read_csv(fc / "series_BC.csv").to_dict("list"),
        "A−B (ascending)": pd.read_csv(fc / "series_AB.csv").to_dict("list"),
        "Null (ascending)": pd.read_csv(fc / "series_NULL.csv").to_dict("list"),
        "A−C (descending, dawn)": desc_series.to_dict("list"),
        "fits": {k: seasonal_amplitude(v) for k, v in (
            ("A−C (ascending, dusk)", committed), ("A−C (descending, dawn)", desc_series))},
    }, title="Aggregated seasonal phase series", group="Charts", status="core",
        description="Zone double-difference series (phaseG). Descending is exploratory.",
        prov=P(TAB / "phaseG_aggregate_series.csv", root=REPO).add(D / "phaseG_aggregate_series_descending.csv", root=D))
    ref = D / "referee"
    null = pd.read_csv(ref / "LT08_null_5000.csv")
    lt08 = pd.read_csv(ref / "LT08_summary.csv")
    check("LT08 null summary matches T07 p (A−C)",
          abs(float(lt08.p_value.iloc[0]) - float(t07.loc["A−C", "p_perm"])) < 0.0006,
          f"LT08 p={lt08.p_value.iloc[0]} vs T07 p_perm={t07.loc['A−C', 'p_perm']}")
    k = int((null.amplitude_mm >= t07.loc["A−C", "amplitude_mm"]).sum())
    check("null exceedances recomputed from LT08_null_5000", None,
          f"{k} of {len(null)} draws ≥ observed; (k+1)/(N+1) = {(k+1)/(len(null)+1):.4f}")
    W.chart("null_distribution", {"amplitude_mm": null.amplitude_mm.round(4).tolist(),
                                  "summary": lt08.to_dict("records")[0]},
            title="Reference-matched null (4,614 draws) vs observed A−C amplitude",
            group="Charts", status="core", prov=P(ref / "LT08_null_5000.csv", ref / "LT08_summary.csv"))
    W.chart("null_superseded_300", pd.read_csv(fc / "nulls_300.csv"),
            title="Superseded 300-draw null (July)", group="Charts", status="superseded",
            description="Kept to show why the p-value changed; never quote.")
    for track, fn in (("ascending", "phaseD_coh_by_zone.csv"), ("descending", "phaseD_coh_by_zone_descending.csv")):
        W.chart(f"coh_by_zone_{track}", pd.read_csv(D / fn), title=f"Per-pair coherence by zone ({track})",
                group="Charts", status="core" if track == "ascending" else "exploratory", prov=P(D / fn))
    # X-037 diurnal table, recomputed from the two QC tables
    pa = pd.read_csv(D / "artifacts/pair_selection.csv")
    pdsc = pd.read_csv(D / "artifacts_descending/pair_selection.csv")
    g = lambda df: df[df.dt_days <= 48].groupby(["season", "dt_days"]).coh_aoi_mean.agg(["mean", "std", "count"])  # noqa: E731
    x37 = g(pa).join(g(pdsc), lsuffix="_asc", rsuffix="_desc").reset_index()
    x37["diff"] = x37.mean_asc - x37.mean_desc
    check("X-037: dusk > dawn in every season × Δt cell", bool((x37["diff"] > 0).all()),
          f"{int((x37['diff'] > 0).sum())}/{len(x37)} cells positive; max gap {x37['diff'].max():.3f}")
    W.chart("diurnal_coherence_table", x37.round(4), title="Coherence by season × Δt, dusk vs dawn (X-037)",
            group="Charts", status="exploratory",
            prov=P(D / "artifacts/pair_selection.csv", D / "artifacts_descending/pair_selection.csv"))
    # 2-LOS decomposition, recomputed (not copied from the ticket)
    from insar_wetlands.geometry import two_los_amplitude_uncertainty, two_los_decompose
    ep = pd.Timestamp(cfg["time"]["start"])
    fa, fd = seasonal_amplitude(committed, epoch=ep), seasonal_amplitude(desc_series, epoch=ep)
    tr = cfg["sentinel1"]["tracks"]
    ga = (tr["ascending"]["incidence_angle_deg"], tr["ascending"]["heading_deg"])
    gd = (tr["descending"]["incidence_angle_deg"], tr["descending"]["heading_deg"])
    ea, va = two_los_decompose(fa["a_cos_mm"], fd["a_cos_mm"], ga, gd)
    eb, vb = two_los_decompose(fa["b_sin_mm"], fd["b_sin_mm"], ga, gd)
    unc = two_los_amplitude_uncertainty(fa, fd, ga, gd, n_trials=20000, rng=np.random.default_rng(0))
    # Geometry actually delivered with the products: HyP3 lv_theta / lv_phi are
    # the elevation and east-referenced orientation of the ground->satellite look
    # vector. Zone-A medians -> equivalent (incidence, heading) for two_los_*.
    from insar_wetlands.stack import load_static_layer
    measured = {}
    for track in ("ascending", "descending"):
        pth = ctx.paths.for_phase("web_atlas", track=track)
        th = on_grid(load_static_layer(pth.cropped, "lv_theta"), "lv_theta").astype(float)
        ph = on_grid(load_static_layer(pth.cropped, "lv_phi"), "lv_phi").astype(float)
        th_a, ph_a = np.nanmedian(th[zmask["A"]]), np.nanmedian(ph[zmask["A"]])
        e_e, e_u = np.cos(th_a) * np.cos(ph_a), np.sin(th_a)
        inc = 90.0 - np.degrees(th_a)
        cos_h = -e_e / np.sin(np.radians(inc))
        h = np.degrees(np.arccos(np.clip(cos_h, -1, 1)))
        h = 360.0 - h   # e_E = -sin(inc)cos(h): asc ~349, desc ~190
        measured[track] = (round(float(inc), 3), round(float(h % 360), 3))
        dh = abs((h - tr[track]["heading_deg"] + 180) % 360 - 180)
        check(f"{track} geometry: config vs HyP3 look vector (Zone A median)",
              abs(inc - tr[track]["incidence_angle_deg"]) < 0.5 and dh < 1.0,
              f"incidence {inc:.2f}° vs config {tr[track]['incidence_angle_deg']}°; "
              f"heading {h % 360:.1f}° vs config {tr[track]['heading_deg']}°; "
              f"unit vector E={e_e:+.4f} U={e_u:.4f}", measured=measured[track])
    unc_meas = two_los_amplitude_uncertainty(fa, fd, measured["ascending"], measured["descending"],
                                             n_trials=20000, rng=np.random.default_rng(0))
    ea2, va2 = two_los_decompose(fa["a_cos_mm"], fd["a_cos_mm"], measured["ascending"], measured["descending"])
    eb2, vb2 = two_los_decompose(fa["b_sin_mm"], fd["b_sin_mm"], measured["ascending"], measured["descending"])
    check("2-LOS with measured geometry (sensitivity)", None,
          f"vertical {np.hypot(va2, vb2):.3f} mm (CI {unc_meas['vertical_amplitude_mm']['ci95']}), "
          f"east {np.hypot(ea2, eb2):.3f} mm (CI {unc_meas['east_amplitude_mm']['ci95']}) "
          f"vs config geometry vertical {np.hypot(va, vb):.3f}, east {np.hypot(ea, eb):.3f}")
    j16 = json.loads((D / "phase16_two_los_decomposition.json").read_text())
    check("phase16 point estimate reproduced", abs(np.hypot(va, vb) - j16["vertical_amplitude_mm"]) < 0.002,
          f"vertical {np.hypot(va, vb):.3f} vs Drive JSON {j16['vertical_amplitude_mm']}; "
          f"east {np.hypot(ea, eb):.3f} vs {j16['east_amplitude_mm']} — a mismatch means Drive's "
          f"phase16 JSON was computed with other geometry than config.yaml (X-042: re-run phase16 on Colab)")
    check("phase16 uncertainty (X-032 CI) recomputed", None,
          f"vertical {unc['vertical_amplitude_mm']}; east {unc['east_amplitude_mm']}")
    # The decomposition as first run on Colab (X-032): its own recorded geometry, uncertainty
    # recomputed here so the page can set it beside the current (measured, X-042) geometry.
    g_run = (tuple(j16["geometry"]["ascending"]), tuple(j16["geometry"]["descending"]))
    unc_run = two_los_amplitude_uncertainty(fa, fd, *g_run, n_trials=20000, rng=np.random.default_rng(0))
    W.chart("two_los", {"fit_ascending": fa, "fit_descending": fd,
                        "first_run": {"geometry": {"ascending": g_run[0], "descending": g_run[1]},
                                      "vertical_amplitude_mm": j16["vertical_amplitude_mm"],
                                      "east_amplitude_mm": j16["east_amplitude_mm"],
                                      "uncertainty": unc_run},
                        "vertical_amplitude_mm": float(np.hypot(va, vb)),
                        "east_amplitude_mm": float(np.hypot(ea, eb)),
                        "uncertainty": unc, "geometry": {"ascending": ga, "descending": gd},
                        "measured_geometry": measured,
                        "measured_geometry_result": {
                            "vertical_amplitude_mm": float(np.hypot(va2, vb2)),
                            "east_amplitude_mm": float(np.hypot(ea2, eb2)),
                            "uncertainty": unc_meas},
                        "pure_vertical_ascending_only_mm": fa["amplitude_mm"] / np.cos(np.radians(ga[0]))},
            title="2-LOS decomposition (X-032)", group="Charts", status="exploratory",
            prov=P(D / "phase16_two_los_decomposition.json"))
    # ERA5 at the site, as the pipeline reads it for anything precipitation-based (X-043)
    from insar_wetlands.hydro import era5_samples_per_day, open_era5
    era = open_era5(D)
    lon, lat = cfg["site"]["centroid"]
    pt = era.sel(latitude=lat, longitude=lon, method="nearest")
    tname = "valid_time" if "valid_time" in pt.coords else "time"
    df = pd.DataFrame({"t2m_c": pt["t2m"].values - 273.15, "tcwv": pt["tcwv"].values,
                       "precip_mm": pt["tp"].values * 1000.0},
                      index=pd.to_datetime(pt[tname].values))
    spd = era5_samples_per_day(era)
    check("ERA5 precipitation sampling", abs(spd - 24) < 0.5,
          f"{spd:.0f} samples/day: tp is a 1-hour accumulation, so a daily sum needs 24 "
          f"(hydro.open_era5 reads era5_hourly/ when present)")
    daily = df.resample("1D").agg({"t2m_c": "mean", "tcwv": "mean", "precip_mm": "sum"})
    daily.index = daily.index.strftime("%Y-%m-%d")
    hourly_dir = D / "era5_hourly"
    W.chart("era5_daily", daily.reset_index(names="date").round(3), title="ERA5 at the site (daily)",
            group="Charts", status="pipeline",
            description="Daily means of t2m and tcwv and daily precipitation summed over all 24 "
                        "hourly accumulations (era5_hourly/, X-043).",
            prov=P(*sorted(hourly_dir.glob("*.nc"))) if hourly_dir.is_dir() else P(D / "era5_rzecin.nc"))
    # every committed table
    tables = {}
    for f in sorted(TAB.glob("T*.csv")):
        tables[f.stem] = pd.read_csv(f, keep_default_na=False, na_values=[""]).astype(object).where(lambda d: pd.notna(d), None).to_dict("records")
    W.chart("results_tables", tables, title="Results tables T01–T16", group="Charts", status="core",
            description="Verbatim copies of results/tables/T*.csv — the numbers' source of truth.")
    ref_tables = {f.stem: pd.read_csv(f, keep_default_na=False, na_values=[""]).astype(object).where(lambda d: pd.notna(d), None).to_dict("records")
                  for f in sorted(ref.glob("*.csv")) if f.stat().st_size < 20000}
    W.chart("robustness_tables", ref_tables, title="Robustness tables (K*, L*, X*)", group="Charts",
            status="supporting")
    pm = REPO / "outputs/phaseM_mechanical_vs_dielectric/mechanical_vs_dielectric_summary.csv"
    if pm.exists():
        W.chart("phaseM_summary", pd.read_csv(pm), title="phaseM mechanical vs dielectric (exploratory)",
                group="Charts", status="exploratory",
                description="Part II of phaseM is synthetic-only. X-027 envelope discrepancy open.",
                prov=P(pm, root=REPO))
    from insar_wetlands.paper_numbers import expected_values
    nums = expected_values(TAB)
    W.chart("numbers", nums, title="Key result numbers", group="Charts", status="core",
            description="paper_numbers.REGISTRY resolved against T*.csv. The only result numbers "
                        "the site prints in text.")
    phases = yaml.safe_load((REPO / "config/phases.yaml").read_text())
    W.chart("phases", phases, title="Declared pipeline (config/phases.yaml)", group="Pipeline", status="core")
    runs = []
    for mf in sorted((D / "runs").glob("*/*/manifest.json")):
        try:
            j = json.loads(mf.read_text())
        except Exception:                                   # noqa: BLE001
            continue
        runs.append({"phase": mf.parts[-3], "run": mf.parts[-2],
                     "git_sha": j.get("git", {}).get("sha") if isinstance(j.get("git"), dict) else j.get("git_sha"),
                     "dirty": j.get("git", {}).get("dirty") if isinstance(j.get("git"), dict) else None,
                     "utc": mf.parts[-2].split("_")[0], "products": list((j.get("products") or {}).keys())})
    W.chart("runs", runs, title="Archived Colab executions", group="Pipeline", status="pipeline")
    tickets = []
    for f in sorted((HUB / "_ledger").glob("[A-Z]-[0-9]*.md")):
        mt = re.match(r"^---\n(.*?)\n---", f.read_text(), re.S)
        if mt:
            try:
                fm = yaml.safe_load(mt.group(1))
            except Exception:                               # noqa: BLE001
                continue
            tickets.append({k: fm.get(k) for k in ("id", "title", "status", "priority", "progress", "start",
                                                   "due", "work_chain", "tags", "dependencies", "sha")})
    W.chart("tickets", tickets, title="Ledger tickets", group="Pipeline", status="pipeline")

    # ------------------------------------------------------------ figure gallery
    gal = out / "figures"
    gal.mkdir()
    listed = set(re.findall(r"`([FS]\d+_[a-z_]+\.png)`", (REPO / "docs/paper/figures_tables.md").read_text()))
    gallery = []
    for f in sorted(FIG.glob("*.png")):
        shutil.copy2(f, gal / f.name)
        gallery.append({"file": f"figures/{f.name}", "source": "docs/paper/figures",
                        "status": "core" if f.name in listed else "superseded"})
    for f in sorted(ref.glob("*.png")):
        shutil.copy2(f, gal / f"robustness_{f.name}")
        gallery.append({"file": f"figures/robustness_{f.name}", "source": "Drive referee/", "status": "supporting"})
    br = subprocess.run(["git", "-C", str(REPO), "for-each-ref", "--format=%(refname)",
                         "refs/remotes/origin/outputs/"], capture_output=True, text=True).stdout.split()
    for b in br:
        files = subprocess.run(["git", "-C", str(REPO), "ls-tree", "-r", "--name-only", b],
                               capture_output=True, text=True).stdout.split()
        for fn in files:
            if fn.startswith("outputs/") and fn.endswith(".png"):
                blob = subprocess.run(["git", "-C", str(REPO), "show", f"{b}:{fn}"], capture_output=True)
                name = "branch_" + fn.replace("outputs/", "").replace("/", "__")
                (gal / name).write_bytes(blob.stdout)
                gallery.append({"file": f"figures/{name}", "source": b.replace("refs/remotes/origin/", ""),
                                "status": "pipeline"})
    if (REPO / "outputs/phaseM_mechanical_vs_dielectric/mechanical_vs_dielectric_dashboard.png").exists():
        shutil.copy2(REPO / "outputs/phaseM_mechanical_vs_dielectric/mechanical_vs_dielectric_dashboard.png",
                     gal / "phaseM_dashboard.png")
        gallery.append({"file": "figures/phaseM_dashboard.png", "source": "phaseM", "status": "exploratory"})
    add_field(W, gallery, gal, tpl, P, drive=D)
    from insar_wetlands import field as _field
    hub_root = _field.field_root().parents[1]
    add_fusion(W, hub_root / "08_deliverables" / "field_fusion_x062", hub_root, P, gallery, gal)
    add_backscatter_x063(W, hub_root, P)
    W.chart("gallery", gallery, title="Figure gallery", group="Figures", status="core")

    for spec in a.attach:
        lid, path = spec.split("=", 1)
        W.chart(lid, json.loads(Path(path).read_text()), title=lid.replace("_", " "),
                group="Pipeline", status="pipeline", description=f"Attached from {Path(path).name}.")

    # 256-entry colour tables from matplotlib for every colormap a layer names, so the
    # site colours pixels exactly as the project figures do (no JS re-implementation).
    import matplotlib
    luts = {}
    for name in sorted({lyr["colormap"] for lyr in W.layers if lyr.get("colormap")}):
        cm = matplotlib.colormaps[name]
        luts[name] = [[int(round(c * 255)) for c in cm(i / 255.0)[:3]] for i in range(256)]
    write_json(out / "colormaps.json", luts)

    write_json(out / "checks.json", CHECKS, indent=1)
    W.write_manifest(title="Rzecin InSAR atlas", git_sha=git_sha(REPO),
                     built_utc=pd.Timestamp.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ"),
                     data_root=str(D), zone_colors=ZONE_COLORS,
                     site={"centroid": cfg["site"]["centroid"], "name": cfg["site"]["name"]},
                     tracks=cfg["sentinel1"]["tracks"])
    total = sum(f.stat().st_size for f in out.rglob("*") if f.is_file())
    print(f"\n{len(W.layers)} layers, {total/1e6:.1f} MB -> {out}")
    fails = [c for c in CHECKS if c["verdict"] == "FAIL"]
    print(f"checks: {len(CHECKS)} ({len(fails)} FAIL)")


if __name__ == "__main__":
    main()

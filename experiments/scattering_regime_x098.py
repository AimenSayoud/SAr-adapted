"""X-098 — scattering regime on the mat: does the consecutive-pair phase read the water level or the surface?

Where water stands among emergent vegetation the return is double bounce and the phase reads the water level
(wdowinski2007, kim2017, fujiwara2019); where it is below the surface the phase reads the surface. Each acquisition of
each pixel is classed from Sentinel-1 backscatter as tsyganskaya2018 does — the Z-score of VV and of the VV/VH ratio
over the pixel's own series (double bounce: both up; open water: VV down) — and the pair phase at the logger plots is
regressed on the water-table change by class. Reading the water level gives about 10·cos θ mm of LOS per cm; reading
the floating surface about g·cos θ (buoyancy, T47), an order of magnitude less.

Pre-registration: ``08_deliverables/insar_full_test/preregistration/X-098.yaml`` (rule and thresholds fixed there).
Reads the cube; writes only to the hub: ``08_deliverables/insar_full_test/scattering_regime/``.

    make run E=scattering_regime_x098
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import p6_validation_x082 as X  # noqa: E402  (consecutive-pair filter, PHASE_TO_MM)

from insar_wetlands.cube import cube_root  # noqa: E402
from insar_wetlands.paths import hub_root  # noqa: E402
from insar_wetlands.scoring import assert_preregistered, increment_stats  # noqa: E402

HUB = hub_root()
CUBE = cube_root()
FT = HUB / "08_deliverables" / "insar_full_test"
OUT = FT / "scattering_regime"
PREREG = FT / "preregistration" / "X-098.yaml"
K_MAIN, K_SENS = 1.5, (1.0, 2.0)
YEARS = (2020, 2021, 2022, 2023, 2024)
CLASSES = ("double_bounce", "surface", "open_water")
WATER_SLOPE, SURFACE_SLOPE, R_MIN, N_MIN = 5.0, 2.5, 0.3, 15
CODE = {0: "surface", 1: "double_bounce", 2: "open_water"}


def regime(track: str, k: float) -> xr.DataArray:
    """Class code per (time, y, x): 2 open-water-like (Z_VV <= -k), 1 double-bounce-like (Z_ratio >= k, Z_VV >= 0),
    0 surface; NaN where backscatter is missing."""
    d = xr.open_dataset(CUBE / "silver" / f"rtc_{track}.nc")
    vv = d.vv_db
    ratio = d.vv_db - d.vh_db
    zv = (vv - vv.mean("time")) / vv.std("time")
    zr = (ratio - ratio.mean("time")) / ratio.std("time")
    c = xr.where(zv <= -k, 2, xr.where((zr >= k) & (zv >= 0), 1, 0)).astype("float32")
    c = c.where(vv.notnull() & d.vh_db.notnull())
    c = c.assign_coords(date=("time", pd.to_datetime(d.time.values).normalize()))
    return c


def pairs() -> pd.DataFrame:
    cols = ["track", "grid", "plot", "row", "col", "pair", "t1", "t2", "revisit_days", "unw_1x1", "ref_hard_unw",
            "d_wtd_plot_cm", "d_laser_los_mm"]
    p = pd.read_csv(CUBE / "gold" / "plots_pair.csv.gz", usecols=cols, low_memory=False)
    p = p[(p.grid == "40m") & p.revisit_days.isin([6, 12])].copy()
    for c in ("t1", "t2"):
        p[c] = pd.to_datetime(p[c], format="ISO8601", utc=True)
    p = p[X.consecutive(p) & p.t1.dt.year.isin(YEARS)].reset_index(drop=True)
    p["inc_mm"] = (p.unw_1x1 - p.ref_hard_unw) * X.PHASE_TO_MM
    return p


def label(p: pd.DataFrame, cls: dict[str, xr.DataArray]) -> pd.Series:
    out = pd.Series("missing", index=p.index)
    for track, c in cls.items():
        dates = pd.DatetimeIndex(c.date.values)
        pos = {d.date(): i for i, d in enumerate(dates)}
        arr = c.values
        s = p.track == track
        for i in p.index[s]:
            a, b = pos.get(p.at[i, "t1"].date()), pos.get(p.at[i, "t2"].date())
            if a is None or b is None:
                continue
            r, cc = int(p.at[i, "row"]), int(p.at[i, "col"])
            ca, cb = arr[a, r, cc], arr[b, r, cc]
            if np.isnan(ca) or np.isnan(cb):
                continue
            out.at[i] = CODE[int(ca)] if ca == cb else "transition"
    return out


def slope_r(x: np.ndarray, y: np.ndarray) -> tuple[float, float, int]:
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if len(x) < 3 or not np.any(x):
        return np.nan, np.nan, int(len(x))
    return float(np.sum(x * y) / np.sum(x * x)), float(np.corrcoef(x, y)[0, 1]), int(len(x))


def read(p: pd.DataFrame, cosi: dict[str, float]) -> pd.DataFrame:
    rows = []
    for (track, cls), g in p.groupby(["track", "regime"]):
        if cls == "missing":
            continue
        sl, r, n = slope_r(g.d_wtd_plot_cm.to_numpy(float), g.inc_mm.to_numpy(float))
        p6 = g[(g["plot"] == "P6") & g.d_laser_los_mm.notna() & g.t1.dt.year.between(2021, 2024)]
        las = increment_stats(p6.d_laser_los_mm, p6.inc_mm) if len(p6) >= 3 else {}
        r6 = las.get("pearson_increments", np.nan)
        if n < N_MIN:
            reading = "too few"
        elif sl >= WATER_SLOPE and r >= R_MIN:
            reading = "water level"
        elif sl <= SURFACE_SLOPE and np.isfinite(r6) and r6 >= R_MIN:
            reading = "surface"
        else:
            reading = "neither"
        rows.append({"track": track, "regime": cls, "pairs": n, "slope_mm_per_cm": sl, "r_vs_wtd": r,
                     "expected_water_level_mm_per_cm": 10 * cosi[track], "p6_pairs": len(p6),
                     "p6_r_vs_laser": r6, "reading": reading})
    return pd.DataFrame(rows)


def main() -> None:
    assert_preregistered(PREREG, [OUT / "readings.csv"])
    OUT.mkdir(parents=True, exist_ok=True)
    st = xr.open_dataset(CUBE / "silver" / "static_40m.nc")
    cosi = {"ascending": float(np.cos(np.radians(np.nanmedian(st.incidence_asc_deg.values)))),
            "descending": float(np.cos(np.radians(np.nanmedian(st.incidence_des_deg.values))))}
    bw = st.boardwalk_share.values > 0
    base = pairs()
    base["boardwalk_cell"] = [bool(bw[int(r), int(c)]) for r, c in zip(base.row, base.col)]
    readings, monthly = [], []
    maps = {}
    for k in (K_MAIN, *K_SENS):
        cls = {t: regime(t, k) for t in X.TRACKS}
        p = base.assign(regime=label(base, cls))
        for subset, q in (("all plot cells", p), ("without boardwalk cells", p[~p.boardwalk_cell])):
            readings.append(read(q, cosi).assign(k=k, cells=subset))
        if k == K_MAIN:
            p.to_csv(OUT / "pairs.csv", index=False)
            mat = st.zone.values == 1
            for t, c in cls.items():
                v = c.values[:, mat]
                months = pd.DatetimeIndex(c.date.values).month
                for m in range(1, 13):
                    sel = v[months == m]
                    if sel.size:
                        monthly.append({"track": t, "month": m, **{CODE[j]: float(np.nanmean(sel == j))
                                                                    for j in CODE}})
                ch = (np.diff(c.values, axis=0) != 0) & np.isfinite(c.values[1:]) & np.isfinite(c.values[:-1])
                maps[f"class_change_rate_{t}"] = ch.mean(axis=0)
                for j, name in CODE.items():
                    maps[f"share_{name}_{t}"] = np.nanmean(c.values == j, axis=0)
    readings = pd.concat(readings, ignore_index=True)
    monthly = pd.DataFrame(monthly)
    readings.to_csv(OUT / "readings.csv", index=False)
    monthly.to_csv(OUT / "mat_class_by_month.csv", index=False)
    xr.Dataset({k: (("y", "x"), v.astype("float32")) for k, v in maps.items()},
               coords={"y": st.y, "x": st.x}).to_netcdf(OUT / "regime_maps.nc")
    (OUT / "README.md").write_text(readme(readings, monthly), encoding="utf-8")
    print(f"wrote {OUT}")


def _md(df: pd.DataFrame, nd: int = 2) -> str:
    d = df.copy()
    for c in d.columns:
        if pd.api.types.is_float_dtype(d[c]):
            d[c] = d[c].map(lambda x: "—" if pd.isna(x) else f"{x:.{nd}f}")
    return "\n".join(["| " + " | ".join(map(str, d.columns)) + " |", "|" + "---|" * len(d.columns)] +
                     ["| " + " | ".join(map(str, r)) + " |" for r in d.itertuples(index=False)])


def readme(readings: pd.DataFrame, monthly: pd.DataFrame) -> str:
    main = readings[(readings.k == K_MAIN) & (readings.cells == "all plot cells")].drop(columns=["k", "cells"])
    sens = readings[readings.regime.isin(CLASSES)].pivot_table(index=["track", "regime", "cells"], columns="k",
                                                                values="slope_mm_per_cm").reset_index()
    return f"""# X-098 — scattering regime and what the phase reads (generated)

Generated by `experiments/scattering_regime_x098.py` (`make run E=scattering_regime_x098`); do not edit.
Pre-registration: `../preregistration/X-098.yaml`. Classes per date from the Z-score of VV and VV/VH over each pixel's
series (tsyganskaya2018), k = {K_MAIN} (our choice; the paper's thresholds came from a reference Rzecin lacks).
Pairs: consecutive 6/12-day, plots P1–P9, 2020–2024, hard-target reference. Reading rule: "water level" if slope >=
{WATER_SLOPE} mm/cm and r >= {R_MIN}; "surface" if slope <= {SURFACE_SLOPE} mm/cm and P6 increments follow the laser
(r >= {R_MIN}); else "neither"; "too few" under {N_MIN} pairs. P6 laser part tried on 2021–2024 only (D-028).

## Readings (k = {K_MAIN}, all plot cells)

{_md(main)}

## Slope (mm per cm) by threshold k and cell set

{_md(sens)}

## Mat pixels by class and month (share of acquisitions)

{_md(monthly)}

Files: `readings.csv` (every k and cell set), `pairs.csv` (each pair with its regime), `mat_class_by_month.csv`,
`regime_maps.nc` (per pixel: share of acquisitions in each class, rate of class change).
"""


if __name__ == "__main__":
    main()

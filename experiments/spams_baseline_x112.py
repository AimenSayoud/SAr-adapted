"""X-112 — SPAMS as a published baseline for the floating mat's prior, at P6, year-blocked.

conroy2023_spams: soft-soil surface height from precipitation minus scaled evapotranspiration over a memory time.
Compared, on the same days and inputs as fusion v3 WP1 (X-074), with buoyancy on the logger water table, the ERA5
bucket water table × buoyancy, and zero change. Rule: ``08_deliverables/insar_full_test/preregistration/X-112.yaml``.
Writes only to the hub: ``08_deliverables/insar_full_test/spams_baseline/``.

    make run E=spams_baseline_x112
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fusion_v3_wp1_x074 as W1  # noqa: E402  (daily table and frozen-window rule of X-074)

from insar_wetlands import fusion as fu  # noqa: E402
from insar_wetlands import spams  # noqa: E402
from insar_wetlands.paths import cube_root, hub_root  # noqa: E402
from insar_wetlands.scoring import Protocol, assert_preregistered, detrended_stats, guard_lockbox  # noqa: E402

HUB = hub_root()
FT = HUB / "08_deliverables" / "insar_full_test"
OUT = FT / "spams_baseline"
PREREG = FT / "preregistration" / "X-112.yaml"
YEARS = (2021, 2022, 2023, 2024, 2025)
LAG = 6
LOGGER_YEARS = (2021, 2022, 2023, 2024)


def table() -> pd.DataFrame:
    d = W1.daily()
    h = pd.read_csv(cube_root() / "silver" / "site_hourly.csv.gz", index_col=0, parse_dates=True, low_memory=False,
                    usecols=["time_utc", "om_era5_et0_fao_evapotranspiration"])
    h.index = pd.DatetimeIndex(h.index, tz="UTC") if h.index.tz is None else h.index
    et0 = h.om_era5_et0_fao_evapotranspiration.resample("1D").sum(min_count=20)
    et0.index = et0.index.tz_convert(None).normalize()
    d["et0_mm"] = et0.reindex(d.index)
    d["pet_hamon"] = fu.hamon_pet(d.air_era5_c.to_numpy(), d.doy.to_numpy(), W1.LAT)
    d["dh_mm"] = (d.laser_cm.shift(-LAG) - d.laser_cm) * 10
    d["dw_cm"] = d.wtd_cm.shift(-LAG) - d.wtd_cm
    d["frozen_win"] = d.frozen.rolling(LAG + 1).max().shift(-LAG).fillna(1) > 0
    return d


def bucket_heldout(d: pd.DataFrame, year: int) -> pd.Series:
    """Bucket water table for ``year``, fitted on the logger of the other years, started from the last logger value."""
    rain, pet, w = d.rain_era5_mm.fillna(0).to_numpy(), d.pet_hamon.to_numpy(), d.wtd_cm.to_numpy()
    other = (d.year != year).to_numpy()
    par = fu.fit_bucket(np.where(other, w, np.nan), rain, pet)
    floor = np.nanpercentile(w[other], 1)
    m = (d.year == year).to_numpy()
    i0 = np.flatnonzero(m)[0]
    prev = w[:i0][np.isfinite(w[:i0])]
    w0 = prev[-1] if len(prev) else par["w0"]
    out = pd.Series(np.nan, index=d.index)
    out.iloc[np.flatnonzero(m)] = fu.simulate_bucket(par, w0, rain[m], pet[m], floor=floor)
    return out


def g_fit(d: pd.DataFrame, train: np.ndarray) -> float:
    ok = train & np.isfinite(d.dh_mm.to_numpy()) & np.isfinite(d.dw_cm.to_numpy())
    x, y = d.dw_cm.to_numpy()[ok], d.dh_mm.to_numpy()[ok]
    return float(np.sum(x * y) / np.sum(x * x))


def main() -> None:
    assert_preregistered(PREREG, [OUT / "verdict.csv"])
    OUT.mkdir(parents=True, exist_ok=True)
    prot = Protocol.load(FT / "protocol.yaml")
    d = table()
    d = d[d.year.isin(YEARS + (2020,))]                  # 2020: logger history for the bucket and buoyancy fits only
    guard_lockbox(d.index[d.year.isin(YEARS)], prot)
    P, E = d.rain_era5_mm.to_numpy(), d.et0_mm.to_numpy()
    valid = (d.dh_mm.notna() & ~d.frozen_win & d.year.isin(YEARS)).to_numpy()
    pred = {k: pd.Series(np.nan, index=d.index) for k in ("zero change", "SPAMS", "SPAMS + irreversible",
                                                          "bucket × g", "logger water table × g")}
    lev = {k: pd.Series(np.nan, index=d.index) for k in ("SPAMS", "bucket × g", "logger water table × g")}
    folds = []
    for y in YEARS:
        te = (d.year == y).to_numpy()
        tr = valid & ~te
        par = spams.fit(P, E, d.dh_mm.to_numpy(), tr, LAG)
        x_i = spams.fit_irreversible(P, E, par, d.dh_mm.to_numpy(), tr, LAG)
        ch = spams.predict_change(P, E, par, LAG)
        R = spams.reversible(P, E, par["x_p"], par["x_e"], par["tau"])
        I = spams.irreversible(R, x_i)
        dI = np.full(len(I), np.nan)
        dI[:-LAG] = I[LAG:] - I[:-LAG]
        g = g_fit(d, tr & d.year.isin(LOGGER_YEARS + (2020,)).to_numpy())
        b = bucket_heldout(d, y)
        pred["zero change"][te] = 0.0
        pred["SPAMS"][te] = ch[te]
        pred["SPAMS + irreversible"][te] = (ch + dI)[te]
        pred["bucket × g"][te] = ((b.shift(-LAG) - b) * g)[te]
        pred["logger water table × g"][te] = (d.dw_cm * g)[te]
        lev["SPAMS"][te] = R[te]
        lev["bucket × g"][te] = (b * g)[te]
        lev["logger water table × g"][te] = (d.wtd_cm * g)[te]
        folds.append({"held_out_year": y, "tau_days": par["tau"], "x_p": par["x_p"], "x_e": par["x_e"], "x_i": x_i,
                      "g_mm_per_cm": g, "n_train_days": par["n"]})
    dh = d.dh_mm
    rows = []
    for name, p in pred.items():
        for y in YEARS + ("pooled",):
            m = valid & p.notna().to_numpy() & (np.ones(len(d), bool) if y == "pooled" else (d.year == y).to_numpy())
            if m.sum() < 10:
                continue
            e, t = p[m].to_numpy(), dh[m].to_numpy()
            rows.append({"candidate": name, "year": y, "n_days": int(m.sum()), "rmse_mm": float(np.sqrt(np.mean((e - t) ** 2))),
                         "pearson": float(np.corrcoef(e, t)[0, 1]) if np.std(e) > 0 else np.nan,
                         "amplitude_ratio": float(np.std(e) / np.std(t))})
    inc = pd.DataFrame(rows)
    # the like-for-like comparison with the logger: the days of 2021–2024 every candidate covers
    common = valid & d.year.isin(LOGGER_YEARS).to_numpy() & np.all([p.notna().to_numpy() for p in pred.values()], 0)
    comm = pd.DataFrame([{"candidate": k, "n_days": int(common.sum()),
                          "rmse_mm": float(np.sqrt(np.mean((p[common].to_numpy() - dh[common].to_numpy()) ** 2)))}
                         for k, p in pred.items()])
    levels = level_views(d, lev, prot)
    verdict = decide(inc, comm)
    for name, df in (("increments", inc), ("common_2021_2024", comm), ("levels_within_season", levels),
                     ("folds", pd.DataFrame(folds)), ("verdict", verdict)):
        df.to_csv(OUT / f"{name}.csv", index=False)
    (OUT / "README.md").write_text(readme(inc, comm, levels, pd.DataFrame(folds), verdict), encoding="utf-8")
    print(verdict.to_string(index=False))
    print(f"wrote {OUT}")


def level_views(d: pd.DataFrame, lev: dict, prot: Protocol) -> pd.DataFrame:
    """Within-season (May–November) segments per year, split at the laser datum breaks: detrended r with the laser."""
    lo, hi = prot.season_months
    seg = d.year.astype(str)
    for a, b in prot.datum_breaks:
        seg = seg.where(~(d.index > pd.Timestamp(b).tz_convert(None)), seg + "B")
    rows = []
    for name, s in lev.items():
        for key, g in d[(d.index.month >= lo) & (d.index.month <= hi) & d.year.isin(YEARS)].groupby(seg):
            m = g.laser_cm.notna() & s.reindex(g.index).notna()
            if m.sum() < 30:
                continue
            st = detrended_stats(g.index[m], g.laser_cm[m].to_numpy() * 10, s.reindex(g.index)[m].to_numpy())
            rows.append({"candidate": name, "segment": key, "days": int(m.sum()),
                         "pearson_detrended": st.get("pearson_detrended", np.nan)})
    return pd.DataFrame(rows)


def decide(inc: pd.DataFrame, comm: pd.DataFrame) -> pd.DataFrame:
    pv = inc[inc.year == "pooled"].set_index("candidate").rmse_mm
    yr = inc[inc.year != "pooled"].pivot(index="year", columns="candidate", values="rmse_mm")
    years_beat_zero = int((yr["SPAMS"] < yr["zero change"]).sum())
    c = comm.set_index("candidate").rmse_mm
    kept = bool(pv["SPAMS"] < pv["bucket × g"] and years_beat_zero >= 3)
    matches = bool(c["SPAMS"] <= 1.10 * c["logger water table × g"])
    return pd.DataFrame([{"pooled_rmse_spams": pv["SPAMS"], "pooled_rmse_bucket": pv["bucket × g"],
                          "years_spams_beats_zero": years_beat_zero, "kept_as_baseline": kept,
                          "rmse_spams_2021_2024": c["SPAMS"], "rmse_logger_2021_2024": c["logger water table × g"],
                          "matches_logger_prior": matches}])


def _md(df: pd.DataFrame, nd: int = 2) -> str:
    d = df.copy()
    for col in d.columns:
        if pd.api.types.is_float_dtype(d[col]):
            d[col] = d[col].map(lambda x: "—" if pd.isna(x) else f"{x:.{nd}f}")
    return "\n".join(["| " + " | ".join(map(str, d.columns)) + " |", "|" + "---|" * len(d.columns)] +
                     ["| " + " | ".join(map(str, r)) + " |" for r in d.itertuples(index=False)])


def readme(inc, comm, levels, folds, verdict) -> str:
    lv = levels.groupby("candidate").pearson_detrended.median().reset_index(name="median_detrended_r") if len(levels) else levels
    return f"""# X-112 — SPAMS as a baseline for the mat's prior (generated)

Generated by `experiments/spams_baseline_x112.py` (`make run E=spams_baseline_x112`); do not edit. Pre-registration:
`../preregistration/X-112.yaml`. P6, 6-day laser change, year-blocked over the design years (each year predicted from
fits on the others), frozen-soil windows excluded, the same daily table as X-074. SPAMS (conroy2023_spams) on ERA5
precipitation and FAO-56 reference evapotranspiration; the bucket and buoyancy as X-074. Verdict rule read as: kept if
SPAMS' pooled RMSE is below the bucket's AND SPAMS beats zero change in at least 3 of the 5 years; "matches the logger
prior" if its RMSE on the common 2021–2024 days is within 10% of buoyancy on the logger.

## Verdict

{_md(verdict)}

## 6-day change, per year and pooled

{_md(inc)}

## Same days, 2021–2024 (every candidate available)

{_md(comm)}

## Within-season levels (detrended r with the laser, median over season segments)

{_md(lv) if len(lv) else "no segment long enough"}

## Fitted per fold

{_md(folds, 4)}
"""


if __name__ == "__main__":
    main()

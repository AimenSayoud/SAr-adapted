"""X-052: the inversion outputs already computed, checked against the laser at P6/CR (no re-run).

The supervisor: "we could also compare the methods on real data". The per-pixel products (SBAS,
EVD phase linking, ISBAS, hybrid network) and the zone aggregation (phaseG A−C) are compared with
the peat-surface position the laser measures at P6, at the ascending overpass times (descending
products were not inverted per pixel, D-016). Three questions, per method and window:

1. Does the series follow the surface? r with trends removed, and on anomalies (annual cycle removed
   from both; circular-shift p), with the slope product / expected.
2. Is the seasonal amplitude right? Amplitude of the product's annual cycle (all its dates) against
   that of the laser over 2022–2024 (snow-free daily means), both in the product's geometry.
3. Does the change between two dates follow the surface change? Pairs ≤ 24 d, as X-048.

Expected signal from the laser: vertical change dh (cm → mm); in LOS dh · cos(incidence) for the
LOS products (MintPy convention: toward the satellite = up = positive). The laser's ~10° mounting
is not corrected (≈ 1.5 %, pending the field team's answer). Outputs only to the private hub.

    PYTHONPATH=src python scripts/field_methods_vs_laser.py
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from insar_wetlands import field  # noqa: E402
from insar_wetlands import field_link as fl  # noqa: E402
from insar_wetlands.aggregate import seasonal_amplitude  # noqa: E402
from insar_wetlands.hydro_link import _detrend  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
INCIDENCE_ASC = 32.26            # config.yaml, measured (X-042)
OVERPASS_ASC = "16:36"
METHODS = {                      # variable → (label, geometry)
    "sbas_los_mm": ("SBAS (MintPy)", "los"),
    "evd_los_mm": ("EVD phase linking", "los"),
    "isbas_los_mm": ("ISBAS", "los"),
    "hybrid_vertical_mm": ("hybrid network", "vertical"),
}


def expected_mm(surface_cm, geometry: str):
    """Laser surface (cm, up +) in the product's geometry (mm)."""
    v = np.asarray(surface_cm, float) * 10
    return fl.los_from_vertical(v, INCIDENCE_ASC) if geometry == "los" else v


def compare(dates: pd.Series, product: np.ndarray, expected: np.ndarray, max_dt: int = 24) -> dict:
    ok = np.isfinite(product) & np.isfinite(expected)
    d, p, e = pd.to_datetime(dates[ok]).reset_index(drop=True), product[ok], expected[ok]
    out = {"n_dates": int(ok.sum())}
    if len(d) < 12:
        return out
    ac = fl.anomaly_correlation(d, e, p)
    t = fl.years_since(d)
    out |= {"r_trend_removed": ac["r_raw"], "r_anom": ac["r_anom"], "p_anom": ac["p_anom"],
            "slope_product_per_expected": float(np.polyfit(_detrend(e, t), _detrend(p, t), 1)[0])}
    i, j = np.triu_indices(len(d), 1)
    dt = (d.to_numpy()[j] - d.to_numpy()[i]) / np.timedelta64(1, "D")
    keep = dt <= max_dt
    dp, de = p[j] - p[i], e[j] - e[i]
    if keep.sum() >= 12:
        r, pv = fl.circular_shift_p(de[keep], dp[keep])
        out |= {"n_pairs": int(keep.sum()), "r_pairs": r, "p_pairs": pv,
                "slope_pairs": float(np.polyfit(de[keep], dp[keep], 1)[0])}
    return out


def main(argv=None):
    hub = field.field_root().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument("--first", default=str(hub / "08_deliverables" / "field_first"))
    ap.add_argument("--out", default=str(hub / "08_deliverables" / "field_methods_laser"))
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    wtd = field.load_wtd_hourly()
    laser = field.load_laser()
    surface = laser.surface_cm.where(~field.snow_mask(wtd["Air_2m"], laser.index))
    ext = pd.read_csv(Path(a.first) / "s1_plot_by_date.csv")
    ext = ext[(ext["plot"] == "P6") & (ext.track == "ascending") & ext.variable.isin(METHODS)]
    agg = pd.read_csv(REPO / "results" / "tables" / "phaseG_aggregate_series.csv")

    # the laser's own seasonal amplitude over the common period, in each geometry
    daily = surface.resample("D").mean()["2022":"2024"].dropna()
    laser_amp = {g: seasonal_amplitude(pd.DataFrame({"date": daily.index.tz_localize(None),
                                                     "disp_mm": expected_mm(daily.to_numpy(), g)}))["amplitude_mm"]
                 for g in ("los", "vertical")}

    rows, series = [], {}
    cases = [(var, win, g) for var, (_, g) in METHODS.items() for win in ("1x1", "3x3", "5x5")]
    cases.append(("zone_A-C", "zone", "los"))
    for var, win, geom in cases:
        if var == "zone_A-C":
            s = agg.rename(columns={"disp_mm": "value"})[["date", "value"]]
            label = "zone aggregation (mat − grassland, phaseG)"
        else:
            s = ext[(ext.variable == var) & (ext.window == win)].rename(columns={"median": "value"})[["date", "value"]]
            label = METHODS[var][0]
        s = s.sort_values("date").reset_index(drop=True)
        times = pd.to_datetime(s.date + f" {OVERPASS_ASC}", utc=True)
        lz = np.array([fl.value_at(surface, t) for t in times])
        exp = expected_mm(lz, geom)
        prod = s.value.to_numpy(float)
        row = {"method": label, "variable": var, "window": win, "geometry": geom,
               "n_product_dates": int(np.isfinite(prod).sum())}
        if np.isfinite(prod).sum() >= 12:
            amp = seasonal_amplitude(pd.DataFrame({"date": pd.to_datetime(s.date), "disp_mm": prod}).dropna())
            row |= {"product_amplitude_mm": amp["amplitude_mm"], "laser_amplitude_mm": laser_amp[geom],
                    "amplitude_ratio": amp["amplitude_mm"] / laser_amp[geom]}
            row |= compare(s.date, prod, exp)
        rows.append(row)
        if win in ("3x3", "zone"):
            series[label] = (pd.to_datetime(s.date), prod, exp)
    res = pd.DataFrame(rows)
    res.to_csv(out / "methods_vs_laser.csv", index=False)

    # figure: series (demeaned) and pair changes, 3×3 + zone
    shown = [k for k, (_, p, _) in series.items() if np.isfinite(p).sum() >= 12]
    fig, ax = plt.subplots(2, len(shown), figsize=(4.2 * len(shown), 7.2), squeeze=False)
    full = daily.reindex(pd.date_range(daily.index[0], daily.index[-1], freq="D"))   # gaps (snow) stay gaps
    ld = full.index.tz_localize(None)
    for c, k in enumerate(shown):
        d, p, e = series[k]
        geom = "vertical" if "hybrid" in k else "los"
        base = expected_mm(full.to_numpy(), geom)
        ax[0, c].xaxis.set_major_locator(mdates.YearLocator())
        ax[0, c].xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
        ax[0, c].plot(ld, base - np.nanmean(base), color="0.6", lw=0.8, label="laser (snow-free, daily)")
        ax[0, c].plot(d, p - np.nanmean(p), "o-", ms=2.5, lw=0.8, color="C3", label="product")
        ax[0, c].set_title(k, fontsize=9)
        ax[0, c].set_ylabel(f"{'vertical' if geom == 'vertical' else 'LOS'} (mm, mean removed)")
        ax[0, c].legend(fontsize=7)
        r = res[(res.method == k) & res.window.isin(["3x3", "zone"])].iloc[0]
        ok = np.isfinite(p) & np.isfinite(e)
        dd, pp, ee = d[ok].to_numpy(), p[ok], e[ok]
        i, j = np.triu_indices(len(dd), 1)
        keep = (dd[j] - dd[i]) / np.timedelta64(1, "D") <= 24
        ax[1, c].scatter(ee[j][keep] - ee[i][keep], pp[j][keep] - pp[i][keep], s=10, alpha=0.6)
        lim = np.nanmax(np.abs(np.r_[ee[j][keep] - ee[i][keep], pp[j][keep] - pp[i][keep], 1]))
        ax[1, c].plot([-lim, lim], [-lim, lim], "k--", lw=0.8, label="1:1 (true motion)")
        ax[1, c].axhline(0, color="0.8", lw=0.5); ax[1, c].axvline(0, color="0.8", lw=0.5)
        ax[1, c].set_xlabel("laser change between the two dates (mm)")
        ax[1, c].set_ylabel("product change (mm)")
        ax[1, c].set_title(f"pairs ≤ 24 d: r = {r.get('r_pairs', np.nan):.2f}, slope {r.get('slope_pairs', np.nan):.2f}", fontsize=9)
        ax[1, c].legend(fontsize=7)
    fig.suptitle("P6/CR, ascending: the inversion outputs against the laser-measured surface", fontsize=10)
    fig.tight_layout()
    fig.savefig(out / "fig_methods_vs_laser.png", dpi=150)
    plt.close(fig)
    pd.set_option("display.width", 220)
    print(res.round(3).to_string(index=False))


if __name__ == "__main__":
    main()

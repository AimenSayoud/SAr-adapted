"""X-054: the aggregated zone series over five years (2020–2024), both tracks.

The zone phase of every interferogram (2022–2024 stacks + the 2020–2021 extension, D-020) was
computed by the X-050 run (`zone_pair_table_2020_2024.csv`: coherence-weighted mean unwrapped phase,
zone minus reference, per pair). Here each network is inverted as one super-pixel
(`aggregate.invert_aggregate`, the phaseG code) and fitted with the seasonal model
(`seasonal_amplitude`), over the five years and year by year. Controls: lake − grassland, mat −
lake, the adjacent null. Three networks: full, ≤ 60 d, and one matched between the periods (12-day
grid, pairs ≤ 24 d). The two stacks share no pair (no bridge pairs, SETTLED.md): the inversion joins
2021 to 2022 with a zero increment, so the five-year fit assumes continuity across the gap and the
per-year fits are the clean comparison. A reproduction check first: the 2022–2024 part alone must match the
series the project already uses. Then the series against the measured water table (plot median at
the overpass, dry wells out): anomaly correlation with circular-shift p.

Outputs only to the private hub (the water-table comparison uses field data).

    PYTHONPATH=src python scripts/field_series_x054.py
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from insar_wetlands import field  # noqa: E402
from insar_wetlands import field_link as fl  # noqa: E402
from insar_wetlands.aggregate import (  # noqa: E402
    PHASE_TO_MM,
    filter_pairs,
    invert_aggregate,
    seasonal_amplitude,
)

OVERPASS = {"ascending": "16:36", "descending": "05:09"}
EPOCH = pd.Timestamp("2020-01-01")   # one phase origin for every fit (comparable phases)


def amp_ci(fit: dict, n: int = 20000, seed: int = 0) -> list[float]:
    ab = np.random.default_rng(seed).multivariate_normal([fit["a_cos_mm"], fit["b_sin_mm"]], fit["cov_a_b"], n)
    a = np.hypot(ab[:, 0], ab[:, 1])
    return [float(np.percentile(a, 2.5)), float(np.percentile(a, 97.5))]


def invert(d: pd.DataFrame) -> pd.DataFrame:
    dd = d.assign(ddphase_rad=d.ddisp_mm / PHASE_TO_MM)[["pair", "dt_days", "ddphase_rad", "ddisp_mm", "weight"]]
    return invert_aggregate(dd.reset_index(drop=True))


def matched_network(d: pd.DataFrame, max_dt: int = 24) -> pd.DataFrame:
    """The same geometry in both periods: dates on the 2022–2024 12-day grid only (drops the 2020–2021
    Sentinel-1B dates, 6 days off it) and pairs ≤ max_dt days. The 2020–2021 extension otherwise has
    6-day pairs and a short-pairs-only network, which the 2022–2024 stack does not."""
    anchor = pd.Timestamp(d.loc[d["stack"] == "2022–2024", "date1"].min())
    on_grid = lambda c: (pd.to_datetime(d[c]) - anchor).dt.days % 12 == 0  # noqa: E731
    return d[on_grid("date1") & on_grid("date2") & (d.dt_days <= max_dt)]


def main(argv=None):
    hub = field.field_root().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", default=str(hub / "08_deliverables" / "field_dew_x050" / "zone_pair_table_2020_2024.csv"))
    ap.add_argument("--committed", default=str(hub / "05_code" / "local" / "results" / "asc_desc" / "asc_desc_web.json"))
    ap.add_argument("--out", default=str(hub / "08_deliverables" / "field_series_x054"))
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    zt = pd.read_csv(a.pairs)
    committed = json.loads(Path(a.committed).read_text())["series"]

    wtd = field.load_wtd_hourly()
    cens = {p: field.censored_flag(wtd[p]) for p in field.PLOTS}

    series, fits, checks, corr = [], [], [], []
    for (track, label), d in zt.groupby(["track", "zones"]):
        d = d.drop_duplicates("pair")
        # reproduction: the 2022–2024 stack alone, inverted the same way
        s22 = invert(d[d["stack"] == "2022–2024"])
        key = f"{track}|{'null' if label == 'adjacent null' else label}|full"
        if key in committed:
            ref = pd.DataFrame(committed[key]).dropna()
            m = s22.assign(date=s22.date.dt.strftime("%Y-%m-%d")).merge(ref, on="date", suffixes=("", "_ref"))
            checks.append({"track": track, "zones": label, "n_dates": len(m),
                           "max_abs_diff_mm": float((m.disp_mm - m.disp_mm_ref).abs().max())})
        for variant, dv in (("full", d), ("<=60d", filter_pairs(d, max_dt_days=60)),
                            ("matched 12-d grid, <=24d", matched_network(d))):
            s = invert(dv)
            series.append(s.assign(track=track, zones=label, variant=variant))
            for period, sp in [("2020–2024", s)] + [(str(y), s[s.date.dt.year == y]) for y in range(2020, 2025)]:
                if len(sp) < 8:
                    continue
                f = seasonal_amplitude(sp, epoch=EPOCH)
                fits.append({"track": track, "zones": label, "variant": variant, "period": period, "n_dates": len(sp),
                             "amplitude_mm": f["amplitude_mm"], "ci95_low": amp_ci(f)[0], "ci95_high": amp_ci(f)[1],
                             "phase_doy": f["phase_doy"], "r2_seasonal": f["r2_seasonal"], "trend_mm_yr": f["trend_mm_yr"]})
        # the full-network series against the measured water table (level), anomalies
        s = invert(d)
        lvl = []
        for t in s.date:
            tt = pd.Timestamp(f"{t.date()} {OVERPASS[track]}", tz="UTC")
            v = [field._interp_at(wtd[p], tt) for p in field.PLOTS if not bool(cens[p].asof(tt))]
            lvl.append(float(np.median(v)) if len(v) >= 5 else np.nan)
        corr.append({"track": track, "zones": label, **fl.anomaly_correlation(s.date, np.array(lvl), s.disp_mm.to_numpy())})

    ser = pd.concat(series, ignore_index=True)
    ser.to_csv(out / "series_2020_2024.csv", index=False)
    fits = pd.DataFrame(fits)
    fits.to_csv(out / "seasonal_fits.csv", index=False)
    pd.DataFrame(checks).to_csv(out / "reproduction_check.csv", index=False)
    pd.DataFrame(corr).to_csv(out / "series_vs_wtd.csv", index=False)

    # figure: the four series per track (full network), with the regional WTD beneath
    fig, ax = plt.subplots(3, 1, figsize=(13, 9), sharex=True)
    for i, track in enumerate(OVERPASS):
        for lab, c in (("A-C", "C3"), ("B-C", "C0"), ("A-B", "C4"), ("adjacent null", "0.5")):
            s = ser[(ser.track == track) & (ser.zones == lab) & (ser.variant == "full")]
            ax[i].plot(s.date, s.disp_mm, lw=1.2 if lab == "A-C" else 0.8, color=c, label=lab)
        ax[i].axvline(pd.Timestamp("2022-01-01"), color="k", ls=":", lw=0.8)
        ax[i].set_ylabel(f"{track}\\nLOS (mm)"); ax[i].legend(fontsize=7, ncol=4)
    daily = wtd[field.PLOTS].where(~pd.DataFrame(cens)).median(axis=1).resample("D").mean()
    ax[2].plot(daily.index.tz_localize(None), daily.values, color="C2", lw=0.9)
    ax[2].set_ylabel("WTD, plot median (cm)")
    fig.suptitle("X-054 — aggregated zone series 2020–2024 (dotted: start of the 2022–2024 stacks)", fontsize=10)
    fig.tight_layout(); fig.savefig(out / "fig_series_2020_2024.png", dpi=140); plt.close(fig)
    pd.set_option("display.width", 220)
    print(pd.DataFrame(checks).round(4).to_string(index=False))
    print(fits[(fits.variant == "full") & (fits.zones == "A-C")].round(3).to_string(index=False))
    print(pd.DataFrame(corr).round(3).to_string(index=False))


if __name__ == "__main__":
    main()

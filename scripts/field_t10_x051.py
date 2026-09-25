"""X-051: T10 (coherence vs water-table change, freeze gain) with the MEASURED water table.

T10 regresses each zone's pair coherence on |Δ water table|, where the water table was a proxy
built from ERA5 precipitation (``stratify.pair_hydro_change`` → ``validation.hydrology_proxy``).
The field delivery measures it. This script

1. re-runs T10 exactly as ``export_figures_en`` does (control: must reproduce the committed
   ``results/tables/T10_hydrology_freeze.csv``);
2. swaps the proxy for the measured |ΔWTD| (median over the 9 plots at the overpass time,
   dry-well stretches out, ≥ 5 plots) — same code (``coherence_vs_hydro``), same pairs where both
   exist — and repeats it on the descending track;
3. adds what T10 lacks: the season-cleaned version (coherence with its annual cycle and
   baseline trend regressed out, ``field_link.season_residual``; p from circular shifts), since
   both |ΔWTD| and coherence are seasonal;
4. checks the proxy itself against the measurement, and the freeze gain with the station's own
   air temperature instead of ERA5.

The committed T10 is not changed. Outputs go ONLY to the private hub (``--out``, default
``<hub>/08_deliverables/field_t10_x051``).

    INSAR_DRIVE_ROOT=<snapshot> PYTHONPATH=src python scripts/field_t10_x051.py
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
from insar_wetlands.bootstrap import start  # noqa: E402
from insar_wetlands.hydro import open_era5  # noqa: E402
from insar_wetlands.stratify import (  # noqa: E402
    coherence_vs_hydro,
    freeze_coherence_gain,
    pair_hydro_change,
)

REPO = Path(__file__).resolve().parents[1]
TABLES = REPO / "results" / "tables"
OVERPASS = {"ascending": "16:36", "descending": "05:09"}
COH_FILES = {"ascending": "phaseD_coh_by_zone.csv", "descending": "phaseD_coh_by_zone_descending.csv"}


def measured_pair_hydro(pairs: list[str], track: str, wtd: pd.DataFrame, cens: dict) -> pd.DataFrame:
    """Per pair: |median ΔWTD over the plots| (cm) at the overpass and the lower of the station's
    daily-mean air temperatures of the two dates (K — the same definition as ``pair_hydro_change``
    with ERA5)."""
    tmin_daily = wtd["Air_2m"].resample("D").mean()
    rows = []
    for p in pairs:
        a, b = p.split("_")
        t1, t2 = (pd.Timestamp(f"{x} {OVERPASS[track]}", tz="UTC") for x in (a, b))
        dw = [field._interp_at(wtd[q], t2) - field._interp_at(wtd[q], t1) for q in field.PLOTS
              if not (bool(cens[q].asof(t1)) or bool(cens[q].asof(t2)))]
        tm = [tmin_daily.get(pd.Timestamp(x, tz="UTC"), np.nan) for x in (a, b)]
        rows.append({"pair": p, "dwtd": abs(float(np.median(dw))) if len(dw) >= 5 else np.nan, "n_plots": len(dw),
                     "tmin_station": min(tm) + 273.15 if np.all(np.isfinite(tm)) else np.nan})
    return pd.DataFrame(rows).set_index("pair")


def season_cleaned(dfz: pd.DataFrame, hyd: pd.DataFrame, max_dt: int | None) -> pd.DataFrame:
    out = []
    for z, g in dfz.groupby("zone"):
        m = g.merge(hyd[["dwtd"]], left_on="pair", right_index=True)
        m = m[np.isfinite(m.dwtd) & np.isfinite(m.mean_coh)]
        if max_dt:
            m = m[m.dt_days <= max_dt]
        m = m.assign(mid=[pd.Timestamp(p[:8]) + (pd.Timestamp(p[9:]) - pd.Timestamp(p[:8])) / 2 for p in m.pair]).sort_values("mid")
        res = fl.season_residual(m.mean_coh, m.mid, m.dt_days)
        r, p = fl.circular_shift_p(m.dwtd.to_numpy(), res)
        out.append({"zone": z, "r_season_cleaned": r, "p_season_cleaned": p, "n": len(m)})
    return pd.DataFrame(out)


def main(argv=None):
    hub = field.field_root().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(hub / "08_deliverables" / "field_t10_x051"))
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    ctx = start("field_t10_x051", mount=False, git=False)
    D, cfg = ctx.paths.drive, ctx.cfg
    lon, lat = cfg["site"]["centroid"]
    era5 = open_era5(D)
    wtd = field.load_wtd_hourly()
    cens = {p: field.censored_flag(wtd[p]) for p in field.PLOTS}

    report, rows, pair_rows = {}, [], []
    for track, f in COH_FILES.items():
        dfz = pd.read_csv(TABLES / f)
        pairs = sorted(dfz.pair.unique())
        proxy = pair_hydro_change(pairs, era5, lon, lat)
        if track == "ascending":   # 1. control: the committed table, same code path
            t10 = coherence_vs_hydro(dfz, proxy).merge(freeze_coherence_gain(dfz, proxy), on="zone", how="outer")
            committed = pd.read_csv(TABLES / "T10_hydrology_freeze.csv")
            report["control_reproduces_committed_T10"] = bool(np.allclose(
                t10[committed.columns[1:]].to_numpy(float), committed[committed.columns[1:]].to_numpy(float), atol=1e-9))
        meas = measured_pair_hydro(pairs, track, wtd, cens)
        both = meas.index[np.isfinite(meas.dwtd) & np.isfinite(proxy.dwtd.reindex(meas.index))]
        dts = pd.Series({p: (pd.Timestamp(p[9:]) - pd.Timestamp(p[:8])).days for p in pairs})
        for variant, hyd in (("proxy (ERA5 precipitation), all pairs", proxy),
                             ("proxy, same pairs as measured", proxy.loc[both]),
                             ("measured WTD", meas.loc[both])):
            for max_dt in (None, 48):
                sub = dfz if max_dt is None else dfz[dfz.dt_days <= max_dt]
                h = hyd if max_dt is None else hyd.loc[[p for p in hyd.index if dts[p] <= max_dt]]
                t = coherence_vs_hydro(sub, h).merge(season_cleaned(sub, h, None), on=["zone", "n"], how="left")
                t.insert(0, "pairs", "all" if max_dt is None else f"≤ {max_dt} d")
                t.insert(0, "water_table", variant)
                t.insert(0, "track", track)
                rows.append(t)
        # 4. the proxy against the measurement, and freeze gain with the station temperature
        j = proxy.join(meas, rsuffix="_measured").loc[both]
        report.setdefault("proxy_vs_measured_abs_change", {})[track] = {
            "r_all_pairs": float(np.corrcoef(j.dwtd, j.dwtd_measured)[0, 1]),
            "r_pairs_le48d": float(np.corrcoef(j.dwtd[dts[j.index] <= 48], j.dwtd_measured[dts[j.index] <= 48])[0, 1]),
            "n": int(len(j))}
        fg_era5 = freeze_coherence_gain(dfz, proxy).assign(temperature="ERA5 t2m (daily mean)")
        st = meas[["tmin_station"]].rename(columns={"tmin_station": "tmin"}).dropna()
        fg_st = freeze_coherence_gain(dfz, st).assign(temperature="station air (daily mean)")
        report.setdefault("freeze_gain", {})[track] = pd.concat([fg_era5, fg_st]).to_dict("records")
        pair_rows.append(j.assign(track=track, dt_days=dts[j.index].values).reset_index())

    tab = pd.concat(rows, ignore_index=True)
    tab.to_csv(out / "t10_proxy_vs_measured.csv", index=False)
    pairs_tab = pd.concat(pair_rows, ignore_index=True).rename(columns={"dwtd": "abs_dproxy", "dwtd_measured": "abs_dwtd_cm"})
    pairs_tab.to_csv(out / "pair_water_table_change.csv", index=False)
    (out / "report.json").write_text(json.dumps(report, indent=1, default=float))

    fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
    asc = pairs_tab[pairs_tab.track == "ascending"]
    ax[0].scatter(asc.abs_dwtd_cm, asc.abs_dproxy, s=6, alpha=0.5, c=np.where(asc.dt_days <= 48, "#2563eb", "#a1a1aa"))
    ax[0].set_xlabel("|ΔWTD| measured, median over plots (cm)"); ax[0].set_ylabel("|Δ proxy| (ERA5 precipitation)")
    ax[0].set_title(f"The T10 proxy against the measurement (ascending pairs; r = "
                    f"{report['proxy_vs_measured_abs_change']['ascending']['r_all_pairs']:.2f})", fontsize=9)
    d = tab[(tab.pairs == "≤ 48 d") & (tab.zone.isin(["A", "B", "C"]))]
    labels = ["proxy, same pairs as measured", "measured WTD"]
    x = np.arange(3)
    for k, (track, lab) in enumerate([(t, v) for t in OVERPASS for v in labels]):
        s = d[(d.track == track) & (d.water_table == lab)].set_index("zone").reindex(["A", "B", "C"])
        ax[1].bar(x + (k - 1.5) * 0.2, s.r_season_cleaned, 0.2, label=f"{track}, {lab.split(',')[0]}",
                  color=["#fbbf24", "#d97706", "#93c5fd", "#2563eb"][k])
        for xi, p in enumerate(s.p_season_cleaned):
            ax[1].text(xi + (k - 1.5) * 0.2, min(0, s.r_season_cleaned.iloc[xi]) - 0.02, "*" if p < 0.05 else "", ha="center")
    ax[1].set_ylim(min(d.r_season_cleaned.min(), 0) - 0.05, 0.02)
    ax[1].axhline(0, c="k", lw=0.6); ax[1].set_xticks(x, ["A mat", "B lake", "C grassland"])
    ax[1].set_ylabel("r (coherence vs |ΔWTD|), season removed"); ax[1].legend(fontsize=7)
    ax[1].set_title("T10 season-cleaned, pairs ≤ 48 d (* p < 0.05, circular shifts)", fontsize=9)
    fig.tight_layout()
    fig.savefig(out / "fig_t10_measured.png", dpi=150)
    plt.close(fig)
    print(json.dumps(report, indent=1, default=float)[:3000])
    print(f"wrote {sorted(p.name for p in out.iterdir())}")


if __name__ == "__main__":
    main()

"""X-049: do vegetation and surface conditions (UAV, LAI) explain why the plots behave differently
in C-band radar?

The supervisor's step 7–8: keep the UAV variables few and physical — NDVI, NDRE, red-edge / NIR,
Altum surface temperature, SunScan LAI — and use them to explain *spatial* differences in
Sentinel-1 coherence and response among the plots, not to correlate every index.

Two analyses, both built so that the seasonal cycle cannot create a relationship:

1. Between plots. Each UAV variable is expressed per campaign relative to that day's mean over the
   plots (removes date, weather and flight-time effects), then summarised per plot over the growing
   season (May–September). Compared with each plot's radar behaviour — mean short-baseline
   coherence, temporal coherence, mean VV / VH, and its sensitivity to the water table (X-047) —
   and its hydrology (mean WTD, distance to the mat edge). Spearman with exact permutation p
   (all 9! orderings). Radar at 3×3 and at the plot's own pixel (1×1), because neighbouring 3×3
   windows overlap.
2. Campaign dates. For each campaign and plot, the nearest ascending acquisition (± 6 d); UAV and
   radar values with the plot mean AND the campaign mean removed (two-way), so what is left is
   "this plot differs from its usual self on this day more than the others do". Permutation p:
   plot labels shuffled within each campaign.

UAV cells the field team marked yellow are excluded (field.load_uav_table). Outputs only to the
private hub.

    PYTHONPATH=src python scripts/field_uav_x049.py
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from insar_wetlands import field  # noqa: E402
from insar_wetlands import field_link as fl  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
GROWING = (5, 9)   # May–September: full canopy, no snow
UAV_VARS = {"ndvi": "NDVI", "ndre": "NDRE", "nir_842": "NIR 842 reflectance", "re_717": "red-edge 717 reflectance",
            "thermal_c": "surface temperature (°C)", "lai": "LAI"}


def uav_indices(u: pd.DataFrame) -> pd.DataFrame:
    nir, red, re = u["REMX_NIR_842_mean"], u["REMX_Red_668_mean"], u["REMX_Red_Edge_717_mean"]
    return pd.DataFrame({"date": u.Date, "plot": u.base_plot, "subplot": u.Plot,
                         "ndvi": (nir - red) / (nir + red), "ndre": (nir - re) / (nir + re), "nir_842": nir, "re_717": re,
                         "thermal_c": u["Altum_Thermal_11um_mean"], "lai": u["LAI"]})


def relative_to_campaign(v: pd.DataFrame) -> pd.DataFrame:
    """Per campaign: plot median over sub-plots, minus that campaign's mean over the plots."""
    per = v.groupby(["date", "plot"])[list(UAV_VARS)].median()
    return per - per.groupby(level="date").transform("mean")


def plot_radar_metrics(first: Path, x047: Path) -> pd.DataFrame:
    j = pd.read_csv(first / "plot_s1_wtd.csv", parse_dates=["date"])
    ext = pd.read_csv(first / "s1_plot_by_date.csv", parse_dates=["date"])
    tc = pd.read_csv(first / "s1_plot_temporal_coherence.csv")
    ppc = pd.read_csv(x047 / "per_pair_correlations.csv")
    pdc = pd.read_csv(x047 / "per_date_correlations.csv")
    gs = lambda d: d[d.date.dt.month.between(*GROWING)]  # noqa: E731
    out = {}
    for win in ("3x3", "1x1"):
        e = gs(ext[(ext.window == win) & (ext.track == "ascending")])
        m = e.pivot_table(index="plot", columns="variable", values="median", aggfunc="mean")
        out[f"coh_{win}"] = m["coh_pairs_le24d"]
        out[f"vv_{win}"] = m["vv_db"]
        out[f"vh_{win}"] = m["vh_db"]
        out[f"tcoh_{win}"] = tc[tc.window == win].set_index("plot")["median"]
    d = gs(j[j.track == "descending"]).groupby("plot").coh_pairs_le24d.mean()
    out["coh_desc_3x3"] = d
    out["coh_sens_to_dwtd"] = ppc[ppc.track == "ascending"].set_index("plot").r_coh_vs_abs_dwtd
    out["vv_sens_to_wtd"] = pdc[(pdc.track == "ascending") & (pdc.s1 == "vv_db") & (pdc.wtd == "wtd_at")].set_index("plot").r_anom
    return pd.DataFrame(out)


def hydro_context() -> pd.DataFrame:
    from pyproj import Transformer
    from shapely.geometry import Point
    from shapely.ops import transform

    from insar_wetlands.aoi import load_aoi
    from insar_wetlands.config import load_config
    wtd = field.load_wtd_hourly()["2022":"2024"]
    cens = pd.DataFrame({p: field.censored_flag(wtd[p]) for p in field.PLOTS})
    w = wtd[field.PLOTS].where(~cens)
    w = w[w.index.month.isin(range(GROWING[0], GROWING[1] + 1))]
    aoi = transform(Transformer.from_crs(4326, 32633, always_xy=True).transform, load_aoi(load_config(REPO / "config" / "config.yaml")))
    plots = field.load_plots().groupby("base_plot")[["x", "y"]].mean()
    return pd.DataFrame({"wtd_mean_gs_cm": w.mean(), "wtd_sd_gs_cm": w.std(),
                         "dist_to_mat_edge_m": [aoi.exterior.distance(Point(r.x, r.y)) for r in plots.reindex(field.PLOTS).itertuples()]},
                        index=field.PLOTS)


def campaign_matched(v: pd.DataFrame, first: Path, max_days: int = 6) -> pd.DataFrame:
    j = pd.read_csv(first / "plot_s1_wtd.csv", parse_dates=["date"])
    s1 = j[j.track == "ascending"][["date", "plot", "coh_pairs_le24d", "vv_db", "vh_db", "wtd_at"]]
    per = v.groupby(["date", "plot"])[list(UAV_VARS)].median().reset_index()
    rows = []
    for (d, p), g in per.groupby(["date", "plot"]):
        c = s1[s1["plot"] == p]
        k = (c.date - d).abs()
        if k.min() <= pd.Timedelta(days=max_days):
            r = c.loc[k.idxmin()]
            rows.append({**g.iloc[0].to_dict(), "s1_date": r.date, "days_apart": (r.date - d).days,
                         "coh": r.coh_pairs_le24d, "vv_db": r.vv_db, "vh_db": r.vh_db, "wtd_at_s1": r.wtd_at})
    return pd.DataFrame(rows)


def two_way_correlation(df: pd.DataFrame, a: str, b: str, n_perm: int = 5000, rng=None) -> dict:
    """Correlation of a and b after removing plot and campaign means; p by shuffling plot labels
    within campaigns."""
    d = df[["date", "plot", a, b]].dropna()
    d = d[d.groupby("date")["plot"].transform("size") >= 3]
    if len(d) < 15:
        return {"n": len(d), "r": np.nan, "p": np.nan}

    date = pd.factorize(d["date"])[0]
    plot = pd.factorize(d["plot"])[0]

    def resid(x):   # remove campaign means, then plot means (numpy, for speed in the permutation loop)
        x = x - (np.bincount(date, x) / np.bincount(date))[date]
        return x - (np.bincount(plot, x) / np.bincount(plot))[plot]

    xa, xb = d[a].to_numpy(float), d[b].to_numpy(float)
    ra = resid(xa)
    r0 = float(np.corrcoef(ra, resid(xb))[0, 1])
    rng = rng or np.random.default_rng(0)
    groups = [np.flatnonzero(date == k) for k in range(date.max() + 1)]
    null = np.empty(n_perm)
    for i in range(n_perm):
        s = xb.copy()
        for g in groups:
            s[g] = s[rng.permutation(g)]
        null[i] = np.corrcoef(ra, resid(s))[0, 1]
    return {"n": len(d), "r": r0, "p": float((np.sum(np.abs(null) >= abs(r0)) + 1) / (n_perm + 1))}


def main(argv=None):
    hub = field.field_root().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument("--first", default=str(hub / "08_deliverables" / "field_first"))
    ap.add_argument("--x047", default=str(hub / "08_deliverables" / "field_x047_x048"))
    ap.add_argument("--out", default=str(hub / "08_deliverables" / "field_uav_x049"))
    a = ap.parse_args(argv)
    out, first = Path(a.out), Path(a.first)
    out.mkdir(parents=True, exist_ok=True)

    u = field.load_uav_table()
    v = uav_indices(u)
    v = v[v.date.dt.year.between(2022, 2024)]
    avail = v.groupby(v.date.dt.strftime("%Y-%m-%d"))[list(UAV_VARS)].count()
    avail.to_csv(out / "uav_availability.csv")

    rel = relative_to_campaign(v).reset_index()
    rel_gs = rel[pd.to_datetime(rel.date).dt.month.between(*GROWING)]
    uav_plot = rel_gs.groupby("plot")[list(UAV_VARS)].median().add_suffix("_rel")
    uav_plot["n_campaigns_gs"] = rel_gs.groupby("plot").ndvi.count()
    radar = plot_radar_metrics(first, Path(a.x047))
    hydro = hydro_context()
    summary = uav_plot.join(radar, how="outer").join(hydro, how="outer").reindex(field.PLOTS)
    summary.index.name = "plot"
    summary.to_csv(out / "plot_summary.csv")

    explanatory = [f"{k}_rel" for k in UAV_VARS] + list(hydro.columns)
    # P8 and P9 share one radar pixel (and window): as radar units they are one. The 8-unit test
    # averages them, so the exact p is not inflated by a duplicated plot.
    units = summary.rename(index={"P8": "P8+P9", "P9": "P8+P9"}).groupby(level=0).mean()
    rows = []
    for ex in explanatory:
        for rm in radar.columns:
            rho, p, n = fl.spearman_exact(summary[ex], summary[rm])
            rho8, p8, n8 = fl.spearman_exact(units[ex], units[rm])
            rows.append({"explanatory": ex, "radar": rm, "n_plots": n, "rho": rho, "p_exact": p,
                         "n_units": n8, "rho_8_units": rho8, "p_exact_8_units": p8})
    sp = pd.DataFrame(rows)
    sp.to_csv(out / "between_plot_spearman.csv", index=False)
    # How much the explanatory variables move together along the transect (the confounding).
    colin = pd.DataFrame([{"a": x, "b": y, "rho_8_units": fl.spearman_exact(units[x], units[y])[0]}
                          for i, x in enumerate(explanatory) for y in explanatory[i + 1:]])
    colin.to_csv(out / "explanatory_intercorrelation.csv", index=False)

    cm = campaign_matched(v, first)
    cm.to_csv(out / "campaign_matched.csv", index=False)
    tw = pd.DataFrame([{"uav": k, "radar": r, **two_way_correlation(cm, k, r)}
                       for k in UAV_VARS for r in ("coh", "vv_db", "vh_db")])
    tw.to_csv(out / "campaign_two_way.csv", index=False)

    # figure: the between-plot relations with the smallest p (3×3 radar), labelled by plot
    top = sp[sp.radar.str.contains("3x3|sens|desc")].nsmallest(6, "p_exact")
    fig, ax = plt.subplots(2, 3, figsize=(13, 7.5))
    for axi, r in zip(ax.ravel(), top.itertuples()):
        axi.scatter(summary[r.explanatory], summary[r.radar], color="C0")
        for p_, x_, y_ in zip(summary.index, summary[r.explanatory], summary[r.radar]):
            axi.annotate(p_, (x_, y_), fontsize=8, xytext=(3, 3), textcoords="offset points")
        axi.set_xlabel(r.explanatory); axi.set_ylabel(r.radar)
        axi.set_title(f"ρ = {r.rho:+.2f}, exact p = {r.p_exact:.3f} (n = {r.n_plots})", fontsize=9)
    fig.suptitle(f"Between plots: the 6 strongest of {len(sp)} tests (≈ {0.05 * len(sp):.0f} expected at p < 0.05 by chance)", fontsize=10)
    fig.tight_layout()
    fig.savefig(out / "fig_uav_between_plots.png", dpi=150)
    plt.close(fig)
    pd.set_option("display.width", 220)
    print("flagged cells masked:", u.attrs.get("n_flagged"))
    print(summary.round(3).to_string())
    print(sp.sort_values("p_exact").head(14).round(3).to_string(index=False))
    print("n tests:", len(sp), "p<0.05 (9 plots):", int((sp.p_exact < 0.05).sum()), "(8 units):", int((sp.p_exact_8_units < 0.05).sum()))
    print(colin[colin.rho_8_units.abs() >= 0.7].round(2).to_string(index=False))
    print(tw.round(3).to_string(index=False))


if __name__ == "__main__":
    main()

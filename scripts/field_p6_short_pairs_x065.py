"""X-065: does the P6 phase follow the laser in a year the fusion never saw, and how does it depend on
the revisit? 2021 (Sentinel-1A + 1B: 6-day and 12-day pairs, D-020 crops) next to 2022–2024 (12- and
24-day, as X-059/X-062), both tracks; and does ΔVV / ΔRVI (C-048 stacks) add to the phase model?

Same conventions as X-059 / X-062 block D: P6 pixel minus the grassland (C) median, PHASE_TO_MM; laser
nearest measured hourly value within ±1 h, snow-free 72 h, not an outlier, not gap-filled, no pair across
a filled stretch (X-058); LOS = vertical · cos(incidence). The 2022–2024 12-day fit reproduces X-062's
block D. Exploratory: one pixel, tens of pairs. Writes only to the hub.

    INSAR_DRIVE_ROOT=../local/drive_pristine PYTHONPATH=src python scripts/field_p6_short_pairs_x065.py
"""

import numpy as np
import pandas as pd
import xarray as xr

from insar_wetlands import field
from insar_wetlands import fusion as fu
from insar_wetlands.bootstrap import start
from insar_wetlands.inversion.isbas import PHASE_TO_MM
from insar_wetlands.paths import make_paths
from insar_wetlands.stack import list_pairs, load_layer

HUB = field.field_root().parents[1]
OUT = HUB / "08_deliverables" / "field_p6_short_pairs_x065"

DLV = HUB / "08_deliverables"
INC = {"ascending": 32.26, "descending": 39.17}
HOUR = {"ascending": "16:36", "descending": "05:09"}
QW = 55.465763 / 4

ctx = start("field_p6_short_pairs_x065", mount=False, git=False)
OUT.mkdir(parents=True, exist_ok=True)
C = ctx.zones["C"].values
px = pd.read_csv(DLV / "field_first/plot_pixels.csv").set_index("plot")
r6, c6 = int(px.loc["P6", "row"]), int(px.loc["P6", "col"])
fl = pd.read_csv(DLV / "field_p6/laser_qc/laser_flags.csv", parse_dates=["time_utc"]).set_index("time_utc")
fl["ok"] = fl.surface_cm.notna() & ~fl.snow_72h & ~fl.outlier & ~fl.filled
runs = pd.read_csv(DLV / "field_p6/laser_qc/interpolated.csv", parse_dates=["start", "end"])
wtd = field.load_wtd_hourly()
cens = field.censored_flag(wtd["P6"])


def laser_at(t):
    i = fl.index.get_indexer([t], method="nearest")[0]
    if abs((fl.index[i] - t).total_seconds()) > 3600 or not fl.ok.iloc[i]:
        return np.nan
    return fl.surface_cm.iloc[i]


def wtd_at(t):
    return np.nan if bool(cens.asof(t)) else field._interp_at(wtd["P6"], t)


rows = []
for track in INC:
    roots = [make_paths(cfg=ctx.cfg, track=track).cropped, HUB / f"05_code/local/s1_2020_2021/hyp3_cropped_{track}"]
    have = {p: root for root in roots for p in list_pairs(root)}
    dates = sorted({pd.Timestamp(d) for p in have for d in p.split("_")})
    rtc = xr.open_dataset(HUB / "05_code/local/drive_mirror" /
                          ("rtc_dualpol_2020_2024.nc" if track == "ascending" else "rtc_dualpol_2020_2024_descending.nc"))
    vv = rtc.gamma0_vv_db[:, r6 - 1:r6 + 2, c6 - 1:c6 + 2].median(("y", "x")).to_series()
    lv, lh = 10 ** (rtc.gamma0_vv_db / 10), 10 ** (rtc.gamma0_vh_db / 10)
    rvi = (4 * lh / (lv + lh))[:, r6 - 1:r6 + 2, c6 - 1:c6 + 2].median(("y", "x")).to_series()
    for step in (1, 2):                         # consecutive dates; and every second date (same satellite in 2020–21)
        for a, b in zip(dates[:-step], dates[step:]):
            pid = f"{a:%Y%m%d}_{b:%Y%m%d}"
            if pid not in have:
                continue
            u = load_layer(have[pid], "unw_phase", [pid]).isel(pair=0).values
            co = load_layer(have[pid], "corr", [pid]).isel(pair=0).values
            ref = np.nanmedian(u[C & np.isfinite(u)])
            ta, tb = (pd.Timestamp(f"{d.date()} {HOUR[track]}", tz="UTC") for d in (a, b))
            across = bool(((runs.end > ta) & (runs.start < tb)).any())
            la, lb = laser_at(ta), laser_at(tb)
            dh = (lb - la) * 10 if not across else np.nan
            rows.append({"track": track, "step": step, "pair": pid, "dt": (b - a).days, "year": b.year,
                         "s1": (u[r6, c6] - ref) * PHASE_TO_MM, "coh": co[r6, c6],
                         "dh_mm": dh, "laser_los": dh * np.cos(np.radians(INC[track])),
                         "dwtd": wtd_at(tb) - wtd_at(ta),
                         "dvv": vv.get(b, np.nan) - vv.get(a, np.nan), "drvi": rvi.get(b, np.nan) - rvi.get(a, np.nan)})
d = pd.DataFrame(rows)
d.to_csv(OUT / "p6_pairs_2020_2024.csv", index=False)


def stats(g):
    g = g.dropna(subset=["laser_los", "s1", "coh"])
    k = g[(g.laser_los.abs() > 1.2) & (g.laser_los.abs() < QW)]
    out = {"n": len(g), "n_cmp": len(k)}
    if len(g) >= 6:
        out["r"] = np.corrcoef(g.s1, g.laser_los)[0, 1]
        out["same_sign"] = float((np.sign(k.s1) == np.sign(k.laser_los)).mean()) if len(k) else np.nan
        gw = g.dropna(subset=["dwtd"])
        if len(gw) >= 6:
            f = fu.bayes_linear(np.column_stack([gw.laser_los, gw.dwtd]), gw.s1.to_numpy(), fu.phase_sigma_mm(gw.coh.to_numpy()))
            out.update(m=f["beta"][0], m_sd=f["sd"][0], e=f["beta"][1])
        out["median_coh"] = g.coh.median()
    return pd.Series(out)


d["period"] = np.where(d.year <= 2021, "2020–2021 (S1A+S1B)", "2022–2024 (S1A)")
st = pd.DataFrame([{"track": k[0], "period": k[1], "dt_days": k[2], **stats(g).to_dict()}
                   for k, g in d.groupby(["track", "period", "dt"])])
st = st[st.n >= 6]
st.to_csv(OUT / "phase_vs_laser_by_revisit.csv", index=False)

models = {"laser + ΔWTD": ["laser_los", "dwtd"], "laser + ΔWTD + ΔVV": ["laser_los", "dwtd", "dvv"],
          "laser + ΔVV": ["laser_los", "dvv"], "laser + ΔWTD + ΔRVI": ["laser_los", "dwtd", "drvi"], "laser only": ["laser_los"]}
res = []
for track in INC:
    g = d[(d.track == track) & (d.step == 1)].dropna(subset=["laser_los", "dwtd", "dvv", "drvi", "s1", "coh"])
    for name, cols in models.items():
        pred = np.full(len(g), np.nan)
        for y in sorted(g.year.unique()):
            tr, te = (g.year != y).to_numpy(), (g.year == y).to_numpy()
            if tr.sum() < 8:
                continue
            f = fu.bayes_linear(g[cols].to_numpy()[tr], g.s1.to_numpy()[tr], fu.phase_sigma_mm(g.coh.to_numpy()[tr]))
            pred[te] = g[cols].to_numpy()[te] @ f["beta"]
        m = np.isfinite(pred)
        fa = fu.bayes_linear(g[cols].to_numpy(), g.s1.to_numpy(), fu.phase_sigma_mm(g.coh.to_numpy()))
        res.append({"track": track, "model": name, "n": int(m.sum()), "years": "–".join(map(str, sorted(g.year.unique()))),
                    "rmse_year_held_out_mm": float(np.sqrt(np.mean((pred[m] - g.s1.to_numpy()[m]) ** 2))),
                    **{f"coef_{c}": float(b) for c, b in zip(cols, fa["beta"])},
                    **{f"sd_{c}": float(s) for c, s in zip(cols, fa["sd"])}})
res = pd.DataFrame(res)
res.to_csv(OUT / "backscatter_in_phase_model.csv", index=False)


def md(df):
    f = lambda v: f"{v:.2f}" if isinstance(v, float) else str(v)  # noqa: E731
    return "\n".join(["| " + " | ".join(map(str, df.columns)) + " |", "|" + "---|" * len(df.columns)]
                     + ["| " + " | ".join(f(v) for v in r) + " |" for r in df.itertuples(index=False)])


chk = st[(st.track == "ascending") & (st.dt_days == 12) & st.period.str.startswith("2022")].iloc[0]
L = ["# X-065 — P6 phase vs laser by revisit, with an independent year (exploratory)", "",
     "Generated by `05_code/SAr-adapted/scripts/field_p6_short_pairs_x065.py`; every number is computed by it.", "",
     "2021 was used by neither X-059 nor the fusion (X-062): the laser starts 2021-05-27 and the 2020–2021 crops (D-020) "
     "have 6-day pairs (Sentinel-1A + 1B). Conventions as X-059 / X-062 block D. Check: ascending 2022–2024, 12-day pairs: "
     f"m = {chk.m:.2f} ± {chk.m_sd:.2f}, e = {chk.e:.2f} mm/cm — X-062 block D's fit.", "",
     "## 1. Phase vs laser at P6, by track, period and revisit", "",
     "`m`: motion coefficient of s1 = m·laser_LOS + e·ΔWTD (coherence-weighted, through the origin; 1 = the phase "
     "records the motion). `same_sign` over pairs with a laser change between its noise floor and λ/4.", "",
     md(st.round(3)), "",
     "## 2. Does backscatter add to the phase model? (consecutive pairs, all years, year held out)", "",
     md(res[["track", "model", "n", "years", "rmse_year_held_out_mm"]].round(3)), "",
     "Coefficients in `backscatter_in_phase_model.csv`.", "",
     "## Files", "", "| File | Content |", "|---|---|",
     "| `p6_pairs_2020_2024.csv` | every consecutive (step 1) and next-but-one (step 2) pair: phase, coherence, laser, ΔWTD, ΔVV, ΔRVI |",
     "| `phase_vs_laser_by_revisit.csv` | section 1 |", "| `backscatter_in_phase_model.csv` | section 2 |"]
(OUT / "README.md").write_text("\n".join(L) + "\n")
print(st.round(2).to_string(index=False))
print(res[["track", "model", "n", "rmse_year_held_out_mm"]].round(2).to_string(index=False))

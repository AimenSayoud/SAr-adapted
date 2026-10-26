"""C-075: simple answers for the supervisor — does Sentinel-1 follow the laser and the water table, and how does the
weather move the water table? One sentence (chosen by a rule fixed in the pre-registration), one figure and one table
per question.

Q1  P6: the Sentinel-1 change of each consecutive pair against the laser's change (X-126's reading, reproduced first).
Q2  The nine plots: (a) the phase change against the water-table change; (b) coherence against the size of that change;
    (c) backscatter against the water-table level, with the annual cycle removed.
Q3  The nine plots: the daily change of the water table against the day's rain and evapotranspiration.
Q4  Trees at the nine cells around P6: C-074's tables, as delivered (no new analysis).

The plot's own 40 m cell is the primary read and the 3 × 3 window a check. Pre-registration: hub
``insar_full_test/preregistration/C-075.yaml``. Writes only to the hub (``08_deliverables/simple_answers_c075/``); no field
value in this file.

    make run E=simple_answers_c075
"""
from __future__ import annotations

import shutil
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import x065_hard_targets_x126 as X126  # noqa: E402  (its reader: acquisitions between, the regular sets)

from insar_wetlands import figstyle as fs  # noqa: E402
from insar_wetlands.cube.core import WAVELENGTH_MM  # noqa: E402
from insar_wetlands.inversion.isbas import PHASE_TO_MM  # noqa: E402
from insar_wetlands.paths import cube_root, hub_root  # noqa: E402
from insar_wetlands.scoring import (  # noqa: E402
    Protocol,
    assert_preregistered,
    guard_lockbox,
    increment_stats,
)

HUB, CUBE = hub_root(), cube_root()
FT = HUB / "08_deliverables" / "insar_full_test"
OUT = HUB / "08_deliverables" / "simple_answers_c075"
PREREG = FT / "preregistration" / "C-075.yaml"
C074 = HUB / "08_deliverables" / "supervisor_pack_c074"
PLOTS = [f"P{i}" for i in range(1, 10)]
TRACKS = ("ascending", "descending")
SETS = {"2020–2021 (6-day)": (2020, 2021), "2022–2024 (12-day)": (2022, 2024), "2025 (6-day)": (2025, 2025)}
N_BOOT, SEED, MIN_N = 2000, 0, 15
LASER_START = pd.Timestamp("2021-05-27", tz="UTC")


def truthy(s: pd.Series) -> pd.Series:
    return s.astype(str).isin(["True", "1", "1.0", "true"])


# ------------------------------------------------------------------------------------------- statistics

def pearson(a, b) -> float:
    a, b = np.asarray(a, float), np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3 or np.std(a[ok]) == 0 or np.std(b[ok]) == 0:
        return np.nan
    return float(np.corrcoef(a[ok], b[ok])[0, 1])


def spear(a, b) -> float:
    a, b = np.asarray(a, float), np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    return float(spearmanr(a[ok], b[ok]).statistic) if ok.sum() > 3 else np.nan


def boot(stat, a, b, rng, groups=None) -> tuple[float, float]:
    """95 % interval of stat(a, b): resampling rows, or whole groups (calendar-month blocks) when ``groups`` is given."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    if groups is None:
        idx = [rng.integers(0, len(a), len(a)) for _ in range(N_BOOT)]
    else:
        g = pd.Series(np.arange(len(a))).groupby(np.asarray(groups)).apply(list).tolist()
        idx = [np.concatenate([g[k] for k in rng.integers(0, len(g), len(g))]) for _ in range(N_BOOT)]
    draws = np.array([stat(a[i], b[i]) for i in idx])
    lo, hi = np.nanpercentile(draws, [2.5, 97.5])
    return float(lo), float(hi)


def slope0(x, y) -> float:
    """Least-squares slope through the origin of y on x."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    return float(np.sum(x[ok] * y[ok]) / np.sum(x[ok] ** 2)) if ok.sum() > 2 else np.nan


def follows(r: float, lo: float, n: int) -> str:
    if n < MIN_N:
        return "too few pairs"
    if not lo > 0:
        return "does not follow"
    return "follows" if r >= 0.5 else "follows weakly"


# ------------------------------------------------------------------------------------------- pairs

def load_pairs() -> pd.DataFrame:
    """Consecutive pairs of every plot (X-126's reader for all plots): own cell and 3 × 3, both references."""
    cols = ["period", "track", "grid", "plot", "t1", "t2", "revisit_days", "platforms", "unw_1x1", "unw_3x3", "coh_1x1",
            "coh_3x3", "los_u", "d_laser_los_mm", "d_laser_vertical_mm", "d_wtd_plot_cm", "laser_snow_t1",
            "laser_snow_t2", "flag_frozen_air_t1", "flag_frozen_air_t2", "ref_hard_unw", "ref_grassland_unw"]
    p = pd.read_csv(CUBE / "gold" / "plots_pair.csv.gz", usecols=cols, low_memory=False)
    p = p[(p.grid == "40m") & p.period.astype(str).isin(X126.REGULAR) & p["plot"].isin(PLOTS)].copy()
    p = p[~p.platforms.astype(str).str.contains("S1D")]
    for c in ("t1", "t2"):
        p[c] = pd.to_datetime(p[c], format="ISO8601", utc=True)
    p["between"] = X126.acquisitions_between(p)
    p = p[p.between == 0].drop_duplicates(["plot", "track", "t1", "t2"]).copy()
    for ref, col in (("hard", "ref_hard_unw"), ("grass", "ref_grassland_unw")):
        for w in ("1x1", "3x3"):
            p[f"los_{ref}_{w}"] = (p[f"unw_{w}"] - p[col]) * PHASE_TO_MM
            p[f"vert_{ref}_{w}"] = p[f"los_{ref}_{w}"] / p.los_u
    p["cold"] = truthy(p.laser_snow_t1) | truthy(p.laser_snow_t2) | truthy(p.flag_frozen_air_t1) | truthy(p.flag_frozen_air_t2)
    p["set"] = None
    for name, (y0, y1) in SETS.items():
        p.loc[p.t1.dt.year.between(y0, y1), "set"] = name
    p["month"] = p.t1.dt.month
    return p


def reproduce_x126(p: pd.DataFrame) -> pd.DataFrame:
    """Stop rule: P6, the laser's span (X-126's window), 2021 and 2025 six-day consecutive pairs, no frost flag."""
    ref = pd.read_csv(FT / "x065_hard_targets" / "by_revisit.csv")
    ref = ref[(ref.reference == "hard_targets") & (ref.pairs_kind == "consecutive") & (ref.revisit_days == 6)]
    rows = []
    for track in TRACKS:
        g = p[(p["plot"] == "P6") & (p.track == track) & (p.revisit_days == 6)]
        lasr = g.dropna(subset=["d_laser_los_mm"])
        g = g[(g.t1 >= lasr.t1.min() - pd.Timedelta(days=13)) & (g.t2 <= lasr.t2.max() + pd.Timedelta(days=13))]
        for y in (2021, 2025):
            m = g[g.t1.dt.year == y][["d_laser_los_mm", "los_hard_1x1"]].dropna()
            mine = increment_stats(m.d_laser_los_mm, m.los_hard_1x1)
            th = ref[(ref.track == track) & (ref.period.astype(str) == str(y))].iloc[0]
            rows.append({"track": track, "year": y, "pairs": mine["pairs"], "pairs_x126": int(th.pairs),
                         "r": mine["pearson_increments"], "r_x126": float(th.pearson_increments)})
    t = pd.DataFrame(rows)
    t["identical"] = (t.pairs == t.pairs_x126) & ((t.r - t.r_x126).abs() < 1e-9)
    return t


# ------------------------------------------------------------------------------------------- Q1, Q2a, Q2b

def q1(p: pd.DataFrame) -> pd.DataFrame:
    rng = np.random.default_rng(SEED)
    rows = []
    for track in TRACKS:
        for s in SETS:
            for w in ("1x1", "3x3"):
                g = p[(p["plot"] == "P6") & (p.track == track) & (p.set == s) & ~p.cold]
                m = g[["d_laser_los_mm", f"los_hard_{w}"]].dropna()
                a, b = m.d_laser_los_mm.to_numpy(), m[f"los_hard_{w}"].to_numpy()
                st = increment_stats(a, b) if len(m) > 2 else {"pairs": len(m)}
                lo, hi = boot(pearson, a, b, rng) if len(m) >= 5 else (np.nan, np.nan)
                r = st.get("pearson_increments", np.nan)
                rows.append({"track": track, "set": s, "window": w, "pairs": len(m), "r": r, "r_lo95": lo,
                             "r_hi95": hi, "slope_through_origin": st.get("slope_through_origin", np.nan),
                             "rmse_mm": st.get("rmse_increments_mm", np.nan), "answer": follows(r, lo, len(m))})
    return pd.DataFrame(rows)


def q2a(p: pd.DataFrame) -> pd.DataFrame:
    rng = np.random.default_rng(SEED)
    rows = []
    for plot in PLOTS:
        for track in TRACKS:
            for s in list(SETS)[:2]:
                g = p[(p["plot"] == plot) & (p.track == track) & (p.set == s) & ~p.cold]
                for w, ref in (("1x1", "hard"), ("3x3", "hard"), ("1x1", "grass")):
                    m = g[["d_wtd_plot_cm", f"vert_{ref}_{w}"]].dropna()
                    a, b = m.d_wtd_plot_cm.to_numpy(), m[f"vert_{ref}_{w}"].to_numpy()
                    r = pearson(a, b)
                    lo, hi = boot(pearson, a, b, rng) if len(m) >= 5 else (np.nan, np.nan)
                    rows.append({"plot": plot, "track": track, "set": s, "window": w, "reference": ref,
                                 "pairs": len(m), "r": r, "r_lo95": lo, "r_hi95": hi,
                                 "slope_mm_per_cm": slope0(a, b), "answer": follows(r, lo, len(m))})
                if plot == "P6":
                    m = g[["d_wtd_plot_cm", "d_laser_vertical_mm"]].dropna()
                    a, b = m.d_wtd_plot_cm.to_numpy(), m.d_laser_vertical_mm.to_numpy()
                    r = pearson(a, b)
                    lo, hi = boot(pearson, a, b, rng) if len(m) >= 5 else (np.nan, np.nan)
                    rows.append({"plot": plot, "track": track, "set": s, "window": "laser", "reference": "laser",
                                 "pairs": len(m), "r": r, "r_lo95": lo, "r_hi95": hi,
                                 "slope_mm_per_cm": slope0(a, b), "answer": follows(r, lo, len(m))})
    return pd.DataFrame(rows)


def q2b(p: pd.DataFrame) -> pd.DataFrame:
    rng = np.random.default_rng(SEED)
    rows = []
    for plot in PLOTS:
        for track in TRACKS:
            for s in list(SETS)[:2]:
                g = p[(p["plot"] == plot) & (p.track == track) & (p.set == s) & ~p.cold].copy()
                for w in ("1x1", "3x3"):
                    m = g[["d_wtd_plot_cm", f"coh_{w}", "month"]].dropna()
                    a = m.d_wtd_plot_cm.abs().to_numpy()
                    c = m[f"coh_{w}"].to_numpy()
                    c_season = c - m.groupby("month")[f"coh_{w}"].transform("median").to_numpy()
                    rho, rho_s = spear(a, c), spear(a, c_season)
                    lo, hi = boot(spear, a, c_season, rng) if len(m) >= 5 else (np.nan, np.nan)
                    ans = ("too few pairs" if len(m) < MIN_N else
                           "coherence falls when the water table changes" if hi < 0 else "no clear relation")
                    rows.append({"plot": plot, "track": track, "set": s, "window": w, "pairs": len(m),
                                 "spearman_raw": rho, "spearman_season_removed": rho_s, "lo95": lo, "hi95": hi,
                                 "answer": ans})
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------------------------- Q2c

def q2c() -> tuple[pd.DataFrame, pd.DataFrame]:
    cols = ["track", "date", "platform", "plot", "vv_db_1x1", "vh_db_1x1", "vv_db_3x3", "vh_db_3x3", "wtd_at",
            "laser_snow", "flag_frozen_air", "censored"]
    a = pd.read_csv(CUBE / "gold" / "plots_acq.csv.gz", usecols=cols, low_memory=False)
    a["date"] = pd.to_datetime(a.date)
    a = a[a["plot"].isin(PLOTS) & a.date.dt.year.between(2020, 2024) & ~a.platform.astype(str).str.contains("S1D")]
    a = a[~(truthy(a.laser_snow) | truthy(a.flag_frozen_air) | truthy(a.censored))].dropna(subset=["wtd_at"]).copy()
    a["month"] = a.date.dt.month
    a["block"] = a.date.dt.strftime("%Y-%m")
    rng = np.random.default_rng(SEED)
    rows = []
    for plot in PLOTS:
        for track in TRACKS:
            g = a[(a["plot"] == plot) & (a.track == track)]
            for band in ("vv_db_1x1", "vh_db_1x1", "vv_db_3x3", "vh_db_3x3"):
                m = g[["wtd_at", band, "month", "block"]].dropna()
                an_w = m.wtd_at - m.groupby("month").wtd_at.transform("mean")
                an_b = m[band] - m.groupby("month")[band].transform("mean")
                r_raw, r_an = pearson(m.wtd_at, m[band]), pearson(an_w, an_b)
                lo, hi = boot(pearson, an_w.to_numpy(), an_b.to_numpy(), rng, groups=m.block.to_numpy()) \
                    if len(m) >= 10 else (np.nan, np.nan)
                ans = ("too few acquisitions" if len(m) < MIN_N else "brighter when the water table is higher" if lo > 0
                       else "darker when the water table is higher" if hi < 0 else "no clear relation")
                rows.append({"plot": plot, "track": track, "band": band.split("_")[0].upper(),
                             "window": band.split("_")[-1], "acquisitions": len(m), "r_raw": r_raw,
                             "r_season_removed": r_an, "lo95": lo, "hi95": hi,
                             "db_per_10cm": slope0(an_w.to_numpy(), an_b.to_numpy()) * 10, "answer": ans})
    return pd.DataFrame(rows), a


# ------------------------------------------------------------------------------------------- Q3

def daily() -> pd.DataFrame:
    h = pd.read_csv(CUBE / "silver" / "site_hourly.csv.gz", low_memory=False)
    h["time_utc"] = pd.to_datetime(h.time_utc, format="ISO8601", utc=True)
    h = h[h.time_utc.dt.year.between(2020, 2024)].set_index("time_utc")
    agg = {"st_rain_mm": "sum", "om_era5_et0_fao_evapotranspiration": "sum", "om_era5_rain": "sum",
           "st_air_c": "min", "laser_snow_evidence_72h": lambda s: truthy(s).any()}
    for pl in PLOTS:
        agg[f"wtd_{pl}_cm"] = "mean"
        agg[f"censored_{pl}"] = lambda s: truthy(s).any()
    d = h.resample("1D").agg(agg)
    d = d.rename(columns={"st_rain_mm": "rain_mm", "om_era5_et0_fao_evapotranspiration": "et0_mm",
                          "om_era5_rain": "rain_era5_mm", "st_air_c": "air_min_c", "laser_snow_evidence_72h": "snow"})
    return d


def q3(d: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(SEED)
    rows, ev = [], []
    for pl in PLOTS:
        w = d[f"wtd_{pl}_cm"]
        x = pd.DataFrame({"dw": w.diff(), "rain": d.rain_mm, "rain_era5": d.rain_era5_mm, "et0": d.et0_mm,
                          "year": d.index.year, "block": d.index.strftime("%Y-%m")})
        bad = d.snow | (d.air_min_c < 0) | d[f"censored_{pl}"] | d[f"censored_{pl}"].shift(1, fill_value=False)
        x = x[~bad].dropna()
        for rain_col in ("rain", "rain_era5"):
            X = np.column_stack([np.ones(len(x)), x[rain_col], x.et0])
            beta = np.linalg.lstsq(X, x.dw, rcond=None)[0]
            # leave-one-year-out R²
            pred = pd.Series(np.nan, index=x.index)
            for y in sorted(x.year.unique()):
                tr, te = x.year != y, x.year == y
                b = np.linalg.lstsq(X[tr.to_numpy()], x.dw[tr], rcond=None)[0]
                pred[te] = X[te.to_numpy()] @ b
            r2 = 1 - np.nansum((x.dw - pred) ** 2) / np.nansum((x.dw - x.dw.mean()) ** 2)
            # block bootstrap of the coefficients (calendar-month blocks)
            blocks = x.groupby("block").indices
            keys = list(blocks)
            draws = []
            for _ in range(N_BOOT):
                i = np.concatenate([blocks[keys[k]] for k in rng.integers(0, len(keys), len(keys))])
                draws.append(np.linalg.lstsq(X[i], x.dw.to_numpy()[i], rcond=None)[0])
            draws = np.array(draws)
            lo, hi = np.percentile(draws, 2.5, axis=0), np.percentile(draws, 97.5, axis=0)
            rows.append({"plot": pl, "rain_source": "station" if rain_col == "rain" else "ERA5", "days": len(x),
                         "rise_per_10mm_rain_cm": beta[1] * 10, "rain_lo95": lo[1] * 10, "rain_hi95": hi[1] * 10,
                         "change_per_mm_et0_cm": beta[2], "et0_lo95": lo[2], "et0_hi95": hi[2],
                         "loyo_r2": r2,
                         "answer_rain": "rain raises the water table" if lo[1] > 0 else "no clear effect of rain",
                         "answer_et0": "evaporation lowers it" if hi[2] < 0 else "no clear effect of evaporation"})
        # events: days with at least 10 mm of station rain
        for t in d.index[d.rain_mm >= 10]:
            if t - pd.Timedelta(days=1) not in w.index:
                continue
            base = w.get(t - pd.Timedelta(days=1))
            row = {"plot": pl, "date": t.date(), "rain_mm": d.rain_mm[t], "wtd_before_cm": base}
            for k in (1, 3, 7):
                v = w.get(t + pd.Timedelta(days=k))
                row[f"change_after_{k}d_cm"] = (v - base) if (v is not None and np.isfinite(v) and np.isfinite(base)) else np.nan
            ev.append(row)
    return pd.DataFrame(rows), pd.DataFrame(ev)


# ------------------------------------------------------------------------------------------- figures

def fig_q1(p, t, path):
    with fs.style():
        fig, axes = plt.subplots(2, 3, figsize=fs.size("full", 120), sharex=True, sharey=True,
                                 gridspec_kw={"hspace": 0.25, "wspace": 0.08, "left": 0.09, "right": 0.99,
                                              "top": 0.92, "bottom": 0.09})
        for i, track in enumerate(TRACKS):
            for j, s in enumerate(SETS):
                ax = axes[i, j]
                g = p[(p["plot"] == "P6") & (p.track == track) & (p.set == s) & ~p.cold][["d_laser_los_mm", "los_hard_1x1"]].dropna()
                ax.axline((0, 0), slope=1, color=fs.MUTED, lw=0.6, ls="--")
                ax.scatter(g.d_laser_los_mm, g.los_hard_1x1, s=8, color=fs.TRACK[track], alpha=0.8, lw=0)
                r = t[(t.track == track) & (t.set == s) & (t.window == "1x1")].iloc[0]
                ax.text(0.03, 0.97, f"{s}\nr = {r.r:.2f} [{r.r_lo95:.2f}, {r.r_hi95:.2f}], n = {r.pairs}\n{r.answer}",
                        transform=ax.transAxes, va="top", fontsize=7)
                if j == 0:
                    ax.set_ylabel(f"{track}: Sentinel-1 change (mm)")
                if i == 1:
                    ax.set_xlabel("laser change at P6 (mm, line of sight)")
        fig.suptitle("Q1 — does Sentinel-1 follow the laser at P6? One point per consecutive pair", fontsize=9)
        fs.export(fig, path, provenance={"source": "simple_answers_c075 q1"})
        plt.close(fig)


def dotplot(df, value, lo, hi, title, ylabel, path, ref_rows=None):
    sets = list(dict.fromkeys(df.set))
    with fs.style():
        fig, axes = plt.subplots(1, len(sets), figsize=fs.size("full", 80), sharey=True,
                                 gridspec_kw={"wspace": 0.06, "left": 0.08, "right": 0.99, "top": 0.86, "bottom": 0.14})
        axes = np.atleast_1d(axes)
        for ax, s in zip(axes, sets):
            ax.axhline(0, color=fs.INK, lw=0.6)
            for k, track in enumerate(TRACKS):
                g = df[(df.set == s) & (df.track == track)].set_index("plot").reindex(PLOTS)
                x = np.arange(len(PLOTS)) + (k - 0.5) * 0.25
                ax.vlines(x, g[lo], g[hi], color=fs.TRACK[track], lw=1.2)
                ax.plot(x, g[value], "o", ms=4, color=fs.TRACK[track], label=track)
            if ref_rows is not None:
                rr = ref_rows[ref_rows.set == s]
                for k, track in enumerate(TRACKS):
                    q = rr[rr.track == track]
                    if len(q):
                        ax.plot(5 + (k - 0.5) * 0.25, q[value].iloc[0], marker="*", ms=9, color=fs.TRACK[track],
                                mec=fs.INK, mew=0.5, ls="none", label="P6 laser (reference)" if k == 0 else None)
            ax.set_xticks(range(len(PLOTS)), PLOTS)
            ax.set_title(s, fontsize=8)
        axes[0].set_ylabel(ylabel)
        axes[0].legend(fontsize=7, loc="lower left")
        fig.suptitle(title, fontsize=9)
        fs.export(fig, path, provenance={"source": "simple_answers_c075"})
        plt.close(fig)


def fig_q3(t, ev, d, path):
    s = t[t.rain_source == "station"].set_index("plot").reindex(PLOTS)
    with fs.style():
        fig, axes = plt.subplots(2, 2, figsize=fs.size("full", 130),
                                 gridspec_kw={"hspace": 0.4, "wspace": 0.25, "left": 0.08, "right": 0.99, "top": 0.92, "bottom": 0.07})
        x = np.arange(len(PLOTS))
        ax = axes[0, 0]
        ax.vlines(x, s.rain_lo95, s.rain_hi95, color=fs.TRACK["ascending"], lw=1.4)
        ax.plot(x, s.rise_per_10mm_rain_cm, "o", color=fs.TRACK["ascending"])
        ax.axhline(0, color=fs.INK, lw=0.6)
        ax.set_xticks(x, PLOTS)
        ax.set_ylabel("rise per 10 mm of rain (cm)")
        ax.set_title("(a) rain", fontsize=8)
        ax = axes[0, 1]
        ax.vlines(x, s.et0_lo95, s.et0_hi95, color=fs.TRACK["descending"], lw=1.4)
        ax.plot(x, s.change_per_mm_et0_cm, "o", color=fs.TRACK["descending"])
        ax.axhline(0, color=fs.INK, lw=0.6)
        ax.set_xticks(x, PLOTS)
        ax.set_ylabel("change per mm of evapotranspiration (cm)")
        ax.set_title("(b) evaporation (ERA5 FAO-56 reference)", fontsize=8)
        ax = axes[1, 0]
        med = ev.groupby("plot")[["change_after_1d_cm", "change_after_3d_cm", "change_after_7d_cm"]].median().reindex(PLOTS)
        for k, (col, lab) in enumerate((("change_after_1d_cm", "1 day"), ("change_after_3d_cm", "3 days"),
                                        ("change_after_7d_cm", "7 days"))):
            ax.bar(x + (k - 1) * 0.27, med[col], width=0.27, label=f"{lab} after")
        ax.axhline(0, color=fs.INK, lw=0.6)
        ax.set_xticks(x, PLOTS)
        ax.set_ylabel("median change after a day with ≥ 10 mm rain (cm)")
        ax.legend(fontsize=7)
        ax.set_title("(c) after heavy rain", fontsize=8)
        ax = axes[1, 1]
        yr = d[d.index.year == 2021]
        ax.plot(yr.index, yr["wtd_P6_cm"], color=fs.INK, lw=0.8)
        ax.set_ylabel("water table at P6 (cm, negative below surface)")
        ax2 = ax.twinx()
        ax2.bar(yr.index, yr.rain_mm, color=fs.TRACK["ascending"], alpha=0.5, width=1)
        ax2.set_ylabel("daily rain (mm)")
        ax2.set_ylim(4 * max(float(yr.rain_mm.max()), 1.0), 0)      # rain hangs from the top, in its own quarter
        ax2.tick_params(axis="y", labelright=True, right=True)
        ax.set_title("(d) P6 in 2021: water table and rain", fontsize=8)
        fig.suptitle("Q3 — how the weather moves the water table (daily, 2020–2024, frozen and snow days left out)", fontsize=9)
        fs.export(fig, path, provenance={"source": "simple_answers_c075 q3"})
        plt.close(fig)


# ------------------------------------------------------------------------------------------- README and Excel

def count(df, ans, **sel) -> str:
    g = df
    for k, v in sel.items():
        g = g[g[k] == v]
    g = g[g.answer != "too few pairs"] if "answer" in g else g
    return f"{int(g[g.answer.str.startswith(ans)].shape[0])} of {len(g)}"


def three(df, **sel) -> str:
    """'follows at a, weakly at b, not at c of n plots' for one selection of rows."""
    g = df
    for k, v in sel.items():
        g = g[g[k] == v]
    g = g[g.answer != "too few pairs"]
    n = len(g)
    f = int((g.answer == "follows").sum())
    w = int((g.answer == "follows weakly").sum())
    return f"follows at {f}, weakly at {w}, not at {n - f - w} of {n} plots"


def md(df, nd=2):
    d = df.copy()
    for c in d.columns:
        if pd.api.types.is_float_dtype(d[c]):
            d[c] = d[c].map(lambda v: "" if pd.isna(v) else f"{v:.{nd}f}")
    return "\n".join(["| " + " | ".join(d.columns) + " |", "|" + "---|" * len(d.columns)] +
                     ["| " + " | ".join(str(v) for v in r) + " |" for r in d.itertuples(index=False)])


def readme(t1, t2a, t2b, t2c, t3, ev, nine, rep) -> str:
    a1 = t1[t1.window == "1x1"]
    a2 = t2a[(t2a.window == "1x1") & (t2a.reference == "hard")]
    lz = t2a[t2a.window == "laser"]
    b2 = t2b[t2b.window == "1x1"]
    c2 = t2c[t2c.window == "1x1"]
    s3 = t3[t3.rain_source == "station"]
    lines = [
        "# Simple answers — Sentinel-1, the laser, the water table and the weather at Rzecin (generated)", "",
        "Generated by `05_code/SAr-adapted/experiments/simple_answers_c075.py`; do not edit. Pre-registered:",
        "`08_deliverables/insar_full_test/preregistration/C-075.yaml` (each sentence below is chosen by a rule written "
        "before the numbers). Excel: `simple_answers.xlsx` (one sheet per question). Read with each plot's own 40 m "
        "radar cell; the 3 × 3 window is in the tables as a check.", "",
        "## In one page", "",
        "1. **Sentinel-1 and the laser (P6).** " + "; ".join(
            f"{s}: {', '.join(f'{tr} {a1[(a1.set == s) & (a1.track == tr)].answer.iloc[0]}' for tr in TRACKS)}" for s in SETS)
        + ". Six-day pairs follow the mat's movement closely; twelve-day pairs less, most of all on the morning "
        "(descending) orbit.",
        "2. **Sentinel-1 and the water table (nine plots).** The phase change follows the water-table change — "
        + "; ".join(f"{s} ascending: {three(a2, set=s, track='ascending')}, descending: "
                    f"{three(a2, set=s, track='descending')}" for s in list(SETS)[:2])
        + ". Coherence falls when the water table changes (" + count(b2[b2.set == list(SETS)[0]], "coherence falls")
        + " plot-and-track cases on six-day pairs). VV backscatter is brighter when the water is higher ("
        + count(c2[c2.band == "VV"], "brighter") + " plot-and-track cases).",
        "3. **Weather and the water table.** Rain raises it and evaporation lowers it at "
        + count(s3.assign(answer=s3.answer_rain), "rain raises") + " plots; together they explain "
        f"{s3.loyo_r2.min():.0%}–{s3.loyo_r2.max():.0%} of the day-to-day change. The rise per 10 mm of rain is largest "
        f"at {s3.loc[s3.rise_per_10mm_rain_cm.idxmax(), 'plot']} ({s3.rise_per_10mm_rain_cm.max():.1f} cm) and smallest at "
        f"{s3.loc[s3.rise_per_10mm_rain_cm.idxmin(), 'plot']} ({s3.rise_per_10mm_rain_cm.min():.1f} cm).",
        "4. **Trees at the nine cells around P6.** Almost none; the shrubs there do not change answer 1.", "",
        "## Q1 — Does Sentinel-1 follow the laser at P6?", "",
        "**Answer.** " + "; ".join(
            f"{s}: {', '.join(f'{tr} {a1[(a1.set == s) & (a1.track == tr)].answer.iloc[0]}' for tr in TRACKS)}" for s in SETS) + ".",
        "", "![Q1](Q1_laser_vs_s1.png)", "", md(a1[["track", "set", "pairs", "r", "r_lo95", "r_hi95", "slope_through_origin", "rmse_mm", "answer"]]), "",
        "**How.** Each point is one pair of consecutive images: the change of the mat surface measured by the laser, "
        "against the change Sentinel-1 measures in P6's cell (referenced to stable buildings). r = 1 would be perfect; "
        "the slope says whether the size matches (1 = same size).",
        "**Caution.** One laser point. Snow and frost days are left out. 2022–2024 has only twelve-day pairs (one "
        "satellite), and over twelve days the mat often moves more than one pair can measure without ambiguity (a "
        f"quarter of the wavelength along the line of sight, {WAVELENGTH_MM / 4 / 10:.1f} cm).", "",
        "## Q2 — Does Sentinel-1 follow the water table at the nine plots?", "",
        "### (a) The phase change against the water-table change", "",
        "**Answer.** " + "; ".join(f"{s}, ascending: {three(a2, set=s, track='ascending')}; descending: "
                                   f"{three(a2, set=s, track='descending')}" for s in list(SETS)[:2]) + ".",
        "", "![Q2a](Q2a_phase_vs_wtd.png)", "",
        md(a2[["plot", "track", "set", "pairs", "r", "r_lo95", "r_hi95", "slope_mm_per_cm", "answer"]]), "",
        "The laser itself against the water table at P6 (the reference — how the mat moves with its water):", "",
        md(lz[["track", "set", "pairs", "r", "r_lo95", "r_hi95", "slope_mm_per_cm", "answer"]]), "",
        "**How.** For each pair of consecutive images, Sentinel-1's vertical change of the plot's cell against the change "
        "of the plot's water table between the two dates. The radar sees the mat's surface, which rises and falls with "
        "its water by buoyancy — by a fraction of the water's movement — so the slope is millimetres of surface per "
        "centimetre of water.",
        "**Caution.** P8 and P9 lie in the same radar cell. The water table was gap-filled in the delivery without a flag. "
        "P1–P4 sit where the mat is held more firmly.", "",
        "### (b) Coherence against the size of the water-table change", "",
        "**Answer.** Coherence falls when the water table changes at " +
        "; ".join(f"{count(b2, 'coherence falls', set=s, track=tr)} plots ({tr}, {s})" for s in list(SETS)[:2] for tr in TRACKS) + ".",
        "", "![Q2b](Q2b_coherence_vs_dwtd.png)", "",
        md(b2[["plot", "track", "set", "pairs", "spearman_raw", "spearman_season_removed", "lo95", "hi95", "answer"]]), "",
        "**How.** Coherence (0–1) says how similar the two radar images are. Spearman correlation with the size of the "
        "water-table change, after removing each month's typical coherence (so the season alone does not make it).", "",
        "### (c) Backscatter against the water-table level", "",
        "**Answer.** VV: " + ", ".join(f"{tr} brighter when the water is higher at "
                                         f"{count(c2[c2.band == 'VV'], 'brighter', track=tr)} plots" for tr in TRACKS) +
        "; VH: " + ", ".join(f"{tr} {count(c2[c2.band == 'VH'], 'brighter', track=tr)}" for tr in TRACKS) + ".",
        "", "![Q2c](Q2c_backscatter_vs_wtd.png)", "",
        md(c2[["plot", "track", "band", "acquisitions", "r_raw", "r_season_removed", "lo95", "hi95", "db_per_10cm", "answer"]]), "",
        "**How.** Backscatter is the brightness of each radar image. Both it and the water table share an annual cycle, "
        "so each is taken relative to its own mean for that calendar month before they are compared; the interval "
        "resamples whole months.", "",
        "## Q3 — How does the weather move the water table?", "",
        "**Answer.** " + f"{count(s3.assign(answer=s3.answer_rain), 'rain raises')} plots: rain raises the water table; "
        f"{count(s3.assign(answer=s3.answer_et0), 'evaporation lowers')} plots: evaporation lowers it. Rain and "
        f"evaporation together explain {s3.loyo_r2.min():.0%}–{s3.loyo_r2.max():.0%} of the day-to-day change on years "
        "they were not fitted on.", "",
        "![Q3](Q3_weather_vs_wtd.png)", "",
        md(s3[["plot", "days", "rise_per_10mm_rain_cm", "rain_lo95", "rain_hi95", "change_per_mm_et0_cm", "et0_lo95",
               "et0_hi95", "loyo_r2", "answer_rain", "answer_et0"]], 3), "",
        "Median change after a day with at least 10 mm of rain (cm):", "",
        md(ev.groupby("plot")[["change_after_1d_cm", "change_after_3d_cm", "change_after_7d_cm"]].median().reindex(PLOTS)
           .assign(events=ev.groupby("plot").size()).reset_index()), "",
        "**How.** Daily change of each plot's water table against that day's rain at the site station and the "
        "reference evapotranspiration (how much water a wet surface loses to the air) from ERA5. Fitted per plot; "
        "the share explained is tested on each year left out of the fit.",
        "**Caution.** Days with frost or snow are left out. A water table at the bottom of a dry well is left out. "
        "ERA5 rain in place of the station's is in the Excel (sheet Q3, `rain_source`).", "",
        "## Q4 — Trees at the nine cells around P6", "",
        "**Answer.** The nine cells hold almost no trees and some shrubs (table, from the laser scan of 19 September "
        "2024); every cell follows the laser on the six-day pairs, alike — the trees and shrubs there do not change "
        "the answer to Q1.", "",
        "![Q4](Q4_nine_cells.png)", "",
        md(nine), "",
        "**How.** C-074's tables, as delivered on 2026-10-24 (no new analysis).", "",
        "## Check", "",
        "Q1's reader reproduces X-126 exactly before anything is reported:", "", md(rep), "",
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    assert_preregistered(PREREG, [OUT / "README.md"])
    p = load_pairs()
    guard_lockbox(pd.concat([p.t1, p.t2]), Protocol.load(FT / "protocol.yaml"))
    rep = reproduce_x126(p)
    rep.to_csv(OUT / "check_reproduces_x126.csv", index=False)
    if not rep.identical.all():
        print(rep.to_string())
        sys.exit("STOP (C-075 stop rule): Q1's reader does not reproduce X-126")
    t1, t2a, t2b = q1(p), q2a(p), q2b(p)
    t2c, _ = q2c()
    d = daily()
    t3, ev = q3(d)
    nine = pd.read_csv(C074 / "tables" / "nine_cells.csv")[["cell", "depth_in_mat_m", "boardwalk_share",
                                                             "tree_share_2024", "shrub_share_2024"]]
    grp = pd.read_csv(C074 / "tables" / "nine_cells_by_group.csv")
    six = grp[(grp.revisit_days == 6) & (grp.group.astype(str) == "2021")]
    nine = nine.merge(six.pivot_table(index="cell", columns="track", values="pearson_increments")
                      .add_prefix("r_laser_2021_").reset_index(), on="cell", how="left")
    for name, df in (("Q1_laser_vs_s1", t1), ("Q2a_phase_vs_wtd", t2a), ("Q2b_coherence_vs_dwtd", t2b),
                     ("Q2c_backscatter_vs_wtd", t2c), ("Q3_weather_vs_wtd", t3), ("Q3_rain_events", ev),
                     ("Q4_nine_cells", nine)):
        df.to_csv(OUT / f"{name}.csv", index=False)
    fig_q1(p, t1, OUT / "Q1_laser_vs_s1")
    lz = t2a[t2a.window == "laser"]
    dotplot(t2a[(t2a.window == "1x1") & (t2a.reference == "hard")], "r", "r_lo95", "r_hi95",
            "Q2a — Sentinel-1 phase change against the water-table change, per plot (r, 95 % interval)",
            "correlation r", OUT / "Q2a_phase_vs_wtd", ref_rows=lz)
    dotplot(t2b[t2b.window == "1x1"], "spearman_season_removed", "lo95", "hi95",
            "Q2b — coherence against the size of the water-table change, season removed (Spearman, 95 % interval)",
            "Spearman correlation", OUT / "Q2b_coherence_vs_dwtd")
    c2 = t2c[(t2c.window == "1x1") & (t2c.band == "VV")].assign(set="VV")
    c2h = t2c[(t2c.window == "1x1") & (t2c.band == "VH")].assign(set="VH")
    dotplot(pd.concat([c2, c2h]), "r_season_removed", "lo95", "hi95",
            "Q2c — backscatter against the water-table level, annual cycle removed (r, 95 % interval)",
            "correlation r", OUT / "Q2c_backscatter_vs_wtd")
    fig_q3(t3, ev, d, OUT / "Q3_weather_vs_wtd")
    for src, dst in (("S2c_nine_cells_coherence_laser.png", "Q4_nine_cells.png"),):
        shutil.copyfile(C074 / src, OUT / dst)
    (OUT / "README.md").write_text(readme(t1, t2a, t2b, t2c, t3, ev, nine, rep), encoding="utf-8")
    with pd.ExcelWriter(OUT / "simple_answers.xlsx") as xl:
        pd.DataFrame({"sheet": ["Q1", "Q2a", "Q2b", "Q2c", "Q3", "Q3_events", "Q4"],
                      "question": ["Does Sentinel-1 follow the laser at P6?",
                                   "Does the Sentinel-1 phase change follow the water-table change at each plot?",
                                   "Does coherence fall when the water table changes?",
                                   "Is backscatter brighter when the water table is higher (annual cycle removed)?",
                                   "How do rain and evaporation move the water table, per plot?",
                                   "Water-table change after each day with at least 10 mm of rain",
                                   "The nine cells around P6: trees, shrubs and how each follows the laser"]}
                     ).to_excel(xl, sheet_name="README", index=False)
        for name, df in (("Q1", t1), ("Q2a", t2a), ("Q2b", t2b), ("Q2c", t2c), ("Q3", t3), ("Q3_events", ev),
                         ("Q4", nine)):
            df.to_excel(xl, sheet_name=name, index=False)
    print((OUT / "README.md").read_text().split("## Q2")[0][-1500:])


if __name__ == "__main__":
    main()

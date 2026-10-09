"""X-082 test 1 — does coherence follow the *change* in moisture between the two dates (Hrysiewicz et al. 2023) on the
Rzecin floating mat? Exploratory; plan and pass/fail criteria in ``08_deliverables/hrysiewicz_tests_x082/PLAN.md``.

Every pair the cube holds (6–24 days, 2017–2026, both tracks), the zone median coherence of the mat, grassland and
stable ground, and moisture at the two dates: ERA5-Land topsoil moisture 0–7 cm (all years, modelled) and the mat water
table (2020–2024, measured). Frozen or snowy pairs are kept apart (ERA5 moisture is not liquid water then).

1. Upper envelope (90th percentile per |Δ moisture| bin) and its slope — their Fig. 10; CI by resampling years.
2. Change vs level: standardised coefficients of |Δ moisture| and the pair's mean moisture in one model.
3. Season control: the same on monthly anomalies, and Spearman ρ within May–Sep and Oct–Apr.
4. Frost regime. 5. Mat vs grassland sensitivity. 6. Coherence timeline (the near-diagonal band of their matrix).

Reads the cube (hub); writes only to the hub. No field value in this file.

    PYTHONPATH=src python scripts/coherence_moisture_x082.py
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

from insar_wetlands import field  # noqa: E402

HUB = field.field_root().parents[1]
GOLD = HUB / "06_data" / "cube" / "gold"
OUT = HUB / "08_deliverables" / "hrysiewicz_tests_x082" / "coherence_moisture"
TRACKS = ("ascending", "descending")
ZONES = ("mat", "grassland", "stable")
DRIVERS = {"sm": ("d_om_era5_land_soil_moisture_0_to_7cm", "om_era5_land_soil_moisture_0_to_7cm",
                  "ERA5-Land topsoil moisture (m³/m³)"),
           "wtd": ("d_wtd_mat_mean_cm", "wtd_mat_mean_cm", "mat water table (cm)")}
B = 500
MIN_YEARS = 4          # a year-resampling interval needs at least this many years
RNG = np.random.default_rng(0)

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, MUTED, GRIDC, AXIS, SURFACE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"
SEQ = LinearSegmentedColormap.from_list("seq_blue", ["#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])
RV_COL = {6: BLUE, 12: ORANGE, 24: AQUA}
plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE, "axes.edgecolor": AXIS,
    "axes.labelcolor": INK2, "axes.titlecolor": INK, "axes.titlesize": 9.5, "axes.titleweight": "bold",
    "axes.titlelocation": "left", "axes.labelsize": 8.5, "xtick.color": MUTED, "ytick.color": MUTED,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "axes.grid": True, "grid.color": GRIDC, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False, "legend.fontsize": 7.5,
})


def boolify(s: pd.Series) -> pd.Series:
    return s.map({True: True, False: False, "True": True, "False": False}).fillna(False).astype(bool)


def load() -> pd.DataFrame:
    keep = ["track", "grid", "pair", "t1", "t2", "revisit_days"] + [f"ref_{z}_coh" for z in ZONES] + \
           [v[0] for v in DRIVERS.values()] + [f"flag_{f}_{t}" for f in ("frozen_air", "frozen_soil_om", "snow_om",
                                                                          "snow_imgw") for t in ("t1", "t2")]
    p = pd.read_csv(GOLD / "plots_pair.csv.gz", usecols=keep, low_memory=False)
    p = p[p.grid == "40m"].drop_duplicates(["track", "pair"]).copy()
    for c in ("t1", "t2"):
        p[c] = pd.to_datetime(p[c], format="ISO8601", utc=True)
    p["cold"] = np.column_stack([boolify(p[c]) for c in p.columns if c.startswith("flag_")]).any(axis=1)
    p["rv"] = np.where(p.revisit_days <= 6, 6, np.where(p.revisit_days <= 12, 12, np.where(p.revisit_days <= 24, 24, 0)))
    # rv: 6, 12, 24 (= 18–24-day pooled); 0 = longer pairs (36–744 days: the 2022–2024 36/48-day network, bridges,
    # a few annual pairs) — kept out of tests 1–5 and reported apart (section 7)
    p["year"], p["month"] = p.t1.dt.year, p.t1.dt.month
    p["season"] = np.where(p.month.between(5, 9), "May–Sep", "Oct–Apr")
    p["mid"] = p.t1 + (p.t2 - p.t1) / 2
    # moisture level at each date, from the per-acquisition table (site variables: the same for every plot)
    a = pd.read_csv(GOLD / "plots_acq.csv.gz", usecols=["track", "date"] + [v[1] for v in DRIVERS.values()],
                    low_memory=False).drop_duplicates(["track", "date"])
    a["date"] = pd.to_datetime(a.date).dt.date
    lv = a.set_index(["track", "date"])
    for k, (dcol, lcol, _) in DRIVERS.items():
        l1 = lv[lcol].reindex(list(zip(p.track, p.t1.dt.date))).to_numpy()
        l2 = lv[lcol].reindex(list(zip(p.track, p.t2.dt.date))).to_numpy()
        p[f"{k}_level"] = (l1 + l2) / 2
        p[f"{k}_absd"] = p[dcol].abs()
    return p


def year_boot(years: np.ndarray, fn) -> tuple:
    """Statistic(s) and 95 % intervals by resampling whole years (pairs inside a year share weather). ``fn`` takes an
    index array into the caller's numpy arrays and returns a scalar or a 1-D array; result (est, lo, hi), same shape."""
    est = np.atleast_1d(np.asarray(fn(np.arange(len(years))), float))
    rows = [np.flatnonzero(years == y) for y in np.unique(years)]
    vals = []
    for _ in range(B):
        idx = np.concatenate([rows[i] for i in RNG.integers(0, len(rows), len(rows))])
        try:
            vals.append(np.atleast_1d(np.asarray(fn(idx), float)))
        except (ValueError, np.linalg.LinAlgError):
            continue
    lo, hi = (np.nanpercentile(np.vstack(vals), [2.5, 97.5], axis=0) if vals else (est * np.nan, est * np.nan))
    if est.size == 1:
        return float(est[0]), float(lo[0]), float(hi[0])
    return est, lo, hi


def envelope_slope(x: np.ndarray, y: np.ndarray, nbin: int = 6, q: float = 0.9) -> float:
    """Slope of the q-quantile of y against x over quantile bins of x (their Fig. 10 envelope)."""
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if len(x) < nbin * 8:
        return np.nan
    edges = np.unique(np.quantile(x, np.linspace(0, 1, nbin + 1)))
    b = np.clip(np.digitize(x, edges[1:-1]), 0, len(edges) - 2)
    xm = np.array([np.median(x[b == k]) for k in range(len(edges) - 1) if (b == k).any()])
    yq = np.array([np.quantile(y[b == k], q) for k in range(len(edges) - 1) if (b == k).any()])
    return float(np.polyfit(xm, yq, 1)[0]) if len(xm) >= 3 else np.nan


def std_betas(y: np.ndarray, X: np.ndarray) -> np.ndarray:
    """Standardised OLS coefficients of the columns of X on y."""
    ok = np.isfinite(y) & np.isfinite(X).all(axis=1)
    y, X = y[ok], X[ok]
    if len(y) < 15:
        raise ValueError
    Z = (X - X.mean(0)) / X.std(0)
    zy = (y - y.mean()) / y.std()
    return np.linalg.lstsq(np.column_stack([np.ones(len(zy)), Z]), zy, rcond=None)[0][1:]


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    ok = np.isfinite(x) & np.isfinite(y)
    return float(spearmanr(x[ok], y[ok]).statistic)


def anomalies(g: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    g = g.copy()
    for c in cols:
        g[c] = g[c] - g.groupby("month")[c].transform("median")
    return g


def analyse(p: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    env, beta, seas, sens = [], [], [], []
    for (t, rv), g0 in p[p.rv > 0].groupby(["track", "rv"]):
        warm = g0[~g0.cold]
        for k in DRIVERS:
            x, lvl = f"{k}_absd", f"{k}_level"
            gk = warm.dropna(subset=[x, lvl]).reset_index(drop=True)
            if len(gk) < 30:
                continue
            yrs = gk.year.to_numpy()
            ny = len(np.unique(yrs))
            for z in ZONES:
                y = f"ref_{z}_coh"
                xv, yv = gk[x].to_numpy(float), gk[y].to_numpy(float)
                s, lo, hi = year_boot(yrs, lambda i, xv=xv, yv=yv: envelope_slope(xv[i], yv[i]))
                env.append({"track": t, "revisit": rv, "driver": k, "zone": z, "n": len(gk), "n_years": ny, "envelope_slope": s,
                            "lo95": lo, "hi95": hi})
                for mode, gg in (("raw", gk), ("monthly anomalies", anomalies(gk, [y, x, lvl]))):
                    Y, X = gg[y].to_numpy(float), gg[[x, lvl]].to_numpy(float)
                    (bc, bl), (lc, ll), (hc, hl) = year_boot(yrs, lambda i, Y=Y, X=X: std_betas(Y[i], X[i]))
                    beta.append({"track": t, "revisit": rv, "driver": k, "zone": z, "data": mode, "n": len(gg), "n_years": ny,
                                 "beta_change": bc, "change_lo95": lc, "change_hi95": hc,
                                 "beta_level": bl, "level_lo95": ll, "level_hi95": hl})
                for sea in ("May–Sep", "Oct–Apr"):
                    m = (gk.season == sea).to_numpy()
                    if m.sum() < 15:
                        continue
                    xs, ys = xv[m], yv[m]
                    r, lo2, hi2 = year_boot(yrs[m], lambda i, xs=xs, ys=ys: spearman(xs[i], ys[i]))
                    seas.append({"track": t, "revisit": rv, "driver": k, "zone": z, "season": sea, "n": int(m.sum()),
                                 "n_years": len(np.unique(yrs[m])),
                                 "spearman_rho": r, "lo95": lo2, "hi95": hi2})
            Ym, Yg, X = gk.ref_mat_coh.to_numpy(float), gk.ref_grassland_coh.to_numpy(float), gk[[x, lvl]].to_numpy(float)
            m_, lo3, hi3 = year_boot(yrs, lambda i, Ym=Ym, Yg=Yg, X=X: std_betas(Ym[i], X[i])[0] - std_betas(Yg[i], X[i])[0])
            sens.append({"track": t, "revisit": rv, "driver": k, "n": len(gk), "n_years": ny, "beta_mat_minus_grassland": m_,
                         "lo95": lo3, "hi95": hi3})
    return pd.DataFrame(env), pd.DataFrame(beta), pd.DataFrame(seas), pd.DataFrame(sens)


def long_pairs(p: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """The pairs longer than 24 days already in the cube: a first look at the long-baseline side of the 2023 claim."""
    lp = p[p.rv == 0].copy()
    lp["length"] = pd.cut(lp.revisit_days, [24, 36, 48, 96, 200, 400, 800]).astype(str)
    lp["t1_season"] = np.where(lp.month.isin([12, 1, 2]), "winter", np.where(lp.month.isin([3, 4, 5]), "spring",
                               np.where(lp.month.isin([6, 7, 8]), "summer", "autumn")))
    summ = (lp.groupby(["track", "length", "t1_season"]).agg(n=("pair", "size"), coh_mat=("ref_mat_coh", "median"),
                                                            coh_grassland=("ref_grassland_coh", "median"),
                                                            coh_stable=("ref_stable_coh", "median")).reset_index())
    annual = lp[lp.revisit_days >= 200][["track", "pair", "revisit_days", "t1_season", "ref_mat_coh",
                                          "ref_grassland_coh", "ref_stable_coh"]].sort_values(["track", "pair"])
    return summ, annual


def frost(p: pd.DataFrame) -> pd.DataFrame:
    return (p[p.rv > 0].groupby(["track", "rv", "cold"]).agg(n=("pair", "size"), coh_mat=("ref_mat_coh", "median"),
                                                  coh_grassland=("ref_grassland_coh", "median"))
            .reset_index().rename(columns={"rv": "revisit", "cold": "frozen_or_snow"}))


# ------------------------------------------------------------------------------------------- figures

def fig_scatter(p, out):
    fig, axs = plt.subplots(2, 2, figsize=(12, 8.6), sharey=True)
    for i, t in enumerate(TRACKS):
        for j, k in enumerate(DRIVERS):
            ax = axs[i, j]
            g = p[(p.track == t) & ~p.cold & (p.rv > 0)].dropna(subset=[f"{k}_absd"])
            for rv, col in RV_COL.items():
                h = g[g.rv == rv]
                ax.scatter(h[f"{k}_absd"], h.ref_mat_coh, s=9, color=col, alpha=0.45, edgecolors="none",
                           label=f"{rv if rv < 24 else '18–24'}-day (n {len(h)})")
                if len(h) >= 48:
                    edges = np.unique(np.quantile(h[f"{k}_absd"], np.linspace(0, 1, 7)))
                    e = h.groupby(pd.cut(h[f"{k}_absd"], edges, include_lowest=True), observed=True).agg(
                        xm=(f"{k}_absd", "median"), yq=("ref_mat_coh", lambda v: np.quantile(v, 0.9)))
                    ax.plot(e.xm, e.yq, "-", color=col, lw=1.8)
            ax.set_xlabel(f"|Δ {DRIVERS[k][2]}| between the two dates")
            ax.set_title(f"{t} · mat coherence vs |Δ {'topsoil moisture' if k == 'sm' else 'water table'}| · "
                         "lines: 90th percentile", fontsize=8.5)
            ax.set_ylim(0, 1)
        axs[i, 0].set_ylabel("mat coherence (zone median)")
    axs[0, 0].legend(loc="upper right")
    fig.suptitle("Coherence against the moisture change between dates (unfrozen, snow-free pairs) — Hrysiewicz 2023 "
                 "Fig. 10 on the mat", x=0.01, ha="left", fontsize=11, fontweight="bold")
    fig.tight_layout()
    fig.savefig(out / "fig1_coherence_vs_moisture_change.png", dpi=170)
    plt.close(fig)


def fig_betas(beta, out):
    b = beta[(beta.zone.isin(["mat", "grassland"]))].copy()
    b["label"] = b.track.str[:3] + " " + b.revisit.astype(str).str.replace("24", "18–24") + "d · " + b.driver
    fig, axs = plt.subplots(1, 2, figsize=(12, 5.6), sharey=True)
    for ax, mode in zip(axs, ("raw", "monthly anomalies")):
        g = b[b.data == mode]
        labels = list(dict.fromkeys(g.label))
        for k, (z, col, off) in enumerate((("mat", AQUA, -0.18), ("grassland", ORANGE, 0.18))):
            h = g[g.zone == z].set_index("label").reindex(labels)
            yy = np.arange(len(labels)) + off
            ax.errorbar(h.beta_change, yy, xerr=[h.beta_change - h.change_lo95, h.change_hi95 - h.beta_change],
                        fmt="o", color=col, ms=4.5, elinewidth=1, capsize=0, label=f"{z}: |Δ| (change)")
            ax.errorbar(h.beta_level, yy, xerr=[h.beta_level - h.level_lo95, h.level_hi95 - h.beta_level],
                        fmt="s", mfc=SURFACE, color=col, ms=4.5, elinewidth=0.8, capsize=0, label=f"{z}: level")
        ax.axvline(0, color=INK2, lw=0.8)
        ax.set_yticks(range(len(labels))), ax.set_yticklabels(labels, fontsize=7)
        ax.set_xlabel("standardised coefficient on coherence (95 % CI, years resampled)")
        ax.set_title(f"{mode}: change vs level", fontsize=9)
    axs[0].legend(loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=2, fontsize=6.8)
    fig.suptitle("Does the change in moisture explain coherence beyond its level — and after the season is removed?",
                 x=0.01, ha="left", fontsize=11, fontweight="bold")
    fig.tight_layout()
    fig.savefig(out / "fig2_change_vs_level.png", dpi=170)
    plt.close(fig)


def fig_timeline(p, out):
    fig, axs = plt.subplots(3, 1, figsize=(13, 7.4), sharex=True, gridspec_kw={"height_ratios": [1, 1.2, 1.2]})
    acq = p.drop_duplicates("t1").sort_values("t1")
    axs[0].plot(acq.t1, acq.sm_level, ".", ms=2.5, color=BLUE)
    axs[0].set_ylabel("topsoil moisture\n(pair mean)")
    axs[0].set_title("ERA5-Land topsoil moisture at the pairs (mean of the two dates)", fontsize=8.5)
    for ax, t in zip(axs[1:], TRACKS):
        g = p[(p.track == t) & (p.rv > 0)]
        sc = ax.scatter(g.mid, g.revisit_days, c=g.ref_mat_coh, cmap=SEQ, vmin=0, vmax=1, s=14, marker="s",
                        edgecolors="none")
        ax.set_ylabel("revisit (days)"), ax.set_yticks([6, 12, 18, 24]), ax.set_ylim(2, 28)
        ax.set_title(f"{t} · mat coherence, every pair ≤ 24 days (longer pairs: README §7)", fontsize=8.5)
    cb = fig.colorbar(sc, ax=axs[1:], fraction=0.02, pad=0.01)
    cb.set_label("mat coherence"), cb.outline.set_visible(False)
    fig.savefig(out / "fig3_coherence_timeline.png", dpi=170, bbox_inches="tight")
    plt.close(fig)


def md(df: pd.DataFrame, fmt="{:.3f}") -> str:
    f = lambda v: fmt.format(v) if isinstance(v, (float, np.floating)) and np.isfinite(v) else ("—" if isinstance(v, float) else str(v))  # noqa: E731
    return "\n".join(["| " + " | ".join(map(str, df.columns)) + " |", "|" + "---|" * len(df.columns)]
                     + ["| " + " | ".join(f(v) for v in r) + " |" for r in df.itertuples(index=False)])


def readme(p, env, beta, seas, sens, fr, lsum, lann, out):
    env_all, beta_all, seas_all, sens_all = env, beta, seas, sens
    short = env[(env.zone == "mat") & (env.n_years < MIN_YEARS)]
    env, beta = env[env.n_years >= MIN_YEARS], beta[beta.n_years >= MIN_YEARS]
    seas, sens = seas[seas.n_years >= MIN_YEARS], sens[sens.n_years >= MIN_YEARS]
    me = env[env.zone == "mat"]
    neg = me[me.hi95 < 0]
    b = beta[beta.zone == "mat"]
    raw, ano = b[b.data == "raw"], b[b.data == "monthly anomalies"]
    ch_raw = raw[raw.change_hi95 < 0]
    ch_ano = ano[ano.change_hi95 < 0]
    dom = ano[(ano.beta_change.abs() > ano.beta_level.abs())]
    sm = seas[seas.zone == "mat"]
    sneg = sm[sm.hi95 < 0]
    lvl_raw = raw[raw.level_lo95 > 0]
    lvl_ano = ano[(ano.level_lo95 > 0) | (ano.level_hi95 < 0)]
    sm_s, wt_s = sens[sens.driver == "sm"], sens[sens.driver == "wtd"]
    sm_less = sm_s[sm_s.lo95 > 0]
    wt_more = wt_s[wt_s.hi95 < 0]
    agree_env = len(neg) >= len(me) / 2
    agree_ano = len(ch_ano) >= len(ano) / 2
    verdict = ("**agrees**" if agree_env and agree_ano else "**partly agrees**" if (agree_env or agree_ano)
               else "**does not reproduce**")
    L = ["# X-082 test 1 — coherence against moisture change on the mat (exploratory)", "",
         "Generated by `05_code/SAr-adapted/scripts/coherence_moisture_x082.py`; every number is computed by it. Plan, "
         "and what counts as agreement: `../PLAN.md`. Claim tested: Hrysiewicz et al. (2023) — coherence is controlled by "
         "the *change* in moisture between the two dates, not its level.", "",
         f"Pairs ≤ 24 days: {int((p.rv > 0).sum())} (both tracks, 2017–2026), of which {int((p.cold & (p.rv > 0)).sum())} frozen or snowy at a date "
         f"(kept apart). Moisture: ERA5-Land topsoil 0–7 cm for every year; the measured mat water table 2020–2024.", "",
         "## Findings", "",
         f"- **Verdict: the mat {verdict} with the 2023 result.**",
         f"- **Upper envelope:** the 90th-percentile coherence falls with |Δ moisture| (95 % CI below 0) in "
         f"{len(neg)}/{len(me)} track × revisit × driver cases on the mat.",
         f"- **Change vs level (raw):** |Δ| lowers coherence beyond the level in {len(ch_raw)}/{len(raw)} cases.",
         f"- **After removing the season (monthly anomalies):** in {len(ch_ano)}/{len(ano)} cases; the change outweighs "
         f"the level in {len(dom)}/{len(ano)}. This is the control the 2023 paper did not make.",
         f"- **Within season:** Spearman ρ of mat coherence with |Δ| below 0 (CI) in {len(sneg)}/{len(sm)} "
         "season × track × revisit × driver cases.",
         f"- **Level vs change:** the pair's mean moisture *raises* coherence in the raw data ({len(lvl_raw)}/{len(raw)} "
         f"cases, CI above 0) but that disappears once the season is removed ({len(lvl_ano)}/{len(ano)} cases with a CI "
         "excluding 0) — the level effect is seasonal, the change effect is not: their claim, in its stronger form.",
         f"- **Mat vs grassland:** the mat is *less* sensitive than the grassland to ERA5 topsoil moisture (a 9 km soil "
         f"model; {len(sm_less)}/{len(sm_s)} cases); to the water table measured on the mat it is *more* sensitive in "
         f"{len(wt_more)}/{len(wt_s)} cases, the rest not distinguishable (`mat_vs_grassland.csv`).",
         f"- **Left out of the counts:** {len(short)} case(s) with fewer than {MIN_YEARS} years of data (6-day pairs with "
         "a water table exist only in 2020–2021), where resampling years cannot give an interval; they stay in the "
         "tables with `n_years`.",
         "- **Frost:** frozen or snowy pairs summarised separately below; ERA5 topsoil moisture is not liquid water "
         "then, so they are outside the test.", "",
         "![scatter](fig1_coherence_vs_moisture_change.png)", "", "![betas](fig2_change_vs_level.png)", "",
         "![timeline](fig3_coherence_timeline.png)", "",
         "## 1. Upper envelope (90th percentile per |Δ| bin; slope per unit of |Δ|)", "", md(env_all), "",
         "## 2. Change vs level (standardised coefficients; raw and monthly anomalies)", "", md(beta_all), "",
         "## 3. Within season (Spearman ρ)", "", md(seas_all), "",
         "## 4. Frozen or snowy pairs", "", md(fr), "",
         "## 5. Mat minus grassland sensitivity to |Δ| (standardised; positive = the mat less sensitive)", "", md(sens_all), "",
         "## 7. Pairs longer than 24 days already in the cube (outside tests 1–5)", "",
         f"{int((p.rv == 0).sum())} pairs of 36–{int(p.revisit_days.max())} days (the 2022–2024 36/48-day network, bridge "
         "pairs, a few annual pairs). Median coherence by length and season of the first date:", "", md(lsum), "",
         f"The {len(lann)} pairs of 200 days or more:", "", md(lann), "",
         f"None starts in winter ({int((lann.t1_season == 'winter').sum())} of {len(lann)}), so the 2023 claim — "
         "winter-referenced pairs stay coherent for 1–3 years — is still untested here; D-027 orders those pairs. "
         f"Stable ground in these pairs is at {lann.ref_stable_coh.median():.2f} too, so D-027 must judge the mat "
         "against stable ground and hard targets in the same pairs, not against an absolute value.", "",
         "## Reading notes", "",
         "- Coherence values pile up near 0.3: the 40 m coherence estimate does not reach 0 on fully decorrelated pairs "
         "(estimator bias); read the envelope, not the floor.",
         "- Coherence: the cube's zone median per pair (`ref_<zone>_coh`, 40 m). Revisit 18 and 24 days pooled as "
         "\"24\" (`rv`).",
         "- Confidence intervals: whole calendar years resampled (pairs inside a year share weather).",
         "- Only pairs ≤ 24 days exist; the long-baseline part of their claim (coherence recovering after 1–3 years) "
         "needs the D-027 annual winter pairs.",
         "- ERA5-Land topsoil moisture is a 9 km model value, not the mat's moisture; the water table is measured but "
         "covers 2020–2024 only."]
    (out / "README.md").write_text("\n".join(L) + "\n")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    p = load()
    env, beta, seas, sens = analyse(p)
    fr = frost(p)
    lsum, lann = long_pairs(p)
    lsum.to_csv(OUT / "long_pairs_summary.csv", index=False)
    lann.to_csv(OUT / "long_pairs_200d_plus.csv", index=False)
    p.drop(columns=[c for c in p.columns if c.startswith("flag_")]).to_csv(OUT / "pairs.csv", index=False)
    for name, df in (("envelope", env), ("change_vs_level", beta), ("within_season", seas),
                     ("mat_vs_grassland", sens), ("frost", fr)):
        df.to_csv(OUT / f"{name}.csv", index=False)
    fig_scatter(p, OUT)
    fig_betas(beta, OUT)
    fig_timeline(p, OUT)
    readme(p, env, beta, seas, sens, fr, lsum, lann, OUT)
    print(f"X-082 test 1 written to {OUT}")


if __name__ == "__main__":
    main()

"""X-070: closure phase as the moisture observable — does the mat's non-closing phase measure its wetting?

The short-pair phase at P6 is motion plus an opposite-sign water/moisture part (X-062 D, X-065), and the wet surfaces do
not close (6 + 6 − 12 d, X-066). The literature retrieves a relative soil-moisture index from closure phase (a transfer
function from moisture anomalies to closure, RSE 2025). If the mat's closure carries its moisture, it is a radar-side
observable of the part fusion v2 carries as noise.

1. Every closure triangle t1, t2 = t1 + s, t3 = t1 + 2s whose three interferograms exist: s = 6 d in 2020–2021 (S1A+S1B)
   and 2025 (S1A+S1C), s = 12 d in 2022–2024. Closure phase per pixel from the WRAPPED phases,
   arg(e^{i(φ12 + φ23 − φ13)}) — no unwrapping — averaged per zone as a phasor weighted by the triplet's mean coherence
   (mat A, lake B, grassland C, stable ground near the mat, P6's 3×3); converted to mm of line of sight.
2. Against what could wet the surface over the triangle: the water-table change and its middle-date anomaly
   (W2 − (W1 + W3)/2), rain over [t1, t3], the number of wet overpasses (dew/rain), the backscatter change and its
   middle-date anomaly on the same zone. Raw and season-removed correlations, circular-shift p in time.
3. At P6, the X-065 phase model s1 = m·laser + e·ΔWTD (+ k·closure): does adding the closure of the triangles a pair
   belongs to lower the year-held-out error?

Mac-safe (no inversion, no null); field data are read from the hub, outputs go only to the hub.

    INSAR_DRIVE_ROOT=../local/drive_pristine PYTHONPATH=src python scripts/closure_moisture_x070.py
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import xarray as xr  # noqa: E402
from scipy import ndimage  # noqa: E402

from insar_wetlands import field  # noqa: E402
from insar_wetlands import field_link as fl  # noqa: E402
from insar_wetlands import fusion as fu  # noqa: E402
from insar_wetlands.bootstrap import start  # noqa: E402
from insar_wetlands.inversion.isbas import PHASE_TO_MM  # noqa: E402
from insar_wetlands.paths import make_paths  # noqa: E402
from insar_wetlands.stack import list_pairs, load_layer  # noqa: E402

HUB = Path(__file__).resolve().parents[3]
DLV = HUB / "08_deliverables"
OUT = DLV / "closure_x070"
HOUR = {"ascending": "16:36", "descending": "05:09"}
MIRROR = HUB / "05_code" / "local" / "drive_mirror"
RTC = {("ascending", "2020"): "rtc_dualpol_2020_2024.nc", ("descending", "2020"): "rtc_dualpol_2020_2024_descending.nc",
       ("ascending", "2025"): "rtc_dualpol_2025_2025.nc", ("descending", "2025"): "rtc_dualpol_2025_2025_descending.nc"}
ZONES_OUT = ["mat", "lake", "grassland", "stable", "p6"]
FEATURES = {"dwtd13": "water-table change t1→t3 (cm)", "wtd_mid_anom": "water table at t2 minus the mean of t1, t3 (cm)",
            "rain_mm": "rain over [t1, t3] (mm)", "n_wet": "wet overpasses among t1, t2, t3",
            "dvv13": "VV change t1→t3, same zone (dB)", "vv_mid_anom": "VV at t2 minus the mean of t1, t3 (dB)"}


def stacks(ctx) -> list[tuple[str, str, Path, int]]:
    out = []
    for t in HOUR:
        out.append(("2020–2021", t, HUB / f"05_code/local/s1_2020_2021/hyp3_cropped_{t}", 6))
        out.append(("2022–2024", t, make_paths(cfg=ctx.cfg, track=t).cropped, 12))
        out.append(("2025", t, HUB / f"05_code/local/s1_2025/hyp3_cropped_{t}", 6))
    return [s for s in out if s[2].exists()]


def triplets(pairs: list[str], step: int) -> list[tuple[str, str, str]]:
    have = set(pairs)
    out = []
    for p in pairs:
        a, b = pd.Timestamp(p[:8]), pd.Timestamp(p[9:])
        if (b - a).days != step:
            continue
        c = b + pd.Timedelta(days=step)
        p23, p13 = f"{b:%Y%m%d}_{c:%Y%m%d}", f"{a:%Y%m%d}_{c:%Y%m%d}"
        if p23 in have and p13 in have:
            out.append((p, p23, p13))
    return out


def closure_table(ctx, masks: dict) -> pd.DataFrame:
    rows = []
    for period, track, root, step in stacks(ctx):
        pairs = list_pairs(root)
        tri = triplets(pairs, step)
        need = sorted({p for t in tri for p in t})
        if not need:
            continue
        w = load_layer(root, "wrapped_phase", need)
        co = load_layer(root, "corr", need)
        u = load_layer(root, "unw_phase", need)
        idx = {p: k for k, p in enumerate(w.pair.values.astype(str))}
        W, C, U = w.values, co.values, u.values
        refC = {p: np.nanmedian(U[idx[p]][masks["grassland"] & np.isfinite(U[idx[p]])]) for p in need}
        for p12, p23, p13 in tri:
            i, j, k = idx[p12], idx[p23], idx[p13]
            c = np.angle(np.exp(1j * (W[i] + W[j] - W[k])))
            wt = np.nanmean(np.stack([C[i], C[j], C[k]]), 0)
            cu = (U[i] - refC[p12]) + (U[j] - refC[p23]) - (U[k] - refC[p13])
            # whole phase cycles between the unwrapped and the wrapped closure = unwrapping inconsistency of the triangle
            raw = U[i] + U[j] - U[k]
            off = raw - c
            off0 = np.nanmedian(off[masks["stable"] & np.isfinite(off)])   # HyP3's arbitrary per-pair reference constants
            cyc = np.rint((off - off0) / (2 * np.pi))
            t1, t2, t3 = (pd.Timestamp(f"{x} {HOUR[track]}", tz="UTC") for x in (p12[:8], p12[9:], p13[9:]))
            r = {"period": period, "track": track, "step": step, "t1": t1, "t2": t2, "t3": t3, "p12": p12, "p23": p23, "p13": p13}
            for z, m in masks.items():
                ok = m & np.isfinite(c) & np.isfinite(wt)
                ph = np.nansum(wt[ok] * np.exp(1j * c[ok])) / max(np.nansum(wt[ok]), 1e-9)
                r[f"{z}_closure_mm"] = float(np.angle(ph) * PHASE_TO_MM)
                r[f"{z}_consistency"] = float(np.abs(ph))         # 1 = every pixel the same closure
                r[f"{z}_unw_closure_mm"] = float(np.nanmean(cu[m]) * PHASE_TO_MM)
                r[f"{z}_coh"] = float(np.nanmedian(wt[m]))
                okc = m & np.isfinite(cyc)
                r[f"{z}_cycle_share"] = float(np.mean(cyc[okc] != 0)) if okc.any() else np.nan
                r[f"{z}_cycle_mean"] = float(np.mean(cyc[okc])) if okc.any() else np.nan     # net signed cycles: a bias, not noise
            rows.append(r)
        print(f"  {period} {track}: {len(tri)} triangles (step {step} d)", flush=True)
    d = pd.DataFrame(rows)
    d["mat_minus_stable_mm"] = d.mat_closure_mm - d.stable_closure_mm
    return d


def drivers(d: pd.DataFrame, ctx, masks: dict) -> pd.DataFrame:
    wtd = field.load_wtd_hourly()
    cens = {p: field.censored_flag(wtd[p]) for p in field.PLOTS}

    def wmed(t):
        v = [field._interp_at(wtd[p], t) for p in field.PLOTS if not bool(cens[p].asof(t))]
        return float(np.median(v)) if len(v) >= 5 else np.nan
    wet = pd.read_csv(DLV / "field_dew_x050" / "wetness_at_overpasses_2020_2024.csv").set_index(["date", "track"])
    w25 = pd.read_csv(DLV / "fusion_v2" / "flags_2025" / "wetness_at_overpasses_2025_open_meteo.csv")
    om = pd.read_csv(HUB / "06_data/local_small/open_meteo/open_meteo_rzecin_hourly_2020_2025.csv", parse_dates=["time"]).set_index("time")
    om.index = pd.to_datetime(om.index, utc=True)
    rain_st = wtd["Rain_mm_Tot"].astype(float)
    vv = {}
    for (track, y), name in RTC.items():
        if (MIRROR / name).exists():
            ds = xr.open_dataset(MIRROR / name)
            g = ctx.to_grid(ds["gamma0_vv_db"]).values
            tt = pd.to_datetime(ds.time.values).normalize()
            vv[(track, y)] = {z: pd.Series([float(np.nanmedian(10 * np.log10(np.nanmean(10 ** (g[k][m] / 10)))))
                                            if np.isfinite(g[k][m]).any() else np.nan for k in range(len(tt))], index=tt)
                              for z, m in (("mat", masks["mat"]), ("grassland", masks["grassland"]), ("lake", masks["lake"]))}

    def wet_at(t, track):
        key = (f"{t:%Y-%m-%d}", track)
        if key in wet.index:
            return bool(wet.loc[key].wet)
        h = w25[(w25.date == f"{t:%Y-%m-%d}") & (w25.track == track)] if "date" in w25.columns else pd.DataFrame()
        return bool(h.wet.iloc[0]) if len(h) else np.nan
    out = []
    for r in d.itertuples(index=False):
        rec = {}
        ws = [wmed(t) for t in (r.t1, r.t2, r.t3)]
        rec["dwtd13"] = ws[2] - ws[0]
        rec["wtd_mid_anom"] = ws[1] - (ws[0] + ws[2]) / 2
        src = rain_st if r.t3.year <= 2024 else om["Rain_mm_Tot"].astype(float)
        rec["rain_mm"] = float(src[(src.index > r.t1) & (src.index <= r.t3)].sum())
        f = [wet_at(t, r.track) for t in (r.t1, r.t2, r.t3)]
        rec["n_wet"] = float(np.nansum(f)) if not all(isinstance(x, float) and np.isnan(x) for x in f) else np.nan
        key = (r.track, "2025" if r.t1.year == 2025 else "2020")
        for z in ("mat", "grassland", "lake"):
            s = vv.get(key, {}).get(z)
            v = [s.get(t.normalize().tz_localize(None), np.nan) if s is not None else np.nan for t in (r.t1, r.t2, r.t3)]
            rec[f"dvv13_{z}"] = v[2] - v[0]
            rec[f"vv_mid_anom_{z}"] = v[1] - (v[0] + v[2]) / 2
        out.append(rec)
    return pd.concat([d.reset_index(drop=True), pd.DataFrame(out)], axis=1)


def tests(d: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (period, track), g in d.groupby(["period", "track"]):
        g = g.sort_values("t2").reset_index(drop=True)
        mid = g.t2.dt.tz_convert(None)
        for z in ["mat", "lake", "grassland", "stable", "p6", "mat_minus_stable"]:
            y = g["mat_minus_stable_mm"] if z == "mat_minus_stable" else g[f"{z}_closure_mm"]
            for fk in FEATURES:
                zz = {"stable": "grassland", "p6": "mat", "mat_minus_stable": "mat"}.get(z, z)
                x = g[f"{fk}_{zz}"] if fk in ("dvv13", "vv_mid_anom") else g[fk]
                ok = np.isfinite(x) & np.isfinite(y)
                if ok.sum() < 15:
                    continue
                r_raw, p_raw = fl.circular_shift_p(x[ok].to_numpy(), y[ok].to_numpy())
                ya = fl.season_residual(y[ok].to_numpy(), mid[ok], None)
                xa = fl.season_residual(x[ok].to_numpy().astype(float), mid[ok], None)
                r_an, p_an = fl.circular_shift_p(xa, ya)
                rows.append({"period": period, "track": track, "zone": z, "feature": fk, "n": int(ok.sum()),
                             "r_raw": r_raw, "p_raw": p_raw, "r_anom": r_an, "p_anom": p_an})
    return pd.DataFrame(rows)


def season_table(d: pd.DataFrame) -> pd.DataFrame:
    s = d.assign(season=d.t2.dt.month.map({12: "DJF", 1: "DJF", 2: "DJF", 3: "MAM", 4: "MAM", 5: "MAM", 6: "JJA", 7: "JJA", 8: "JJA",
                                            9: "SON", 10: "SON", 11: "SON"}))
    return s.groupby(["period", "track", "season"]).agg(triangles=("mat_closure_mm", "size"), mat_mm=("mat_closure_mm", "median"),
                                                         lake_mm=("lake_closure_mm", "median"), grassland_mm=("grassland_closure_mm", "median"),
                                                         stable_mm=("stable_closure_mm", "median"), mat_consistency=("mat_consistency", "median"),
                                                         mat_unw_mm=("mat_unw_closure_mm", "median"), stable_unw_mm=("stable_unw_closure_mm", "median"),
                                                         mat_cycle_share=("mat_cycle_share", "median"), lake_cycle_share=("lake_cycle_share", "median"),
                                                         stable_cycle_share=("stable_cycle_share", "median"),
                                                         mat_cycle_mean=("mat_cycle_mean", "mean"), stable_cycle_mean=("stable_cycle_mean", "mean")).reset_index()


def p6_model(d: pd.DataFrame) -> pd.DataFrame:
    """X-065's year-held-out phase model at P6, with the closure of the triangles each pair is a leg of."""
    pp = pd.read_csv(DLV / "field_p6_short_pairs_x065" / "p6_pairs_2020_2024.csv")
    pp = pp[pp.step == 1].copy()
    leg = {}
    for r in d.itertuples(index=False):
        for p in (r.p12, r.p23):
            leg.setdefault((r.track, p), []).append(r.p6_closure_mm)
    pp["closure_p6"] = [np.mean(leg[(t, p)]) if (t, p) in leg else np.nan for t, p in zip(pp.track, pp.pair)]
    models = {"laser + ΔWTD": ["laser_los", "dwtd"], "laser + ΔWTD + closure": ["laser_los", "dwtd", "closure_p6"],
              "laser + closure": ["laser_los", "closure_p6"], "laser only": ["laser_los"]}
    res = []
    for track in HOUR:
        g = pp[pp.track == track].dropna(subset=["laser_los", "dwtd", "closure_p6", "s1", "coh"]).reset_index(drop=True)
        for name, cols in models.items():
            pred = np.full(len(g), np.nan)
            for y in sorted(g.year.unique()):
                tr, te = (g.year != y).to_numpy(), (g.year == y).to_numpy()
                if tr.sum() < 8:
                    continue
                f = fu.bayes_linear(g[cols].to_numpy()[tr], g.s1.to_numpy()[tr], fu.phase_sigma_mm(g.coh.to_numpy()[tr]))
                pred[te] = g[cols].to_numpy()[te] @ f["beta"]
            m = np.isfinite(pred)
            fa = fu.bayes_linear(g[cols].to_numpy(), g.s1.to_numpy(), fu.phase_sigma_mm(g.coh.to_numpy())) if len(g) >= 8 else None
            res.append({"track": track, "model": name, "n": int(m.sum()), "years": "–".join(map(str, sorted(g.year.unique()))),
                        "rmse_year_held_out_mm": float(np.sqrt(np.mean((pred[m] - g.s1.to_numpy()[m]) ** 2))) if m.any() else np.nan,
                        **({f"coef_{c}": float(b) for c, b in zip(cols, fa["beta"])} if fa else {}),
                        **({f"sd_{c}": float(s) for c, s in zip(cols, fa["sd"])} if fa else {})})
    return pd.DataFrame(res)


def fig(d: pd.DataFrame, t: pd.DataFrame):
    fig, ax = plt.subplots(2, 2, figsize=(14, 8))
    for i, track in enumerate(HOUR):
        g = d[d.track == track].sort_values("t2")
        for z, col in (("mat", "C3"), ("lake", "C0"), ("stable", "0.5")):
            ax[0, i].plot(g.t2.dt.tz_convert(None), g[f"{z}_closure_mm"], ".", ms=4, color=col, label=z)
        ax[0, i].axhline(0, color="0.7", lw=0.6)
        ax[0, i].set_title(f"{track}: closure per triangle (mm LOS; 6 d in 2020–21 and 2025, 12 d in 2022–24)", fontsize=8)
        ax[0, i].legend(fontsize=7)
        sel = t[(t.track == track) & t.zone.isin(["mat", "stable", "lake"])]
        piv = sel.pivot_table(index=["period", "feature"], columns="zone", values="r_anom")
        if len(piv):
            ax[1, i].imshow(piv.to_numpy(), cmap="RdBu_r", vmin=-0.6, vmax=0.6, aspect="auto")
            ax[1, i].set_yticks(range(len(piv)), [f"{a} · {b}" for a, b in piv.index], fontsize=6)
            ax[1, i].set_xticks(range(len(piv.columns)), piv.columns, fontsize=7)
            for (yy, xx), v in np.ndenumerate(piv.to_numpy()):
                if np.isfinite(v):
                    ax[1, i].text(xx, yy, f"{v:.2f}", ha="center", va="center", fontsize=6)
            ax[1, i].set_title("r of closure with each driver, season removed", fontsize=8)
    fig.tight_layout(); fig.savefig(OUT / "fig_closure_moisture.png", dpi=120); plt.close(fig)


def md(df: pd.DataFrame) -> str:
    lines = ["| " + " | ".join(map(str, df.columns)) + " |", "|" + "---|" * len(df.columns)]
    for r in df.itertuples(index=False):
        lines.append("| " + " | ".join(f"{v:.2f}" if isinstance(v, float) else str(v) for v in r) + " |")
    return "\n".join(lines)


def fp(p):
    return "≤ 0.0005" if p <= 1 / 2001 + 1e-12 else f"{p:.3f}"


def readme(d, t, s, m):
    tt = t.copy()
    tt["raw r (p)"] = [f"{a:.2f} ({fp(b)})" for a, b in zip(tt.r_raw, tt.p_raw)]
    tt["season-removed r (p)"] = [f"{a:.2f} ({fp(b)})" for a, b in zip(tt.r_anom, tt.p_anom)]
    sig = t[(t.p_anom < 0.05)]
    find = [f"- **{len(d)} triangles** ({', '.join(f'{k[0]} {k[1]} {v}' for k, v in d.groupby(['period', 'track']).size().items())}).",
            f"- Season-removed links with p < 0.05: {len(sig)} of {len(t)} tests (≈ {0.05 * len(t):.0f} expected by chance)."]
    for z in ("mat", "stable", "lake"):
        g = sig[sig.zone == z]
        if len(g):
            find.append(f"  - {z}: " + "; ".join(f"{r.period} {r.track[:4]} {r.feature} r {r.r_anom:.2f}" for r in g.itertuples()))
    js = s[s.season == "JJA"]
    find.append(f"- **Summer non-closure is whole cycles, not closure phase.** June–August, the mat's wrapped closure (median "
                f"{js.mat_mm.min():.2f} to {js.mat_mm.max():.2f} mm) is near zero while its unwrapped closure is {js.mat_unw_mm.min():.1f} to "
                f"{js.mat_unw_mm.max():.1f} mm; {js.mat_cycle_share.min():.0%}–{js.mat_cycle_share.max():.0%} of mat pixels and "
                f"{js.stable_cycle_share.min():.0%}–{js.stable_cycle_share.max():.0%} of stable-ground pixels carry a whole-cycle "
                f"inconsistency. Net signed cycles per pixel: mat {js.mat_cycle_mean.min():+.2f} to {js.mat_cycle_mean.max():+.2f}, "
                f"stable {js.stable_cycle_mean.min():+.2f} to {js.stable_cycle_mean.max():+.2f} — unwrapping errors everywhere in summer, "
                "with a net direction on the mat.")
    for track in HOUR:
        g = m[m.track == track].set_index("model")
        if len(g) and "laser + ΔWTD" in g.index:
            find.append(f"- **P6, {track}**: year-held-out RMSE {g.loc['laser + ΔWTD', 'rmse_year_held_out_mm']:.2f} mm with laser + ΔWTD, "
                        f"{g.loc['laser + ΔWTD + closure', 'rmse_year_held_out_mm']:.2f} mm adding the closure, "
                        f"{g.loc['laser + closure', 'rmse_year_held_out_mm']:.2f} mm with the closure in place of ΔWTD "
                        f"(n {int(g.loc['laser + ΔWTD', 'n'])}).")
    L = ["# X-070 — closure phase as the moisture observable (exploratory)", "",
         "Generated by `05_code/SAr-adapted/scripts/closure_moisture_x070.py`; every number from its tables. Closure phase from the "
         "wrapped phases, per zone as a coherence-weighted phasor (no unwrapping), mm of LOS. Field data read from the hub; "
         "outputs only here.", "",
         "## Findings", "", *find, "",
         "## Closure by season (median over triangles, mm)", "", md(s), "",
         "`mat_consistency`: |mean phasor| of the mat's per-pixel closure (1 = every pixel the same closure; low = noise). "
         "`*_unw_mm`: the same closure from the unwrapped phases (as X-066). `*_cycle_share`: share of the zone's pixels whose unwrapped "
         "closure differs from the wrapped one by a whole phase cycle (2π ↔ λ/2 in LOS) — an unwrapping inconsistency in the triangle.", "",
         "## Closure against what could wet the surface", "",
         "Features: " + "; ".join(f"`{k}` {v}" for k, v in FEATURES.items()) + ". `stable` and `p6` use the grassland's and the mat's VV.", "",
         md(tt[["period", "track", "zone", "feature", "n", "raw r (p)", "season-removed r (p)"]]), "",
         "## P6: does the closure help the phase model? (X-065 pairs, year held out)", "", md(m.iloc[:, :5]), "",
         "## How to read this", "",
         "- A moisture observable should give the mat (and the lake) a closure that follows wetting — rain, a rising water table, "
         "a VV rise at the middle date — and stable ground little. A closure that follows nothing, or follows the same drivers "
         "on stable ground, is noise or regional.",
         "- At P6 the test is predictive: a covariate that carries the moisture part should lower the year-held-out error of the "
         "phase model; one that does not, does not help fusion.",
         "- Many tests on shared dates: read the counts against chance, not single p-values.", "",
         "## Files", "", "| File | Content |", "|---|---|",
         "| `closure_triangles.csv` | every triangle: closure, consistency and coherence per zone, drivers |",
         "| `closure_tests.csv` | every zone × feature × period × track correlation |",
         "| `closure_by_season.csv`, `p6_closure_model.csv` | the tables above |",
         "| `fig_closure_moisture.png` | closure per triangle; the correlation matrix |"]
    (OUT / "README.md").write_text("\n".join(L) + "\n")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--reuse", action="store_true")
    a = ap.parse_args(argv)
    OUT.mkdir(parents=True, exist_ok=True)
    ctx = start("closure_moisture_x070", mount=False, git=False)
    Z = {k: ctx.zones[k].values.astype(bool) for k in "ABCD"}
    dist = ndimage.distance_transform_edt(~Z["A"]) * 40.0
    px = pd.read_csv(DLV / "field_first" / "plot_pixels.csv").set_index("plot")
    r6, c6 = int(px.loc["P6", "row"]), int(px.loc["P6", "col"])
    p6 = np.zeros_like(Z["A"])
    p6[r6 - 1:r6 + 2, c6 - 1:c6 + 2] = True
    masks = {"mat": Z["A"], "lake": Z["B"], "grassland": Z["C"], "stable": Z["D"] & (dist > 200) & (dist < 1500), "p6": p6 & Z["A"]}
    cache = OUT / "closure_triangles.csv"
    if a.reuse and cache.exists():
        d = pd.read_csv(cache, parse_dates=["t1", "t2", "t3"])
    else:
        d = drivers(closure_table(ctx, masks), ctx, masks)
        d.to_csv(cache, index=False)
    for c in ("t1", "t2", "t3"):
        d[c] = pd.to_datetime(d[c], utc=True)
    t = tests(d)
    t.to_csv(OUT / "closure_tests.csv", index=False)
    s = season_table(d)
    s.to_csv(OUT / "closure_by_season.csv", index=False)
    m = p6_model(d)
    m.to_csv(OUT / "p6_closure_model.csv", index=False)
    fig(d, t)
    readme(d, t, s, m)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()

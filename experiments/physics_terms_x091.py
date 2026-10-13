"""X-091(a) — the gas and frost terms the buoyancy prior lacks: do they predict the P6 mat level better, year-blocked?

The prior that ties fusion and the segment-and-tie tests to physics is plain buoyancy on the water table (T47). Two
terms it never had: ebullition after atmospheric-pressure drops (terry2017: falls trigger sudden gas release, so a
floating mat should step down) and frost heave (kharanzhevskaya2026: freezing lifts the surface and breaks the
water-table relation). Each is added to the prior and kept only by the pre-registered rule
(``08_deliverables/insar_full_test/preregistration/X-091.yaml``): leave-one-year-out RMSE of the daily level change.
A superposed-epoch composite of the pressure-fall events shows the gas term's raw signal.

Reads the cube's hourly site table (laser, P6 water table, ERA5 pressure and soil temperature); writes only to the hub:
``08_deliverables/insar_full_test/physics_terms/``. Tried on the laser in hand (D-028).

    make run E=physics_terms_x091
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from insar_wetlands.cube import cube_root
from insar_wetlands.paths import hub_root
from insar_wetlands.scoring import assert_preregistered

HUB = hub_root()
CUBE = cube_root()
FT = HUB / "08_deliverables" / "insar_full_test"
OUT = FT / "physics_terms"
PREREG = FT / "preregistration" / "X-091.yaml"
YEARS = (2021, 2022, 2023, 2024)
RELEVEL = (pd.Timestamp("2022-07-05", tz="UTC"), pd.Timestamp("2022-08-08", tz="UTC"))
WINTER = {11, 12, 1, 2, 3, 4}
MIN_HOURS = 12
EVENT_Q, EVENT_GAP_D = 0.05, 5
N_PERM = 1000
KEEP_GAIN = 0.05


def load() -> pd.DataFrame:
    cols = ["time_utc", "laser_surface_cm", "laser_ok", "wtd_P6_cm", "om_era5_surface_pressure",
            "om_era5_land_soil_temperature_0_to_7cm"]
    s = pd.read_csv(CUBE / "silver" / "site_hourly.csv.gz", index_col=0, usecols=cols, low_memory=False)
    s.index = pd.to_datetime(s.index, format="ISO8601", utc=True)
    ok = s.laser_ok.map({True: True, "True": True}).fillna(False).astype(bool)
    s["laser_mm"] = (s.laser_surface_cm * 10.0).where(ok)
    s["p"] = s.om_era5_surface_pressure
    s["dp24"] = s.p - s.p.shift(24, freq="h").reindex(s.index)
    s["tsoil"] = s.om_era5_land_soil_temperature_0_to_7cm
    return s[(s.index.year >= min(YEARS)) & (s.index.year <= max(YEARS))]


def daily(s: pd.DataFrame) -> pd.DataFrame:
    g = s.resample("D")
    d = pd.DataFrame({"laser_mm": g.laser_mm.mean(), "n_laser": g.laser_mm.count(), "wtd_cm": g.wtd_P6_cm.mean(),
                      "gas": g.dp24.apply(lambda v: np.nansum(np.minimum(v.to_numpy(), 0.0)) / 24.0),
                      "tsoil": g.tsoil.mean()})
    d.loc[d.n_laser < MIN_HOURS, "laser_mm"] = np.nan
    d["frost"] = np.maximum(0.0, -d.tsoil).rolling(7, min_periods=7).sum()
    step = pd.DataFrame({"d_level": d.laser_mm.diff(), "d_wtd": d.wtd_cm.diff(), "gas": d.gas, "d_frost": d.frost.diff()},
                        index=d.index)
    span = (d.index - pd.Timedelta(days=1) < RELEVEL[1]) & (d.index > RELEVEL[0])
    step = step[~span].dropna()
    step["year"] = step.index.year
    step["winter"] = step.index.month.isin(list(WINTER))
    return step


def fit(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    return np.linalg.lstsq(x, y, rcond=None)[0]


def loyo(step: pd.DataFrame, cols: list[str], score_mask: np.ndarray) -> tuple[float, list[dict]]:
    """Pooled held-out RMSE (on ``score_mask`` days) and per-fold coefficients, leave one year out."""
    err, coefs = [], []
    for y in YEARS:
        tr, te = step.year != y, (step.year == y) & score_mask
        if te.sum() == 0 or tr.sum() < 10:
            continue
        b = fit(step.loc[tr, cols].to_numpy(), step.loc[tr, "d_level"].to_numpy())
        pred = step.loc[te, cols].to_numpy() @ b
        err.append(step.loc[te, "d_level"].to_numpy() - pred)
        coefs.append({"held_out": y, "days_scored": int(te.sum()), **dict(zip(cols, b))})
    e = np.concatenate(err) if err else np.array([np.nan])
    return float(np.sqrt(np.mean(e ** 2))), coefs


def events(s: pd.DataFrame) -> pd.DatetimeIndex:
    dp = s.dp24.dropna()
    low = dp <= dp.quantile(EVENT_Q)
    starts = dp.index[low & ~low.shift(1, fill_value=False)]
    out, last = [], None
    for t in starts:
        if last is None or t - last >= pd.Timedelta(days=EVENT_GAP_D):
            out.append(t)
            last = t
    return pd.DatetimeIndex(out)


def response(s: pd.DataFrame, ts, g: float) -> pd.DataFrame:
    """Residual (laser minus g·ΔWTD) change 0–24 h and 0–72 h after each time in ``ts``."""
    lev, w = s.laser_mm, s.wtd_P6_cm
    rows = []
    for t in pd.DatetimeIndex(ts):
        for h in (24, 72):
            a, b = lev.get(t), lev.get(t + pd.Timedelta(hours=h))
            wa, wb = w.get(t), w.get(t + pd.Timedelta(hours=h))
            if pd.notna(a) and pd.notna(b) and pd.notna(wa) and pd.notna(wb):
                rows.append({"t0": t, "hours": h, "laser_change_mm": b - a, "wtd_change_cm": wb - wa,
                             "residual_mm": (b - a) - g * (wb - wa)})
    return pd.DataFrame(rows, columns=["t0", "hours", "laser_change_mm", "wtd_change_cm", "residual_mm"])


def composite(s: pd.DataFrame, ev: pd.DatetimeIndex, g: float) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Responses after the events, and after control hours (noon, same months, dP24 above the median, far from any
    event)."""
    dp = s.dp24
    far = pd.Series(True, index=s.index)
    for t in ev:
        far[(s.index > t - pd.Timedelta(days=EVENT_GAP_D)) & (s.index < t + pd.Timedelta(days=EVENT_GAP_D))] = False
    ctrl = s.index[far & (dp > dp.median()) & s.index.month.isin(list(set(ev.month))) & (s.index.hour == 12)]
    return response(s, ev, g), response(s, ctrl, g)


def perm_p(ev_r: pd.DataFrame, s: pd.DataFrame, n_ev: int, g: float, rng: np.random.Generator) -> dict[int, float]:
    """Block permutation: event times redrawn uniformly within the same years; p of a mean residual at least as low."""
    out = {}
    valid = s.index[s.laser_mm.notna() & s.wtd_P6_cm.notna()]
    for h in (24, 72):
        obs = ev_r[ev_r.hours == h].residual_mm.mean()
        if not np.isfinite(obs):
            out[h] = np.nan
            continue
        null = []
        for _ in range(N_PERM):
            ts = pd.DatetimeIndex(rng.choice(valid, size=n_ev, replace=False))
            r = response(s, ts, g)
            r = r[r.hours == h]
            null.append(r.residual_mm.mean() if len(r) else np.nan)
        null = np.array(null)
        null = null[np.isfinite(null)]
        out[h] = float((1 + np.sum(null <= obs)) / (1 + len(null)))
    return out


def main() -> None:
    assert_preregistered(PREREG, [OUT / "verdict.csv"])
    OUT.mkdir(parents=True, exist_ok=True)
    s = load()
    step = daily(s)
    allm, winm = np.ones(len(step), bool), step.winter.to_numpy()
    base_all, cb_all = loyo(step, ["d_wtd"], allm)
    gas_all, cg = loyo(step, ["d_wtd", "gas"], allm)
    base_win, cb_win = loyo(step, ["d_wtd"], winm)
    frost_win, cf = loyo(step, ["d_wtd", "d_frost"], winm)
    verdict = pd.DataFrame([
        {"term": "gas (pressure falls)", "days": "all", "rmse_prior_mm": base_all, "rmse_prior_plus_term_mm": gas_all,
         "gain": 1 - gas_all / base_all, "sign_ok_every_fold": all(c["gas"] > 0 for c in cg),
         "folds": len(cg)},
        {"term": "frost (soil freezing degree-days)", "days": "winter (Nov–Apr)", "rmse_prior_mm": base_win,
         "rmse_prior_plus_term_mm": frost_win, "gain": 1 - frost_win / base_win,
         "sign_ok_every_fold": all(c["d_frost"] > 0 for c in cf), "folds": len(cf)},
    ])
    verdict["kept"] = (verdict.gain >= KEEP_GAIN) & verdict.sign_ok_every_fold
    coefs = pd.concat([pd.DataFrame(cb_all).assign(model="prior", days="all"),
                       pd.DataFrame(cg).assign(model="prior + gas", days="all"),
                       pd.DataFrame(cb_win).assign(model="prior", days="winter"),
                       pd.DataFrame(cf).assign(model="prior + frost", days="winter")], ignore_index=True)
    g_all = float(fit(step[["d_wtd"]].to_numpy(), step.d_level.to_numpy())[0]) / 1.0   # mm per cm (daily fit)
    ev = events(s)
    ev_r, ct_r = composite(s, ev, g_all)
    p = perm_p(ev_r, s, max(1, ev_r[ev_r.hours == 24].t0.nunique()), g_all, np.random.default_rng(91))
    comp = pd.concat([ev_r.assign(kind="event"), ct_r.assign(kind="control")], ignore_index=True)
    summ = comp.groupby(["kind", "hours"]).agg(n=("residual_mm", "size"), residual_mean_mm=("residual_mm", "mean"),
                                               laser_mean_mm=("laser_change_mm", "mean"),
                                               wtd_mean_cm=("wtd_change_cm", "mean")).reset_index()
    summ["perm_p_lower"] = summ.apply(lambda r: p.get(r.hours) if r.kind == "event" else np.nan, axis=1)
    verdict.to_csv(OUT / "verdict.csv", index=False)
    coefs.to_csv(OUT / "coefficients.csv", index=False)
    comp.to_csv(OUT / "pressure_events.csv", index=False)
    summ.to_csv(OUT / "pressure_composite.csv", index=False)
    (OUT / "README.md").write_text(readme(verdict, coefs, summ, len(step), int(step.winter.sum()), len(ev)),
                                   encoding="utf-8")
    print(f"wrote {OUT}")


def _md(df: pd.DataFrame, nd: int = 3) -> str:
    d = df.copy()
    for c in d.columns:
        if pd.api.types.is_float_dtype(d[c]):
            d[c] = d[c].map(lambda x: "—" if pd.isna(x) else f"{x:.{nd}f}")
    return "\n".join(["| " + " | ".join(map(str, d.columns)) + " |", "|" + "---|" * len(d.columns)] +
                     ["| " + " | ".join(map(str, r)) + " |" for r in d.itertuples(index=False)])


def readme(verdict, coefs, summ, n_days, n_winter, n_ev) -> str:
    return f"""# X-091(a) — gas and frost terms in the buoyancy prior (generated)

Generated by `experiments/physics_terms_x091.py` (`make run E=physics_terms_x091`); do not edit. Pre-registration:
`../preregistration/X-091.yaml`. P6 only, laser in hand (2021–2024, years with the water-table logger; D-028): tried on
one point and four years, to be solidified with the longer laser record and a second surface point (D-031).

Daily level changes used: {n_days} (winter: {n_winter}). Pressure-fall events (lowest {EVENT_Q:.0%} of 24-h changes, at
least {EVENT_GAP_D} days apart): {n_ev}.

## Verdict (kept if held-out RMSE falls by >= {KEEP_GAIN:.0%} and the sign is right in every fold)

{_md(verdict)}

## Coefficients per held-out year

{_md(coefs)}

## Pressure-fall composite: residual change after events vs controls (mm; p = share of {N_PERM} random event sets
at least as low, floor 1/{N_PERM + 1})

{_md(summ)}

Files: `verdict.csv`, `coefficients.csv`, `pressure_events.csv` (every event and control window), `pressure_composite.csv`.
"""


if __name__ == "__main__":
    main()

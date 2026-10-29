"""X-106 — freeze the mat method v4 and report it on the design years.

v4 (selected by X-119, T78): the chain of consecutive HyP3 pairs at P6, read at the half-cell-corrected position (T75,
T77), hard-target reference, vertical, each winter's level step set by SPAMS. This script (1) fits SPAMS once on all
valid 6-day laser changes of the design years, (2) writes the frozen parameters with the code sha to the hub
(``08_deliverables/insar_full_test/v4/v4_frozen_params.yaml`` — field-derived, never in this public repo), (3) reports v4
and the plain chain on the design years (in-sample, for reference). ``build_v4`` is what X-107 runs, once, on the held-
out years with the frozen file and no refit. Rule: ``preregistration/X-106.yaml``.

    make run E=freeze_v4_x106
"""
from __future__ import annotations

import subprocess
import sys
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import p6_validation_x082 as X  # noqa: E402  (consecutive-pair rule, laser readers)
import spams_baseline_x112 as SB  # noqa: E402  (the daily weather/laser table)

from insar_wetlands import spams, v4  # noqa: E402
from insar_wetlands.cube import cube_root  # noqa: E402
from insar_wetlands.cube.core import rad_to_los_mm  # noqa: E402
from insar_wetlands.paths import hub_root  # noqa: E402
from insar_wetlands.scoring import Protocol, assert_preregistered, chain, detrended_stats, level_views, winter_steps  # noqa: E402

HUB = hub_root()
FT = HUB / "08_deliverables" / "insar_full_test"
V4DIR = FT / "v4"
PARAMS = V4DIR / "v4_frozen_params.yaml"
OUT = V4DIR / "design_report"
PREREG = FT / "preregistration" / "X-106.yaml"
DESIGN = (2021, 2022, 2023, 2024, 2025)
IFG_SETS = ("2017_2019", "2020_2021", "2022_2024", "2025", "2026")
LAG = 6


# ------------------------------------------------------------------------------ the frozen parts

def fit_and_freeze() -> dict:
    d = SB.table()
    valid = (d.dh_mm.notna() & ~d.frozen_win & d.year.isin(DESIGN)).to_numpy()
    par = spams.fit(d.rain_era5_mm.to_numpy(), d.et0_mm.to_numpy(), d.dh_mm.to_numpy(), valid, LAG)
    shift = pd.read_csv(FT / "hyp3_placement" / "decisive" / "verdict.csv").set_index("track")
    px = pd.read_csv(HUB / "08_deliverables" / "field_first" / "plot_pixels.csv").set_index("plot")
    sha = subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True,
                         cwd=Path(__file__).resolve().parents[1]).stdout.strip()
    frozen = {"method": "v4 = chain of consecutive HyP3 pairs at P6 (corrected position), hard-target reference, "
                        "vertical, each winter step set by SPAMS (X-106, selected by X-119)",
              "frozen_on": str(date.today()), "code_sha": sha, "design_years": list(DESIGN),
              "read_position": {t: {"x": float(px.loc["P6", "E"]) - float(shift.loc[t, "median_dx_m"]),
                                    "y": float(px.loc["P6", "N"]) - float(shift.loc[t, "median_dy_m"])} for t in X.TRACKS},
              "pairs": {"consecutive": True, "exclude_platforms": ["S1D"], "reference": "ref_hard_unw"},
              "spams": {"x_p": par["x_p"], "x_e": par["x_e"], "tau_days": par["tau"], "fit_days": par["n"],
                        "inputs": "ERA5 daily precipitation, FAO-56 reference evapotranspiration (cube site_hourly)"},
              "winters": "protocol.yaml season_months / spring_window; same laser datum for scoring"}
    PARAMS.write_text(yaml.safe_dump(frozen, sort_keys=False, allow_unicode=True), encoding="utf-8")
    return frozen


def unw_at(track: str, x: float, y: float) -> dict[str, float]:
    out = {}
    for s in IFG_SETS:
        f = cube_root() / "silver" / f"ifg_{s}_{track}.nc"
        if not f.exists():
            continue
        with xr.open_dataset(f) as ds:
            ix = int(np.argmin(np.abs(ds.x.values - x)))
            iy = int(np.argmin(np.abs(ds.y.values - y)))
            sub = ds.unw.isel(x=slice(max(ix - 2, 0), ix + 3), y=slice(max(iy - 2, 0), iy + 3))
            v = sub.interp(x=x, y=y, method="linear").values
            k = pd.to_datetime(ds.t1.values).strftime("%Y%m%d") + "_" + pd.to_datetime(ds.t2.values).strftime("%Y%m%d")
            out.update({kk: float(vv) for kk, vv in zip(k, v)})
    return out


def spams_level(frozen: dict) -> pd.Series:
    d = SB.table()
    s = frozen["spams"]
    R = spams.reversible(d.rain_era5_mm.to_numpy(), d.et0_mm.to_numpy(), s["x_p"], s["x_e"], int(s["tau_days"]))
    return pd.Series(R, index=d.index)


def build_v4(track: str, frozen: dict, years: tuple[int, ...], relevel: bool = True) -> pd.DataFrame:
    """The v4 level series (vertical mm, uplift positive) of one track over ``years``: time, seg, est."""
    g = pd.read_csv(cube_root() / "gold" / "plots_pair.csv.gz", low_memory=False,
                    usecols=["track", "grid", "plot", "t1", "t2", "platforms", "ref_hard_unw", "los_u"])
    g = g[(g["plot"] == "P6") & (g.grid == "40m") & (g.track == track)].copy()
    for c in ("t1", "t2"):
        g[c] = pd.to_datetime(g[c], format="ISO8601", utc=True)
    g = g[X.consecutive(g) & g.t1.dt.year.isin(years) & g.t2.dt.year.isin(years)]
    g = g[~g.platforms.astype(str).str.contains("|".join(frozen["pairs"]["exclude_platforms"]))]
    pos = frozen["read_position"][track]
    u = unw_at(track, pos["x"], pos["y"])
    key = g.t1.dt.strftime("%Y%m%d") + "_" + g.t2.dt.strftime("%Y%m%d")
    g["y_vert_mm"] = -rad_to_los_mm(key.map(u).to_numpy() - g.ref_hard_unw.to_numpy()) / g.los_u.to_numpy()
    g = g.dropna(subset=["y_vert_mm"]).sort_values("t1")
    c = chain(g.assign(est_mm=g.y_vert_mm), "est_mm")
    if relevel:
        R = spams_level(frozen)
        at = lambda t: float(R.reindex([pd.Timestamp(t).tz_convert(None).normalize()]).iloc[0])
        c = v4.relevel(c, at, v4.winter_bounds(c.time))
    return c


# ------------------------------------------------------------------------------ the design-year report

def report(frozen: dict, prot: Protocol) -> tuple[pd.DataFrame, pd.DataFrame]:
    lev = X.laser()[1]
    rows, wins = [], []
    for track in X.TRACKS:
        for name, rl in (("v4", True), ("plain chain (corrected position)", False)):
            c = build_v4(track, frozen, DESIGN, relevel=rl)
            c["truth"] = lev.reindex(pd.DatetimeIndex(c.time).round("h")).to_numpy()
            for view, dd in level_views(c, prot):
                if view.startswith("season") and len(dd) >= 4:
                    rows.append({"track": track, "estimate": name, "view": view, **detrended_stats(dd.time, dd.truth, dd.est)})
            z = c.dropna(subset=["truth"]).copy()
            for _, cs in z.groupby("seg"):
                z.loc[cs.index, "truth"] -= cs.truth.iloc[0]
                z.loc[cs.index, "est"] -= cs.est.iloc[0]
            wins.append(winter_steps(z, prot).assign(track=track, estimate=name))
    return pd.DataFrame(rows), pd.concat(wins, ignore_index=True)


def attribution(frozen: dict, W: pd.DataFrame) -> pd.DataFrame:
    """Why v4's in-sample winter steps differ from X-119's A: per scored winter, SPAMS's step between the scored (laser)
    epochs with the frozen fit and with X-119's fit (without the winter's two years), and the part of v4's step the
    chain carries outside the re-levelled span (v4 re-levels between its own epochs, which cannot depend on the laser)."""
    d, R = SB.table(), spams_level(frozen)
    out = []
    for _, r in W[W.estimate == "v4"].iterrows():
        f, t, y = pd.Timestamp(r["from"]), pd.Timestamp(r["to"]), int(r.winter[:4])
        valid = (d.dh_mm.notna() & ~d.frozen_win & d.year.isin(DESIGN) & ~d.year.isin({y, y + 1})).to_numpy()
        P, E = d.rain_era5_mm.to_numpy(), d.et0_mm.to_numpy()
        par = spams.fit(P, E, d.dh_mm.to_numpy(), valid, LAG)
        Rx = pd.Series(spams.reversible(P, E, par["x_p"], par["x_e"], par["tau"]), index=d.index)
        out.append({"track": r.track, "winter": r.winter, "laser_change_mm": r.laser_change_mm,
                    "v4_change_mm": r.insar_change_mm, "spams_frozen_scored_mm": R[t] - R[f],
                    "spams_without_winter_years_mm": Rx[t] - Rx[f],
                    "chain_outside_relevel_mm": r.insar_change_mm - (R[t] - R[f])})
    return pd.DataFrame(out)


def main() -> None:
    assert_preregistered(PREREG, [OUT / "summary.csv"])
    OUT.mkdir(parents=True, exist_ok=True)
    prot = Protocol.load(FT / "protocol.yaml")
    frozen = fit_and_freeze()
    S, W = report(frozen, prot)
    summ = (S.groupby(["estimate", "track"]).rmse_detrended_mm.median().rename("season_rmse_median_mm").to_frame()
            .join(W.assign(ae=W.difference_mm.abs()).groupby(["estimate", "track"]).ae.mean().rename("winter_mae_mm"))
            .reset_index())
    A = attribution(frozen, W)
    for name, df in (("seasons", S), ("winters", W), ("summary", summ), ("attribution", A)):
        df.to_csv(OUT / f"{name}.csv", index=False)
    md = lambda df: "\n".join(["| " + " | ".join(map(str, df.columns)) + " |", "|" + "---|" * len(df.columns)] +
                              ["| " + " | ".join(f"{x:.2f}" if isinstance(x, float) else str(x) for x in r) + " |"
                               for r in df.itertuples(index=False)])
    (OUT / "README.md").write_text(f"""# v4 frozen — design-year report (generated, in-sample)

Generated by `experiments/freeze_v4_x106.py`; do not edit. Frozen parameters: `../v4_frozen_params.yaml` (code sha
{frozen['code_sha']}). In-sample on the design years — for reference only; the scored comparison is X-107 on the held-out
years. Vertical mm against the laser.

## Summary

{md(summ)}

## Winter steps

{md(W[['track', 'estimate', 'winter', 'laser_change_mm', 'insar_change_mm', 'difference_mm']])}

## Attribution against X-119's A

X-119's A re-levelled between the scored (laser) epochs with a SPAMS fit made without the winter's two years. The frozen
v4 re-levels between its own chain epochs (the protocol's winter, read on the chain: it cannot depend on the laser) with
one SPAMS fit on all design years. `chain_outside_relevel_mm` is what the chain adds between the scored and the
re-levelled epochs; the two SPAMS columns show how much the anchor's step depends on the years it is fitted on.

{md(A)}
""", encoding="utf-8")
    print(summ.to_string(index=False))
    print(f"frozen: {PARAMS}")


if __name__ == "__main__":
    main()

"""X-143 — do the skip-one pairs at P6 carry a phase that does not close, common to both tracks? (R4 M2)

Triplets of three consecutive six-day acquisitions (2021, 2025) at P6: closure c = wrap(w02 − w01 − w12) from the
wrapped phase alone. Motion, and any whole cycle of it, cancels in c exactly, so what is left is the non-closing phase
(and noise). The ascending closure is matched with the descending triplet 1.5 ± 0.6 days later at all three dates
(X-136's rule); a common non-closing term shows as a correlation of the two.

Pre-registration: hub ``08_deliverables/insar_full_test/preregistration/X-143.yaml`` (checked before the run).
Reads the cube and the hub; writes only to the hub: ``08_deliverables/insar_full_test/closure_x143/``. No field value
in this file.

    make run E=closure_x143
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import p6_validation_x082 as X  # noqa: E402  (PHASE_TO_MM, the acquisition list, the tracks)

from insar_wetlands.paths import cube_root, hub_root  # noqa: E402
from insar_wetlands.scoring import assert_preregistered  # noqa: E402

HUB, CUBE = hub_root(), cube_root()
FT = HUB / "08_deliverables" / "insar_full_test"
OUT = FT / "closure_x143"
PREREG = FT / "preregistration" / "X-143.yaml"
OFFSET_D, TOL_D = 1.5, 0.6
YEARS = (2021, 2025)
LASER_START = pd.Timestamp("2021-05-27", tz="UTC")
N_BOOT, SEED, BLOCK = 2000, 142, 3
COLS = ["track", "plot", "grid", "t1", "t2", "platforms", "wrapped_1x1", "ref_hard_wrapped", "laser_snow_t1",
        "laser_snow_t2", "d_laser_los_mm", "d_wtd_P6_cm", "d_om_era5_land_soil_moisture_0_to_7cm", "unw_1x1",
        "ref_hard_unw"]


def wrap(x):
    return np.angle(np.exp(1j * np.asarray(x, float)))


def truthy(s: pd.Series) -> pd.Series:
    return s.map({True: True, False: False, "True": True, "False": False}).fillna(False).astype(bool)


def load() -> pd.DataFrame:
    p = pd.read_csv(CUBE / "gold" / "plots_pair.csv.gz", usecols=COLS, low_memory=False)
    p = p[(p["plot"] == "P6") & (p.grid == "40m") & ~p.platforms.astype(str).str.contains("S1D")].copy()
    for c in ("t1", "t2"):
        p[c] = pd.to_datetime(p[c], format="ISO8601", utc=True)
    p["snow"] = truthy(p.laser_snow_t1) | truthy(p.laser_snow_t2)
    return p.drop_duplicates(["track", "t1", "t2"])


def triplets(p: pd.DataFrame) -> pd.DataFrame:
    """Every t0 < t1 < t2 of consecutive six-day acquisitions with its three pairs present, closure in mm."""
    acq = pd.read_csv(X.ASF, parse_dates=["date"])
    rows = []
    for track in X.TRACKS:
        g = p[p.track == track].set_index(["t1", "t2"])
        days = np.sort(acq[acq.track == track].date.dt.normalize().unique())
        times = pd.DatetimeIndex(sorted(set(g.index.get_level_values(0)) | set(g.index.get_level_values(1))))
        by_day = {t.tz_convert(None).normalize(): t for t in times}
        for i in range(len(days) - 2):
            d0, d1, d2 = (pd.Timestamp(days[i + j]) for j in range(3))
            if (d1 - d0).days != 6 or (d2 - d1).days != 6 or d0.year not in YEARS:
                continue
            t0, t1, t2 = by_day.get(d0), by_day.get(d1), by_day.get(d2)
            if t0 is None or t1 is None or t2 is None or t0 < LASER_START:
                continue
            try:
                a, b, s = g.loc[(t0, t1)], g.loc[(t1, t2)], g.loc[(t0, t2)]
            except KeyError:
                continue
            rows.append({
                "track": track, "t0": t0, "t1": t1, "t2": t2, "year": t0.year,
                "snow": bool(a.snow or b.snow or s.snow),
                "closure_mm": float(wrap(s.wrapped_1x1 - a.wrapped_1x1 - b.wrapped_1x1) * X.PHASE_TO_MM),
                "closure_ref_mm": float(wrap((s.wrapped_1x1 - s.ref_hard_wrapped) - (a.wrapped_1x1 - a.ref_hard_wrapped)
                                             - (b.wrapped_1x1 - b.ref_hard_wrapped)) * X.PHASE_TO_MM),
                "skip_residual_mm": float((s.unw_1x1 - s.ref_hard_unw) * X.PHASE_TO_MM - s.d_laser_los_mm),
                "d_laser_los_mm": s.d_laser_los_mm, "d_wtd_cm": s.d_wtd_P6_cm,
                "d_soil_moisture": s.d_om_era5_land_soil_moisture_0_to_7cm})
    return pd.DataFrame(rows)


def match(t: pd.DataFrame) -> pd.DataFrame:
    a = t[t.track == "ascending"].add_prefix("a_")
    d = t[t.track == "descending"].add_prefix("d_")
    m = a.merge(d, how="cross")
    for k in ("t0", "t1", "t2"):
        dt = (m[f"d_{k}"] - m[f"a_{k}"]).dt.total_seconds() / 86400
        m = m[(dt - OFFSET_D).abs().le(TOL_D)]
        if not len(m):
            break
    return m.sort_values("a_t0").drop_duplicates("a_t0").drop_duplicates("d_t0").reset_index(drop=True)


def block_boot(n: int, rng: np.random.Generator) -> np.ndarray:
    starts = rng.integers(0, max(n - BLOCK + 1, 1), size=int(np.ceil(n / BLOCK)))
    return np.concatenate([np.arange(s, min(s + BLOCK, n)) for s in starts])[:n]


def r_ci(x: np.ndarray, y: np.ndarray, rng: np.random.Generator) -> tuple[float, float, float, int]:
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if len(x) < 4:
        return np.nan, np.nan, np.nan, len(x)
    r = float(np.corrcoef(x, y)[0, 1])
    bs = []
    for _ in range(N_BOOT):
        i = block_boot(len(x), rng)
        if np.std(x[i]) > 0 and np.std(y[i]) > 0:
            bs.append(np.corrcoef(x[i], y[i])[0, 1])
    lo, hi = np.percentile(bs, [2.5, 97.5])
    return r, float(lo), float(hi), len(x)


def circ_mean_ci(c_mm: np.ndarray, rng: np.random.Generator) -> tuple[float, float, float]:
    ph = c_mm / X.PHASE_TO_MM
    f = lambda v: float(np.angle(np.mean(np.exp(1j * v))) * X.PHASE_TO_MM)  # noqa: E731
    bs = [f(ph[block_boot(len(ph), rng)]) for _ in range(N_BOOT)]
    lo, hi = np.percentile(bs, [2.5, 97.5])
    return f(ph), float(lo), float(hi)


def main() -> None:
    assert_preregistered(PREREG, [OUT / "verdict.json"])
    rng = np.random.default_rng(SEED)
    t = triplets(load())
    ok = t[~t.snow].sort_values(["track", "t0"]).reset_index(drop=True)
    counts = ok.groupby("track").size().reindex(list(X.TRACKS), fill_value=0)
    OUT.mkdir(parents=True, exist_ok=True)
    t.to_csv(OUT / "triplets.csv", index=False)
    if (counts < 10).any():
        v = {"verdict": "stopped: fewer than 10 snow-free triplets on a track", "counts": counts.to_dict()}
        (OUT / "verdict.json").write_text(json.dumps(v, indent=1), encoding="utf-8")
        print(v)
        return
    m = match(ok)
    m.to_csv(OUT / "matched.csv", index=False)
    r, lo, hi, n = r_ci(m.a_closure_mm.to_numpy(), m.d_closure_mm.to_numpy(), rng)
    if n < 20:
        verdict = "(c) no non-closing phase common to both tracks is detected; with fewer than 20 matched triplets the test is weak"
    elif lo > 0:
        verdict = "(a) the skip-one pairs carry a non-closing phase common to both tracks"
    else:
        verdict = "(b) no non-closing phase common to both tracks is detected in the skip-one pairs"
    rr, rlo, rhi, rn = r_ci(m.a_closure_ref_mm.to_numpy(), m.d_closure_ref_mm.to_numpy(), rng)
    drv, per = [], []
    for track, g in ok.groupby("track"):
        c = g.closure_mm.to_numpy()
        mu, mlo, mhi = circ_mean_ci(c, rng)
        res = g.skip_residual_mm.to_numpy()
        mres = np.isfinite(res)
        share = float(np.var(c[mres]) / np.var(res[mres])) if mres.sum() > 3 else np.nan
        per.append({"track": track, "triplets": len(g), "circular_mean_mm": mu, "ci_lo": mlo, "ci_hi": mhi,
                    "sd_mm": float(np.std(c, ddof=1)), "var_share_of_skip_residual": share,
                    "n_with_laser": int(mres.sum())})
        for name, col in (("water table (12-day change)", "d_wtd_cm"),
                          ("ERA5-Land soil moisture 0–7 cm (12-day change)", "d_soil_moisture"),
                          ("laser (12-day change, LOS)", "d_laser_los_mm")):
            a, b, c2, k = r_ci(c, g[col].to_numpy(float), rng)
            drv.append({"track": track, "driver": name, "r": a, "ci_lo": b, "ci_hi": c2, "n": k})
        for y, gy in g.groupby("year"):
            per.append({"track": track, "year": int(y), "triplets": len(gy),
                        "circular_mean_mm": circ_mean_ci(gy.closure_mm.to_numpy(), rng)[0],
                        "sd_mm": float(np.std(gy.closure_mm, ddof=1))})
    pd.DataFrame(drv).to_csv(OUT / "drivers.csv", index=False)
    pd.DataFrame(per).to_csv(OUT / "per_track.csv", index=False)
    v = {"verdict": verdict, "r_cross_track": r, "ci": [lo, hi], "matched_triplets": n,
         "r_cross_track_hard_ref": rr, "ci_hard_ref": [rlo, rhi], "snow_free_triplets": counts.to_dict(),
         "matched_by_year": m.a_year.value_counts().sort_index().to_dict()}
    (OUT / "verdict.json").write_text(json.dumps(v, indent=1, default=int), encoding="utf-8")
    (OUT / "README.md").write_text(readme(v, pd.DataFrame(per), pd.DataFrame(drv)), encoding="utf-8")
    print(json.dumps(v, indent=1, default=int))


def md(df: pd.DataFrame) -> str:
    d = df.copy().round(3)
    head = "| " + " | ".join(map(str, d.columns)) + " |\n|" + "---|" * len(d.columns) + "\n"
    return head + "\n".join("| " + " | ".join("" if pd.isna(x) else str(x) for x in r) + " |"
                            for r in d.itertuples(index=False))


def readme(v: dict, per: pd.DataFrame, drv: pd.DataFrame) -> str:
    return f"""# X-143 — closure of six-day triplets at P6 (R4 M2)

Generated by `05_code/SAr-adapted/experiments/closure_x143.py`; pre-registered in `preregistration/X-143.yaml`.
Closure c = wrap(w02 − w01 − w12), wrapped phase at P6 alone, mm of line of sight; motion and whole cycles cancel in it.

## Verdict

**{v["verdict"]}.** Cross-track correlation of the closures {v["r_cross_track"]:.2f} (95 % block-bootstrap interval
{v["ci"][0]:.2f} to {v["ci"][1]:.2f}; {v["matched_triplets"]} matched triplets); on the hard-target-referenced phase
{v["r_cross_track_hard_ref"]:.2f} ({v["ci_hard_ref"][0]:.2f} to {v["ci_hard_ref"][1]:.2f}).

## Per track

{md(per)}

## Closure against drivers

{md(drv)}
"""


if __name__ == "__main__":
    main()

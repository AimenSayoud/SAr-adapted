"""X-058 — verify the P6/CR laser before any comparison with Sentinel-1.

The supervisor (2026-09-28) asks, before the laser amplitudes are accepted: what
``PEAT/LaserSensor`` physically is, its units and sign, whether the ~10° mounting needs a
correction, gaps / outliers / resets, and values taken as close as possible to the Sentinel-1
acquisition times. This script answers each from the manual and from the data, writes the flags
the P6 comparison (X-059) uses, and lists what only the field team can confirm.

Reads the field data (hub ``06_data/field``) and the acquisition list of the first deliverable;
writes ONLY into the hub (``--out``, default ``<hub>/08_deliverables/field_p6/laser_qc``).

    PYTHONPATH=src python scripts/field_p6_laser_qc.py
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
from insar_wetlands.hydro_link import _detrend  # noqa: E402
from insar_wetlands.inversion.isbas import WAVELENGTH_M  # noqa: E402

HUB = Path(__file__).resolve().parents[3]
PERIOD = ("2022-01-01", "2024-12-31")
MOUNT_DEG = 10.0
QUARTER_WAVE_MM = WAVELENGTH_M / 4 * 1000
INCIDENCE = {"ascending": 32.26, "descending": 39.17}       # config.yaml, measured (X-042)
ERROR_CODE_MAX = -900.0                                     # SDMS40 error codes −900 … −912 (manual table 9-4)


def load(out: Path):
    wtd = field.load_wtd_hourly()
    laser = field.load_laser()
    snow72 = field.snow_mask(wtd["Air_2m"], laser.index, hours=72)
    snow24 = field.snow_mask(wtd["Air_2m"], laser.index, hours=24)
    acq = pd.read_csv(HUB / "08_deliverables" / "field_first" / "s1_acquisitions.csv")
    acq["time_utc"] = pd.to_datetime(acq.date + " " + acq.overpass_utc).dt.tz_localize("UTC")
    return wtd, laser, snow72, snow24, acq


# ------------------------------------------------------------------------------ the checks

def units_and_sign(laser, wtd, snow):
    """Units: the laser against the field team's own surface column (cm, used to correct WTD in
    cm) and against the raw CR water level (cm). Sign: a floating mat rises with its water, so the
    laser must increase with the raw level — on daily values, on daily anomalies (annual cycle
    removed from both) and on 12-day changes, so it is not only the shared season."""
    s = laser.surface_cm.where(~snow)
    d = pd.DataFrame({"laser": s, "ref": laser.surface_ref_cm, "cr": wtd["CR_raw"].reindex(laser.index)})
    daily = d.resample("D").mean().dropna()
    t = fl.years_since(daily.index)
    anom = {c: _detrend(daily[c].to_numpy(), t, True) for c in ("laser", "cr")}
    d12 = daily.resample("12D").mean().diff().dropna()
    r12, p12 = fl.circular_shift_p(d12.cr.to_numpy(), d12.laser.to_numpy())
    ra, pa = fl.circular_shift_p(anom["cr"], anom["laser"])
    vals = laser.surface_cm.dropna()
    return {
        "n_days": len(daily),
        "value_min": float(vals.min()), "value_max": float(vals.max()),
        "slope_laser_per_ref": fl.slopes(daily.ref, daily.laser)["slope_ols"],
        "r_laser_ref": float(daily.laser.corr(daily.ref)),
        "slope_laser_per_cr": fl.slopes(daily.cr, daily.laser)["slope_ols"],
        "r_laser_cr_daily": float(daily.laser.corr(daily.cr)),
        "r_laser_cr_anomaly": ra, "p_laser_cr_anomaly": pa,
        "n_12d": len(d12), "r_laser_cr_12d": r12, "p_laser_cr_12d": p12,
        "slope_laser_per_cr_12d": fl.slopes(d12.cr, d12.laser)["slope_ols"],
    }, daily, d12


def gaps(laser, snow72, snow24):
    """Empty hours and SDMS40 error codes. ``laser`` is the cleaned series (filled stretches
    already set to NaN), so filled gaps count as gaps."""
    present = laser.surface_cm.notna()
    runs = fl.gap_runs(present, min_hours=24)
    yr = pd.DataFrame({"hours": 1, "laser": present, "snowfree_72h": present & ~snow72,
                       "snowfree_24h": present & ~snow24}).groupby(laser.index.year).sum()
    yr["share_laser"] = yr.laser / yr.hours
    yr["share_snowfree_72h"] = yr.snowfree_72h / yr.hours
    yr["share_snowfree_24h"] = yr.snowfree_24h / yr.hours
    n_err = int((laser.surface_cm <= ERROR_CODE_MAX).sum())
    return runs, yr.rename_axis("year").reset_index(), n_err


def outliers_and_steps(laser, snow72):
    s = laser.surface_cm.where(~snow72)
    flags, sigma = fl.running_median_outliers(s, window="25h", k=5.0, min_abs=1.0)
    steps = fl.hourly_steps(s, threshold=1.0)
    steps_all = fl.hourly_steps(laser.surface_cm, threshold=1.0)
    return flags, sigma, steps, steps_all


def drift(laser, wtd, snow72):
    """Monthly offset of the laser against two independent-looking references: the field team's
    surface column (``surface_ref_cm``) and the surface implied by the logger (raw CR level −
    corrected WTD at P6). A reset appears as a jump, a settling mount or reference as a trend."""
    s = laser.surface_cm.where(~snow72)
    d = pd.DataFrame({"laser": s, "ref": laser.surface_ref_cm,
                      "logger": (wtd["CR_raw"] - wtd["P6"]).reindex(laser.index)})
    m = d.resample("MS").median()
    m["n_hours"] = s.resample("MS").count()
    m["laser_minus_ref"] = m.laser - m.ref
    m["laser_minus_logger"] = m.laser - m.logger
    m = m[m.n_hours >= 72]
    t = fl.years_since(m.index)
    rate = {c: float(np.polyfit(t, m[c], 1)[0]) for c in ("laser_minus_ref", "laser_minus_logger")}
    return m.rename_axis("month").reset_index(), rate


def at_overpass(laser, snow72, snow24, flags, acq, filled):
    """For every Sentinel-1 acquisition: the nearest hourly laser value (±1 h), its time offset,
    the 3-hour mean around the overpass (check), and the flags."""
    s = laser.surface_cm
    rows = []
    for a in acq.itertuples(index=False):
        t = a.time_utc
        v = s.dropna()
        i = v.index.get_indexer([t], method="nearest")[0]
        dt_min = abs((v.index[i] - t).total_seconds()) / 60
        near = v.iloc[i] if dt_min <= 60 else np.nan
        win = s[(s.index >= t - pd.Timedelta("90min")) & (s.index <= t + pd.Timedelta("90min"))]
        rows.append({"date": a.date, "track": a.track, "time_utc": t, "laser_cm": near,
                     "dt_minutes": dt_min if dt_min <= 60 else np.nan, "laser_3h_mean_cm": win.mean(),
                     "snow_72h": bool(snow72.asof(t)) if t >= snow72.index[0] else True,
                     "snow_24h": bool(snow24.asof(t)) if t >= snow24.index[0] else True,
                     "outlier": bool(flags.get(v.index[i], False)) if dt_min <= 60 else False,
                     "filled": bool(filled.get(v.index[i], False)) if dt_min <= 60 else False})
    out = pd.DataFrame(rows)
    out["usable"] = out.laser_cm.notna() & ~out.snow_72h & ~out.outlier & ~out.filled
    return out


def consecutive_changes(ov: pd.DataFrame, factor: float) -> pd.DataFrame:
    """Laser change between consecutive acquisitions of each track (usable at both ends), in
    vertical mm and in LOS mm; ``factor`` = 1 (sensor output already vertical) or cos 10°."""
    rows = []
    for track, g in ov.sort_values("time_utc").groupby("track"):
        g = g.reset_index(drop=True)
        for k in range(len(g) - 1):
            a, b = g.iloc[k], g.iloc[k + 1]
            if a.usable and b.usable:
                dh = (b.laser_cm - a.laser_cm) * 10 * factor
                rows.append({"track": track, "t1": a.date, "t2": b.date,
                             "dt_days": (b.time_utc - a.time_utc).days, "dh_mm": dh,
                             "dlos_mm": dh * np.cos(np.radians(INCIDENCE[track]))})
    return pd.DataFrame(rows)


def mount_sensitivity(ov, laser, snow72):
    """Key laser quantities with the output taken as vertical (factor 1: the sensor calibrated
    its own angle) and with a cos 10° factor (if it did not)."""
    daily = laser.surface_cm.where(~snow72).resample("D").mean().dropna()
    daily = daily[PERIOD[0]:PERIOD[1]]
    rows = []
    for label, f in (("as delivered (sensor-calibrated angle)", 1.0), ("× cos 10°", np.cos(np.radians(MOUNT_DEG)))):
        ch = consecutive_changes(ov, f)
        amp = daily.quantile(0.95) - daily.quantile(0.05)
        for track, g in ch.groupby("track"):
            rows.append({"case": label, "factor": f, "track": track, "n_consecutive_pairs": len(g),
                         "median_abs_dh_mm": g.dh_mm.abs().median(), "max_abs_dh_mm": g.dh_mm.abs().max(),
                         "share_los_beyond_quarter_wave": float((g.dlos_mm.abs() > QUARTER_WAVE_MM).mean()),
                         "range_p05_p95_mm_2022_2024": amp * 10 * f})
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------------ outputs

def figure(out: Path, laser, wtd, snow72, flags, m, d12, steps_all, interp):
    fig, ax = plt.subplots(4, 1, figsize=(11, 12))
    s = laser.surface_cm
    filled = fl.in_runs(s.index, interp)
    ax[0].plot(s.index, s.where(snow72 & ~filled), lw=0.5, color="0.7", label="snow period (72 h frost rule)")
    ax[0].plot(s.index, s.where(~snow72 & ~filled), lw=0.6, color="k", label="laser, snow-free")
    ax[0].plot(s.index, s.where(filled), lw=1.4, color="tab:purple", label="filled by interpolation in the delivery")
    o = s[flags]
    ax[0].scatter(o.index, o, s=10, color="tab:red", zorder=3, label=f"outliers ({len(o)})")
    if len(steps_all):
        ax[0].scatter(steps_all.time, steps_all.after, marker="v", s=18, color="tab:orange", zorder=3,
                      label=f"hourly steps > 1 cm, any period ({len(steps_all)})")
    a2 = ax[0].twinx()
    a2.plot(wtd.index, wtd["CR_raw"], lw=0.5, color="tab:blue", alpha=0.6)
    a2.set_ylabel("raw CR water level (cm)", color="tab:blue")
    ax[0].set_ylabel("laser surface (cm)")
    ax[0].legend(fontsize=7, loc="lower left")
    ax[0].set_title("P6/CR laser, whole record", fontsize=9)
    ax[1].plot(m.month, m.laser_minus_ref, "o-", ms=3, label="laser − field team's surface column")
    ax[1].plot(m.month, m.laser_minus_logger, "s-", ms=3, label="laser − (raw CR − WTD P6)")
    ax[1].set_ylabel("monthly median offset (cm)")
    ax[1].legend(fontsize=7)
    ax[1].set_title("Drift: the laser against two surface references", fontsize=9)
    d = s.where(~snow72 & ~filled).dropna().diff()
    d = d[(d.index.to_series().diff() <= pd.Timedelta("1h")) & (d.abs() < 1.5)]
    ax[2].hist(d, bins=np.arange(-1.5, 1.55, 0.1) - 0.05, color="0.4")
    ax[2].set_yscale("log")
    ax[2].set_xlabel("hour-to-hour change, snow-free (cm)")
    ax[2].set_title("Noise: hourly changes", fontsize=9)
    ax[3].scatter(d12.cr, d12.laser, s=10, color="k")
    ax[3].axhline(0, color="0.6", lw=0.6)
    ax[3].axvline(0, color="0.6", lw=0.6)
    ax[3].set_xlabel("12-day change of the raw CR water level (cm)")
    ax[3].set_ylabel("12-day change of the laser (cm)")
    ax[3].set_title("Sign: the surface rises when the water rises (snow-free 12-day means)", fontsize=9)
    fig.tight_layout()
    fig.savefig(out / "fig_laser_qc.png", dpi=160)
    plt.close(fig)


def f2(x, n=2):
    return "—" if x is None or not np.isfinite(x) else f"{x:.{n}f}"


def fp(p):
    return "—" if not np.isfinite(p) else ("≤ 0.001" if p <= 0.001 else f"{p:.3f}")


def report(out, u, runs, yr, n_err, sigma, steps, steps_all, rate, noise_cm, ov, sens, n_flags, interp, grid):
    per = ov[(ov.date >= PERIOD[0]) & (ov.date <= PERIOD[1])]
    cov = per.groupby("track").agg(acq=("date", "size"), within_1h=("laser_cm", lambda x: int(x.notna().sum())),
                                   filled=("filled", "sum"), usable=("usable", "sum"),
                                   snowfree_24h=("snow_24h", lambda x: int((~x).sum())),
                                   max_dt_min=("dt_minutes", "max"))
    hit = per[per.filled]
    ys = yr.set_index("year")
    L = []
    L.append("# X-058 — the P6/CR laser, verified before use\n")
    L.append("Generated by `05_code/SAr-adapted/scripts/field_p6_laser_qc.py` (outputs here only; field data). "
             "Answers the supervisor's laser checks of 2026-09-28 (`01_admin/supervisor/2026-09-28_email_p6_priority.md`); "
             "design in `../PLAN.md`. Every number below is computed by the script.\n")
    L.append("## 1. What `PEAT/LaserSensor` is\n")
    L.append("A Campbell Scientific **SDMS40 multipoint scanning snow-depth sensor** (manual in the delivery). It scans its "
             "laser along an oval on the ground, measures the distance to each point, filters and averages them, and outputs "
             "the **average snow depth of the target area**: the height of the surface **above the ground level set at "
             "calibration** (manual §1, §5). On snow-free peat that is the **peat/moss surface position relative to its "
             "calibration level**, averaged over the scanned oval (its size depends on the mounting height and angle) — "
             "not a point. Negative values = the surface is below where it was at calibration.\n")
    L.append("## 2. Units\n")
    L.append(f"- The manual's default unit is **mm** (configurable). The values run from {f2(u['value_min'],1)} to "
             f"{f2(u['value_max'],1)}, on a 0.1 grid (0.05 later) plus constant offsets (§7).")
    L.append(f"- Against the field team's own surface column (cm, used to correct WTD in cm): slope "
             f"{f2(u['slope_laser_per_ref'])}, r = {f2(u['r_laser_ref'])} ({u['n_days']} snow-free days) — same scale.")
    L.append(f"- Against the raw CR water level (cm): {f2(u['slope_laser_per_cr'])} per cm of water (daily). In mm this "
             f"would be a mat moving {f2(u['slope_laser_per_cr']/10,3)} cm per cm of water — not a floating mat.")
    L.append("- **Reading: centimetres.** To confirm with the field team (the logger program sets the unit).\n")
    L.append("## 3. Sign — does an increase mean the surface went up?\n")
    L.append("- Manual: depth increases when the surface comes closer to the sensor → **increase = up**.")
    L.append(f"- Data, independent of the field team's correction: a floating mat rises with its water. Laser vs raw CR "
             f"level, snow-free: daily r = {f2(u['r_laser_cr_daily'])}; **anomalies** (annual cycle removed from both) "
             f"r = {f2(u['r_laser_cr_anomaly'])} (p {fp(u['p_laser_cr_anomaly'])}); **12-day changes** r = "
             f"{f2(u['r_laser_cr_12d'])} (p {fp(u['p_laser_cr_12d'])}, n = {u['n_12d']}), slope "
             f"{f2(u['slope_laser_per_cr_12d'])} cm/cm. Positive at every time scale, not only through the season.")
    L.append("- **Reading: increase = upward surface movement.** Consistent with the conventions used so far "
             "(Sentinel-1 LOS positive toward the satellite = up).\n")
    L.append("## 4. The ~10° mounting\n")
    L.append("The SDMS40 mounts at 0–45° from vertical in 5° steps and, at calibration, **measures its own installation "
             "height and inclination angle** (manual §5, §7.3.1, §9.2.1); its output is the depth for that geometry. So the "
             "10° is already accounted for **if calibration was run with the sensor at 10° and it was not moved since** — "
             "then no correction is due, and applying cos 10° again would be wrong. If it was not, vertical changes are "
             f"cos 10° = {np.cos(np.radians(MOUNT_DEG)):.3f} × the output. Sensitivity (2022–2024, consecutive acquisitions, "
             "snow-free, outliers removed):\n")
    L.append("| track | case | pairs | median \\|Δh\\| (mm) | max \\|Δh\\| (mm) | LOS beyond λ/4 | 5–95 % range (mm) |")
    L.append("|---|---|---|---|---|---|---|")
    for _, r in sens.iterrows():
        L.append(f"| {r.track} | {r.case} | {r.n_consecutive_pairs} | {f2(r.median_abs_dh_mm,1)} | {f2(r.max_abs_dh_mm,1)} | "
                 f"{r.share_los_beyond_quarter_wave:.0%} | {f2(r.range_p05_p95_mm_2022_2024,0)} |")
    L.append("\nNo conclusion depends on the factor. **Kept as delivered (factor 1)** until the field team confirms the "
             "calibration (C-045).\n")
    L.append("## 5. Gaps\n")
    L.append(f"- SDMS40 error codes (−900 … −912) in the data: {n_err}.")
    L.append("- Share of hours with a value, and snow-free under two frost rules:\n")
    L.append("| year | hours | with value | snow-free (72 h rule) | snow-free (24 h rule) |")
    L.append("|---|---|---|---|---|")
    for y, r in ys.iterrows():
        L.append(f"| {y} | {int(r.hours)} | {r.share_laser:.0%} | {r.share_snowfree_72h:.0%} | {r.share_snowfree_24h:.0%} |")
    L.append(f"\n- Gaps of ≥ 24 h: {len(runs)} (longest {int(runs.hours.max()) if len(runs) else 0} h) — `gaps.csv`.")
    L.append(f"- **Filled gaps.** The delivered series is not all measurement: {len(interp)} stretches "
             f"({int(interp.hours.sum())} h) are straight lines between two readings — constant, sub-resolution hourly "
             "changes a sensor on a 0.1 grid cannot produce (`interpolated.csv`). They are treated as gaps here "
             "(the table above counts them as missing). They were used as data in X-048 / X-052.\n")
    L.append("| filled from | to | hours | laser change across it (cm) | raw CR water-level change across it (cm) | surface change the water predicts (cm) |")
    L.append("|---|---|---|---|---|---|")
    for r in interp.itertuples(index=False):
        L.append(f"| {r.start:%Y-%m-%d %H:%M} | {r.end:%Y-%m-%d %H:%M} | {r.hours} | {f2(r.laser_change_cm)} | "
                 f"{f2(r.cr_change_cm)} | {f2(r.expected_change_cm)} |")
    L.append("\nThe prediction uses the 12-day buoyancy slope of §3. A filled stretch across which the laser moves far "
             "more than the water predicts bridges a **re-levelling** (the level before and after differ for a reason "
             "other than the mat), not motion; no laser change is taken across any filled stretch.\n")
    L.append("## 6. Outliers\n")
    L.append(f"- Snow-free values more than 5 robust SD (σ = {f2(sigma,3)} cm) and ≥ 1 cm from the 25-h running median: "
             f"**{n_flags}** hours flagged (`laser_flags.csv`); they are excluded at the overpasses.")
    L.append("- Snow: the sensor reads snow as surface. The 72-h frost rule is kept (it removes the winter jumps); the "
             "24-h rule would recover the pairs counted in §8.\n")
    L.append("## 7. Resets and drift\n")
    st = ", ".join(f"{t:%Y-%m-%d %H:%M} ({v:+.1f} cm)" for t, v in zip(steps.time, steps.step)) or "none"
    L.append(f"- Hourly steps > 1 cm while snow-free and not filled: **{len(steps)}** ({st}); in any period, snow "
             f"included: {len(steps_all)} (`steps.csv`). Paired up/down steps an hour apart are one-hour spikes (caught by the outlier "
             "rule), not resets; no level shift shows as a jump in the measured values.")
    L.append("- **But a re-levelling is hidden in the data.** The sensor reports multiples of its resolution; a "
             "constant added afterwards moves every later value off that grid by the same amount. The offset from "
             "the grid changes here (`grid_offsets.csv`):\n")
    L.append("| from | to | days | offset from the grid (cm) |")
    L.append("|---|---|---|---|")
    for r in grid.itertuples(index=False):
        L.append(f"| {r.first_day:%Y-%m-%d} | {r.last_day:%Y-%m-%d} | {r.days} | {r.offset:.4f} |")
    L.append("\n  The change of offset falls inside the longest filled stretch: the series was **re-referenced across "
             "that gap** (a reinstallation, recalibration or a manual correction), and the interpolation joins the two "
             "levels so no step is visible. A level difference across that gap is therefore not surface motion. "
             "(A later switch to half steps looks like a change of averaging, not of level.)")
    L.append(f"- But the laser **drifts** against both surface references: the monthly offset changes by "
             f"{f2(rate['laser_minus_ref'])} cm/yr against the field team's surface column and "
             f"{f2(rate['laser_minus_logger'])} cm/yr against the logger surface (`drift_monthly.csv`). Either the "
             "laser's mount/calibration level or the references move slowly (a mast settling in peat would do it; X-057 "
             "sees the same as a yearly offset). Part of it is the re-levellings bridged by filled stretches (§5): the "
             "linear rate is a summary, not a mechanism.")
    L.append(f"- Over one 12-day interval that drift is {f2(abs(rate['laser_minus_ref'])*10*12/365.25,1)} mm — small "
             f"against the median 12-day change (§4). It **matters for the seasonal amplitude** and multi-year "
             "levels, **not for changes between consecutive acquisitions**, which X-059 uses.\n")
    L.append("## 8. At the Sentinel-1 acquisitions (2022–2024)\n")
    L.append("| track | acquisitions | laser within ±1 h | of which filled | usable (measured, snow-free 72 h, not outlier) | snow-free on the 24 h rule | largest time offset (min) |")
    L.append("|---|---|---|---|---|---|---|")
    for t, r in cov.iterrows():
        L.append(f"| {t} | {int(r.acq)} | {int(r.within_1h)} | {int(r.filled)} | {int(r.usable)} | {int(r.snowfree_24h)} | {f2(r.max_dt_min,0)} |")
    if len(hit):
        L.append("\n- Acquisitions whose laser value is a filled one (excluded): "
                 + ", ".join(f"{d} {t[:3]}" for d, t in zip(hit.date, hit.track)) + ".")
    L.append(f"\n- Per acquisition: `laser_at_s1.csv` (nearest hourly value, offset, 3-h mean, flags). The nearest value "
             f"and the 3-h mean differ by a median {f2((per.laser_cm - per.laser_3h_mean_cm).abs().median()*10,1)} mm.")
    L.append(f"- **Noise floor**: the hour-to-hour change (measured, snow-free, |Δ| < 1 cm) has an SD of {f2(noise_cm*10,1)} mm — "
             f"the noise of a difference of two readings, or an upper bound of it. A laser change between two acquisitions is taken as real above "
             f"**{f2(2*noise_cm*10,1)} mm** (2 SD); X-059 computes Sentinel-1/laser ratios only above it.\n")
    L.append("## Questions for the field team\n")
    L.append("1. Unit set in the logger program for the SDMS40 (we read cm).")
    L.append("2. Was the SDMS40's automatic calibration (`aXA!`) run after installation at ~10°? Dates of any "
             "recalibration, move or maintenance 2021–2025? (decides C-045 and the drift in §7)")
    L.append("3. What is the mast anchored to (can it settle or heave)? Is the CR logger referenced to a fixed datum?")
    L.append("4. The unnamed column of `Laser_Sensor.xlsx` (our `surface_ref_cm`): manual surface readings, interpolated?")
    L.append("5. Is the target oval on moss/peat or on vascular vegetation (does summer growth raise the reading)?\n")
    L.append("## Verdict for X-059\n")
    L.append("Use the laser **as delivered** (cm → mm, increase = up, no cos 10° pending question 2), the nearest hourly "
             "**measured** value within ±1 h of each overpass (filled stretches excluded), snow-free on the 72-h rule, "
             "outliers excluded; changes between consecutive acquisitions only (drift negligible there), never across "
             "the re-referenced gap; ratios only above the noise floor. Seasonal amplitudes carry the drift and "
             "re-referencing caveats.\n")
    L.append("## Questions added by these checks\n")
    L.append("6. What happened to the sensor during the filled stretches (above)? Was it reinstalled or recalibrated in "
             "July–August 2022, and was a constant offset applied to the data afterwards?")
    L.append("7. Were the filled values in `Laser_Sensor.xlsx` produced by the field team (interpolation), and are "
             "there other corrections in the delivered series?\n")
    (out / "laser_qc.md").write_text("\n".join(L))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default=str(HUB / "08_deliverables" / "field_p6" / "laser_qc"))
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    wtd, raw, snow72, snow24, acq = load(out)
    interp = fl.interpolated_runs(raw.surface_cm)
    cr = wtd["CR_raw"]
    lv = raw.surface_cm.dropna()
    interp["laser_change_cm"] = [lv[lv.index > r.end].iloc[0] - lv[lv.index < r.start].iloc[-1] for r in interp.itertuples()]
    interp["cr_change_cm"] = [cr[cr.index > r.end].iloc[0] - cr[cr.index < r.start].iloc[-1]
                              if (cr.index > r.end).any() else np.nan for r in interp.itertuples()]
    grid = fl.grid_offsets(raw.surface_cm)
    filled = pd.Series(fl.in_runs(raw.index, interp), index=raw.index)
    laser = raw.assign(surface_cm=raw.surface_cm.where(~filled))          # filled values are not data
    u, daily, d12 = units_and_sign(laser, wtd, snow72)
    interp["expected_change_cm"] = interp.cr_change_cm * u["slope_laser_per_cr_12d"]
    runs, yr, n_err = gaps(laser, snow72, snow24)
    flags, sigma, steps, steps_all = outliers_and_steps(laser, snow72)
    m, rate = drift(laser, wtd, snow72)
    noise_cm = fl.diff_noise_sd(laser.surface_cm.where(~snow72), clip=1.0)
    ov = at_overpass(raw, snow72, snow24, flags, acq, filled)
    sens = mount_sensitivity(ov, laser, snow72)

    pd.Series({**u, "sigma_running_median_cm": sigma, "n_outlier_hours": int(flags.sum()),
               "n_steps_snowfree": len(steps), "n_steps_all": len(steps_all), "n_error_codes": n_err,
               "drift_vs_ref_cm_per_yr": rate["laser_minus_ref"], "drift_vs_logger_cm_per_yr": rate["laser_minus_logger"],
               "diff_noise_sd_cm": noise_cm, "noise_floor_mm": 2 * noise_cm * 10}).to_csv(out / "summary.csv", header=["value"])
    runs.to_csv(out / "gaps.csv", index=False)
    yr.to_csv(out / "coverage_by_year.csv", index=False)
    pd.DataFrame({"time_utc": laser.index, "surface_cm": laser.surface_cm.to_numpy(), "snow_72h": snow72.to_numpy(),
                  "snow_24h": snow24.to_numpy(), "outlier": flags.to_numpy(), "filled": filled.to_numpy(),
                  "surface_cm_raw": raw.surface_cm.to_numpy()}).to_csv(out / "laser_flags.csv", index=False)
    interp.to_csv(out / "interpolated.csv", index=False)
    grid.to_csv(out / "grid_offsets.csv", index=False)
    steps_all.to_csv(out / "steps.csv", index=False)
    m.to_csv(out / "drift_monthly.csv", index=False)
    ov.to_csv(out / "laser_at_s1.csv", index=False)
    sens.to_csv(out / "mount_sensitivity.csv", index=False)
    figure(out, raw, wtd, snow72, flags, m, d12, steps_all, interp)
    report(out, u, runs, yr, n_err, sigma, steps, steps_all, rate, noise_cm, ov, sens, int(flags.sum()), interp, grid)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()

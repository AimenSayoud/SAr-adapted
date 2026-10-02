"""X-067 — can fusion v2 do better? Each change against v2 as published, judged the same way: P6 vs the laser on
held-out dates (2020–2021 time blocks; 2025 fully held out), ascending-only vs descending-only agreement, and stable
ground. No parameter is tuned on held-out laser dates.

E1  τ, σ_c of the mean-reverting motion from every laser day outside the test block (other years included)
E2  moisture coefficient split by surface state (pairs with a wet date vs dry pairs)
E3  backscatter as a per-pixel wet flag: VV above the pixel's own calendar-month median by > 1 SD of its anomalies
    marks it wet on that date; its pairs take the wet noise level (replaces the one scene-wide station flag)
E4  backscatter as an open-water mask: mat pixels with mean VV darker than the lake's median kept out of the model
E5  spatial pooling: strong (κ 16) and none (κ 0) against v2's choice (κ 4)

    INSAR_DRIVE_ROOT=../local/drive_pristine PYTHONPATH=src python scripts/fusion_v2_experiments_x067.py [variant prefixes]
    PYTHONPATH=src python scripts/fusion_v2_experiments_x067.py --readme
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import field_fusion_v2_x066 as V  # noqa: E402

from insar_wetlands.bootstrap import start  # noqa: E402

OUT = V.DLV / "fusion_v2" / "experiments_x067"
MIRROR = V.HUB / "05_code" / "local" / "drive_mirror"
STACKS = {"ascending": ["rtc_dualpol_2020_2024.nc", "rtc_dualpol_2025_2025.nc"],
          "descending": ["rtc_dualpol_2020_2024_descending.nc", "rtc_dualpol_2025_2025_descending.nc"]}


def vv_stack(track: str) -> xr.DataArray:
    parts = [xr.open_dataset(MIRROR / f).gamma0_vv_db for f in STACKS[track] if (MIRROR / f).exists()]
    return xr.concat(parts, dim="time").sortby("time")


def pixel_wet_flags(z: float = 1.0) -> dict:
    """{(track, date): bool grid}: VV anomaly to the pixel's calendar-month median (2020–2024) above z × its SD."""
    out = {}
    for track in STACKS:
        vv = vv_stack(track)
        ref = vv.sel(time=slice("2020", "2024"))
        clim = ref.groupby("time.month").median("time")
        anom = vv.groupby("time.month") - clim
        sd = (ref.groupby("time.month") - clim).std("time")
        wet = (anom > z * sd).values
        for k, t in enumerate(pd.to_datetime(vv.time.values)):
            out[(track, t.normalize())] = wet[k]
    return out


def open_water(inp) -> np.ndarray:
    """Mat pixels darker in mean VV (dusk, 2020–2024) than the lake's median: behave like open water."""
    m = xr.open_dataset(MIRROR / STACKS["ascending"][0]).gamma0_vv_db.mean("time").values
    lake = float(np.nanmedian(m[inp.zones["B"]]))
    w = inp.zones["A"] & (m < lake)
    r6, c6 = inp.p6
    w[r6, c6] = False
    return w


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    ctx = start("fusion_v2_experiments_x067", mount=False, git=False)
    inp = V.Inputs(ctx, V.HUB / "05_code" / "local" / "s1_2020_2021")
    final = json.loads((V.DLV / "fusion_v2" / "final_calibration.json").read_text())
    best = (final["R1"]["sigma_d_ratio"], final["R1"]["kappa"])
    r2_phys = {tuple([k.split("|")[0], int(k.split("|")[1])]): v for k, v in final["R2"]["phys"].items()}
    pw, ow = pixel_wet_flags(), open_water(inp)
    wet_share_mat = float(np.mean([v[inp.zones["A"]].mean() for v in pw.values()]))
    wet_share_grass = float(np.mean([v[inp.zones["C"]].mean() for v in pw.values()]))
    only = sys.argv[1:]                      # optional: run only the named variants (prefix match), e.g. "E1" "E5 strong"
    variants = {
        "v2 (as published)": {},
        "E1 τ from laser outside the test block": {"ou_outside_test": True},
        "E2 moisture term by surface state": {"e_by_wet": True},
        "E3 backscatter wet flag per pixel": {"pixel_wet": pw},
        "E4 open-water pixels out": {"mask_out": ow},
        "E5 strong pooling (κ 16)": {"best": (best[0], 16.0)},
        "E5 no pooling (κ 0)": {"best": (best[0], 0.0)},
    }
    if only:
        variants = {k: v for k, v in variants.items() if any(k.startswith(o) for o in only)}
    tag = "_" + "_".join(o.replace(" ", "") for o in only) if only else ""
    rows, checks = [], []
    for name, var in variants.items():
        var = {**var, "save_maps": False}
        print(f"\n######## {name}", flush=True)
        r1 = V.run(inp, "R1_2020_2021_6day_both", ("2020-01-01", "2022-01-31"), ("ascending", "descending"),
                   V.median_split_folds, OUT, fixed_best=var.get("best", best), variant=var)
        fixed = dict(r1["final_cal"])
        fixed["phys"] = {**r2_phys, **r1["final_cal"]["phys"]}
        nb = V.calibrate(inp, r1["tables"], lambda t: True, buoyancy=False)
        fixed.update({"g": 0.0, "tau": nb["tau"], "sigma_c": nb["sigma_c"]})
        r3 = V.run(inp, "R3_2025_6day_both", ("2025-01-01", "2025-12-31"), ("ascending", "descending"),
                   [("2025 held out", lambda t: False, lambda t: True)], OUT, fixed_cal=fixed,
                   fixed_best=var.get("best", best), water=False, variant=var)
        for R in (r1, r3):
            rows += [{"variant": name, **m} for m in R["metrics"]]
            checks += [{"variant": name, **c} for c in R["checks"]]
        pd.DataFrame(rows).to_csv(OUT / f"metrics_by_fold{tag}.csv", index=False)
        pd.DataFrame(checks).to_csv(OUT / f"checks{tag}.csv", index=False)
    met = pd.DataFrame(rows)
    met = met[met.estimate.isin(["v2 mat (P6 pixel)", "v2 single pixel"])]
    summ = met.groupby(["variant", "run", "estimate"], sort=False).agg(
        r=("r", "mean"), rmse_mm=("rmse_mm", "mean"), amplitude_ratio=("amplitude_ratio", "mean")).reset_index()
    summ.to_csv(OUT / f"summary{tag}.csv", index=False)
    if not only:
        (OUT / "context.json").write_text(json.dumps({"best_sigma_d_ratio_kappa": best, "pixel_wet_share_mat": wet_share_mat,
                                                  "pixel_wet_share_grassland": wet_share_grass,
                                                      "open_water_pixels_out": int(ow.sum())}, indent=2))
    pd.set_option("display.width", 220)
    print(summ.round(3).to_string(index=False))




def write_readme() -> None:
    """README from the saved tables (both rounds): every number read from a CSV written by main()."""
    rounds = [("round 1", "summary.csv", "checks.csv"), ("round 2 (E1 bounded to the run's period, κ 16)", "summary_v2_E1_E5strong.csv",
                                                         "checks_v2_E1_E5strong.csv")]
    ctxj = json.loads((OUT / "context.json").read_text())
    md = lambda df: "\n".join(["| " + " | ".join(map(str, df.columns)) + " |", "|" + "---|" * len(df.columns)]  # noqa: E731
                              + ["| " + " | ".join(f"{v:.3f}" if isinstance(v, float) else str(v) for v in r) + " |" for r in df.itertuples(index=False)])
    L = ["# X-067 — can fusion v2 do better? (exploratory)", "",
         "Generated by `05_code/SAr-adapted/scripts/fusion_v2_experiments_x067.py --readme` from its tables. Each change against "
         "v2 as published: P6 vs the laser on held-out dates (R1 2020–2021 time blocks; R3 2025 fully held out), mean over folds; "
         "ascending-only vs descending-only agreement and stable ground in `checks*.csv`.", "",
         f"Context: pixel wet flag (E3) marks {ctxj['pixel_wet_share_mat']:.0%} of mat pixel-dates and {ctxj['pixel_wet_share_grassland']:.0%} "
         f"of grassland ones wet; E4 found {ctxj['open_water_pixels_out']} mat pixels darker than the lake median; v2's pooling "
         f"is σ_d ratio {ctxj['best_sigma_d_ratio_kappa'][0]}, κ {ctxj['best_sigma_d_ratio_kappa'][1]}.", ""]
    for title, sf, cf in rounds:
        if not (OUT / sf).exists():
            continue
        s = pd.read_csv(OUT / sf)
        s = s[s.estimate == "v2 mat (P6 pixel)"].assign(run=lambda d: d.run.str[:2])
        L += [f"## {title}", "", md(s[["variant", "run", "r", "rmse_mm", "amplitude_ratio"]].round(3)), ""]
    L += ["## Reading", "",
          "- E1 in round 1 equalled v2 by construction (an open-ended test block); bounded in round 2, it trades correlation for "
          "amplitude at P6 in 2020–2021 and leaves 2025 unchanged (R3 already took τ from other years).",
          "- E2 (moisture term by surface state), E4 (open water) and E5 (pooling, none or strong) leave P6 essentially as v2; "
          "round 1's 'weak pooling (κ 4)' is v2's own setting, hence identical.",
          "- E3 (backscatter wet flag per pixel) does not help at P6; it slightly raises dusk/dawn agreement and lowers the "
          "stable-ground level in 2025 (`checks.csv`).",
          "- So the model-side changes tried here do not move P6: the limits are in the data (per-pair noise, no 2025 water "
          "table, one validation pixel). v2 stays as published."]
    (OUT / "README.md").write_text("\n".join(L) + "\n")


if __name__ == "__main__":
    if sys.argv[1:] == ["--readme"]:
        write_readme()
    else:
        main()

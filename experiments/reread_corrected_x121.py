"""X-121 — re-run the cell- and plot-scale tests unchanged on the HyP3-corrected inputs (X-120), keeping the originals.

For every entry, in an order where each sees its corrected upstream: (1) back up all the entries' declared output
folders (hub ``05_code/local/x121_backup/<UTC>/``, outside git); (2) run each entry through the register with the
relocating crop reader on (``RZECIN_HYP3_RELOCATE=1``) — the cube is already corrected; (3) copy every new or changed
output file to ``insar_full_test/reread_corrected/<entry>/``; (4) restore the originals from the backup (files the runs
created are moved, not deleted). Writes outside the declared folders are listed. Rule: ``preregistration/X-121.yaml``.

    INSAR_DRIVE_ROOT=<hub>/05_code/local/drive_mirror make run E=reread_corrected_x121
"""
from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

from insar_wetlands.paths import hub_root
from insar_wetlands.register import load_register as load
from insar_wetlands.scoring import assert_preregistered

HUB = hub_root()
FT = HUB / "08_deliverables" / "insar_full_test"
OUT = FT / "reread_corrected"
PREREG = FT / "preregistration" / "X-121.yaml"
REPO = Path(__file__).resolve().parents[1]
ORDER = ["field_first_deliverables", "field_p6_validation", "field_p6_short_pairs_x065", "field_boardwalk_x060",
         "field_plots_x061", "field_p6_cells_x079", "baseline_x099", "slc_p6_x102", "fusion_v3_wp3_x076",
         "loss_of_lock_x096", "matched_wt_x094", "bias_precheck_x097"]


def sha(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def files_under(d: Path) -> dict[str, str]:
    return {str(f.relative_to(HUB)): sha(f) for f in d.rglob("*") if f.is_file()} if d.exists() else {}


def hub_dirs(entry) -> list[Path]:
    out = []
    for o in entry.outputs:
        if str(o).startswith("hub:"):
            p = HUB / str(o)[4:]
            out.append(p if p.suffix == "" else p.parent)
    return out


def git_status() -> set[str]:
    r = subprocess.run(["git", "status", "--porcelain"], cwd=HUB, capture_output=True, text=True)
    return {line[3:] for line in r.stdout.splitlines()}


SNAPSHOT = HUB / "06_data" / "cube_snapshots" / "as_delivered_2026-10-15"


def main() -> None:
    control = "--control" in sys.argv          # amendment: today's code on the as-delivered inputs
    out_root = OUT / "_control_as_delivered" if control else OUT
    assert_preregistered(PREREG, [out_root / "runs.json"])
    reg = load()
    entries = {n: reg[n] for n in ORDER}
    dirs = sorted({d for e in entries.values() for d in hub_dirs(e)}, key=lambda p: len(p.parts))
    dirs = [d for d in dirs if not any(d != o and o in d.parents for o in dirs)]       # outermost folders only
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    backup = HUB / "05_code" / "local" / "x121_backup" / stamp
    before = {}
    for d in dirs:
        if d.exists():
            shutil.copytree(d, backup / d.relative_to(HUB))
        before.update(files_under(d))
    git0 = git_status()
    env = dict(os.environ, RZECIN_HYP3_RELOCATE="0" if control else "1", PYTHONPATH=str(REPO / "src"))
    if control:
        env["RZECIN_CUBE_ROOT"] = str(SNAPSHOT)
    runs = []
    for name in ORDER:
        t0 = datetime.now(timezone.utc).isoformat()
        r = subprocess.run([sys.executable, "-m", "insar_wetlands.register", "run", name], cwd=REPO, env=env,
                           capture_output=True, text=True)
        runs.append({"entry": name, "started": t0, "exit": r.returncode, "tail": (r.stdout + r.stderr)[-1500:]})
        print(name, "exit", r.returncode, flush=True)
    moved = []
    for d in dirs:
        after = files_under(d)
        for rel, h in after.items():
            if before.get(rel) == h:
                continue
            src = HUB / rel
            owner = next((n for n in ORDER for dd in hub_dirs(entries[n]) if dd in src.parents or dd == src.parent), "other")
            dst = out_root / owner / src.relative_to(next(dd for dd in dirs if dd in src.parents))
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dst)
            if rel in before:
                shutil.copy2(backup / rel, src)          # the original back in place
            else:
                src.unlink()                             # created by the re-run: now lives only in reread_corrected
            moved.append({"file": rel, "to": str(dst.relative_to(HUB)), "was_new": rel not in before})
    outside = sorted(p for p in git_status() - git0 if not p.startswith(str(OUT.relative_to(HUB))))
    out_root.mkdir(parents=True, exist_ok=True)
    (out_root / "runs.json").write_text(json.dumps({"backup": str(backup.relative_to(HUB)), "runs": runs, "moved": moved,
                                               "changed_outside_declared": outside}, indent=1), encoding="utf-8")
    print(f"moved {len(moved)} files; outside the declared folders: {outside}")


# ------------------------------------------------------------------------------ the comparison (control vs corrected)

SPECS = [  # (test, entry, file under the entry's output folder, keys, values)
    ("X-059", "field_p6_validation", "p6_stats.csv", ["track", "window", "subset"],
     ["n_pairs", "r_s1_vs_laser", "p_s1_vs_laser", "slope_ols", "slope_tls"]),
    ("X-065", "field_p6_short_pairs_x065", "phase_vs_laser_by_revisit.csv", ["track", "period", "dt_days"],
     ["n", "r", "same_sign", "m", "m_sd"]),
    ("X-060", "field_boardwalk_x060", "p6_windows.csv", ["track", "window", "subset"],
     ["n_pairs", "r_s1_vs_laser", "p_s1_vs_laser", "slope_tls"]),
    ("X-060", "field_boardwalk_x060", "boardwalk_tests_summary.csv", ["variable", "track", "res_m"],
     ["median_anomaly_boardwalk", "share_beyond_null_90", "p_mannwhitney"]),
    ("X-061", "field_plots_x061", "plots_groups.csv", ["track", "period", "extraction", "subset", "group"],
     ["median_r_s1_vs_dwtd", "n_p05_s1_vs_dwtd", "median_slope_mm_per_cm", "median_r_coh_vs_abs_dwtd",
      "n_p05_coh_vs_abs_dwtd"]),
    ("X-079", "field_p6_cells_x079", "cells.csv", ["track", "period", "cell"],
     ["n", "r", "p_circ", "slope_tls", "tls_lo95", "tls_hi95", "m"]),
    ("X-079", "field_p6_cells_x079", "drivers.csv", ["track", "period", "slope", "against"], ["r", "p_perm"]),
    ("X-076", "fusion_v3_wp3_x076", "wp3/e31_shape_vs_raft.csv",
     ["period", "track", "grid", "revisit_days", "shape", "subset"], ["slope", "r", "perm_p", "above_patch_null"]),
    ("X-094", "matched_wt_x094", "verdict.csv", ["track"], ["matched_share_mat", "unmatched_share_mat", "keeps_phase"]),
    ("X-097", "bias_precheck_x097", "verdicts.csv", ["correction", "track", "class"], ["verdict"]),
    ("X-099", "baseline_x099", "levels.csv", ["track", "reference", "view"],
     ["pearson_detrended", "rmse_detrended_mm", "amplitude_ratio"]),
    ("X-099", "baseline_x099", "increments.csv", ["track", "reference"],
     ["pearson_increments", "rmse_increments_mm", "slope_through_origin", "within_quarter_wave"]),
    ("X-099", "baseline_x099", "winters.csv", ["track", "reference", "winter"], ["difference_mm"]),
    ("X-102", "slc_p6_x102", "verdict.csv", ["variant"], ["rmse_change_asc", "rmse_change_desc", "result"]),
    ("X-102", "slc_p6_x102", "scores.csv", ["track", "variant"], ["rmse_increments_mm", "pearson_increments"]),
]


def _version(entry: str, rel: str, which: str):
    import pandas as pd
    reg = load()
    orig = hub_dirs(reg[entry])[0] / rel
    cands = ([OUT / "_control_as_delivered" / entry / rel] if which == "control" else [OUT / entry / rel]) + [orig]
    f = next((c for c in cands if c.exists()), None)
    return (pd.read_csv(f), f) if f is not None else (None, None)


def _flag(stat: str, a, b) -> str:
    import math
    if isinstance(a, str) or isinstance(b, str) or isinstance(a, bool) or isinstance(b, bool):
        return "changed" if str(a) != str(b) else ""
    if a is None or b is None or (isinstance(a, float) and math.isnan(a)) or (isinstance(b, float) and math.isnan(b)):
        return "" if (a is None or a != a) and (b is None or b != b) else "one side missing"
    if stat.startswith("p_") or stat in ("perm_p", "p_perm", "p_circ"):
        return "crosses 0.05" if (a < 0.05) != (b < 0.05) else ""
    if stat in ("r", "slope", "slope_tls", "slope_ols", "m", "r_s1_vs_laser", "pearson_increments", "difference_mm",
                "median_r_s1_vs_dwtd", "median_r_coh_vs_abs_dwtd", "median_slope_mm_per_cm", "pearson_detrended",
                "median_anomaly_boardwalk") and a * b < 0 and max(abs(a), abs(b)) > 0.05:
        return "sign flips"
    return ""


def compare() -> None:
    import pandas as pd
    rows = []
    for test, entry, rel, keys, vals in SPECS:
        c, fc = _version(entry, rel, "control")
        k, fk = _version(entry, rel, "corrected")
        if c is None or k is None:
            continue
        m = c.merge(k, on=keys, how="outer", suffixes=("__c", "__k"))
        for _, r in m.iterrows():
            for v in vals:
                a, b = r.get(v + "__c"), r.get(v + "__k")
                rows.append({"test": test, "table": f"{entry}/{rel}", "key": " · ".join(str(r[x]) for x in keys),
                             "stat": v, "control": a, "corrected": b,
                             "delta": (b - a) if isinstance(a, (int, float)) and isinstance(b, (int, float))
                             and not isinstance(a, bool) else None, "flag": _flag(v, a, b)})
    D = pd.DataFrame(rows)
    D.to_csv(OUT / "comparison.csv", index=False)
    S = (D.assign(flagged=D.flag != "").groupby(["test", "table"]).agg(values=("stat", "size"),
                                                                        flagged=("flagged", "sum")).reset_index())
    F = D[D.flag != ""]
    fmt = lambda x: f"{x:.3f}" if isinstance(x, float) else str(x)  # noqa: E731
    md = lambda df: "\n".join(["| " + " | ".join(map(str, df.columns)) + " |", "|" + "---|" * len(df.columns)] +  # noqa: E731
                              ["| " + " | ".join(fmt(x) for x in r) + " |" for r in df.itertuples(index=False)])
    (OUT / "COMPARISON.md").write_text(f"""# X-121 — the cell- and plot-scale tests, as delivered vs corrected (generated)

Generated by `experiments/reread_corrected_x121.py --compare`; do not edit. Rule: `../preregistration/X-121.yaml`
(with its amendment). **control** = today's code on the as-delivered inputs (`_control_as_delivered/`, or the committed
deliverable where the control reproduced it exactly); **corrected** = the same code on the HyP3-corrected inputs.
Every value: `comparison.csv`. A flag marks a sign that flips, a p-value that crosses 0.05, a verdict or boolean that
changes, or a value present on one side only — the reading of each test is in the ticket (X-121) and the try.

## Values and flags per table

{md(S)}

## Flagged values

{md(F[["test", "table", "key", "stat", "control", "corrected", "flag"]]) if len(F) else "None."}
""", encoding="utf-8")
    print(S.to_string(index=False))
    print(F.groupby(["test", "flag"]).size().to_string())


if __name__ == "__main__":
    if "--compare" in sys.argv:
        compare()
    else:
        main()

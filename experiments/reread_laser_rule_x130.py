"""X-130 — re-read the article's laser comparisons and the supervisor tables with the cube's P6 laser on the
snow-evidence rule (X-129), beside the committed outputs (rule A, the 72-h frost rule).

As X-121: (1) back up the entries' declared output folders (hub ``05_code/local/x130_backup/<UTC>/``, outside git);
(2) run each entry unchanged through the register, in order; X-099's outputs stand in for
``reread_corrected/baseline_x099`` while the later entries run (X-126 checks against them); (3) copy every new or
changed file to ``insar_full_test/reread_laser_rule_b/<entry>/``; (4) restore the originals. ``--control`` runs the same
on a rule-A cube into ``reread_laser_rule_b/_control_rule_A/``; ``--compare`` writes the comparison (tables side by
side, every article number that moves). Rule: ``preregistration/X-130.yaml``.

    make run E=reread_laser_rule_x130 ARGS=--control     # on the rule-A cube (before the switch)
    make run E=reread_laser_rule_x130                     # on the rule-B cube
    make run E=reread_laser_rule_x130 ARGS=--compare
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from reread_corrected_x121 import files_under, git_status, hub_dirs  # noqa: E402
from insar_wetlands.paths import cube_root, hub_root  # noqa: E402
from insar_wetlands.register import load_register as load  # noqa: E402
from insar_wetlands.scoring import assert_preregistered  # noqa: E402

HUB = hub_root()
FT = HUB / "08_deliverables" / "insar_full_test"
OUT = FT / "reread_laser_rule_b"
CONTROL = OUT / "_control_rule_A"
PREREG = FT / "preregistration" / "X-130.yaml"
REPO = Path(__file__).resolve().parents[1]
ORDER = ["baseline_x099", "x065_hard_targets_x126", "skip_pairs_x127", "article_fig_f2", "article_fig_f4",
         "article_fig_f6", "supervisor_table_c071"]
STAND_IN = (FT / "baseline", FT / "reread_corrected" / "baseline_x099")   # X-099's output → where its readers look
NUMBERS = HUB / "08_deliverables" / "article_plan" / "numbers.yaml"


def cube_rule() -> str:
    s = pd.read_csv(cube_root() / "silver" / "site_hourly.csv.gz", usecols=["laser_snow_72h_source"], low_memory=False)
    return "snow_evidence_72h" if (s.laser_snow_72h_source == "snow_evidence").any() else "frost_72h"


def run_pass(control: bool) -> None:
    out_root = CONTROL if control else OUT
    want = "frost_72h" if control else "snow_evidence_72h"
    if cube_rule() != want:
        sys.exit(f"the cube's laser rule is {cube_rule()}, this pass needs {want}")
    assert_preregistered(PREREG, [out_root / "runs.json"])
    reg = load()
    entries = {n: reg[n] for n in ORDER}
    dirs = sorted({d for e in entries.values() for d in hub_dirs(e)} | {STAND_IN[1]}, key=lambda p: len(p.parts))
    dirs = [d for d in dirs if not any(d != o and o in d.parents for o in dirs)]
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    backup = HUB / "05_code" / "local" / "x130_backup" / stamp
    before = {}
    for d in dirs:
        if d.exists():
            shutil.copytree(d, backup / d.relative_to(HUB))
        before.update(files_under(d))
    git0 = git_status()
    env = dict(os.environ, PYTHONPATH=str(REPO / "src"))
    runs = []
    for name in ORDER:
        t0 = datetime.now(timezone.utc).isoformat()
        r = subprocess.run([sys.executable, "-m", "insar_wetlands.register", "run", name], cwd=REPO, env=env,
                           capture_output=True, text=True)
        runs.append({"entry": name, "started": t0, "exit": r.returncode, "tail": (r.stdout + r.stderr)[-1500:]})
        print(name, "exit", r.returncode, flush=True)
        if name == "baseline_x099":
            for f in STAND_IN[0].glob("*.csv"):
                shutil.copy2(f, STAND_IN[1] / f.name)
    moved = []
    for d in dirs:
        for rel, h in files_under(d).items():
            if before.get(rel) == h:
                continue
            src = HUB / rel
            if STAND_IN[1] in src.parents:                       # the stand-in copy: restore only
                shutil.copy2(backup / rel, src) if rel in before else src.unlink()
                continue
            owner = next((n for n in ORDER for dd in hub_dirs(entries[n]) if dd in src.parents or dd == src.parent),
                         "other")
            dst = out_root / owner / src.relative_to(next(dd for dd in dirs if dd in src.parents))
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dst)
            if rel in before:
                shutil.copy2(backup / rel, src)
            else:
                src.unlink()
            moved.append({"file": rel, "to": str(dst.relative_to(HUB)), "was_new": rel not in before})
    outside = sorted(p for p in git_status() - git0 if not p.startswith(str(OUT.relative_to(HUB))))
    out_root.mkdir(parents=True, exist_ok=True)
    (out_root / "runs.json").write_text(json.dumps({"cube_rule": want, "backup": str(backup.relative_to(HUB)),
                                                    "runs": runs, "moved": moved,
                                                    "changed_outside_declared": outside}, indent=1), encoding="utf-8")
    print(f"moved {len(moved)} files; outside the declared folders: {outside}")


# ------------------------------------------------------------------------------ the comparison

def committed_of(rel: str) -> Path:
    """The committed file a re-read file stands beside (X-099's corrected outputs live in reread_corrected)."""
    p = HUB / rel
    return STAND_IN[1] / p.name if STAND_IN[0] == p.parent else p


def table_diff(a: pd.DataFrame, b: pd.DataFrame) -> dict:
    same_shape = a.shape == b.shape and list(a.columns) == list(b.columns)
    out = {"rows_A": len(a), "rows_B": len(b), "same_columns": list(a.columns) == list(b.columns)}
    if same_shape:
        num = [c for c in a.columns if pd.api.types.is_numeric_dtype(a[c]) and pd.api.types.is_numeric_dtype(b[c])]
        x, y = a[num].to_numpy(float), b[num].to_numpy(float)
        same = (np.isnan(x) & np.isnan(y)) | (np.abs(x - y) <= 1e-9)
        oth = [c for c in a.columns if c not in num]
        txt = lambda d: d[oth].astype(object).where(d[oth].notna(), "<NA>").astype(str).to_numpy()   # NaN = NaN
        so = (txt(a) == txt(b)) if oth else np.ones((len(a), 0), bool)
        out["cells_changed"] = int((~same).sum() + (~so).sum())
    else:
        out["cells_changed"] = None
    return out


def resolve(entry: dict, path: Path) -> float | None:
    d = pd.read_csv(path)
    if "where" in entry:
        for k, v in entry["where"].items():
            d = d[d[k].astype(str) == str(v)]
        if len(d) != 1:
            return None
        v = d[entry["column"]].iloc[0]
    else:
        v = d[entry["column"]].iloc[int(entry["row"])]
    try:
        v = float(v) * float(entry.get("scale", 1)) * (100 if entry.get("percent") else 1)
    except (TypeError, ValueError):
        return None
    return round(v, int(entry.get("decimals", 2)))


def compare() -> None:
    def moved(root):
        f = root / "runs.json"
        return json.loads(f.read_text())["moved"] if f.exists() else []
    lines = ["# X-130 — the P6 laser on the snow-evidence rule: what moves (generated)", "",
             "Generated by SAr-adapted `experiments/reread_laser_rule_x130.py --compare`; do not edit. Rule:",
             "`../preregistration/X-130.yaml`. A = the committed outputs (72-h frost rule); B = the re-read on the",
             "snow-evidence cube (`<entry>/`); control = today's code on the rule-A cube (`_control_rule_A/`).", ""]
    ctrl = []
    for m in moved(CONTROL):
        if not m["file"].endswith((".csv", ".csv.gz")):
            continue
        a = committed_of(m["file"])
        ctrl.append({"file": str(a.relative_to(HUB)), **table_diff(pd.read_csv(a), pd.read_csv(HUB / m["to"]))}
                    if a.exists() else {"file": m["file"], "rows_A": None, "rows_B": None, "cells_changed": "new"})
    c = pd.DataFrame(ctrl)
    lines += ["## Control: today's code on the rule-A cube against the committed tables", "",
              "Only files that differ from the committed ones were kept; a CSV listed with 0 changed cells differs in "
              "formatting only.", "", md(c) if len(c) else "(no CSV differs: the control reproduces every committed "
                                                            "table)", ""]
    rows, remap = [], {}
    for m in moved(OUT):
        a = committed_of(m["file"])
        remap[str(a.relative_to(HUB))] = HUB / m["to"]
        if m["file"].endswith((".csv", ".csv.gz")) and a.exists():
            rows.append({"file": str(a.relative_to(HUB)), **table_diff(pd.read_csv(a), pd.read_csv(HUB / m["to"]))})
    t = pd.DataFrame(rows)
    lines += ["## Tables that change under B", "", md(t) if len(t) else "(none)", ""]
    nums = yaml.safe_load(NUMBERS.read_text())["numbers"]
    nr = []
    for n in nums:
        if n["csv"] not in remap:
            continue
        va, vb = resolve(n, HUB / n["csv"]), resolve(n, remap[n["csv"]])
        if va != vb:
            nr.append({"id": n["id"], "A": va, "B": vb, "table": Path(n["csv"]).name, "note": n.get("note", "")[:90]})
    nt = pd.DataFrame(nr)
    lines += ["## Article numbers that move (numbers.yaml)", "",
              f"{len(nt)} of the {sum(n['csv'] in remap for n in nums)} numbers read from a re-read table.", "",
              md(nt) if len(nt) else "(none)", ""]
    for name, keys in (("x065_hard_targets/verdict.csv", ["track"]),
                       ("x065_hard_targets/by_revisit.csv", ["track", "period", "revisit_days", "pairs_kind", "reference"]),
                       ("reread_corrected/baseline_x099/years.csv", None),
                       ("reread_corrected/baseline_x099/levels.csv", None)):
        a = FT / name
        b = remap.get(str(a.relative_to(HUB)))
        if b is None:
            continue
        A, B = pd.read_csv(a), pd.read_csv(b)
        keys = keys or [c for c in A.columns if not pd.api.types.is_float_dtype(A[c])]
        j = A.merge(B, on=keys, how="outer", suffixes=(" A", " B"))
        lines += [f"## {name}, side by side", "", md(j), ""]
    (OUT / "COMPARISON.md").write_text("\n".join(lines), encoding="utf-8")
    nt.to_csv(OUT / "numbers_moved.csv", index=False)
    t.to_csv(OUT / "tables_changed.csv", index=False)
    print(f"wrote {OUT / 'COMPARISON.md'}: {len(t)} tables change, {len(nt)} numbers move")


def md(df: pd.DataFrame, nd: int = 3) -> str:
    d = df.copy()
    for c in d.columns:
        if pd.api.types.is_float_dtype(d[c]):
            d[c] = d[c].map(lambda x: "" if pd.isna(x) else f"{x:.{nd}f}")
    return "\n".join(["| " + " | ".join(map(str, d.columns)) + " |", "|" + "---|" * len(d.columns)] +
                     ["| " + " | ".join(str(x) for x in r) + " |" for r in d.itertuples(index=False)])


def main() -> None:
    if "--compare" in sys.argv:
        compare()
    else:
        run_pass(control="--control" in sys.argv)


if __name__ == "__main__":
    main()

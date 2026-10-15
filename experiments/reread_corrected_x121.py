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


if __name__ == "__main__":
    main()

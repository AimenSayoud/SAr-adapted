"""X-123 — re-run the cell- and plot-level backscatter results unchanged with the RTC put back under its labels (X-122),
against a control pass on the RTC as delivered, keeping the originals.

Both passes read the HyP3 layers corrected (``RZECIN_HYP3_RELOCATE=1``, X-120). ``--control`` runs first, on the cube
before its RTC layer is rebuilt, with ``RZECIN_RTC_RELOCATE=0``; the main pass runs after the rebuild with
``RZECIN_RTC_RELOCATE=1``. Each pass backs up the declared output folders (hub ``05_code/local/x123_backup/<UTC>/``,
outside git), runs the entries, moves new or changed files to ``insar_full_test/reread_rtc/[_control_rtc_as_delivered/]
<entry>/`` and restores the originals (``reread_corrected_x121``'s procedure). ``--compare`` matches every table of the
two passes. Rule: ``preregistration/X-123.yaml``.

    INSAR_DRIVE_ROOT=<hub>/05_code/local/drive_mirror make run E=reread_rtc_x123 ARGS="--control"
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from reread_corrected_x121 import files_under, git_status, hub_dirs  # noqa: E402

from insar_wetlands.paths import hub_root  # noqa: E402
from insar_wetlands.register import load_register as load  # noqa: E402
from insar_wetlands.scoring import assert_preregistered  # noqa: E402

HUB = hub_root()
FT = HUB / "08_deliverables" / "insar_full_test"
OUT = FT / "reread_rtc"
CONTROL = OUT / "_control_rtc_as_delivered"
PREREG = FT / "preregistration" / "X-123.yaml"
REPO = Path(__file__).resolve().parents[1]
ORDER = ["backscatter_x063", "field_first_deliverables", "field_pixel_x053", "field_boardwalk_x060",
         "field_fusion_x062", "field_p6_short_pairs_x065", "fusion_v2_experiments_x067", "closure_moisture_x070",
         "cube_analysis_x081", "scattering_regime_x098"]


def run_pass(control: bool) -> None:
    out_root = CONTROL if control else OUT
    assert_preregistered(PREREG, [out_root / "runs.json"])
    reg = load()
    entries = {n: reg[n] for n in ORDER}
    dirs = sorted({d for e in entries.values() for d in hub_dirs(e)}, key=lambda p: len(p.parts))
    dirs = [d for d in dirs if not any(d != o and o in d.parents for o in dirs)]       # outermost folders only
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    backup = HUB / "05_code" / "local" / "x123_backup" / stamp
    before = {}
    for d in dirs:
        if d.exists():
            shutil.copytree(d, backup / d.relative_to(HUB))
        before.update(files_under(d))
    git0 = git_status()
    env = dict(os.environ, RZECIN_HYP3_RELOCATE="1", RZECIN_RTC_RELOCATE="0" if control else "1",
               PYTHONPATH=str(REPO / "src"))
    runs = []
    for name in ORDER:
        t0 = datetime.now(timezone.utc).isoformat()
        r = subprocess.run([sys.executable, "-m", "insar_wetlands.register", "run", name], cwd=REPO, env=env,
                           capture_output=True, text=True)
        runs.append({"entry": name, "started": t0, "exit": r.returncode, "tail": (r.stdout + r.stderr)[-1500:]})
        print(name, "exit", r.returncode, flush=True)
    moved = []
    for d in dirs:
        for rel, h in files_under(d).items():
            if before.get(rel) == h:
                continue
            src = HUB / rel
            owner = next((n for n in ORDER for dd in hub_dirs(entries[n]) if dd in src.parents or dd == src.parent),
                         "other")
            dst = out_root / owner / src.relative_to(next(dd for dd in dirs if dd in src.parents))
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dst)
            if rel in before:
                shutil.copy2(backup / rel, src)          # the original back in place
            else:
                src.unlink()                             # created by the re-run: now lives only in reread_rtc
            moved.append({"file": rel, "to": str(dst.relative_to(HUB)), "was_new": rel not in before})
    outside = sorted(p for p in git_status() - git0 if not p.startswith(str(OUT.relative_to(HUB))))
    out_root.mkdir(parents=True, exist_ok=True)
    (out_root / "runs.json").write_text(json.dumps({"backup": str(backup.relative_to(HUB)), "runs": runs,
                                                    "moved": moved, "changed_outside_declared": outside}, indent=1),
                                        encoding="utf-8")
    print(f"moved {len(moved)} files; outside the declared folders: {outside}")


# ------------------------------------------------------------------------------ the comparison (control vs corrected)

def _flag(col: str, a, b) -> str:
    import math
    if isinstance(a, (str, bool)) or isinstance(b, (str, bool)):
        return "changed" if str(a) != str(b) else ""
    fa = a is None or (isinstance(a, float) and math.isnan(a))
    fb = b is None or (isinstance(b, float) and math.isnan(b))
    if fa or fb:
        return "" if fa and fb else "one side missing"
    c = col.lower()
    if c.startswith("p_") or c.endswith("_p") or c in ("p", "pval", "p_value", "perm_p", "p_perm", "p_circ") \
            or "pvalue" in c or c.startswith("p_perm"):
        return "crosses 0.05" if (a < 0.05) != (b < 0.05) else ""
    if a * b < 0 and max(abs(a), abs(b)) > 0.05:
        return "sign flips"
    return ""


def _pair(entry: str, rel: Path):
    import pandas as pd
    c, k = CONTROL / entry / rel, OUT / entry / rel
    orig = hub_dirs(load()[entry])[0] / rel
    c = c if c.exists() else orig                          # a file the control left unchanged equals the original
    k = k if k.exists() else orig
    if not (c.exists() and k.exists()):
        return None, None
    try:
        return pd.read_csv(c), pd.read_csv(k)
    except Exception:  # noqa: BLE001 — not a table
        return None, None


def compare() -> None:
    import numpy as np
    import pandas as pd
    rels = set()
    for root in (CONTROL, OUT):
        for entry in ORDER:
            d = root / entry
            if d.exists():
                rels |= {(entry, f.relative_to(d)) for f in d.rglob("*.csv")}
    rows = []
    for entry, rel in sorted(rels):
        if str(rel).startswith("_") or rel.name.startswith("comparison.csv"):
            continue
        c, k = _pair(entry, rel)
        if c is None:
            continue
        num = [x for x in c.columns if x in k.columns and pd.api.types.is_float_dtype(c[x])]
        keys = [x for x in c.columns if x in k.columns and x not in num]
        if keys and not c.duplicated(keys).any() and not k.duplicated(keys).any():
            m = c.merge(k, on=keys, how="outer", suffixes=("__c", "__k"))
            label = lambda r, keys=keys: " · ".join(str(r[x]) for x in keys)  # noqa: E731
        else:
            m = c.add_suffix("__c").join(k.add_suffix("__k"), how="outer")
            keys = []
            label = lambda r: f"row {r.name}"  # noqa: E731
        vals = num if keys else [x for x in c.columns if x in k.columns]
        for _, r in m.iterrows():
            for v in vals:
                a, b = r.get(v + "__c"), r.get(v + "__k")
                if isinstance(a, np.generic):
                    a = a.item()
                if isinstance(b, np.generic):
                    b = b.item()
                num_ab = isinstance(a, (int, float)) and isinstance(b, (int, float)) and not isinstance(a, bool)
                rows.append({"entry": entry, "table": str(rel), "key": label(r), "stat": v, "control": a,
                             "corrected": b, "delta": (b - a) if num_ab else None, "flag": _flag(v, a, b)})
    D = pd.DataFrame(rows)
    D.to_csv(OUT / "comparison.csv.gz", index=False)          # gzipped: the plain table passes 50 MB
    D["changed"] = D.delta.fillna(0).abs() > 1e-9
    S = (D.assign(flagged=D.flag != "").groupby(["entry", "table"])
         .agg(values=("stat", "size"), changed=("changed", "sum"), flagged=("flagged", "sum")).reset_index())
    F = D[D.flag != ""]
    fmt = lambda x: f"{x:.3f}" if isinstance(x, float) else str(x)  # noqa: E731
    md = lambda df: "\n".join(["| " + " | ".join(map(str, df.columns)) + " |", "|" + "---|" * len(df.columns)] +  # noqa: E731
                              ["| " + " | ".join(fmt(x) for x in r) + " |" for r in df.itertuples(index=False)])
    (OUT / "COMPARISON.md").write_text(f"""# X-123 — the backscatter results, RTC as delivered vs corrected (generated)

Generated by `experiments/reread_rtc_x123.py --compare`; do not edit. Rule: `../preregistration/X-123.yaml`.
**control** = the RTC as delivered (`_control_rtc_as_delivered/`, or the committed file where the control reproduced
it exactly); **corrected** = the same code with the RTC under its labels. Both read the HyP3 layers corrected (X-120).
Every value: `comparison.csv.gz`. A flag marks a sign that flips (|value| > 0.05), a p-value that crosses 0.05, a
boolean or text that changes, or a value present on one side only; the reading per test is in `README.md`.

## Values, changed values and flags per table

{md(S)}

## Flagged values

{md(F[["entry", "table", "key", "stat", "control", "corrected", "flag"]]) if len(F) else "None."}
""", encoding="utf-8")
    print(S.to_string(index=False))
    print(F.groupby(["entry", "flag"]).size().to_string() if len(F) else "no flags")


def smoothing() -> None:
    """Amendment diagnostic: does the bilinear relocation of the nearest-resampled RTC smooth speckle? The as-delivered
    layer (the snapshot), moved by the measured shift and by the mirrored shift: the same interpolation weights in the
    wrong direction, so what both change is smoothing, not placement. Amplitude dispersion and the dB temporal spread,
    median per zone, 2020–2024."""
    import numpy as np
    import pandas as pd
    import xarray as xr

    from insar_wetlands import placement as P
    snap = HUB / "06_data" / "cube_snapshots" / "rtc_as_delivered_2026-10-16" / "silver"
    zone = xr.open_dataset(snap / "static_40m.nc").zone.values
    rows = []
    for t in ("ascending", "descending"):
        r = xr.open_dataset(snap / f"rtc_{t}.nc").sel(time=slice("2020", "2024"))
        dx, dy = P.RTC_CONTENT_SHIFT_M[t]
        lin = 10 ** (r.vv_db.values / 10)
        for name, v in (("as delivered", lin), ("moved by the shift", P.relocate(lin, r.x.values, r.y.values, dx, dy)),
                        ("moved by the mirrored shift", P.relocate(lin, r.x.values, r.y.values, -dx, -dy))):
            a = np.sqrt(v)
            with np.errstate(all="ignore"):
                da, sd = np.nanstd(a, 0) / np.nanmean(a, 0), np.nanstd(10 * np.log10(v), 0)
            for z, zn in ((1, "mat"), (3, "grassland")):
                rows.append({"track": t, "layer": name, "zone": zn, "amplitude_dispersion": float(np.nanmedian(da[zone == z])),
                             "vv_temporal_sd_db": float(np.nanmedian(sd[zone == z]))})
    d = pd.DataFrame(rows)
    out = OUT / "_bilinear_withdrawn"
    d.to_csv(out / "smoothing.csv", index=False)
    md = "\n".join(["| " + " | ".join(d.columns) + " |", "|" + "---|" * len(d.columns)] +
                   ["| " + " | ".join(f"{x:.3f}" if isinstance(x, float) else str(x) for x in r) + " |"
                    for r in d.itertuples(index=False)])
    (out / "SMOOTHING.md").write_text(f"""# X-123 amendment — the bilinear relocation smooths the RTC (generated)

Generated by `experiments/reread_rtc_x123.py --smoothing`; do not edit. The RTC on the 40 m grid is nearest-resampled
(one 10 m sample per cell); moving it by bilinear interpolation averages neighbouring samples. Moved by the mirrored
shift (same weights, wrong direction), the layer changes as much as when moved by the measured shift: the change is
smoothing, not placement. The first corrected pass (this folder) is therefore withdrawn; the correction is remade by
re-sampling the 10 m source at the corrected positions (`pipelines/rtc_download.py --placed`).

{md}
""", encoding="utf-8")
    print(d.to_string(index=False))


if __name__ == "__main__":
    if "--smoothing" in sys.argv:
        smoothing()
    elif "--compare" in sys.argv:
        compare()
    else:
        run_pass(control="--control" in sys.argv)

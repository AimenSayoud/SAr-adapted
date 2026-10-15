"""X-109 decisive test — is the HyP3 interferogram grid placed right? Per pair, against our right-placed SLC coherence.

For every consecutive 2021 pair HyP3 and the X-101 CSLC stack share, per track: the shift of HyP3's coherence (5 m
steps, +-3 cells) that maximises its correlation with our 40 m cell coherence (unfiltered). Offset if the median shift
over pairs exceeds 20 m on an axis and its bootstrap interval excludes zero. Rule: X-109's amendment of 2026-10-15
(``08_deliverables/insar_full_test/preregistration/X-109.yaml``). Reads the CSLCs (``INSAR_DRIVE_ROOT``) and the cube;
writes only to the hub: ``08_deliverables/insar_full_test/hyp3_placement/decisive/``.

    INSAR_DRIVE_ROOT=<hub>/05_code/local/drive_mirror make run E=placement_decisive_x109
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parent))
import slc_qc_x101 as Q  # noqa: E402

from insar_wetlands.cube import cube_root  # noqa: E402
from insar_wetlands.paths import hub_root, resolve_drive  # noqa: E402
from insar_wetlands.scoring import assert_preregistered  # noqa: E402

HUB = hub_root()
FT = HUB / "08_deliverables" / "insar_full_test"
OUT = FT / "hyp3_placement" / "decisive"
PREREG = FT / "preregistration" / "X-109.yaml"
STEP, MAXS, CELL, TOL = 5.0, 120.0, 40.0, 20.0
SHIFTS = np.arange(-MAXS, MAXS + STEP / 2, STEP)


def best_shift(hy: np.ndarray, ours: np.ndarray) -> tuple[float, float, float, float]:
    """(dx east m, dy north m, r at best, r unshifted): the shift applied to HyP3's image that best matches ours."""
    v = np.isfinite(hy)
    img = np.where(v, hy, 0.0)
    best = (0.0, 0.0, -np.inf)
    r0 = np.nan
    for dy in SHIFTS:
        for dx in SHIFTS:
            sh = ndimage.shift(img, (-dy / CELL, dx / CELL), order=1, mode="constant", cval=0.0)
            sv = ndimage.shift(v.astype(float), (-dy / CELL, dx / CELL), order=1, mode="constant", cval=0.0) > 0.999
            m = sv & np.isfinite(ours)
            if m.sum() < 200:
                continue
            r = np.corrcoef(sh[m], ours[m])[0, 1]
            if dx == 0 and dy == 0:
                r0 = r
            if r > best[2]:
                best = (float(dx), float(dy), float(r))
    return best[0], best[1], best[2], float(r0)


def main() -> None:
    assert_preregistered(PREREG, [OUT / "verdict.csv"])
    OUT.mkdir(parents=True, exist_ok=True)
    st = xr.open_dataset(cube_root() / "silver" / "static_40m.nc")
    gx, gy = st.x.values, st.y.values
    shape = (len(gy), len(gx))
    rows = []
    for track in ("ascending", "descending"):
        folder = resolve_drive() / "slc_pilot" / "2021" / "cslc" / track
        files = {Q.DATE_RE.search(f.name).group(1): f for f in sorted(folder.glob("*.h5"))
                 if "static" not in f.name.lower() and Q.DATE_RE.search(f.name)}
        dates = sorted(files)
        hy = Q.hyp3_pairs(track)
        for a, b in zip(dates[:-1], dates[1:]):
            if (a, b) not in hy:
                continue
            s1, x, y = Q.cslc(files[a])
            s2, _, _ = Q.cslc(files[b])
            rr, cc = Q.axis_index(x, y, gx, gy)
            _, coh = Q.look(s1, s2, rr, cc, shape)
            dx, dy, r, r0 = best_shift(hy[(a, b)][1], coh)
            rows.append({"track": track, "t1": a, "t2": b, "dx_m": dx, "dy_m": dy, "r_best": r, "r_unshifted": r0})
            print(track, a, b, dx, dy, round(r, 3), flush=True)
    p = pd.DataFrame(rows)
    rng = np.random.default_rng(0)
    v = []
    for track, g in p.groupby("track"):
        rec = {"track": track, "pairs": len(g)}
        off = False
        for ax in ("dx_m", "dy_m"):
            x = g[ax].to_numpy()
            med = float(np.median(x))
            boot = np.median(rng.choice(x, (2000, len(x)), replace=True), 1)
            lo, hi = np.percentile(boot, [2.5, 97.5])
            rec.update({f"median_{ax}": med, f"ci_lo_{ax}": float(lo), f"ci_hi_{ax}": float(hi)})
            off |= abs(med) > TOL and not (lo <= 0 <= hi)
        rec["verdict"] = "offset" if off else "placed right"
        v.append(rec)
    v = pd.DataFrame(v)
    p.to_csv(OUT / "pairs.csv", index=False)
    v.to_csv(OUT / "verdict.csv", index=False)
    (OUT / "README.md").write_text(readme(p, v), encoding="utf-8")
    print(v.to_string(index=False))


def _md(df: pd.DataFrame, nd: int = 1) -> str:
    d = df.copy()
    for col in d.columns:
        if pd.api.types.is_float_dtype(d[col]):
            d[col] = d[col].map(lambda x: "—" if pd.isna(x) else f"{x:.{nd}f}")
    return "\n".join(["| " + " | ".join(map(str, d.columns)) + " |", "|" + "---|" * len(d.columns)] +
                     ["| " + " | ".join(map(str, r)) + " |" for r in d.itertuples(index=False)])


def readme(p, v) -> str:
    s = p.groupby("track")[["r_unshifted", "r_best"]].median().reset_index()
    return f"""# X-109 decisive test — HyP3 interferogram grid against our right-placed SLC coherence (generated)

Generated by `experiments/placement_decisive_x109.py`; do not edit. Rule: the amendment of 2026-10-15 in
`../../preregistration/X-109.yaml`. Per consecutive 2021 pair shared by HyP3 and the X-101 CSLCs: the shift of HyP3's
coherence (5 m steps, +-3 cells; dx east, dy north) best matching our 40 m unfiltered coherence; offset if the median
over pairs exceeds 20 m on an axis and its bootstrap interval (2000 resamples) excludes zero.

## Verdict

{_md(v)}

## Median correlation, unshifted vs at the best shift

{_md(s, 3)}
"""


if __name__ == "__main__":
    main()

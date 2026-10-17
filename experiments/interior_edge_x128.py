"""X-128 — do six-day pairs keep phase more often in the mat's interior than at its edge once the area around the plot
transect and its boardwalk is left out?

The article's §5.2 compares the median per-cell share of six-day pairs above each pair's open-water floor in the mat's
interior and at its edge (``tables/short_pairs_by_depth.csv``), and notes that much of the contrast lies near the plot
transect. Here the same per-cell shares (the C-067 computation behind that table and F5) are compared outside a zone
fixed in advance — mat cells within three cells (120 m, Chebyshev) of a plot or of the boardwalk — with a block
bootstrap over the cube's spatial blocks, since neighbouring cells are not independent. Radar only: no laser is read.

Pre-registration: hub ``08_deliverables/insar_full_test/preregistration/X-128.yaml`` (checked before the run). Reads the
cube; writes only to the hub (``08_deliverables/insar_full_test/interior_edge/``). No field value in this file.

    make run E=interior_edge_x128
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import maximum_filter

sys.path.insert(0, str(Path(__file__).resolve().parent))
import article_figures_c067 as C67  # noqa: E402  (static(), pair_table(): the computation behind Table F5 and §5.2)
from insar_wetlands.paths import hub_root  # noqa: E402
from insar_wetlands.scoring import assert_preregistered  # noqa: E402

HUB = hub_root()
FT = HUB / "08_deliverables" / "insar_full_test"
OUT = FT / "interior_edge"
PREREG = FT / "preregistration" / "X-128.yaml"
EXCLUDE_CELLS = {"primary (120 m)": 3, "sensitivity (200 m)": 5}
N_BOOT, SEED, MIN_CELLS = 2000, 0, 30


def near(mask: np.ndarray, cells: int) -> np.ndarray:
    """Cells within ``cells`` (Chebyshev distance) of any True cell of ``mask``."""
    return maximum_filter(mask.astype(np.uint8), size=2 * cells + 1, mode="constant") > 0


def contrast(share, interior, edge, blocks, rng) -> dict:
    """Median interior minus median edge, with a block-bootstrap 95 % interval over spatial blocks."""
    diff = float(np.nanmedian(share[interior]) - np.nanmedian(share[edge]))
    ids = np.unique(blocks[interior | edge])
    by_block = {b: (share[interior & (blocks == b)], share[edge & (blocks == b)]) for b in ids}
    draws = []
    for _ in range(N_BOOT):
        pick = rng.choice(ids, size=len(ids), replace=True)
        si = np.concatenate([by_block[b][0] for b in pick])
        se = np.concatenate([by_block[b][1] for b in pick])
        if np.isfinite(si).any() and np.isfinite(se).any():
            draws.append(np.nanmedian(si) - np.nanmedian(se))
    lo, hi = np.nanpercentile(draws, [2.5, 97.5])
    return {"interior_median": float(np.nanmedian(share[interior])), "edge_median": float(np.nanmedian(share[edge])),
            "difference": diff, "diff_lo95": float(lo), "diff_hi95": float(hi), "blocks": len(ids)}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    assert_preregistered(PREREG, [OUT / "contrast.csv", OUT / "README.md"])
    st, masks, water = C67.static()
    mat = masks["mat"]
    depth = st.depth_in_mat_m.values
    split = float(np.nanmedian(depth[mat]))                     # the split of short_pairs_by_depth.csv (whole mat)
    blocks = st.spatial_block.values
    anchors = (st.plot_cell.values > 0) | (np.nan_to_num(st.boardwalk_share.values) > 0)
    rng = np.random.default_rng(SEED)
    rows = []
    for track in ("ascending", "descending"):
        _, extra = C67.pair_table(track, masks, water)
        share = extra["above6"] / max(extra["n6"], 1)
        for label, cells in [("whole mat", 0)] + list(EXCLUDE_CELLS.items()):
            keep = mat & ~(near(anchors, cells) if cells else np.zeros_like(mat))
            interior, edge = keep & (depth > split), keep & (depth <= split)
            row = {"track": track, "area": label, "split_depth_m": split, "six_day_pairs": extra["n6"],
                   "mat_cells_kept": int(keep.sum()), "interior_cells": int(interior.sum()), "edge_cells": int(edge.sum())}
            if min(interior.sum(), edge.sum()) < MIN_CELLS:
                rows.append({**row, "verdict": "not judged (too few cells)"})
                continue
            c = contrast(share, interior, edge, blocks, rng)
            holds = c["difference"] > 0 and c["diff_lo95"] > 0
            rows.append({**row, **c, "verdict": "holds" if holds else "does not hold"})
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "contrast.csv", index=False)
    d = t.copy()
    for c in d.columns:
        if pd.api.types.is_float_dtype(d[c]):
            d[c] = d[c].map(lambda x: f"{x:.2f}")
    table = "\n".join(["| " + " | ".join(d.columns) + " |", "|" + "---|" * len(d.columns)] +
                      ["| " + " | ".join(str(x) for x in row) + " |" for row in d.itertuples(index=False)])
    (OUT / "README.md").write_text(f"""# X-128 — interior vs edge, with and without the transect area (generated)

Generated by `experiments/interior_edge_x128.py` (`make run E=interior_edge_x128`); do not edit. Pre-registration:
`../preregistration/X-128.yaml`. Per-cell share of six-day pairs above the pair's open-water floor (C-067's
`pair_table`, 2020–2021 set, corrected cube); interior = deeper than the whole mat's median depth; excluded = mat cells
within the given distance of a plot cell or a boardwalk cell; 95 % interval of the difference by a block bootstrap over
the cube's spatial blocks ({N_BOOT} resamples, seed {SEED}). The primary is the 120 m exclusion.

{table}
""", encoding="utf-8")
    print(t[["track", "area", "difference", "diff_lo95", "diff_hi95", "verdict"]].to_string())


if __name__ == "__main__":
    main()

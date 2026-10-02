"""Finer products onto the 40 m template (D-024: Sentinel-1 2026 at 5x1 looks, 20 m pixel spacing).

Every analysis of the project runs on one 40 m template grid. A finer product is brought onto it by averaging,
per template cell, the fine pixels whose centres fall inside that cell — so it works whether or not the fine grid
is aligned with the template. Coherence and unwrapped phase are averaged as numbers; wrapped phase as a phasor
(the angle of the mean of e^{iφ}), so that values near ±π do not cancel to 0.
"""
from __future__ import annotations

import numpy as np


def cell_index(fine_x: np.ndarray, fine_y: np.ndarray, tpl_x: np.ndarray, tpl_y: np.ndarray):
    """Row and column of the template cell containing each fine pixel centre (−1 outside the template)."""
    dx, dy = tpl_x[1] - tpl_x[0], tpl_y[1] - tpl_y[0]
    col = np.floor((fine_x - (tpl_x[0] - dx / 2)) / dx).astype(int)
    row = np.floor((fine_y - (tpl_y[0] - dy / 2)) / dy).astype(int)
    col[(col < 0) | (col >= len(tpl_x))] = -1
    row[(row < 0) | (row >= len(tpl_y))] = -1
    return row, col


def to_template(values: np.ndarray, fine_x, fine_y, tpl_x, tpl_y, kind: str = "mean", min_count: int = 1) -> np.ndarray:
    """Average a fine 2-D field (y, x) onto the template. ``kind``: "mean" (coherence, unwrapped phase) or
    "phase" (wrapped phase, circular mean). Cells with fewer than ``min_count`` finite fine pixels are NaN."""
    v = np.asarray(values, float)
    r, c = cell_index(np.asarray(fine_x, float), np.asarray(fine_y, float), np.asarray(tpl_x, float), np.asarray(tpl_y, float))
    R, C = np.meshgrid(r, c, indexing="ij")
    ok = (R >= 0) & (C >= 0) & np.isfinite(v)
    flat = R[ok] * len(tpl_x) + C[ok]
    n_cells = len(tpl_y) * len(tpl_x)
    count = np.bincount(flat, minlength=n_cells).astype(float)
    if kind == "phase":
        re = np.bincount(flat, weights=np.cos(v[ok]), minlength=n_cells)
        im = np.bincount(flat, weights=np.sin(v[ok]), minlength=n_cells)
        out = np.arctan2(im, re)
    elif kind == "mean":
        out = np.bincount(flat, weights=v[ok], minlength=n_cells) / np.where(count > 0, count, 1)
    else:
        raise ValueError(kind)
    out[count < min_count] = np.nan
    return out.reshape(len(tpl_y), len(tpl_x))


KIND = {"wrapped_phase": "phase", "corr": "mean", "unw_phase": "mean"}


def stack_to_template(da, tpl_x, tpl_y, kind: str | None = None, min_count: int = 2):
    """A (pair, y, x) stack from ``stack.load_layer`` on a finer grid → the same stack on the template grid
    (``kind`` from the layer name when not given: wrapped phase as a phasor, coherence and unwrapped phase as means).
    A grid that already is the template is returned unchanged."""
    import xarray as xr

    if da.sizes["x"] == len(tpl_x) and da.sizes["y"] == len(tpl_y) and np.allclose(da.x, tpl_x) and np.allclose(da.y, tpl_y):
        return da
    k = kind or KIND.get(str(da.name), "mean")
    vals = np.stack([to_template(da.values[i], da.x.values, da.y.values, tpl_x, tpl_y, kind=k, min_count=min_count)
                     for i in range(da.sizes["pair"])])
    out = xr.DataArray(vals, dims=("pair", "y", "x"), coords={"pair": da.pair.values, "y": tpl_y, "x": tpl_x}, name=da.name)
    for c in ("ref_date", "sec_date"):
        if c in da.coords:
            out = out.assign_coords({c: ("pair", da[c].values)})
    return out

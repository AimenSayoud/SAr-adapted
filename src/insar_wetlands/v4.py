"""Mat method v4 (X-106, selected by X-119): the chain of consecutive short pairs, each winter's level step set by an
independent anchor (SPAMS, weather only).

Pure functions; the data access and the frozen parameters (field-derived, kept in the hub) live in the experiment
``freeze_v4_x106`` and in ``08_deliverables/insar_full_test/v4/v4_frozen_params.yaml``.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def winter_bounds(times: pd.Series, season_months=(5, 11), spring_window=(3, 6), datum=None) -> list[tuple[int, int]]:
    """Index pairs (last in-season epoch of year y, first spring-window epoch of y+1); with ``datum`` (a Series of
    datum labels aligned with ``times``) a winter spanning two datums is skipped."""
    t = pd.to_datetime(pd.Series(times).reset_index(drop=True))
    out = []
    for y in sorted(t.dt.year.unique()):
        a = t.index[(t.dt.year == y) & t.dt.month.between(*season_months)]
        b = t.index[(t.dt.year == y + 1) & t.dt.month.between(*spring_window)]
        if len(a) and len(b):
            if datum is not None and datum.iloc[a[-1]] != datum.iloc[b[0]]:
                continue
            out.append((int(a[-1]), int(b[0])))
    return out


def relevel(c: pd.DataFrame, anchor_at, bounds: list[tuple[int, int]]) -> pd.DataFrame:
    """``c``: a chain (time, seg, est) in time order. Each winter step (a → b) is set to ``anchor_at(time_b) −
    anchor_at(time_a)``; later epochs shift by the difference; the segments either side of the winter are joined."""
    c = c.reset_index(drop=True).copy()
    for a, b in bounds:
        step = anchor_at(c.time.iloc[b]) - anchor_at(c.time.iloc[a])
        if not np.isfinite(step):
            continue
        delta = step - (c.est.iloc[b] - c.est.iloc[a])
        c.loc[b:, "est"] += delta
        old = c.seg.iloc[b]
        c.loc[(c.index >= b) & (c.seg == old), "seg"] = c.seg.iloc[a]
    return c


def anchor_only(c: pd.DataFrame, anchor_at) -> pd.DataFrame:
    """The anchor alone, read at the chain's epochs (X-107 amendment_4: the weather model with no radar): time, seg,
    est with ``est = anchor_at(time)`` and one segment, since nothing breaks a weather series at a gap or a winter;
    epochs where the anchor is undefined are dropped."""
    out = pd.DataFrame({"time": c.time.reset_index(drop=True), "seg": 0})
    out["est"] = [anchor_at(t) for t in out.time]
    return out[np.isfinite(out.est.to_numpy(float))].reset_index(drop=True)


def block_cells(x: np.ndarray, y: np.ndarray, xq: float, yq: float, size: int) -> tuple[np.ndarray, np.ndarray]:
    """Row and column indices of the cell block read by the neighbour-average variant (X-125) on a regular grid with
    cell-centre coordinates ``x`` (ascending) and ``y`` (descending). ``size`` 2: the four cells whose centres surround
    (xq, yq); ``size`` 3: the cell containing (xq, yq) and its eight neighbours."""
    dx, dy = float(x[1] - x[0]), float(y[1] - y[0])
    fc, fr = (xq - x[0]) / dx, (yq - y[0]) / dy                   # fractional column and row of the point
    if size == 2:
        c0, r0 = int(np.floor(fc)), int(np.floor(fr))
        cols, rows = np.array([c0, c0 + 1]), np.array([r0, r0 + 1])
    elif size == 3:
        c0, r0 = int(np.rint(fc)), int(np.rint(fr))
        cols, rows = np.arange(c0 - 1, c0 + 2), np.arange(r0 - 1, r0 + 2)
    else:
        raise ValueError("size must be 2 or 3")
    if cols.min() < 0 or rows.min() < 0 or cols.max() >= len(x) or rows.max() >= len(y):
        raise ValueError("block leaves the grid")
    return rows, cols


def block_mean(values: np.ndarray, rows: np.ndarray, cols: np.ndarray, min_share: float = 0.5) -> np.ndarray:
    """Mean over the block of ``values`` (..., y, x), NaN where fewer than ``min_share`` of its cells are finite."""
    b = values[..., rows[:, None], cols[None, :]].reshape(*values.shape[:-2], -1)
    ok = np.isfinite(b)
    with np.errstate(invalid="ignore"):
        m = np.nanmean(np.where(ok, b, np.nan), axis=-1)
    return np.where(ok.mean(axis=-1) >= min_share, m, np.nan)

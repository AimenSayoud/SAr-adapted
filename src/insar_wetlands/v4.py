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

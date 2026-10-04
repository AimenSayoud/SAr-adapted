"""Data contracts of the cube's gold tables: keys that must be unique, columns that must exist, and physical ranges.

``validate`` returns a list of problems (empty = the table honours its contract); ``check`` in ``gold.py`` fails the build
on any. Ranges are physical plausibility bounds, not expectations: a value outside them is an ingestion error.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

PI = float(np.pi) + 1e-6

CONTRACTS: dict[str, dict] = {
    "plots_acq": {
        "keys": ["track", "date", "plot"],
        "required": ["time_utc", "period", "unit", "row", "col", "air_c", "flag_wet_any", "flag_frozen_air"],
        "ranges": {"vv_db_3x3": (-40, 15), "vh_db_3x3": (-45, 10), "rvi_3x3": (0, 4), "air_c": (-40, 45),
                   "wtd_at": (-200, 80), "s2_ndvi_3x3": (-1, 1), "s2_ndmi_3x3": (-1, 1), "lst_c_3x3": (-40, 70),
                   "uav_ndvi": (-1, 1), "uav_lai": (0, 15), "s2_lag_days": (-10, 10), "lst_lag_days": (-16, 16),
                   "uav_lag_days": (-30, 30), "eco_lst_c_3x3": (-40, 70)},
    },
    "plots_pair": {
        "keys": ["period", "track", "grid", "pair", "plot"],
        "required": ["t1", "t2", "revisit_days", "wrapped_1x1", "unw_1x1", "coh_1x1", "ref_stable_unw",
                     "ref_grassland_unw", "fold_year"],
        "ranges": {"revisit_days": (1, 800), "coh_1x1": (0, 1), "coh_3x3": (0, 1), "wrapped_1x1": (-PI, PI),
                   "wrapped_3x3": (-PI, PI), "consistency_3x3": (0, 1 + 1e-9), "rain_sum_mm": (0, 3000)},
    },
    "triangles": {
        "keys": ["period", "track", "grid", "p12", "p23", "p13", "unit"],
        "required": ["closure_wrapped", "closure_consistency", "cycle_share"],
        "ranges": {"closure_wrapped": (-PI, PI), "closure_consistency": (0, 1 + 1e-9), "cycle_share": (0, 1)},
    },
    "truth_laser_daily": {
        "keys": ["date"],
        "required": ["surface_cm", "n_ok_hours"],
        "ranges": {"surface_cm": (-100, 100), "n_ok_hours": (0, 24)},
    },
}


def validate(df: pd.DataFrame, name: str) -> list[str]:
    c = CONTRACTS[name]
    out = []
    miss = [k for k in c["keys"] + c["required"] if k not in df.columns]
    if miss:
        out.append(f"{name}: missing columns {miss}")
        return out
    dup = int(df.duplicated(c["keys"]).sum())
    if dup:
        out.append(f"{name}: {dup} duplicated keys {c['keys']}")
    for col, (lo, hi) in c["ranges"].items():
        if col in df.columns:
            v = pd.to_numeric(df[col], errors="coerce")
            bad = int(((v < lo) | (v > hi)).sum())
            if bad:
                out.append(f"{name}: {bad} values of {col} outside [{lo}, {hi}]")
    return out

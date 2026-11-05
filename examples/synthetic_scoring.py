"""A deterministic scoring example; no field data, files, credentials or cloud jobs."""
from __future__ import annotations

import numpy as np
import pandas as pd

from insar_wetlands.scoring import increment_stats


def synthetic_series() -> tuple[pd.DatetimeIndex, np.ndarray]:
    """Return the demo's deterministic epochs and displacement (mm)."""
    epochs = pd.date_range("2021-05-01", periods=30, freq="6D", tz="UTC")
    truth_mm = 10 * np.sin(np.linspace(0, 2 * np.pi, len(epochs)))
    return epochs, truth_mm


def main() -> None:
    epochs, truth_mm = synthetic_series()
    increments = np.diff(truth_mm)
    for name, estimate in (("Perfect estimate", increments), ("Half-sized estimate", increments / 2)):
        metrics = increment_stats(increments, estimate)
        print(f"{name}: pairs={metrics['pairs']}, "
              f"slope={metrics['slope_through_origin']:.3f}, "
              f"RMSE={metrics['rmse_increments_mm']:.3f} mm")


if __name__ == "__main__":
    main()

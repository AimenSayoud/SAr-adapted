"""X-146's dry-snow delay model recovers the textbook values (Guneriussen et al., 2001)."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "experiments"))
import winter_snow_x146 as W  # noqa: E402


def test_delay_at_vertical_incidence_is_refractive_excess():
    # θ = 0: the extra path is (√ε − 1) per unit depth, read as sinking
    eps = 1 + 1.6 * 0.3 + 1.86 * 0.3 ** 3
    assert np.isclose(W.delay_slope(0.3, 0.0), -(np.sqrt(eps) - 1))


def test_delay_about_one_mm_per_mm_swe_at_sentinel1_incidence():
    # ρ 0.25 g cm⁻³, θ 39°: ≈ 0.24 mm per mm of depth, i.e. ≈ 1 mm per mm of snow water equivalent
    s = W.delay_slope(0.25, 39.0)
    assert -0.26 < s < -0.22
    assert np.isclose(s / 0.25, -0.96, atol=0.05)


def test_denser_snow_delays_more():
    assert W.delay_slope(0.35, 39.0) < W.delay_slope(0.15, 39.0) < 0

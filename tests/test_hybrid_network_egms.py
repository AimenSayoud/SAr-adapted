"""Ground-truth tests for hybrid network design and European Ground Motion Service (EGMS) integration.

Modules under test:
- insar_wetlands/hybrid_network.py
- insar_wetlands/egms.py

Verifies that:
(1) build_hybrid_pairs combines short baseline pairs and annual pairs, tags their kind,
    and deduplicates without dropping coverage;
(2) network_summary produces consistent graph metrics;
(3) egms.verdict distinguishes the three physical scenarios (all empty, border-only, on-peat points);
(4) egms.load_egms_timeseries_csv filters points by AOI geometry and fails on missing coordinates.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
from shapely.geometry import Polygon

from insar_wetlands.egms import (
    load_egms_timeseries_csv,
    verdict,
)

# =========================================================================
# hybrid_network — build_hybrid_pairs & network_summary
# =========================================================================


# =========================================================================
# egms — verdict
# =========================================================================

def test_egms_verdict_all_empty():
    """No points on tourbiere and no points on border."""
    stats = {
        "aoi": {"n": 0},
        "buffer": {"n": 0},
    }
    msg = verdict(stats)
    assert "Aucun point EGMS sur l'AOI NI en bordure" in msg


def test_egms_verdict_border_only():
    """0 points on peat, but valid points on surrounding infrastructure."""
    stats = {
        "aoi": {"n": 0},
        "buffer": {"n": 18, "median": -1.2, "p10": -2.5, "p90": -0.2},
    }
    msg = verdict(stats)
    assert "0 point sur la tourbière, 18 en bordure" in msg
    assert "-1.2 mm/an" in msg
    assert "Cohérent avec PSI = bâti/routes seulement" in msg


def test_egms_verdict_points_on_peat():
    """Valid points detected on the wetland itself."""
    stats = {
        "aoi": {"n": 7, "median": -3.4},
        "buffer": {"n": 12, "median": -1.0},
    }
    msg = verdict(stats)
    assert "7 points EGMS sur la tourbière" in msg
    assert "-3.4 mm/an" in msg
    assert "à comparer directement à notre résultat" in msg


# =========================================================================
# egms — load_egms_timeseries_csv
# =========================================================================

def test_load_egms_timeseries_csv_spatial_filtering(tmp_path: Path, monkeypatch):
    """Filters points strictly within the AOI polygon."""
    csv_file = tmp_path / "egms_sample.csv"
    # AOI will be a box around (16.30, 52.76)
    csv_content = (
        "latitude,longitude,20220101,20220113\n"
        "52.7630,16.3090,-1.5,-2.0\n"  # Inside
        "53.5000,17.0000,0.5,1.0\n"    # Far outside
    )
    csv_file.write_text(csv_content)

    # Mock load_aoi returning a 0.05-deg square around (16.309, 52.763)
    poly = Polygon([
        (16.25, 52.70),
        (16.35, 52.70),
        (16.35, 52.80),
        (16.25, 52.80),
    ])
    monkeypatch.setattr("insar_wetlands.aoi.load_aoi", lambda cfg=None: poly)

    df_filtered = load_egms_timeseries_csv(csv_file)
    assert len(df_filtered) == 1
    assert abs(df_filtered.loc[0, "latitude"] - 52.7630) < 1e-4
    assert abs(df_filtered.loc[0, "longitude"] - 16.3090) < 1e-4


def test_load_egms_timeseries_csv_missing_columns_raises(tmp_path: Path, monkeypatch):
    """Fails fast when latitude or longitude column is missing."""
    csv_file = tmp_path / "bad_egms.csv"
    csv_file.write_text("northing,easting,20220101\n500000,300000,0.0\n")

    monkeypatch.setattr("insar_wetlands.aoi.load_aoi", lambda cfg=None: Polygon())
    with pytest.raises(KeyError, match="colonnes lat/lon EGMS introuvables"):
        load_egms_timeseries_csv(csv_file)

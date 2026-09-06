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
from insar_wetlands.hybrid_network import (
    build_hybrid_pairs,
    network_summary,
)

# =========================================================================
# hybrid_network — build_hybrid_pairs & network_summary
# =========================================================================

def test_build_hybrid_pairs_combines_short_and_annual_without_duplicates():
    """Hybrid network must contain both short and annual pairs, properly tagged

    and strictly deduplicated.
    """
    # 3 years of synthetic 12-day Sentinel-1 dates covering April targets
    dates = pd.date_range("2021-03-15", "2023-05-15", freq="12D")
    dates_series = pd.Series(dates)

    hybrid = build_hybrid_pairs(
        dates_series,
        short_max_days=48,
        max_pairs_per_date=3,
        annual_target_month=4,
        annual_target_day=15,
        annual_window_days=20,
        max_gap_years=2,
    )

    assert isinstance(hybrid, pd.DataFrame)
    for col in ["ref_date", "sec_date", "dt_days", "pair", "kind"]:
        assert col in hybrid.columns

    # Check deduplication
    assert len(hybrid) == len(hybrid.drop_duplicates("pair"))

    # Both types must be present
    kinds = set(hybrid["kind"].unique())
    assert kinds == {"short", "annual"}

    # Short pairs must be within short_max_days
    short_pairs = hybrid[hybrid["kind"] == "short"]
    assert (short_pairs["dt_days"] <= 48).all()
    assert len(short_pairs) > 0

    # Annual pairs must span across years (~365 days or ~730 days)
    annual_pairs = hybrid[hybrid["kind"] == "annual"]
    assert (annual_pairs["dt_days"] >= 300).all()
    assert len(annual_pairs) > 0


def test_network_summary_computes_consistent_invariants():
    """network_summary accurately counts short/annual pairs and temporal span."""
    pairs_data = [
        {"pair": "20220101_20220113", "dt_days": 12, "kind": "short"},
        {"pair": "20220113_20220125", "dt_days": 12, "kind": "short"},
        {"pair": "20220101_20220125", "dt_days": 24, "kind": "short"},
        {"pair": "20220415_20230416", "dt_days": 366, "kind": "annual"},
    ]
    df = pd.DataFrame(pairs_data)

    summary = network_summary(df)

    assert summary["n_pairs_total"] == 4
    assert summary["n_short"] == 3
    assert summary["n_annual"] == 1
    assert summary["dt_days_min"] == 12
    assert summary["dt_days_max"] == 366
    assert summary["n_dates"] == 5  # 20220101, 20220113, 20220125, 20220415, 20230416
    assert summary["span_days"] == (pd.Timestamp("2023-04-16") - pd.Timestamp("2022-01-01")).days


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

"""The viewing geometry in config.yaml must be physically possible for the bursts used.

X-042: the descending track carried incidence 34.10° / heading 193.6° while its burst sees
the AOI at 38.99–39.36° — every two-geometry result inherited the error. These invariants
need no data and would have failed on that config.
"""
from pathlib import Path

import yaml

CFG = yaml.safe_load((Path(__file__).resolve().parents[1] / "config" / "config.yaml").read_text())
TRACKS = CFG["sentinel1"]["tracks"]

# Sentinel-1 IW sub-swath incidence envelopes (ESA user guide, rounded outward).
SUBSWATH_INCIDENCE = {"IW1": (29.0, 36.5), "IW2": (34.5, 42.5), "IW3": (39.5, 46.5)}
# Heading of a sun-synchronous orbit at ~52° N, clockwise from north.
HEADING = {"ASCENDING": (340.0, 356.0), "DESCENDING": (184.0, 200.0)}


def test_incidence_lies_inside_the_measured_aoi_range():
    for name, t in TRACKS.items():
        lo, hi = t["aoi_incidence_range_deg"]
        assert lo <= t["incidence_angle_deg"] <= hi, f"{name}: {t['incidence_angle_deg']} not in [{lo}, {hi}]"


def test_incidence_is_possible_for_the_subswath():
    for name, t in TRACKS.items():
        lo, hi = SUBSWATH_INCIDENCE[t["subswath"]]
        assert lo <= t["incidence_angle_deg"] <= hi, f"{name} {t['subswath']}: {t['incidence_angle_deg']}"
        assert t["burst_id"].endswith(t["subswath"])


def test_heading_matches_the_flight_direction():
    for name, t in TRACKS.items():
        lo, hi = HEADING[t["flight_direction"]]
        assert lo <= t["heading_deg"] <= hi, f"{name}: heading {t['heading_deg']}"

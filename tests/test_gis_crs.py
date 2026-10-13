"""Every GeoJSON shipped in data/gis has coordinates inside the valid area of the CRS it declares (C-042: the AOI
outline declared UTM 33N while holding lon/lat, so GIS tools placed it near the equator)."""
import json
from pathlib import Path

import pytest

GIS = Path(__file__).resolve().parents[1] / "data" / "gis"
# (x range, y range) each CRS can hold for a site in Poland
VALID = {"urn:ogc:def:crs:OGC:1.3:CRS84": ((14.0, 24.5), (49.0, 55.0)),
         "urn:ogc:def:crs:EPSG::4326": ((14.0, 24.5), (49.0, 55.0)),
         "urn:ogc:def:crs:EPSG::32633": ((166_000.0, 834_000.0), (5_400_000.0, 6_150_000.0))}


def _points(c):
    if isinstance(c[0], (int, float)):
        yield c
    else:
        for x in c:
            yield from _points(x)


@pytest.mark.parametrize("path", sorted(GIS.glob("*.geojson")), ids=lambda p: p.name)
def test_coordinates_match_declared_crs(path):
    d = json.loads(path.read_text(encoding="utf-8"))
    crs = d.get("crs", {}).get("properties", {}).get("name", "urn:ogc:def:crs:OGC:1.3:CRS84")
    assert crs in VALID, f"{path.name}: undeclared CRS {crs}"
    (x0, x1), (y0, y1) = VALID[crs]
    pts = [p for f in d["features"] if f.get("geometry") for p in _points(f["geometry"]["coordinates"])]
    assert pts
    assert all(x0 <= p[0] <= x1 and y0 <= p[1] <= y1 for p in pts), f"{path.name}: coordinates outside {crs}"

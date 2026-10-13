"""X-110 helpers: per-cell sums of a geocoded raster land in the right cube cells; the chain of a phase series returns
its wrapped consecutive differences in mm."""
import sys
from pathlib import Path

import numpy as np
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "experiments"))
import slc_phaselink_x110 as P  # noqa: E402


def test_cell_sums_of_a_constant_phase_raster(tmp_path):
    import rasterio
    from rasterio.transform import from_origin
    gx = 585720.0 + 40 * np.arange(4)
    gy = 5849440.0 - 40 * np.arange(3)
    st = xr.Dataset(coords={"x": gx, "y": gy})
    a = np.full((12, 32), np.exp(1j * 0.4), "complex64")            # 10 m rows, 5 m cols over the 3 x 4 cells
    with rasterio.open(tmp_path / "l.tif", "w", driver="GTiff", height=12, width=32, count=1, dtype="complex64",
                       crs="EPSG:32633", transform=from_origin(gx[0] - 20, gy[0] + 20, 5, 10)) as w:
        w.write(a, 1)
    s = P.cell_sums(tmp_path / "l.tif", st)
    assert s.shape == (3, 4) and np.allclose(np.angle(s), 0.4) and np.allclose(np.abs(s), 32)


def test_chain_from_phase_is_wrapped_differences():
    import p6_validation_x082 as X
    ph = np.array([0.0, 0.5, 1.2, -3.0])
    q = P.chain_from_phase(["20210101", "20210107", "20210113", "20210119"], ph)
    d = np.angle(np.exp(1j * np.diff(ph)))
    assert np.allclose(q.est_mm.to_numpy(), d * X.PHASE_TO_MM)


def test_linked_files_picks_dates_only(tmp_path):
    d = tmp_path / "t175_374052_iw1" / "linked_phase"
    (d / "20210101_20210326").mkdir(parents=True)
    for n in ("20210107.slc.tif", "20210101.slc.tif", "compressed_20210101_20210401_20210624.tif",
              "temporal_coherence_20210101_20210326.tif"):
        (d / n).write_bytes(b"")
    (d / "20210101_20210326" / "20210113.slc.tif").write_bytes(b"")
    assert [f.name for f in P.linked_files(tmp_path)] == ["20210101.slc.tif", "20210107.slc.tif"]

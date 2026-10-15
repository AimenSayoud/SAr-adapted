"""The crop reader leaves HyP3 crops as delivered by default and, with RZECIN_HYP3_RELOCATE=1, puts their content back
under the labels by the X-120 shift of the file's track (exact on a linear field)."""
import numpy as np
import rasterio
from rasterio.transform import from_origin

from insar_wetlands import placement as P
from insar_wetlands import stack


def _write(root, pair, burst, field):
    d = root / pair
    d.mkdir(parents=True)
    t = from_origin(585720.0, 5849280.0, 40.0, 40.0)
    with rasterio.open(d / f"S1_{burst}_IW1_{pair}_VV_INT40_0000_corr.tif", "w", driver="GTiff", height=field.shape[0],
                       width=field.shape[1], count=1, dtype="float32", crs="EPSG:32633", transform=t) as f:
        f.write(field.astype("float32"), 1)
    (d / "ok").touch()          # a finished crop, as products.download_and_crop marks it


def test_reader_relocates_only_when_asked(tmp_path, monkeypatch):
    ny, nx = 30, 40
    x = 585740.0 + 40 * np.arange(nx)
    y = 5849260.0 - 40 * np.arange(ny)
    X, Y = np.meshgrid(x, y)
    dx, dy = P.HYP3_CONTENT_SHIFT_M[40]["ascending"]
    truth = 1e-3 * (X - 585000) + 2e-3 * (Y - 5848000)
    delivered = 1e-3 * (X + dx - 585000) + 2e-3 * (Y + dy - 5848000)    # content of ground at L + (dx, dy)
    _write(tmp_path, "20210101_20210107", "374052", delivered)
    monkeypatch.delenv(stack.RELOCATE_ENV, raising=False)
    a = stack.load_layer(tmp_path, "corr").isel(pair=0).values
    assert np.allclose(a, delivered)
    monkeypatch.setenv(stack.RELOCATE_ENV, "1")
    b = stack.load_layer(tmp_path, "corr").isel(pair=0).values
    inner = np.s_[2:-2, 2:-2]
    assert np.allclose(b[inner], truth[inner], atol=1e-6)


def test_track_from_file_name():
    assert stack.hyp3_track_of("S1_374052_IW1_20210101_20210107_VV_INT40_0000_corr.tif") == "ascending"
    assert stack.hyp3_track_of("S1_045867_IW2_20210101_20210107_VV_INT40_0000_corr.tif") == "descending"

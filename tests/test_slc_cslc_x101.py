"""X-101 helpers that run on Colab only, checked here: the burst id COMPASS looks up, and the AOI burst database
answering COMPASS's own query (compass.utils.helpers.burst_bbox_from_db), so only the AOI is geocoded."""
import sqlite3
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "experiments"))
import slc_cslc_x101 as M  # noqa: E402

COMPASS_QUERY = "SELECT epsg, xmin, ymin, xmax, ymax FROM burst_id_map WHERE burst_id_jpl = ?"


def test_jpl_id_matches_s1reader_form():
    assert M.jpl_id("175_374052_IW1") == "t175_374052_iw1"
    assert M.jpl_id("022_045867_IW2") == "t022_045867_iw2"
    assert M.jpl_id("22_045867_IW2") == "t022_045867_iw2"


def test_burst_db_answers_compass_query(tmp_path):
    box = (585000.0, 5843000.0, 592000.0, 5850000.0)
    db = M.burst_db(tmp_path / "b.sqlite3", ["t175_374052_iw1", "t022_045867_iw2"], box)
    con = sqlite3.connect(db)
    con.row_factory = sqlite3.Row
    r = con.execute(COMPASS_QUERY, ("t022_045867_iw2",)).fetchone()
    assert r["epsg"] == M.EPSG and (r["xmin"], r["ymin"], r["xmax"], r["ymax"]) == box
    rows = con.execute("SELECT * FROM burst_id_map WHERE burst_id_jpl IN (?, ?)",
                       ("t175_374052_iw1", "t022_045867_iw2")).fetchall()
    assert len(rows) == 2


def test_dates_on_drive_ignore_static_layers(tmp_path):
    for n in ("t175_374052_iw1_20210107.h5", "t175_374052_iw1_20210113_v1.h5", "static_layers_t175_374052_iw1.h5"):
        (tmp_path / n).write_bytes(b"")
    assert M.done_dates(tmp_path) == {"20210107", "20210113"}


def test_launcher_installs_everything_with_conda():
    from insar_wetlands.colab import launcher
    nb = launcher("slc_cslc_x101", conda="compass isce3 s1reader burst2safe sentineleof sardem gdal",
                  secrets="EARTHDATA_USERNAME EARTHDATA_PASSWORD")
    src = "\n".join("".join(c["source"]) for c in nb["cells"])
    assert "condacolab.install()" in src
    assert "mamba install" in src and "burst2safe" in src and "xarray" in src
    assert "bash colab/setup.sh conda" in src and "pip install -q burst2safe" not in src
    assert "os.environ[k] = userdata.get(k)" in src                 # secrets reach the subprocess as env vars
    assert src.index("userdata.get") < src.index("register run")


def test_keep_puts_the_date_in_the_name(tmp_path):
    work, dest = tmp_path / "w", tmp_path / "d"
    (work / "t175_374052_iw1" / "20210107").mkdir(parents=True)
    (work / "t175_374052_iw1" / "20210107" / "t175_374052_iw1.h5").write_bytes(b"x")
    (work / "t175_374052_iw1" / "20210113").mkdir(parents=True)
    (work / "t175_374052_iw1" / "20210113" / "t175_374052_iw1_20210113.h5").write_bytes(b"x")
    (work / "scratch").mkdir()
    (work / "scratch" / "tmp_20210107.h5").write_bytes(b"x")
    dest.mkdir()
    rows = M.keep(work, dest, "ascending", True)
    assert M.done_dates(dest) == {"20210107", "20210113"}
    assert len(rows) == 2 and not (dest / "tmp_20210107.h5").exists()


def _tiny_dem(path):
    import numpy as np
    import rasterio
    from rasterio.transform import from_origin
    with rasterio.open(path, "w", driver="GTiff", height=3, width=4, count=1, dtype="float32", crs="EPSG:4326",
                       transform=from_origin(16.30, 52.77, 0.001, 0.001)) as w:
        w.write(np.full((3, 4), 94.0, "float32"), 1)


def test_dem_gets_the_geoid_added(tmp_path, monkeypatch):
    import numpy as np
    import rasterio
    _tiny_dem(tmp_path / "egm.tif")
    monkeypatch.setattr(M, "geoid_undulation", lambda lon, lat: np.full(np.shape(lon), 34.5))
    M.to_ellipsoid(tmp_path / "egm.tif", tmp_path / "ell.tif")
    with rasterio.open(tmp_path / "ell.tif") as r:
        assert np.allclose(r.read(1), 128.5)
        assert r.crs.to_epsg() == 4326          # plain 2-D lon/lat, not the compound geoid CRS


def test_missing_geoid_grid_is_refused(tmp_path, monkeypatch):
    import numpy as np
    import pytest
    _tiny_dem(tmp_path / "egm.tif")
    monkeypatch.setattr(M, "geoid_undulation", lambda lon, lat: np.zeros(np.shape(lon)))   # PROJ without the grid
    with pytest.raises(RuntimeError, match="no EGM2008 grid"):
        M.to_ellipsoid(tmp_path / "egm.tif", tmp_path / "ell.tif")


def test_runconfig_polarization_list_becomes_string(tmp_path):
    import yaml
    f = tmp_path / "geo_runconfig_20210101_t175_374052_iw1.yaml"
    f.write_text(yaml.safe_dump({"runconfig": {"name": "cslc_s1_workflow_default",
                                               "groups": {"processing": {"polarization": ["co-pol"], "x": 1}}}}))
    M.fix_runconfig(f)
    d = yaml.safe_load(f.read_text())
    assert d["runconfig"]["groups"]["processing"] == {"polarization": "co-pol", "x": 1}
    assert d["runconfig"]["groups"]["output"]["chunk_size"] == [M.CHUNK, M.CHUNK]
    assert d["runconfig"]["name"] == "cslc_s1_workflow_default"


def test_dem_check_wants_the_copernicus_heights(tmp_path):
    """The chain geolocates correctly only on the delivered (orthometric) heights: equal passes, +N is refused."""
    import numpy as np
    import pytest
    import rasterio
    import xarray as xr
    from rasterio.transform import from_origin
    from insar_wetlands.cube import cube_root
    st = xr.open_dataset(cube_root() / "silver" / "static_40m.nc")
    ortho = float(np.nanmedian(st.dem_glo30_m.where(st.zone == 1)))
    for name, val, ok in (("ortho.tif", ortho, True), ("ell.tif", ortho + 34.5, False)):
        with rasterio.open(tmp_path / name, "w", driver="GTiff", height=200, width=200, count=1, dtype="float32",
                           crs="EPSG:4326", transform=from_origin(16.25, 52.80, 0.001, 0.0005)) as w:
            w.write(np.full((200, 200), val, "float32"), 1)
        if ok:
            assert abs(M.check_dem(tmp_path / name)) < M.DEM_TOL_M
        else:
            with pytest.raises(RuntimeError, match="does not hold the Copernicus heights"):
                M.check_dem(tmp_path / name)


def test_vh_goes_to_its_own_folder(monkeypatch, tmp_path):
    import insar_wetlands.paths as paths
    monkeypatch.setattr(paths, "resolve_drive", lambda *a, **k: tmp_path)
    assert M.drive_dest("ascending", "VV") == tmp_path / "slc_pilot" / "2021" / "cslc" / "ascending"
    assert M.drive_dest("ascending", "VH") == tmp_path / "slc_pilot" / "2021" / "cslc_vh" / "ascending"


def test_build_safes_passes_a_dates_granules_together(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(M, "sh", lambda cmd, **kw: calls.append([str(c) for c in cmd]))
    M.build_safes([["S1_374052_IW1_20210101T163611_VV_3243-BURST", "S1_374052_IW1_20210101T163611_VH_3243-BURST"]],
                  tmp_path)
    assert len(calls) == 1
    cmd = calls[0]
    assert "burst2safe" in cmd and cmd.count("--all-anns") == 1
    assert any("_VV_" in c for c in cmd) and any("_VH_" in c for c in cmd)

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

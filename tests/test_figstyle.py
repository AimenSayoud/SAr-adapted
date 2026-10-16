"""figstyle: a figure drawn at a journal width is exported at exactly that width, with its text at least 7 pt and its
fonts embedded as TrueType, and with a provenance record beside it."""
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pytest  # noqa: E402

from insar_wetlands import figstyle as fs  # noqa: E402


def test_size_is_the_publisher_width():
    assert fs.size("single", 60)[0] * 25.4 == pytest.approx(90.0)
    assert fs.size("full", 60)[0] * 25.4 == pytest.approx(190.0)
    assert fs.size(120, 45)[1] * 25.4 == pytest.approx(45.0)


def test_text_is_never_below_the_minimum():
    with fs.style():
        for k in ("font.size", "axes.labelsize", "xtick.labelsize", "ytick.labelsize", "legend.fontsize"):
            assert float(matplotlib.rcParams[k]) >= fs.MIN_FONT_PT


def test_export_keeps_the_page_size_and_embeds_truetype(tmp_path):
    pypdf = pytest.importorskip("pypdf")
    with fs.style():
        fig, ax = plt.subplots(figsize=fs.size("single", 60))
        ax.plot([0, 1], [0, 1])
        ax.set_xlabel("x (unit)")
        rec = fs.export(fig, tmp_path / "F0", provenance={"data": "synthetic"})
    plt.close(fig)
    page = pypdf.PdfReader(tmp_path / "F0.pdf").pages[0]
    assert float(page.mediabox.width) / 72 * 25.4 == pytest.approx(90.0, abs=0.1)     # not cropped
    assert float(page.mediabox.height) / 72 * 25.4 == pytest.approx(60.0, abs=0.1)
    fonts = {f["/Subtype"] for f in page["/Resources"]["/Font"].get_object().values()
             for f in [f.get_object()]}
    assert fonts <= {"/TrueType", "/Type0"}                                            # no Type 3 bitmaps
    assert (tmp_path / "F0.png").exists()
    meta = json.loads((tmp_path / "F0.json").read_text())
    assert meta["width_mm"] == 90.0 and meta["data"] == "synthetic" and rec["figure"] == "F0"

"""Tests for the manuscript assembly pipeline.

Checks that: (1) sections are concatenated in MANUSCRIPT order, not alphabetical
or filesystem order; (2) missing files are skipped without crashing but are
visible in the report; (3) unresolved image references are REPORTED by name
rather than silently dropped — a manuscript that builds while quietly losing a
figure is worse than one that fails; (4) the real docs/paper/ tree assembles and
all its image references resolve.

Run: python tests/test_paper_build.py
"""

import tempfile
import zipfile
from pathlib import Path

from insar_wetlands.paper_build import (SECTION_ORDER, add_page_setup,
                                        polish_tables,
                                        assemble_markdown, build_data_appendix,
                                        collect_sections, csv_to_markdown,
                                        find_missing_images,
                                        repair_content_types)

CT = ('<?xml version="1.0" encoding="UTF-8"?><Types xmlns="http://schemas.'
      'openxmlformats.org/package/2006/content-types"><Default Extension="xml"'
      ' ContentType="application/xml"/></Types>')


def _fake_docx(path: Path, body: str, media=("word/media/rId1.png",),
               content_types: str = CT):
    doc = ('<?xml version="1.0" encoding="UTF-8"?><w:document xmlns:w="http://'
           'schemas.openxmlformats.org/wordprocessingml/2006/main"><w:body>'
           f"{body}</w:body></w:document>")
    with zipfile.ZipFile(path, "w") as z:
        z.writestr("[Content_Types].xml", content_types)
        z.writestr("word/document.xml", doc)
        for m in media:
            z.writestr(m, b"\x89PNG\r\n\x1a\n")
    return path


def _fake_paper(tmp: Path):
    (tmp / "figures").mkdir(parents=True, exist_ok=True)
    (tmp / "figures" / "F01_x.png").write_bytes(b"\x89PNG\r\n\x1a\n" + b"0" * 50)
    # deliberately out of alphabetical order relative to manuscript order
    (tmp / "04_results.md").write_text(
        "## 4. Results\n\n![Fig 1](figures/F01_x.png)\n\n![Gone](figures/NOPE.png)\n")
    (tmp / "00_title_abstract.md").write_text("# Title\n\n## Abstract\n\nText.\n")
    (tmp / "01_introduction.md").write_text("## 1. Introduction\n\nText.\n")
    return tmp


def test_sections_are_ordered_by_manuscript_not_filesystem():
    with tempfile.TemporaryDirectory() as d:
        tmp = _fake_paper(Path(d))
        got = [f.name for f in collect_sections(tmp)]
        assert got == ["00_title_abstract.md", "01_introduction.md",
                       "04_results.md"], got
        # absent sections are skipped, not fatal
        assert "05_discussion.md" not in got


def test_assembly_keeps_order_and_reports_missing_images():
    with tempfile.TemporaryDirectory() as d:
        tmp = _fake_paper(Path(d))
        rep = assemble_markdown(tmp, tmp / "_m.md", title="T")
        text = (tmp / "_m.md").read_text()
        assert rep["n_sections"] == 3
        assert text.index("Abstract") < text.index("Introduction") < text.index("Results")
        # the resolvable image is kept, the broken one is NAMED
        assert "figures/F01_x.png" in rep["images"]
        assert rep["missing_images"] == ["figures/NOPE.png"], rep["missing_images"]


def test_find_missing_images_ignores_urls():
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        md = "![a](https://example.com/x.png)\n![b](figures/absent.png)\n"
        assert find_missing_images(md, tmp) == ["figures/absent.png"]


def test_empty_directory_raises():
    with tempfile.TemporaryDirectory() as d:
        try:
            assemble_markdown(Path(d), Path(d) / "out.md")
        except FileNotFoundError:
            return
        raise AssertionError("an empty paper directory must raise")


def test_real_paper_tree_assembles():
    """The actual manuscript must assemble, and its figure references must be
    consistent with whatever has been generated.

    Two distinct states, deliberately treated differently:
      - `figures/` empty  -> figures not exported yet. Pending, not a defect.
      - `figures/` populated but a reference does not resolve -> BROKEN link,
        which must fail: a manuscript that builds while quietly dropping a
        figure is worse than one that fails."""
    paper = Path(__file__).resolve().parents[1] / "docs" / "paper"
    if not paper.exists():
        print("  (docs/paper absent -> skipped)")
        return
    with tempfile.TemporaryDirectory() as d:
        rep = assemble_markdown(paper, Path(d) / "_m.md")
        print(f"  {rep['n_sections']} sections, {rep['chars']} chars, "
              f"{len(rep['images'])} image refs")
        assert rep["n_sections"] >= 5, rep["files"]
        assert rep["chars"] > 20000, "manuscript looks truncated"
        assert rep["images"], "no figure is referenced in the manuscript"

        generated = list((paper / "figures").glob("*.png")) \
            if (paper / "figures").exists() else []
        if not generated:
            print(f"  ({len(rep['missing_images'])} figures pending export "
                  f"-> run export_figures_en.ipynb)")
        else:
            assert not rep["missing_images"], (
                f"figures/ is populated ({len(generated)} PNG) but these "
                f"references are broken: {rep['missing_images']}")


def test_csv_to_markdown_escapes_pipes_and_flags_truncation():
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "T99_x.csv"
        p.write_text("zone,note\nA,a|b\n" + "".join(f"Z{i},v\n" for i in range(50)))
        md = csv_to_markdown(p, max_rows=10)
        rows = [ln for ln in md.splitlines() if ln.startswith("|")]
        # header + separator + exactly max_rows body rows, no more
        assert len(rows) == 12, len(rows)
        # a literal pipe inside a cell must not create a phantom column
        assert r"a\|b" in md
        # the reader must be told the table was cut
        assert "41 further rows omitted" in md, md[-200:]


def test_data_appendix_is_regenerated_from_the_csvs():
    with tempfile.TemporaryDirectory() as d:
        figs = Path(d) / "figures"
        figs.mkdir()
        (figs / "T01_zones.csv").write_text("zone,px\nA,499\nC,398\n")
        (figs / "T02_temporal_coherence.csv").write_text("zone,tcoh\nA,0.58\n")
        out = Path(d) / "09_appendix_data.md"
        rep = build_data_appendix(figs, out)
        text = out.read_text()
        assert rep["n_tables"] == 2, rep
        # values come from the CSV, so the document cannot drift from the data
        assert "499" in text and "0.58" in text
        # known ids get their real caption, so the appendix reads as a document
        assert "Zone definition, pixel counts and areas" in text
        assert text.index("T01") < text.index("T02")


def test_data_appendix_without_csvs_says_so():
    """An empty appendix must announce itself, not look like a finished one."""
    with tempfile.TemporaryDirectory() as d:
        figs = Path(d) / "figures"
        figs.mkdir()
        out = Path(d) / "09_appendix_data.md"
        rep = build_data_appendix(figs, out)
        assert rep["n_tables"] == 0
        assert "No table has been exported yet" in out.read_text()


def test_data_appendix_is_in_the_section_order():
    """Generated or not, the appendix must have a defined place in the
    manuscript; collect_sections skips it while it is absent."""
    assert SECTION_ORDER[-1] == "09_appendix_data.md"
    with tempfile.TemporaryDirectory() as d:
        tmp = _fake_paper(Path(d))
        assert "09_appendix_data.md" not in [f.name for f in collect_sections(tmp)]
        (tmp / "09_appendix_data.md").write_text("## Appendix B\n\nx\n")
        got = [f.name for f in collect_sections(tmp)]
        assert got[-1] == "09_appendix_data.md", got


def test_repair_content_types_declares_embedded_media():
    """pandoc 3.1.x embeds images without declaring their extension, which makes
    the OPC package invalid. The repair must add it, and must be idempotent."""
    with tempfile.TemporaryDirectory() as d:
        p = _fake_docx(Path(d) / "a.docx", "<w:p/>",
                       media=("word/media/rId1.png", "word/media/rId2.jpeg"))
        added = repair_content_types(p)
        assert set(added) == {"png", "jpeg"}, added
        ct = zipfile.ZipFile(p).read("[Content_Types].xml").decode()
        assert 'Extension="png" ContentType="image/png"' in ct
        assert 'Extension="jpeg" ContentType="image/jpeg"' in ct
        # the pre-existing declaration must survive
        assert 'Extension="xml"' in ct
        # running twice must change nothing
        assert repair_content_types(p) == []


def test_repair_content_types_noop_without_media():
    with tempfile.TemporaryDirectory() as d:
        p = _fake_docx(Path(d) / "a.docx", "<w:p/>", media=())
        assert repair_content_types(p) == []


def test_add_page_setup_replaces_the_empty_pandoc_section():
    """pandoc emits a self-closing <w:sectPr/>: it must be replaced wholesale,
    with lnNumType in its schema-mandated position (after pgMar, before cols)."""
    with tempfile.TemporaryDirectory() as d:
        p = _fake_docx(Path(d) / "a.docx", "<w:p/><w:sectPr />")
        assert add_page_setup(p, line_numbers=True) is True
        doc = zipfile.ZipFile(p).read("word/document.xml").decode()
        assert "<w:sectPr />" not in doc
        assert doc.index("<w:pgMar") < doc.index("<w:lnNumType") \
            < doc.index("<w:cols")
        # sectPr must remain the last child of body
        assert doc.index("</w:sectPr>") < doc.index("</w:body>")
        assert add_page_setup(p, line_numbers=True) is False  # idempotent


def test_add_page_setup_appends_when_no_section_exists():
    with tempfile.TemporaryDirectory() as d:
        p = _fake_docx(Path(d) / "a.docx", "<w:p/>")
        assert add_page_setup(p, line_numbers=True) is True
        doc = zipfile.ZipFile(p).read("word/document.xml").decode()
        assert doc.index("<w:sectPr>") > doc.index("<w:p/>")
        assert "<w:lnNumType" in doc


def test_add_page_setup_keeps_a_populated_section():
    """A real section must not be overwritten — only line numbering added."""
    with tempfile.TemporaryDirectory() as d:
        p = _fake_docx(Path(d) / "a.docx",
                       '<w:sectPr><w:pgSz w:w="99" w:h="88"/>'
                       '<w:cols w:space="720"/></w:sectPr>')
        assert add_page_setup(p, line_numbers=True) is True
        doc = zipfile.ZipFile(p).read("word/document.xml").decode()
        assert 'w:w="99"' in doc, "existing page size must be preserved"
        assert doc.index("<w:lnNumType") < doc.index("<w:cols")


def test_add_page_setup_can_be_declined():
    with tempfile.TemporaryDirectory() as d:
        p = _fake_docx(Path(d) / "a.docx",
                       '<w:sectPr><w:cols w:space="720"/></w:sectPr>')
        assert add_page_setup(p, line_numbers=False) is False
        assert "<w:lnNumType" not in \
            zipfile.ZipFile(p).read("word/document.xml").decode()


def test_line_numbers_are_off_unless_asked_for():
    """Continuous numbering is what a journal wants at submission, but it
    clutters a document being read. It must be opt-in."""
    import inspect
    from insar_wetlands.paper_build import build_manuscript
    for fn in (build_manuscript, add_page_setup):
        assert inspect.signature(fn).parameters["line_numbers"].default is False, fn
    with tempfile.TemporaryDirectory() as d:
        p = _fake_docx(Path(d) / "a.docx", "<w:p/><w:sectPr />")
        add_page_setup(p)                      # default: no numbering
        doc = zipfile.ZipFile(p).read("word/document.xml").decode()
        assert "<w:lnNumType" not in doc
        assert "<w:pgSz" in doc, "page geometry must still be set"


def test_polish_tables_forces_a_consistent_width():
    """pandoc sizes each table to its content, so neighbouring tables come out
    different widths and none lines up with the text column."""
    with tempfile.TemporaryDirectory() as d:
        body = ('<w:tbl><w:tblPr><w:tblW w:type="auto" w:w="0" /></w:tblPr></w:tbl>'
                '<w:tbl><w:tblPr><w:tblW w:type="auto" w:w="0"/></w:tblPr></w:tbl>')
        p = _fake_docx(Path(d) / "a.docx", body)
        assert polish_tables(p) == 2
        doc = zipfile.ZipFile(p).read("word/document.xml").decode()
        assert doc.count('w:type="pct" w:w="5000"') == 2, doc
        assert 'w:type="auto"' not in doc
        # both spellings pandoc emits must be handled, not just the spaced one


def test_polish_tables_is_a_noop_without_tables():
    with tempfile.TemporaryDirectory() as d:
        p = _fake_docx(Path(d) / "a.docx", "<w:p/>")
        assert polish_tables(p) == 0


if __name__ == "__main__":
    test_sections_are_ordered_by_manuscript_not_filesystem()
    test_assembly_keeps_order_and_reports_missing_images()
    test_find_missing_images_ignores_urls()
    test_empty_directory_raises()
    test_csv_to_markdown_escapes_pipes_and_flags_truncation()
    test_data_appendix_is_regenerated_from_the_csvs()
    test_data_appendix_without_csvs_says_so()
    test_data_appendix_is_in_the_section_order()
    test_repair_content_types_declares_embedded_media()
    test_repair_content_types_noop_without_media()
    test_add_page_setup_replaces_the_empty_pandoc_section()
    test_add_page_setup_appends_when_no_section_exists()
    test_add_page_setup_keeps_a_populated_section()
    test_add_page_setup_can_be_declined()
    test_line_numbers_are_off_unless_asked_for()
    test_polish_tables_forces_a_consistent_width()
    test_polish_tables_is_a_noop_without_tables()
    test_real_paper_tree_assembles()
    print("ALL PAPER-BUILD TESTS PASSED")

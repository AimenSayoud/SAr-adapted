"""Tests for the manuscript/data transcription check.

A referee found seven quantities where the running text disagreed with the
exported tables. The code was correct in every case; only the hand-typed
transcription was stale. These tests pin the guard that makes that class of
error visible at build time.

Run: python tests/test_paper_numbers.py
"""

import tempfile
from pathlib import Path

from insar_wetlands.paper_numbers import (
    GENERATED_SECTIONS,
    REGISTRY,
    SUPERSEDED,
    check_manuscript_numbers,
    expected_values,
    fixed,
    hand_written_text,
    scientific,
)


def test_formatters_match_manuscript_typography():
    # U+2212 MINUS SIGN, not an ASCII hyphen: the manuscript uses the former,
    # so a check built on the latter would never match.
    assert fixed(3)(-0.08085) == "−0.081"
    assert "-" not in fixed(3)(-0.08085)
    assert fixed(0, percent=True)(0.8932) == "89 %"
    assert fixed(1, percent=True)(0.647) == "64.7 %"
    assert fixed(2)(3.286) == "3.29"
    assert scientific(2)(4.838844e-46) == "4.84 × 10⁻⁴⁶"
    assert scientific(2)(2.2e-49) == "2.20 × 10⁻⁴⁹"


def _paper(tmp: Path, csv_rows: str, text: str):
    (tmp / "figures").mkdir(parents=True, exist_ok=True)
    (tmp / "figures" / "T07_seasonal_amplitudes.csv").write_text(csv_rows)
    (tmp / "04_results.md").write_text(text)
    return tmp


def test_a_stale_number_is_reported():
    with tempfile.TemporaryDirectory() as d:
        tmp = _paper(Path(d), "series,amplitude_mm\nA−C,3.286\n",
                     "The amplitude is 9.99 mm.\n")
        bad = check_manuscript_numbers(tmp)
        names = [b["name"] for b in bad]
        assert "seasonal amplitude A - C" in names, names
        assert any(b["expected"] == "3.29" for b in bad), bad


def test_a_correct_number_passes():
    with tempfile.TemporaryDirectory() as d:
        tmp = _paper(Path(d), "series,amplitude_mm\nA−C,3.286\n",
                     "The amplitude is 3.29 mm.\n")
        assert [b["name"] for b in check_manuscript_numbers(tmp)
                if b["name"] == "seasonal amplitude A - C"] == []


def test_missing_csv_is_pending_not_a_failure():
    """Before the first export there is nothing to check against. That is a
    pending state and must not be reported as a stale number."""
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        (tmp / "figures").mkdir()
        (tmp / "04_results.md").write_text("no numbers here\n")
        assert check_manuscript_numbers(tmp) == []
        vals = expected_values(tmp / "figures")
        assert all(v["expected"] is None for v in vals)
        assert all("not exported" in v.get("note", "") for v in vals)


def test_generated_appendix_is_excluded():
    """Appendix B is generated FROM the CSVs, so it agrees with them by
    construction. Including it would make the whole check vacuous."""
    assert "09_appendix_data.md" in GENERATED_SECTIONS
    with tempfile.TemporaryDirectory() as d:
        tmp = _paper(Path(d), "series,amplitude_mm\nA−C,3.286\n",
                     "The amplitude is 9.99 mm.\n")
        # the generated appendix carries the true value ...
        (tmp / "09_appendix_data.md").write_text("| A−C | 3.29 |\n")
        assert "3.29" not in hand_written_text(tmp), "appendix leaked into the check"
        # ... and must NOT rescue the stale prose
        names = [b["name"] for b in check_manuscript_numbers(tmp)]
        assert "seasonal amplitude A - C" in names, names


def test_registry_entries_are_wellformed():
    seen = set()
    for name, csv_name, column, where, fmt in REGISTRY:
        assert name not in seen, f"duplicate registry name {name!r}"
        seen.add(name)
        assert csv_name.endswith(".csv"), csv_name
        assert callable(fmt), name
        assert where is None or isinstance(where, dict), name


def test_real_manuscript_matches_its_exported_data():
    """The actual manuscript. If figures/ has not been exported the registry
    resolves to nothing and this is vacuously true, which is the correct
    pending behaviour."""
    paper = Path(__file__).resolve().parents[1] / "docs" / "paper"
    if not paper.exists():
        print("  (docs/paper absent -> skipped)")
        return
    from insar_wetlands.paper_numbers import format_report
    bad = check_manuscript_numbers(paper)
    resolved = [v for v in expected_values(paper / "figures")
                if v["expected"] is not None]
    print(f"  {len(resolved)}/{len(REGISTRY)} registered numbers resolved "
          f"against exported CSVs")
    assert not bad, "\n" + format_report(bad)


def test_a_superseded_value_is_caught_even_when_partially_updated():
    """The presence check alone proves a number appears somewhere, not that
    every occurrence was updated. These figures repeat across the abstract,
    results, conclusions and appendix, so fixing four of five sites still
    passes. Superseded values are therefore checked for absence."""
    with tempfile.TemporaryDirectory() as d:
        tmp = _paper(Path(d), "series,amplitude_mm\nA−C,3.286\n",
                     "The amplitude is 3.29 mm.\n")
        # correct value present -> clean
        assert check_manuscript_numbers(tmp) == []
        # a superseded figure left behind in one place -> caught
        (tmp / "01_introduction.md").write_text(
            "Earlier we reported Δ = −0.069 here.\n")
        bad = check_manuscript_numbers(tmp)
        assert len(bad) == 1, bad
        assert bad[0]["superseded"] == "−0.069", bad


def test_superseded_list_does_not_collide_with_current_values():
    """A superseded string must never be a substring of a value the registry
    currently expects, or the manuscript could never pass."""
    paper = Path(__file__).resolve().parents[1] / "docs" / "paper"
    if not paper.exists():
        return
    current = {v["expected"] for v in expected_values(paper / "figures")
               if v["expected"]}
    for old in SUPERSEDED:
        assert not any(old in c for c in current), \
            f"superseded {old!r} collides with a current value"


if __name__ == "__main__":
    test_formatters_match_manuscript_typography()
    test_a_stale_number_is_reported()
    test_a_correct_number_passes()
    test_missing_csv_is_pending_not_a_failure()
    test_generated_appendix_is_excluded()
    test_registry_entries_are_wellformed()
    test_a_superseded_value_is_caught_even_when_partially_updated()
    test_superseded_list_does_not_collide_with_current_values()
    test_real_manuscript_matches_its_exported_data()
    print("ALL PAPER-NUMBER TESTS PASSED")

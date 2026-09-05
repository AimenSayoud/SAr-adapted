"""Check that the numbers written in the manuscript match the exported data.

Why this module exists. A referee found seven quantities where the running text
and Appendix B disagreed — a stale Wilcoxon *p*, a paired Δ that matched neither
the mean nor the median, a jackknife range that did not overlap the exported
one. Every one came from the same cause: a number was computed once, typed into
the prose by hand, and never revisited when the pipeline was re-run. Unit tests
did not catch it because the code was right; only the *transcription* was stale.

So the check is deliberately shaped around transcription, not computation. Each
entry below names a quantity, says where its truth lives in the exported CSVs,
and says how it is written in the text. The check then asserts that the string
actually appears in the assembled manuscript. It cannot verify prose, and it
does not try; it verifies the handful of load-bearing figures that a referee
will cross-check first, which is exactly the set that went stale.

`REGISTRY` is meant to grow. Adding a number here costs one line and buys
permanent protection against it drifting.
"""

from __future__ import annotations

import csv
import math
import re
from pathlib import Path

# Characters the manuscript uses that differ from what Python prints.
_MINUS = "−"          # U+2212 MINUS SIGN, not ASCII hyphen
_SUPERSCRIPT = str.maketrans("0123456789-", "⁰¹²³⁴"
                                            "⁵⁶⁷⁸⁹"
                                            "⁻")


def _read_csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def _cell(rows: list[dict], column: str, where: dict | None = None) -> float:
    """One value from a CSV: `column`, optionally on the row matching `where`."""
    if where:
        for r in rows:
            if all(str(r.get(k, "")).strip() == str(v) for k, v in where.items()):
                return float(r[column])
        raise KeyError(f"no row matching {where}")
    if len(rows) != 1:
        raise KeyError(f"{len(rows)} rows but no row selector given")
    return float(rows[0][column])


# ---------------------------------------------------------------- formatters
def fixed(dp: int, percent: bool = False, signed: bool = False):
    """Plain decimal, in the manuscript's typography (U+2212 for minus)."""
    def fmt(v: float) -> str:
        x = v * 100 if percent else v
        s = f"{abs(x):.{dp}f}"
        if x < 0:
            s = _MINUS + s
        elif signed:
            s = "+" + s
        return s + (" %" if percent else "")
    return fmt


def scientific(dp: int = 2):
    """e.g. 4.84e-46 -> '4.84 × 10⁻⁴⁶', matching the manuscript's rendering."""
    def fmt(v: float) -> str:
        exp = int(math.floor(math.log10(abs(v))))
        mant = v / (10 ** exp)
        return f"{mant:.{dp}f} × 10{str(exp).translate(_SUPERSCRIPT)}"
    return fmt


# Values that were published once and are now wrong. Checked for ABSENCE.
#
# The presence check alone proves a number appears somewhere, not that every
# occurrence was updated: these figures are repeated across the abstract,
# results, conclusions and appendix, so correcting four of five sites still
# passes. Listing the superseded values catches exactly that partial update,
# and also catches a copy-paste from the previous draft.
SUPERSEDED = {
    "2.2 × 10⁻⁴⁹": "old Wilcoxon p, superseded by 4.84 × 10⁻⁴⁶",
    "−0.069": "old paired Δ; the mean is −0.081 and the median −0.050",
    "−0.0705": "old date-jackknife minimum, superseded by −0.0842",
    "−0.0652": "old date-jackknife maximum, superseded by −0.0774",
    "23.1 %": "old zone-D usable fraction, superseded by 23.2 %",
    "+0.077)": "grassland freeze gain truncated instead of rounded; 0.07775 -> 0.078",
}


# ------------------------------------------------------------------ registry
# (name, csv, column, row selector, formatter)
REGISTRY = [
    ("Wilcoxon p, A vs C", "T05_paired_test.csv", "wilcoxon_p", None,
     scientific(2)),
    ("paired delta mean, A - C", "T05_paired_test.csv", "delta_mean", None,
     fixed(3)),
    ("date-jackknife minimum", "T05_paired_test.csv", "date_jackknife_min", None,
     fixed(4)),
    ("date-jackknife maximum", "T05_paired_test.csv", "date_jackknife_max", None,
     fixed(4)),
    ("fraction of pairs with A lower", "T05_paired_test.csv", "frac_a_lower",
     None, fixed(0, percent=True)),
    # Usable fractions come from the multi-threshold sweep, which is the one
    # code path that drops invalid pixels correctly; T02 is derived from it.
    ("usable fraction A at 0.7", "T03_multithreshold.csv", "A",
     {"threshold": "0.7"}, fixed(1, percent=True)),
    ("usable fraction C at 0.7", "T03_multithreshold.csv", "C",
     {"threshold": "0.7"}, fixed(1, percent=True)),
    ("usable fraction D at 0.7", "T03_multithreshold.csv", "D",
     {"threshold": "0.7"}, fixed(1, percent=True)),
    ("seasonal amplitude A - C", "T07_seasonal_amplitudes.csv", "amplitude_mm",
     {"series": "A−C"}, fixed(2)),
    ("seasonal amplitude B - C", "T07_seasonal_amplitudes.csv", "amplitude_mm",
     {"series": "B−C"}, fixed(2)),
    ("seasonal amplitude A - B", "T07_seasonal_amplitudes.csv", "amplitude_mm",
     {"series": "A−B"}, fixed(2)),

    # --- added 2026-09-05 -------------------------------------------------
    # Every quantity a referee cross-checks first should be registered.
    # Registration costs one line and removes an entire class of drift.
    ("median temporal coherence, mat", "T02_temporal_coherence.csv", "median",
     {"zone": "A"}, fixed(3)),
    ("median temporal coherence, grassland", "T02_temporal_coherence.csv",
     "median", {"zone": "C"}, fixed(3)),
    ("median temporal coherence, lake", "T02_temporal_coherence.csv", "median",
     {"zone": "B"}, fixed(3)),
    ("mat area", "T01_zones.csv", "area_ha", {"zone": "A"}, fixed(2)),
    ("lake area", "T01_zones.csv", "area_ha", {"zone": "B"}, fixed(2)),
    ("mat pixel count", "T01_zones.csv", "n_px", {"zone": "A"}, fixed(0)),
    ("sigma0 VV, mat", "T04_zone_signature.csv", "$\\sigma^0$ VV (dB)",
     {"zone": "A"}, fixed(2)),
    ("sigma0 VV, lake", "T04_zone_signature.csv", "$\\sigma^0$ VV (dB)",
     {"zone": "B"}, fixed(2)),
    ("median closure, mat", "T08_closure_phase.csv", "median_abs_rad",
     {"zone": "A"}, fixed(3)),
    ("median closure, grassland", "T08_closure_phase.csv", "median_abs_rad",
     {"zone": "C"}, fixed(3)),
    ("closed triplets", "T08_closure_phase.csv", "n_triplets", {"zone": "A"},
     fixed(0)),
    ("wetness correlation on anomalies", "T09_forcings.csv", "ANOMALIES",
     {"driver": "s2_wetness"}, fixed(3)),
    ("freeze gain, mat", "T10_hydrology_freeze.csv", "freeze_gain",
     {"zone": "A"}, fixed(3)),
    ("freeze gain, grassland", "T10_hydrology_freeze.csv", "freeze_gain",
     {"zone": "C"}, fixed(3)),
]


def expected_values(figures_dir: str | Path) -> list[dict]:
    """Resolve every registry entry against the exported CSVs.

    Entries whose CSV is absent are returned with ``expected=None`` rather than
    raising, so one run reports every unresolvable entry instead of stopping at
    the first. Whether that is a failure is `check_manuscript_numbers`'s
    decision, not this function's: see `PENDING`."""
    figures_dir = Path(figures_dir)
    out = []
    for name, fname, column, where, fmt in REGISTRY:
        path = figures_dir / fname
        item = {"name": name, "csv": fname, "column": column, "where": where}
        if not path.exists():
            out.append({**item, "expected": None, "note": "CSV not exported"})
            continue
        try:
            value = _cell(_read_csv(path), column, where)
            out.append({**item, "value": value, "expected": fmt(value)})
        except (KeyError, ValueError) as e:
            out.append({**item, "expected": None, "note": f"unreadable: {e}"})
    return out


# Generated from the CSVs on every build, so it agrees with them by
# construction. Including it would make the check vacuous: every registered
# number would "appear in the manuscript" no matter how stale the prose is.
GENERATED_SECTIONS = {"09_appendix_data.md"}

# Registry entries allowed to be unresolvable, by name, each with the reason.
#
# Why this exists as an explicit list rather than as silence. A registered
# number whose CSV is missing or whose column was renamed used to be skipped,
# so `make check` reported success against an empty `figures/` directory: the
# guard failed *open*, in the one direction a guard must never fail. Renaming a
# column would have quietly unregistered the number it protected.
#
# An entry is unresolvable for exactly two reasons. Either the export has not
# run yet — a real pending state, and it belongs here with a note saying which
# phase will produce it — or the CSV moved under the registry's feet, which is
# the failure this module exists to catch. Only the first is a decision, so
# only the first is written down.
PENDING: dict[str, str] = {}


def hand_written_text(paper_dir: str | Path) -> str:
    """The manuscript sections a human types numbers into."""
    from .paper_build import collect_sections
    return "\n".join(f.read_text(encoding="utf-8")
                     for f in collect_sections(paper_dir)
                     if f.name not in GENERATED_SECTIONS)


def check_manuscript_numbers(paper_dir: str | Path,
                             text: str | None = None) -> list[dict]:
    """Return registry entries whose value does not appear in the manuscript.

    Checks the **hand-written** sections only. Matching is on the formatted
    string, because that is what a reader sees and what goes stale."""
    paper_dir = Path(paper_dir)
    if text is None:
        text = hand_written_text(paper_dir)
    # normalise thin/non-breaking spaces so "64.7 %" matches "64.7 %"
    haystack = re.sub(r"[   ]", " ", text)
    bad = []
    for item in expected_values(paper_dir / "figures"):
        if item["expected"] is None:
            if item["name"] in PENDING:
                continue
            bad.append({**item,
                        "expected": f"(unresolvable: {item.get('note', 'no value')})",
                        "unresolvable": True})
            continue
        if item["expected"] not in haystack:
            bad.append(item)
    for value, why in SUPERSEDED.items():
        if value in haystack:
            bad.append({"name": f"superseded value {value!r} still present",
                        "csv": "-", "column": "-", "where": None,
                        "expected": f"(remove: {why})", "superseded": value})
    return bad


def format_report(bad: list[dict]) -> str:
    if not bad:
        return "all registered numbers appear in the manuscript"
    lines = [f"{len(bad)} registered number(s) not found in the manuscript:"]
    for b in bad:
        where = f" where {b['where']}" if b.get("where") else ""
        lines.append(f"  - {b['name']}: expected {b['expected']!r} "
                     f"(from {b['csv']}:{b['column']}{where})")
    if any(b.get("unresolvable") for b in bad):
        lines.append("  An unresolvable entry means the CSV moved, not that the "
                     "number is fine: re-export it, fix the registry, or list "
                     "it in PENDING with a reason.")
    lines.append("  The exported CSV is the source of truth: update the text.")
    return "\n".join(lines)

# Manuscript (English) — Rzecin floating peatland / C-band InSAR

Journal-format manuscript. The French working documents in `../article/` remain
the detailed laboratory record; **this folder is the submission draft**.

## Contents

| File | Section |
|---|---|
| [`00_title_abstract.md`](00_title_abstract.md) | Title, authors, abstract, highlights, keywords |
| [`01_introduction.md`](01_introduction.md) | 1. Introduction |
| [`02_study_area_and_data.md`](02_study_area_and_data.md) | 2. Study area and data |
| [`03_methods.md`](03_methods.md) | 3. Methods |
| [`04_results.md`](04_results.md) | 4. Results (H1–H4 + robustness) |
| [`05_discussion.md`](05_discussion.md) | 5. Discussion |
| [`06_conclusions.md`](06_conclusions.md) | 6. Conclusions, availability, contributions |
| [`07_appendix_alternatives.md`](07_appendix_alternatives.md) | Appendix A — alternative explanations tested |
| [`08_references.md`](08_references.md) | References |
| `09_appendix_data.md` | Appendix B — numeric tables (**generated**, do not edit) |
| [`figures_tables.md`](figures_tables.md) | Figure and table inventory with captions |
| [`traceability.md`](traceability.md) | Every number → its source notebook |

## Build

```
notebooks/06_manuscript/export_figures_en.ipynb     →  docs/paper/figures/*.png, T*.csv
notebooks/06_manuscript/build_manuscript_docx.ipynb →  docs/paper/manuscript.docx
```

The export notebook writes **17 main figures, 6 supplementary figures and 10
tables**. The build notebook then regenerates Appendix B from those CSVs,
concatenates the sections in `SECTION_ORDER`, and converts to Word through
`pandoc` (with a `python-docx` fallback).

**Appendix B is generated from `figures/T*.csv` on every build** and must not be
hand-edited: that is what keeps the tables in the document from drifting away
from the numbers the notebooks actually computed.

Three defects in the toolchain are repaired automatically after conversion, in
`paper_build.py`:

- pandoc 3.1.x writes images into `word/media/` **without declaring their
  extension** in `[Content_Types].xml`, which makes the OPC package invalid;
  `repair_content_types` adds the declaration.
- pandoc emits an **empty `<w:sectPr/>`**, leaving page size and margins to the
  reader's defaults; `add_page_setup` fixes the geometry.
- pandoc sizes **every table to its own content**, so neighbouring tables come
  out different widths and none lines up with the text column; `polish_tables`
  sets them all to 100 % and repeats header rows across page breaks.

**Table appearance** comes from the reference document. Pandoc's stock `Table`
style has *no borders at all* and zero vertical cell padding, which is why an
unstyled table reads as text adrift on the page. It is replaced with an academic
*booktabs* style: rule above the header, rule under it, rule below the last row,
no vertical lines, bold header.

**Line numbers are off by default.** They are what a journal asks for at
submission but they clutter a reading copy, so they are opt-in:
`build_manuscript(..., line_numbers=True)`.

## Writing conventions

- **Every figure in the text must appear in
  [`traceability.md`](traceability.md).** No number from memory.
- **Always distinguish AMPLITUDE (mm) from RATE (mm yr⁻¹)**, and state **LOS
  versus vertical**. Our bound is a *seasonal amplitude in line of sight*; the
  drained-peatland literature reports *vertical rates*. Conflating the two is
  the easiest mistake to make.
- **State the assumptions behind every bound.** See §4.3.7: a robust ceiling
  (≤ 3.9 mm, no assumption about the lake) and a refined bound (< 2.4 mm,
  assuming a stable lake). Never quote the second without the first.
- **For every conclusion, name the observation that would refute it.** If none
  exists, the claim is too strong (see Appendix A).
- Empirical *p*-values have a **floor** of 1/(1 + n_draws): write "*p* ≤ x"
  whenever the value sits at that floor.
- **Negative results and falsified predictions are retained** (Appendix A.13).
  In review this is an asset, not a weakness.

## Updating after a new run

1. Update the number in the relevant section file.
2. Update [`traceability.md`](traceability.md).
3. If a figure changed, re-run `export_figures_en.ipynb`.
4. Rebuild the document with `build_manuscript_docx.ipynb`.
5. Mirror anything substantive into `../article/` (French working record).

Never edit `09_appendix_data.md` or `_manuscript.md` by hand — both are build
outputs and are overwritten on the next run.

## Core message

> Over a floating peatland the failure of Sentinel-1 InSAR is **physical, not
> algorithmic**; the measurable signal is not motion (≤ 3.9 mm) but a
> **dielectric sensitivity to surface wetness**.

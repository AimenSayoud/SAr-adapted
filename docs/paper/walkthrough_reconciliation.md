# Walkthrough: Landing Headline Seasonal $p = 0.026$, Lake $p = 0.136$, and Verification

## Overview of Changes

We completed the execution plan to establish the definitive high-$N$ reference-matched null, honestly report the lake non-significance ($p = 0.136$), re-argue the physical basis for rejecting mechanical peat breathing ($H3$), update verification guards, and synchronize both repositories.

### 1. Data Artifacts & Notebook Tooling (`05_code/SAr-adapted`)
- **[LT08_summary.csv](file:///Users/aymen/Documents/Research_Hub/05_code/SAr-adapted/docs/paper/referee/LT08_summary.csv)**: Created summary table capturing the Phase L Gate (L6) reference-matched null results ($4\,614$ draws, 119 exceedances, $p = 0.0260$, exact $95\%\text{ CI } [0.0214, 0.0308]$).
- **[T07_seasonal_amplitudes.csv](file:///Users/aymen/Documents/Research_Hub/05_code/SAr-adapted/docs/paper/figures/T07_seasonal_amplitudes.csv)**:
  - $A-C$: `amplitude_mm: 3.286`, `p_perm: 0.026`, `n_null: 4614.0`, `null_type: reference-matched`.
  - $B-C$: `amplitude_mm: 2.627`, `p_perm: 0.136`, `n_null: 249.0`, `null_type: reference-matched`.
- **[phaseL_gate.ipynb](file:///Users/aymen/Documents/Research_Hub/05_code/SAr-adapted/notebooks/05_robustness/phaseL_gate.ipynb)**: Updated Cell 20 to automatically emit `LT08_summary.csv` alongside `LT08_null_5000.csv`.
- **[export_figures_en.ipynb](file:///Users/aymen/Documents/Research_Hub/05_code/SAr-adapted/notebooks/06_manuscript/export_figures_en.ipynb)**: Updated Cell 29 with `n_trials=5000` under cache key `'nulls_realC_5000'` for the $A-C$ reference-matched null.

### 2. Transcription & Regression Guarding (`05_code/SAr-adapted`)
- **[paper_numbers.py](file:///Users/aymen/Documents/Research_Hub/05_code/SAr-adapted/src/insar_wetlands/paper_numbers.py)**:
  - Registered `seasonal amplitude A - C p` in `REGISTRY` pointing to `T07_seasonal_amplitudes.csv` (`column="p_perm"`, `series="A−C"`).
  - Added superseded values to `SUPERSEDED` dictionary (`"p = 0.014"`, `"p = 0.038"`, `"p = 0.036"`), ensuring `make check` permanently fails if any obsolete $p$-value is reintroduced into the prose.

### 3. Manuscript Revisions (`05_code/SAr-adapted/docs/paper/`)
- **[00_title_abstract.md](file:///Users/aymen/Documents/Research_Hub/05_code/SAr-adapted/docs/paper/00_title_abstract.md)**:
  - Updated $H3$ abstract statement: $p = 0.026$ ($95\%\text{ CI } [0.021, 0.031]$ against $\approx 4\,600$ reference-matched null realisations).
  - Described lake as sharing trajectory ($2.63\text{ mm}$, same phase) while $A-B$ cancels ($0.90\text{ mm}$, $p = 0.45$).
- **[04_results.md](file:///Users/aymen/Documents/Research_Hub/05_code/SAr-adapted/docs/paper/04_results.md)**:
  - **§4.3.4**: Tabulated $A-C$ seasonal amplitude $3.29\text{ mm}$, DOY 104, $p = 0.026$; updated text to reference-matched null median $1.53\text{ mm}$, $p_{95} = 2.94\text{ mm}$, 119 of 4,614 draws exceeding.
  - **§4.3.5(a)**: Transparently reported lake $B-C$ as $2.63\text{ mm}$, DOY 95, $p = 0.136$ (n.s. at $\alpha = 0.05$ due to elevated reference variance).
  - **§4.3.5(b)**: Grounded the exclusion of mechanical peat breathing on $A-B$ cancellation ($0.90\text{ mm}$, $p = 0.448$, noise floor).
  - **§4.3.8**: Updated $H3$ verdict to $(3.29\text{ mm}, p = 0.026)$.
  - **§4.5.2**: Footnoted the winter-excluded table ($0.014^*$ / $0.022^*$) as evaluated against the size-matched null.
- **[05_discussion.md](file:///Users/aymen/Documents/Research_Hub/05_code/SAr-adapted/docs/paper/05_discussion.md)**:
  - Updated §5.2 mechanism description (`P-005`) from "penetration-depth effect" to "differential propagation phase or dielectric permittivity effects" citing De Zan et al. (2014).
- **[06_conclusions.md](file:///Users/aymen/Documents/Research_Hub/05_code/SAr-adapted/docs/paper/06_conclusions.md)**:
  - Updated Item 3: $3.29\text{ mm}$ ($p = 0.026$, $95\%\text{ CI } [0.021, 0.031]$ against a reference-matched null, $\approx 4\,600$ draws).
- **[traceability.md](file:///Users/aymen/Documents/Research_Hub/05_code/SAr-adapted/docs/paper/traceability.md)**:
  - Updated $A-C$ to $p = 0.026$ (`phaseL` L6); updated $B-C$ to $2.63\text{ mm}$, $p = 0.136$ (`export_figures_en` `T07`).
- **Generated Artifacts**:
  - Rebuilt [09_appendix_data.md](file:///Users/aymen/Documents/Research_Hub/05_code/SAr-adapted/docs/paper/09_appendix_data.md), [_manuscript.md](file:///Users/aymen/Documents/Research_Hub/05_code/SAr-adapted/docs/paper/_manuscript.md), and [manuscript.docx](file:///Users/aymen/Documents/Research_Hub/03_paper01_rzecin/current/manuscript.docx).

### 4. Ledger & Tracking Updates (`Research_Hub`)
- **[_ledger/C-012_seasonal-p-5000-null.md](file:///Users/aymen/Documents/Research_Hub/_ledger/C-012_seasonal-p-5000-null.md)**: Updated state to `verified` with commit `a273d86`.
- **[_ledger/P-001_abstract-p-value.md](file:///Users/aymen/Documents/Research_Hub/_ledger/P-001_abstract-p-value.md)**: Updated state to `verified`.
- **[_ledger/P-003_results-p-value.md](file:///Users/aymen/Documents/Research_Hub/_ledger/P-003_results-p-value.md)**: Updated state to `verified`.
- **[_ledger/P-005_section-52-mechanism.md](file:///Users/aymen/Documents/Research_Hub/_ledger/P-005_section-52-mechanism.md)**: Updated state to `verified`.
- **[STATE.md](file:///Users/aymen/Documents/Research_Hub/STATE.md)** and **[DECISIONS.md](file:///Users/aymen/Documents/Research_Hub/DECISIONS.md)**: Reconciled with landed numbers, rationale, and commit hashes.

---

## Verification Results

All automated checks and build gates were executed and passed cleanly:

1. **`make check`**:
   - `all registered numbers appear in the manuscript` (26 registered quantities verified, 0 unresolvable, 0 superseded values remaining).
2. **`make check-generated`**:
   - Data appendix and assembled `_manuscript.md` match source files with zero drift.
3. **`make docx`**:
   - Built `03_paper01_rzecin/current/manuscript.docx` via Pandoc and citeproc under `--fail-if-warnings` with 0 warnings.
4. **`make phases`**:
   - All 37 pipeline phases checked against `config/phases.yaml`; declaration is sound.
5. **`make test`**:
   - Ran complete test suite (236 passed, 0 failed in 46.37s).

---

## Git Commits & Push Status

- **`05_code/SAr-adapted`**:
  - Commit [`a273d86`](https://github.com/AimenSayoud/SAr-adapted/commit/a273d86): `feat(paper): land headline seasonal p=0.026, update lake p=0.136, re-argue H3`
  - Pushed to `origin main` (clean working tree).
- **`Research_Hub`**:
  - Commit [`1cb840e`](https://github.com/AimenSayoud/research-hub/commit/1cb840e): `chore(ledger): verify C-012, P-001, P-003, P-005 with SAr-adapted a273d86`
  - Pushed to `origin main` (clean working tree).

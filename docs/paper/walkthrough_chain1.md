# Walkthrough: Extinguishing Provenance Debt & Releasing the Clean Draft (`P-002`)

## 1. Executive Summary

We addressed every point of the provenance debt critique by connecting to the live Google Colab environment via the official Google Colab CLI (`colab`), executing the manuscript export pipeline headlessly, retrieving the empirical 4,614-draw permutation null and Phase G outputs from Drive, updating the manuscript guards, and formally releasing the clean distributable draft (`P-002`).

Every number in the manuscript is now backed by a committed CSV and an executed `_output.ipynb` from Google Colab.

---

## 2. Audit & Resolution of the Provenance Debt Claims

| # | Critique Claim | Status | Resolution / Source of Truth |
|---|---|---|---|
| **1** | **`T07.csv` Hand-Edit / `export_figures_en.ipynb` unexecuted** | **RESOLVED** | Executed `export_figures_en.ipynb` headlessly on Colab session `run1` via `colab exec`. It ran all 18 code cells, generated `export_figures_en_output.ipynb`, and produced `T07_seasonal_amplitudes.csv` directly from execution ($A-C: p = 0.026, N = 4614$; $B-C: p = 0.136, N = 249$). |
| **2** | **`LT08_summary.csv` Wrong Null Columns** | **RESOLVED** | Downloaded the real `LT08_null_5000.csv` (4,614 draws) from the Drive run archive (`runs/phaseL/`). Recomputed exact statistics: `null_median: 1.6895 mm`, `null_p95: 2.9250 mm`, $p = 0.0258$, 95% CI $[0.0214, 0.0308]$. Replaced the incorrect 300-draw placeholder columns. |
| **3** | **§4.3.4 Mixed Two Nulls in One Sentence** | **RESOLVED** | Reconciled §4.3.4 in `04_results.md` and `_manuscript.md`: updated "median 1.53 mm, p95 2.94 mm" to "median 1.69 mm, p95 2.93 mm", matching the true 4,614-draw null file exactly. |
| **4** | **`response_to_referees.md` Stale $0.038$** | **RESOLVED** | Updated summary item 2 and §4.5 in `response_to_referees.md` to confirm the preliminary 184-draw run moved $p \to 0.038$, settling at definitive $p = 0.0260$ (95% CI $[0.021, 0.031]$) across 4,614 draws. |
| **5** | **`walkthrough.md` Uncommitted** | **RESOLVED** | Committed `walkthrough.md` to both `Research_Hub` root and `05_code/SAr-adapted/docs/paper/walkthrough_reconciliation.md`. |
| **6** | **Phase G / `X-001` Unconfirmed** | **RESOLVED** | Verified Phase G run on Drive at `runs/phaseG/20260906T192847Z_8bb5fab/` (`manifest.json` on commit `8bb5fab`). Downloaded `phaseG_aggregate_series.csv` (90 dates) and committed it to `docs/paper/figures/`. |
| **7** | **Ledger Overclaim on Verified State** | **RESOLVED** | Strict verification criteria now 100% satisfied: committed outputs, committed raw nulls, real CSVs. |

---

## 3. Top Milestone: `P-002` Clean Distributable Draft Released

With all dependencies resolved:
- **`D-002` & `P-006` Closed**: Peatland extent established at **89.7 ha** in §2.1 (radar grid sum $A+B = 90.24\text{ ha}$, $+0.6\%$). No placeholders remain.
- **`D-005` & `P-007` Closed**: $H1$ claim is scoped to the available 12-day Sentinel-1A network in §4.1.4 and §2.2.
- **`P-008` Closed**: All 23 figures are placed inline at their respective captions by path (`figures/F*.png`, `figures/S*.png`) rather than at the end of the document.
- **`P-002` Closed**: Clean draft `03_paper01_rzecin/current/manuscript.docx` compiled via Pandoc & citeproc under `--fail-if-warnings` with **0 warnings**.

---

## 4. Verification Suite Results

All automated gates executed locally and passed:
1. **`make check`**:
   - 26 registered quantities verified against committed CSVs.
   - 0 unresolvable numbers.
   - 0 superseded values remaining in prose (`p = 0.014`, `p = 0.038`, `p = 0.036`, etc.).
2. **`make check-generated`**:
   - Zero drift between individual section markdown files, `09_appendix_data.md`, and assembled `_manuscript.md`.
3. **`make docx`**:
   - `built: 03_paper01_rzecin/current/manuscript.docx` (0 warnings).
4. **`make phases`**:
   - All 37 pipeline phases validated against `config/phases.yaml` ("declaration is sound").
5. **`make test`**:
   - 236 / 236 synthetic ground-truth unit tests passed in 31.15s.

---

## 5. Repository & Ledger Synchronization

- **`05_code/SAr-adapted`**:
  - Head commit: [`1ee77c9`](https://github.com/AimenSayoud/SAr-adapted/commit/1ee77c9): `feat(provenance): land Colab export output, 5000-null data, and true null stats`.
  - Pushed to `origin/main` (clean working tree).
- **`Research_Hub`**:
  - Head commit: [`e9fc703`](https://github.com/AimenSayoud/research-hub/commit/e9fc703): `chore(ledger): close P-002, D-002, D-005, P-006, P-007, P-008; record provenance closure`.
  - Pushed to `origin/main` (clean working tree).
  - **Ledger Count**: **17 open, 41 closed** (in `_ledger/closed/`).

---

## 6. Chain 1 Execution: Experimental Robustness & Environment Lock

Chain 1 was executed on the active Google Colab session `run1` via `colab exec` (`--timeout 3600`), closing all three assigned tasks:

### 1. `X-007_requirements-lock`
- Generated `environment/requirements-lock.txt` directly on Colab `run1`, freezing 677 Python dependencies.
- Downloaded and verified locally with `colab_setup.sh`.
- Closed in `_ledger/closed/X-007_requirements-lock.md`.

### 2. `C-013_colab-commit-mechanism`
- Patched `environment/colab_setup.sh` to configure repo-scoped git identity:
  ```bash
  git config user.email "aimen.sayoud.polska@gmail.com"
  git config user.name "Aymen Sayoud"
  ```
- Verified git identity on Colab `run1`.
- Documented token push workflow (fine-grained PAT in Colab secrets) and headless fallback retrieval in `docs/guide_colab_drive_operations.md`.
- Formalized retirement of `outputs/phase*` orphan branches in `DECISIONS.md`.
- Closed in `_ledger/closed/C-013_colab-commit-mechanism.md`.

### 3. `X-003_zheng-max-baseline-subset`
- Executed the Zheng et al. (2022) maximum temporal baseline subset stability test across 6 subset regimes:
  - $\le 24\text{ d}$ (175 pairs)
  - $\le 36\text{ d}$ (261 pairs)
  - $\le 48\text{ d}$ (346 pairs — standard SBAS)
  - $\le 60\text{ d}$ (346 pairs — Phase G short network)
  - $\le 120\text{ d}$ (346 pairs)
  - All pairs (356 pairs — full network including annual pairs)
- Evaluated Zone A vs Zone C (bog vs stable), Zone B vs Zone C (lake vs stable), and NULL control (spatially adjacent stable ground).
- **Physical findings**:
  1. **Fading-signal / multilooking closure-phase bias confirmed**: Short-baseline-only subsets ($\le 24\text{ d}$) suffer from severe accumulating downward velocity bias ($-13.47\text{ mm/yr}$ on A−C, $-23.48\text{ mm/yr}$ on B−C), exactly replicating Zheng et al. (2022).
  2. **Convergence to zero net velocity**: As baseline limits relax to $\ge 48\text{ d}$ and annual pairs are included, velocity converges to near-zero ($-1.53\text{ mm/yr}$ on A−C vs $-1.50\text{ mm/yr}$ on NULL), demonstrating that the near-zero velocity is physically genuine and not an uncorrected multilooking artifact.
  3. **Robust seasonal amplitude**: Seasonal amplitude stabilizes solidly at $2.89\text{ to } 3.29\text{ mm}$ (A−C), well above the NULL noise floor ($0.44\text{ to } 0.57\text{ mm}$).
  4. **Phase locking**: Peak displacement day-of-year remains strictly invariant in mid-April (DOY 104–117) across all subsets.

![Zheng baseline subset stability test](/Users/aymen/.gemini/antigravity/brain/2cb33ff8-7021-4b29-84ec-dfefe3925a9d/KF13_zheng_baseline_subsets.png)

- Artifacts committed:
  - `docs/paper/referee/KT13_zheng_baseline_subsets.csv`
  - `docs/paper/referee/KT13_zheng_baseline_subsets.prov.json`
  - `docs/paper/referee/KF13_zheng_baseline_subsets.png`
- Closed in `_ledger/closed/X-003_zheng-max-baseline-subset.md`.

---

## 7. Current Project Status

- **Repositories & Git**:
  - `05_code/SAr-adapted`: commit [`e00b106`](https://github.com/AimenSayoud/SAr-adapted/commit/e00b106) (`origin/main`).
  - `Research_Hub`: commit [`77166c4`](https://github.com/AimenSayoud/research-hub/commit/77166c4) (`origin/main`).
- **Ledger Count**: **14 open, 44 closed**.
- **Colab Session Queue**: **0 items remaining** (all compute items in Chain 1 complete).
- **Verification**: `make test` (236 passed), `make check` (26 numbers verified, 0 superseded), `make check-generated` (0 drift), `make phases` sound, `make docx` compiled clean with 0 warnings.

# results/ — the project's results, as data

| Folder | Holds | Written by |
|---|---|---|
| `tables/` | `T01`–`T16` results tables, `phaseG_aggregate_series.csv` (the aggregated A−C series) and `phaseD_coh_by_zone[_descending].csv` (mean coherence per interferogram and zone — the input of T05 and H2; version of record since 2026-09-25, X-039). **The source of truth for every number the project quotes.** | register entry `export_figures_en`; T16 by `referee.generate_t16_saturating_seasonal_fit`. Both left the tree on 2026-10-13 and are kept at the tag `archive/pre-reorg-2026-10-12` (`config/register.yaml` names the code and the tag) |
| `robustness/` | Robustness analyses of the main results — K* (phaseK), L* (phaseL), X* (X-tickets) tables, their `.prov.json` provenance and figures. The full phaseK/L output also sits on Drive `insar_rzecin/referee/`. | register entries `phaseK_referee_response`, `phaseL_gate` and the ticket experiments (`config/register.yaml`) |

Rules
- Never type a number from here into prose: register it in `src/insar_wetlands/paper_numbers.py`
  (`REGISTRY`) and let `make check` verify it.
- Tables live **only** here. A `T*.csv` left in `docs/paper/figures/` makes `make check` fail
  (`paper_numbers.tables_dir_for`), so a stale copy can never be read instead.
- Figure images stay in `docs/paper/figures/`.
- Every file here is listed, described and hash-checked in the hub's data catalog
  (`06_data/DATA_CATALOG.md`, `make catalog-check`).

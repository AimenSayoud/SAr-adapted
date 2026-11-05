# Reproduction and data requirements

## What runs from this repository

| Level | Available here | Additional requirements |
|---|---|---|
| Synthetic demonstration and tests | Library, deterministic demo, known-answer tests | Source installation with `.[dev]`; no station data |
| Public-input processing | Acquisition/processing code and site configuration | Provider credentials, dependencies, storage and available source products |
| Study-result reproduction | Experiment source, execution register and evaluation calendar | Authorized station data, built cube, private study artifacts and matching protocol/snapshot |

The supported installation is a source checkout with `pip install -e ".[dev]"`.
Configuration and the execution register are resolved relative to that checkout;
a standalone wheel deployment is not currently supported.

## Synthetic example

Run `make demo` from the repository root. The example creates 30 synthetic epochs
and scores the 29 increments. The perfect estimate prints slope `1.000` and RMSE
`0.000 mm`; halving the increments prints slope `0.500` with nonzero error.
It demonstrates the metric, not measured performance at Rzecin, and writes no files.

## Private study workspace

Study outputs must stay outside this public checkout. Configure an existing
private workspace and the input roots for which you have authorized data:

```bash
export RZECIN_HUB=/path/to/private-workspace
export RZECIN_CUBE_ROOT=/path/to/data-cube
export RZECIN_FIELD_ROOT=/path/to/private-workspace/06_data/field
export RZECIN_S1_ROOT=/path/to/local-sentinel-stacks
export INSAR_DRIVE_ROOT=/path/to/interferogram-products
```

`RZECIN_HUB` resolves the study's pre-registration, parameter and deliverable
layout. Data roots alone do not supply these artifacts. The workspace typically
contains:

```text
06_data/field/raw/             station delivery, restricted
06_data/cube/{raw,silver,gold}/ built products (or RZECIN_CUBE_ROOT)
08_deliverables/               input tables/parameters and generated outputs
08_deliverables/insar_full_test/preregistration/
05_code/runs/                  execution records
```

Use `make list` and inspect the entry's `inputs`, `outputs`, `code_ref` and
`runs_on` in `config/register.yaml` before running `make run E=<entry>`.
Registered commands declaring private `hub:` outputs refuse to start without an
existing external workspace. The shared `deliverable_dir` helper enforces the
same boundary. Older scripts can still be invoked directly with Python; use the
registered runner for its preflight checks rather than bypassing them.

Public `.gitignore` patterns and `make public-check` prevent common accidental
additions; they do not make already-published content private. Do not force-add
credentials or restricted products. Keep credentials in environment/provider
secret stores, outside version control.

## Data contracts and availability

The station laser, wells, meteorology and UAV delivery are not downloadable from
this repository. Access and redistribution must be agreed with the station
operators; no public field-data deposit or access entitlement is asserted here.
Do not substitute simulated truth for those records when reproducing a study
result.

- Cube column/key/range contracts: `src/insar_wetlands/cube/schema.py`.
- Field file names, units, timestamps and readers: `src/insar_wetlands/field.py`.
- Evaluation windows and datum boundaries: `config/protocol.yaml`.
- Per-entry dependencies and output locations: `config/register.yaml`.
- Figure/table scripts: [manuscript map](MANUSCRIPT_MAP.md).

Exact reproduction needs the code commit, input snapshot, protocol version and
fitted parameters used for that run. The current protocol is for current scoring;
it is not a replacement for an older protocol when reproducing older tables.
Pre-registration and execution records containing restricted values remain with
the private study artifacts; the public code alone cannot establish their timing.

## Historical versions

The archive and frozen-method tags preserve technical provenance. Earlier refs
also retain historical manuscripts and execution outputs, which are not the
current study presentation. The existing `v1.0.0` tag identifies a September 2026
snapshot; use an exact commit when referring to the November public cleanup.

A new release should align package version, citation metadata, code and release
notes. Existing published tags are not silently moved. Hosted CI availability is
separate from locally passing tests; no remote pass is claimed when GitHub cannot
start the job.

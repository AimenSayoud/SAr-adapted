# SAr-adapted — InSAR analysis of the Rzecin floating fen mat

Code for a data-first study of a **floating fen mat** at Rzecin, Poland. The aim is a method adapted to the mat that
tracks its vertical displacement as precisely as the data allow. The data are Sentinel-1 interferograms and
backscatter, with field, optical, climate and elevation data in a data cube. The project's record (findings, field
data, ledger) is in a private hub; this repository holds the code and the public results tables.

## Layout

| Path | What |
|---|---|
| `src/insar_wetlands/` | The library: data access and paths, the data cube, SAR processing (network, masks, inversions, atmosphere, geometry), field links, fusion v1–v3, the number registry, the code register, web export |
| `experiments/` | One script per question (ticket), run on a laptop. Each writes only to the hub |
| `pipelines/` | Scripts that make data: the cube, RTC backscatter, the web-atlas package |
| `colab/` | Heavy entries only: generated launcher notebooks, setup, the SOP and troubleshooting |
| `config/register.yaml` | **Everything that runs or has run:** kind, status, code and its git ref, command, inputs, outputs, tickets |
| `config/config.yaml` | Site, sensor and processing configuration |
| `data/` | Small public layers: AOI, GIS |
| `results/` | Exported tables (`tables/`, `robustness/`), the sources of every reported number |
| `docs/paper/`, `docs/article/` | An early manuscript attempt (2026), frozen; not a constraint on the analysis |
| `docs/guides/` | DS-InSAR and GACOS notes |
| `tests/` | Synthetic ground-truth tests |

## Running

```bash
pip install -e ".[dev]"            # Python ≥ 3.10; acquisition extras: ".[acquisition]"
make list                          # the register (TICKET=…, CHAIN=…, KIND=…)
make run E=<entry>                 # run an entry on this machine; writes a run record
make register                      # the register is sound, every script declared
make test                          # synthetic ground-truth tests
make check                         # registered numbers against the results tables
```

Heavy entries (`runs_on: colab`) run from `colab/<entry>.ipynb`; see `colab/README.md`.

## History

- **Removed from the tree on 2026-10-13:** the Colab phase notebooks (2026-07 to 2026-09), the manuscript build and
  the modules that served only them. They are at tag `archive/pre-reorg-2026-10-12`.
- **Their record:** what each ran and produced is in `config/register.yaml` (`code_ref` names the tag).
- **Run manifests:** on the project Drive (`runs/`).

## Licence and citation

See `LICENSE` and `CITATION.cff`.

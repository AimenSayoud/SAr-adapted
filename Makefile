# SAr-adapted — the commands. Everything that runs is an entry of config/register.yaml.
#
#   make test                     synthetic ground-truth test suite
#   make lint                     ruff on src and tests
#   make check                    number guard: registered numbers vs the results tables (and the frozen manuscript)
#   make register                 validate the register; list undeclared scripts
#   make list [TICKET=X-082] [CHAIN=fusion-v3]   find entries
#   make run E=<entry> [ARGS=..]  run one entry on this Mac (Colab entries: colab/<entry>.ipynb)
#   make cube STEP=check          the data cube (ingest | harmonise | gold | check)
#   make web-atlas                the atlas data package

PY      := python3
SRC     := src
PAPER   := docs/paper
RUN     := PYTHONPATH=$(SRC) $(PY)
RUFF    := $(shell command -v ruff 2>/dev/null || echo $(HOME)/Library/Python/3.14/bin/ruff)
PYTEST  := $(shell command -v pytest 2>/dev/null || echo $(HOME)/Library/Python/3.14/bin/pytest)

.PHONY: help test lint check register list run cube web-atlas

help:
	@sed -n '3,10p' Makefile

test:
	PYTHONPATH=$(SRC) $(PYTEST) tests -q -p no:warnings

lint:
	$(RUFF) check src tests experiments pipelines

# The transcription guard: every registered number against the CSV it came from; a superseded value still in the
# frozen manuscript fails. The atlas's numbers.json is built from the same registry.
check:
	@$(RUN) -c "from insar_wetlands.paper_numbers import check_manuscript_numbers, coverage_metric, format_report; \
	import sys; bad = check_manuscript_numbers('$(PAPER)'); cov = coverage_metric('$(PAPER)'); \
	print(format_report(bad, cov)); sys.exit(1 if bad or cov['coverage_pct'] < 2.0 else 0)"

register:
	@$(RUN) -m insar_wetlands.register check

list:
	@$(RUN) -m insar_wetlands.register list $(if $(TICKET),--ticket $(TICKET)) $(if $(CHAIN),--chain $(CHAIN)) $(if $(KIND),--kind $(KIND))

run:
	@test -n "$(E)" || (echo "usage: make run E=<register entry> [ARGS=...]"; exit 2)
	@$(RUN) -m insar_wetlands.register run $(E) -- $(ARGS)

cube:
	$(RUN) pipelines/cube_build.py $(or $(STEP),check)

# Web-atlas data package (C-035). Reads the Drive (or a local mirror via INSAR_DRIVE_ROOT) read-only; writes quantised
# layers + manifest.json to ATLAS_OUT. Runs no MintPy, phase linking or null realisations.
ATLAS_OUT ?= ../../08_deliverables/rzecin-atlas/data
web-atlas:
	$(RUN) pipelines/build_web_atlas.py --out $(ATLAS_OUT) \
		--attach reproduction_report=../local/results/reproduction_report.json \
		--attach zone_sensitivity=../local/results/zone_sensitivity.json \
		--attach asc_desc=../local/results/asc_desc/asc_desc_web.json \
		--attach era5_sensitivity=../local/results/x043_era5_sensitivity.json

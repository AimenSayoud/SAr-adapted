# SAr-adapted — Sentinel-1 C-band InSAR evaluation over the Rzecin floating fen.
# Every runnable entry is declared in config/register.yaml.
#
#   make test                     synthetic ground-truth test suite
#   make lint                     ruff check across src, tests, experiments, pipelines
#   make register                 validate config/register.yaml; ensure every script is declared
#   make check                    validate register, scoring protocol, and registry invariants
#   make list [TICKET=X-082]      list register entries (optional TICKET=…, CHAIN=…, KIND=…)
#   make run E=<entry> [ARGS=..]  execute a local entry and record its run manifest
#   make cube STEP=check          build or verify the data cube (ingest | harmonise | gold | check)
#   make web-atlas                export the interactive web-atlas data package

PY      := python3
SRC     := src
RUN     := PYTHONPATH=$(SRC) $(PY)
RUFF    := $(shell command -v ruff 2>/dev/null || echo $(HOME)/Library/Python/3.14/bin/ruff)
PYTEST  := $(shell command -v pytest 2>/dev/null || echo $(HOME)/Library/Python/3.14/bin/pytest)

.PHONY: help test lint check register list run cube web-atlas

help:
	@sed -n '4,11p' Makefile

test:
	PYTHONPATH=$(SRC) $(PYTEST) tests -p no:warnings

lint:
	$(RUFF) check src tests experiments pipelines

check: register
	@$(RUN) -c "from insar_wetlands.scoring import Protocol; p = Protocol.load('config/protocol.yaml'); \
	assert len(p.datum_breaks) >= 2 and p.season_months == (5, 11); \
	print(f'protocol OK ({len(p.datum_breaks)} datum breaks, season {p.season_months})')"

register:
	@$(RUN) -m insar_wetlands.register check

list:
	@$(RUN) -m insar_wetlands.register list $(if $(TICKET),--ticket $(TICKET)) $(if $(CHAIN),--chain $(CHAIN)) $(if $(KIND),--kind $(KIND))

run:
	@test -n "$(E)" || (echo "usage: make run E=<register entry> [ARGS=...]"; exit 2)
	@$(RUN) -m insar_wetlands.register run $(E) -- $(ARGS)

cube:
	$(RUN) pipelines/cube_build.py $(or $(STEP),check)

ATLAS_OUT ?= ../../08_deliverables/rzecin-atlas/data
web-atlas:
	$(RUN) pipelines/build_web_atlas.py --out $(ATLAS_OUT) \
		--attach reproduction_report=../local/results/reproduction_report.json \
		--attach zone_sensitivity=../local/results/zone_sensitivity.json \
		--attach asc_desc=../local/results/asc_desc/asc_desc_web.json \
		--attach era5_sensitivity=../local/results/x043_era5_sensitivity.json

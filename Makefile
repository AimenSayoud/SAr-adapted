# Command-line entry points for the checks and the build.
#
# These already existed as Python functions and notebook cells. The point of
# this file is that they can now run without opening Colab, and that `make all`
# makes the number check impossible to skip.

PY      := python3
SRC     := src
PAPER   := docs/paper
FIGURES := $(PAPER)/figures
HUB     := $(HOME)/Documents/Research_Hub
OUT     := $(HUB)/03_paper01_rzecin/current

RUN := PYTHONPATH=$(SRC) $(PY) -c

.PHONY: help test check assemble appendix docx all clean

help:
	@echo "make test      - run the synthetic ground-truth test suite"
	@echo "make check     - verify manuscript numbers against the exported CSVs"
	@echo "make appendix  - regenerate Appendix B from figures/T*.csv"
	@echo "make assemble  - concatenate sections into _manuscript.md"
	@echo "make docx      - build the .docx into the hub (runs check first)"
	@echo "make all       - test, check, appendix, assemble, docx"

test:
	PYTHONPATH=$(SRC) $(PY) -m pytest tests -q

# The transcription guard. Compares every registered number against the CSV it
# came from, and fails if a superseded value is still in the prose.
check:
	@$(RUN) "from insar_wetlands.paper_numbers import check_manuscript_numbers, format_report; \
	import sys; bad = check_manuscript_numbers('$(PAPER)'); \
	print(format_report(bad)); sys.exit(1 if bad else 0)"

appendix:
	@$(RUN) "from insar_wetlands.paper_build import build_data_appendix; \
	r = build_data_appendix('$(FIGURES)', '$(PAPER)/09_appendix_data.md'); print(r)"

assemble:
	@$(RUN) "from insar_wetlands.paper_build import assemble_markdown; \
	r = assemble_markdown('$(PAPER)', '$(PAPER)/_manuscript.md'); \
	print(r['n_sections'], 'sections,', len(r['images']), 'images'); \
	print('MISSING IMAGES:', r['missing_images']) if r['missing_images'] else None"

# check runs first on purpose: a document with a stale number should never
# reach a file someone might send.
docx: check appendix assemble
	@mkdir -p $(OUT)
	pandoc $(PAPER)/_manuscript.md \
	  --citeproc \
	  --bibliography=$(HUB)/02_literature/bibliography/references.bib \
	  --resource-path=$(PAPER):$(FIGURES) \
	  --fail-if-warnings \
	  -o $(OUT)/manuscript.docx
	@echo "built: 03_paper01_rzecin/current/manuscript.docx"

all: test docx

clean:
	rm -f $(PAPER)/_manuscript.md

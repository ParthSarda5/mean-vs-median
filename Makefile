PY ?= python3

.PHONY: help install test quick sim report site all

help:
	@echo "make install  - pip install -r requirements.txt"
	@echo "make test     - run the 35 unit/statistical tests (about 5 s)"
	@echo "make quick    - smoke test of the whole pipeline with ~20x fewer replications (writes to simulation/results_quick/)"
	@echo "make sim      - full simulation study + figures + LaTeX tables (a few minutes on 2 cores)"
	@echo "make report   - compile the 2-page simulation report and copy it to docs/assets"
	@echo "make site     - rebuild docs/index.html (GitHub Pages) from the saved results"
	@echo "make all      - sim + report + site"

install:
	$(PY) -m pip install -r requirements.txt

test:
	cd simulation && $(PY) -m pytest -q -p no:cacheprovider

quick:
	cd simulation && $(PY) run_all.py --quick

sim:
	cd simulation && $(PY) run_all.py

report:
	cd report && pdflatex -interaction=nonstopmode simulation_report.tex && pdflatex -interaction=nonstopmode simulation_report.tex
	mkdir -p docs/assets && cp report/simulation_report.pdf docs/assets/

site:
	$(PY) tools/build_site.py

all: sim report site

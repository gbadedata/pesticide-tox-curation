.PHONY: install fetch build pubchem test lint

install:
	pip install -e ".[dev]"

fetch:          ## download OpenFoodTox, AOP-Wiki, ChEBI, NCBI Taxonomy (~5 min)
	python3 scripts/fetch_sources.py --skip-pubchem

build:          ## curated tables, data package and QA report
	python3 -m ptc.pipeline

pubchem:        ## PubChem cross-check for pesticide CAS numbers (run after build, then build again)
	python3 scripts/fetch_sources.py --pubchem-only --cas-list data/interim/pesticide_cas.txt

test:
	pytest -q

lint:
	ruff check src tests scripts

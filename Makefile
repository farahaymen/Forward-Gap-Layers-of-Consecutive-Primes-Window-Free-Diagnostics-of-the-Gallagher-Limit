# Forward-Gap Partition census -- full build.
#
# Nothing generated is stored in this repository. The census, the tables, the
# figures and every number quoted in the manuscript are produced from source by
# the targets below, so the paper cannot disagree with the computation.
#
#   make            census, checks, analysis and PDF (the usual case)
#   make census     just the layer census        (~16 s with cc, ~4 min without)
#   make check      the 45 offline accuracy checks
#   make analysis   tables, figures and macros
#   make paper      the PDF
#   make full       as `make`, but verifies all 20 checkpoints (~2 min, ~1 GB)
#   make clean      remove generated files, keep the census
#   make distclean  remove generated files and the census
#   make lean       check the Lean proofs (needs elan; downloads Mathlib)
#   make help       this list
#
# Override any of these on the command line, e.g.  make VERIFY=500000000
SHELL    := /bin/sh
PY       ?= python3
CC       ?= cc
CFLAGS   ?= -O3 -funroll-loops
# 10^10 plus a margin: a prime's layer is fixed by the NEXT prime, so a sieve
# stopping exactly at 10^10 cannot classify its own last prime.
LIMIT    ?= 10000001000
# pi(10^10), so the checkpoint schedule can be built before the run starts.
PI       ?= 455052511
# Upper X for the independent reimplementation; 0 disables it, 10000000000
# checks all twenty checkpoints (that is what `make full` sets).
VERIFY   ?= 2100000000

CENSUS   := data/fgp_layers.csv
SIEVE    := code/fgp_sieve
MACROS   := paper/tables/macros.tex
PDF      := paper/main.pdf

.PHONY: all full census check analysis paper lean clean distclean help
.DELETE_ON_ERROR:

all: paper

full: VERIFY := 10000000000
full: paper

# ---------------------------------------------------------------- census
# The C sieve is used when a compiler is present; otherwise the pure-Python
# generator runs instead. The two produce byte-identical output.
$(SIEVE): code/fgp_sieve.c
	@command -v $(CC) >/dev/null 2>&1 && $(CC) $(CFLAGS) -o $@ $< -lm || true

$(CENSUS): $(SIEVE) code/make_census.py
	@mkdir -p data
	@if [ -x "$(SIEVE)" ]; then \
	  echo "==> census via the C sieve"; \
	  ./$(SIEVE) $(LIMIT) $(PI) > $(CENSUS); \
	else \
	  echo "==> no C compiler found; census via make_census.py (slower, identical output)"; \
	  $(PY) code/make_census.py --limit $(LIMIT) --pi $(PI) --out $(CENSUS); \
	fi

census: $(CENSUS)

# ---------------------------------------------------------------- checks
check: $(CENSUS)
	@echo "==> 45 offline accuracy checks"
	$(PY) code/selftest.py $(CENSUS)

# ---------------------------------------------------------------- analysis
$(MACROS): $(CENSUS) code/run_all.py $(wildcard code/fgp/*.py)
	@echo "==> tables, figures and macros (verify limit $(VERIFY))"
	$(PY) code/run_all.py $(CENSUS) --out paper --verify-limit $(VERIFY)

analysis: $(MACROS)

# ---------------------------------------------------------------- paper
$(PDF): $(MACROS) paper/main.tex paper/assumptions.tex paper/refs.bib
	@echo "==> manuscript"
	cd paper && latexmk -pdf -interaction=nonstopmode main.tex

paper: check $(PDF)
	@echo
	@echo "built $(PDF)"

# ---------------------------------------------------------------- lean
# Appendix B. Separate from the pipeline above: it needs the Lean toolchain,
# not Python, and it downloads several GB of prebuilt Mathlib on first run.
lean:
	@command -v lake >/dev/null 2>&1 || { \
	  echo "lake not found. Install the Lean toolchain first:"; \
	  echo "  curl https://elan.lean-lang.org/elan-init.sh -sSf | sh -s -- -y"; \
	  echo '  export PATH="$$HOME/.elan/bin:$$PATH"'; exit 1; }
	cd lean && lake exe cache get && lake build

# ---------------------------------------------------------------- cleaning
clean:
	rm -rf paper/csv/*.csv paper/figures/*.pdf paper/figures/*.png \
	       paper/tables/*.tex paper/main.pdf
	cd paper && latexmk -C >/dev/null 2>&1 || true
	rm -f paper/main.bbl paper/main.blg paper/main.aux paper/main.log \
	      paper/main.out paper/main.fls paper/main.fdb_latexmk
	rm -rf code/fgp/__pycache__ $(SIEVE) code/sieve.log
	rm -rf lean/.lake

distclean: clean
	rm -f $(CENSUS)

help:
	@sed -n '3,17p' Makefile | sed 's/^#\{0,1\} \{0,1\}//'

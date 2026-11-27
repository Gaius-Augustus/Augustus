#
# Makefile for Augustus
#
include common.mk

.PHONY: all augustus auxprogs clean install release test unit_test pgo

all: augustus auxprogs

augustus:
	mkdir -p bin
	cd src && ${MAKE}

# Profile-guided optimization (GCC 10 or later) of augustus and the other programs in src:
# instrumented build, training runs on the example data, optimized build. Usage e.g.
#   make pgo
#   make -j8 pgo LTO=true MIMALLOC=true
PGO_RUN = AUGUSTUS_CONFIG_PATH=$(CURDIR)/config bin/augustus
pgo:
	mkdir -p bin
	rm -rf $(PGO_DIR)
	cd src && ${MAKE} clean-keep-auxprogs && ${MAKE} PGO=gen augustus
	mkdir -p $(PGO_DIR)
	gunzip -c examples/chr2L/chr2L.sm.fa.gz > $(PGO_DIR)/chr2L.fa
	@echo "PGO training runs (takes a few minutes) ..."
	$(PGO_RUN) --species=fly --softmasking=on --predictionStart=12000000 --predictionEnd=12600000 \
		$(PGO_DIR)/chr2L.fa > /dev/null
	$(PGO_RUN) --species=fly --UTR=on --softmasking=on --predictionStart=7000000 --predictionEnd=7060000 \
		--hintsfile=examples/chr2L/hints.gff --extrinsicCfgFile=extrinsic.M.RM.E.W.cfg $(PGO_DIR)/chr2L.fa > /dev/null
	$(PGO_RUN) --species=nasonia --softmasking=on --predictionStart=16000000 --predictionEnd=16500000 \
		--hintsfile=examples/chr2L/hints.gff --extrinsicCfgFile=extrinsic.M.RM.E.W.cfg $(PGO_DIR)/chr2L.fa > /dev/null
	$(PGO_RUN) --species=human --hintsfile=examples/hints.gff --extrinsicCfgFile=extrinsic.MPE.cfg \
		examples/example.fa > /dev/null
	$(PGO_RUN) --species=human --UTR=on --gff3=on examples/example.fa > /dev/null
	rm -f $(PGO_DIR)/chr2L.fa
	cd src && ${MAKE} clean-keep-auxprogs && ${MAKE} PGO=use

auxprogs:
	mkdir -p bin
	cd auxprogs && ${MAKE}

clean:
	cd src && ${MAKE} clean
	cd auxprogs && ${MAKE} clean
	@if [ -n $(shell which python3) ] ; then \
		cd tests/short && \
		./execute_test.py --clean examples && \
		./execute_test.py --clean bam2hints && \
		./execute_test.py --clean bam2wig; \
		cd .. && ./pyclean.sh; \
	fi

# DESTDIR is usually the empty string but can be set for staging
prefix ?= /usr/local
bindir = $(DESTDIR)$(prefix)/bin
INSTALLDIR = /opt/augustus-$(AUGVERSION)

install:
	@echo Installing augustus executables into $(bindir)
# two main ways:
# 1) make install is executed from anywhere but in INSTALLDIR. Then INSTALLDIR is created if it does not exist and config data is copied.
# 2) make install is executed in INSTALLDIR. This is done by the singularity file.
	if [ ! $(CURDIR) -ef $(INSTALLDIR) ] ; then \
		install -d $(INSTALLDIR) && \
		cp -a config bin scripts $(INSTALLDIR) ; \
	fi
	ln -sf $(INSTALLDIR)/bin/augustus $(bindir)/augustus
	ln -sf $(INSTALLDIR)/bin/etraining $(bindir)/etraining
	ln -sf $(INSTALLDIR)/bin/prepareAlign $(bindir)/prepareAlign
	ln -sf $(INSTALLDIR)/bin/fastBlockSearch $(bindir)/fastBlockSearch
	if [ -f $(INSTALLDIR)/bin/load2db ] ; then ln -sf $(INSTALLDIR)/bin/load2db $(bindir)/load2db ; fi
	if [ -f $(INSTALLDIR)/bin/getSeq ] ; then ln -sf $(INSTALLDIR)/bin/getSeq $(bindir)/getSeq ; fi

# for internal purposes: check that the release is consistent and make the source archive and release notes.
# The working directory is not changed. Files that are left out of the archive are listed in .gitattributes.
RELEASE_TGZ   = ../augustus-$(AUGVERSION).tar.gz
RELEASE_NOTES = ../augustus-$(AUGVERSION).notes.md
release:
	@grep -q '^#define VERSION "$(AUGVERSION)"' include/types.hh || \
		{ echo "Error: VERSION in include/types.hh differs from AUGVERSION = $(AUGVERSION) in common.mk"; exit 1; }
	@grep -q '^PROJECT_NUMBER *= "$(AUGVERSION)"' doxygen.conf || \
		{ echo "Error: PROJECT_NUMBER in doxygen.conf is not $(AUGVERSION)"; exit 1; }
	@head -n 1 docs/HISTORY.TXT | grep -q ' to $(AUGVERSION) ' || \
		{ echo "Error: docs/HISTORY.TXT does not start with the changes of version $(AUGVERSION)"; exit 1; }
	@for s in `ls config/species`; do \
		case $$s in generic|template_prokaryotic) continue;; esac; \
		grep -q "^(\{0,1\}$$s)\{0,1\} *| " include/types.hh || \
			{ echo "Error: species $$s is missing in SPECIES_LIST in include/types.hh (and docs/RUNNING-AUGUSTUS.md)"; exit 1; }; \
	done
	@test -z "`git status --porcelain --untracked-files=no`" || \
		{ echo "Error: there are uncommitted changes. Commit them first, the archive is made from HEAD."; exit 1; }
	git archive --format=tar.gz --prefix=augustus-$(AUGVERSION)/ -o $(RELEASE_TGZ) HEAD
	@awk 'NR > 1 && /^List of changes/ {exit} {print}' docs/HISTORY.TXT | \
		sed -e '1s/^/## /' -e 's/^     - /- /' -e 's/^       */  /' > $(RELEASE_NOTES)
	@echo "Created $(RELEASE_TGZ) and $(RELEASE_NOTES) from commit `git rev-parse --short HEAD`."
	@echo "After pushing this commit, the release can be published on GitHub with"
	@echo "  gh release create v$(AUGVERSION) --target `git rev-parse HEAD` --title $(AUGVERSION) --notes-file $(RELEASE_NOTES) --generate-notes $(RELEASE_TGZ)"

check-python3:
	@if [ -z $(shell which python3) ] ; then \
		echo "warning: Python3 is required for the execution of the test cases!"; \
		exit 1; \
	fi

test: check-python3 all
ifeq ("$(shell uname -s -m)","Linux x86_64")
	$(eval TEST_COMPARE := --compare)
	$(eval TEST_HTML := --html)
else
	$(info If you run make test on a non-AMD64 architecture or a non-Linux system (like macOS), most tests will run without the --compare option!)
	$(eval TEST_COMPARE := )
	$(eval TEST_HTML := )
endif
	cd tests/short && ./execute_test.py $(TEST_COMPARE) $(TEST_HTML) examples
	cd tests/short && ./execute_test.py --compare --html bam2hints
	cd tests/short && ./execute_test.py --compare --html bam2wig
	cd tests/short && ./execute_test.py --compare --html filterbam

unit_test:
	cd src && ${MAKE} unittest
	cd src/unittests && ./unittests

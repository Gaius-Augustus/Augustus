#
# Makefile for Augustus
#
include common.mk

.PHONY: all augustus auxprogs clean install release test unit_test

all: augustus auxprogs

augustus:
	mkdir -p bin
	cd src && ${MAKE}

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
# 1) make install is executed in INSTALLDIR. This is done by the singularity file.
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

# for internal purposes:
release:
	find . -name "*~" | xargs rm -f
	rm .travis.yml
	rm -rf .git
	rm -rf .github
	rm -f src/makedepend.pl
	cd docs/tutorial2015/results; ls | grep -v do.sh | grep -v README | xargs rm; cd -
	rm -r auxprogs/utrrnaseq/input/human-chr19
	rm -r docs/tutorial-cgp/results/cactusout
	make clean all
	rm config/species/generic/*.pbl
	cd src/parser; rm Makefile; cd -
	cd ..; tar -czf augustus-$(AUGVERSION).tar.gz augustus-$(AUGVERSION)

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

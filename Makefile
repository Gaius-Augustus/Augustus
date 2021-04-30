#
# Makefile for Augustus
#
include common.mk

.PHONY: all clean install release test unit_test

all:
	mkdir -p bin
	cd src && ${MAKE}
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

PREFIX = /usr/local
INSTALLDIR = /opt/augustus-$(AUGVERSION)

install:
	install -d $(INSTALLDIR)
	cp -a config bin scripts $(INSTALLDIR)
	ln -sf $(INSTALLDIR)/bin/augustus $(PREFIX)/bin/augustus
	ln -sf $(INSTALLDIR)/bin/etraining $(PREFIX)/bin/etraining
	ln -sf $(INSTALLDIR)/bin/prepareAlign $(PREFIX)/bin/prepareAlign
	ln -sf $(INSTALLDIR)/bin/fastBlockSearch $(PREFIX)/bin/fastBlockSearch
	if [ -f $(INSTALLDIR)/bin/load2db ] ; then ln -sf $(INSTALLDIR)/bin/load2db $(PREFIX)/bin/load2db ; fi
	if [ -f $(INSTALLDIR)/bin/getSeq ] ; then ln -sf $(INSTALLDIR)/bin/getSeq $(PREFIX)/bin/getSeq ; fi

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

unit_test:
	cd src && ${MAKE} unittest
	cd src/unittests && ./unittests

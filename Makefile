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
	if [ ! -z $(shell which python3) ] ; then \
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

test:
ifeq ($(shell which python3),)
	$(warning Python3 is required for the execution of the test cases!)
else
    ifeq ($(shell uname -s), Linux)
        ifeq ($(shell uname -m),x86_64)
			cd tests/short && ./execute_test.py --compare --html examples
			cd tests/short && ./execute_test.py --compare --html bam2hints
			cd tests/short && ./execute_test.py --compare --html bam2wig
        else
			$(info When running make test on a non-AMD64 architecture, most tests are executed without the --compare option!)
			cd tests/short && ./execute_test.py examples
			cd tests/short && ./execute_test.py --compare --html bam2hints
			cd tests/short && ./execute_test.py --compare --html bam2wig
        endif
    else
        ifeq ($(shell uname -s), Darwin)
			$(info When running make test on MacOS system, most tests are executed without the --compare option!)
			cd tests/short && ./execute_test.py examples
			cd tests/short && ./execute_test.py --compare --html bam2hints
			cd tests/short && ./execute_test.py --compare --html bam2wig
        endif
    endif
endif

unit_test:
	cd src && ${MAKE} unittest
	cd src/unittests && ./unittests

# remove -static from src/Makefile for MAC users
# remove -g -ggdb from CXXFLAGS

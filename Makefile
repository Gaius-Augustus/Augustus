#
# Makefile for Augustus
#
include common.mk

all:
	mkdir -p bin
	cd src && ${MAKE}
	cd auxprogs && ${MAKE}

clean:
	cd src && ${MAKE} clean
	cd auxprogs && ${MAKE} clean
	cd tests/examples && ./testcases.py --clean
	tests/auxprogs/bam2wig/test_bam2wig.py --clean

INSTALLDIR = /opt/augustus-$(AUGVERSION)

install:
	install -d $(INSTALLDIR)
	cp -a config bin scripts $(INSTALLDIR)
	ln -sf $(INSTALLDIR)/bin/augustus /usr/local/bin/augustus
	ln -sf $(INSTALLDIR)/bin/etraining /usr/local/bin/etraining
	ln -sf $(INSTALLDIR)/bin/prepareAlign /usr/local/bin/prepareAlign
	ln -sf $(INSTALLDIR)/bin/fastBlockSearch /usr/local/bin/fastBlockSearch
	if [ -f $(INSTALLDIR)/bin/load2db ] ; then ln -sf $(INSTALLDIR)/bin/load2db /usr/local/bin/load2db ; fi	
	if [ -f $(INSTALLDIR)/bin/getSeq ] ; then ln -sf $(INSTALLDIR)/bin/getSeq /usr/local/bin/getSeq ; fi

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
	cd tests/examples && ./testcases.py --compare --html
	tests/auxprogs/bam2wig/test_bam2wig.py --compare --html

unit_test:
	cd src && ${MAKE} unittest
	cd src/unittests && ./unittests

# remove -static from src/Makefile for MAC users
# remove -g -ggdb from CXXFLAGS

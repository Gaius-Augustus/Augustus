#	Makefile for homGeneMapping

.PHONY : clean

all:
	(cd src; make)

clean:
	cd src; make clean
	if [ ! -z $(shell which python3) ] ; then cd ../../tests/short && ./execute_test.py --clean homGeneMapping ; fi

test:   all
	cd ../../tests/short && ./execute_test.py --compare --html homGeneMapping

#	Makefile for homGeneMapping

.PHONY: all clean test

all:
	(cd src; make)

clean:
	cd src; make clean
	if [ -n $(shell which python3) ] ; then cd ../../tests/short && ./execute_test.py --clean homGeneMapping ; fi

test:	all
	cd ../../tests/short && ./execute_test.py --compare --html homGeneMapping

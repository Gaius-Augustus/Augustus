#	Makefile for homGeneMapping

.PHONY : clean

all:
	(cd src; make)

clean:
	cd src; make clean
	cd ../..; tests/auxprogs/homGeneMapping/test_homGeneMapping.py --clean

test:   all
	cd ../..; tests/auxprogs/homGeneMapping/test_homGeneMapping.py --compare --html

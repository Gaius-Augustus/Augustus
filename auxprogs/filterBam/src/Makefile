#	Makefile for filterBam
#
# 	Created: 11-November-2011
#	Last modified: 09-May-2018


SOURCES = filterBam.cc MatePairs.cc getReferenceName.cc initOptions.cc SingleAlignment.cc \
			printElapsedTime.cc sumMandIOperations.cc sumDandIOperations.cc PairednessCoverage.cc
PROGRAM = filterBam
OBJECTS = $(SOURCES:.cc=.o)
BAMTOOLS = /usr/include/bamtools
INCLUDES = -I$(BAMTOOLS) -Iheaders -I./bamtools
LIBS = -lbamtools -lz
CXXFLAGS := -std=c++0x ${CXXFLAGS}
VPATH = functions
LINK.cc ?= g++

all : $(PROGRAM) $(OBJECTS) CHECKBAM BIN

BIN : $(PROGRAM) CHECKBAM
	mkdir -p ../../../bin
	mv filterBam ../../../bin/filterBam

CHECKBAM:
	@if [ -z "${BAMTOOLS}" ]; then \
		echo '[Makefile]: $${BAMTOOLS} has not been set.'; \
		exit 101; \
	elif [ "${BAMTOOLS}" ]; then \
		echo "filterBam compiled with BAMTOOLS=${BAMTOOLS}"; \
	fi

$(PROGRAM) : $(OBJECTS)
	$(LINK.cc) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) $^ -o $@ $(LIBS)

$(OBJECTS) : %.o:%.cc
	$(LINK.cc) $(CXXFLAGS) $(CPPFLAGS) -c $^ -o $@ $(INCLUDES)


clean:
	rm -f *~ $(OBJECTS)
	rm -f filterBam
	rm -f ../../../bin/filterBam

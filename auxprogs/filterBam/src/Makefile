#	Makefile for filterBam
#
# 	Created: 11-November-2011
#	Last modified: 09-May-2018


include ../../../common.mk

SOURCES = filterBam.cc MatePairs.cc getReferenceName.cc initOptions.cc SingleAlignment.cc \
			printElapsedTime.cc sumMandIOperations.cc sumDandIOperations.cc PairednessCoverage.cc
PROGRAM = filterBam
OBJECTS = $(SOURCES:.cc=.o)
#BAMTOOLS = /usr/include/bamtools
INCLUDES = -I$(BAMTOOL_INC) -Iheaders -I./bamtools
LIBS = $(BAMTOOL_LIBS) -lz
CFLAGS = -std=c++0x
VPATH = functions

all : $(PROGRAM) $(OBJECTS) CHECKBAM BIN

BIN : $(PROGRAM) CHECKBAM
	mv filterBam ../../../bin/filterBam

CHECKBAM:
	@if [ -z "${BAMTOOL}" ]; then \
		echo '[Makefile]: $${BAMTOOL} has not been set.'; \
		exit 101; \
	elif [ "${BAMTOOL}" ]; then \
		echo "filterBam compiled with BAMTOOL=${BAMTOOL}"; \
	fi

$(PROGRAM) : $(OBJECTS)
	$(LINK.cc) $(CFLAGS) $(LDFLAGS) $^ -o $@ $(LIBS) 

$(OBJECTS) : %.o:%.cc
	$(LINK.cc) $(CFLAGS) $(CPPFLAGS) -c $^ -o $@ $(INCLUDES)


clean:
	rm -f *~ $(OBJECTS)
	rm -f ../../../bin/filterBam

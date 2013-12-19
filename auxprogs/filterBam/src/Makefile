#	Makefile for filterBam
#	
# 	Created: 11-November-2011
#	Last modified: 4-March-2012


SOURCES = filterBam.cc MatePairs.cc getReferenceName.cc initOptions.cc SingleAlignment.cc \
			printElapsedTime.cc sumMandIOperations.cc sumDandIOperations.cc PairednessCoverage.cc
PROGRAM = filterBam
OBJECTS = $(SOURCES:.cc=.o)
CFLAGS = $(INCLUDES) -std=c++0x 
INCLUDES = -I$(BAMTOOLS)/include -Iheaders -I$(BAMTOOLS)/src/toolkit
VPATH = functions 
LIBS = $(BAMTOOLS)/lib/libbamtools.a -lz 
INSTDIR = ../bin

all : $(PROGRAM) $(OBJECTS) CHECKBAM BIN

BIN : $(PROGRAM) CHECKBAM
	mkdir -p $(INSTDIR);\
	mv filterBam $(INSTDIR);

CHECKBAM: 
	@if [ -z "${BAMTOOLS}" ]; then \
		echo '[Makefile]: $${BAMTOOLS} has not been set.'; \
		exit 101; \
	elif [ "${BAMTOOLS}" ]; then \
		echo "filterBam compiled with BAMTOOLS=${BAMTOOLS}"; \
	fi

$(PROGRAM) : $(OBJECTS) 
	$(LINK.cc) $(CFLAGS) $^ -o $@ $(LIBS) 

$(OBJECTS) : %.o:%.cc
	$(LINK.cc) $(CFLAGS) -c $^ -o $@ 


clean: 
	rm -f *~ $(OBJECTS) 

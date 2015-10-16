#	Makefile for filterBam
#	
# 	Created: 11-November-2011
#	Last modified: 4-March-2012


SOURCES = filterBam.cc MatePairs.cc getReferenceName.cc initOptions.cc SingleAlignment.cc \
			printElapsedTime.cc sumMandIOperations.cc sumDandIOperations.cc PairednessCoverage.cc
PROGRAM = filterBam
OBJECTS = $(SOURCES:.cc=.o)
BAMTOOLS = /usr/include/bamtools
INCLUDES = -I$(BAMTOOLS) -Iheaders -I./bamtools
LIBS = -lbamtools -lz
CFLAGS = -std=c++0x
VPATH = functions 
INSTDIR = ../bin

all : $(PROGRAM) $(OBJECTS) CHECKBAM BIN

BIN : $(PROGRAM) CHECKBAM
	mkdir -p $(INSTDIR);\
	mv filterBam $(INSTDIR);
	cp $(INSTDIR)/filterBam ../../../bin/filterBam

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
	$(LINK.cc) $(CFLAGS) -c $^ -o $@ $(INCLUDES) 


clean: 
	rm -f *~ $(OBJECTS) 

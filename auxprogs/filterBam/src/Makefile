#	Makefile for filterBam
#


SOURCES = filterBam.cc MatePairs.cc initOptions.cc \
			printElapsedTime.cc sumDandIOperations.cc PairednessCoverage.cc
PROGRAM = filterBam
OBJECTS = $(SOURCES:.cc=.o)
SEQLIB_INSTALL_DIR = /builds/nweeks/repos/SeqLib
HTSLIB_INSTALL_DIR = /builds/nweeks/repos/SeqLib/htslib
INCLUDES = -I$(SEQLIB_INSTALL_DIR)/ -I${HTSLIB_INSTALL_DIR}/ -Iheaders
LDFLAGS=-L${SEQLIB_INSTALL_DIR}/src -L${HTSLIB_INSTALL_DIR}
LIBS = -lseqlib -lz -lhts
CXXFLAGS := -std=c++11 ${CXXFLAGS}
VPATH = functions
LINK.cc ?= g++

all : $(PROGRAM) $(OBJECTS) CHECKBAM BIN

BIN : $(PROGRAM) CHECKBAM
	mkdir -p ../../../bin
	mv filterBam ../../../bin/filterBam

CHECKBAM:
	@if [ -z "${SEQLIB_INSTALL_DIR}" ]; then \
		echo '[Makefile]: $${SEQLIB_INSTALL_DIR} has not been set.'; \
		exit 101; \
	elif [ "${SEQLIB_INSTALL_DIR}" ]; then \
		echo "filterBam compiled with SEQLIB=${SEQLIB_INSTALL_DIR}"; \
	fi

$(PROGRAM) : $(OBJECTS)
	$(LINK.cc) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) $^ -o $@ $(LIBS)

$(OBJECTS) : %.o:%.cc
	$(LINK.cc) $(CXXFLAGS) $(CPPFLAGS) -c $^ -o $@ $(INCLUDES)


clean:
	rm -f *~ $(OBJECTS)
	rm -f filterBam
	rm -f ../../../bin/filterBam

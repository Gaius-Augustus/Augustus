#	Makefile for filterBam
#


SOURCES = filterBam.cc MatePairs.cc getReferenceName.cc initOptions.cc SingleAlignment.cc \
			printElapsedTime.cc sumMandIOperations.cc sumDandIOperations.cc PairednessCoverage.cc
PROGRAM = filterBam
OBJECTS = $(SOURCES:.cc=.o)
ifeq ($(shell uname -s), Darwin) #path for default homebrew installation of bamtools
	BAMTOOLS_INSTALL_DIR = /usr/local/opt/bamtools
else
	BAMTOOLS_INSTALL_DIR = /usr
endif
INCLUDES = -I$(BAMTOOLS_INSTALL_DIR)/include/bamtools -Iheaders -I./bamtools
LIBS = -lbamtools -lz
CXXFLAGS := -std=c++0x ${CXXFLAGS}
VPATH = functions
LINK.cc ?= g++

all : $(PROGRAM) $(OBJECTS) CHECKBAM BIN

BIN : $(PROGRAM) CHECKBAM
	mkdir -p ../../../bin
	mv filterBam ../../../bin/filterBam

CHECKBAM:
	@if [ -z "${BAMTOOLS_INSTALL_DIR}" ]; then \
		echo '[Makefile]: $${BAMTOOLS_INSTALL_DIR} has not been set.'; \
		exit 101; \
	elif [ "${BAMTOOLS_INSTALL_DIR}" ]; then \
		echo "filterBam compiled with BAMTOOLS=${BAMTOOLS_INSTALL_DIR}/include/bamtools"; \
	fi

$(PROGRAM) : $(OBJECTS)
	$(LINK.cc) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) $^ -o $@ $(LIBS)

$(OBJECTS) : %.o:%.cc
	$(LINK.cc) $(CXXFLAGS) $(CPPFLAGS) -c $^ -o $@ $(INCLUDES)


clean:
	rm -f *~ $(OBJECTS)
	rm -f filterBam
	rm -f ../../../bin/filterBam

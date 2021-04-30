#	Makefile for filterBam

CXX ?= g++
CXXFLAGS := -std=c++0x ${CXXFLAGS}

SOURCES = filterBam.cc MatePairs.cc getReferenceName.cc initOptions.cc SingleAlignment.cc \
          printElapsedTime.cc sumMandIOperations.cc sumDandIOperations.cc PairednessCoverage.cc
OBJS = $(SOURCES:.cc=.o)
ifeq ($(shell uname -s), Darwin) #path for default homebrew installation of bamtools
	BAMTOOLS_INSTALL_DIR = /usr/local/opt/bamtools
else
	BAMTOOLS_INSTALL_DIR = /usr
endif
INCLS += -I$(BAMTOOLS_INSTALL_DIR)/include/bamtools -I./headers -I./bamtools
LDLIBS += -lbamtools -lz
VPATH = functions

.PHONY: all BIN CHECKBAM clean

all: CHECKBAM filterBam $(OBJS) BIN

BIN: CHECKBAM filterBam
	mkdir -p ../../../bin
	mv filterBam ../../../bin/filterBam

CHECKBAM:
	@if [ -z "${BAMTOOLS_INSTALL_DIR}" ]; then \
		echo '[Makefile]: $$(BAMTOOLS_INSTALL_DIR) has not been set.'; \
		exit 101; \
	elif [ "${BAMTOOLS_INSTALL_DIR}" ]; then \
		echo "filterBam compiled with BAMTOOLS=$(BAMTOOLS_INSTALL_DIR)/include/bamtools"; \
	fi

filterBam: $(OBJS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@  $^ $(LDLIBS)

$(OBJS): %.o: %.cc
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLS) -c $^


clean:
	rm -f *~ $(OBJS)
	rm -f filterBam
	rm -f ../../../bin/filterBam

#	Makefile for filterBam

include ../../../common.mk

CXX ?= g++
CXXFLAGS := -std=c++0x ${CXXFLAGS}

SOURCES = filterBam.cc MatePairs.cc getReferenceName.cc initOptions.cc SingleAlignment.cc \
          printElapsedTime.cc sumMandIOperations.cc sumDandIOperations.cc PairednessCoverage.cc
OBJS = $(SOURCES:.cc=.o)

# keep compatibility with old versions still using BAMTOOLS_INSTALL_DIR
BAMTOOLS_INSTALL_DIR ?= /usr

# sets the default locations if it has not already been set in common.mk
INCLUDE_PATH_BAMTOOLS ?= -I$(BAMTOOLS_INSTALL_DIR)/include/bamtools
LIBRARY_PATH_BAMTOOLS ?= -L$(BAMTOOLS_INSTALL_DIR)/lib -Wl,-rpath,$(BAMTOOLS_INSTALL_DIR)/lib

INCLS += -I./headers -I./bamtools
INCLS += $(INCLUDE_PATH_BAMTOOLS)
LDFLAGS += $(LIBRARY_PATH_BAMTOOLS)
LDLIBS += -lbamtools
INCLS += $(INCLUDE_PATH_ZLIB)
LDFLAGS += $(LIBRARY_PATH_ZLIB)
LDLIBS += -lz
VPATH = functions

.PHONY: all clean

all: filterBam
	mkdir -p ../../../bin
	mv filterBam ../../../bin/filterBam

filterBam: $(OBJS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@  $^ $(LDLIBS)

$(OBJS): %.o: %.cc
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLS) -c $^


clean:
	rm -f *~ $(OBJS)
	rm -f filterBam
	rm -f ../../../bin/filterBam

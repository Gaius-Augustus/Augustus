#	Makefile corresponding to BAM2HINTS
#
#	NOTE:
#	a) "-g" option allows for debugging
#	b) libbamtools.a and lz library are used during the linking process
#

include ../../common.mk

CXX ?= g++
CXXFLAGS := -Wall -O3  ${CXXFLAGS} # -g -p -g -ggdb

# keep compatibility with old versions still using BAMTOOLS_INSTALL_DIR
BAMTOOLS_INSTALL_DIR ?= /usr

# sets the default locations if it has not already been set in common.mk
INCLUDE_PATH_BAMTOOLS ?= -I$(BAMTOOLS_INSTALL_DIR)/include/bamtools
LIBRARY_PATH_BAMTOOLS ?= -L$(BAMTOOLS_INSTALL_DIR)/lib -Wl,-rpath,$(BAMTOOLS_INSTALL_DIR)/lib

INCLS += $(INCLUDE_PATH_BAMTOOLS)
LDFLAGS += $(LIBRARY_PATH_BAMTOOLS)
LDLIBS += -lbamtools
INCLS += $(INCLUDE_PATH_ZLIB)
LDFLAGS += $(LIBRARY_PATH_ZLIB)
LDLIBS += -lz

.PHONY: all clean

all: bam2hints

bam2hints: bam2hints.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) $(INCLS) -o $@ $^ $(LDLIBS)
	mkdir -p ../../bin
	cp -f bam2hints ../../bin

bam2hints.o: bam2hints.cc
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLS) -c $^

clean:
	rm -f bam2hints.o bam2hints ../../bin/bam2hints

#	Makefile corresponding to BAM2HINTS
#
#	NOTE:
#	a) "-g" option allows for debugging
#	b) libbamtools.a and lz library are used during the linking process
#
CXX ?= g++
CXXFLAGS := -Wall -O3  ${CXXFLAGS} # -g -p -g -ggdb

# Variable definition
ifeq ($(shell uname -s), Darwin) #path for default homebrew installation of bamtools
	BAMTOOLS_INSTALL_DIR = /usr/local/opt/bamtools
else
	BAMTOOLS_INSTALL_DIR = /usr
endif

INCLS += -I$(BAMTOOLS_INSTALL_DIR)/include/bamtools
LDLIBS += -lbamtools -lz

.PHONY: all clean

all: bam2hints

bam2hints: bam2hints.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)
	mkdir -p ../../bin
	cp -f bam2hints ../../bin

bam2hints.o: bam2hints.cc
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLS) -c $^

clean:
	rm -f bam2hints.o bam2hints ../../bin/bam2hints

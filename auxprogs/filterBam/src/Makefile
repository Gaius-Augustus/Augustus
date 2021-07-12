#	Makefile for filterBam

include ../../../common.mk

CXX ?= g++
CXXFLAGS := -std=c++0x ${CXXFLAGS}

SOURCES = filterBam.cc MatePairs.cc initOptions.cc \
          printElapsedTime.cc sumDandIOperations.cc PairednessCoverage.cc
OBJS = $(SOURCES:.cc=.o)

INCLS += -I./headers
INCLS += $(shell pkg-config --cflags  jsoncpp htslib libseqlib)
LDFLAGS += $(shell pkg-config --ldflags  jsoncpp htslib libseqlib)
LDLIBS += -lpthread $(shell pkg-config --libs jsoncpp htslib libseqlib) -lz
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

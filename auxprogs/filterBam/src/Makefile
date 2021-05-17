#	Makefile for filterBam
#


SOURCES = filterBam.cc MatePairs.cc initOptions.cc \
			printElapsedTime.cc sumDandIOperations.cc PairednessCoverage.cc
PROGRAM = filterBam
OBJECTS = $(SOURCES:.cc=.o)
PKG_CONFIG=$(shell command -v pkg-config 2>/dev/null)
ifneq ($(PKG_CONFIG),)
CXXFLAGS += -Iheaders $(shell pkg-config --cflags  jsoncpp htslib libseqlib)
LIBS += -lpthread $(shell pkg-config --libs jsoncpp htslib libseqlib)
else
CXXFLAGS += -Iheaders
LIBS += -ljsoncpp -lhts -lseqlib -lpthread
endif
CXXFLAGS += -std=c++11 -O2
VPATH = functions
LINK.cc ?= g++

all : $(PROGRAM) $(OBJECTS) BIN

BIN : $(PROGRAM)
	mkdir -p ../../../bin
	mv filterBam ../../../bin/filterBam

$(PROGRAM) : $(OBJECTS)
	$(LINK.cc) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) $^ -o $@ $(LIBS)

$(OBJECTS) : %.o:%.cc
	$(LINK.cc) $(CXXFLAGS) $(CPPFLAGS) -c $^ -o $@ $(INCLUDES)


clean:
	rm -f *~ $(OBJECTS)
	rm -f filterBam
	rm -f ../../../bin/filterBam

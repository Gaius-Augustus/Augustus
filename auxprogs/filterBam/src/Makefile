#	Makefile for filterBam

include ../../../common.mk

# set to true if the SeqLib library should used to handle bam files
#USE_SEQLIB ?= true
# set to true if the seqlib library version is 1.1.1 (various bugs and api changes are taken into account)
#SEQLIB_1_1_1 ?= true
# set to true if the seqlib library version is 1.1.2 (various bugs and api changes are taken into account)
#SEQLIB_1_1_2 ?= true

CXX ?= g++
CXXFLAGS := -Wall -Wextra -Wpedantic -std=c++11 ${CXXFLAGS}

SOURCES = filterBam.cc initOptions.cc printElapsedTime.cc MatePairs.cc PairednessCoverage.cc
INCLS += -I./headers

ifneq (,$(findstring $(USE_SEQLIB),1 true True TRUE)) # if USE_SEQLIB is defined and contains 1, true, True or TRUE
	# use SeqLib library
	SOURCES += BamSeqLibAccess.cc
	# default locations:
	INCLUDE_PATH_SEQLIB ?= -I /usr/include/SeqLib -I/usr/include/htslib -I/usr/include/jsoncpp -I/usr/local/include/SeqLib -I/usr/local/include/htslib -I/usr/include/jsoncpp
	INCLS += $(INCLUDE_PATH_SEQLIB)
	LDFLAGS += $(LIBRARY_PATH_SEQLIB)
	LDLIBS += -lseqlib -lfml -lbwa -lhts -lm -llzma -lpthread
	CPPFLAGS += -DUSE_SEQLIB
	
	ifneq (,$(findstring $(SEQLIB_1_1_1),1 true True TRUE)) # if SEQLIB_1_1_1 is defined and contains 1, true, True or TRUE
		LDLIBS += -lssw -ljsoncpp
		CPPFLAGS += -DSEQLIB_1_1_1
	else ifneq (,$(findstring $(SEQLIB_1_1_2),1 true True TRUE)) # if SEQLIB_1_1_2 is defined and contains 1, true, True or TRUE
		LDLIBS += -lssw -ljsoncpp
		CPPFLAGS += -DSEQLIB_1_1_2
	else
		LDLIBS += -lbz2
		CPPFLAGS += -DSEQLIB_1_2
	endif
else # use BamTools library
	SOURCES += BamToolsAccess.cc
	# keep compatibility with old versions still using BAMTOOLS_INSTALL_DIR
	BAMTOOLS_INSTALL_DIR ?= /usr

	# sets the default locations if it has not already been set in common.mk
	INCLUDE_PATH_BAMTOOLS ?= -I$(BAMTOOLS_INSTALL_DIR)/include/bamtools
	LIBRARY_PATH_BAMTOOLS ?= -L$(BAMTOOLS_INSTALL_DIR)/lib -Wl,-rpath,$(BAMTOOLS_INSTALL_DIR)/lib

	INCLS += $(INCLUDE_PATH_BAMTOOLS)
	LDFLAGS += $(LIBRARY_PATH_BAMTOOLS)
	LDLIBS += -lbamtools
endif

OBJS = $(SOURCES:.cc=.o)
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

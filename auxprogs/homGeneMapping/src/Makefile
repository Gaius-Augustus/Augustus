#
# Makefile for homGeneMapping

include ../../../common.mk

# set SQLITE = false in common.mk to disable the usage of the SQLite library, option --dbaccess will not be available
# set BOOST = false in common.mk to disable the usage of the Boost graph library, the option --printHomologs will not be available
BOOST ?= true

ifeq (,$(findstring $(COMPGENEPRED),0 false False FALSE)) # if COMPGENEPRED is not defined or not contains 0, false, False or FALSE
# comment this flag to disable database access for retrieval of hints, the option --dbaccess will not be available
	SQLITE ?= true
endif

CXX   ?= g++

# Notes: - "-Wno-sign-compare" eliminates a high number of warnings (see footnote below). Please adopt
#          a strict signed-only usage strategy to avoid mistakes since we are not warned about this.
#        - The order of object files in $(OBJS) IS IMPORTANT (see lldouble.hh)
CXXFLAGS := -Wall -Wno-sign-compare -ansi -pedantic -std=c++0x -pthread -O2 ${CXXFLAGS}

INCLS	= -I../include
LIBS	=
OBJS	= gene.o genome.o

ifeq (,$(findstring $(BOOST),0 false False FALSE)) # if BOOST is not defined or not contains 0, false, False or FALSE
	CPPFLAGS += -DBOOST
	INCLS +=  # add the boost include path here, if boost ist not installed system-wide
endif

ifneq (,$(findstring $(SQLITE),1 true True TRUE)) # if SQLITE is defined and contains 1, true, True or TRUE
	LIBS    += -lsqlite3 # add the sqlite library path here, if sqlite is not install system-wide
	OBJS    += sqliteDB.o
	CPPFLAGS  += -DM_SQLITE
endif

all: homGeneMapping

.SUFFIXES:
.SUFFIXES: .cc .o .so

.cc.o:
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) -o $@ $< $(INCLS)

homGeneMapping: main.cc $(OBJS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) -o $@ $^ $(INCLS) $(LIBS)
	mkdir -p ../../../bin/
	cp homGeneMapping ../../../bin/homGeneMapping

clean:
	rm -f homGeneMapping $(OBJS)
	rm -f ../../../bin/homGeneMapping

main.o : \
	../include/genome.hh \
	../include/projectio.hh \
	main.cc

gene.o : \
	../include/gene.hh \
	../include/projectio.hh \
	gene.cc

genome.o : \
	../include/genome.hh \
	../include/bitmasking.hh \
	../include/gene.hh \
	genome.cc

sqliteDB.o: \
        ../include/sqliteDB.hh \
        sqliteDB.cc

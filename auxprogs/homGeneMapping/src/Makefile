#
# Makefile for homGeneMapping
#
# comment this flag to disable the usage of the Boost graph library, the option --printHomologs will not be available
BOOST = true

CC      = g++

# Notes: - "-Wno-sign-compare" eliminates a high number of warnings (see footnote below). Please adopt
#          a strict signed-only usage strategy to avoid mistakes since we are not warned about this.
#        - The order of object files in $(OBJS) IS IMPORTANT (see lldouble.hh)
CFLAGS := -Wall -Wno-sign-compare -ansi -pedantic -std=c++0x -pthread -O2 ${CFLAGS} # -DDEBUG -g -ggdb -pg

INCLS	= -I../include
LIBS	= -lsqlite3 # add the sqlite library path here, if sqlite is not install system-wide
OBJS	= gene.o genome.o sqliteDB.o

ifdef BOOST
	CFLAGS += -DBOOST
	INCLS +=  # add the boost include path here, if boost ist not installed system-wide
endif

all: homGeneMapping

.SUFFIXES:
.SUFFIXES: .cc .o .so

.cc.o:
	$(CC) -c $(CFLAGS) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) -o $@ $< $(INCLS)

homGeneMapping: main.cc $(OBJS)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $^ $(INCLS)  $(LIBS)
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

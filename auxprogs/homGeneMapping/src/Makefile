#
# Makefile for homGeneMapping
#
# BOOST = true

# set this flag to enable
# database access for retrieval of hints
# SQLITE = true

CC      = g++

# Notes: - "-Wno-sign-compare" eliminates a high number of warnings (see footnote below). Please adopt
#          a strict signed-only usage strategy to avoid mistakes since we are not warned about this.
#        - The order of object files in $(OBJS) IS IMPORTANT (see lldouble.hh)
CFLAGS := -Wall -Wno-sign-compare -ansi -pedantic -std=c++0x -pthread -O2 ${CFLAGS} # -DDEBUG -g -ggdb -pg

INCLS	= -I../include
LIBS	=
OBJS	= gene.o genome.o

ifdef BOOST
	CFLAGS += -DBOOST
	INCLS +=  # add the boost include path here, if boost ist not installed system-wide
endif

ifdef SQLITE
INCLS	+=           # add the sqlite include path here, if sqlite is not installed system-wide
LIBS    += -lsqlite3 # add the sqlite library path here, if sqlite is not install system-wide
OBJS    += sqliteDB.o
CFLAGS  += -DSQLITE
endif

all: homGeneMapping

.SUFFIXES:
.SUFFIXES: .cc .o .so

.cc.o:
	$(CC) -c $(CFLAGS) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) -o $@ $< $(INCLS)

homGeneMapping: main.cc $(OBJS)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o $@ $^ $(INCLS)  $(LIBS)
	cp homGeneMapping ../../../bin/homGeneMapping

clean:
	rm -f homGeneMapping $(OBJS)

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

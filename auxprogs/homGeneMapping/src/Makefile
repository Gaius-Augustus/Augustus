#
# Makefile for homGeneMapping
#
# BOOST = true

CC      = g++

# Notes: - "-Wno-sign-compare" eliminates a high number of warnings (see footnote below). Please adopt 
#          a strict signed-only usage strategy to avoid mistakes since we are not warned about this.
#        - The order of object files in $(OBJS) IS IMPORTANT (see lldouble.hh)
CFLAGS := -Wall -Wno-sign-compare -ansi -pedantic -std=c++0x -pthread -g -ggdb -O2 ${CFLAGS} # -DDEBUG -g -ggdb -pg

INCLS	= -I../include

OBJS	= gene.o genome.o 

ifdef BOOST
	CFLAGS += -DBOOST
	INCLS +=  # add the boost include path here, if boost ist not installed system-wide
endif

all: homGeneMapping

.SUFFIXES:
.SUFFIXES: .cc .o .so

.cc.o:
	$(CC) -c $(CFLAGS) -o $@ $< $(INCLS)

homGeneMapping: main.cc $(OBJS)
	 $(CC) $(CFLAGS) -o $@ $^ $(INCLS)
	mv homGeneMapping ../bin/
	cp ../bin/homGeneMapping ../../../bin/homGeneMapping
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

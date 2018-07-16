# Makefile of bam2wig
#
# NOTE: Modify the variable TOOLDIR according to where your samtools/htslib/bcftools/tabix software has been installed
#
# Created: 12-June-2012
# Last modified: 10-NOV-2017
#

ifndef TOOLDIR
	TOOLDIR=$(HOME)/tools
# replace with your own parent directory of samtools, htslib, bcftools, tabix
# if TOOLDIR is not specified as environment variable
endif

PROGRAM = bam2wig
SOURCES = $(PROGRAM)
OBJECTS = $(SOURCES:.c=.o)
SAMTOOLS=$(TOOLDIR)/samtools/
HTSLIB=$(TOOLDIR)/htslib/
BCFTOOLS=$(TOOLDIR)/bcftools/
TABIX=$(TOOLDIR)/tabix/
INCLUDES=-I$(SAMTOOLS) -I. -I$(HTSLIB) -I$(BCFTOOLS) -I$(TABIX)
VPATH=$(SAMTOOLS)
LIBS=$(SAMTOOLS)/libbam.a $(HTSLIB)/libhts.a -lcurses -lm -lz -lpthread -lcurl -lssl -lcrypto
CFLAGS=-Wall -O2 $(INCLUDES)
CC=gcc

$(PROGRAM) : bam2wig.o
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS) -lbz2 -llzma
	cp bam2wig ../../bin/bam2wig

bam2wig.o : bam2wig.c
	$(CC) $(CFLAGS) -c $^ -o $@

clean:
	rm -f $(OBJECTS) bam2wig.o
	rm -f ../../bin/bam2wig

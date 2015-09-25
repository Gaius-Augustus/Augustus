# Makefile of bam2wig
#
# NOTE: Modify the variable SAMTOOLS according to where your samtools software has been installed
#
# Created: 12-June-2012
# Last modified: 07-April-2014
#

PROGRAM = bam2wig
SOURCES = $(PROGRAM) bam_index.o 
OBJECTS = $(SOURCES:.c=.o)
TOOLDIR=/home/mario/tools
SAMTOOLS=$(TOOLDIR)/samtools/
HTSLIB=$(TOOLDIR)/htslib/
BCFTOOLS=$(TOOLDIR)/bcftools/
TABIX=$(TOOLDIR)/tabix/
INCLUDES=-I$(SAMTOOLS) -I. -I$(HTSLIB) -I$(BCFTOOLS) -I$(TABIX)
VPATH=$(SAMTOOLS)
LIBS=$(SAMTOOLS)/libbam.a $(HTSLIB)/libhts.a -lcurses -lm -lz -lpthread
CFLAGS= -g -Wall -O2 $(INCLUDES)
CC=gcc

$(PROGRAM) : bam2wig.o bam_index.o
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

bam2wig.o : bam2wig.c 
	$(CC) $(CFLAGS) -c $^ -o $@ $(INCLUDES)

clean: 
	rm -f $(OBJECTS) bam2wig.o

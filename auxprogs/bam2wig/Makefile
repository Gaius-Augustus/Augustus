# Makefile of bam2wig
#
# NOTE: Modify the variable SAMTOOLS according to where your samtools software has been installed
#
# Created: 12-June-2012
# Last modified: 19-June-2012
#

PROGRAM = bam2wig
SOURCES = $(PROGRAM) bam_index.o 
OBJECTS = $(SOURCES:.c=.o) 
SAMTOOLS=/fs1.cm/katharina/tools/samtools/
INCLUDES=-I$(SAMTOOLS) -I.
VPATH=$(SAMTOOLS)
LIBS=$(SAMTOOLS)/libbam.a -lcurses -lm -lz -lpthread
CFLAGS= -g -Wall -O2 $(INCLUDES)
CC=gcc

$(PROGRAM) : bam2wig.o bam_index.o
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

bam2wig.o : bam2wig.c 
	$(CC) $(CFLAGS) -c $^ -o $@ $(INCLUDES)

clean: 
	rm -f $(OBJECTS) bam2wig.o


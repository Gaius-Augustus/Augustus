# Makefile of bam2wig
#
# NOTE: Modify the variable SAMTOOLS according to where your samtools software has been installed
#
# Created: 12-June-2012
# Last modified: 13-June-2012
#

PROGRAM = bam2wig
SOURCES = $(PROGRAM)  		
OBJECTS = $(SOURCES:.c=.o) 
CFLAGS=-g -Wall -O2
SAMTOOLS=/home/tonatiuh/Documents/git/samtools
INCLUDES=-I$(SAMTOOLS) -I.
VPATH=$(SAMTOOLS)
LIBS=$(SAMTOOLS)/libbam.a -lcurses -lm -lz 
CC=gcc

# RECIPES (just to remember how compiling and linking are done)
all : bam2wig.o bam_index.o
	$(CC) $(CFLAGS) bam2wig.o $(SAMTOOLS)/bam_index.o -o bam2wig $(LIBS)

bam2wig.o : bam2wig.c 
	$(CC) $(CFLAGS) -c bam2wig.c -o bam2wig.o -I$(SAMTOOLS) -I.


clean :
	rm -f *~ $(OBJECTS) *.o bam2wig
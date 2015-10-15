# Makefile of bam_index.o

# Created: October 15th 2015
# Last modified: October 15th by Katharina J. Hoff

SAMTOOLS=/usr/include/samtools
HTSLIB=/usr/lib/x86_64-linux-gnu

INCLUDES=-I$(SAMTOOLS) -I. -I$(HTSLIB)
LIBS=$(HTSLIB)/libhts.a -lcurses -lm -lz -lpthread
CFLAGS= -g -Wall -O2 $(INCLUDES)
CC=gcc

bam_index.o : bam_index.c
	$(CC) $(CFLAGS) -c $^ -o $@ $(INCLUDES)

clean:
	rm -f bam_index.o

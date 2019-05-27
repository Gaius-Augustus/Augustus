# Makefile of bam2wig
#
# NOTE: Modify the variable TOOLDIR according to where your samtools/htslib/bcftools software has been installed
# 
# License: Artistic License, see file LICENSE.TXT or 
#          https://opensource.org/licenses/artistic-license-1.0
#

include ../../common.mk

#ifndef TOOLDIR
#	TOOLDIR=$(HOME)/tools
# replace with your own parent directory of samtools, htslib, bcftools
# if TOOLDIR is not specified as environment variable
#endif

PROGRAM = bam2wig
SOURCES = $(PROGRAM)
OBJECTS = $(SOURCES:.c=.o)
INCLUDES=-I$(SAMTOOLS_INC) -I. -I$(HTSLIB_INC) -I$(BCFTOOLS_INC)
# VPATH=$(SAMTOOLS)
LIBS=-lcurses -lm -lz -lpthread -lcurl -lssl -lcrypto $(SAMTOOLS_LIBS) $(HTSLIB_LIBS)
CFLAGS=-Wall -O2 $(INCLUDES)
CC=gcc

$(PROGRAM) : bam2wig.o
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS) -lbz2 -llzma
	mkdir -p ../../bin
	cp bam2wig ../../bin/bam2wig

bam2wig.o : bam2wig.c
	$(CC) $(CFLAGS) -c $^ -o $@

clean:
	rm -f $(OBJECTS) bam2wig.o
	rm -f ../../bin/bam2wig

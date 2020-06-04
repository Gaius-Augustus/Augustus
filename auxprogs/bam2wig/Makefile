# Makefile of bam2wig
#
# NOTE: Modify the variable TOOLDIR according to where your htslib software has been installed
# 
# License: Artistic License, see file LICENSE.TXT or 
#          https://opensource.org/licenses/artistic-license-1.0
#

ifndef TOOLDIR
	TOOLDIR=$(HOME)/tools
# replace with your own parent directory of htslib
# if TOOLDIR is not specified as environment variable
endif

PROGRAM = bam2wig
SOURCES = $(PROGRAM)
OBJECTS = $(SOURCES:.c=.o)
HTSLIB=$(TOOLDIR)/htslib
INCLUDES=-I. -I$(HTSLIB) -lbz2 -llzma
HTSLIBS?=$(HTSLIB)/libhts.a -lm -lz -lpthread -lcurl -lssl -lcrypto
CFLAGS:=-Wall -O2 $(INCLUDES) $(CFLAGS)
CC?=gcc

$(PROGRAM) : bam2wig.o
	$(CC) $(LDFLAGS) $^ -o $@ $(HTSLIBS)
	mkdir -p ../../bin
	cp bam2wig ../../bin/bam2wig

bam2wig.o : bam2wig.c
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $^ -o $@

clean:
	rm -f $(OBJECTS) bam2wig.o bam2wig
	rm -f ../../bin/bam2wig

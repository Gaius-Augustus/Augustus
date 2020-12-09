# Makefile of bam2wig
#
# NOTE: install libhts-dev >= 1.10 or modify the variable TOOLDIR according to where your htslib software has been installed
# 
# License: Artistic License, see file LICENSE.TXT or 
#          https://opensource.org/licenses/artistic-license-1.0
#

ifdef TOOLDIR
	DEF_TOOLDIR=1
else
	TOOLDIR=$(HOME)/tools
	# replace with your own parent directory of htslib
	# if TOOLDIR is not specified as environment variable
	# and libhts-dev is not installed
endif

PROGRAM = bam2wig
SOURCES = $(PROGRAM)
OBJECTS = $(SOURCES:.c=.o)

ifneq ($(wildcard ${TOOLDIR}/htslib/htslib/.),) # if TOOLDIR exists and contains htslib
	INCLUDES=-I$(TOOLDIR)/htslib -I$(TOOLDIR)/htslib/htslib
	HTSLIBS=-L$(TOOLDIR)/htslib -Wl,-rpath,$(TOOLDIR)/htslib -lhts -lcurses -lm -lz -lpthread -lcurl -lssl -lcrypto -lbz2 -llzma
	HTSERRHINT="Check if the installed HTSlib is of version 1.10 or higher - see README.txt"
else
	INCLUDES=-I/usr/include/htslib
	HTSLIBS=-lhts -lz -lpthread
	ifdef DEF_TOOLDIR
		HTSERRHINT="There is no htslib folder in directory $(TOOLDIR) \n"
	endif
	HTSERRHINT:=$(HTSERRHINT)"Check if HTSlib is installed and of version 1.10 or higher - see README.txt"
endif

CFLAGS:=-Wall -O2 $(CFLAGS)
CC?=gcc

.PHONY : clean test	

$(PROGRAM) : bam2wig.o
	$(CC) $(CFLAGS)	$(CPPFLAGS) $(LDFLAGS) $^ -o $@ $(HTSLIBS)
	mkdir -p ../../bin
	cp bam2wig ../../bin/bam2wig

bam2wig.o : bam2wig.c
	-$(CC) $(CFLAGS) $(CPPFLAGS) -c $^ -o $@ $(INCLUDES)
	@if [ ! -f $@ ]; then echo $(HTSERRHINT); exit 1; fi

clean:
	rm -f $(OBJECTS) bam2wig.o bam2wig
	rm -f ../../bin/bam2wig
	cd ../..; tests/auxprogs/bam2wig/test_bam2wig.py --clean
	
test:   $(PROGRAM)
	cd ../..; tests/auxprogs/bam2wig/test_bam2wig.py --compare --html

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

ifneq ($(wildcard ${TOOLDIR}/include/htslib/.),) # if TOOLDIR exists and contains htslib
	INCLUDES=-I$(TOOLDIR)/include/htslib
	HTSLIBS=-L$(TOOLDIR)/lib -Wl,-rpath,$(TOOLDIR)/lib -lhts -lcurses -lm -lz -lpthread -lcurl -lssl -lcrypto -lbz2 -llzma
	HTSERRHINT="Check if the installed HTSlib is of version 1.10 or higher - see README.md"
else
    ifeq ($(shell uname -s), Darwin) #path for default homebrew installation of htslib
		INCLUDES=-I/usr/local/opt/htslib/include/htslib
    else
		INCLUDES=-I/usr/include/htslib
    endif
	HTSLIBS=-lhts -lz -lpthread
	ifdef DEF_TOOLDIR
		HTSERRHINT="There is no htslib folder in directory $(TOOLDIR)/include \n"
	endif
	HTSERRHINT:=$(HTSERRHINT)"Check if HTSlib is installed and of version 1.10 or higher - see README.md"
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
	if [ ! -z $(shell which python3) ] ; then cd ../../tests/short && ./execute_test.py --clean bam2wig ; fi

test:   $(PROGRAM)
	cd ../../tests/short && ./execute_test.py --compare --html bam2wig

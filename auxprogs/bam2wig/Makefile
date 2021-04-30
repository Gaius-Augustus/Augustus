# Makefile of bam2wig
#
# NOTE: install libhts-dev >= 1.10 or modify the variable TOOLDIR according to where your htslib software has been installed
#
# License: Artistic License, see file LICENSE.TXT or
#          https://opensource.org/licenses/artistic-license-1.0
#

CC ?= gcc
CFLAGS := -Wall -O2 $(CFLAGS)

ifdef TOOLDIR
	DEF_TOOLDIR=1
else
	TOOLDIR=$(HOME)/tools
	# replace with your own parent directory of htslib
	# if TOOLDIR is not specified as environment variable
	# and libhts-dev is not installed
endif

ifneq ($(wildcard ${TOOLDIR}/include/htslib/.),) # if TOOLDIR exists and contains htslib
	INCLS += -I$(TOOLDIR)/include/htslib
	LDLIBS += -L$(TOOLDIR)/lib -Wl,-rpath,$(TOOLDIR)/lib -lhts -lcurses -lm -lz -lpthread -lcurl -lssl -lcrypto -lbz2 -llzma
	HTSERRHINT="Check if the installed HTSlib is of version 1.10 or higher - see README.md"
else
    ifeq ($(shell uname -s), Darwin) #path for default homebrew installation of htslib
		INCLS += -I/usr/local/opt/htslib/include/htslib
    else
		INCLS += -I/usr/include/htslib
    endif
	LDLIBS += -lhts -lz -lpthread
	ifdef DEF_TOOLDIR
		HTSERRHINT="There is no htslib folder in directory $(TOOLDIR)/include \n"
	endif
	HTSERRHINT:=$(HTSERRHINT)"Check if HTSlib is installed and of version 1.10 or higher - see README.md"
endif

.PHONY: all clean test

all: bam2wig

bam2wig: bam2wig.o
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)
	mkdir -p ../../bin
	cp -f bam2wig ../../bin/bam2wig

bam2wig.o: bam2wig.c
	-$(CC) $(CPPFLAGS) $(CFLAGS) $(INCLS) -c $^
	@if [ ! -f $@ ]; then echo $(HTSERRHINT); exit 1; fi

clean:
	rm -f bam2wig.o bam2wig
	rm -f ../../bin/bam2wig
	if [ ! -z $(shell which python3) ] ; then cd ../../tests/short && ./execute_test.py --clean bam2wig ; fi

test:   bam2wig
	cd ../../tests/short && ./execute_test.py --compare --html bam2wig

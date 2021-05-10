# Makefile of bam2wig
#
# NOTE: install libhts-dev >= 1.10
#  or modify the variable HTSLIB_INSTALL_DIR according to where your htslib software has been installed
#  or set variables INCLUDE_PATH_HTSLIB and LIBRARY_PATH_HTSLIB that point to hts header files and hts library files respectively.
#
# License: Artistic License, see file LICENSE.TXT or
#          https://opensource.org/licenses/artistic-license-1.0
#

include ../../common.mk

CC ?= gcc
CFLAGS := -Wall -O2 $(CFLAGS)
	
HTSERRHINT=""
ifndef INCLUDE_PATH_HTSLIB
	# set INCLUDE_PATH_HTSLIB and LIBRARY_PATH_HTSLIB
	# or HTSLIB_INSTALL_DIR if neither of these is specified as environment variable
	# and libhts-dev is not installed
	# keep compatibility with old versions still using TOOLDIR
  ifdef TOOLDIR
	HTSLIB_INSTALL_DIR = $(TOOLDIR)/htslib
  endif
	HTSLIB_INSTALL_DIR ?= $(HOME)/tools/htslib
  ifneq ($(wildcard ${HTSLIB_INSTALL_DIR}/include/htslib/.),) # if HTSLIB_INSTALL_DIR exists and contains include/htslib
	INCLUDE_PATH_HTSLIB = -I$(HTSLIB_INSTALL_DIR)/include/htslib
	LIBRARY_PATH_HTSLIB = -L$(HTSLIB_INSTALL_DIR)/lib -Wl,-rpath,$(HTSLIB_INSTALL_DIR)/lib
  else
    ifneq ($(wildcard ${HTSLIB_INSTALL_DIR}/htslib/.),) # if HTSLIB_INSTALL_DIR exists and contains htslib
	INCLUDE_PATH_HTSLIB = -I$(HTSLIB_INSTALL_DIR)/htslib
	LIBRARY_PATH_HTSLIB = -L$(HTSLIB_INSTALL_DIR) -Wl,-rpath,$(HTSLIB_INSTALL_DIR)
    else
	HTSERRHINT="There is no htslib folder in directory $(HTSLIB_INSTALL_DIR) - "
    endif
  endif
endif

# set default include paths:
INCLUDE_PATH_HTSLIB ?= -I/usr/include/htslib -I/usr/local/include/htslib/
HTSERRHINT:=$(HTSERRHINT)"Check if HTSlib is installed and of version 1.10 or higher - see README.md"

INCLS += $(INCLUDE_PATH_HTSLIB)
LDFLAGS += $(LIBRARY_PATH_HTSLIB)
LDLIBS += -lhts
INCLS += $(INCLUDE_PATH_ZLIB)
LDFLAGS += $(LIBRARY_PATH_ZLIB)
LDLIBS += -lz -lpthread
# still needed? : LDLIBS +=  -lm -lcurses -lcurl -lssl -lcrypto -lbz2 -llzma

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

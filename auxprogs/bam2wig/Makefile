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

all: bam2wig

bam2wig: bam2wig
	cp -f bam2wig ../../bin/bam2wig

test:   bam2wig
	cd ../../tests/short && ./execute_test.py --compare --html bam2wig

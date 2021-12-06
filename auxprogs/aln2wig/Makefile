# Makefile of aln2wig

CC ?= gcc
CFLAGS := -Wall -Wno-sign-compare -ansi -pedantic -O2 $(CFLAGS)

.PHONY: all clean

aln2wig: aln2wig.o
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) -o aln2wig $^
	mkdir -p ../../bin
	cp -f aln2wig ../../bin

aln2wig.o: aln2wig.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $^

all: aln2wig

clean:
	rm -f aln2wig aln2wig.o ../../bin/aln2wig

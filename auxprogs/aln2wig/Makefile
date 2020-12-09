CFLAGS := -Wall -Wno-sign-compare -ansi -pedantic -O2 ${CFLAGS}
CC?=gcc

aln2wig : aln2wig.o
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) -o aln2wig aln2wig.o
	mkdir -p ../../bin
	cp aln2wig ../../bin
aln2wig.o : aln2wig.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c aln2wig.c

all : aln2wig

clean :
	rm -f aln2wig aln2wig.o; rm -f ../../bin/aln2wig 

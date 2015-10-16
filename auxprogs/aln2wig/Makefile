CFLAGS := -Wall -Wno-sign-compare -ansi -pedantic -O2 -ggdb

psl2wig : aln2wig.o
	gcc $(CFLAGS) -o aln2wig aln2wig.o;
	cp aln2wig ../../bin
psl2wig.o : aln2wig.c
	gcc $(CFLAGS) -c aln2wig.c

all : psl2wig

clean : 
	rm -f aln2wig aln2wig.o 

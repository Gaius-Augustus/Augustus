CFLAGS := -Wall -Wno-sign-compare -ansi -pedantic -O2 -ggdb ${CFLAGS}

psl2wig : aln2wig.o
	gcc $(CFLAGS) -o aln2wig aln2wig.o;
	cp aln2wig ../../bin
psl2wig.o : aln2wig.c
	gcc $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) -c aln2wig.c

all : psl2wig

clean :
	rm -f aln2wig aln2wig.o; rm -f ../../bin/aln2wig 

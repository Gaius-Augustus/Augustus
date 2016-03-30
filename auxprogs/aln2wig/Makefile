CFLAGS := -Wall -Wno-sign-compare -ansi -pedantic -O2 -ggdb ${CFLAGS}

aln2wig : aln2wig.o
	gcc $(CFLAGS) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o aln2wig aln2wig.o;
	cp aln2wig ../../bin
aln2wig.o : aln2wig.c
	gcc $(CPPFLAGS) $(CXXFLAGS) $(CFLAGS) $(LDFLAGS) -c aln2wig.c

all : aln2wig

clean :
	rm -f aln2wig aln2wig.o; rm -f ../../bin/aln2wig 

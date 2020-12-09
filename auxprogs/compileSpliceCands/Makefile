CC?=gcc

compileSpliceCands : compileSpliceCands.o list.h list.o
	$(CC) $(CFLAGS)	$(CPPFLAGS) $(LDFLAGS) -o compileSpliceCands compileSpliceCands.o list.o
	cp compileSpliceCands ../../bin/compileSpliceCands

compileSpliceCands.o : compileSpliceCands.c 
	$(CC) -Wall -pedantic -ansi $(CFLAGS) $(CPPFLAGS) -c compileSpliceCands.c

all : compileSpliceCands

clean : 
	rm -f compileSpliceCands compileSpliceCands.o ../../bin/compileSpliceCands list.o

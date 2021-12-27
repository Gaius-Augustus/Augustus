# Makefile of compileSpliceCands

CC ?= gcc
CFLAGS := -Wall -pedantic -ansi ${CFLAGS}

.PHONY: all clean

all: compileSpliceCands

compileSpliceCands: compileSpliceCands.o list.h list.o
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) -o compileSpliceCands $(filter %.o,$^)
	cp -f compileSpliceCands ../../bin/compileSpliceCands

compileSpliceCands.o: compileSpliceCands.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $^

clean:
	rm -f compileSpliceCands *.o ../../bin/compileSpliceCands

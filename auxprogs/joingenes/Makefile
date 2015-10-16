CC=g++
CFLAGS=-c -Wall -std=gnu++0x

all: joingenes

joingenes: joingenes.o jg_transcript.o jg_ios.o jg_transcript.h jg_ios.h
	$(CC) joingenes.o jg_transcript.o jg_ios.o -o joingenes
	cp joingenes ../../bin/

joingenes.o: joingenes.cpp jg_transcript.h jg_ios.h
	$(CC) $(CFLAGS) joingenes.cpp

jg_transcript.o: jg_transcript.cpp jg_transcript.h
	$(CC) $(CFLAGS) jg_transcript.cpp

jg_ios.o: jg_ios.cpp jg_ios.h jg_transcript.h
	$(CC) $(CFLAGS) jg_ios.cpp

clean:
	rm -rf *o joingenes; rm -rf ../../bin/joingenes

CC=g++
CFLAGS=-c -Wall -std=gnu++0x

all: gp_call

gp_call: gp_call.o joingenes.o gp_ios.o joingenes.h gp_ios.h
	$(CC) gp_call.o joingenes.o gp_ios.o -o gp_call

gp_call.o: gp_call.cpp joingenes.h gp_ios.h
	$(CC) $(CFLAGS) gp_call.cpp

joingenes.o: joingenes.cpp joingenes.h
	$(CC) $(CFLAGS) joingenes.cpp

gp_ios.o: gp_ios.cpp gp_ios.h joingenes.h
	$(CC) $(CFLAGS) gp_ios.cpp

clean:
	rm -rf *o gp_call

# Makefile for joingenes

CXX ?= g++
CXXFLAGS := -Wall -std=gnu++0x -O2 ${CXXFLAGS}

.PHONY: all clean

all: joingenes

joingenes: joingenes.o jg_transcript.o jg_ios.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^
	mkdir -p ../../bin/
	cp -f joingenes ../../bin/

joingenes.o: joingenes.cpp jg_transcript.h jg_ios.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c joingenes.cpp

jg_transcript.o: jg_transcript.cpp jg_transcript.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c jg_transcript.cpp

jg_ios.o: jg_ios.cpp jg_ios.h jg_transcript.h
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c jg_ios.cpp

clean:
	rm -f *.o joingenes ../../bin/joingenes

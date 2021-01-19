CXX?=g++
CXXFLAGS+=-Wall -std=gnu++0x -O2

all: joingenes

joingenes: joingenes.o jg_transcript.o jg_ios.o jg_transcript.h jg_ios.h
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) joingenes.o jg_transcript.o jg_ios.o -o joingenes
	mkdir -p ../../bin/
	cp joingenes ../../bin/

joingenes.o: joingenes.cpp jg_transcript.h jg_ios.h
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) joingenes.cpp

jg_transcript.o: jg_transcript.cpp jg_transcript.h
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) jg_transcript.cpp

jg_ios.o: jg_ios.cpp jg_ios.h jg_transcript.h
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) jg_ios.cpp

clean:
	rm -rf *.o joingenes; rm -rf ../../bin/joingenes

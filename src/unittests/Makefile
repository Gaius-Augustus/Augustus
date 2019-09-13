include ../../common.mk
AINC = ../../include
CXX = g++
CXXFLAGS := -std=c++11 -Wall -pedantic -pthread ${CXXFLAGS}
LDFLAGS = -L../googletest -lgtest -lgtest_main
INCLS	+= -I$(AINC) -isystem ../googletest/include/

# Link all tests to a single test binary
unittests: lldouble_unittest.o fasta_unittest.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) $(INCLS) lldouble_unittest.o fasta_unittest.o ../lldouble.o ../fasta.o -o unittests

# Here follows each individual Unit Test
lldouble_unittest.o: lldouble_unittest.cc $(AINC)/lldouble.hh
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) $(INCLS) -c lldouble_unittest.cc

fasta_unittest.o: fasta_unittest.cc $(AINC)/fasta.hh
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) $(INCLS) -c fasta_unittest.cc

# execute all tests: run test binary
test: unittests
	./unittests

clean:
	rm -f *.o unittests

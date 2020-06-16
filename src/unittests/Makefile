include ../../common.mk
AINC = ../../include
CXX ?= g++
CXXFLAGS := -std=c++11 -Wall -pedantic -pthread ${CXXFLAGS}
LDFLAGS = -L../googletest -lgtest -lgtest_main -lgsl -lgslcblas -lmysqlclient -lmysqlpp -lsqlite3
INCLS	+= -I$(AINC) -isystem ../googletest/include/ -I/usr/include/mysql -I/usr/include/mysql++

# Link all tests to a single test binary
unittests: lldouble_unittest.o fasta_unittest.o evo_unittest.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLS) lldouble_unittest.o fasta_unittest.o evo_unittest.o \
	../lldouble.o ../fasta.o ../contTimeMC.o ../properties.o ../types.o ../geneticcode.o \
	$(LDFLAGS)  -o unittests
#evo_unittest.o \
# Here follows each individual Unit Test
lldouble_unittest.o: lldouble_unittest.cc $(AINC)/lldouble.hh
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLS) -c lldouble_unittest.cc

fasta_unittest.o: fasta_unittest.cc $(AINC)/fasta.hh
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLS) -c fasta_unittest.cc

evo_unittest.o: evo_unittest.cc $(AINC)/contTimeMC.hh
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLS) -c evo_unittest.cc

# execute all tests: run test binary
test: unittests
	./unittests

clean:
	rm -f *.o unittests

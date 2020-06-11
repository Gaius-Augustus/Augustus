include ../../common.mk
AINC = ../../include
CXX ?= g++
CXXFLAGS := -std=c++11 -Wall -pedantic -pthread ${CXXFLAGS}
LDFLAGS = -L../googletest -lgtest -lgtest_main -lgsl -lgslcblas -lmysqlclient -lmysqlpp -lsqlite3
INCLS	+= -I$(AINC) -isystem ../googletest/include/ -I/usr/include/mysql -I/usr/include/mysql++

# Link all tests to a single test binary
unittests: lldouble_unittest.o fasta_unittest.o evo_unittest.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) $(INCLS) lldouble_unittest.o fasta_unittest.o ../lldouble.o ../fasta.o ../contTimeMC.o ../genbank.o ../properties.o ../pp_profile.o ../pp_hitseq.o ../pp_scoring.o ../statemodel.o ../namgene.o \
 ../types.o ../gene.o ../evaluation.o ../motif.o ../geneticcode.o ../hints.o ../extrinsicinfo.o ../projectio.o \
 ../intronmodel.o ../exonmodel.o ../igenicmodel.o ../utrmodel.o ../merkmal.o ../vitmatrix.o ../lldouble.o ../mea.o ../graph.o \
 ../meaPath.o ../exoncand.o ../randseqaccess.o ../fasta.o ../ncmodel.o -o unittests

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

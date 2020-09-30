include ../../common.mk
AINC = ../../include
CXX ?= g++
CXXFLAGS := -std=c++11 -Wall -pedantic -pthread ${CXXFLAGS}
LDFLAGS = -L../googletest -lgtest -lgtest_main -lgsl -lgslcblas -lmysqlclient -lmysqlpp \
	-lsqlite3 -llpsolve55 -lcolamd -ldl -lsqlite3 -lboost_iostreams
INCLS	+= -I$(AINC) -isystem ../googletest/include/ -I/usr/include/mysql -I/usr/include/mysql++
OBJS    = ../lldouble.o ../fasta.o ../contTimeMC.o ../codonevo.o \
	../properties.o ../types.o \
	../geneticcode.o ../phylotree.o ../geneMSA.o ../alignment.o \
	../orthoexon.o ../exoncand.o ../intronmodel.o ../introntrain.o \
	../merkmal.o ../motif.o ../namgene.o ../gene.o ../pp_profile.o \
	../pp_scoring.o ../pp_hitseq.cc ../vitmatrix.o \
	../extrinsicinfo.o ../hints.o ../igenicmodel.o ../igenictrain.o \
	../statemodel.o ../utrmodel.o ../utrtrain.o ../exontrain.o ../exonmodel.o \
	../commontrain.o ../ncmodel.o ../evaluation.o ../projectio.o \
	../mea.o ../meaPath.o ../speciesgraph.o ../graph.o ../randseqaccess.o \
	../sqliteDB.o ../genbank.o \
	../scanner/lex.o ../parser/parse.o \

# Link all tests to a single test binary
unittests: lldouble_unittest.o fasta_unittest.o evo_unittest.o codonevo_unittest.o alignment_unittest.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLS) \
	lldouble_unittest.o fasta_unittest.o evo_unittest.o codonevo_unittest.o alignment_unittest.o \
	$(OBJS) $(LDFLAGS) -o unittests
#evo_unittest.o \
# Here follows each individual Unit Test
lldouble_unittest.o: lldouble_unittest.cc $(AINC)/lldouble.hh
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLS) -c lldouble_unittest.cc

fasta_unittest.o: fasta_unittest.cc $(AINC)/fasta.hh
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLS) -c fasta_unittest.cc

evo_unittest.o: evo_unittest.cc $(AINC)/contTimeMC.hh
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLS) -c evo_unittest.cc

codonevo_unittest.o: codonevo_unittest.cc $(AINC)/contTimeMC.hh $(AINC)/codonevo.hh 
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLS) -c codonevo_unittest.cc

alignment_unittest.o: alignment_unittest.cc $(AINC)/alignment.hh
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLS) -c alignment_unittest.cc

# execute all tests: run test binary
test: unittests
	./unittests

clean:
	rm -f *.o unittests

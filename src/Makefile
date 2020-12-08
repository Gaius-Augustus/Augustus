#
# Makefile for Augustus
#
include ../common.mk

# set TESTING to false if you don't need testing for the comparative gene prediction mode (CGP) or
# the required libraries 
# libboost-serialization-dev 
# are not available
TESTING = true

CXX?=g++

# Notes: - "-Wno-sign-compare" eliminates a high number of warnings (see footnote below). Please adopt
#          a strict signed-only usage strategy to avoid mistakes since we are not warned about this.
#        - The order of object files in $(OBJS) IS IMPORTANT (see lldouble.hh)
CXXFLAGS := -Wall -Wno-sign-compare -pedantic -O3 -std=c++11 ${CXXFLAGS}#-DDEBUG -g -ggdb -pg -DDEBUG_STATES

ifeq (,$(findstring $(ZIPINPUT),0 false False FALSE))  # if ZIPINPUT is not defined or is something else than 0, false, False or FALSE
	CPPFLAGS += -DZIPINPUT
	LIBS    = -lboost_iostreams
endif

INCLS	= -I../include

OBJS	= genbank.o properties.o pp_profile.o pp_hitseq.o pp_scoring.o statemodel.o namgene.o \
 types.o gene.o evaluation.o motif.o geneticcode.o hints.o extrinsicinfo.o projectio.o \
 intronmodel.o exonmodel.o igenicmodel.o utrmodel.o merkmal.o vitmatrix.o lldouble.o mea.o graph.o \
 meaPath.o exoncand.o randseqaccess.o fasta.o ncmodel.o
ifeq (,$(findstring $(COMPGENEPRED),0 false False FALSE))  # if COMPGENEPRED is not defined or is something else than 0, false, False or FALSE
	SQLITE   ?= true
	MYSQL    ?= true
	CPPFLAGS += -DCOMPGENEPRED
	OBJS += parser/parse.o scanner/lex.o genomicMSA.o geneMSA.o contTimeMC.o codonevo.o compgenepred.o phylotree.o orthograph.o orthoexon.o alignment.o speciesgraph.o codonMSA.o train_logReg_param.o tokenizer.o
	LIBS += -lgsl -lgslcblas # for matrix exponentiation that is required in comparative gene finding
	LIBS += -llpsolve55 -lcolamd -ldl # for mixed integer linear programming (alignment.cc)
	INCLS += -I/usr/include/lpsolve
endif

ifeq (,$(findstring $(TESTING),0 false False FALSE))  # if TESTING is not defined or is something else than 0, false, False or FALSE
	CPPFLAGS += -DTESTING
	LIBS += -lboost_serialization # for serialization/deserialization
endif

DUMOBJS = dummy.o
TOBJS	= commontrain.o igenictrain.o introntrain.o exontrain.o utrtrain.o # contentmodel.o baumwelch.o
PROGR	= augustus etraining prepareAlign fastBlockSearch pp_simScore
INFO    = cxxflags

ifneq (,$(findstring $(MYSQL),1 true True TRUE)) # if MYSQL is defined and contains 1, true, True or TRUE
	# MySQL support
	# replace RUNTIME_LIBPATH with your own run-time libpath if mysql++ is not installed system-wide, e.g.
	# RUNTIME_LIBPATH = -Wl,-rpath=~/augustus/trunks/mysql++/mysql++-3.1.0/lib
	RUNTIME_LIBPATH =
	INCLS    += -I/usr/include/mysql -I/usr/include/mysql++ # the path to mysql++ may have to be adjusted
	LIBS     += -lmysqlclient -lmysqlpp ${RUNTIME_LIBPATH} # -L/usr/local/lib/libmysqlpp.so.3.1.0
	CPPFLAGS += -DM_MYSQL
	PROGR    += load2db getSeq
endif

ifneq (,$(findstring $(SQLITE),1 true True TRUE)) # if SQLITE is defined and contains 1, true, True or TRUE
	# SQLite support
	LIBS      += -lsqlite3
	OBJS      += sqliteDB.o
	CPPFLAGS  += -DM_SQLITE
	PROGR     += load2sqlitedb getSeq
endif

all: $(OBJS) $(TOBJS) $(DUMOBJS) $(PROGR) info

.SUFFIXES:
.SUFFIXES: .cc .o .so
.PHONY: googletest

.cc.o:
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) -o $@ $< $(INCLS)

augustus: augustus.cc $(OBJS) $(DUMOBJS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(INCLS) $(LIBS)
	cp augustus ../bin/

etraining: etraining.cc $(TOBJS) $(OBJS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(INCLS) $(LIBS)
	cp etraining ../bin/

evaluate: evaluate.cc $(OBJS) $(DUMOBJS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(INCLS) $(LIBS)

fastBlockSearch: fastBlockSearch.cc pp_fastBlockSearcher.o \
	types.o properties.o geneticcode.o pp_profile.o lldouble.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(INCLS) $(LIBS)
	cp fastBlockSearch ../bin/

pp_simScore: pp_simscore.cc fasta.o types.o properties.o geneticcode.o pp_profile.o lldouble.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(INCLS) $(LIBS)
	cp pp_simScore ../bin/

prepareAlign: pp_prepare_align.cc
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^
	cp prepareAlign ../bin/

load2db: load2db.cc fasta.o hints.o types.o properties.o lldouble.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(INCLS) $(LIBS)
	cp load2db ../bin/

load2sqlitedb: load2sqlitedb.cc fasta.o hints.o types.o properties.o lldouble.o sqliteDB.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(INCLS) $(LIBS)
	cp load2sqlitedb ../bin/

getSeq: getSeq.cc $(OBJS) $(DUMOBJS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(INCLS) $(LIBS)
	cp getSeq ../bin/

espoca: espoca.cc $(OBJS) $(DUMOBJS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(INCLS) $(LIBS)
	cp espoca ../bin/

info:
	echo "$(CXXFLAGS) $(CPPFLAGS)" > $(INFO)

clean:
	rm -f ../bin/* $(PROGR) $(OBJS) $(DUMOBJS) $(TOBJS) pp_fastBlockSearcher.o $(INFO)
	cd googletest && ${MAKE} clean
	cd unittests && ${MAKE} clean

tidy: clean
	rm -f *~ *.o *.rej *.orig ../include/*~ ../include/*.orig ../include/*.rej $(INFO)
	cd googletest && ${MAKE} tidy

googletest:
	cd googletest && ${MAKE}

unittest: googletest
	cd unittests && ${MAKE}

augustus.o : \
	../include/evaluation.hh \
	../include/extrinsicinfo.hh \
	../include/genbank.hh \
	../include/gene.hh \
	../include/geneticcode.hh \
	../include/hints.hh \
	../include/lldouble.hh \
	../include/matrix.hh \
	../include/motif.hh \
	../include/namgene.hh \
	../include/pp_hitseq.hh \
	../include/pp_profile.hh \
	../include/pp_scoring.hh \
	../include/properties.hh \
	../include/statemodel.hh \
	../include/types.hh \
	../include/codonMSA.hh \
	augustus.cc

dummy.o : \
	../include/exonmodel.hh \
	../include/igenicmodel.hh \
	../include/intronmodel.hh \
	../include/lldouble.hh \
	../include/matrix.hh \
	../include/statemodel.hh \
	../include/types.hh \
	../include/utrmodel.hh \
	../include/vitmatrix.hh \
	dummy.cc

etraining.o : \
	../include/extrinsicinfo.hh \
	../include/genbank.hh \
	../include/gene.hh \
	../include/geneticcode.hh \
	../include/hints.hh \
	../include/lldouble.hh \
	../include/matrix.hh \
	../include/merkmal.hh \
	../include/motif.hh \
	../include/pp_hitseq.hh \
	../include/pp_profile.hh \
	../include/pp_scoring.hh \
	../include/properties.hh \
	../include/statemodel.hh \
	../include/types.hh \
	../include/vitmatrix.hh \
	etraining.cc

evaluate.o : \
	../include/evaluation.hh \
	../include/genbank.hh \
	../include/gene.hh \
	../include/geneticcode.hh \
	../include/lldouble.hh \
	../include/matrix.hh \
	../include/motif.hh \
	../include/pp_hitseq.hh \
	../include/pp_profile.hh \
	../include/pp_scoring.hh \
	../include/types.hh \
	../include/vitmatrix.hh \
	evaluate.cc

evaluation.o : \
	../include/evaluation.hh \
	../include/extrinsicinfo.hh \
	../include/gene.hh \
	../include/geneticcode.hh \
	../include/hints.hh \
	../include/lldouble.hh \
	../include/matrix.hh \
	../include/motif.hh \
	../include/namgene.hh \
	../include/pp_hitseq.hh \
	../include/pp_profile.hh \
	../include/pp_scoring.hh \
	../include/properties.hh \
	../include/types.hh \
	../include/vitmatrix.hh \
	evaluation.cc

exonmodel.o : \
	../include/exonmodel.hh \
	../include/extrinsicinfo.hh \
	../include/geneticcode.hh \
	../include/hints.hh \
	../include/lldouble.hh \
	../include/matrix.hh \
	../include/merkmal.hh \
	../include/motif.hh \
	../include/pp_hitseq.hh \
	../include/pp_profile.hh \
	../include/pp_scoring.hh \
	../include/projectio.hh \
	../include/properties.hh \
	../include/statemodel.hh \
	../include/types.hh \
	../include/vitmatrix.hh \
	exonmodel.cc

exontrain.o : \
	../include/commontrain.hh \
	../include/exonmodel.hh \
	../include/gene.hh \
	../include/geneticcode.hh \
	../include/lldouble.hh \
	../include/matrix.hh \
	../include/merkmal.hh \
	../include/motif.hh \
	../include/pp_hitseq.hh \
	../include/pp_profile.hh \
	../include/pp_scoring.hh \
	../include/projectio.hh \
	../include/properties.hh \
	../include/statemodel.hh \
	../include/types.hh \
	../include/vitmatrix.hh \
	exontrain.cc

extrinsicinfo.o : \
	../include/extrinsicinfo.hh \
	../include/gene.hh \
	../include/geneticcode.hh \
	../include/hints.hh \
	../include/lldouble.hh \
	../include/matrix.hh \
	../include/motif.hh \
	../include/pp_hitseq.hh \
	../include/pp_profile.hh \
	../include/pp_scoring.hh \
	../include/projectio.hh \
	../include/properties.hh \
	../include/types.hh \
	../include/vitmatrix.hh \
	extrinsicinfo.cc

genbank.o : \
	../include/genbank.hh \
	../include/fasta.hh \
	../include/gene.hh \
	../include/geneticcode.hh \
	../include/lldouble.hh \
	../include/matrix.hh \
	../include/motif.hh \
	../include/pp_hitseq.hh \
	../include/pp_profile.hh \
	../include/pp_scoring.hh \
	../include/projectio.hh \
	../include/properties.hh \
	../include/types.hh \
	../include/vitmatrix.hh \
	genbank.cc

gene.o : \
	../include/gene.hh \
	../include/geneticcode.hh \
	../include/hints.hh \
	../include/intronmodel.hh \
	../include/lldouble.hh \
	../include/matrix.hh \
	../include/motif.hh \
	../include/pp_hitseq.hh \
	../include/pp_profile.hh \
	../include/pp_scoring.hh \
	../include/statemodel.hh \
	../include/types.hh \
	../include/vitmatrix.hh \
	gene.cc

geneticcode.o : \
	../include/geneticcode.hh \
	../include/types.hh \
	geneticcode.cc

hints.o : \
	../include/hints.hh \
	../include/lldouble.hh \
	../include/types.hh \
	hints.cc

igenicmodel.o : \
	../include/extrinsicinfo.hh \
	../include/geneticcode.hh \
	../include/hints.hh \
	../include/igenicmodel.hh \
	../include/lldouble.hh \
	../include/matrix.hh \
	../include/merkmal.hh \
	../include/projectio.hh \
	../include/properties.hh \
	../include/statemodel.hh \
	../include/types.hh \
	../include/vitmatrix.hh \
	igenicmodel.cc

igenictrain.o : \
	../include/gene.hh \
	../include/geneticcode.hh \
	../include/igenicmodel.hh \
	../include/lldouble.hh \
	../include/matrix.hh \
	../include/merkmal.hh \
	../include/motif.hh \
	../include/pp_hitseq.hh \
	../include/pp_profile.hh \
	../include/pp_scoring.hh \
	../include/properties.hh \
	../include/statemodel.hh \
	../include/types.hh \
	../include/vitmatrix.hh \
	igenictrain.cc

intronmodel.o : \
	../include/extrinsicinfo.hh \
	../include/geneticcode.hh \
	../include/hints.hh \
	../include/intronmodel.hh \
	../include/lldouble.hh \
	../include/matrix.hh \
	../include/merkmal.hh \
	../include/motif.hh \
	../include/pp_hitseq.hh \
	../include/pp_profile.hh \
	../include/pp_scoring.hh \
	../include/projectio.hh \
	../include/properties.hh \
	../include/statemodel.hh \
	../include/types.hh \
	../include/vitmatrix.hh \
	intronmodel.cc

introntrain.o : \
	../include/commontrain.hh \
	../include/gene.hh \
	../include/geneticcode.hh \
	../include/intronmodel.hh \
	../include/lldouble.hh \
	../include/matrix.hh \
	../include/merkmal.hh \
	../include/motif.hh \
	../include/pp_hitseq.hh \
	../include/pp_profile.hh \
	../include/pp_scoring.hh \
	../include/projectio.hh \
	../include/properties.hh \
	../include/statemodel.hh \
	../include/types.hh \
	../include/vitmatrix.hh \
	introntrain.cc

lldouble.o : \
	../include/lldouble.hh \
	lldouble.cc

logdoubletest.o : \
	../include/lldouble.hh \
	logdoubletest.cc

motif.o : \
	../include/geneticcode.hh \
	../include/lldouble.hh \
	../include/matrix.hh \
	../include/motif.hh \
	../include/projectio.hh \
	../include/properties.hh \
	../include/types.hh \
	motif.cc

namgene.o : \
	../include/extrinsicinfo.hh \
	../include/gene.hh \
	../include/geneticcode.hh \
	../include/hints.hh \
	../include/lldouble.hh \
	../include/matrix.hh \
	../include/mea.hh \
	../include/motif.hh \
	../include/namgene.hh \
	../include/pp_hitseq.hh \
	../include/pp_profile.hh \
	../include/pp_scoring.hh \
	../include/projectio.hh \
	../include/properties.hh \
	../include/statemodel.hh \
	../include/types.hh \
	../include/vitmatrix.hh \
	../include/exoncand.hh \
	namgene.cc

projectio.o : \
	../include/projectio.hh \
	projectio.cc

properties.o : \
	../include/lldouble.hh \
	../include/properties.hh \
	../include/types.hh \
	properties.cc

statemodel.o : \
	../include/exonmodel.hh \
	../include/extrinsicinfo.hh \
	../include/geneticcode.hh \
	../include/hints.hh \
	../include/igenicmodel.hh \
	../include/intronmodel.hh \
	../include/lldouble.hh \
	../include/matrix.hh \
	../include/pp_hitseq.hh \
	../include/pp_profile.hh \
	../include/pp_scoring.hh \
	../include/properties.hh \
	../include/statemodel.hh \
	../include/types.hh \
	../include/utrmodel.hh \
	../include/vitmatrix.hh \
	statemodel.cc

types.o : \
	../include/lldouble.hh \
	../include/properties.hh \
	../include/types.hh \
	types.cc

utrmodel.o : \
	../include/extrinsicinfo.hh \
	../include/gene.hh \
	../include/geneticcode.hh \
	../include/hints.hh \
	../include/intronmodel.hh \
	../include/lldouble.hh \
	../include/matrix.hh \
	../include/merkmal.hh \
	../include/motif.hh \
	../include/pp_hitseq.hh \
	../include/pp_profile.hh \
	../include/pp_scoring.hh \
	../include/projectio.hh \
	../include/properties.hh \
	../include/statemodel.hh \
	../include/types.hh \
	../include/utrmodel.hh \
	../include/vitmatrix.hh \
	utrmodel.cc

utrtrain.o : \
	../include/commontrain.hh \
	../include/gene.hh \
	../include/geneticcode.hh \
	../include/lldouble.hh \
	../include/matrix.hh \
	../include/motif.hh \
	../include/pp_hitseq.hh \
	../include/pp_profile.hh \
	../include/pp_scoring.hh \
	../include/properties.hh \
	../include/statemodel.hh \
	../include/types.hh \
	../include/utrmodel.hh \
	../include/vitmatrix.hh \
	utrtrain.cc

ncmodel.o : \
	../include/extrinsicinfo.hh \
	../include/gene.hh \
	../include/geneticcode.hh \
	../include/hints.hh \
	../include/intronmodel.hh \
	../include/utrmodel.hh \
	../include/lldouble.hh \
	../include/matrix.hh \
	../include/merkmal.hh \
	../include/motif.hh \
	../include/projectio.hh \
	../include/properties.hh \
	../include/statemodel.hh \
	../include/types.hh \
	../include/utrmodel.hh \
	../include/vitmatrix.hh \
	ncmodel.cc

merkmal.o : \
	../include/evaluation.hh \
	../include/extrinsicinfo.hh \
	../include/gene.hh \
	../include/geneticcode.hh \
	../include/hints.hh \
	../include/lldouble.hh \
	../include/matrix.hh \
	../include/merkmal.hh \
	../include/motif.hh \
	../include/namgene.hh \
	../include/pp_hitseq.hh \
	../include/pp_profile.hh \
	../include/pp_scoring.hh \
	../include/projectio.hh \
	../include/properties.hh \
	../include/statemodel.hh \
	../include/types.hh \
	../include/vitmatrix.hh \
	merkmal.cc

pp_profile.o : \
	../include/geneticcode.hh \
	../include/lldouble.hh \
	../include/pp_profile.hh \
	../include/properties.hh \
	../include/types.hh \
	pp_profile.cc

pp_hitseq.o : \
	../include/pp_hitseq.hh \
	../include/geneticcode.hh \
	../include/lldouble.hh \
	../include/pp_profile.hh \
	../include/types.hh \
	pp_hitseq.cc

pp_scoring.o : \
	../include/pp_scoring.hh \
	../include/geneticcode.hh \
	../include/lldouble.hh \
	../include/pp_hitseq.hh \
	../include/pp_profile.hh \
	../include/properties.hh \
	../include/types.hh \
	../include/vitmatrix.hh \
	pp_scoring.cc

pp_fastBlockSearcher.o : \
	../include/pp_fastBlockSearcher.hh \
	../include/geneticcode.hh \
	../include/lldouble.hh \
	../include/pp_profile.hh \
	../include/types.hh \
	../include/vitmatrix.hh \
	pp_fastBlockSearcher.cc

vitmatrix.o : \
	../include/vitmatrix.hh \
	../include/lldouble.hh \
	../include/types.hh \
	vitmatrix.cc

commontrain.o : \
	../include/commontrain.hh \
	../include/lldouble.hh \
	../include/types.hh \
	commontrain.cc

mea.o : \
	../include/mea.hh \
	../include/gene.hh \
	../include/graph.hh \
	../include/meaPath.hh \
	mea.cc

graph.o : \
	../include/graph.hh \
	graph.cc

meaPath.o : \
	../include/meaPath.hh \

exoncand.o : \
       ../include/exoncand.hh \
       ../include/intronmodel.hh \
       ../include/types.hh \
       ../include/exonmodel.hh \
       exoncand.cc

compgenepred.o : \
	../include/compgenepred.hh \
	../include/contTimeMC.hh \
	../include/orthograph.hh \
	../include/mea.hh \
        ../include/genomicMSA.hh \
        ../include/geneMSA.hh \
	../include/orthoexon.hh \
	../include/namgene.hh \
	../include/randseqaccess.hh \
	../include/liftover.hh \
	../include/train_logReg_param.hh \
	../include/tokenizer.hh \
       compgenepred.cc

randseqaccess.o : \
       ../include/randseqaccess.hh \
       ../include/genbank.hh \
       ../include/table_structure.h \
       randseqaccess.cc

orthoexon.o : \
	../include/orthoexon.hh \
	../include/graph.hh \
	../include/exoncand.hh \
	../include/projectio.hh \
	../include/types.hh \
	../include/phylotree.hh \
	orthoexon.cc

orthograph.o : \
	../include/orthograph.hh \
	../include/contTimeMC.hh \
	../include/mea.hh \
	../include/speciesgraph.hh \
	../include/phylotree.hh \
	../include/orthogene.hh \
	../include/geneMSA.hh \
	orthograph.cc

genomicMSA.o : \
	../include/genomicMSA.hh \
	../include/exoncand.hh \
	../include/randseqaccess.hh \
	../include/alignment.hh \
	../include/geneMSA.hh \
	genomicMSA.cc

geneMSA.o : \
	../include/geneMSA.hh \
	../include/genomicMSA.hh \
	../include/exoncand.hh \
	../include/orthoexon.hh \
	../include/intronmodel.hh \
	../include/namgene.hh \
	../include/randseqaccess.hh \
	../include/alignment.hh \
	../include/orthograph.hh \
	../include/geneticcode.hh \
	geneMSA.cc

parser/parse.o: \
	parser/parser.ih \
	parser/parse.cc

scanner/lex.o: \
	scanner/scanner.ih \
	scanner/lex.cc

phylotree.o: \
	../include/phylotree.hh \
	parser/parser.h \
	../include/properties.hh \
	../include/contTimeMC.hh \
	phylotree.cc

speciesgraph.o: \
	../include/speciesgraph.hh \
	../include/graph.hh \
	speciesgraph.cc

fasta.o: \
	../include/fasta.hh \
	fasta.cc

contTimeMC.o: \
	../include/contTimeMC.hh \
	contTimeMC.cc

codonevo.o: \
	../include/codonevo.hh \
	../include/contTimeMC.hh \
	codonevo.cc

alignment.o: \
	../include/alignment.hh \
	alignment.cc

sqliteDB.o: \
	../include/sqliteDB.hh \
	sqliteDB.cc

codonMSA.o: \
	../include/codonMSA.hh \
	codonMSA.cc

train_logReg_param.o: \
	../include/train_logReg_param.hh \
	train_logReg_param.cc

tokenizer.o: \
	../include/tokenizer.hh \
	tokenizer.cc

#    1) signed-unsigned comparisons are needed because the methods
#       std::vector<T>::size or std::string::length unfortunately
#       return unsigned types. Using a negative integer in a
#       comparison i<size() would yield false by default (which in
#       most cases is actually a desired behaviour since negative
#       indices are invalid). The same behaviour would be achieved by
#       (unsigned)i<size().  However, to remove the need of writing
#       this in every single loop header, we eliminate the warning
#       (which does not warn in every dangerous case anyway:
#       e.g. consider an expression like "u < size()-1" which compiles
#       without warning if u is unsigned but will still be true if
#       size() equals 0).

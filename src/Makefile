#
# Makefile for Augustus
#
include ../common.mk
CC	= g++

# Notes: - "-Wno-sign-compare" eliminates a high number of warnings (see footnote below). Please adopt 
#          a strict signed-only usage strategy to avoid mistakes since we are not warned about this.
#        - The order of object files in $(OBJS) IS IMPORTANT (see lldouble.hh)
CFLAGS := -Wall -Wno-sign-compare -ansi -pedantic -g -ggdb -O2 ${CFLAGS} # -DDEBUG -g -ggdb -pg

ifdef ZIPINPUT
	CFLAGS += -DZIPINPUT
	LIBS    = -lboost_iostreams
endif

ifdef COMPGENEPRED
	CFLAGS += -std=c++0x -DCOMPGENEPRED
endif
INCLS	= -I../include

OBJS	= genbank.o properties.o pp_profile.o pp_hitseq.o pp_scoring.o statemodel.o namgene.o \
 types.o gene.o evaluation.o motif.o geneticcode.o hints.o extrinsicinfo.o projectio.o \
 intronmodel.o exonmodel.o igenicmodel.o utrmodel.o merkmal.o vitmatrix.o lldouble.o mea.o graph.o \
 meaPath.o exoncand.o randseqaccess.o speciesgraph.o\
  fasta.o
ifdef COMPGENEPRED
	OBJS += parser/parse.o scanner/lex.o genomicMSA.o geneMSA.o contTimeMC.o compgenepred.o phylotree.o orthograph.o orthoexon.o alignment.o
	LIBS += -lgsl -lgslcblas # for matrix exponentiation that is required in comparative gene finding
endif
DUMOBJS = dummy.o
TOBJS	= commontrain.o igenictrain.o introntrain.o exontrain.o utrtrain.o # contentmodel.o baumwelch.o 
PROGR	= augustus etraining consensusFinder curve2hints prepareAlign fastBlockSearch
INFO    = cflags

ifdef MYSQL
# replace RUNTIME_LIBPATH with your own run-time libpath if mysql++ is not installed system-wide, e.g.
# RUNTIME_LIBPATH = -Wl,-rpath=/home/mario/augustus/trunks/mysql++/mysql++-3.1.0/lib 
RUNTIME_LIBPATH = 
INCLS	+= -I/usr/include/mysql -I/usr/include/mysql++ # the path to mysql++ may have to be adjusted
LIBS 	+=  -lmysqlclient -lmysqlpp ${RUNTIME_LIBPATH} # -L/usr/local/lib/libmysqlpp.so.3.1.0
CFLAGS  += -DAMYSQL
PROGR   += load2db getSeq
endif

ifdef SQLITE
LIBS 	+= -lsqlite3
OBJS	+= sqliteDB.o
CFLAGS	+= -DSQLITE
PROGR   += load2sqlitedb getSeq
endif

all: $(OBJS) $(TOBJS) $(DUMOBJS) $(PROGR) info

.SUFFIXES:
.SUFFIXES: .cc .o .so

.cc.o:
	$(CC) -c $(CFLAGS) -o $@ $< $(INCLS)

augustus: augustus.cc $(OBJS) $(DUMOBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(INCLS) $(LIBS)
	cp augustus ../bin/

etraining: etraining.cc $(TOBJS) $(OBJS) 
	$(CC) $(CFLAGS) -o $@ $^ $(INCLS) $(LIBS)
	cp etraining ../bin/

evaluate: evaluate.cc $(OBJS) $(DUMOBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(INCLS) $(LIBS)

consensusFinder: consensusFinder.cc consensus.o 
	$(CC) $(CFLAGS) -o $@ $^ $(INCLS) $(LIBS)
	cp consensusFinder ../bin/

curve2hints: curve2hints.cc exon_seg.o 
	$(CC) $(CFLAGS) -o $@ $^ $(INCLS) $(LIBS)
	cp curve2hints ../bin/

fastBlockSearch: fastBlockSearch.cc pp_fastBlockSearcher.o \
	types.o properties.o geneticcode.o pp_profile.o lldouble.o
	$(CC) $(CFLAGS) -o $@ $^ $(INCLS) $(LIBS)
	cp fastBlockSearch ../bin/

prepareAlign: pp_prepare_align.cc
	$(CC) $(CFLAGS) -o $@ $^
	cp prepareAlign ../bin/

load2db: load2db.cc fasta.o hints.o types.o properties.o lldouble.o
	$(CC) $(CFLAGS) -o $@ $^ $(INCLS) $(LIBS)
	cp load2db ../bin/

load2sqlitedb: load2sqlitedb.cc fasta.o hints.o types.o properties.o lldouble.o sqliteDB.o
	$(CC) $(CFLAGS) -o $@ $^ $(INCLS) $(LIBS)
	cp load2sqlitedb ../bin/

getSeq: getSeq.cc $(OBJS) $(DUMOBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(INCLS) $(LIBS)
	cp getSeq ../bin/

info:
	echo "$(CFLAGS)" > $(INFO)

clean:
	rm -f $(PROGR) $(OBJS) $(DUMOBJS) $(TOBJS) consensus.o exon_seg.o pp_fastBlockSearcher.o  $(INFO) 

tidy: clean
	rm -f *~ *.o *.rej *.orig ../include/*~ ../include/*.orig ../include/*.rej $(INFO)

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

consensus.o : \
	../include/consensus.hh \
	../include/geneticcode.hh \
	../include/types.hh \
	consensus.cc

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
	../include/orthograph.hh \
	../include/mea.hh \
        ../include/genomicMSA.hh \
        ../include/geneMSA.hh \
	../include/orthoexon.hh \
	../include/namgene.hh \
	../include/randseqaccess.hh \
       compgenepred.cc

randseqaccess.o : \
       ../include/randseqaccess.hh \
       ../include/genbank.hh \
       ../include/table_structure.h \
       randseqaccess.cc

orthoexon.o : \
	../include/orthoexon.hh \
	../include/exoncand.hh \
	../include/projectio.hh \
	../include/orthograph.hh \
	../include/types.hh \
	orthoexon.cc

orthograph.o : \
	../include/orthograph.hh \
	../include/mea.hh \
	orthograph.cc

genomicMSA.o : \
	../include/genomicMSA.hh \
	../include/exoncand.hh \
	../include/orthoexon.hh \
	../include/randseqaccess.hh \
	../include/alignment.hh \
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
	../include/properties.hh\
	phylotree.cc

speciesgraph.o: \
	../include/speciesgraph.hh \
	speciesgraph.cc

fasta.o: \
	../include/fasta.hh \
	fasta.cc

contTimeMC.o: \
	../include/contTimeMC.hh \
	contTimeMC.cc

alignment.o: \
	../include/alignment.hh \
	alignment.cc

sqliteDB.o: \
	../include/sqliteDB.hh \
	sqliteDB.cc

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

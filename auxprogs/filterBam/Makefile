#	Makefile corresponding to computeMatePairs
#	
# 	NOTE:
#	a) "-g" option allows for debugging
#	b) libbamtools.a and lz library are used during the linking process
#
# 	Created: 11-November-2011
#	Last modified: 24-February-2012


# test: ttview checkBamEnv

# checkBamEnv: 
# 	@if [ -z "${BAMTOOLS}" ]; then \
# 		echo '[Makefile]: $${BAMTOOLS} has not been set.'; \
# 		exit 101; \
# 	elif [ "${BAMTOOLS}" ]; then \
# 		echo '[Makefile]: variable $${BAMTOOLS} is defined.'; \
# 		echo "Compiling with BAMTOOLS=${BAMTOOLS}"; \
# 	fi

# ttview: checkBamEnv
# 	@echo $(BAMTOOLS)/filterBam.cc 


# # Variable definition
# all:build pdf checkBamEnv

# checkBamEnv: 
# 	@if [ -z "${BAMTOOLS}" ]; then \
# 		echo '[Makefile]: $${BAMTOOLS} has not been set.'; \
# 		exit 101; \
# 	elif [ "${BAMTOOLS}" ]; then \
# 		echo '[Makefile]: variable $${BAMTOOLS} is defined.'; \
# 		echo "Compiling with BAMTOOLS=${BAMTOOLS}"; \
# 	fi

BAMTOOLS = /home/tonatiuh/web-home/bamtools
INCLUDES = -I$(BAMTOOLS)/include -Iheaders -I$(BAMTOOLS)/src/toolkit
VPATH = functions # Required to find files in functions/foo.cc
LIBS = $(BAMTOOLS)/lib/libbamtools.a -lz 
SOURCES = filterBam.cc MatePairs.cc getReferenceName.cc initOptions.cc SingleAlignment.cc \
			printElapsedTime.cc sumMandIOperations.cc sumDandIOperations.cc PairednessCoverage.cc
PROGRAM = filterBam
OBJECTS = $(SOURCES:.cc=.o)
CFLAGS = -g $(INCLUDES) -std=c++0x -pg 
#LINK.cc = g++


######################################################
# Compilation of filterBam
######################################################
# Recipes
# $@: full target name of current target. 
# $<: .cc file of target. 
build : $(PROGRAM) $(OBJECTS)
$(PROGRAM) : $(OBJECTS) 
	$(LINK.cc) $(CFLAGS) $^ -o $@ $(LIBS) 

# This recipe is equivalent to specifying:
# $^: is a list of all the prerequisites of the rule
# % contains the stem-name of the file
$(OBJECTS) : %.o:%.cc
	$(LINK.cc) $(CFLAGS) -c $^ -o $@ 


######################################################
# Generation of PDF document
######################################################

pdf : $(PROGRAM).tex
	pdflatex $^
	pdflatex $^ 


######################################################
# Generation of Dummy datasets
######################################################
SAMFILES = example_1.sam
SCRIPT = sam2bam.sh

bam : $(SAMFILES)
	sh $(SCRIPT) $^

######################################################
# Cleaning routine 
######################################################

PDFOBJECTS = $(PROGRAM).aux $(PROGRAM).blg $(PROGRAM).bbl $(PROGRAM).log $(PROGRAM).out

clean: 
	rm -f *~ $(OBJECTS) $(PDFOBJECTS)

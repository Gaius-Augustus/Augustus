#	Makefile corresponding to computeMatePairs
#	
# 	NOTE:
#	a) "-g" option allows for debugging
#	b) libbamtools.a and lz library are used during the linking process
#
# 	Created: 11-November-2011
#	Last modified: 2-February-2012

# Variable definition
BAMTOOLS_ROOT := $(shell find ~ -name bamtools -type d) # Change "~" by "/" ?
BAMTOOLS_ROOT := $(shell echo $(BAMTOOLS_ROOT)| sed -e "s/\([a-zA-Z,\/]*\) .*/\1/" )
#BAMTOOLS_ROOT = /home/tonatiuh/web-home/bamtools
INCLUDES = -I$(BAMTOOLS_ROOT)/include -Iheaders -I$(BAMTOOLS_ROOT)/src/toolkit
VPATH = functions # Required to find files in functions/foo.cc
LIBS = $(BAMTOOLS_ROOT)/lib/libbamtools.a -lz 
SOURCES = filterBam.cc MatePairs.cc getReferenceName.cc initOptions.cc SingleRead.cc \
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
$(PROGRAM) : $(OBJECTS) 
	$(LINK.cc) $(CFLAGS) $^ -o $@ $(LIBS) 

# This recipe is equivalent to specifying:
# $^: is a list of all the prerequisites of the rule
# % contains the stem-name of the file
all : $(OBJECTS) # Compiles each foo.cc into foo.o
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

PDFOBJECTS = $(PROGRAM).aux $(PROGRAM).blg $(PROGRAM).bbl $(PROGRAM).log

clean: 
	rm -f *~ $(OBJECTS) $(PDFOBJECTS)

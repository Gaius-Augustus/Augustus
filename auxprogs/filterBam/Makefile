#	Makefile corresponding to computeMatePairs
#	
# 	NOTE:
#	a) "-g" option allows for debugging
#	b) libbamtools.a and lz library are used during the linking process
#
# 	Created: 11-November-2011
#	Last modified: 2-February-2012

# Variable definition
HOME = /home/tonatiuh/web-home/
#BT1 = $(HOME)/pezmaster31-bamtools-12aca98
#BT2 = $(HOME)/bamtools-1.0.2
BT3 = $(HOME)/bamtools
CODE = $(HOME)/filterBam/bamtoolsCode
INCLUDES = -I$(BT3)/include -Iheaders -I$(BT3)/src/toolkit
VPATH = functions # Required to find files in functions/foo.cc
LIBS = $(BT3)/lib/libbamtools.a -lz 
SOURCES = filterBam.cc MatePairs.cc getReferenceName.cc initOptions.cc SingleRead.cc \
			printElapsedTime.cc sumMandIOperations.cc sumDandIOperations.cc PairednessCoverage.cc
PROGRAM = filterBam
OBJECTS = $(SOURCES:.cc=.o)
CFLAGS = -g $(INCLUDES) -std=c++0x 
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
build : $(OBJECTS) # Compiles each foo.cc into foo.o
$(OBJECTS) : %.o:%.cc
	$(LINK.cc) $(CFLAGS) -c $^ -o $@ 


######################################################
# Generation of PDF document
######################################################

PDFOBJECTS = $(PROGRAM).aux $(PROGRAM).blg $(PROGRAM).bbl $(PROGRAM).log

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

clean: 
	rm -f *~ $(OBJECTS) $(PDFOBJECTS)

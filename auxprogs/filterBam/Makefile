#	Makefile corresponding to computeMatePairs
#	
# 	NOTE:
#	a) "-g" option allows for debugging
#	b) libbamtools.a and lz library are used during the linking process
#
# 	Created: 11-November-2011
#	Last modified: 25-January-2012

# Variable definition
HOME = /home/tonatiuh/
BT1 = $(HOME)/pezmaster31-bamtools-12aca98
BT2 = $(HOME)/bamtools-1.0.2
BT3 = $(HOME)/bamtools
CODE = $(HOME)/filterBam/bamtoolsCode
INCLUDES = -I$(BT3)/include -Iheaders -I$(BT3)/src/toolkit
VPATH = functions # Required to find files in functions/foo.cc
LIBS = $(BT3)/lib/libbamtools.a -lz -lboost_regex
SOURCES = filterBam.cc MatePairs.cc getReferenceName.cc initOptions.cc \
			printElapsedTime.cc sumMandIOperations.cc sumDandIOperations.cc
PROGRAM = filterBam
OBJECTS = $(SOURCES:.cc=.o)
CFLAGS = -g $(INCLUDES) -std=c++0x 
#LINK.cc = g++


# Recipe(s)
# $@: full target name of current target. 
# $<: .cc file of target. 
$(PROGRAM) : $(OBJECTS) 
	$(LINK.cc) $(CFLAGS) $^ -o $@ $(LIBS) 

# This recipe is equivalent to specifying:
# $^: is a list of all the prerequisites of the rule
# % contains the stem-name of the file
all:$(OBJECTS) # Compiles each foo.cc into foo.o
$(OBJECTS) : %.o:%.cc
	$(LINK.cc) $(CFLAGS) -c $^ -o $@ 


clean: 
	rm -f *~ $(OBJECTS)

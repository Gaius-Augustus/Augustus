#	Makefile corresponding to printCoverage
#	
# 	NOTE:
#	a) "-g" option allows for debugging
#	b) libbamtools.a and lz library are used during the linking process
#
# 	Created: 25-Aug-11
#	Last modified: 4-July-2012

# Variable definition
HOME = /home/tonatiuh/Documents/git
BT3 = $(HOME)/bamtools
CODE = ../
INCLUDES = -I$(BT3)/include -I$(CODE)/headers -I$(BT3)/src/toolkit
VPATH = $(CODE)/functions # Required to find files in functions/foo.cc
LIBS = $(BT3)/lib/libbamtools.a -lz 
SOURCES = chunkSortedBam.cc getReferenceName.cc SingleAlignment.cc printElapsedTime.cc initOptions.cc
PROGRAM = chunkSortedBam
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
all : $(OBJECTS) # Compiles each foo.cc into foo.o
$(OBJECTS) : %.o:%.cc
	$(LINK.cc) $(CFLAGS) -c $^ -o $@ 


clean: 
	rm -f *~ $(OBJECTS)

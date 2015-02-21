#	Makefile corresponding to BAM2HINTS
#	
# 	NOTE:
#	a) "-g" option allows for debugging
#	b) libbamtools.a and lz library are used during the linking process
#
# 	Created: 14-December-2011
#	Last modified:  4-May-2012

# Variable definition
# on greifweb1
BAMTOOLS = /home/mario/tools/bamtools
INCLUDES = $(BAMTOOLS)/include
LIBS = $(BAMTOOLS)/lib/libbamtools.a -lz
SOURCES = bam2hints.cc 
OBJECTS = $(SOURCES:.cc=.o)
CFLAGS = -Wall -O0 # -g -p

LINK.cc = g++

# Recipe(s)
# $@: full target name of current target. 
# $<: .c file of target. 
bam2hints : $(OBJECTS) 
	$(LINK.cc) $(CFLAGS) -o $@ $(OBJECTS) $(LIBS) 

all:$(OBJECTS) # Compiles each foo.cc into foo.o

bam2hints.o : $(SOURCES)
	$(LINK.cc) $(CFLAGS) -c $< -o $@ -I$(INCLUDES) 

clean: 
	rm -f $(OBJECTS)

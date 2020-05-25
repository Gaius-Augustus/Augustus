#	Makefile corresponding to BAM2HINTS
#	
# 	NOTE:
#	a) "-g" option allows for debugging
#	b) libbamtools.a and lz library are used during the linking process
#
# 	Created: 14-December-2011
#	Last modified:  09-October-2015 by Katharina J. Hoff

# Variable definition
INCLUDES = -I/usr/include/bamtools
LIBS = -lbamtools -lz
SOURCES = bam2hints.cc 
OBJECTS = $(SOURCES:.cc=.o)
CXXFLAGS += -Wall -O2 # -g -p -g -ggdb 

LINK.cc ?= g++

# Recipe(s)
# $@: full target name of current target. 
# $<: .c file of target. 
bam2hints : $(OBJECTS) 
	$(LINK.cc) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) -o $@ $(OBJECTS) $(LIBS)
	mkdir -p ../../bin
	cp bam2hints ../../bin

all:$(OBJECTS) # Compiles each foo.cc into foo.o

bam2hints.o : $(SOURCES)
	$(LINK.cc) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@ $(INCLUDES) 

clean:
	rm -f bam2hints.o bam2hints ../../bin/bam2hints

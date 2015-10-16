#	Makefile corresponding to BAM2HINTS
#	
# 	NOTE:
#	a) "-g" option allows for debugging
#	b) libbamtools.a and lz library are used during the linking process
#
# 	Created: 14-December-2011
#	Last modified:  09-October-2015 by Katharina J. Hoff

# Variable definition
INCLUDES = /usr/include/bamtools
LIBS = -lbamtools -lz
SOURCES = bam2hints.cc 
OBJECTS = $(SOURCES:.cc=.o)
CFLAGS = -Wall -O0 -g -ggdb # -g -p

LINK.cc = g++

# Recipe(s)
# $@: full target name of current target. 
# $<: .c file of target. 
bam2hints : $(OBJECTS) 
	$(LINK.cc) $(CFLAGS) -o $@ $(OBJECTS) $(LIBS) 
	cp bam2hints ../../bin

all:$(OBJECTS) # Compiles each foo.cc into foo.o

bam2hints.o : $(SOURCES)
	$(LINK.cc) $(CFLAGS) -c $< -o $@ -I$(INCLUDES) 

clean: 
	rm -f bam2hints.o bam2hints ../../bin/bam2hints

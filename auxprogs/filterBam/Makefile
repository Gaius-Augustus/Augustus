#	Makefile corresponding to FILTERBAM
#	
# 	NOTE:
#	a) "-g" option allows for debugging
#	b) libbamtools.a and lz library are used during the linking process
#
# 	Created: 28-September-2011
#	Last modified: 14-October-2011

# Variable definition
HOME = .
BAMTOOLS = /home/tonatiuh/web-home/pezmaster31-bamtools-12aca98/
INCLUDES = $(BAMTOOLS)/include
LIBS = $(BAMTOOLS)/lib/libbamtools.a -lz
SOURCES = filterBam.cc 
OBJECTS = $(SOURCES:.cc=.o)
CFLAGS = -g -std=c++0x # -std=c++0x is required for support of unordered maps
#LINK.cc = g++

# Recipe(s)
# $@: full target name of current target. 
# $<: .c file of target. 
filterBam : $(OBJECTS) 
	$(LINK.cc) $(CFLAGS) -o $@ $(OBJECTS) $(LIBS) 

filterBam.o : $(SOURCES)
	$(LINK.cc) $(CFLAGS) -c $< -o $@ -I$(INCLUDES) 

clean: 
	rm -f *~ $(OBJECTS)

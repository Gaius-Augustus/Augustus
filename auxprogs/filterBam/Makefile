#	Makefile corresponding to FILTERBAM
#	
# 	NOTE:
#	a) "-g" option allows for debugging
#	b) libbamtools.a and lz library are used during the linking process
#
# 	Created: 28-September-2011
#	Last modified: 15-November-2011

# Variable definition
HOME = /home/tonatiuh/web-home
BT1 = $(HOME)/pezmaster31-bamtools-12aca98
BT2 = $(HOME)/bamtools-1.0.2
CODE = $(HOME)/filterBam/bamtoolsCode
INCLUDES = -I$(BT2)/include -I$(CODE)/headers
VPATH = $(CODE)/functions # Required to find files in functions/foo.cc
LIBS = $(BT2)/lib/libbamtools.a -lz 
SOURCES = filterBam.cc 
PROGRAM = filterBam
OBJECTS = $(SOURCES:.cc=.o)
CFLAGS = -g $(INCLUDES) -std=c++0x # -std=c++0x is required for support of unordered maps
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
	$(LINK.cc) $(CFLAGS) -c $^ -o $@ $(LIBS) 

# Create accompanying pdf file
pdf:filterBam.tex
	pdflatex filterBam.tex

clean: 
	rm -f *~ $(OBJECTS) *.log *.aux

#	Makefile for filterBam
#	
# 	Created: 11-November-2011
#	Last modified: 4-March-2012

SOURCES = filterBam.cc MatePairs.cc getReferenceName.cc initOptions.cc SingleAlignment.cc \
			printElapsedTime.cc sumMandIOperations.cc sumDandIOperations.cc PairednessCoverage.cc
PROGRAM = filterBam
OBJECTS = $(SOURCES:.cc=.o)
CFLAGS = -g $(INCLUDES) -std=c++0x -pg 
#LINK.cc = g++
# #BAMTOOLS = /home/tonatiuh/web-home/bamtools
INCLUDES = -I$(BAMTOOLS)/include -Iheaders -I$(BAMTOOLS)/src/toolkit
VPATH = functions # Required to find files in functions/foo.cc
LIBS = $(BAMTOOLS)/lib/libbamtools.a -lz 

# Variable definition
all : $(PROGRAM) $(OBJECTS) CHECKBAM

# Define destination directory
CHECKBAM: 
	@if [ -z "${BAMTOOLS}" ]; then \
		echo '[Makefile]: $${BAMTOOLS} has not been set.'; \
		exit 101; \
	elif [ "${BAMTOOLS}" ]; then \
		echo "filterBam compiled with BAMTOOLS=${BAMTOOLS}"; \
	fi

$(PROGRAM) : $(OBJECTS) 
	$(LINK.cc) $(CFLAGS) $^ -o $@ $(LIBS) 

# This recipe is equivalent to specifying:
# $^: is a list of all the prerequisites of the rule
# % contains the stem-name of the file
$(OBJECTS) : %.o:%.cc
	$(LINK.cc) $(CFLAGS) -c $^ -o $@ 


# ######################################################
# # Generation of Dummy datasets
# ######################################################
# SAMFILES = example_1.sam
# SCRIPT = sam2bam.sh

# bam : $(SAMFILES)
# 	sh $(SCRIPT) $^

######################################################
# Cleaning routine 
######################################################


clean: 
	rm -f *~ $(OBJECTS) 

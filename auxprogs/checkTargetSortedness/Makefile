# Tests whether BAM file has been sorted by target reference (RefID)
# 
# NOTE: 
# Depending on the version of the compiler,  the call to "-lcurses" might have 
# to be replaced to "-lncurses"
#
#
# Created: 11-April-2012
# Last modified:  9-July-2012
# 

SOURCES = checkTargetSortedness.c bam_index.o ### Keep in mind that "bam_index.o" needs to be included!!!
PROGRAM = checkTargetSortedness
OBJECTS = $(SOURCES:.c=.o) 
CFLAGS=-g -Wall -O2
# Change variable below so that it points to your SAMTOOLS folder...
SAMTOOLS=/home/tonatiuh/Documents/git/samtools
INCLUDES=-I$(SAMTOOLS) -I.
VPATH=$(SAMTOOLS)
LIBS=$(SAMTOOLS)/libbam.a -lcurses -lm -lz -lpthread
CC=gcc

all:$(PROGRAM) $(OBJECTS)

$(PROGRAM):$(OBJECTS)
	$(CC) $(CFLAGS)	$^ -o $@ $(LIBS) 

$(OBJECTS):%.o:%.c
	$(CC) $(CFLAGS) -c $^ -o $@ $(INCLUDES)

clean: 
	rm *.o checkTargetSortedness

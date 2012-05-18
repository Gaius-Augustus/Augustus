# Compiles an example program that makes use of the SAMtools API
# 
# Created: 11-April-2012
# Last modified: 18-May-2012
# 

SOURCES = checkTargetSortedness.c bam_index.o ### Keep in mind that "bam_index.o" needs to be included!!!
PROGRAM = checkTargetSortedness
OBJECTS = $(SOURCES:.c=.o) 
CFLAGS=-g -Wall -O2
SAMTOOLS=/home/tonatiuh/Documents/git/samtools
INCLUDES=-I$(SAMTOOLS) -I.
VPATH=$(SAMTOOLS)
LIBS=$(SAMTOOLS)/libbam.a -lcurses -lm -lz 
CC=gcc

all:$(PROGRAM) $(OBJECTS)

$(PROGRAM):$(OBJECTS)
	$(CC) $(CFLAGS)	$^ -o $@ $(LIBS) 

$(OBJECTS):%.o:%.c
	$(CC) $(CFLAGS) -c $^ -o $@ $(INCLUDES)

clean: 
	rm *.o 
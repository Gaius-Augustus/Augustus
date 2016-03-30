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

SOURCES = checkTargetSortedness.c
PROGRAM = checkTargetSortedness
OBJECTS = $(SOURCES:.c=.o) 
CFLAGS := -g -Wall -O2 ${CFLAGS}
# Change variable below so that it points to your SAMTOOLS folder...
SAMTOOLS = /usr/include/samtools
INCLUDES = -I$(SAMTOOLS) -I.
# replace -lbam with the following in case you have your own samtools library: $(SAMTOOLS)/libbam.a
LIBS=-lbam -lcurses -lm -lz -lpthread
CC=gcc

all:$(PROGRAM) $(OBJECTS)

$(PROGRAM):$(OBJECTS)
	$(CC) $(CFLAGS)	$(CPPFLAGS) $(LDFLAGS) $^ -o $@ $(LIBS) 

$(OBJECTS):%.o:%.c
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) -c $^ -o $@ $(INCLUDES)

clean:
	rm -f *.o checkTargetSortedness

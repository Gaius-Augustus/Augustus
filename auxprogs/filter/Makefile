# Export the variable BAMTOOLS with value pointing 
# to bamtools
#	
# Created:  9-December-2012
# Last modified: 11-April-2013

# Variable definition 
BAMTOOLS=/home/tonatiuh/Documents/git/bamtools
JSON=/home/tonatiuh/Documents/git/Bamtools # change accordingly
INCLUDES = -I$(BAMTOOLS)/include -I$(BAMTOOLS)/src/toolkit -I$(BAMTOOLS)/src -I$(BAMTOOLS)/src/utils -I$(BAMTOOLS)/src/third_party -I.
CFLAGS = -O3 -DNDEBUG -g $(INCLUDES) -L$(BAMTOOLS)/lib
LIBS = -lz -lbamtools -Wl,-rpath,$(BAMTOOLS)/lib $(BAMTOOLS)/lib/libjsoncpp.so.1.0.0 $(BAMTOOLS)/lib/libbamtools-utils.so.2.2.0
VPATH=$(BAMTOOLS)/src/toolkit
#OBJECTS=bamtools_convert.o bamtools_count.o bamtools_coverage.o bamtools_filter.o bamtools_header.o bamtools_index.o bamtools_merge.o bamtools_random.o bamtools_revert.o bamtools_sort.o bamtools_split.o bamtools_stats.o bamtools_illumina.o filter.o
OBJECTS=bamtools_illumina.o filter.o

all : $(OBJECTS)
	$(LINK.cc) -o filter $(OBJECTS) $(LIBS) $(CFLAGS)

# bamtools_stats.o : bamtools_stats.cpp
# 	$(LINK.cc) $(CFLAGS) -c $(BAMTOOLS)/src/toolkit/bamtools_stats.cpp -o bamtools_stats.o

bamtools_illumina.o : bamtools_illumina.cpp
	$(LINK.cc) $(CFLAGS) -c bamtools_illumina.cpp -o bamtools_illumina.o

filter.o : filter.cpp
	$(LINK.cc) $(CFLAGS) -c filter.cpp -o filter.o

clean: 
	rm -f filter *.o


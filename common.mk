# Definitions common to all Makefiles
# This file is included from other Makefiles in the augustus project.
AUGVERSION = 3.3.2

# make this a comment if you get compilation errors about the boost iostreams library
# cannot install this library and do not need to input gzipped input genome files
#ZIPINPUT = true

# uncomment this line to enable comparative gene finding (requires compiler which supports C++11 standard)
COMPGENEPRED = true

# uncomment this line when you need MySQL access to sequences (most users don't)
# MYSQL = true

# uncomment this line to enable access to SQLite databases that store
# file offsets of sequence data in flat files and hints
# SQLITE = true

# base path for all other tools

TOOLDIR = $(HOME)/tools/

BAMTOOL=$(TOOLDIR)/bamtools
BAMTOOL_INC = $(BAMTOOL)/include
BAMTOOL_LIBS = -lbamtools

SAMTOOLS = $(TOOLDIR)/samtools
SAMTOOLS_LIBS = -lbam -L$(SAMTOOLS)
SAMTOOLS_INC = $(SAMTOOLS)/include

HTSLIB = $(TOOLDIR)/htslib
HTSLIB_LIBS=-lhts -L$(HTSLIB)/libs
HTSLIB_INC=$(HTSLIB)/include

LPSOLVE_INC = /usr/include/lpsolve

BCFTOOLS_INC = $(TOOLDIR)/bcftools

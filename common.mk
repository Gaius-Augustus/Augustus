# Definitions common to all Makefiles
# This file is included from other Makefiles in the augustus project.
AUGVERSION = 3.4.0

# set ZIPINPUT to false if you do not require input gzipped input genome files,
# get compilation errors about the boost iostreams library or
# the required libraries libboost-iostreams-dev and lib1g-dev are not available
ZIPINPUT = true

# set COMPGENEPRED to false if you do not require the comparative gene prediction mode (CGP) or
# the required libraries
# libgsl-dev, libboost-all-dev, libsuitesparse-dev, liblpsolve55-dev, libmysql++-dev and libsqlite3-dev
# are not available
COMPGENEPRED = true

# set EBONY to true to enable ClaMSA inference server queries for the ebony model (splice site scores) (only for CGP mode)
EBONY = false

# set these paths to the correct locations if you have installed the corresponding packages in non-default locations:
#INCLUDE_PATH_ZLIB        := -I/usr/include
#LIBRARY_PATH_ZLIB        := -L/usr/lib/x86_64-linux-gnu -Wl,-rpath,/usr/lib/x86_64-linux-gnu
#INCLUDE_PATH_BOOST       := -I/usr/include
#LIBRARY_PATH_BOOST       := -L/usr/lib/x86_64-linux-gnu -Wl,-rpath,/usr/lib/x86_64-linux-gnu
#INCLUDE_PATH_LPSOLVE     := -I/usr/include/lpsolve
#LIBRARY_PATH_LPSOLVE     := -L/usr/lib/lp_solve/ -Wl,-rpath,/usr/l ib/lp_solve/
#INCLUDE_PATH_SUITESPARSE := -I/usr/include/suitesparse/
#LIBRARY_PATH_SUITESPARSE := -L/usr/lib/x86_64-linux-gnu -Wl,-rpath,/usr/lib/x86_64-linux-gnu
#INCLUDE_PATH_GSL         := -I/usr/include
#LIBRARY_PATH_GSL         := -L/usr/lib/x86_64-linux-gnu -Wl,-rpath,/usr/lib/x86_64-linux-gnu
#INCLUDE_PATH_MYSQL       := -I/usr/include -I/usr/include/mysql      # the path to mysql++ may have to be adjusted
#LIBRARY_PATH_MYSQL       := -L/usr/lib -L/usr/lib/x86_64-linux-gnu -Wl,-rpath,/usr/lib -Wl,-rpath,/usr/lib/x86_64-linux-gnu
#INCLUDE_PATH_SQLITE      := -I/usr/include
#LIBRARY_PATH_SQLITE      := -L/usr/lib/x86_64-linux-gnu -Wl,-rpath,/usr/lib/x86_64-linux-gnu
#INCLUDE_PATH_BAMTOOLS    := -I/usr/include/bamtools
#LIBRARY_PATH_BAMTOOLS    := -L/usr/lib/x86_64-linux-gnu -Wl,-rpath,/usr/lib/x86_64-linux-gnu
#INCLUDE_PATH_HTSLIB      := -I/usr/include/htslib
#LIBRARY_PATH_HTSLIB      := -L/usr/lib/x86_64-linux-gnu -Wl,-rpath,/usr/lib/x86_64-linux-gnu
#INCLUDE_PATH_SEQLIB      := -I /usr/include/SeqLib -I/usr/include/htslib -I/usr/include/jsoncpp
#LIBRARY_PATH_SEQLIB      := -L/usr/lib/x86_64-linux-gnu -Wl,-rpath,/usr/lib/x86_64-linux-gnu

# alternatively add paths with header files to INCLS and paths with library files to LDFLAGS

ifeq ($(shell uname -s), Darwin)
	# path for default homebrew installation of lp_solve
	INCLUDE_PATH_LPSOLVE = -I/usr/local/opt/lp_solve/include
	# path for default homebrew installation of mysql and mysql++
	INCLUDE_PATH_MYSQL = -I/usr/local/opt/mysql/include/mysql -I/usr/local/opt/mysql++/include/mysql
	# path for default homebrew installation of bamtools
	INCLUDE_PATH_BAMTOOLS = -I/usr/local/opt/bamtools/include/bamtools
	# path for default homebrew installation of htslib
	INCLUDE_PATH_HTSLIB = -I/usr/local/opt/htslib/include/htslib
endif

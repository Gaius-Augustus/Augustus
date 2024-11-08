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
EBONY = true

# set these paths to the correct locations if you have installed the corresponding packages in non-default locations:
INCLUDE_PATH_ZLIB        := -I/home/saenkos/zlib/zlib_install/include
LIBRARY_PATH_ZLIB        := -L/home/saenkos/zlib/zlib_install/lib -Wl,-rpath,/home/saenkos/zlib/zlib_install/lib
INCLUDE_PATH_BOOST       := -I/home/saenkos/boost/boost_build/boost_1_76_0/include/
LIBRARY_PATH_BOOST       := -L/home/saenkos/boost/boost_build/boost_1_76_0/lib/ -Wl,-rpath,/home/saenkos/boost/boost_build/boost_1_76_0/lib/
INCLUDE_PATH_LPSOLVE     := -I/home/saenkos/lpsolve/lpsolve_install/include
LIBRARY_PATH_LPSOLVE     := -L/home/saenkos/lpsolve/lpsolve_install/lib -Wl,-rpath,/home/saenkos/lpsolve/lpsolve_install/lib
INCLUDE_PATH_SUITESPARSE := -I/home/saenkos/suitesparse/include
LIBRARY_PATH_SUITESPARSE := -L//home/saenkos/suitesparse/lib -Wl,-rpath,/home/saenkos/suitesparse/lib
INCLUDE_PATH_GSL         := -I/home/saenkos/gsl/gsl_install/include
LIBRARY_PATH_GSL         := -L/home/saenkos/gsl/gsl_install/lib -Wl,-rpath,/home/saenkos/gsl/gsl_install/lib
INCLUDE_PATH_MYSQL       := -I/home/saenkos/mysql/mysql_install/usr/include -I/home/saenkos/mysql/mysql_install/usr/include/mysql -I/home/saenkos/mysql/mysql_install/usr/include/mysql++     # the path to mysql++ may have to be adjusted
LIBRARY_PATH_MYSQL       := -L/home/saenkos/mysql/mysql_install/usr/lib/x86_64-linux-gnu -Wl,-rpath,/home/saenkos/mysql/mysql_install/usr/lib/x86_64-linux-gnu
INCLUDE_PATH_BAMTOOLS    := -I/home/saenkos/bamtools/bamtools_install/include/bamtools
LIBRARY_PATH_BAMTOOLS    := -L/home/saenkos/bamtools/bamtools_install/lib -Wl,-rpath,/home/saenkos/bamtools/bamtools_install/lib
INCLUDE_PATH_HTSLIB      := -I/home/saenkos/htslib/include/htslib
LIBRARY_PATH_HTSLIB      := -L/home/saenkos/htslib/lib -Wl,-rpath,/home/saenkos/htslib/lib
INCLUDE_PATH_SEQLIB      := -I /home/saenkos/SeqLib -I/home/saenkos/SeqLib/htslib
LIBRARY_PATH_SEQLIB      := -L//home/saenkos/SeqLib/lib -Wl,-rpath,/home/saenkos/SeqLib/lib
OPT = -g
CXXFLAGS = -DDEBUG -g -ggdb -pg -DDEBUG_STATES


INCLUDE_PATH_SQLITE := -I/home/saenkos/sqlite/sqlite3_install/include
LIBRARY_PATH_SQLITE := -L/home/saenkos/sqlite/sqlite3_install/lib -Wl,-rpath,/home/saenkos/sqlite/sqlite3_install/lib
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

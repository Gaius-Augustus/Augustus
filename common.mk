# Definitions common to all Makefiles
# This file is included from other Makefiles in the augustus project.
AUGVERSION = 3.3.3

# set ZIPINPUT to false if you do not require input gzipped input genome files,
# get compilation errors about the boost iostreams library or
# the required libraries libboost-iostreams-dev and lib1g-dev are not available
ZIPINPUT = true

<<<<<<< HEAD
# set COMPGENEPRED to false if you do not require the comparative gene prediction mode (CGP) or
# the required libraries libgsl-dev, libboost-all-dev, libsuitesparse-dev and liblpsolve55-dev are not available
COMPGENEPRED = true
=======
# uncomment this line to enable comparative gene finding (requires compiler which supports C++11 standard)
COMPGENEPRED = true

# uncomment this line when you need MySQL access to sequences (most users don't)
#MYSQL = true

# uncomment this line to enable access to SQLite databases that store
# file offsets of sequence data in flat files and hints
SQLITE = true
>>>>>>> 0aaefc96961df96de4c4269dd9a2b2d5eb028909

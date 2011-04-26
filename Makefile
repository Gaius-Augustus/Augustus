#
# Makefile for Augustus
#

all: 
	cd src && ${MAKE}
	cd scripts && ${MAKE}

clean: 
	cd src && ${MAKE} clean
	cd scripts && ${MAKE} clean

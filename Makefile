#
# Makefile for Augustus
#
include common.mk

all: 
	cd src && ${MAKE}
	cd scripts && ${MAKE}

clean: 
	cd src && ${MAKE} clean
	cd scripts && ${MAKE} clean

# for internal purposes:
release: 
	find . -name .svn | xargs rm -rf
	find . -name "*~" | xargs rm -f
	rm -rf scripts/compileSpliceCands
	make clean all
	make clean
	rm -f src/.kdbgrc*      
	rm -f src/makedepend.pl
	cd auxprogs/filterBam/; make clean all; cd -
	cd auxprogs/bam2hints; make clean; make ; cd -
	cd config/species; rm -rf tobacco xeno1 bombus_terrestris{1,3}
	tar -czf ../augustus.2.7.tar.gz .

# remove -static from src/Makefile for MAC users

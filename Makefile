#
# Makefile for Augustus
#

all: 
	cd src && ${MAKE}
	cd scripts && ${MAKE}

clean: 
	cd src && ${MAKE} clean
	cd scripts && ${MAKE} clean

release: 
	find . -name .svn | xargs rm -rf
	find . -name "*~" | xargs rm -f
	rm -rf scripts/bam2wig
	rm -rf scripts/fixed2variableStep
	rm -rf scripts/compileSpliceCands
	make clean all
	cd src; cp augustus etraining fastBlockSearch prepareAlign ../bin/; cd ..
	make clean
	rm -f src/.kdbgrc*      
	rm -f src/makedepend.pl
	rm -f auxprogs/filterBam/*.o
	cd auxprogs/bam2hints; make clean; cd -
	cd config/species; rm -rf tobacco xeno1 bombus_terrestris{1,3}
	tar -czf ../augustus.2.6.1.tar.gz .

# remove -static from Makefile for MAC users

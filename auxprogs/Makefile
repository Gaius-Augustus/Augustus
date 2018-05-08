# Makefile for auxprogs
#
# Created on: 16-October-2016
# Last modified on: 16-October-2016 by Katharina Hoff (katharina.hoff@gmail.com)

all:
	cd bam2hints; make;
	cd compileSpliceCands; make;
	cd filterBam; make;
	cd homGeneMapping; make;
	cd joingenes; make;
	cd bam2wig; make;
	cd utrrnaseq/Debug; make -B all;
clean:
	cd bam2hints; make clean;
	cd compileSpliceCands; make clean;
	cd filterBam; make clean;
	cd homGeneMapping; make clean;
	cd joingenes; make clean;
	cd bam2wig; make clean;
	cd utrrnaseq/Debug; make clean;

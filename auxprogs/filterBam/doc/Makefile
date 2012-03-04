# Generates documentation of filterBAM
#
# Created: 28-February-2012
# Last modified: 4-March-2012

######################################################
# Generation of PDF document
######################################################


all : filterBam.pdf

filterBam.pdf : filterBam.dvi
	dvipdf $<

filterBam.dvi:filterBam.tex
	latex $< > filterBam.log
	latex $< > filterBam.log
	@echo "filterBam.tex compilation results being stored in filterBam.log"

clean : 
	rm -f *~ *.aux *.blg *.bbl *.log *.out *.dvi filterBam.log

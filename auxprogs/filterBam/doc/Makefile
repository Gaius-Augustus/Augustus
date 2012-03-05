# Generates documentation of filterBAM
#
# Created: 28-February-2012
# Last modified: 4-March-2012

all: compile

have_pdf := $(wildcard filterBam.pdf)

ifeq ($(strip $(have_pdf)),)

compile : filterBam.pdf
	@echo "filterBam.pdf does not exist. Compiling...";

filterBam.pdf : filterBam.dvi
	dvipdf $<

filterBam.dvi:filterBam.tex
	latex $< #> filterBam.log
	latex $< #> filterBam.log
	@echo "filterBam.tex compilation results being stored in filterBam.log"

else 

compile : 
	@echo "filterBam.pdf exists.";

endif


clean : 
	rm -f *~ *.aux *.blg *.bbl *.log *.out *.dvi filterBam.log

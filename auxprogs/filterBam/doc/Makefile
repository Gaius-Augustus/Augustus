# Generates documentation of filterBAM
#
# Created: 28-February-2012
# Last modified: 28-February-2012

######################################################
# Generation of PDF document
######################################################

PROGRAM=filterBam

doc : $(PROGRAM).tex
	pdflatex $^
	pdflatex $^ 

PDFOBJECTS = $(PROGRAM).aux $(PROGRAM).blg $(PROGRAM).bbl $(PROGRAM).log $(PROGRAM).out


clean : 
	rm -f *~ $(PDFOBJECTS)
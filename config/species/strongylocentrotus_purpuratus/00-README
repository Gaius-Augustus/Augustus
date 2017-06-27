ABOUT:
This training set was created by Joseph Ryan using version 3.1 of the Strongylocentrotus purpuratus gene models from Echinobase using the procedures outlined here: http://www.molecularevolution.org/molevolfiles/exercises/augustus/training.html

INSTALL:
move the strongylocentrotus_purpuratus directory into the augustus/config/species directory

RUN: 
once installed, add --species=strongylocentrotus_purpuratus to the augustus cmd line

# notes below from training
SEE: http://bioinf.uni-greifswald.de/augustus/binaries/tutorial/scipio.html
genes.gb contains 88 gene structures for training. Website says: As a minimum for good performance you should have 200 gene structures for training

# these files were downloaded from Echinobase

cp -R $AUGUSTUS_CONFIG_PATH/species/generic $AUGUSTUS_CONFIG_PATH/species/strongylocentrotus_purpuratus

mkdir 01-TRAIN
cd 01-TRAIN/
mv ~/Spur3.1.fa ~/Transcriptome.gff3 .
/usr/local/augustus-3.0.3/scripts/gff2gbSmallDNA.pl Transcriptome.gff3 Spur3.1.fa 1000 genes.raw.gb
/usr/local/augustus-3.0.3/bin/etraining --species=strongylocentrotus_purpuratus --stopCodonExcludedFromCDS=true genes.raw.gb 2> train.err
cat train.err | perl -pe 's/.*in sequence (\S+): .*/$1/' > badgenes.lst
cat badgenes.lst
locate filterGenes
/usr/local/augustus-3.0.3/scripts/filterGenes.pl badgenes.lst genes.raw.gb > genes.gb
grep -c "LOCUS" genes.raw.gb genes.gb

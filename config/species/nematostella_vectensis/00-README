ABOUT:
This training set was created by Joseph Ryan using version 1 of the Nematostella gene models using the procedures outlined here: http://www.molecularevolution.org/molevolfiles/exercises/augustus/training.html

INSTALL:
move the nematostella_vectensis directory into the augustus/config/species directory

RUN: 
once installed, add --species=nematostella_vectensis to the augustus cmd line

# notes from the training
lwp-download ftp://ftp.jgi-psf.org/pub/JGI_data/Nematostella_vectensis/v1.0/annotation/transcripts.Nemve1FilteredModels1.fasta.gz
lwp-download ftp://ftp.jgi-psf.org/pub/JGI_data/Nematostella_vectensis/v1.0/annotation/Nemve1.FilteredModels1.gff.gz

gzip -d transcripts.Nemve1FilteredModels1.fasta.gz

gzip -d Nemve1.FilteredModels1.gff.gz

/usr/local/augustus-3.0.3/scripts/gff2gbSmallDNA.pl  Nemve1.FilteredModels1.gff transcripts.Nemve1FilteredModels1.fasta 1000 genes.raw.gb

/usr/local/augustus-3.0.3/bin/etraining --species=Nematostella_vectensis --stopCodonExcludedFromCDS=true genes.raw.gb 2> train.err

cp /usr/local/augustus-3.0.3/config/species/generic/generic_parameters.cfg /usr/local/augustus-3.0.3/config/species/nematostella_vectensis/nematostella_vectensis.cfg

cp /usr/local/augustus-3.0.3/config/species/generic/generic_weightmatrix.txt /usr/local/augustus-3.0.3/config/species/nematostella_vectensis/nematostella_vectensis_weightmatrix.txt

mv /usr/local/augustus-3.0.3/config/species/nematostella_vectensis/nematostella_vectensis.cfg /usr/local/augustus-3.0.3/config/species/nematostella_vectensis/nematostella_vectensis_parameters.cfg

/usr/local/augustus-3.0.3/scripts/gff2gbSmallDNA.pl Nemve1.FilteredModels1.gff /data1/jfryan/00-DATA/NvecGenome 1000 genes.raw.gb

/usr/local/augustus-3.0.3/bin/etraining --species=nematostella_vectensis --stopCodonExcludedFromCDS=true genes.raw.gb 2> train.err

/usr/local/augustus-3.0.3/scripts/filterGenes.pl badgenes.lst genes.raw.gb > genes.gb






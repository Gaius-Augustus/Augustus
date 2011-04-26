This directory contains toy data to test the pipeline autoAug.pl for training and running
AUGUSTUS.

genome.fa          region 2000001..3000000 of chromosome I from C.elegans
cdna.fa            sample cDNA sequences (ESTs)
traingenes.gff     coordinates of nGASP training genes falling in region of genome.fa
                   coding regions (CDS) only
traingenes.gb      genes and sequences together in genbank format
                   this file was created from genome.fa and traingenes.gff with the command
		   gff2gbSmallDNA.pl traingenes.gff genome.fa 1000 traingenes.gb
cdna.psl           EST alignments in a format produced by BLAT and GMAP

1) Try whether AUGUSTUS and BLAT are installed correctly by executing


autoAug.pl --species=xyz --genome=genome.fa --trainingset=traingenes.gb --cdna=cdna.fa
	   --singleCPU -v -v -v --opt=0 --useexisting --noutr

This will run the training and prediction pipeline on a single machine without meta parameter optimization
and without training and predicting untranslated regions.


2) For a first test, whether you have installed everything correctly, including PASA, do

autoAug.pl --species=test --genome=genome.fa --cdna=cdna.fa --singleCPU -v -v -v --pasa --useexisting --opt=0

This will run the complete pipeline, including training set generation with PASA.
If something fails, you can fix the problem, delete any incomplete files and rerun the same command to resume.
This command takes about 40m on my computer to finish.


3) Here is an example for how to predict genes genome-wide using a cluster:

autoAugPred.pl --species=elegans -g genome.fa --hints=hints.E.gff

Then follow the commands, and submit the job(s) on your cluster. E.g. do

cd x/y/z/autoAugPred_hints/shells/
qsub aug1

When they are finished, then rerun

autoAugPred.pl --species=elegans -g genome.fa --hints=hints.E.gff --useexisting

from the original directory.
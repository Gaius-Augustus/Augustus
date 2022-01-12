# Katharina J. Hoff
# November 25th 2020
# Species specific AUGUSTUS parameters for Paramecium tetraurelia

Augustus parameter set name: Paramecium_tetraurelia

Lineage: Eukaryota; Sar; Alveolata; Ciliophora; Intramacronucleata; Oligohymenophorea; Peniculida; Parameciidae; Paramecium

This parameter set was created with BRAKER on the basis of genome assembly https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/165/425/GCF_000165425.1_ASM16542v1/GCF_000165425.1_ASM16542v1_genomic.fna.gz and OrthoDb proteins https://v100.orthodb.org/download/odb10_protozoa_fasta.tar.gz

Call:

braker.pl --genome=genome.fa --prot_seq=proteins.fasta --softmasking --cores=8 --species=Paramecium_tetraurelia --epmode

Flanking region size: 688
Number of generated training genes: 2750 nonredundant genes, 200 set aside for testing, 200 set aside for optimize_augustus.pl time consuming steps

Initial accuracy on test set: 0.282073333333333
Accuracy on test set after optimize_augustus.pl: 0.451706666666667

With this parameter set, BRAKER predicted 9896 transcripts (Augustus with hints in the target genome). For comparison, the RefSeq annotation contains 39580 transcripts.

This parameter set was generated in a fully automated fashion for supporing the POMPU project (http://www.pompu-project.de/). It probably heavily underpredicts genes.

If using this parameter set, please cite:

Bruna, T., Hoff, K.J., Stanke, M., Lomsadze, A., & Borodovsky, M. (2021). BRAKER2: Automatic Eukaryotic Genome Annotation with GeneMark-EP+ and AUGUSTUS Supported by a Protein Database. NAR genomics and bioinformatics, lqaa108.

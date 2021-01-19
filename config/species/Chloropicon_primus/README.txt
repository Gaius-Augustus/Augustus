# Katharina J. Hoff
# November 25th 2020
# Species specific AUGUSTUS parameters for Chloropicon primus

Augustus parameter set name: Chloropicon_primus

Lineage: Eukaryota; Viridiplantae; Chlorophyta; Chloropicophyceae; Chloropicales; Chloropicaceae; Chloropicon

This parameter set was created with BRAKER on the basis of genome assembly https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/007/859/695/GCA_007859695.1_ASM785969v1/GCA_007859695.1_ASM785969v1_genomic.fna.gz and OrthoDb proteins https://v100.orthodb.org/download/odb10_plants_fasta.tar.gz

Call:

braker.pl --genome=genome.fa --prot_seq=proteins.fasta --softmasking --cores=8 --species=Chloropicon_primus --epmode

Flanking region size: 801
Number of generated training genes: 1213 nonredundant genes, 300 set aside for testing, 300 set aside for optimize_augustus.pl time consuming steps

Initial accuracy on test set: 0.842
Accuracy on test set after optimize_augustus.pl: 0.858466666666667

With this parameter set, BRAKER predicted 9582 genes (Augustus with hints in the target genome). For comparison, the GenBank annotation contains 8639 genes.

This parameter set was generated in a fully automated fashion for supporing the POMPU project (http://www.pompu-project.de/).

If using this parameter set, please cite:

Bruna, T., Hoff, K., Stanke, M., Lomsadze, A., & Borodovsky, M. (2020). BRAKER2: Automatic Eukaryotic Genome Annotation with GeneMark-EP+ and AUGUSTUS Supported by a Protein Database. bioRxiv. https://www.biorxiv.org/content/10.1101/2020.08.10.245134v1.abstract

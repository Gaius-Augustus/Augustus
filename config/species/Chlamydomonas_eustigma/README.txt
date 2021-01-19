# Katharina J. Hoff
# November 25th 2020
# Species specific AUGUSTUS parameters for Chlamydomonas eustigma

Augustus parameter set name: Chlamydomonas_eustigma

Lineage: Eukaryota; Viridiplantae; Chlorophyta; core chlorophytes; Chlorophyceae; Chlamydomonadales; Chlamydomonadaceae; Chlamydomonas

This parameter set was created with BRAKER on the basis of genome assembly https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/335/675/GCA_002335675.1_C.eustigma_genome_v1.0/GCA_002335675.1_C.eustigma_genome_v1.0_genomic.fna.gz and OrthoDb proteins https://v100.orthodb.org/download/odb10_plants_fasta.tar.gz

Call:

braker.pl --genome=genome.fa --prot_seq=proteins.fasta --softmasking --cores=8 --species=Chlamydomonas_eustigma --epmode

Flanking region size: 1937
Number of generated training genes: 5115 nonredundant genes, 300 set aside for testing, 300 set aside for optimize_augustus.pl time consuming steps

Initial accuracy on test set: 0.8192
Accuracy on test set after optimize_augustus.pl: 0.8076 (? it slightly decreased)

With this parameter set, BRAKER predicted 15942 genes (Augustus with hints in the target genome). For comparison, the DDBJ annotation contains 14161 genes.

This parameter set was generated in a fully automated fashion for supporing the POMPU project (http://www.pompu-project.de/).

If using this parameter set, please cite:

Bruna, T., Hoff, K., Stanke, M., Lomsadze, A., & Borodovsky, M. (2020). BRAKER2: Automatic Eukaryotic Genome Annotation with GeneMark-EP+ and AUGUSTUS Supported by a Protein Database. bioRxiv. https://www.biorxiv.org/content/10.1101/2020.08.10.245134v1.abstract

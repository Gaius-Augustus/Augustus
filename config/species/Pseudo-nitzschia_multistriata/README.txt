# Katharina J. Hoff
# November 25th 2020
# Species specific AUGUSTUS parameters for Pseudo-nitzschia multistriata

Augustus parameter set name: Pseudo-nitzschia_multistriata

Lineage: Eukaryota; Sar; Stramenopiles; Ochrophyta; Bacillariophyta; Bacillariophyceae; Bacillariophycidae; Bacillariales; Bacillariaceae; Pseudo-nitzschia

This parameter set was created with BRAKER on the basis of genome assembly https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/660/405/GCA_900660405.1_ASM90066040v1/GCA_900660405.1_ASM90066040v1_genomic.fna.gz and OrthoDb proteins https://v100.orthodb.org/download/odb10_plants_fasta.tar.gz

Call:

braker.pl --genome=genome.fa --prot_seq=proteins.fasta --softmasking --cores=8 --species=Pseudo-nitzschia_multistriata --epmode

Flanking region size: 986
Number of generated training genes: 2435 nonredundant genes, 300 set aside for testing, 300 set aside for optimize_augustus.pl time consuming steps

Initial accuracy on test set: 0.428133333333333
Accuracy on test set after optimize_augustus.pl: 0.4622

With this parameter set, BRAKER predicted 13455 genes (Augustus with hints in the target genome). For comparison, the EMBL annotation contains 12039 genes.

This parameter set was generated in a fully automated fashion for supporing the POMPU project (http://www.pompu-project.de/).

If using this parameter set, please cite:

Bruna, T., Hoff, K., Stanke, M., Lomsadze, A., & Borodovsky, M. (2020). BRAKER2: Automatic Eukaryotic Genome Annotation with GeneMark-EP+ and AUGUSTUS Supported by a Protein Database. bioRxiv. https://www.biorxiv.org/content/10.1101/2020.08.10.245134v1.abstract

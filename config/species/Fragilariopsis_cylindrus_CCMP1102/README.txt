# Katharina J. Hoff
# November 25th 2020
# Species specific AUGUSTUS parameters for Fragilariopsis cylindrus CCMP1102

Augustus parameter set name: Fragilariopsis_cylindrus_CCMP1102

Lineage: Eukaryota; Sar; Stramenopiles; Ochrophyta; Bacillariophyta; Bacillariophyceae; Bacillariophycidae; Bacillariales; Bacillariaceae; Fragilariopsis

This parameter set was created with BRAKER on the basis of genome assembly https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/750/085/GCA_001750085.1_Fracy1/GCA_001750085.1_Fracy1_genomic.fna.gz and OrthoDb proteins https://v100.orthodb.org/download/odb10_plants_fasta.tar.gz

Call:

braker.pl --genome=genome.fa --prot_seq=proteins.fasta --softmasking --cores=8 --species=Fragilariopsis_cylindrus_CCMP1102 --epmode

Flanking region size: 895
Number of generated training genes: 1233 nonredundant genes, 300 set aside for testing, 300 set aside for optimize_augustus.pl time consuming steps

Initial accuracy on test set: 0.874733333333333
Accuracy on test set after optimize_augustus.pl: 0.940866666666667

With this parameter set, BRAKER predicted 26889 genes (Augustus with hints in the target genome). For comparison, the GenBank annotation contains 18111 genes.

This parameter set was generated in a fully automated fashion for supporing the POMPU project (http://www.pompu-project.de/).

If using this parameter set, please cite:

Bruna, T., Hoff, K., Stanke, M., Lomsadze, A., & Borodovsky, M. (2020). BRAKER2: Automatic Eukaryotic Genome Annotation with GeneMark-EP+ and AUGUSTUS Supported by a Protein Database. bioRxiv. https://www.biorxiv.org/content/10.1101/2020.08.10.245134v1.abstract

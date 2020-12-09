# Katharina J. Hoff
# November 25th 2020
# Species specific AUGUSTUS parameters for Thalassiosira pseudonana CCMP1335

Augustus parameter set name: Thalassiosira_pseudonana_CCMP1335

Lineage: Eukaryota; Sar; Stramenopiles; Ochrophyta; Bacillariophyta; Coscinodiscophyceae; Thalassiosirophycidae; Thalassiosirales; Thalassiosiraceae; Thalassiosira

This parameter set was created with BRAKER on the basis of genome assembly https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/149/405/GCF_000149405.2_ASM14940v2/GCF_000149405.2_ASM14940v2_genomic.fna.gz and OrthoDb proteins https://v100.orthodb.org/download/odb10_plants_fasta.tar.gz

Call:

braker.pl --genome=genome.fa --prot_seq=proteins.fasta --softmasking --cores=8 --species=Thalassiosira_pseudonana_CCMP1335 --epmode

Flanking region size: 935
Number of generated training genes: 1861 nonredundant genes, 300 set aside for testing, 300 set aside for optimize_augustus.pl time consuming steps

Initial accuracy on test set: 0.840333333333333
Accuracy on test set after optimize_augustus.pl: 0.844666666666667

With this parameter set, BRAKER predicted 13073 genes (Augustus with hints in the target genome). For comparison, the GenBank annotation contains 11673 genes.

This parameter set was generated in a fully automated fashion for supporing the POMPU project (http://www.pompu-project.de/).

If using this parameter set, please cite:

Bruna, T., Hoff, K., Stanke, M., Lomsadze, A., & Borodovsky, M. (2020). BRAKER2: Automatic Eukaryotic Genome Annotation with GeneMark-EP+ and AUGUSTUS Supported by a Protein Database. bioRxiv. https://www.biorxiv.org/content/10.1101/2020.08.10.245134v1.abstract

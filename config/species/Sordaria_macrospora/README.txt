# Katharina J. Hoff
# November 25th 2020
# Species specific AUGUSTUS parameters for Sordaria macrospora

Augustus parameter set name: Sordaria_macrospora

Lineage: 

This parameter set was created with BRAKER on the basis of genome assembly https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/182/805/GCF_000182805.2_ASM18280v2/GCF_000182805.2_ASM18280v2_genomic.fna.gz and OrthoDb proteins https://v100.orthodb.org/download/odb10_fungi_fasta.tar.gz

Call:

braker.pl --genome=genome.fa --prot_seq=proteins.fasta --softmasking --cores=8 --species=Sordaria macrospora --epmode

Flanking region size: 878
Number of generated training genes: 3315 nonredundant genes, 300 set aside for testing, 300 set aside for optimize_augustus.pl time consuming steps

Initial accuracy on test set: 0.852066666666667
Accuracy on test set after optimize_augustus.pl: 0.851933333333333

With this parameter set, BRAKER predicted 10837 genes (Augustus with hints in the target genome). For comparison, the RefSeq annotation contains 9896 genes.

This parameter set was generated in a fully automated fashion for supporing the POMPU project (http://www.pompu-project.de/).

If using this parameter set, please cite:

Bruna, T., Hoff, K., Stanke, M., Lomsadze, A., & Borodovsky, M. (2020). BRAKER2: Automatic Eukaryotic Genome Annotation with GeneMark-EP+ and AUGUSTUS Supported by a Protein Database. bioRxiv. https://www.biorxiv.org/content/10.1101/2020.08.10.245134v1.abstract

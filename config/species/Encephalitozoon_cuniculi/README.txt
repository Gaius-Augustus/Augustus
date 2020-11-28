# Katharina J. Hoff
# November 28th 2020
# Species specific AUGUSTUS parameters for Encephalitozoon cuniculi

Augustus parameter set name: Encephalitozoon_cuniculi

Lineage: Eukaryota; Opisthokonta; Fungi; Fungi incertae sedis; Microsporidia; Apansporoblastina; Unikaryonidae; Encephalitozoon


This parameter set was created with BRAKER on the basis of genome assembly https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/091/225/GCF_000091225.1_ASM9122v1/GCF_000091225.1_ASM9122v1_genomic.fna.gz and protein file https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/091/225/GCF_000091225.1_ASM9122v1/GCF_000091225.1_ASM9122v1_protein.faa.gz

Call:

braker.pl --species=Encephalitozoon_cuniculi --genome=genome.fa --prot_seq=proteins.fa --prg=gth --trainFromGth --softmasking --cores=4

Flanking region size: 544
Number of generated training genes: 1694 nonredundant genes, 300 set aside for testing, 300 set aside for optimize_augustus.pl time consuming steps

Initial accuracy on test set: 0.726933333333333
Accuracy on test set after optimize_augustus.pl: 0.805466666666667

With this parameter set, BRAKER predicted 2310 genes (Augustus with hints in the target genome). For comparison, the RefSeq annotation contains 1996 genes.

This parameter set was generated in a fully automated fashion for supporing the POMPU project (http://www.pompu-project.de/).

If using this parameter set, please cite:

Bruna, T., Hoff, K., Stanke, M., Lomsadze, A., & Borodovsky, M. (2020). BRAKER2: Automatic Eukaryotic Genome Annotation with GeneMark-EP+ and AUGUSTUS Supported by a Protein Database. bioRxiv. https://www.biorxiv.org/content/10.1101/2020.08.10.245134v1.abstract

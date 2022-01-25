# Katharina J. Hoff
# January 25th 2022
# Species specific AUGUSTUS parameters for Berviolum minutum

Augustus parameter set name: Berviolum_minutum

Lineage: Eukaryota; Sar; Alveolata; Dinophyceae; Suessiales; Symbiodiniaceae; Breviolum

This parameter set was created with BRAKER on the basis of genome https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/507/305/GCA_000507305.1_ASM50730v1/GCA_000507305.1_ASM50730v1_genomic.fna.gz and OrthoDb proteins https://v100.orthodb.org/download/odb10_protozoa_fasta.tar.gz

Call:

braker.pl --genome=genome.fa --prot_seq=proteins.fasta --softmasking --cores=8 --species=Berviolum_minutum --epmode

Flanking region size: 2279
Number of generated training genes: 3475 nonredundant genes, 300 set aside for testing, 300 set aside for optimize_augustus.pl time consuming steps

Initial accuracy on test set: 0.47988
Accuracy on test set after optimize_augustus.pl: 0.510793333333333

Note: we have not generated an annotation for Berviolum minutum on the basis of this parameter set, yet.

This parameter set was generated in a fully automated fashion for supporing the POMPU project (http://www.pompu-project.de/).

If using this parameter set, please cite:

Bruna, T., Hoff, K.J., Stanke, M., Lomsadze, A., & Borodovsky, M. (2021). BRAKER2: Automatic Eukaryotic Genome Annotation with GeneMark-EP+ and AUGUSTUS Supported by a Protein Database. NAR genomics and bioinformatics, lqaa108.

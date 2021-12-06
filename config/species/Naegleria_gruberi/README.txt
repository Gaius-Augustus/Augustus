# Katharina J. Hoff
# August 25th 2021
# Species specific AUGUSTUS parameters for Naegleria gruberi

Augustus parameter set name: Naegleria_gruberi

Lineage: Eukaryota; Discoba; Heterolobosea; Tetramitia; Eutetramitia; Vahlkampfiidae; Naegleria

This parameter set was created with BRAKER on the basis of genome assembly https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/985/GCF_000004985.1_V1.0/GCF_000004985.1_V1.0_genomic.fna.gz and OrthoDb proteins https://v100.orthodb.org/download/odb10_protozoa_fasta.tar.gz

Call:

braker.pl --genome=genome.fa --prot_seq=proteins.fasta --softmasking --cores=8 --species=Naegleria_gruberi --epmode

Flanking region size: 816
Number of generated training genes: 3278 nonredundant genes, 200 set aside for testing, 200 set aside for optimize_augustus.pl time consuming steps

Initial accuracy on test set: 0.8546
Accuracy on test set after optimize_augustus.pl: 0.910333333333333

With this parameter set, BRAKER predicted 19688 transcripts (Augustus with hints in the target genome). For comparison, the RefSeq annotation contains 15711 transcripts.

This parameter set was generated in a fully automated fashion for supporing the POMPU project (http://www.pompu-project.de/).

If using this parameter set, please cite:

Bruna, T., Hoff, K.J., Stanke, M., Lomsadze, A., & Borodovsky, M. (2021). BRAKER2: Automatic Eukaryotic Genome Annotation with GeneMark-EP+ and AUGUSTUS Supported by a Protein Database. NAR genomics and bioinformatics, lqaa108.

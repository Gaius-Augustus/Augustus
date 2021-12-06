# Katharina J. Hoff
# August 25th 2021
# Species specific AUGUSTUS parameters for Monosiga brevicollis

Augustus parameter set name: Monosiga_brevicollis

Lineage: Eukaryota; Opisthokonta; Choanoflagellata; Craspedida; Salpingoecidae; Monosiga

This parameter set was created with BRAKER on the basis of genome assembly https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/865/GCF_000002865.3_V1.0/GCF_000002865.3_V1.0_genomic.fna.gz and OrthoDb proteins https://v100.orthodb.org/download/odb10_plants_fasta.tar.gz

Call:

braker.pl --genome=genome.fa --prot_seq=proteins.fasta --softmasking --cores=8 --species=Monosiga_brevicollis --epmode

Flanking region size: 1405
Number of generated training genes: 3553 nonredundant genes, 200 set aside for testing, 200 set aside for optimize_augustus.pl time consuming steps

Initial accuracy on test set: 0.825
Accuracy on test set after optimize_augustus.pl: 0.831

With this parameter set, BRAKER predicted 11381 transcripts (Augustus with hints in the target genome). For comparison, the RefSeq annotation contains 9202 transcripts.

This parameter set was generated in a fully automated fashion for supporing the POMPU project (http://www.pompu-project.de/).

If using this parameter set, please cite:

Bruna, T., Hoff, K.J., Stanke, M., Lomsadze, A., & Borodovsky, M. (2021). BRAKER2: Automatic Eukaryotic Genome Annotation with GeneMark-EP+ and AUGUSTUS Supported by a Protein Database. NAR genomics and bioinformatics, lqaa108.

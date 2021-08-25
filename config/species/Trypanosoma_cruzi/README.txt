# Katharina J. Hoff
# August 25th 2021
# Species specific AUGUSTUS parameters for Trypanosoma cruzi

Augustus parameter set name: Trypanosoma_cruzi

Lineage: Eukaryota; Discoba; Euglenozoa; Kinetoplastea; Metakinetoplastina; Trypanosomatida; Trypanosomatidae; Trypanosoma; Schizotrypanum

This parameter set was created with BRAKER on the basis of genome assembly https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/209/065/GCF_000209065.1_ASM20906v1/GCF_000209065.1_ASM20906v1_genomic.fna.gz and OrthoDb proteins https://v100.orthodb.org/download/odb10_protozoa_fasta.tar.gz

Call:

braker.pl --genome=genome.fa --prot_seq=proteins.fasta --softmasking --cores=8 --species=Trypanosoma_cruzi --epmode

Flanking region size: 777
Number of generated training genes: 1705 nonredundant genes, 200 set aside for testing, 200 set aside for optimize_augustus.pl time consuming steps

Initial accuracy on test set: 0.849933333333334
Accuracy on test set after optimize_augustus.pl: 0.905

With this parameter set, BRAKER predicted 43221 transcripts (Augustus with hints in the target genome). For comparison, the RefSeq annotation contains 19607 transcripts.

This parameter set was generated in a fully automated fashion for supporing the POMPU project (http://www.pompu-project.de/).

If using this parameter set, please cite:

Bruna, T., Hoff, K.J., Stanke, M., Lomsadze, A., & Borodovsky, M. (2021). BRAKER2: Automatic Eukaryotic Genome Annotation with GeneMark-EP+ and AUGUSTUS Supported by a Protein Database. NAR genomics and bioinformatics, lqaa108.

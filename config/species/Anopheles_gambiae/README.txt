# Katharina J. Hoff
# August 26th 2021
# Species specific AUGUSTUS parameters for Anopheles gambiae

Augustus parameter set name: Anopheles_gambiae

Lineage: 

This parameter set was created with BRAKER on the basis of genome assembly https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/575/GCF_000005575.2_AgamP3/GCF_000005575.2_AgamP3_genomic.fna.gz and OrthoDb proteins https://v100.orthodb.org/download/odb10_arthropoda_fasta.tar.gz

Call:

braker.pl --genome=genome.fa --prot_seq=proteins.fasta --softmasking --cores=8 --species=Anopheles_gambiae --epmode

Flanking region size: 1889
Number of generated training genes: 5302 nonredundant genes, 300 set aside for testing, 300 set aside for optimize_augustus.pl time consuming steps

Initial accuracy on test set: 0.859866666666667
Accuracy on test set after optimize_augustus.pl: 0.8722

With this parameter set, BRAKER predicted 16607 transcripts (Augustus with hints in the target genome). For comparison, the RefSeq assembly has 14102 transcripts.

This parameter set was generated in a fully automated fashion for supporing the POMPU project (http://www.pompu-project.de/).

If using this parameter set, please cite:

Bruna, T., Hoff, K.J., Stanke, M., Lomsadze, A., & Borodovsky, M. (2021). BRAKER2: Automatic Eukaryotic Genome Annotation with GeneMark-EP+ and AUGUSTUS Supported by a Protein Database. NAR genomics and bioinformatics, lqaa108.

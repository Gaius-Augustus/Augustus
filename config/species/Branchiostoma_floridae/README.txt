# Katharina J. Hoff
# August 26th 2021
# Species specific AUGUSTUS parameters for Branchiostoma floridae

Augustus parameter set name: Branchiostoma_floridae

Lineage: Eukaryota; Opisthokonta; Metazoa; Eumetazoa; Bilateria; Deuterostomia; Chordata; Cephalochordata; Leptocardii; Amphioxiformes; Branchiostomidae; Branchiostoma

This parameter set was created with BRAKER on the basis of genome assembly https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/003/815/GCF_000003815.2_Bfl_VNyyK/GCF_000003815.2_Bfl_VNyyK_genomic.fna.gz and OrthoDb proteins https://v100.orthodb.org/download/odb10_vertebrata_fasta.tar.gz

Call:

braker.pl --genome=genome.fa --prot_seq=proteins.fasta --softmasking --cores=8 --species=Branchiostoma_floridae --epmode

Flanking region size: 3040
Number of generated training genes: 3085 nonredundant genes, 300 set aside for testing, 300 set aside for optimize_augustus.pl time consuming steps

Initial accuracy on test set: 0.7474
Accuracy on test set after optimize_augustus.pl: 0.761133333333333

With this parameter set, BRAKER predicted 38641 transcripts (Augustus with hints in the target genome). For comparison, the RefSeq assembly has 43041 transcripts.

This parameter set was generated in a fully automated fashion for supporing the POMPU project (http://www.pompu-project.de/).

If using this parameter set, please cite:

Bruna, T., Hoff, K.J., Stanke, M., Lomsadze, A., & Borodovsky, M. (2021). BRAKER2: Automatic Eukaryotic Genome Annotation with GeneMark-EP+ and AUGUSTUS Supported by a Protein Database. NAR genomics and bioinformatics, lqaa108.

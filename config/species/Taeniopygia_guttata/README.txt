# Katharina J. Hoff
# August 26th 2021
# Species specific AUGUSTUS parameters for Taeniopygia guttata

Augustus parameter set name: Taeniopygia_guttata

Lineage: Eukaryota; Opisthokonta; Metazoa; Eumetazoa; Bilateria; Deuterostomia; Chordata; Craniata; Vertebrata; Gnathostomata; Teleostomi; Euteleostomi; Sarcopterygii; Dipnotetrapodomorpha; Tetrapoda; Amniota; Sauropsida; Sauria; Archelosauria; Archosauria; Dinosauria; Saurischia; Theropoda; Coelurosauria; Aves; Neognathae; Passeriformes; Passeroidea; Estrildidae; Estrildinae; Taeniopygia

This parameter set was created with BRAKER on the basis of genome assembly https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/957/565/GCF_003957565.2_bTaeGut1.4.pri/GCF_003957565.2_bTaeGut1.4.pri_genomic.fna.gz and OrthoDb proteins https://v100.orthodb.org/download/odb10_vertebrata_fasta.tar.gz

Call:

braker.pl --genome=genome.fa --prot_seq=proteins.fasta --softmasking --cores=8 --species=Taeniopygia_guttata --epmode

Flanking region size: 2632
Number of generated training genes: 2853 nonredundant genes, 300 set aside for testing, 300 set aside for optimize_augustus.pl time consuming steps

Initial accuracy on test set: 0.362106666666667
Accuracy on test set after optimize_augustus.pl: 0.377813333333333

With this parameter set, BRAKER predicted 39928 transcripts (Augustus with hints in the target genome). For comparison, the RefSeq assembly has 41214 transcripts.

This parameter set was generated in a fully automated fashion for supporing the POMPU project (http://www.pompu-project.de/).

If using this parameter set, please cite:

Bruna, T., Hoff, K.J., Stanke, M., Lomsadze, A., & Borodovsky, M. (2021). BRAKER2: Automatic Eukaryotic Genome Annotation with GeneMark-EP+ and AUGUSTUS Supported by a Protein Database. NAR genomics and bioinformatics, lqaa108.

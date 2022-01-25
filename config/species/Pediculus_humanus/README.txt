# Katharina J. Hoff
# January 25th 2022
# Species specific AUGUSTUS parameters for Pediculus humanus

Augustus parameter set name: Pediculus_humanus

Lineage: Eukaryota; Opisthokonta; Metazoa; Eumetazoa; Bilateria; Protostomia; Ecdysozoa; Panarthropoda; Arthropoda; Mandibulata; Pancrustacea; Hexapoda; Insecta; Dicondylia; Pterygota; Neoptera; Paraneoptera; Psocodea; Phthiraptera; Anoplura; Pediculidae; Pediculus

This parameter set was created with BRAKER on the basis of genome https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/295/GCF_000006295.1_JCVI_LOUSE_1.0/GCF_000006295.1_JCVI_LOUSE_1.0_genomic.fna.gz and OrthoDb proteins https://v100.orthodb.org/download/odb10_arthropoda_fasta.tar.gz

Call:

braker.pl --genome=genome.fa --prot_seq=proteins.fasta --softmasking --cores=8 --species=Pediculus_humanus --epmode

Flanking region size: 1863
Number of generated training genes: 3535 nonredundant genes, 300 set aside for testing, 300 set aside for optimize_augustus.pl time consuming steps

Initial accuracy on test set: 0.871466666666667
Accuracy on test set after optimize_augustus.pl: 0.890666666666667

With this parameter set, BRAKER predicted 13268 transcripts (Augustus with hints in the target genome). For comparison, the RefSeq annotation contains 10775 transcripts.

This parameter set was generated in a fully automated fashion for supporing the POMPU project (http://www.pompu-project.de/).

If using this parameter set, please cite:

Bruna, T., Hoff, K.J., Stanke, M., Lomsadze, A., & Borodovsky, M. (2021). BRAKER2: Automatic Eukaryotic Genome Annotation with GeneMark-EP+ and AUGUSTUS Supported by a Protein Database. NAR genomics and bioinformatics, lqaa108.

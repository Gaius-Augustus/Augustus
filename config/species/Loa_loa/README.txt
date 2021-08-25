# Katharina J. Hoff
# August 25th 2021
# Species specific AUGUSTUS parameters for Loa loa

Augustus parameter set name: Loa_loa

Lineage: Eukaryota; Opisthokonta; Metazoa; Eumetazoa; Bilateria; Protostomia; Ecdysozoa; Nematoda; Chromadorea; Rhabditida; Spirurina; Spiruromorpha; Filarioidea; Onchocercidae; Loa

This parameter set was created with BRAKER on the basis of genome assembly https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/183/805/GCF_000183805.2_Loa_loa_V3.1/GCF_000183805.2_Loa_loa_V3.1_genomic.fna.gz and OrthoDb proteins https://v100.orthodb.org/download/odb10_metazoa_fasta.tar.gz

Call:

braker.pl --genome=genome.fa --prot_seq=proteins.fasta --softmasking --cores=8 --species=Loa_loa --epmode

Flanking region size: 2179
Number of generated training genes: 3551 nonredundant genes, 200 set aside for testing, 200 set aside for optimize_augustus.pl time consuming steps

Initial accuracy on test set: 0.7862
Accuracy on test set after optimize_augustus.pl: 0.7996

With this parameter set, BRAKER predicted 14614 transcripts (Augustus with hints in the target genome). For comparison, the RefSeq annotation contains 15440 transcripts.

This parameter set was generated in a fully automated fashion for supporing the POMPU project (http://www.pompu-project.de/).

If using this parameter set, please cite:

Bruna, T., Hoff, K.J., Stanke, M., Lomsadze, A., & Borodovsky, M. (2021). BRAKER2: Automatic Eukaryotic Genome Annotation with GeneMark-EP+ and AUGUSTUS Supported by a Protein Database. NAR genomics and bioinformatics, lqaa108.

# Katharina J. Hoff
# August 26th 2021
# Species specific AUGUSTUS parameters for Cryptosporidium hominis

Augustus parameter set name: Cryptosporidium_hominis

Lineage: Eukaryota; Sar; Alveolata; Apicomplexa; Conoidasida; Coccidia; Eucoccidiorida; Eimeriorina; Cryptosporidiidae; Cryptosporidium

This parameter set was created with BRAKER on the basis of genome assembly https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/425/GCF_000006425.1_ASM642v1/GCF_000006425.1_ASM642v1_genomic.fna.gz and OrthoDb proteins https://v100.orthodb.org/download/odb10_protozoa_fasta.tar.gz

Call:

braker.pl --genome=genome.fa --prot_seq=proteins.fasta --softmasking --cores=8 --species=Cryptosporidium_hominis --epmode

Flanking region size: 775
Number of generated training genes: 1243 nonredundant genes, 300 set aside for testing, 300 set aside for optimize_augustus.pl time consuming steps

Initial accuracy on test set: 0.559733333333333
Accuracy on test set after optimize_augustus.pl: 0.5596

With this parameter set, BRAKER predicted 4735 transcripts (Augustus with hints in the target genome). For comparison, the RefSeq assembly has 3885 transcripts.

This parameter set was generated in a fully automated fashion for supporing the POMPU project (http://www.pompu-project.de/).

If using this parameter set, please cite:

Bruna, T., Hoff, K.J., Stanke, M., Lomsadze, A., & Borodovsky, M. (2021). BRAKER2: Automatic Eukaryotic Genome Annotation with GeneMark-EP+ and AUGUSTUS Supported by a Protein Database. NAR genomics and bioinformatics, lqaa108.

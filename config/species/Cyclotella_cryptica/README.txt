# Katharina J. Hoff
# August 25th 2021
# Species specific AUGUSTUS parameters for Cyclotella cryptica

Augustus parameter set name: Cyclotella_cryptica

Lineage: Eukaryota; Sar; Stramenopiles; Ochrophyta; Bacillariophyta; Coscinodiscophyceae; Thalassiosirophycidae; Thalassiosirales; Stephanodiscaceae; Cyclotella

This parameter set was created with BRAKER on the basis of genome assembly https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/013/187/285/GCA_013187285.1_ASM1318728v1/GCA_013187285.1_ASM1318728v1_genomic.fna.gz and OrthoDb proteins https://v100.orthodb.org/download/odb10_plants_fasta.tar.gz

Call:

braker.pl --genome=genome.fa --prot_seq=proteins.fasta --softmasking --cores=8 --species=Sclerotinia_sclerotiorum --epmode

Flanking region size: 988
Number of generated training genes: 2258 nonredundant genes, 200 set aside for testing, 200 set aside for optimize_augustus.pl time consuming steps

Initial accuracy on test set: 0.5982
Accuracy on test set after optimize_augustus.pl: 0.6696

With this parameter set, BRAKER predicted 28883 transcripts (Augustus with hints in the target genome). To our knowledge, this might be the first annotation. Annotation is available upon request from katharina.hoff@uni-greifswald.de

This parameter set was generated in a fully automated fashion for supporing the POMPU project (http://www.pompu-project.de/).

If using this parameter set, please cite:

Bruna, T., Hoff, K.J., Stanke, M., Lomsadze, A., & Borodovsky, M. (2021). BRAKER2: Automatic Eukaryotic Genome Annotation with GeneMark-EP+ and AUGUSTUS Supported by a Protein Database. NAR genomics and bioinformatics, lqaa108.

# Katharina J. Hoff
# January 25th 2022
# Species specific AUGUSTUS parameters for Skeletonema costatum

Augustus parameter set name: Skeletonema_costatum

Lineage: Eukaryota; Sar; Stramenopiles; Ochrophyta; Bacillariophyta; Coscinodiscophyceae; Thalassiosirophycidae; Thalassiosirales; Skeletonemataceae; Skeletonema

This parameter set was created with BRAKER on the basis of genome https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/018/806/925/GCA_018806925.1_FSU_Scostatum_1.0/GCA_018806925.1_FSU_Scostatum_1.0_genomic.fna.gz and OrthoDb proteins https://v100.orthodb.org/download/odb10_plants_fasta.tar.gz

Call:

braker.pl --genome=genome.fa --prot_seq=proteins.fasta --softmasking --cores=8 --species=Skeletonema_costatum --epmode

Flanking region size: 716
Number of generated training genes: 1876 nonredundant genes, 300 set aside for testing, 300 set aside for optimize_augustus.pl time consuming steps

Initial accuracy on test set: 0.599733333333333
Accuracy on test set after optimize_augustus.pl: 0.6378

With this parameter set, BRAKER predicted 18379 transcripts (Augustus with hints, first BRAKER iteration, in the target genome). This may to our knowlege be the first gene set for this organism. It is available upon request from katharina.hoff@uni-greifswald.de.

This parameter set was generated in a fully automated fashion for supporing the POMPU project (http://www.pompu-project.de/).

If using this parameter set, please cite:

Bruna, T., Hoff, K.J., Stanke, M., Lomsadze, A., & Borodovsky, M. (2021). BRAKER2: Automatic Eukaryotic Genome Annotation with GeneMark-EP+ and AUGUSTUS Supported by a Protein Database. NAR genomics and bioinformatics, lqaa108.

# Katharina J. Hoff
# January 25th 2022
# Species specific AUGUSTUS parameters for Hydra vulgaris

Augustus parameter set name: Hydra_vulgaris

Lineage: Eukaryota; Opisthokonta; Metazoa; Eumetazoa; Cnidaria; Hydrozoa; Hydroidolina; Anthoathecata; Aplanulata; Hydridae; Hydra

This parameter set was created with BRAKER on the basis of genome https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/095/GCF_000004095.1_Hydra_RP_1.0/GCF_000004095.1_Hydra_RP_1.0_genomic.fna.gz and OrthoDb proteins https://v100.orthodb.org/download/odb10_metazoa_fasta.tar.gz

Call:

braker.pl --genome=genome.fa --prot_seq=proteins.fasta --softmasking --cores=8 --species=Hydra_vulgaris --epmode

Flanking region size: 1598
Number of generated training genes: 2128 nonredundant genes, 300 set aside for testing, 300 set aside for optimize_augustus.pl time consuming steps

Initial accuracy on test set: 0.733066666666667
Accuracy on test set after optimize_augustus.pl: 0.774466666666667

With this parameter set, BRAKER predicted 108397 transcripts/~52109 genes (Augustus with hints). For comparison, the reference annotation contains 21990 genes. Warning: these parameter may lead to split gene predictions/tends to overpredict.

This parameter set was generated in a fully automated fashion for supporing the POMPU project (http://www.pompu-project.de/).

If using this parameter set, please cite:

Bruna, T., Hoff, K.J., Stanke, M., Lomsadze, A., & Borodovsky, M. (2021). BRAKER2: Automatic Eukaryotic Genome Annotation with GeneMark-EP+ and AUGUSTUS Supported by a Protein Database. NAR genomics and bioinformatics, lqaa108.

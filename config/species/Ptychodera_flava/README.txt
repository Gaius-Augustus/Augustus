# Katharina J. Hoff
# December 2nd 2020
# Species specific AUGUSTUS parameters for Ptychodera flava

Warning: This parameter set is suspected to lead to falsely split gene predictions.

Augustus parameter set name: Ptychodera_flava

Lineage: Eukaryota; Metazoa; Hemichordata; Enteropneusta; Ptychoderidae; Ptychodera

This parameter set was created with BRAKER on the basis of genome assembly https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/465/055/GCA_001465055.1_ptychodera_flava_version_1.0.14/GCA_001465055.1_ptychodera_flava_version_1.0.14_genomic.fna.gz and OrthoDb proteins https://v100.orthodb.org/download/odb10_metazoa_fasta.tar.gz

Call:

braker.pl --genome=genome.fa --prot_seq=proteins.fasta --softmasking --cores=8 --species=Ptychodera_flava --epmode

Flanking region size: 2673
Number of generated training genes: 2617 nonredundant genes, 300 set aside for testing, 300 set aside for optimize_augustus.pl time consuming steps

Initial accuracy on test set: 0.837933333333333
Accuracy on test set after optimize_augustus.pl: 0.863266666666667

With this parameter set, BRAKER predicted 72062 genes (Augustus with hints in the target genome). To our knowledge, this might be the first gene set for this species. Warning: the gene set either contains a huge number of repeat proteins, or the parameter set causes falsely split gene predictions. The gene set is available upon request from katharina.hoff@uni-greifswald.de

This parameter set was generated in a fully automated fashion for supporing the POMPU project (http://www.pompu-project.de/).

If using this parameter set, please cite:

Bruna, T., Hoff, K., Stanke, M., Lomsadze, A., & Borodovsky, M. (2020). BRAKER2: Automatic Eukaryotic Genome Annotation with GeneMark-EP+ and AUGUSTUS Supported by a Protein Database. bioRxiv. https://www.biorxiv.org/content/10.1101/2020.08.10.245134v1.abstract

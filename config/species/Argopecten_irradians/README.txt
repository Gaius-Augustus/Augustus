# Katharina J. Hoff
# December 2nd 2020
# Species specific AUGUSTUS parameters for Argopecten irradians

Warning: This parameter set is suspected to lead to falsely split gene predictions.

Augustus parameter set name: Argopecten_irradians

Lineage: Eukaryota; Metazoa; Spiralia; Lophotrochozoa; Mollusca; Bivalvia; Autobranchia; Pteriomorphia; Pectinida; Pectinoidea; Pectinidae; Argopecten

This parameter set was created with BRAKER on the basis of genome assembly https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/004/382/745/GCA_004382745.1_QAU_Airr_1.1/GCA_004382745.1_QAU_Airr_1.1_genomic.fna.gz and OrthoDb proteins https://v100.orthodb.org/download/odb10_metazoa_fasta.tar.gz

Call:

braker.pl --genome=genome.fa --prot_seq=proteins.fasta --softmasking --cores=8 --species=Argopecten_irradians --epmode

Flanking region size: 2542
Number of generated training genes: 4577 nonredundant genes, 300 set aside for testing, 300 set aside for optimize_augustus.pl time consuming steps

Initial accuracy on test set: 0.839866666666667
Accuracy on test set after optimize_augustus.pl: 0.8646

With this parameter set, BRAKER predicted 57682 genes (Augustus with hints in the target genome). To our knowledge, this might be the first gene set for this species. Warning: the gene set either contains a huge number of repeat proteins, or the parameter set causes falsely split gene predictions. The gene set is available upon request from katharina.hoff@uni-greifswald.de

This parameter set was generated in a fully automated fashion for supporing the POMPU project (http://www.pompu-project.de/).

If using this parameter set, please cite:

Bruna, T., Hoff, K., Stanke, M., Lomsadze, A., & Borodovsky, M. (2020). BRAKER2: Automatic Eukaryotic Genome Annotation with GeneMark-EP+ and AUGUSTUS Supported by a Protein Database. bioRxiv. https://www.biorxiv.org/content/10.1101/2020.08.10.245134v1.abstract

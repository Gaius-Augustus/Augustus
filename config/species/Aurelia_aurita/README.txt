# Katharina J. Hoff
# November 29th 2020
# Species specific AUGUSTUS parameters for Aurelia aurita

Augustus parameter set name: Aurelia_aurita

Lineage: Eukaryota; Opisthokonta; Metazoa; Eumetazoa; Cnidaria; Scyphozoa; Semaeostomeae; Ulmaridae; Aurelia

This parameter set was created with BRAKER on the basis of genome assembly https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/004/194/415/GCA_004194415.1_ABSv1/GCA_004194415.1_ABSv1_genomic.fna.gz and OrthoDb proteins https://v100.orthodb.org/download/odb10_metazoa_fasta.tar.gz

Call:

braker.pl --genome=genome.fa --prot_seq=proteins.fasta --softmasking --cores=8 --species=Aurelia_aurita --epmode

Flanking region size: 3440
Number of generated training genes: 2734 nonredundant genes, 300 set aside for testing, 300 set aside for optimize_augustus.pl time consuming steps

Initial accuracy on test set: 0.7782
Accuracy on test set after optimize_augustus.pl: 0.793933333333333

With this parameter set, BRAKER predicted 29078 genes (Augustus with hints in the target genome). To our knowledge, this might be the first gene set for this species. The gene set is available upon request from katharina.hoff@uni-greifswald.de.

This parameter set was generated in a fully automated fashion for supporing the POMPU project (http://www.pompu-project.de/).

If using this parameter set, please cite:

Bruna, T., Hoff, K., Stanke, M., Lomsadze, A., & Borodovsky, M. (2020). BRAKER2: Automatic Eukaryotic Genome Annotation with GeneMark-EP+ and AUGUSTUS Supported by a Protein Database. bioRxiv. https://www.biorxiv.org/content/10.1101/2020.08.10.245134v1.abstract

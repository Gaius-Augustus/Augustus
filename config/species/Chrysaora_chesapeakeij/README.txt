# Katharina J. Hoff
# November 30th 2020
# Species specific AUGUSTUS parameters for Chrysaora chesapeakei

# WARNING: this parameter set has low quality and will likely lead to poor results!

Augustus parameter set name: Chrysaora_chesapeakei

Lineage: Eukaryota; Opisthokonta; Metazoa; Eumetazoa; Cnidaria; Scyphozoa; Semaeostomeae; Pelagiidae; Chrysaora

This parameter set was created with BRAKER on the basis of genome assembly https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/011/763/395/GCA_011763395.1_ASM1176339v1/GCA_011763395.1_ASM1176339v1_genomic.fna.gz and OrthoDb proteins https://v100.orthodb.org/download/odb10_metazoa_fasta.tar.gz

Call:

braker.pl --genome=genome.fa --prot_seq=proteins.fasta --softmasking --cores=8 --species=Chrysaora_chesapeakei --epmode

Flanking region size: 1246
Number of generated training genes: 1712 nonredundant genes, 300 set aside for testing, 300 set aside for optimize_augustus.pl time consuming steps

Initial accuracy on test set: 0.27872
Accuracy on test set after optimize_augustus.pl: 0.310626666666667

With this parameter set, BRAKER predicted 126426 genes (Augustus with hints in the target genome). To our knowledge, this might be the first gene set for this species. The gene set is available upon request from katharina.hoff@uni-greifswald.de. However, the parameter set and the gene number indicate that this parameter set does not work very well for this species!

This parameter set was generated in a fully automated fashion for supporing the POMPU project (http://www.pompu-project.de/).

If using this parameter set, please cite:

Bruna, T., Hoff, K., Stanke, M., Lomsadze, A., & Borodovsky, M. (2020). BRAKER2: Automatic Eukaryotic Genome Annotation with GeneMark-EP+ and AUGUSTUS Supported by a Protein Database. bioRxiv. https://www.biorxiv.org/content/10.1101/2020.08.10.245134v1.abstract

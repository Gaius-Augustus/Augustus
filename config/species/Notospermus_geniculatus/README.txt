# Katharina J. Hoff
# November 30th 2020
# Species specific AUGUSTUS parameters for Notospermus geniculatus

Augustus parameter set name: Notospermus_geniculatus

Lineage: Eukaryota; Opisthokonta; Metazoa; Eumetazoa; Bilateria; Protostomia; Spiralia; Lophotrochozoa; Nemertea; Pilidiophora; Heteronemertea; Lineidae; Notospermus

This parameter set was created with BRAKER on the basis of genome assembly https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/633/025/GCA_002633025.1_ASM263302v1/GCA_002633025.1_ASM263302v1_genomic.fna.gz and OrthoDb proteins https://v100.orthodb.org/download/odb10_metazoa_fasta.tar.gz

Call:

braker.pl --genome=genome.fa --prot_seq=proteins.fasta --softmasking --cores=8 --species=Notospermus_geniculatus --epmode

Flanking region size: 1981
Number of generated training genes: 2660 nonredundant genes, 300 set aside for testing, 300 set aside for optimize_augustus.pl time consuming steps

Initial accuracy on test set: 0.811066666666667
Accuracy on test set after optimize_augustus.pl: 0.843

With this parameter set, BRAKER predicted 63141 genes (Augustus with hints in the target genome). We assume that many of these genes are falsely split. To our knowledge, this might be the first gene set for this genome assembly. The gene set is available upon request from katharina.hoff@uni-greifswald.de.

This parameter set was generated in a fully automated fashion for supporing the POMPU project (http://www.pompu-project.de/).

If using this parameter set, please cite:

Bruna, T., Hoff, K., Stanke, M., Lomsadze, A., & Borodovsky, M. (2020). BRAKER2: Automatic Eukaryotic Genome Annotation with GeneMark-EP+ and AUGUSTUS Supported by a Protein Database. bioRxiv. https://www.biorxiv.org/content/10.1101/2020.08.10.245134v1.abstract

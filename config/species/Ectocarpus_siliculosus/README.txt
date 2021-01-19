# Katharina J. Hoff
# November 29th 2020
# Species specific AUGUSTUS parameters for Ectocarpus siliculosus

Augustus parameter set name: Ectocarpus_siliculosus

Lineage: Eukaryota; Sar; Stramenopiles; Ochrophyta; PX clade; Phaeophyceae; Ectocarpales; Ectocarpaceae; Ectocarpus

This parameter set was created with BRAKER on the basis of genome assembly https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/310/025/GCA_000310025.1_ASM31002v1/GCA_000310025.1_ASM31002v1_genomic.fna.gz and OrthoDb proteins https://v100.orthodb.org/download/odb10_metazoa_fasta.tar.gz

Call:

braker.pl --genome=genome.fa --prot_seq=proteins.fasta --softmasking --cores=8 --species=Ectocarpus_siliculosus --epmode

Flanking region size: 2936
Number of generated training genes: 3546 nonredundant genes, 300 set aside for testing, 300 set aside for optimize_augustus.pl time consuming steps

Initial accuracy on test set: 0.7918
Accuracy on test set after optimize_augustus.pl: 0.803066666666667

With this parameter set, BRAKER predicted 18129 genes (Augustus with hints in the target genome). For comparison, the EMBL annotation contains 16269 genes.

This parameter set was generated in a fully automated fashion for supporing the POMPU project (http://www.pompu-project.de/).

If using this parameter set, please cite:

Bruna, T., Hoff, K., Stanke, M., Lomsadze, A., & Borodovsky, M. (2020). BRAKER2: Automatic Eukaryotic Genome Annotation with GeneMark-EP+ and AUGUSTUS Supported by a Protein Database. bioRxiv. https://www.biorxiv.org/content/10.1101/2020.08.10.245134v1.abstract

# Katharina J. Hoff
# November 25th 2020
# Species specific AUGUSTUS parameters for Ostreococcus tauri

Augustus parameter set name: Ostreococcus_tauri

Lineage: Eukaryota; Viridiplantae; Chlorophyta; Mamiellophyceae; Mamiellales; Bathycoccaceae; Ostreococcus

This parameter set was created with BRAKER on the basis of genome assembly https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/214/015/GCF_000214015.3_version_140606/GCF_000214015.3_version_140606_genomic.fna.gz and OrthoDb proteins https://v100.orthodb.org/download/odb10_plants_fasta.tar.gz

Call:

braker.pl --genome=genome.fa --prot_seq=proteins.fasta --softmasking --cores=8 --species=Ostreococcus_tauri --epmode

Flanking region size: 755
Number of generated training genes: 958 nonredundant genes, 200 set aside for testing, 200 set aside for optimize_augustus.pl time consuming steps

Initial accuracy on test set: 0.839733333333333
Accuracy on test set after optimize_augustus.pl: 0.849666666666667

With this parameter set, BRAKER predicted 9815 genes (Augustus with hints in the target genome). For comparison, the RefSeq annotation contains 7766 genes.

This parameter set was generated in a fully automated fashion for supporing the POMPU project (http://www.pompu-project.de/).

If using this parameter set, please cite:

Bruna, T., Hoff, K., Stanke, M., Lomsadze, A., & Borodovsky, M. (2020). BRAKER2: Automatic Eukaryotic Genome Annotation with GeneMark-EP+ and AUGUSTUS Supported by a Protein Database. bioRxiv. https://www.biorxiv.org/content/10.1101/2020.08.10.245134v1.abstract

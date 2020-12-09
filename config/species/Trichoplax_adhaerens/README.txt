# Katharina J. Hoff
# November 29th 2020
# Species specific AUGUSTUS parameters for Trichoplax adhaerens

Augustus parameter set name: Trichoplax_adhaerens

Lineage: Eukaryota; Opisthokonta; Metazoa; Eumetazoa; Placozoa; Trichoplacidae; Trichoplax

This parameter set was created with BRAKER on the basis of genome assembly https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/150/275/GCF_000150275.1_v1.0/GCF_000150275.1_v1.0_genomic.fna.gz and OrthoDb proteins https://v100.orthodb.org/download/odb10_metazoa_fasta.tar.gz

Call:

braker.pl --genome=genome.fa --prot_seq=proteins.fasta --softmasking --cores=8 --species=Trichoplax_adherens --epmode

Flanking region size: 1927
Number of generated training genes: 3326 nonredundant genes, 300 set aside for testing, 300 set aside for optimize_augustus.pl time consuming steps

Initial accuracy on test set: 0.8682
Accuracy on test set after optimize_augustus.pl: 0.884133333333333

With this parameter set, BRAKER predicted 14293 genes (Augustus with hints in the target genome). For comparison, the RefSeq annotation contains 11520 genes.

This parameter set was generated in a fully automated fashion for supporing the POMPU project (http://www.pompu-project.de/).

If using this parameter set, please cite:

Bruna, T., Hoff, K., Stanke, M., Lomsadze, A., & Borodovsky, M. (2020). BRAKER2: Automatic Eukaryotic Genome Annotation with GeneMark-EP+ and AUGUSTUS Supported by a Protein Database. bioRxiv. https://www.biorxiv.org/content/10.1101/2020.08.10.245134v1.abstract

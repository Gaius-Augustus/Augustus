# Katharina J. Hoff
# November 27th 2020
# Species specific AUGUSTUS parameters for Vitrella brassicaformis

Augustus parameter set name: Vitrella_brassicaformis

Lineage: Eukaryota; Sar; Alveolata; Colpodellida; Vitrellaceae; Vitrella

This parameter set was created with BRAKER on the basis of genome assembly https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/179/505/GCA_001179505.1_Vbrassicaformis/GCA_001179505.1_Vbrassicaformis_genomic.fna.gz and the annotated proteins https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/179/505/GCA_001179505.1_Vbrassicaformis/GCA_001179505.1_Vbrassicaformis_protein.faa.gz

Call:

braker.pl --species=Vitrella_brassicaformis --genome=genome.fa --prot_seq=proteins.fa --prg=gth --trainFromGth --softmasking --cores=2

Flanking region size: 1140
Number of generated training genes: 7847 nonredundant genes, 300 set aside for testing, 300 set aside for optimize_augustus.pl time consuming steps

Initial accuracy on test set: 0.7878
Accuracy on test set after optimize_augustus.pl: 0.787

With this parameter set, BRAKER predicted 24148 genes (Augustus with hints in the target genome). For comparison, the EMBL annotation contains 23034 genes.

This parameter set was generated in a fully automated fashion for supporing the POMPU project (http://www.pompu-project.de/).

If using this parameter set, please cite:

Bruna, T., Hoff, K., Stanke, M., Lomsadze, A., & Borodovsky, M. (2020). BRAKER2: Automatic Eukaryotic Genome Annotation with GeneMark-EP+ and AUGUSTUS Supported by a Protein Database. bioRxiv. https://www.biorxiv.org/content/10.1101/2020.08.10.245134v1.abstract

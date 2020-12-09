# Katharina J. Hoff
# November 25th 2020
# Species specific AUGUSTUS parameters for Gonapodya prolifera

Augustus parameter set name: Gonapodya_prolifera

Lineage: Eukaryota; Opisthokonta; Fungi; Fungi incertae sedis; Chytridiomycota; Chytridiomycota incertae sedis; Monoblepharidomycetes; Monoblepharidales; Gonapodyaceae; Gonapodya

This parameter set was created with BRAKER on the basis of genome assembly https://www.ncbi.nlm.nih.gov/genome/18227?genome_assembly_id=266567 and OrthoDb proteins https://v100.orthodb.org/download/odb10_fungi_fasta.tar.gz

Call:

braker.pl --genome=genome.fa --prot_seq=proteins.fasta --softmasking --cores=8 --species=Gonapodya_prolifera --epmode

Flanking region size: 1046
Number of generated training genes: 4501 nonredundant genes, 300 set aside for testing, 300 set aside for optimize_augustus.pl time consuming steps

Initial accuracy on test set: 0.825266666666667
Accuracy on test set after optimize_augustus.pl: 0.8356

With this parameter set, BRAKER predicted 17500 genes (Augustus with hints in the target genome). For comparison, the GenBank annotation contains 13831 genes.

This parameter set was generated in a fully automated fashion for supporing the POMPU project (http://www.pompu-project.de/).

If using this parameter set, please cite:

Bruna, T., Hoff, K., Stanke, M., Lomsadze, A., & Borodovsky, M. (2020). BRAKER2: Automatic Eukaryotic Genome Annotation with GeneMark-EP+ and AUGUSTUS Supported by a Protein Database. bioRxiv. https://www.biorxiv.org/content/10.1101/2020.08.10.245134v1.abstract

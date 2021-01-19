# Katharina J. Hoff
# December 2nd 2020
# Species specific AUGUSTUS parameters for Nemopilema nomurai

Augustus parameter set name: Nemopilema_nomurai

Lineage: Eukaryota; Metazoa; Cnidaria; Scyphozoa; Rhizostomeae; Rhizostomatidae; Nemopilema

This parameter set was created with BRAKER on the basis of genome assembly https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/864/495/GCA_003864495.1_NemNom1.0/GCA_003864495.1_NemNom1.0_genomic.fna.gz and OrthoDb proteins https://v100.orthodb.org/download/odb10_metazoa_fasta.tar.gz

Call:

braker.pl --genome=genome.fa --prot_seq=proteins.fasta --softmasking --cores=8 --species=Nemopilema_nomurai --epmode

Flanking region size: 2432
Number of generated training genes: 2724 nonredundant genes, 300 set aside for testing, 300 set aside for optimize_augustus.pl time consuming steps

Initial accuracy on test set: 0.811266666666667
Accuracy on test set after optimize_augustus.pl: 0.828666666666667

With this parameter set, BRAKER predicted 21756 genes (Augustus with hints in the target genome). To our knowledge, this might be the first gene set for this species. The gene set is available upon request from katharina.hoff@uni-greifswald.de

This parameter set was generated in a fully automated fashion for supporing the POMPU project (http://www.pompu-project.de/).

If using this parameter set, please cite:

Bruna, T., Hoff, K., Stanke, M., Lomsadze, A., & Borodovsky, M. (2020). BRAKER2: Automatic Eukaryotic Genome Annotation with GeneMark-EP+ and AUGUSTUS Supported by a Protein Database. bioRxiv. https://www.biorxiv.org/content/10.1101/2020.08.10.245134v1.abstract

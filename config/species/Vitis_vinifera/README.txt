# Katharina J. Hoff
# January 25th 2022
# Species specific AUGUSTUS parameters for Vitis vinifera

Augustus parameter set name: Vitis_vinifera

Lineage: Eukaryota; Viridiplantae; Streptophyta; Streptophytina; Embryophyta; Tracheophyta; Euphyllophyta; Spermatophyta; Magnoliopsida; Mesangiospermae; eudicotyledons; Gunneridae; Pentapetalae; rosids; rosids incertae sedis; Vitales; Vitaceae; Viteae; Vitis

This parameter set was created with BRAKER on the basis of genome https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/003/745/GCF_000003745.3_12X/GCF_000003745.3_12X_genomic.fna.gz and OrthoDb proteins https://v100.orthodb.org/download/odb10_plants_fasta.tar.gz

Call:

braker.pl --genome=genome.fa --prot_seq=proteins.fasta --softmasking --cores=8 --species=Vitis_vinivera --epmode

Flanking region size: 1337
Number of generated training genes: 7243 nonredundant genes, 300 set aside for testing, 300 set aside for optimize_augustus.pl time consuming steps

Initial accuracy on test set: 0.803066666666667
Accuracy on test set after optimize_augustus.pl: 0.8334

Note: we have not generated an annotation for Vitis vinivera on the basis of this parameter set, yet.

This parameter set was generated in a fully automated fashion for supporing the POMPU project (http://www.pompu-project.de/).

If using this parameter set, please cite:

Bruna, T., Hoff, K.J., Stanke, M., Lomsadze, A., & Borodovsky, M. (2021). BRAKER2: Automatic Eukaryotic Genome Annotation with GeneMark-EP+ and AUGUSTUS Supported by a Protein Database. NAR genomics and bioinformatics, lqaa108.

# Katharina J. Hoff
# August 26th 2021
# Species specific AUGUSTUS parameters for Ricinus communis

Augustus parameter set name: Ricinus_communis

Lineage: Eukaryota; Viridiplantae; Streptophyta; Streptophytina; Embryophyta; Tracheophyta; Euphyllophyta; Spermatophyta; Magnoliopsida; Mesangiospermae; eudicotyledons; Gunneridae; Pentapetalae; rosids; fabids; Malpighiales; Euphorbiaceae; Acalyphoideae; Acalypheae; Ricinus

This parameter set was created with BRAKER on the basis of genome assembly https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/151/685/GCF_000151685.1_JCVI_RCG_1.1/GCF_000151685.1_JCVI_RCG_1.1_genomic.fna.gz and OrthoDb proteins https://v100.orthodb.org/download/odb10_plants_fasta.tar.gz

Call:

braker.pl --genome=genome.fa --prot_seq=proteins.fasta --softmasking --cores=8 --species=Ricinus_communis --epmode

Flanking region size: 1127
Number of generated training genes: 8000 nonredundant genes, 300 set aside for testing, 300 set aside for optimize_augustus.pl time consuming steps

Initial accuracy on test set: 0.8158
Accuracy on test set after optimize_augustus.pl: 0.840466666666667

With this parameter set, BRAKER predicted 53951 transcripts (Augustus with hints in the target genome). For comparison, the RefSeq annotation contains 28583 transcripts.

This parameter set was generated in a fully automated fashion for supporing the POMPU project (http://www.pompu-project.de/). This parameter set might heavily overpredict!

If using this parameter set, please cite:

Bruna, T., Hoff, K.J., Stanke, M., Lomsadze, A., & Borodovsky, M. (2021). BRAKER2: Automatic Eukaryotic Genome Annotation with GeneMark-EP+ and AUGUSTUS Supported by a Protein Database. NAR genomics and bioinformatics, lqaa108.

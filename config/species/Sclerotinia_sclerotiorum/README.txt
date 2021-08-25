# Katharina J. Hoff
# August 25th 2021
# Species specific AUGUSTUS parameters for Sclerotinia sclerotiorum

Augustus parameter set name: Sclerotinia_sclerotiorum

Lineage: Eukaryota; Opisthokonta; Fungi; Dikarya; Ascomycota; saccharomyceta; Pezizomycotina; leotiomyceta; sordariomyceta; Leotiomycetes; Helotiales; Sclerotiniaceae; Sclerotinia

This parameter set was created with BRAKER on the basis of genome assembly https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/945/GCF_000146945.2_ASM14694v2/GCF_000146945.2_ASM14694v2_genomic.fna.gz and OrthoDb proteins https://v100.orthodb.org/download/odb10_fungi_fasta.tar.gz

Call:

braker.pl --genome=genome.fa --prot_seq=proteins.fasta --softmasking --cores=8 --species=Sclerotinia_sclerotiorum --epmode

Flanking region size: 814
Number of generated training genes: 4766 nonredundant genes, 200 set aside for testing, 200 set aside for optimize_augustus.pl time consuming steps

Initial accuracy on test set: 0.9038
Accuracy on test set after optimize_augustus.pl: 0.912666666666667

With this parameter set, BRAKER predicted 11432 transcripts (Augustus with hints in the target genome). For comparison, the RefSeq annotation contains 14490 transcripts.

This parameter set was generated in a fully automated fashion for supporing the POMPU project (http://www.pompu-project.de/).

If using this parameter set, please cite:

Bruna, T., Hoff, K.J., Stanke, M., Lomsadze, A., & Borodovsky, M. (2021). BRAKER2: Automatic Eukaryotic Genome Annotation with GeneMark-EP+ and AUGUSTUS Supported by a Protein Database. NAR genomics and bioinformatics, lqaa108.

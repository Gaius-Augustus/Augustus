# Katharina J. Hoff
# November 25th 2020
# Species specific AUGUSTUS parameters for Sphaceloma murrayae

Augustus parameter set name: Sphaceloma_murrayae

Lineage: Eukaryota; Opisthokonta; Fungi; Dikarya; Ascomycota; saccharomyceta; Pezizomycotina; leotiomyceta; dothideomyceta; Dothideomycetes; Dothideomycetidae; Myriangiales; Elsinoaceae; Sphaceloma

This parameter set was created with BRAKER on the basis of genome assembly https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/895/985/GCA_002895985.1_ASM289598v1/GCA_002895985.1_ASM289598v1_genomic.fna.gz and OrthoDb proteins https://v100.orthodb.org/download/odb10_fungi_fasta.tar.gz

Call:

braker.pl --genome=genome.fa --prot_seq=proteins.fasta --softmasking --cores=8 --species=Sphaceloma_murrayae --epmode

Flanking region size: 829
Number of generated training genes: 3274 nonredundant genes, 300 set aside for testing, 300 set aside for optimize_augustus.pl time consuming steps

Initial accuracy on test set: 0.849266666666667
Accuracy on test set after optimize_augustus.pl: 0.852733333333333

This parameter set was generated in a fully automated fashion for supporing the POMPU project (http://www.pompu-project.de/).

If using this parameter set, please cite:

Bruna, T., Hoff, K., Stanke, M., Lomsadze, A., & Borodovsky, M. (2020). BRAKER2: Automatic Eukaryotic Genome Annotation with GeneMark-EP+ and AUGUSTUS Supported by a Protein Database. bioRxiv. https://www.biorxiv.org/content/10.1101/2020.08.10.245134v1.abstract

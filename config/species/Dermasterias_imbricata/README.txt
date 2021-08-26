# Katharina J. Hoff
# August 26th 2021
# Species specific AUGUSTUS parameters for Dermasterias imbricata

Augustus parameter set name: Dermasterias_imbricata

Lineage: Eukaryota; Opisthokonta; Metazoa; Eumetazoa; Bilateria; Deuterostomia; Echinodermata; Eleutherozoa; Asterozoa; Asteroidea; Valvatacea; Valvatida; Poraniidae; Dermasterias

This parameter set was created in a fully automated fashion using BRAKER1 on a confidential genome assembly & private RNA-Seq data.

Call:

braker.pl --genome=genome.fa --bam=rnaseq.bam --softmasking --cores=20 --species=Dermasterias_imbricata

Flanking region size: 4213
Number of generated training genes: 4474 nonredundant genes, 300 set aside for testing, 300 set aside for optimize_augustus.pl time consuming steps

Initial accuracy on test set: 0.864066666666667
Accuracy on test set after optimize_augustus.pl: 0.8778

With this parameter set, BRAKER predicted 27375 transcripts (Augustus with hints in the target genome).

BUSCO scores with metazoa_odb10: C:95.2%[S:81.2%,D:14.0%],F:3.5%,M:1.3%,n:954

To our knowledge, this might be the first annotation. Annotation is available upon request from katharina.hoff@uni-greifswald.de in case the consortium agrees to the individual request.

The consortium for this project consists of:

Melissa DeBiasse, Katharina J. Hoff, Lars Gabriel, Lauren Schiebelhut & Michael Dawson.

If using this parameter set, please cite:

Hoff, K.J., S. Lange, Lomsadze, A., Borodovsky, M. & Stanke M. (2016). BRAKER1: unsupervised RNA-Seq-based genome annotation with GeneMark-ET and AUGUSTUS. Bioinformatics 32 (5): 767-769.

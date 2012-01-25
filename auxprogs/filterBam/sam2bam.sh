# Script to make conversion SAM->BAM with samtools,
# assuming "samtools" is in Tonatiuh's drive

SAMTOOLS=/home/tonatiuh/samtools;
$SAMTOOLS/samtools view -bS example_1.sam -o example_1.bam
# Script to make conversion SAM->BAM with samtools,
# assuming "samtools" is in Tonatiuh's drive

SAMTOOLS=/home/tonatiuh/samtools;
SAM=$1;
BAM=`echo $SAM|sed 's/.sam/.bam/'` 

$SAMTOOLS/samtools view -bS $SAM -o $BAM
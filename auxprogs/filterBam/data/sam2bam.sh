# Script to make conversion SAM->BAM with samtools.
#
# NOTE: 
# Modify the variable SAMTOOLS to point to the
# directory where SAMTOOLS has been installed. 
#
# Created: 25-January-2012

SAMTOOLS=/home/tonatiuh/samtools;
SAM=$1;
BAM=`echo $SAM|sed 's/.sam/.bam/'` 

echo "----------------------------------------------------";
echo "Converting: $SAM into $BAM with samtools";
echo "Check samtools documentation for futher reference";
echo "----------------------------------------------------";

$SAMTOOLS/samtools view -bS $SAM -o $BAM
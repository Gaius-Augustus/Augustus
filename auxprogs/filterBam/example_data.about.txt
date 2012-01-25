# Creation of test file for filterBam demo.
#
# NOTE:
# Unless implied by an environment variable, the working path 
# assumed is the one containing the "filterBam" source code
#
# Created: 24-January-2012
# Last modified: 24-January-2012
#

########################################################
# First merge files with bamtools.
#
# The source files are files that were previously subsampled with the  
# "bamtools random" utility
#
# NOTE:
# The name of each environment variables refers to files 
# with alignments that have a feature that will be filtered 
# out by "filterBam". 
# E.g. "$COVERAGE/Example.coverage.bam" is a file with alignments
# that will be filtered out by the deafault value of "coverage" 
# the coverage criterion of filterBam. 
#
########################################################

COVERAGE=/home/tonatiuh/web-home/filterBam/bamtoolsCode/filterExample/OUT/randomSamples;
PERCID=/home/tonatiuh/web-home/filterBam/bamtoolsCode/filterExample/OUT/randomSamples
BEST=/home/tonatiuh/web-home/alignments/ExampleBam/filteredFiles/randomSamples
REFID=/home/tonatiuh/web-home/alignments/ExampleBam/filteredFiles/randomSamples

bamtools merge -in $COVERAGE/Example.coverage.samples.bam -in $PERCID/Example.percId.samples.bam -in $BEST/Example.best.bam -in $REFID/Example.refID.bam -out test.bam 

# Converting file from BAM->SAM
samtools view test.bam > test.sam


# Converting from BAM->SAM, in order to modify query names
samtools view -h test.bam > test.sam


########################################################
# In emacs
########################################################

# Change query names to "r[number]/1", so they constitute "forward" queries
^chr17_\([0-9,a-z,:,_]+\) -> r\#/1

# Rename file:
test.sam -> example_1.sam

# Rename (manually) reads: 
# r22/1, r24/1, r6/1, r20/1, r10/1, r18/1, r23/1, r14/1 
# by r/2
# Note: these are the reads that will be detected as "best" or "unique"

########################################################
# Intron file 
# INTRON=/home/tonatiuh/web-home/filterBam/bamtoolsCode/filterExample/OUT;
########################################################

# Add noIntrons file, eg1.introns.sam
cd $HOME/web-home/filterBam/bamtoolsCode/filterExample/OUT
samtools view eg1.introns.bam > eg1.introns.sam # conversion without header
cat $FILTERBAMPATH/example_1.sam eg1.introns.sam > $FILTERBAMPATH/temp.sam
cd $FILTERBAMPATH
mv temp.sam example_1.sam
# Change reference name: gi|9626243|ref|NC_001416.1| -> chr17 and 
# modify sequence number


########################################################
#
# Including reads that will be passed by the filter when 
# in modality of single reads
#
########################################################


# Append into example_1.sam the first five readings of: 
# /home/tonatiuh/web-home/filterBam/bamtoolsCode/filterExample/IN/eg1.f.bam
# into the file and rename reads to respect sequence of example_1.sam

# Convert file from SAM->BAM
samtools view -bS example_1.sam -o example_1.bam


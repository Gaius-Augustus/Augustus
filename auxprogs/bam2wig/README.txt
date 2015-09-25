Installation instructions:

0. Install dependencies:
   * bcftools
   * htslib
   * samtools
   * tabix
  You find those tools at https://github.com/samtools . bam2wig has been tested with samtools 0.2.0-rc7

1. Edit Makefile and modify the value of SAMTOOLS variable
2. Type

 make

on command prompt

Examples:
A file test.bam has ben included. You can try bam2wig by trying out
the following examples,

a) ./bam2wig test.bam 
b) ./bam2wig -t "my_specified_track" -r chr3L test.s.bam 

Example (b) can only be done if an index file for test.s.bam exists. 
Do "samtools index test.s.bam" and a file "test.s.bam.bai" will be generated.

Tonatiuh Pena
17-June-2012

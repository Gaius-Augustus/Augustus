Installation instructions:

0. Install dependencies:
   * bcftools
   * htslib
   * samtools
   * tabix
  You find those tools at https://github.com/samtools . bam2wig has been tested with samtools 0.2.0-rc7

  We recommend that you install all four tools into one directory that we will further refer to as the
  TOOLDIR. For example, you can install all tools into your home directory, which will result into
  four folders in your home directory: bcftools, htslib, samtools and tabix.

1. Export TOOLDIR path

  export TOOLDIR=/path/where/four/dependencies/reside

2. Type

 make

on command prompt

Examples:
A file test.bam has been included. You can try bam2wig by trying out
the following examples,

a) ./bam2wig test.bam 
b) ./bam2wig -t "my_specified_track" -r chr3L test.s.bam 

Example (b) can only be done if an index file for test.s.bam exists. 
Do "samtools index test.s.bam" and a file "test.s.bam.bai" will be generated.

Tonatiuh Pena Centeno, Katharina Hoff
29 November 2017

Installation instructions:

0. Install dependencies:
   * bcftools
   * htslib
   * samtools
  You find those tools at https://github.com/samtools . bam2wig has been tested with samtools 0.2.0-rc7

  We recommend that you install all three tools into one directory that we will further refer to as the
  TOOLDIR. For example, you can install all tools into your home directory, which will result into
  three folders in your home directory: bcftools, htslib and samtools.

  Installation example:

  git clone https://github.com/samtools/htslib.git
  cd htslib
  autoheader
  autoconf
  ./configure
  make
  sudo make install
  cd ..
  git clone https://github.com/samtools/bcftools.git
  cd bcftools
  autoheader
  autoconf
  ./configure
  make
  sudo make install
  cd ..
  git clone https://github.com/samtools/samtools.git
  cd samtools
  autoheader
  autoconf -Wno-syntax
  ./configure
  make
  sudo make install
  cd ..

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

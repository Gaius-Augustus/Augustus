Installation instructions:

0. Install dependencies:
   * htslib
  You find this library at https://github.com/samtools/htslib .

  We recommend that you install htslib in one directory that we will further refer to as the
  TOOLDIR.

  Installation example:

  git clone https://github.com/samtools/htslib.git
  cd htslib
  autoheader
  autoconf
  ./configure
  make
  sudo make install

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

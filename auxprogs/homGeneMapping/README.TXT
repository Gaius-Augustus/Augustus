# 6.2.2015
# Stefanie Koenig

homGeneMapping manual
=====================

homGeneMapping takes a set of gene predictions of different genomes and a hal
alignment of the genomes and prints a summary for each gene, e.g.

- how many of its exons/introns are in agreement with genes of other genomes
- how many of its exons/introns are supported by extrinsic evidence from any of the genomes
- a list of geneids of homologous genes


     INSTALLATION
=====================

### Requirements ###
* gcc 4.6 or newer
* HAL Tools (available at https://github.com/ComparativeGenomicsToolkit/hal)
  with this dependencies:
    * package libhdf5-dev or hdf5 from http://www.hdfgroup.org (version 1.10)
    * sonLib from https://github.com/benedictpaten/sonLib
  install: https://github.com/ComparativeGenomicsToolkit/hal#installing-dependencies or see Dockerfile for an example
* Boost graph library for option --printHomologs that clusters transcripts into groups of homologs
  Install the package libboost-all-dev with your package manager
  (optional: add BOOST = false to common.mk if the Boost library is unavailable and option --printHomologs is not required)
* SQLite library for option --dbaccess that enables the retrieval hints from a database:
  Install the package libsqlite3-dev with your package manager or 
  download the SQLite source code from http://www.sqlite.org/download.html 
  (tested with  SQLite 3.8.5 ) and install as follows:

   > tar zxvf sqlite-autoconf-3080500.tar.gz
   > cd sqlite-autoconf-3080500
   > ./configure
   > sudo make
   > sudo make install

   If you encounter an "SQLite header and source version mismatch" error, try

   > ./configure --disable-dynamic-extensions --enable-static --disable-shared
  (optional: add SQLITE = false to common.mk if the SQLite library is unavailable and option --dbaccess is not required)
  
From within the auxprogs/homGeneMapping directory, type "make"


        USAGE
=====================

homGeneMapping [Options] --gtfs=gffilenames.tbl --halfile=aln.hal

ARGUMENTS:
--halfile=aln.hal             input hal file
--gtfs=gtffilenames.tbl       a text file containing the locations of the input gene files
                              and optionally the hints files (both in GTF format).
                              The file is formatted as follows:

                              name_of_genome_1  path/to/genefile/of/genome_1  path/to/hintsfile/of/genome_1
                              name_of_genome_2  path/to/genefile/of/genome_2  path/to/hintsfile/of/genome_2
                              ...
                              name_of_genome_N  path/to/genefile/of/genome_N  path/to/hintsfile/of/genome_N

OPTIONS:
--help                        print this usage info
--cpus=N                      N is the number of CPUs to use (default: 1)
--noDupes                     do not map between duplications in hal graph. (default: off)
--details                     print detailed output (default: off)
--halLiftover_exec_dir=DIR    Directory that contains the executable halLiftover from the HAL Tools package
                              If not specified it must be in $PATH environment variable.
--unmapped                    print a GTF attribute with a list of all genomes, that are not aligned to the
                              corresponding gene feature, e.g. hgm_unmapped "1,4,5"; (default; off)
--tmpdir=DIR                  a temporary file directory that stores lifted over files. (default 'tmp/' in current directory)
--outdir=DIR                  file directory that stores output gene files. (default: current directory)
--printHomologs=FILE          prints disjunct sets of homologous transcripts to FILE, e.g.
                              # 0     dana
                              # 1     dere
                              # 2     dgri
                              # 3     dmel
                              # 4     dmoj
                              # 5     dper
                              (0,jg4139.t1) (0,jg4140.t1) (1,jg7797.t1) (2,jg3247.t1) (4,jg6720.t1) (5,jg313.t1)
                              (1,jg14269.t1) (3,jg89.t1) (5,jg290.t1)
                              ...
                              Two transcripts are in the same set, if all their exons/introns are homologs and their are
                              no additional exons/introns.
                              This option requires the Boost C++ Library
--dbaccess=db                 retrieve hints from an SQLite database. In order to set up a database and populate it with hints
                              a separate tool 'load2sqlitedb' is provided. For more information, see the documentation in
                              README-cgp.txt (section 8a+b) in the Augustus package. If both a database and hint files in 'gtffilenames.tbl'
                              are specified, hints are retrieved from both sources.

examples:
homGeneMapping --noDupes --halLiftOver_exec_dir=~/tools/progressiveCactus/submodules/hal/bin --gtfs=gtffilenames.tbl --halfile=msca.hal
homGeneMapping --gtfs=gtffilenames.tbl --halfile=aln.hal --outdir=outdir --cpus=4

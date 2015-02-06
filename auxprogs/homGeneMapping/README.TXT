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

### Requirements
* gcc 4.6 or newer
* HAL Tools (available at https://github.com/glennhickey/hal)


From within the src directory, type "make"



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
--halLiftover_exec_dir=DIR    Directory that contains the executable halLiftover from the HAL Tools package
                              If not specified it must be in $PATH environment variable.
--tmpdir=DIR                  a temporary file directory that stores lifted over files. (default 'tmp/' in current directory)
--outdir=DIR                  file direcory that stores output gene files. (default: current directory)

example:
homGeneMapping --noDupes --halLiftOver_exec_dir=~/tools/progressiveCactus/submodules/hal/bin --gtfs=gtffilenames.tbl --halfile=msca.hal

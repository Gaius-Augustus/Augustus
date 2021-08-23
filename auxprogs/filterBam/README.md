
# filterBam: a tool for cleaning alignment files in BAM format

[INTRODUCTION](#introduction)  
[DEPENDENCIES](#dependencies)
[INSTALLATION](#installation)  
[RUNNING](running)

# INTRODUCTION

filterBAM is a C++ routine for filtering BAM alignment files; it is based on [filterPSL.pl](../../scripts/filterPSL.pl), a Perl script written
for the cleaning of PSL alignment files.


filterBam only accepts as input BAM files that have been previously sorted by 'query name'.

# DEPENDENCIES

1. Bamtools' API is required for compilation. Bamtools allows the sorting, sampling, filtering, among other things, of BAM files. 

  * On Ubunutu or Debian

        apt-get install bamtools libbamtools-dev

     See [docs/INSTALL.md](../../docs/INSTALL.md#Bamtools) (Bamtools section) for further details, especially for manual 
     installation from sources or installation without root rights.

2. Although not required by filterBam, it may also be convenient to have a copy of the software Samtools.
   See [docs/INSTALL.md](../../docs/INSTALL.md#SAMtools) (SAMtools section) for further details.

   NOTE: Some examples on how to use Samtools and Bamtools utility kits are provided in [data/example_data.about.txt](data/example_data.about.txt).

# INSTALLATION

  1. Download and extract the latest version of filterBam from the AUGUSTUS repository.

  2. Compile filterBam by typing

        make

  3. A binary file 'filterBam' should be stored in the folder 'bin'

# RUNNING

  Get the 'help' menu by typing

        ./filterBam --help

  Some toy data sets are stored in the folder 'data' in case you want to see how the filter works.


[![Build Status](https://travis-ci.org/Gaius-Augustus/Augustus.svg?branch=master)](https://travis-ci.org/Gaius-Augustus/Augustus)
[![Test examples](https://github.com/Gaius-Augustus/Augustus/workflows/Test%20examples/badge.svg)](https://github.com/Gaius-Augustus/Augustus/actions?query=workflow%3A"Test+examples")

# Gene Prediction with AUGUSTUS

[INTRODUCTION](#introduction)  
[INSTALLATION](#installation)  
[RUNNING AUGUSTUS](docs/RUNNING-AUGUSTUS.md)  
[WEB-SERVER](#web-server)  
[COMPARATIVE GENE PREDICTION](docs/README-cgp.md)  
[AUTHORS AND CONTACT](docs/CONTACT.md)  
[REFERENCES](#references-and-documentation)  
[LICENSES](#licenses)  

# INTRODUCTION

AUGUSTUS is a program to find genes and their structures in one or more genomes. [More ...](docs/ABOUT.md)

# INSTALLATION

## Windows
Windows users can use the Windows Subsystem for Linux (WSL) to install AUGUSTUS exactly as described below for Linux. How to set up the WSL for AUGUSTUS is described [here](docs/AUGUSTUS-ON-WINDOWS.md).

## Ubuntu 18.04 or later
This currently installs only a single-genome version without comparative gene prediction capability:
```
sudo apt install augustus augustus-data augustus-doc
```

## Clone from GitHUB

To obtain the most recent complete version, first, clone the repository

```
git clone https://github.com/Gaius-Augustus/Augustus.git
```
or, alternatively, [download](http://bioinf.uni-greifswald.de/augustus/binaries/augustus.current.tar.gz) and unpack the AUGUSTUS source package with
```
tar -xzf augustus.current.tar.gz
```

## Docker

After retrieving Augustus change into the main directory containing the 
[Dockerfile](Dockerfile) and issue

```
docker build -t augustus .
```

## Install dependencies

The following dependencies are required for AUGUSTUS:
- For gzip compressed input:
 (set ZIPINPUT = false in [common.mk](common.mk) if this feature is not required or the required libraries are not available)
  - libboost-iostreams-dev
  - zlib1g-dev
- For [comparative AUGUSTUS](docs/README-cgp.md) (multi-species, CGP):
  (set COMPGENEPRED = false in [common.mk](common.mk) if the libraries required by the CGP version are not available. Augustus can then only be run in single-genome mode, which is what most users need.)
  - libgsl-dev
  - libboost-all-dev
  - libsuitesparse-dev
  - liblpsolve55-dev
  - libsqlite3-dev (add SQLITE = false to [common.mk](common.mk) if this feature is not required or the required library is not available)
  - libmysql++-dev (add MYSQL = false to [common.mk](common.mk) if this feature is not required or the required library is not available)
- For compiling bam2hints and filterBam:
  - libbamtools-dev
- For compiling utrrnaseq:
  - libboost-all-dev (version must be >Boost_1_49_0)
- For compiling bam2wig:
  - Follow [these instructions](./auxprogs/bam2wig/README.md). Note that it shouldn't be a problem to compile AUGUSTUS without bam2wig. In practice, you can simply use `bamToWig.py` to accomplish the same task.
- For compiling homgenemapping
  (set BOOST = FALSE in [./auxprogs/homgenemapping/src/Makefile](./auxprogs/homgenemapping/src/Makefile) if the option --printHomologs is not required or the required libraries are not available)
  - libboost-all-dev

## Compile from sources

Once all dependencies are available, you can compile AUGUSTUS using make.

```
make
```

After compilation has finished, the command bin/augustus should be executable and print a usage message.

## Install locally

As a normal user, add the directory of the executables to the PATH environment variable, for example:

```
export PATH=$PATH:~/augustus/bin:~/augustus/scripts
```

## Install globally

You can install AUGUSTUS globally, if you have root privileges, for example: 

```
sudo make install
```

Alternatively, you can exectue similar commands to those in the "install" section of the top-level Makefile to customize the global installation. 

## Optional: set environment variable AUGUSTUS_CONFIG_PATH

If the environment variable AUGUSTUS_CONFIG_PATH is set, augustus and etraining will look there for the config directory that contains the configuration and parameter files, e.g. '~/augustus/config'. You may want to add this line to a startup script (like ~/.bashrc).

```
export AUGUSTUS_CONFIG_PATH=/my_path_to_AUGUSTUS/augustus/config/
```

If this environment variable is not set, then the programs will look in the path ../config relative to the directory in which the executable lies. As a third alternative, you can specify this directory on the command line when you run augustus:
--AUGUSTUS_CONFIG_PATH=/my_path_to_AUGUSTUS/augustus/config/

# WEB-SERVER

AUGUSTUS can also be run through a web-interface at http://bioinf.uni-greifswald.de/augustus/ and a web service at http://bioinf.uni-greifswald.de/webaugustus/index.gsp.

# Scripts

Many scripts require Perl.

Running the python3 script bamToWig.py has the following software dependencies:
  - Python3
  - twoBitInfo and faToTwoBit from http://hgdownload.soe.ucsc.edu/admin/exe . bamToWig.py will automatically download these tools to the working directory during execution	if they	are not	in your	$PATH.
  - samtools (available e.g. at https://github.com/samtools/samtools or via package managers)

# REFERENCES AND DOCUMENTATION

Mario Stanke, Mark Diekhans, Robert Baertsch, David Haussler (2008).
[Using native and syntenically mapped cDNA alignments to improve de novo gene finding](https://academic.oup.com/bioinformatics/article/24/5/637/202844). Bioinformatics, 24(5), pages 637–644, doi: 10.1093/bioinformatics/btn013

For further references see [docs/REFERENCES.md](docs/REFERENCES.md)

[3 book chapters with command line walkthroughs](https://math-inf.uni-greifswald.de/en/department/about-us/employees/prof-dr-mario-stanke-english/publications/#c302071)

# LICENSES

All source code, i.e.
  - the AUGUSTUS source code (src/*.cc, include/*.hh)
  - the scripts (scripts/*)
  - the auxiliary programs (auxprogs/)
  - the tree-parser (src/scanner, src/parser)
  - the unit tests (src/unittests)
  
is under the [Artistic License](src/LICENSE.TXT).

[![Build Status](https://travis-ci.org/Gaius-Augustus/Augustus.svg?branch=master)](https://travis-ci.org/Gaius-Augustus/Augustus)

# Gene Prediction with AUGUSTUS

[INTRODUCTION](#introduction)  
[INSTALLATION](#installation)  
[RUNNING AUGUSTUS](docs/RUNNING-AUGUSTUS.md)  
[WEB-SERVER](#web-server)  
[AUTHORS AND CONTACT](docs/CONTACT.md)  
[REFERENCES](#references)  
[LICENCES](#licenses)  

# INTRODUCTION

AUGUSTUS is a program to find genes and their structures in one or more genomes. [More ...](docs/ABOUT.md)

# INSTALLATION

## Ubuntu 18.04 or later
This currently installs only a single-genome version without comparative gene prediction capability:
```
sudo apt install augustus augustus-data augustus-doc
```

## Docker
Download [Dockerfile](Dockerfile) into some new directory, change into this directory and issue

```
docker build -t augustus https://github.com/Gaius-Augustus/Augustus.git
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

# Install dependencies

The following dependencies may be required for AUGUSTUS:
- Optional libraries for gzip compressed input (uncomment ZIPINPUT = TRUE in common.mk):
  - libboost-iostreams-dev
  - zlib1g-dev
- Optional for comparative (multi-species, CGP) AUGUSTUS with SQLITE (uncomment COMPGENEPRED = true and SQLITE = true in common.mk):
  - libgsl-dev
  - libboost-graph-dev
  - libsuitesparse-dev
  - liblpsolve55-dev
  - libsqlite3-dev
- Optional for comparative (mutli-species, CGP) AUGUSTUS with MySQL (uncomment COMPGENEPRED = true and MYSQL = true in common.mk):
  - libgsl-dev
  - libboost-graph-dev
  - libsuitesparse-dev
  - liblpsolve55-dev
  - libmysql++-dev
- For compiling bam2hints and filterBam:
  - libbamtools-dev
- For compiling bam2wig:
  - Follow [these instructions](./auxprogs/bam2wig/README.txt)
- For compiling utrrnaseq:
  - libboost-all-dev (version must be >Boost_1_49_0)

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

# REFERENCES

Mario Stanke, Mark Diekhans, Robert Baertsch, David Haussler (2008).
[Using native and syntenically mapped cDNA alignments to improve de novo gene finding](https://academic.oup.com/bioinformatics/article/24/5/637/202844). Bioinformatics, 24(5), pages 637–644, doi: 10.1093/bioinformatics/btn013

For further references see [docs/REFERENCES.md](docs/REFERENCES.md)

# LICENSES

All source code, i.e.
  - the AUGUSTUS source code (src/*.cc, include/*.hh)
  - the scripts (scripts/*.pl)
  - the auxiliary programs (auxprogs/)
  - the tree-parser (src/scanner,src/parser)
  
is under the [Artistic License](src/LICENSE.TXT).

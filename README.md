[![Build and test](https://github.com/Gaius-Augustus/Augustus/workflows/Build%20and%20test/badge.svg)](https://github.com/Gaius-Augustus/Augustus/actions?query=workflow%3A"Build+and+test")

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

## Ubuntu 18.04, Debian 9 or later
Until Ubuntu 21.04 and Debian 11 only as single-genome version, since then with capability for comparative gene prediction.

    sudo apt install augustus augustus-data augustus-doc

## Docker

Create a docker image from [Dockerfile](Dockerfile) using:

    git clone https://github.com/Gaius-Augustus/Augustus.git
    cd Augustus
    docker build -t augustus .

## Singularity

Create a Singularity Image File from the [Singularity Definition File](Singularity.def) using

    git clone https://github.com/Gaius-Augustus/Augustus.git
    cd Augustus
    singularity build augustus.sif Singularity.def

## Building AUGUSTUS from source

See [INSTALL.md](docs/INSTALL.md) for details.

Download source code from github and compile:

    git clone https://github.com/Gaius-Augustus/Augustus.git
    cd Augustus
    make augustus

After compilation has finished, the command bin/augustus should be executable and print a usage message.

For utilities use

    make auxprogs


### Install locally

As a normal user, add the directory of the executables to the PATH environment variable, for example:

    export PATH=~/augustus/bin:~/augustus/scripts:$PATH

### Install globally

You can install AUGUSTUS globally, if you have root privileges, for example:

    sudo make install

Alternatively, you can exectue similar commands to those in the "install" section of the top-level Makefile to customize the global installation.

### Optional: set environment variable AUGUSTUS_CONFIG_PATH

If the environment variable AUGUSTUS_CONFIG_PATH is set, augustus and etraining will look there for the config directory that contains the configuration and parameter files, e.g. '~/augustus/config'. You may want to add this line to a startup script (like ~/.bashrc).

    export AUGUSTUS_CONFIG_PATH=/my_path_to_AUGUSTUS/augustus/config/

If this environment variable is not set, then the programs will look in the path ../config relative to the directory in which the executable lies. As a third alternative, you can specify this directory on the command line when you run augustus:
--AUGUSTUS_CONFIG_PATH=/my_path_to_AUGUSTUS/augustus/config/

# WEB-SERVER

AUGUSTUS can also be run through a web-interface at http://bioinf.uni-greifswald.de/augustus/ and a web service at http://bioinf.uni-greifswald.de/webaugustus/.

# REFERENCES AND DOCUMENTATION

Mario Stanke, Mark Diekhans, Robert Baertsch, David Haussler (2008).
[Using native and syntenically mapped cDNA alignments to improve de novo gene finding](https://academic.oup.com/bioinformatics/article/24/5/637/202844). Bioinformatics, 24(5), pages 637–644, doi: 10.1093/bioinformatics/btn013

For further references see [docs/REFERENCES.md](docs/REFERENCES.md)

[3 book chapters with command line walkthroughs](https://math-inf.uni-greifswald.de/en/department/about-us/employees/prof-dr-mario-stanke-english/publications/#c302071)

# LICENSES

All source code, i.e.
  - the AUGUSTUS source code (`src/*.cc`, `include/*.hh`)
  - the scripts (`scripts/*`)
  - the auxiliary programs (`auxprogs/`)
  - the tree-parser (`src/scanner`, `src/parser`)
  - the unit tests (`src/unittests`)

is under the [Artistic License](src/LICENSE.TXT).

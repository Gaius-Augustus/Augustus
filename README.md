[![Build and test](https://github.com/Gaius-Augustus/Augustus/workflows/Build%20and%20test/badge.svg)](https://github.com/Gaius-Augustus/Augustus/actions?query=workflow%3A"Build+and+test") ![GitHub all releases](https://img.shields.io/github/downloads/gaius-augustus/augustus/total) 
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/augustus/README.html)
[![European Galaxy server](https://img.shields.io/badge/usegalaxy-.eu-brightgreen?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABgAAAASCAYAAABB7B6eAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAACXBIWXMAAAsTAAALEwEAmpwYAAACC2lUWHRYTUw6Y29tLmFkb2JlLnhtcAAAAAAAPHg6eG1wbWV0YSB4bWxuczp4PSJhZG9iZTpuczptZXRhLyIgeDp4bXB0az0iWE1QIENvcmUgNS40LjAiPgogICA8cmRmOlJERiB4bWxuczpyZGY9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkvMDIvMjItcmRmLXN5bnRheC1ucyMiPgogICAgICA8cmRmOkRlc2NyaXB0aW9uIHJkZjphYm91dD0iIgogICAgICAgICAgICB4bWxuczp0aWZmPSJodHRwOi8vbnMuYWRvYmUuY29tL3RpZmYvMS4wLyI+CiAgICAgICAgIDx0aWZmOlJlc29sdXRpb25Vbml0PjI8L3RpZmY6UmVzb2x1dGlvblVuaXQ+CiAgICAgICAgIDx0aWZmOkNvbXByZXNzaW9uPjE8L3RpZmY6Q29tcHJlc3Npb24+CiAgICAgICAgIDx0aWZmOk9yaWVudGF0aW9uPjE8L3RpZmY6T3JpZW50YXRpb24+CiAgICAgICAgIDx0aWZmOlBob3RvbWV0cmljSW50ZXJwcmV0YXRpb24+MjwvdGlmZjpQaG90b21ldHJpY0ludGVycHJldGF0aW9uPgogICAgICA8L3JkZjpEZXNjcmlwdGlvbj4KICAgPC9yZGY6UkRGPgo8L3g6eG1wbWV0YT4KD0UqkwAAAn9JREFUOBGlVEuLE0EQruqZiftwDz4QYT1IYM8eFkHFw/4HYX+GB3/B4l/YP+CP8OBNTwpCwFMQXAQPKtnsg5nJZpKdni6/6kzHvAYDFtRUT71f3UwAEbkLch9ogQxcBwRKMfAnM1/CBwgrbxkgPAYqlBOy1jfovlaPsEiWPROZmqmZKKzOYCJb/AbdYLso9/9B6GppBRqCrjSYYaquZq20EUKAzVpjo1FzWRDVrNay6C/HDxT92wXrAVCH3ASqq5VqEtv1WZ13Mdwf8LFyyKECNbgHHAObWhScf4Wnj9CbQpPzWYU3UFoX3qkhlG8AY2BTQt5/EA7qaEPQsgGLWied0A8VKrHAsCC1eJ6EFoUd1v6GoPOaRAtDPViUr/wPzkIFV9AaAZGtYB568VyJfijV+ZBzlVZJ3W7XHB2RESGe4opXIGzRTdjcAupOK09RA6kzr1NTrTj7V1ugM4VgPGWEw+e39CxO6JUw5XhhKihmaDacU2GiR0Ohcc4cZ+Kq3AjlEnEeRSazLs6/9b/kh4eTC+hngE3QQD7Yyclxsrf3cpxsPXn+cFdenF9aqlBXMXaDiEyfyfawBz2RqC/O9WF1ysacOpytlUSoqNrtfbS642+4D4CS9V3xb4u8P/ACI4O810efRu6KsC0QnjHJGaq4IOGUjWTo/YDZDB3xSIxcGyNlWcTucb4T3in/3IaueNrZyX0lGOrWndstOr+w21UlVFokILjJLFhPukbVY8OmwNQ3nZgNJNmKDccusSb4UIe+gtkI+9/bSLJDjqn763f5CQ5TLApmICkqwR0QnUPKZFIUnoozWcQuRbC0Km02knj0tPYx63furGs3x/iPnz83zJDVNtdP3QAAAABJRU5ErkJggg==)](https://usegalaxy.eu/root?tool_id=augustus)

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

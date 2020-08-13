[![Build Status](https://travis-ci.org/Gaius-Augustus/Augustus.svg?branch=master)](https://travis-ci.org/Gaius-Augustus/Augustus)

# Gene Prediction with AUGUSTUS

[INTRODUCTION](#introduction)  
[DATASET](#dataset)  
[BASIC](#basic)  
[CAVEAT](#caveat)  
[LICENSES](#licenses)  

# INTRODUCTION

The present version of Augustus allows to run prediction over a minimal dataset to assess how changes in the code may have affected accuracy. To use testing related function, please set TESTING = true within [../src/Makefile](../src/Makefile). If either this feature is not required or the required libraries are not available, leave/set TESTING = false (default value).

# DATASET

The present distribution comes with an unbiased dataset containing 352 genes from chromosome 1 of human genome (hg38). 

# BASIC

To make testing available, make sure TESTING is enabled. Under such condition normal prediction can be still run instead of testing, invoking Augustus command line as usual.

## Dependencies

Please make sure boost-serialization-dev is installed. If it is not the case, try and run:

```
sudo apt-get update
sudo apt-get install libboost-serialization-dev
```
In order to have accuracy returned, [EVAL](https://mblab.wustl.edu/software/download/eval-2.2.8.tar.gz) should be installed on your machine. Make sure it is the case

Also, joingenes is required component: please build it running make from within auxprogs/joingenes directory.

Running the python3 script executeTestCGP.py has the following software dependency:
  - Python3

## First use

Run the following script located within scripts sub-directory:
```
python3 executeTestCGP.py --chunks 27 30 47 54 57 80 86 101 118 --run --augustusDir=my_path_to_binaries --workingDir=my_path_to_dataSet
```
replacing my_path_to_binaries with the path to the directory named **bin**, containing comparative Augustus executable on your machine and my_path_to_dataSet with the path to data set used for testing (the dataset can be downloaded from [DATASET](http://bioinf.uni-greifswald.de/bioinf/downloads/data/aug-test/cgp12way.tgz). When done:
```
python3 executeTestCGP.py --chunks 27 30 47 54 57 80 86 101 118 --eval --evalDir=my_path_to_Eval
```
replacing my_path_to_Eval with the path to Eval script. In --eval mode, the script computes accuracy for predictions obtained over the data-set during the execution of --run mode. Results are stored in example/cgp12way/ACCURACY. 

# CAVEAT

Anyone who plans to modify executeTestCGP.py or use part of it, is recommended not to discard the following option, unless really aware of what she or he is doing:
```
--stopCodonExcludedFromCDS=true
```

# LICENSES

All source code, i.e.
  - the AUGUSTUS source code (src/*.cc, include/*.hh)
  - the scripts (scripts/*.pl)
  - the auxiliary programs (auxprogs/)
  - the tree-parser (src/scanner,src/parser)
  - include latest config files which return the best accuracy so far obtained (pending)
  
is under the [Artistic License](src/LICENSE.TXT).

[![Build Status](https://travis-ci.org/Gaius-Augustus/Augustus.svg?branch=master)](https://travis-ci.org/Gaius-Augustus/Augustus)

# Gene Prediction with AUGUSTUS in comparative gene prediction (cgp) mode
  genes are predicted simultaneously in several aligned genomes

[INTRODUCTION](#introduction)  
[INSTALLATION](#installation)  
[RUNNING AUGUSTUS IN CGP MODE](RUNNING-AUGUSTUS-IN-CGP-MODE.md)  
[BUILDING THE NEWICK PARSER FROM SCRATCH](CHANGE-TREE-PARSER.md)  
   (not needed unless you run into compiler errors related to 'parse.cc' or 'lex.cc')  
[AUTHORS AND CONTACT](CONTACT.md)  
[REFERENCES](../README.md#references-and-documentation)  
[LICENSES](../README.md#licenses)  

# INTRODUCTION

   The cgp mode is an extension to AUGUSTUS that takes an alignment of two or more genomes
   and simultaneously predicts genes in all of them.
   Beside the genomes and the alignment, a phylogenetic tree of the species is required input.
   AUGUSTUS-cgp can either be used

   - de novo, or
   - with extrinsic evidence for any subset of species  
     Such evidence includes for example already existing and trusted gene structures 
     or hints from RNA-Seq alignments.

   Both genomes and extrinsic evidence can either be read in from a flat file or 
   alternatively retrieved from a MySQL or SQLite database.  
   All three approaches are described below in more detail.

   This manual assumes that you are already familiar with AUGUSTUS
   and that you know how to use AUGUSTUS for gene prediction in a single genome.

# INSTALLATION

## Install dependencies

   See [these instructions](../README.md#install-dependencies) for a complete overview.
   
## Compile from sources

   Once all dependencies are available, you can compile AUGUSTUS using make.

    make

   In case you had previously compiled AUGUSTUS with disabled cgp mode first you 
   have to call

    make clean

   After compilation has finished, the command bin/augustus should be executable 
   and print a usage message.

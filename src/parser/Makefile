#!/bin/bash

bisonc++ -V -l grammar
cd ../scanner
flex lexer
#quick and dirty fix to some strange flex bug
sed -i 's/extern "C" int isatty (int );/extern "C" int isatty (int ) throw();/' scanner.cc

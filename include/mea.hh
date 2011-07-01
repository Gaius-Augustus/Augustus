#ifndef _MEA_HH
#define _MEA_HH

#include "gene.hh"

list<Gene>* getMEAtranscripts(Gene **sampledGeneStructures, int n, int strlength);
list<Gene>* getMEAtranscripts(list<Gene> *alltranscripts, int strlength);

Gene* listToGene(list<Gene> *genelist);
list<Gene>* geneToList(Gene *genes);

#endif

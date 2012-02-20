#ifndef _MEA7_HH
#define _MEA7_HH

#include "gene.hh"
#include "graph.hh"

void getMEAtranscripts7(list<Gene> *meaGenes, list<Gene> *alltranscripts, int strlength);
void getMeaGenelist7(list<Node*> meaPath, list<Gene> *meaGenes);


#endif

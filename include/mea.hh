#ifndef _MEA_HH
#define _MEA_HH

#include "gene.hh"
#include "graph.hh"

void getMEAtranscripts(list<Gene> *MEAtranscripts, Gene **sampledGeneStructures, int n, int strlength);
void getMEAtranscripts(list<Gene> *meaGenes, list<Gene> *alltranscripts, int strlength);

void addToList(State *st, Statename name, list<Status> *slist);
bool compareStatus(Status first, Status second);
bool compareGenes(Gene first, Gene second);
void printStatelist(list<Status> *list);
void getMeaGenelist(list<Node*> meaPath, list<Gene> *meaGenes);
void addExonToGene(Gene *gene, State *exon);
void addIntronToGene(Gene *gene, Node *predExon, Node *succExon);
StateType getIntronStateType(State *exon1, State *exon2);
void setGeneProperties(Gene *gene);
#endif

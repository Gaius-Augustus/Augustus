#ifndef _MEA_HH
#define _MEA_HH

#include "gene.hh"
#include "graph.hh"

/*
 * interface to AUGUSTUS, getMEAtranscripts() is called in findGenes() : namgene.cc
 */

void getMEAtranscripts(list<Gene> *MEAtranscripts, Gene **sampledGeneStructures, int n, const char* dna);
list<Transcript*> &getMEAtranscripts(list<Transcript*> &alltranscripts, const char* dna);
void buildStatusList(list<Transcript*> &alltranscripts, bool utr, list<Status> &stlist);
void addToList(State *st, Statename name, list<Status> *slist);
bool compareStatus(Status first, Status second);
bool compareGenes(Gene first, Gene second);
void printStatelist(list<Status> *list);
void getMeaGenelist(list<Node*> meaPath, list<Transcript*> *meaGenes);
void addExonToGene(Transcript *tx, State *exon);
void addIntronToGene(Transcript *tx, Node *predExon, Node *succExon);
void addIntronToGene(Transcript *tx, State *intr);
StateType getIntronStateType(State *exon1, State *exon2);
void setGeneProperties(Transcript *tx);
#endif

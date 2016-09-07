#include <stdlib.h>
#include <string.h>
#include <cstring>
#include <iostream>
#include <fstream>
#include <list>
#include <unordered_map>
#include <algorithm>
#include <map>
#include <sstream>

using namespace std;

map<string,int> getFileNames (string fofn);
void load(unordered_map<string,Gene*> &geneMap, string &filename, int &priority, unordered_map<string,bool> &taxaMap, Properties &properties);
void saveOverlap(list<Transcript*> &overlap, string outFileName, Properties &properties);
void saveNew(list<Transcript*> &overlap, string outFileName, Properties &properties);

string nextFreeGeneID(Properties &properties, unordered_map<string,Gene*>* addGeneMap);
string nextFreeTxID(Gene* gene, Properties &properties, unordered_map<string,Transcript*>* addTranscriptMap);
//void mergeGenes(Properties &properties, Gene* acc, Gene* don);

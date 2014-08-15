#include <stdlib.h>
#include <string.h>
#include <cstring>
#include <iostream>
#include <fstream>
#include <list>
#include <unordered_map>
#include <algorithm>
#include <unordered_map>

using namespace std;

void load(unordered_map<string,Gene> &gene_map, list<Transcript> &transcript_list, string &filename, int &priority);
void saveOverlap(list<Transcript*> &overlap, string &filename);
void save(unordered_map<string,Gene> &gene_map, list<Transcript> &transcript_list, string &filename);

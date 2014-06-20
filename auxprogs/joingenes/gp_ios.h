#include <stdlib.h>
#include <string.h>
#include <vector>
#include <cstring>
#include <iostream>
#include <vector>
#include <fstream>
#include <stack>
#include <list>
#include <unordered_map>
#include <sstream>
#include <algorithm>

using namespace std;

void load(unordered_map<string,Gene> &gene_map, list<Transcript> &transcript_list, string &filename, int &priority);
void save_overlap(list<Transcript*> &overlap, string &filename);
void save(unordered_map<string,Gene> &gene_map, list<Transcript> &transcript_list, string &filename);

/**********************************************************************
 * file:    orthoexon.hh
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  maintains orthologous exons for comparative gene prediction
 * authors: Stefanie König
 *
 *********************************************************************/

#ifndef _ORTHOEXON_HH
#define _ORTHOEXON_HH

#include "graph.hh"
#include "orthograph.hh"
#include <map>
#include <vector>
#include <string>

//forward declarations:
class ExonCandidate;

class OrthoExon{

public:
  OrthoExon();
  ~OrthoExon();
  //copy constructor
  OrthoExon(const OrthoExon& other);
  //copy with permutation of vector entries
  OrthoExon(const OrthoExon& other, const vector<size_t> &permutation);

  vector<ExonCandidate*> orthoex;
  string labelpattern;  //changes dynamically and has to be updated after every optimization step
  string getKey(const OrthoGraph &orthograph); //determines the key of an orthoex for the map labelscore
  
};

/*
 * read and write functions for orthologous exons
 */
list<OrthoExon> readOrthoExons(string filename); //read list of orthologous exons from a file
void writeOrthoExons(const list<OrthoExon> &all_orthoex);


ostream& operator<<(ostream& ostrm, const OrthoExon &ex_tuple);
istream& operator>>(istream& istrm, OrthoExon& ex_tuple);

struct Score{
  double treescore;   //stores the score of a label pattern
  int count;          //counts the number of exontuples which have that specific pattern
};

/*
 * hashfunction storing all label patterns and their Score
 * key: string over alphabet {0,1,2}^k, k = # species^k
 * 0 codes for exon in graph, but not part of the maximum weight path
 * 1 codes for exon in graph and part of the maximum weight path
 * 2 codes for exon not in graph: Status* = NULL
 * first digit in string refers to first species in PhyloTree::species, second digit refers to second species, ...
 */

namespace cache{

  extern map<string, Score> labelscore;
  /*
   * cache functions
   */
  bool inHash(string key);
  void resetCounter();
  void addToHash(string key, double score);
  double getScore(string key);
  void incrementCounter(string key);
  void printCache(list<OrthoExon> &ortho);
}

#endif

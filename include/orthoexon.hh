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
#include <map>
#include <vector>
#include <string>

using namespace std;

class OrthoExon{

public:
  OrthoExon();
  ~OrthoExon();

  vector<Status*> orthoex;
  static vector<string> species;
  static size_t getVectorPositionSpecies(string name);
  /*
   * cache functions to store score of already calculated label patterns
   */
  string getKey();
  bool inHash();
  void addToHash(double score);
  double getScore();
  void incrementCounter();

};
/*
 * read and write functions for orthologous exons
 */
list<OrthoExon> readOrthoExons(string filename); //read list of orthologous exons from a file
ostream& operator<<(ostream& ostrm, OrthoExon& ex_tuple);
istream& operator>>(istream& istrm, OrthoExon& ex_tuple);

struct Score{
  double treescore;   //stores the score of a label pattern
  int count;          //counts the number of exontuples which have that specific pattern
};
/*
 * hashfunction storing all label patterns and their Score
 * key: string over alphabet {0,1,2}^k, k = # species
 * 0 codes for exon in graph, but not part of the best gene structure
 * 1 codes for exon in graph and part of the best gene structure
 * 2 codes for exon not in graph: Status* = NULL
 * first digit in string refers to first species in vector<string> species, second digit refers to second species, ...
 */
extern map<string, Score> labelscore;


#endif

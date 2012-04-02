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

class OrthoExon{

public:
  OrthoExon(){
    orthoex.resize(OrthoExon::species.size());
  }
  ~OrthoExon(){
    for(size_t pos = 0; pos < orthoex.size(); pos++){
      delete orthoex[pos];
    }
  }
  //copy constructor
  OrthoExon(const OrthoExon& other){
    orthoex.resize(other.orthoex.size());
    for(size_t pos = 0; pos < orthoex.size(); pos++){
      if (other.orthoex[pos]){
	orthoex[pos] = new State(*other.orthoex[pos]);
      }
    }
  }

  static vector<string> species;
  vector<State*> orthoex;
  static size_t getVectorPositionSpecies(string name);
  string getKey(const OrthoGraph &orthograph) const; //determines the key of an orthoex for the map labelscore

};

/*
 * read and write functions for orthologous exons
 */
map< vector<string>, list<OrthoExon> > readOrthoExons(string filename); //read list of orthologous exons from a file
void writeOrthoExons(const map< vector<string>, list<OrthoExon> > &all_orthoex);
ostream& operator<<(ostream& ostrm, const OrthoExon &ex_tuple);
istream& operator>>(istream& istrm, OrthoExon& ex_tuple);


struct Score{
  double treescore;   //stores the score of a label pattern
  int count;          //counts the number of exontuples which have that specific pattern
};

/*
 * hashfunction storing all label patterns and their Score
 * key: string over alphabet {0,1,2}^k, k = # species^k
 * 0 codes for exon in graph, but not part of the best gene structure
 * 1 codes for exon in graph and part of the best gene structure
 * 2 codes for exon not in graph: Status* = NULL
 * first digit in string refers to first species in OrthoExon::species, second digit refers to second species, ...
 */

namespace cache{

  extern map<string, Score> labelscore;
  /*
   * cache functions
   */
  bool inHash(string key);
  void addToHash(string key, double score);
  double getScore(string key);
  void incrementCounter(string key);
}

#endif

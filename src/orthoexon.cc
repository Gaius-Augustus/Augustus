/**********************************************************************
 * file:    orthoexon.cc
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  maintains orthologous exons for comparative gene prediction
 * authors: Stefanie König
 *
 * date    |   author      |  changes
 * --------|---------------|------------------------------------------
 * 09.03.12|Stefanie König | creation of the file
 **********************************************************************/


#include "orthoexon.hh"
#include "exoncand.hh"
#include "projectio.hh"
#include "types.hh"
#include <fstream>
#include <iostream>


OrthoExon::OrthoExon(){
  orthoex.resize(OrthoGraph::tree->species.size());
}
OrthoExon::~OrthoExon(){
  for(size_t pos = 0; pos < orthoex.size(); pos++){
    delete orthoex[pos];
  }
}
//copy constructor
OrthoExon::OrthoExon(const OrthoExon& other){
  orthoex.resize(other.orthoex.size());
  for(size_t pos = 0; pos < orthoex.size(); pos++){
    if (other.orthoex[pos]){
      orthoex[pos] = new ExonCandidate(*other.orthoex[pos]);
    }
  }
}

//copy with permutation of vector entries
OrthoExon::OrthoExon(const OrthoExon& other, const vector<size_t> &permutation){
  orthoex.resize(other.orthoex.size());
  for(size_t pos = 0; pos < orthoex.size(); pos++){
    if (other.orthoex[pos]){
      orthoex[permutation[pos]] = new ExonCandidate(*other.orthoex[pos]);
    }
  }
}

list<OrthoExon> readOrthoExons(string filename){

  list<OrthoExon> all_orthoex;

  ifstream istrm; 
  istrm.open(filename.c_str(), ifstream::in);
  if (istrm) {
    int nspecies;
    string species;
    vector<size_t> permutation;
    istrm >> goto_line_after( "[SPECIES]");
    istrm >> comment >> nspecies;
    if (nspecies != OrthoGraph::tree->species.size()){
      throw ProjectError("readOrthoExons: number of species in " + filename + 
			 " is " + itoa(nspecies) + ". Number of species in treefile is " + itoa(OrthoGraph::tree->species.size()));
    }
    istrm >> comment;
    for (int i = 0; i < nspecies; i++){
      istrm >> species;
      size_t pos = OrthoGraph::tree->getVectorPositionSpecies(species);
      if (pos == OrthoGraph::tree->species.size()){
	throw ProjectError("readOrthoExons: species name in " + filename + 
			   " is not a species name in treefile.");
      }
      permutation.push_back(pos);
    }
    vector<string> chr(nspecies);
    while(istrm){
      istrm >> goto_line_after( "[CHR]") >> comment;
      for (int i = 0; i < nspecies; i++){
	istrm >> chr[permutation[i]];
      } 
      cout << endl;
      istrm >> goto_line_after( "[ORTHOEX]");
      istrm >> comment;
       while( istrm >> comment >> ws, istrm && istrm.peek() != '[' ){
	 OrthoExon ex_tuple;
	 istrm >> ex_tuple;
	 all_orthoex.push_back(OrthoExon(ex_tuple, permutation));
       }
    } 
  }
  else
    throw ProjectError("readOrthoExons: Could not open this file!");

  return all_orthoex;
}

void writeOrthoExons(const list<OrthoExon> &all_orthoex){
  cout << "# orthologous exons\n" << "#\n" <<"[SPECIES]\n" << "# number of species" << endl;
  cout << OrthoGraph::tree->species.size() << endl;
  cout << "# species names" << endl;
  for (size_t i = 0; i < OrthoGraph::tree->species.size(); i++){
    cout << OrthoGraph::tree->species[i] << "\t";
  }
  cout << endl;
  cout << "#[ORTHOEX]" << endl;
  for(list<OrthoExon>::const_iterator it = all_orthoex.begin(); it != all_orthoex.end(); it++){
    cout << *it << endl;
  }
}

ostream& operator<<(ostream& ostrm, const OrthoExon &ex_tuple){

  int j = 0;
  while ( ex_tuple.orthoex.at(j) == NULL ){
    j++;
  }
  if (j < ex_tuple.orthoex.size()){
    ostrm << stateTypeIdentifiers[(ex_tuple.orthoex.at(j)->getStateType())];
    for (int i = 0; i < ex_tuple.orthoex.size(); i++){
      if (ex_tuple.orthoex.at(i) == NULL){
	ostrm << "\t" << 0 << "\t" << 0;
      }
      else{
	ostrm << "\t" << ex_tuple.orthoex.at(i)->begin+1 << "\t" << ex_tuple.orthoex.at(i)->end - ex_tuple.orthoex.at(i)->begin + 1;
      }
    }
  }
  else{
    cerr<<"Error in writing orthoexon. vector<State*> orthoex only containts null pointers"<<endl;
  }
  return ostrm;
}

istream& operator>>(istream& istrm, OrthoExon& ex_tuple){

  string exontype;
  long int begin, length;

  istrm >> exontype;
  for (int i = 0; i < OrthoGraph::tree->species.size(); i++){
    istrm >> begin >> length;
    if (begin != 0 && length != 0){
      ExonCandidate *exoncand = new ExonCandidate(toExonType(exontype.c_str()), begin-1, begin+length-2);
      ex_tuple.orthoex[i] = exoncand;
    }
  }
  return istrm;
}

map<string, Score> cache::labelscore; //stores score of prunning algorithm for each pattern (leaf labelling)

string OrthoExon::getKey(const OrthoGraph &orthograph){

  labelpattern = "";
  for (size_t i = 0; i < this->orthoex.size(); i++){
    if (this->orthoex.at(i) == NULL){
      labelpattern += "2";
    }
    else {
      map<string, Node*>::iterator it = orthograph.graphs.at(i)->existingNodes.find(orthograph.graphs.at(i)->getKey(this->orthoex.at(i)));
      if (it != orthograph.graphs.at(i)->existingNodes.end()){
	bool label = it->second->label;
	if (label == 1){
	  labelpattern += "1";
	}
	else if (label == 0){
	  labelpattern += "0";
	}
      }
      else
      throw ProjectError("Error in OrthoExon::getKey: exon " + orthograph.graphs.at(i)->getKey(this->orthoex.at(i)) +" not in graph!");
      }
  }
  return labelpattern;
}

bool cache::inHash(string key){
  return ( labelscore.find(key) != labelscore.end() );
}

void cache::addToHash(string key, double score){
  Score s;
  s.treescore = score;
  labelscore[key] = s;

}

void cache::resetCounter(){

  for(map<string, Score>::iterator it = labelscore.begin(); it != labelscore.end(); it++){
    it->second.count = 0;
  }

}

void cache::printCache(list<OrthoExon> &ortho){

  resetCounter();

  cout << "*************************************************************************" << endl;
  cout << "--- orthologous exons + labelpattern ---" << endl;
  for(list<OrthoExon>::iterator it = ortho.begin(); it != ortho.end(); it++){
    incrementCounter(it->labelpattern);
    cout << *it << "\t" << it->labelpattern << endl;
  }
  cout << "\n--- cache summary ---" << endl;
  cout.width(4); cout << "";
  cout.width(12); cout << "score";
  cout.width(4); cout << "#" << endl;

  for(map<string, Score>::iterator it = labelscore.begin(); it != labelscore.end(); it++){

    cout.width(4); cout << it->first;
    cout.width(12); cout << it->second.treescore;
    cout.width(4); cout << it->second.count << endl;

  }
  cout << "*************************************************************************" << endl;
}

double cache::getScore(string key){
  return labelscore[key].treescore;
}
void cache::incrementCounter(string key){
  labelscore[key].count++;
}

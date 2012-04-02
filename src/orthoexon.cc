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
#include "projectio.hh"
#include "types.hh"
#include <fstream>
#include <iostream>

vector<string> OrthoExon::species;

size_t OrthoExon::getVectorPositionSpecies(string name) {
  for(size_t pos=0; pos < species.size(); pos++){
    if (species.at(pos) == name){
      return pos;
    }
  }
  return species.size();
}

map< vector<string>, list<OrthoExon> > readOrthoExons(string filename){

  map< vector<string>, list<OrthoExon> > all_orthoex;

  ifstream istrm; 
  istrm.open(filename.c_str(), ifstream::in);
  if (istrm) {
    int nspecies;
    istrm >> goto_line_after( "[SPECIES]");
    istrm >> comment >> nspecies;
    OrthoExon::species.resize(nspecies);
    istrm >> comment;
    for (int i = 0; i < OrthoExon::species.size(); i++){
      istrm >> OrthoExon::species[i];
    }
    vector<string> chr(nspecies);
    while(istrm){
      istrm >> goto_line_after( "[CHR]") >> comment;
      for (int i = 0; i < nspecies; i++){
	istrm >> chr[i];
      }
      map< vector<string>, list<OrthoExon> >::iterator it = all_orthoex.find(chr);
      if (it == all_orthoex.end()){
	list<OrthoExon> empty_list;
	all_orthoex[chr] = empty_list;
      }  
      cout << endl;
      istrm >> goto_line_after( "[ORTHOEX]");
      istrm >> comment;
       while( istrm >> comment >> ws, istrm && istrm.peek() != '[' ){
	 OrthoExon ex_tuple;
	 istrm >> ex_tuple;
	 all_orthoex[chr].push_back(ex_tuple);
       }
    } 
  }
  else
    throw ProjectError("readOrthoExons: Could not open this file!");

  return all_orthoex;
}

void writeOrthoExons(const map< vector<string>, list<OrthoExon> > &all_orthoex){
  cout << "# orthologous exons\n" << "#\n" <<"[SPECIES]\n" << "# number of species" << endl;
  cout << OrthoExon::species.size() << endl;
  cout << "# species names" << endl;
  for (size_t i = 0; i < OrthoExon::species.size(); i++){
    cout << OrthoExon::species[i] << "\t";
  }
  cout << endl;
  for(map< vector<string>, list<OrthoExon> >::const_iterator it = all_orthoex.begin(); it != all_orthoex.end(); it++){
    cout << "#[CHR]" << endl;
    for (size_t i = 0; i < OrthoExon::species.size(); i++){
      cout << it->first[i] << "\t";
      }
    cout << endl;
    cout << "#[ORTHOEX]" << endl;
    for(list<OrthoExon>::const_iterator ortho = it->second.begin(); ortho != it->second.end(); ortho++){
      cout << *ortho << endl;
    }
  }
}

ostream& operator<<(ostream& ostrm, const OrthoExon &ex_tuple){

  int j = 0;
  while ( ex_tuple.orthoex.at(j) == NULL ){
    j++;
  }
  if (j < ex_tuple.orthoex.size()){
    ostrm << stateTypeIdentifiers[(ex_tuple.orthoex.at(j)->type)];
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
  int begin, length;

  istrm >> exontype;
  for (int i = 0; i < OrthoExon::species.size(); i++){
    istrm >> begin >> length;
    if (begin != 0 && length != 0){
      State *state = new State(begin-1, begin+length-2, toStateType(exontype.c_str()));
      ex_tuple.orthoex.push_back(state);
    }
    else{
     ex_tuple.orthoex.push_back(NULL);
    }
  }
  return istrm;
}

map<string, Score> cache::labelscore; //stores score of prunning algorithm for each pattern (leaf labelling)

string OrthoExon::getKey(const OrthoGraph &orthograph) const{

  string key = "";
  for (size_t i = 0; i < this->orthoex.size(); i++){
    if (this->orthoex.at(i) == NULL){
      key += "2";
    }
    else {
      map<string, Node*>::iterator it = orthograph.graphs.at(i)->existingNodes.find(orthograph.graphs.at(i)->getKey(this->orthoex.at(i)));
      if (it != orthograph.graphs.at(i)->existingNodes.end()){
	bool label = it->second->label;
	if (label == 1){
	  key += "1";
	}
	else if (label == 0){
	  key += "0";
	}
      }
      else
	throw ProjectError("Error in OrthoExon::getKey: exon not in graph!");
      }
  }
  return key;
}

bool cache::inHash(string key){
  return ( labelscore.find(key) != labelscore.end() );
}

void cache::addToHash(string key, double score){
  Score s;
  s.treescore = score;
  s.count = 1;
  labelscore[key] = s;

}

double cache::getScore(string key){
  return labelscore[key].treescore;
}
void cache::incrementCounter(string key){
  labelscore[key].count++;
}

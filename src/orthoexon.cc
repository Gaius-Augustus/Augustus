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

OrthoExon::OrthoExon(){
}

OrthoExon::~OrthoExon(){
}

size_t OrthoExon::getVectorPositionSpecies(string name){
  for(size_t pos=0; pos < species.size(); pos++){
    if (species.at(pos) == name){
      return pos;
    }
  }
  return species.size();
}

list<OrthoExon> readOrthoExons(string filename){

  list<OrthoExon> all_orthoex;

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
      cout << endl;
      istrm >> goto_line_after( "[ORTHOEX]");
      istrm >> comment;
       while( istrm >> comment >> ws, istrm && istrm.peek() != '[' ){
	 OrthoExon ex_tuple;
	 istrm >> ex_tuple;
	 all_orthoex.push_back(ex_tuple);
       }
    } 
  }
  else
    throw ProjectError("readOrthoExons: Could not open this file!");

  return all_orthoex;
}

ostream& operator<<(ostream& ostrm, OrthoExon& ex_tuple){

  int predictionStart;
  try {
    predictionStart = Properties::getIntProperty( "predictionStart" ) - 1;
  } catch (...) {
    predictionStart = 0;
  }

  int j = 0;
  while ( ex_tuple.orthoex.at(j) == NULL ){
    j++;
  }
  if (j < ex_tuple.orthoex.size()){
    ostrm << stateTypeIdentifiers[(((State*)ex_tuple.orthoex.at(j)->item)->type)];
    for (int i = 0; i < ex_tuple.orthoex.size(); i++){
      if (ex_tuple.orthoex.at(i) == NULL){
	ostrm << "\t" << 0 << "\t" << 0;
      }
      else{
	ostrm << "\t" << ex_tuple.orthoex.at(i)->begin+1+predictionStart << "\t" << ex_tuple.orthoex.at(i)->end - ex_tuple.orthoex.at(i)->begin + 1;
      }
    }
  }
  else{
    cerr<<"Error in writing orthoexon. vector<State*> orthoex only containts null pointers"<<endl;
  }
  return ostrm;
}

istream& operator>>(istream& istrm, OrthoExon& ex_tuple){

  int predictionStart;
  try {
    predictionStart = Properties::getIntProperty( "predictionStart" ) - 1;
  } catch (...) {
    predictionStart = 0;
  }

  string exontype;
  int begin, length;

  istrm >> exontype;
  for (int i = 0; i < OrthoExon::species.size(); i++){
    istrm >> begin >> length;
    if (begin != 0 && length != 0){
      State *st = new State(begin-1-predictionStart, begin+length-2-predictionStart, toStateType(exontype.c_str()));
      Status *state = new Status(CDS, begin-1-predictionStart, begin+length-2-predictionStart, 0.0, st);
      ex_tuple.orthoex.push_back(state);
    }
    else{
     ex_tuple.orthoex.push_back(NULL);
    }
  }
  return istrm;
}

map<string, Score> labelscore; //stores score of prunning algorithm for each pattern (leaf labelling)

string OrthoExon::getKey(){

  string key = "";
  for (size_t i=0; i < this->orthoex.size(); i++){
    if (this->orthoex.at(i) == NULL){
      key += "2";
    }
    //TODO:
    /*else if(this->orthoex.at(i)->label == 1){
      key += "1";
    }
    else if(this->orthoex.at(i)->label == 0){
      key += "0";
      }*/
  }
  return key;
}

bool OrthoExon::inHash(){
  return ( labelscore.find(this->getKey()) != labelscore.end() );
}

void OrthoExon::addToHash(double score){
  Score s;
  s.treescore = score;
  s.count = 1;
  labelscore[this->getKey()] = s;

}

double OrthoExon::getScore(){
  return labelscore[this->getKey()].treescore;
}
void OrthoExon::incrementCounter(){
  labelscore[this->getKey()].count++;
}

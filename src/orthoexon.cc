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
#include <fstream>
#include <iostream>

vector<string> OrthoExon::species;

OrthoExon::OrthoExon(){
}

size_t OrthoExon::getVectorPositionSpecies(string name){
  for(size_t pos=0; pos < species.size(); pos++){
    if (species.at(pos) == name){
      return pos;
    }
  }
  return species.size();
}

void readOrthoExons(string filename){
  ifstream istrm; 
  istrm.open(filename.c_str(), ifstream::in);
  if (istrm) {
    int nspecies;
    string exontype;
    int begin, length;
    istrm >> goto_line_after( "[SPECIES]");
    istrm >> comment >> nspecies;
    vector<string> species(nspecies);
    istrm >> comment;
    for (int i = 0; i < species.size(); i++){
      istrm >> species[i];
    }
    vector<string> chr(nspecies);
    while(istrm){
      istrm >> goto_line_after( "[CHR]") >> comment;
      for (int i = 0; i < species.size(); i++){
	istrm >> chr[i];
      }
      for (int i = 0; i < species.size(); i++){
	cout<<chr[i]<<"\t";
      } 
      cout << endl;
      istrm >> goto_line_after( "[ORTHOEX]");
      istrm >> comment;
       while( istrm >> comment >> ws, istrm && istrm.peek() != '[' ){
	istrm >> exontype;
	for (int i = 0; i < nspecies; i++){
	  istrm >> begin >> length;
	  if (begin!=0 && length!=0){
	    string key = itoa(begin-1) + ":" + itoa(begin+length-2) + ":" + exontype;
	    cout<<key<<endl; 
	    //search map existingNodes for key. in case that the key is not in
	    //in the map -> cerr key. might be a typo in the input file
	    //else set pointer to node
	  }
	  else{
	    //set pointer to NULL
	  }
	}
      }
    }
    
  }
  else
  throw ProjectError("Could not open input file " + filename);
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

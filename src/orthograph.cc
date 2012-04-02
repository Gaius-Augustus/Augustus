/**********************************************************************
 * file:    orthograph.cc
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  orthologous graphs for comparative gene prediction
 * authors: Stefanie König
 *
 * date    |   author      |  changes
 * --------|---------------|------------------------------------------
 * 23.03.12|Stefanie König | creation of the file
 **********************************************************************/


#include "orthograph.hh"
#include "orthoexon.hh"

OrthoGraph::OrthoGraph(){
  graphs.resize(OrthoExon::species.size());
}

OrthoGraph::~OrthoGraph(){
  for(int i = 0; i < graphs.size(); i++){
    delete graphs[i];
  }
}

void OrthoGraph::addSingleGraph(string species, AugustusGraph* singleGraph){
  size_t pos = OrthoExon::getVectorPositionSpecies(species);
  if (pos < this->graphs.size()){
    this->graphs[pos] = singleGraph;
  }
  else
    cerr << "species names in Orthograph and OrthoExon don't match" << endl;
}

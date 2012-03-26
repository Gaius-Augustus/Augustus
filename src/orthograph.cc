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

OrthoGraph::OrthoGraph(list<Status> *states, int dnalength, list<Status> *additionalExons){
    AugustusGraph *graph = new AugustusGraph(states, dnalength);
    graph->buildGraph(additionalExons);
    graphs.push_back(graph);
}

OrthoGraph::~OrthoGraph(){
  for(int i = 0; i < graphs.size(); i++){
    delete graphs[i];
  }
}

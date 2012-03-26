/**********************************************************************
 * file:    orthograph.hh
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  orthologous graphs for comparative gene prediction
 * authors: Stefanie König
 *
 *********************************************************************/

#ifndef _ORTHOGRAPH_HH
#define _ORTHOGRAPH_HH

#include "graph.hh"

using namespace std;

class OrthoGraph{

public:
  OrthoGraph(list<Status> *states, int dnalength, list<Status> *additionalExons);
  ~OrthoGraph();
  vector<AugustusGraph*> graphs;
};

#endif

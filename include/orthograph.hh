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
#include "phylotree.hh"
#include "orthoexon.hh"
#include "randseqaccess.hh"

using namespace std;

class OrthoGraph{

public:
  OrthoGraph(RandSeqAccess *rsa);
  ~OrthoGraph();

  vector<AugustusGraph*> graphs;
  vector<AnnoSequence*> orthoSeqRanges;
  static PhyloTree *tree;
  static vector<ofstream*> filestreams;
  list<OrthoExon> all_orthoex;
  vector< list<Gene> *> ptrs_to_alltranscripts; //stores pointers to alltranscripts until they can be deleted (destructor of OrthoGraph)
  double score;


  static void initOutputFiles();
  static void closeOutputFiles();
  void optimizeLabels();
  void outputGenes(Strand strand);
  inline void storePtrsToAlltranscripts(list<Gene> *alltranscripts){
    this->ptrs_to_alltranscripts.push_back(alltranscripts);
  }
};
#endif

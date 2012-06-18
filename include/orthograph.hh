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

#include "speciesgraph.hh"
#include "phylotree.hh"
#include "orthoexon.hh"

using namespace std;

class OrthoGraph{

public:
  OrthoGraph();
  ~OrthoGraph();

  static size_t numSpecies;               //the number of species
  vector<SpeciesGraph*> graphs;
  vector<AnnoSequence*> orthoSeqRanges;
  static PhyloTree *tree;
  static vector<ofstream*> filestreams;
  list<OrthoExon> all_orthoex;
  vector< list<Gene> *> ptrs_to_alltranscripts; //stores pointers to alltranscripts until they can be deleted (destructor of OrthoGraph)

  //functions to redirect filestreams
  static void initOutputFiles();
  static void closeOutputFiles();


  //optimization functions
  void optimize();
  void localMove(vector<MoveObject*> &orthomove);
  double calculateScoreDiff(vector<MoveObject*> &orthomove);
  inline double pruningAlgor(){                    // pruning Algor. for all OrthoExons
    return pruningAlgor(all_orthoex);
  }                          
  double pruningAlgor(list<OrthoExon> &orthoex);   // pruning Algor. for a list of OrthoExons in a range
  list<OrthoExon> orthoExInRange(vector<MoveObject*> &orthomove); //determine all OrthoExons in a range
  void addOrthoIntrons(vector<MoveObject*> &orthomove, list<OrthoExon> &local_orthoexons);

  void outputGenes(Strand strand);
  inline void storePtrsToAlltranscripts(list<Gene> *alltranscripts){
    this->ptrs_to_alltranscripts.push_back(alltranscripts);
  }
};

#endif

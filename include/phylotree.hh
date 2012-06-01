/**********************************************************************
 * file:    phylotree.hh
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  representation of the phylogenetic relationship of the species
 *          in comparative gene prediction
 * authors: Stefanie König
 *
 *********************************************************************/

#ifndef _PHYLOTREE_HH
#define _PHYLOTREE_HH

#include <list>
#include <string>
#include <vector>
#include "graph.hh"



//forward declarations
class Treenode;
class PhyloTree;
class OrthoGraph;
class OrthoExon;
class ExonEvoModel;


class Treenode{

private:
  string species;
  double dist_to_parent;
  Treenode *parent;
  list<Treenode*> children;
  vector<double> alpha; //DB variable of Prunning Algorithm
 
public:  
  Treenode(string s="interior", double t=0.0, Treenode *p=NULL):
    species(s),
    dist_to_parent(t),
    parent(p)
  {
    alpha.resize(2);
  }
  ~Treenode(){}

  inline void addDistance(double dist){
    this->dist_to_parent=dist;
  }

  inline void addChild(Treenode *child){
    child->parent=this;
    this->children.push_back(child);
  }

  inline bool isLeaf() const {
    return (this->children.empty());
  }

  inline bool isRoot() const {
    return (this->parent == NULL);
  }

  void printNode() const;
  double calculateAlphaScore(bool label, ExonEvoModel &evo);

  friend class PhyloTree;
  
};

class ExonEvoModel{

public:
  ExonEvoModel();
  ~ExonEvoModel() {}
  double P(bool label1, bool label2, double dist) const;  //calculates probability of transition from label1 to label2 in time dist
  double getEquilibriumFreq(bool label) const;  //equilibrium frequency of a label

  //private:
  double mu;  // rate for exon loss
  double lambda; // rate for exon gain

};

class PhyloTree{

private:
  list<Treenode*> treenodes; // leaf to root order!

public:
  PhyloTree(string filename);
  ~PhyloTree();

  vector<string> species;
  ExonEvoModel evo;
  size_t getVectorPositionSpecies(string name);
  void printTree() const;
  void printWithGraphviz(string filename) const;
  double pruningAlgor(OrthoExon &orthoex, const OrthoGraph &orthograph);

};

#endif

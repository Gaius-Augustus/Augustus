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
#include "orthoexon.hh"

using namespace std;

//forward declarations
class Treenode;
class PhyloTree;


class Treenode{
 
public:  
  Treenode(string s="interior", double t=0.0, Treenode *p=NULL):
    species(s),
    dist_to_parent(t),
    parent(p)
  {
    alpha.resize(2);
  }
  string species;
  double dist_to_parent;
  Treenode *parent;
  list<Treenode*> children;
  vector<double> alpha; //DB variable of Prunning Algorithm

  inline void addDistance(double dist){
    this->dist_to_parent=dist;
  }

  inline void addChild(Treenode *child){
    child->parent=this;
    this->children.push_back(child);
  }

  inline bool isLeaf(){
    return (this->children.empty());
  }

  inline bool isRoot(){
    return (this->parent == NULL);
  }
  void printNode();
  
};

class PhyloTree{

public:
  PhyloTree(){}
  ~PhyloTree();
  list<Treenode*> treenodes; // leaf to root order!
  void printTree();
  void printWithGraphviz(string filename);
  double pruningAlgor(OrthoExon &orthoex);
  double P(bool label1, bool label2, double dist) const;  //calculates probability of transition from label1 to label2 in time dist
  double getAlphaScore(Treenode* node, bool label);

};

PhyloTree parseTree(string filename);

#endif

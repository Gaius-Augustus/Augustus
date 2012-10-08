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


//forward declarations
class Treenode;
class PhyloTree;
class ExonEvoModel;


class Treenode{

private:
    std::string species;             // speciesname (only leaf nodes)
    double dist_to_parent;           // branch length to parent
    Treenode *parent;                // parent node
    std::list<Treenode*> children;   // child nodes
    std::vector<double> alpha;       /* stores dynamic programing variable of the pruning algorithm:
				      * alpha[i] is the conditional probability of the label pattern 
				      * in the subtree rooted in this node given this node has label i,
				      * where i is from the state space { 0,1 }
				      */
 
public:  
    Treenode(std::string s="interior", double t=0.0, Treenode *p=NULL):
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
    double calculateAlphaScore(bool label, ExonEvoModel &evo); // calculates the dynamic programing variables

    friend class PhyloTree;
  
};

/*
 * model of exon evolution
 */

class ExonEvoModel{

public:
    ExonEvoModel();
    ~ExonEvoModel() {}
    double P(bool label1, bool label2, double dist) const;  //returns probability of transition from label1 to label2 along a branch of length dist
    double getEquilibriumFreq(bool label) const;  //equilibrium frequency of a label
    inline double getPhyloFactor() const{
	return phylo_factor;
    }
    inline double getMu() const{
	return mu;
    }
    inline double getLambda() const{
	return lambda;
    }
private:
    double mu;           // rate for exon loss
    double lambda;       // rate for exon gain
    double phylo_factor; // tuning parameter to fit phylogenetic score to scores in the graph


};

/*
 * the phylogenetic tree
 */

class PhyloTree{

private:
    std::list<Treenode*> treenodes; // leaf to root order!

public:
    PhyloTree(std::string filename);
    ~PhyloTree();

    std::vector<std::string> species;  //vector of species names
    ExonEvoModel evo;
    size_t getVectorPositionSpecies(std::string name);
    void printTree() const;
    void printWithGraphviz(std::string filename) const;
    double pruningAlgor(std::string labelpattern);   //executes pruning algorithm for the tree and a the given model of exon evolution

};

#endif

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

// project includes
#include "contTimeMC.hh"

//forward declarations
class Treenode;
class PhyloTree;
class Evo;
class ExonEvo;

typedef std::vector<bool> bit_vector;

class Treenode{

private:
    //tree topology
    std::string species;             // speciesname (only leaf nodes)
    double distance;                 // branch length to parent
    int idx;                         // index of speciesname as in RandSeqAccess:speciesIndex (interior nodes: -1)
    Treenode *parent;                // parent node
    std::list<Treenode*> children;   // child nodes 
    bool label;
    gsl_matrix *logP;                // log probability matrix for branch length t='distance'    

    /*
     * tables needed for storing recursion variables
     * during pruning algorithm and MAP inference
     *
     * i-th entry is the best assignment given that the parent is assigned to i
     * only needed for MAP inference
     */
    std::vector<int> bestAssign;
    /* 
     * recursion variables:
     * Felsenstein pruning algorithm: i-the entry is the probability of the subtree rooted in this node 
     * given that the root is assigned to i
     * MAP inference: i-the entry is the probability of the most likely assignment of all nodes of the subtree 
     * rooted in this node given that the root is assigned to i
     */
    std::vector<double> table;
 
public:  
    Treenode(std::string s="i", double t=0.0, int i=-1, Treenode *p=NULL):
	species(s),
	distance(t),
	idx(i),
	parent(p),
	label(0),
	logP(NULL)
    {}
    ~Treenode(){}
    void addDistance(double dist){this->distance=dist;}
    void addChild(Treenode *child); 
    void removeChild(Treenode *child){children.remove(child);}
    void makeRoot();
    bool isLeaf() const {return (this->children.empty());}
    bool isRoot() const {return (this->parent == NULL);}
    void setTable(int pos, double value){table.at(pos)=value;}
    double getTable(int pos) const{return table.at(pos);}
    double getDist() const {return distance;}
    int getIdx() const {return idx;}
    void setIdx(int i) {this->idx = i;}
    gsl_matrix* getLogP() const {return logP;}
    void setLogP(gsl_matrix *m) {logP = m;}
    void resizeTable(int size, double val=0);
    std::string getSpecies() const {return species;}
    Treenode *getParent() const {return this->parent;}
    void printNode() const;

    friend class PhyloTree;
};

/*
 * the phylogenetic tree
 */
class PhyloTree{

private:
    std::list<Treenode*> treenodes; // leaf to root order!
    int numSp; // number of species
public:
    PhyloTree(std::string filename);
    //create star-like tree with unit branch length b from a set of species identifiers
    PhyloTree(const std::vector<std::string> &speciesnames, double b=1);
    //copy constructor
    PhyloTree(const PhyloTree &other);
    ~PhyloTree();

    /*
     * general look-up functions
     */
    void getBranchLengths(std::vector<double> &branchset) const;
    void scaleTree(double factor); // scale tree by multiplying each branch length by a factor 
    double getDiversity(); // sum of branch lengths
    void getSpeciesNames(std::vector<std::string> &speciesnames);
    int numSpecies(){return numSp;}
    Treenode *getLeaf(std::string species) const;
    void printTree() const;
    void printNewick(std::string filename) const; // output tree in Newick format
    void recursiveNWK(std::ofstream &file, Treenode *node) const; // subroutine of printNewick()

    /*
     * functions to remove terminal branches and possibly interior nodes
     * drop removes a terminal branch given a leaf node or the species name
     * if after removing a terminal branch the junction node has only one child left, the node
     * will be collapsed and the branch length will be added to its child
     */
    void drop(std::string species);  
    void drop(Treenode *node, Evo *evo);  
    void collapse(Treenode *node, Evo *evo);
    /*
     * prune all leaf nodes of species that are not present in the tree as indicated by a bit vector
     * (i-th bit in the vector is 1 if species i is present and 0 if species i is absent)
     */
    void prune(bit_vector &bv, Evo *evo); 
    /*
     * Felsensteins Pruning Algorithm
     */
    double pruningAlgor(std::string labelpattern, Evo *evo, int u=0);
    double pruningAlgor(std::vector<int> &tuple, Evo *evo, int u=0);
    void printRecursionTable() const;

    /*
     * MAP inference given a set of weights for leaf nodes
     * if the flag 'fixLeafLabels' is turned on MAP inference is carried out
     * with the leaf labels fixed to the corresponding labels in the hect
     * (needed in 'optimization via dual decomposition' to make the solutions
     * of the two subproblems consistent)
     */
    double MAP(std::vector<int> &labels, std::vector<double> &weights, Evo *evo, double k=1.0, bool fixLeafLabels=false);
    void MAPbacktrack(std::vector<int> &labels, Treenode* root, int bestAssign, bool fixLeafLabels);
};


/*
 * Locus tree
 * Charlotte Janas playground
 */
class LocusTree {

};

#endif

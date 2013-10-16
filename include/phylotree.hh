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
class Evo;
class ExonEvo;
class OrthoExon;
class RandSeqAccess;

class Treenode{

private:
    //tree topology
    std::string species;             // speciesname (only leaf nodes)
    double distance;                 // branch length to parent
    Treenode *parent;                // parent node
    std::list<Treenode*> children;   // child nodes 
    bool label;

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
    Treenode(std::string s="i", double t=0.0, Treenode *p=NULL):
	species(s),
	distance(t),
	parent(p),
	label(0)
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
    std::vector<std::string> species;  //vector of species names
    static RandSeqAccess *rsa;
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
    static void setRSA(RandSeqAccess *r){rsa=r;}
    void getBranchLengths(std::vector<double> &branchset) const;
    void getSpeciesNames(std::vector<std::string> &speciesnames) const {speciesnames=this->species;}
    int numSpecies() const {return species.size();}
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
    void drop(Treenode *node, Evo *evo=NULL);  
    void collapse(Treenode *node, Evo *evo=NULL);

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

    // calculate diversity (sum of branch lengths) of the subtree induced by a HECT
    double sumBranches(OrthoExon &oe);
};

#endif

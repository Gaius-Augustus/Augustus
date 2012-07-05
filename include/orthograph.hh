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

class OrthoExon;

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

    string getLabelpattern(OrthoExon &ex); //determines the key of an orthoex for the map labelscore

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


    //functions to create different types of moves:
  
    /* 
     * majorityRule: if # of 0's in labelpattern of an OrthoExon is smaller than # of 1's,
     * all nodes with the labels 0 are made to 1 and vice versa
     */
    vector<MoveObject*> majorityRuleMove(OrthoExon *orthoex);
    vector<MoveObject*> allToOne(OrthoExon *orthoex);

    void outputGenes(Strand strand);
    inline void storePtrsToAlltranscripts(list<Gene> *alltranscripts){
	this->ptrs_to_alltranscripts.push_back(alltranscripts);
    }
};

struct Score{
    double treescore;   //stores the score of a label pattern
 
    int count;          //counts the number of exontuples which have that specific pattern

    Score() : treescore(0), count(0) {}
    ~Score() {}
};

/*
 * hashfunction storing all label patterns and their Score
 * key: string over alphabet {0,1,2}^k, k = # species^k
 * 0 codes for exon in graph, but not part of the maximum weight path
 * 1 codes for exon in graph and part of the maximum weight path
 * 2 codes for exon not in graph: Status* = NULL
 * first digit in string refers to first species in PhyloTree::species, second digit refers to second species, ...
 */

namespace cache{

    extern map<string, Score> labelscore;
    /*
     * cache functions
     */
    bool inHash(string key);
    void resetCounter();
    void addToHash(string key, double score);
    double getScore(string key);
    void incrementCounter(string key);
    void printCache(list<OrthoExon> &ortho);
}


#endif

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

//forward declarations
class OrthoExon;

class OrthoGraph{

public:
    OrthoGraph(){
	graphs.resize(numSpecies);
	ptrs_to_alltranscripts.resize(numSpecies);
	print_change = false;

	try {
	    maxIterations = Properties::getIntProperty("/CompPred/maxIterations");
	} catch (...) {
	    maxIterations = 0;
	}
	try {
	    shift_size = Properties::getIntProperty("/CompPred/shift_size");
	} catch (...) {
	    shift_size = 1;
	}

    }
    ~OrthoGraph(){
	for(int i = 0; i < numSpecies; i++){
	    delete graphs[i];
	    delete ptrs_to_alltranscripts[i];
	}
    }

    static size_t numSpecies; //the number of species
    vector<SpeciesGraph*> graphs;
    static PhyloTree *tree;
    list<OrthoExon> all_orthoex;
    vector< list<Gene> *> ptrs_to_alltranscripts; //stores pointers to alltranscripts until they can be deleted (destructor of OrthoGraph)
    int maxIterations;  // max number a move can be repeated
    int shift_size;     // number of exons local_head/local_tail is shifted to the left/right on the current path

    bool print_change; //only temporary flag

    void addScoreSelectivePressure(); //adds score for selective pressure to all nodes representing orthologous exons
    void globalPathSearch(); //determine initial labeling of the graphs
    string getLabelpattern(OrthoExon &ex); //determines the current labelpattern of an orthoex
    void printHTMLgBrowse(OrthoExon &ex);  //temp: html output for gBrowse
    inline void printCache(){
	printCache(all_orthoex);
    }
    void printCache(list<OrthoExon> &ortho);

    //optimization functions
    void optimize();                                 // create MoveObjects
    void localMove(vector<Move*> &orthomove);  // do the local move for a MoveObject
    inline double pruningAlgor(){                    // pruning Algor. for all OrthoExons
	return pruningAlgor(all_orthoex);
    }               
    double pruningAlgor(list<OrthoExon> &orthoex);   // pruning Algor. for a list of OrthoExons in a range
    list<OrthoExon> orthoExInRange(vector<Move*> &orthomove); //determine all OrthoExons in a range
   


    //functions to create different types of moves:
  
    /* 
     * majorityRule: if # of 0's in labelpattern of an OrthoExon is smaller than # of 1's,
     * all nodes with the labels 0 are made to 1 and vice versa
     */
    vector<Move*> majorityRuleMove(OrthoExon &orthoex);

    void outputGenes(vector<ofstream*> filestreams, vector<int> &geneid);
    inline void storePtrsToAlltranscripts(list<Gene> *alltranscripts){
	this->ptrs_to_alltranscripts.push_back(alltranscripts);
    }
};

struct Score{
    double treescore;   //stores the score of a label pattern
 
    int count;          //counts the number of exon candidate tuples which have that specific pattern

    Score() : treescore(0), count(0) {}
    ~Score() {}
};

/*
 * hashfunction storing all label patterns and their Score
 * labelpattern: string over alphabet {0,1,2}^k, k = # species^k
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
    bool inHash(string labelpattern);
    void resetCounter();
    void addToHash(string labelpattern, double score);
    double getScore(string labelpattern);
    void incrementCounter(string labelpattern);
}

//functions to redirect filestreams
vector<ofstream*> initOutputFiles(string extension = string());
void closeOutputFiles(vector<ofstream*> filestreams);


#endif

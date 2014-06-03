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

//project includes
#include "speciesgraph.hh"
#include "orthoexon.hh"
#include "phylotree.hh"

class OrthoGraph{

public:
    OrthoGraph(){
	graphs.resize(numSpecies);
	ptrs_to_alltranscripts.resize(numSpecies);
	sfcs.resize(numSpecies);
	print_change = false;
    }
    ~OrthoGraph(){
	for(int i = 0; i < numSpecies; i++){
	    delete graphs[i];
	    delete ptrs_to_alltranscripts[i];
	    delete sfcs[i];
	}
    }

    static size_t numSpecies; //the number of species
    vector<SpeciesGraph*> graphs;
    static PhyloTree *tree;
    list<OrthoExon> all_orthoex;
    vector< list<Gene> *> ptrs_to_alltranscripts; //stores pointers to alltranscripts until they can be deleted (destructor of OrthoGraph)
    vector< SequenceFeatureCollection* > sfcs;  // stores extrinsic evidence for each species

    bool print_change; // only temporary flag
    void linkToOEs(list<OrthoExon> &oes); // link ECs in HECTS to nodes in orthograph 
    void addScoreSelectivePressure(); //const. reward for orthologous exons and const. penalty for non-orthologous exons. Only temporary until PAML is integrated.

    /*
     * optimization by making small local changes called moves
     */
    void optimize(ExonEvo &evo);// main routine, in which different moves are created
    void localMove(vector<Move*> &orthomove, ExonEvo &evo, int shift_size);  // execution a single move
    /*
     * currently, only a single type of move is implemented:
     * if # of 0's in labelpattern of an OrthoExon is smaller than # of 1's,
     * all nodes with the labels 0 are made to 1 and vice versa
     */
    vector<Move*> majorityRuleMove(OrthoExon &orthoex, int shift_size);
    double pruningAlgor(ExonEvo &evo){                    // pruning algor. for all OrthoExons
	return pruningAlgor(all_orthoex, evo);
    }               
    double pruningAlgor(list<OrthoExon> &orthoex, ExonEvo &evo);   // pruning algor. for a list of OrthoExons
    list<OrthoExon> orthoExInRange(vector<Move*> &orthomove); //determine all OrthoExons in a range
    string getLabelpattern(OrthoExon &ex); //determines the current labelpattern of an orthoex
    void printHTMLgBrowse(OrthoExon &ex);  //temp: html output for gBrowse
    inline void printCache(){
	printCache(all_orthoex);
    }
    void printCache(list<OrthoExon> &ortho);

    /*
     * optimization via dual decomposition
     * problem is decomposed into two subproblems whith can be solved efficiently:
     * horizontal problem: DAG longest path
     * verical problem: MAP inference on a set of disjoint phylogenetic trees whose
     * leaf nodes are assigned to weights.
     */
    double dualdecomp(ExonEvo &evo,vector< list<Gene> *> &genelist, int gr_ID, int T, double c);  //main routine
    double treeMAPInf(ExonEvo &evo, int &numInconsistent);  //vertical problem
    double globalPathSearch(); // horizontal problem
    double getStepSize(double c,int t, int v);    // specifies a sequence of steps
    double makeConsistent(ExonEvo &evo);
    void printSummary();

    // transform graph labeling into list of genes + filter + output
    void buildGeneList(vector< list<Gene>* > &genelist);
    void filterGeneList(vector< list<Gene> *> &genelist, vector<ofstream*> &filestreams, vector<int> &geneid);
    void outputGenes(vector<ofstream*> &filestreams, vector<int> &geneid){
	vector< list<Gene> *> genelist(numSpecies);
	buildGeneList(genelist);
	filterGeneList(genelist,filestreams, geneid);
    }
};

struct Score{
    double treescore;   //stores the score of a label pattern
 
    int count;          //counts the number of exon candidate tuples which have that specific pattern

    Score() : treescore(0), count(0) {}
    ~Score() {}
};

/*
 * hashfunction storing all label patterns and their score
 * labelpattern: string over alphabet {0,1,2}^k, k = # species^k
 * the i-th character in a label pattern is
 * 0 if the exon in the i-th species has label 0
 * 1 if the exon in the i-th species has label 1
 * 2 if exon in the i-th species does not exist
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
vector<ofstream*> initOutputFiles(string outdir, string extension = string());
void closeOutputFiles(vector<ofstream*> filestreams);


#endif

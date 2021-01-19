/*
 * compgenepred.hh
 *
 * License: Artistic License, see file LICENSE.TXT or 
 *          https://opensource.org/licenses/artistic-license-1.0
 */

#ifndef _COMPGENEPRED_HH
#define _COMPGENEPRED_HH

// project includes
#include "extrinsicinfo.hh"
#include "randseqaccess.hh"
#include "phylotree.hh"
#include "geneMSA.hh"
#include "genomicMSA.hh"
#include "namgene.hh"



#include <stack>


// addeded for inetrval merging of generanges belonging to the same species 
#ifdef TESTING
void mergeInterval(list<tuple<string,int,int> >& interval);
void mergeIntervals(vector<string>& speciesNames, vector<list<tuple<string,int,int> > >& intervals);
bool sortInterval(const tuple<string,int,int>& a, const tuple<string,int,int>& b);
void writeIntervals(string dirname, vector<string>& speciesNames, vector<list<tuple<string,int,int> > >& intervals);
bool shiftGFF(string filename);
#endif

/**
 * @brief comparative gene prediction on multiple species
 * 
 * @author Mario Stanke
 * @author Alexander Gebauer
 * @author Stefanie Koenig
 */
class CompGenePred {
    public:
        CompGenePred();
        ~CompGenePred() { delete rsa;}

        void start();    
        void runPredictionOrTest();

        #ifdef TESTING
        // helpers for testing
        void postprocTest();
        bool readInterval(string filename, list<tuple<string,int,int> >& grlist);
        #endif
        
        RandSeqAccess *rsa;
        PhyloTree tree;

 

};

#endif  // _COMPGENEPRED_HH



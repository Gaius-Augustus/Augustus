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

    #ifdef TESTING
    // helpers for testing
    void prepareTest();   
    void runTest();
    void postprocTest();
    bool readInterval(string filename, list<tuple<string,int,int> >& grlist);
    #endif

    void start();
    RandSeqAccess *rsa;
    PhyloTree tree;
};

#endif  // _COMPGENEPRED_HH

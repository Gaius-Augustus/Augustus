/*
 * compgenepred.hh
 *
 * License: Artistic License, see file LICENSE.TXT or 
 *          https://opensource.org/licenses/artistic-license-1.0
 */

#ifndef _COMPGENEPRED
#define _COMPGENEPRED

// project includes
#include "extrinsicinfo.hh"
#include "randseqaccess.hh"
#include "phylotree.hh"


/*
*   added by Giovanna Migliorelli 10.10.2019 
*   headers for functions in charge of printing MSA associated with some OE
*/
void ec2Chr(int ecStart, int ecEnd, int offset, int chrLen, Strand strand, int& chrStart, int& chrEnd, bool referToPlus);
void printMSA_printEC_helper(ostream& os, GeneMSA* geneRange, int sp, vector<string>& speciesNames, RandSeqAccess* rsa, int chrStart, int chrEnd, AnnoSequence* as, bool expandUnaligned);
void printMSA_printEC(ostream& os, vector<string>& speciesNames, RandSeqAccess* rsa, GeneMSA* geneRange, int sp, ExonCandidate* ec, bool downcase, bool expandUnaligned);
void printMSA(ostream& os, vector<string>& speciesNames, RandSeqAccess* rsa, GeneMSA* geneRange, OrthoExon& oe, bool extraInfo, bool expandUnaligned);


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
    RandSeqAccess *rsa;
    PhyloTree tree;
};

#endif  // _COMPGENEPRED


/**********************************************************************
 * file:    compgenepred.hh
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  comparative gene prediction on multiple species
 * authors: Mario Stanke
 *
 *********************************************************************/

#ifndef _COMPGENEPRED
#define _COMPGENEPRED

// project includes
#include "extrinsicinfo.hh"
#include "randseqaccess.hh"
#include "phylotree.hh"

class CompGenePred {
public:
    CompGenePred();
    ~CompGenePred() { delete rsa;}

    void start();
    RandSeqAccess *rsa;
    PhyloTree tree;
};

#endif  // _COMPGENEPRED

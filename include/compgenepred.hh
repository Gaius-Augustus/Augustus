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

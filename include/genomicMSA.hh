/**********************************************************************
 * file:    genomicMSA.hh
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  multiple sequence alignment of genomes for comparative gene prediction
 * authors: Mario Stanke, Alexander Gebauer
 *
 *********************************************************************/

#ifndef _GENOMICMSA
#define _GENOMICMSA

#include "exoncand.hh"
#include "geneMSA.hh"
#include "randseqaccess.hh"

class GenomicMSA {
public:
    GenomicMSA(RandSeqAccess *rsa) {this->rsa = rsa;}
   ~GenomicMSA(){}

    void readAlignment(string alignFilename); // reads a multiple species alignment from a *.maf file
    void printAlignment(string outFname); // print alignment in .maf format, to stdout if outFname is empty string
    int numAlignments() { return alignment.size(); }
 
    GeneMSA *getNextGene();
    /** 
     * chromosome lengths as specified in the maf file
     * for each species on hash with sequence names as keys and lengths as values
     */
    vector<map<string,int>> chrLen;
    /**
     * merges pairs of alignments in order to reduce the alignment number in trivial cases
     * without doing any potentially false mergers
     */
    void compactify();
    /**
     * changes alignment list, so that afterwards, each alignment contains a gene range:
     * a single alignment that may contain one or more gene, usually the merger of many .maf alignments
     */
    void findGeneRanges();
private:
    list<Alignment*> alignment;
    int numSpecies;
    RandSeqAccess *rsa; // the actual data is manages in CompGenePred
};

#endif  // _GENOMICMSA

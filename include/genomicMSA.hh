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
   /*
     * merges the merged alignment parts so that a gene is possibly in this segment
     * returns NULL if alignment is empty
     */
    void prepareExons() {
        mergeAlignment (6,0.6); // make this reasonable after experience with the data
    }
    GeneMSA *getNextGene();
    /** 
     * chromosome lengths as specified in the maf file
     * for each species on hash with sequence names as keys and lengths as values
     */
    vector<map<string,int>> chrLen; 
private:
    list<Alignment*> alignment;
    /*
     * merges neighboring AlignmentBlocks in alignment if all species are at most "maxGapLen" bases apart and
     * at least 'percentSpeciesAligned' percent of the species are contained in the alignment parts
     * Purpose: An ortho exon will be assumed to be contained in a single (merged) alignment block.
     */
    void mergeAlignment(int maxGapLen, float percentSpeciesAligned);
    int numSpecies;
    RandSeqAccess *rsa; // the actual data is manages in CompGenePred
};

#endif  // _GENOMICMSA

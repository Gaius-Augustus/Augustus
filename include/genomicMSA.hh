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


class GenomicMSA {
public:
    GenomicMSA() {}
   ~GenomicMSA(){}

   void readAlignment(string alignFilename, vector<string> speciesnames); // reads a multiple species alignment from a *.maf file
   void printAlignment(string outFname); // print alignment in .maf format, to stdout if outFname is empty string
   /*
     * merges the merged alignment parts so that a gene is possibly in this segment
     * returns NULL if alignment is empty
     */
    void prepareExons() {
        mergeAlignment (6,0.6); // make this reasonable after experience with the data
    }
    GeneMSA *getNextGene();

private:
    list<AlignmentBlock*> alignment;
    /*
     * merges multiple alignment parts if all species are less than "maxGapLen" bases apart and
     * 'percentSpeciesAligned' percent of the species are contained in the alignment parts
     */
    void mergeAlignment(int maxGapLen, float percentSpeciesAligned);

};

#endif  // _GENOMICMSA

/**********************************************************************
 * file:    geomicMSA.hh
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  multiple sequence alignment of genomes for comparative gene prediction
 * authors: Mario Stanke, Alexander Gebauer
 *
 *********************************************************************/

#ifndef _GENEMSA
#define _GENEMSA

// project includes
#include "exoncand.hh"
#include "orthoexon.hh"

struct AlignmentBlock {
    vector<AlignSeq*> alignSpeciesTupel;
};

class GeneMSA {
public:
    static ofstream *exonCands_outfile, *orthoExons_outfile; // pointer to the output files
    list<AlignmentBlock*> alignment;		// list of the alignment parts which possibly belong to a gene segment
    vector< list<ExonCandidate*>* > exoncands;		// exon candidates found in the different species in a gene segment
    vector< map<string, ExonCandidate*>* > existingCandidates;		// stores the keys of the exon candidates for the different species
    list<OrthoExon> orthoExonsList;		// list of orthologue exons found in a gene segment

    GeneMSA() {};
    ~GeneMSA();

    int getStart(int speciesIdx);
    int getEnd(int speciesIdx);
    string getChr(int speciesIdx);
    Strand getStrand(int speciesIdx);

    map<string,ExonCandidate*>* addToHash(list<ExonCandidate*> *ec); // adds the keys to the map function

    /*
     * createExonCands: get all exon candidates
     * assqthresh, dssqthresh are between 0 and 1 and thresholds for the inclusion of
     * acceptor/donor splice sites based on the pattern probability
     * assqthresh=0.05 means that only acceptor ss are considered
     * that have a pattern, such that 5% of true splice site patterns have lower probability.
     * The default threshold of 0 means that all splice site patterns are considered.
     */
    //void createExonCands(const char *dna, float assqthresh, float dssqthresh);
    void createExonCands(const char *dna, float assqthresh, float dssqthresh, float motifqthresh);

    pair<int,int> getAlignedPosition(AlignSeq* ptr, int pos);	// computes the aligned position of a base in an alignment and the 'block' where the base is found
    int getRealPosition(AlignSeq* ptr, int pos, int idx);	// computes the real position of a base dependent on its position in the alignment
    void createOrthoExons(vector<int> offsets);	// searches for the orthologue exons of the exon candidates of the reference species
    list<ExonCandidate*>* getExonCands(int speciesIdx);
    list<OrthoExon> getOrthoExons();
    void openOutputFiles();
    void printExonCands(vector<int> offsets);
    void printOrthoExons(vector<int> offsets);
    void closeOutputFiles();
};

#endif  // _GENEMSA

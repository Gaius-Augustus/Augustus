/**********************************************************************
 * file:    genomicMSA.hh
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
#include "phylotree.hh"
#include "randseqaccess.hh"


struct AlignmentBlock {
    vector<AlignSeq*> alignSpeciesTuple;
};

class GeneMSA {
public:
    static PhyloTree *tree;
    static int utr_range;
    static int orthoExonID; // stores an ID for exons of different species which are orthologous to each other
    static int geneRangeID; // stores an ID for the possible gene ranges of the different species which belong together
    static vector<int> exonCandID; // stores the IDs for exon candidates of different species
    static ofstream *pamlFile;
    static vector< ofstream* > exonCands_outfiles, orthoExons_outfiles, geneRanges_outfiles; // pointers to the output files
    list<AlignmentBlock*> alignment;		// list of the alignment parts which possibly belong to a gene segment
    vector< list<ExonCandidate*>* > exoncands;		// exon candidates found in the different species in a gene segment
    vector< map<string, ExonCandidate*>* > existingCandidates;		// stores the keys of the exon candidates for the different species
    list<OrthoExon> orthoExonsList;		// list of orthologue exons found in a gene segment

    GeneMSA() {};
    ~GeneMSA(){
        if (!alignment.empty()) {
            alignment.clear();
        }
        for (int i=0; i<exoncands.size(); i++) {
            if (exoncands[i]!=NULL) {
                for(list<ExonCandidate*>::iterator it = exoncands[i]->begin(); it != exoncands[i]->end(); it++){
                    delete *it;
                }
                exoncands.at(i)->clear();
                delete exoncands[i];
            }
        }
        for (int i=0; i<existingCandidates.size(); i++) {
            if (existingCandidates[i]!=NULL) {
                existingCandidates.at(i)->clear();
                delete existingCandidates[i];
            }
        }
    };

    string getName(int speciesIdx);
    long int getSeqIDLength(int speciesIdx);
    string getSeqID(int speciesIdx);
    Strand getStrand(int speciesIdx);
    int getStart(int speciesIdx);
    int getEnd(int speciesIdx);
    int getGFF3FrameForExon(ExonCandidate *ec);

    string reverseString(string text);
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
    void createExonCands(const char *dna, double assmotifqthresh=0.15, double assqthresh=0.3, double dssqthresh=0.7);

    pair<int,int> getAlignedPosition(AlignSeq* ptr, int pos);	// computes the aligned position of a base in an alignment and the 'block' where the base is found
    int getRealPosition(AlignSeq* ptr, int pos, int idx);	// computes the real position of a base dependent on its position in the alignment
    void createOrthoExons(RandSeqAccess *rsa, vector<int> offsets);	// searches for the orthologue exons of the exon candidates of the reference species
    list<ExonCandidate*>* getExonCands(int speciesIdx);
    list<OrthoExon> getOrthoExons();
    vector<ExonCandidate*> cutIncompleteCodons(vector<ExonCandidate*> orthoex);
    string getAlignedOrthoExon(AlignSeq *as_ptr, ExonCandidate* ec, string seq, vector<int> offsets, int speciesIdx);
    static void openOutputFiles();
    void printExonsForPamlInput(RandSeqAccess *rsa, vector<int> offsets);
    void printGeneRanges();
    void printExonCands(vector<int> offsets);
    void printOrthoExons(vector<int> offsets);
    void printSingleOrthoExon(OrthoExon &oe, vector<int> offsets);
    static void closeOutputFiles();
};

#endif  // _GENEMSA

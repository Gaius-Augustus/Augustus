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
#include "contTimeMC.hh"

//forward declarations
class OrthoGraph;

class AlignmentBlock {
public:
    AlignmentBlock(size_t n) : rows(n, NULL){} // initialize with NULL, which stand for missing AlignSeqs
    ~AlignmentBlock(){
	// Steffi: this causes a segmentation fault for more than two species. I don't know why.
	// for (int i=0; i<rows.size(); i++) 
	    //delete rows.at(i);	
    }
    friend bool mergeable (AlignmentBlock *b1, AlignmentBlock *b2, int maxGapLen, float mergeableFrac);
public: // should rather be private
    vector<AlignSeq*> rows;
};

/*
 * b1 and b2 can be merged in that order because they are very similar and right next to each other.
 * In at least 'mergeableFrac' of the alignment block rows the aligned sequenes are
 * present, refer to the same terget sequence, are on the same strand and satisfy 0 <= gaplen <= maxGapLen
 */
bool mergeable (AlignmentBlock *b1, AlignmentBlock *b2, int maxGapLen, float mergeableFrac);


class GeneMSA {
public:
    static int utr_range;
    static int orthoExonID; // stores an ID for exons of different species which are orthologous to each other
    static int geneRangeID; // stores an ID for the possible gene ranges of the different species which belong together
    static vector<int> exonCandID; // stores the IDs for exon candidates of different species
    static ofstream *pamlFile;
    static vector< ofstream* > exonCands_outfiles, orthoExons_outfiles, geneRanges_outfiles, omega_outfiles; // pointers to the output files
    list<AlignmentBlock*> alignment;            // list of the alignment parts which possibly belong to a gene segment
    vector< list<ExonCandidate*>* > exoncands;  // exon candidates found in the different species in a gene segment
    vector< map<string, ExonCandidate*>* > existingCandidates; // stores the keys of the exon candidates for the different species
    list<OrthoExon> orthoExonsList;		// list of ortholog exons found in a gene segment
    list<OrthoExon> orthoExonsWithOmega;        // list of ortholog exons with a computed omega=dN/dS ratio

    GeneMSA() {};
    static void setTree(PhyloTree *t){tree = t;}
    static void setCodonEvo(CodonEvo *c){codonevo = c;}
    ~GeneMSA(){
        if (!alignment.empty()) {
	    // for (list<AlignmentBlock*>::iterator it = alignment.begin(); it != alignment.end(); it++) {
	    // 	delete *it;
	    // }
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
     * assmotifqthresh, assqthresh, dssqthresh are between 0 and 1 and thresholds for the inclusion of
     * acceptor/donor splice sites based on the pattern probability
     * assqthresh=0.05 means that only acceptor ss are considered
     * that have a pattern, such that 5% of true splice site patterns have lower probability.
     * The default threshold of 0 means that all splice site patterns are considered.
     */
    void createExonCands(const char *dna, double assmotifqthresh=0.15, double assqthresh=0.3, double dssqthresh=0.7); // get all exon candidates
    Double computeSpliceSiteScore(Double exonScore, Double minProb, Double maxProb); //computes the score for the splice sites of an exon candidate
    pair<int,int> getAlignedPosition(AlignSeq* ptr, int pos);	// computes the aligned position of a base in an alignment and the 'block' where the base is found
    int getRealPosition(AlignSeq* ptr, int pos, int idx);	// computes the real position of a base dependent on its position in the alignment
    void createOrthoExons(vector<int> offsets, OrthoGraph &orthograph);	// searches for the orthologue exons of the exon candidates of the reference species
    list<ExonCandidate*>* getExonCands(int speciesIdx);
    list<OrthoExon> getOrthoExons();
    void cutIncompleteCodons(vector<ExonCandidate*> &orthoex);
    void readOmega(string file);
    string getAlignedOrthoExon(AlignSeq *as_ptr, ExonCandidate* ec, string seq, int offset);
    vector <string> getSeqForPaml(AlignmentBlock *it_ab, vector<ExonCandidate*> oe, vector<string> seq, vector<int> offsets, vector<int> speciesIdx);
    static void openOutputFiles();
    void printGeneRanges();
    void printExonCands(vector<int> offsets);
    void printOrthoExons(RandSeqAccess *rsa, vector<int> offsets);
    void printSingleOrthoExon(OrthoExon &oe, vector<int> offsets, bool files = true, double omega=-1, int numSub=-1);
    void printExonWithOmega(vector<int> offsets);
    void printExonsForPamlInput(RandSeqAccess *rsa,  OrthoExon &oe,  vector<int> offsets);
    static void closeOutputFiles();
private:
    static PhyloTree *tree;
    static CodonEvo *codonevo;
};

#endif  // _GENEMSA

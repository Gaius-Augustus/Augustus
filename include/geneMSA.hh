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
#include "alignment.hh"
#include "exoncand.hh"
#include "orthoexon.hh"
#include "phylotree.hh"
#include "randseqaccess.hh"
#include "contTimeMC.hh"

//forward declarations
class OrthoGraph;

class GeneMSA {
public:
    GeneMSA(RandSeqAccess *rsa, Alignment *a);
    ~GeneMSA();
    static int padding; // add this many bases to the region before and after the aligned region
    static int orthoExonID; // stores an ID for exons of different species which are orthologous to each other
    static int geneRangeID; // stores an ID for the possible gene ranges of the different species which belong together
    static vector<int> exonCandID; // stores the IDs for exon candidates of different species
    // pointers to the output files
    static vector< ofstream* > exonCands_outfiles, orthoExons_outfiles, geneRanges_outfiles, omega_outfiles; 
    Alignment* alignment;            // alignment of regions which possibly belong to a gene
    vector< list<ExonCandidate*>* > exoncands;  // exon candidates found in the different species in a gene segment
    list<OrthoExon> orthoExonsList;		// list of ortholog exons found in a gene segment

    static void setTree(PhyloTree *t){tree = t;}
    static void setCodonEvo(CodonEvo *c){ codonevo = c; }
    static int numSpecies(){ return tree->numSpecies(); }

    string getSeqID(int speciesIdx) {
	if (alignment->rows[speciesIdx])
	    return alignment->rows[speciesIdx]->seqID;
	else 
	    return "";
    }
    Strand getStrand(int speciesIdx){
	if (alignment->rows[speciesIdx])
	    return alignment->rows[speciesIdx]->strand;
	else 
	    return STRAND_UNKNOWN;
    }
    // start and end position (0-based, inclusive) of the gene range for the given species
    // start <= end refer to the FORWARD strand of the sequence (different from class 'alignment')
    // start = end = -1 if the species is absent from the alignment
    int getStart(int speciesIdx){ return starts[speciesIdx]; }
    int getEnd(int speciesIdx){ return ends[speciesIdx]; }
    int getGFF3FrameForExon(ExonCandidate *ec);

    string reverseString(string text);
    map<string,ExonCandidate*>* getECHash(list<ExonCandidate*> *ec); // adds the keys to the map function

    /*
     * assmotifqthresh, assqthresh, dssqthresh are between 0 and 1 and thresholds for the inclusion of
     * acceptor/donor splice sites based on the pattern probability
     * assqthresh=0.05 means that only acceptor ss are considered
     * that have a pattern, such that 5% of true splice site patterns have lower probability.
     * The default threshold of 0 means that all splice site patterns are considered.
     */
    void createExonCands(int s, const char *dna, double assmotifqthresh=0.15, double assqthresh=0.3, double dssqthresh=0.7);

    //computes the score for the splice sites of an exon candidate
    Double computeSpliceSiteScore(Double exonScore, Double minProb, Double maxProb); 
    // computes the aligned position of a base in an alignment and the 'block' where the base is found
    pair<int,int> getAlignedPosition(AlignmentRow* row, int pos);  

    /**
     * find all ortholog exon candidates, that are present in at least max(2, consThresh * m)
     * where m <= numSpecies is the number of species that are present in the alignment
     * ortholog exon candidates:
     * - both splice sites align exactly
     * - the exon candidate types agrees (single, rsingle, internal0, ...)
     * - the phases at both boundaries agree (i.e. exon candidate types and length modulo 3)
     * EC coordinates are region-based, as they are used in the OrthoGraph
     */
    void createOrthoExons(float consThres = 0.5);
    list<ExonCandidate*>* getExonCands(int speciesIdx){ return exoncands.at(speciesIdx); }
    list<OrthoExon> getOrthoExons() { return orthoExonsList; }
    string getAlignedOrthoExon(AlignmentRow *as_ptr, ExonCandidate* ec, string seq, int offset);
    void printStats(); // to stdout
    void printGeneRanges();
    void printExonCands();
    void printOrthoExons(RandSeqAccess *rsa);
    void printSingleOrthoExon(OrthoExon &oe, bool files = true, double omega=-1, int numSub=-1);
    void computeOmegas(vector<AnnoSequence> const &seqRanges);
    static void openOutputFiles();
    static void closeOutputFiles();

private:
    vector<string> getCodonAlignment(OrthoExon const &oe, vector<AnnoSequence> const &seqRanges);
    void cutIncompleteCodons(OrthoExon &oe);
    static PhyloTree *tree;
    static CodonEvo *codonevo;
    vector<int> starts, ends; // gene ranges for each species
    vector<int> offsets; // this many bases are upstream from the region
    RandSeqAccess *rsa;
};

#endif  // _GENEMSA

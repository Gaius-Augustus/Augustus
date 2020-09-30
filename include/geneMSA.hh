/*
 * geneMSA.hh
 *
 * License: Artistic License, see file LICENSE.TXT or 
 *          https://opensource.org/licenses/artistic-license-1.0
 */

#ifndef _GENEMSA
#define _GENEMSA

// project includes
#include "alignment.hh"
#include "exoncand.hh"
#include "orthoexon.hh"
#include "randseqaccess.hh"
#include "speciesgraph.hh"
#include "contTimeMC.hh"
#include<unordered_map>
#include<boost/functional/hash.hpp>

//forward declarations

/**
 * 
 * @author Mario Stanke
 * @author Alexander Gebauer
 */
class OrthoGraph;
string printRFC(vector<int>);

/* store cumulative sums of log likelihoods for each omega (optional: cumulative sum of number of substitutions) 
 * we create an object of cumValues for each bit_vector and reading frame combination  
 */
struct cumValues{
  vector<double> logliks;
  int numSubs;
  int id;
  cumValues(int i, int s=-1):numSubs(s), id(i){};
  void addLogliks(vector<double>* ll){
    if(logliks.size() == 0){
      logliks.resize(ll->size(),0.0);
      if(logliks.size() == 0){
	cerr<<"logliks still empty!"<<endl;
      }
    }
    for(int u = 0; u < ll->size(); u++){
      logliks[u] += (*ll)[u];
    }
  }
  void addNumSubs(int subs){
    numSubs += subs;
   }
};

/**
 * @brief multiple sequence alignment of genomes for comparative gene prediction
 * 
 * @author Mario Stanke
 * @author Alexander Gebauer
 */
class GeneMSA {
public:
    GeneMSA(RandSeqAccess *rsa, Alignment *a);
    ~GeneMSA();

    /*
     *  get and set functions
     */
    // start and end position (0-based, inclusive) of the gene range for the given species
    // start <= end refer to the FORWARD strand of the sequence (different from class 'alignment')
    // start = end = -1 if the species is absent from the alignment
    string getSeqID(int speciesIdx);
    Strand getStrand(int speciesIdx);
    int getStart(int speciesIdx){ return starts[speciesIdx]; }
    int getEnd(int speciesIdx){ return ends[speciesIdx]; }
    list<ExonCandidate*>* getExonCands(int speciesIdx){ return exoncands.at(speciesIdx); }
    //list<OrthoExon> getOrthoExons() { return orthoExonsList; }
    Alignment * getAlignment() {return alignment;}
    vector<int> getOffsets() {return offsets;}
 
    map<string,ExonCandidate*>* getECHash(list<ExonCandidate*> *ec); // adds the keys to the map function

    /*
     * assmotifqthresh, assqthresh, dssqthresh are between 0 and 1 and thresholds for the inclusion of
     * acceptor/donor splice sites based on the pattern probability
     * assqthresh=0.05 means that only acceptor ss are considered
     * that have a pattern, such that 5% of true splice site patterns have lower probability.
     * The default threshold of 0 means that all splice site patterns are considered.
     */
    void createExonCands(int s, const char *dna, map<int_fast64_t, ExonCandidate*> &ecs, map<int_fast64_t, ExonCandidate*> &addECs);
    void setExonCands(vector<map<int_fast64_t, ExonCandidate*> > &ecs);

    /**
     * find all ortholog exon candidates, that are present in at least max(2, consThresh * m)
     * where m <= numSpecies is the number of species that are present in the alignment
     * Only report OrthoExons oe with at least 'minAvLen' as average length of the exon candidates in oe.
     * ortholog exon candidates:
     * - both splice sites align exactly
     * - the exon candidate types agrees (single, rsingle, internal0, ...)
     * - the phases at both boundaries agree (i.e. exon candidate types and length modulo 3)
     * EC coordinates are region-based, as they are used in the OrthoGraph
     */
    void createOrthoExons(list<OrthoExon> &orthoExonsList, map<int_fast64_t, list<pair<int,ExonCandidate*> > > &alignedECs, Evo *evo, float consThres = 0.0, int minAvLen = 0);
    void printStats(); // to stdout
    void printGeneRanges();
    void printExonCands();
    void printOrthoExons(list<OrthoExon> &orthoExonsList);
    void computeOmegas(list<OrthoExon> &orthoExonsList, vector<AnnoSequence*> const &seqRanges, PhyloTree *ctree);
    void computeOmegasEff(list<OrthoExon> &orthoExonsList, vector<AnnoSequence*> const &seqRanges, PhyloTree *ctree, ofstream *codonAli);
    vector<string> pruneToBV(vector<string> *cs, bit_vector bv); // prune codon strings to bit_vector
    vector<int> pruneToBV(vector<int> *rfc, bit_vector bv); // prune RFC to bit_vector
    double omegaForCodonTuple(vector<double> *loglik);
    void printOmegaForCodon(string outdir);
    void printCumOmega();
    void comparativeSignalScoring(list<OrthoExon> &orthoExonsList);
    // Charlotte Janas playground
    LocusTree *constructTree(); // creates, stores are returns the locus tree for the sequence tuple

    // calculate a columnwise conservation score and output it (for each species) in wiggle format
    void calcConsScore(list<OrthoExon> &orthoExonsList, vector<AnnoSequence*> const &seqRanges, string outdir);
    double calcColumnScore(int a, int c, int t, int g); // input: number of a,c,t and g's in one alignment column 
    void consToWig(vector<double> &consScore, string outdir);

    // static functions
    static void setTree(PhyloTree *t){tree = t;}
    static void setCodonEvo(CodonEvo *c){ codonevo = c; }
    static int numSpecies(){ return tree->numSpecies(); }
    static void openOutputFiles(string outdir);
    static void closeOutputFiles();

    // static data members
    static int padding; // add this many bases to the region before and after the aligned region
    static int orthoExonID; // ID for exons of different species which are orthologous to each other
    static int geneRangeID; // stores an ID for the possible gene ranges of the different species which belong together
    static vector<int> exonCandID; // IDs for exon candidates of different species
    static unordered_map< bit_vector, PhyloTree*, boost::hash<bit_vector>> topologies;
    // pointers to the output files
    static vector< ofstream* > exonCands_outfiles, orthoExons_outfiles, geneRanges_outfiles_bed, geneRanges_outfiles_gff, omega_outfiles; 
    // map that stored all codon combinations on which fitch and pruning algorithm already have been run (for calculation of omega and number of substitutions)
    static map<vector<string>, pair<vector<double>, int> > computedCumValues;

    void printSingleOrthoExon(OrthoExon &oe, bool files = true);
    void collect_features(int species, list<OrthoExon> *hects, SpeciesGraph *speciesgraph);

    /**
     * getAllOEMsas obtains and prints multiple sequence alignments (MSAs) and their label y=0,1, whether
     * it constitutes a real CDS or not in the reference species.
     *
     * @param[in] species is the reference species index
     * @param[in] hects ortho exons to get alignments from
     * @param[in] ref_class assignment of CDS keys as 1 (= in reference anno) or 0.
     * @param[in] seqRanges the aligned sequences
     */
    void getAllOEMsas(int species, list<OrthoExon> *hects, unordered_map<string,int> *ref_class,
                         vector<AnnoSequence*> const &seqRanges);

    /**
     * getMsa obtains the alignment in string format, including
     * unaligned insertions that are not part of any fragment.
     *
     * @return msa vector of alignment rows
     * @param[in] oe the OrthoExon whose alignment is sought
     * @param[in] seqRanges contains the alignment of the larger region tuple
     * @param[in] flanking allows to add padding on both sides, e.g. for sequence signals
     */
    StringAlignment getMsa(OrthoExon const &oe, vector<AnnoSequence*> const &seqRanges, size_t flanking = 0);
private:
  vector<string> getCodonAlignment(OrthoExon const &oe, vector<AnnoSequence*> const &seqRanges, const vector<vector<fragment>::const_iterator > &froms, map<unsigned, vector<int> > *alignedCodons = NULL, bool generateString=true, vector<vector<int> > *posStoredCodons = NULL, ofstream *codonAli = NULL);
    void cutIncompleteCodons(OrthoExon &oe);
    cumValues* findCumValues(bit_vector bv, vector<int> rfc);
    static PhyloTree *tree;
    LocusTree *ltree;
    static CodonEvo *codonevo;
    vector<int> starts, ends; // gene ranges for each species
    vector<int> offsets; // this many bases are upstream from the region
    RandSeqAccess *rsa;
    Alignment* alignment;            // alignment of regions which possibly belong to a gene
    vector< list<ExonCandidate*>* > exoncands;  // exon candidates found in the different species in a gene segment
    //list<OrthoExon> orthoExonsList;		// Steffi: due to multiple copying, I remove this as attribute of any class. Instead it can be passed from compgenepred.cc by reference, whenever necessary.
    unordered_map<bit_vector, vector<pair<vector<int>, cumValues> >, boost::hash<bit_vector> > cumOmega; // stores cumulative omega values for every reading frame combination and every bitvector that exist 
    map<bit_vector, map<vector<int>, vector<double> > > codonOmega;

    #ifdef TESTING
    /*
    *   added by Giovanna Migliorelli 14.05.2020 
    *   code responsible for serialization
    *   aknowledgement : https://www.boost.org/doc/libs/1_70_0/libs/serialization/doc/tutorial.html
    */

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & starts;
        ar & alignment;

    }

    /*  added by Giovanna Migliorelli 14.05.2020
    *   default consttruct to allow for serialization
    */

    GeneMSA(){}; 

    #endif
};



#endif  // _GENEMSA

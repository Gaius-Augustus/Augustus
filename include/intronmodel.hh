/*****************************************************************************\
 * Filename : intronmodel.hh
 * Authors   : Emmanouil Stafilarakis, Mario Stanke
 *
 * Description: Intron Model Header File
 *
 *
 * Date       |   Author              |  Changes
 *------------|-----------------------|----------------------------------------
 * 10.11.2001 | Stafilarakis Emm.     | Creation of the file
 * 31.05.2002 | Mario Stanke          | static data members
 \******************************************************************************/

#ifndef _INTRONMODEL_HH
#define _INTRONMODEL_HH

#include "statemodel.hh"


/**
 * @memo    The intron model class.
 *
 * @doc     
 *
 * @authors  Emmanouil Stafilarakis, Mario Stanke
 */
 
class IntronModel : public StateModel {
public:
    IntronModel();
    ~IntronModel();
    
    StateType getStateType () const {
	return itype;
    }
    
    /**
     * Build all needed probabilities from the given
     * gene set.
     *
     * @param   annoseq A single linked list annotated sequences
     */
    void buildModel( const AnnoSequence* annoseq, int parIndex );
    void registerPars( Parameters* parameters);
    /**
     * Print the wanted probabilities in files.
     * The names of the files should be given with
     * the Properties object /IntronModel/outfile
     */
    void printProbabilities( int parIndex, BaseCount *bc, const char* suffix = NULL );
    void initAlgorithms(Matrix<Double>&, int);
    virtual void updateToLocalGCEach( Matrix<Double>& trans, int cur);
    void viterbiForwardAndSampling(ViterbiMatrixType&, ViterbiMatrixType&, int, int, 
				   AlgorithmVariant, OptionListItem&);
    Double emiProbUnderModel  (int begin, int end) const;
    Double seqProb            (int left, int right) const;
    static Double dSSProb     (int base, bool forwardStrand);
    static Double aSSProb     (int base, bool forwardStrand);
    static void init();
    static void resetPars() {
	initSnippetProbs();
	initAlgorithmsCalled = false;
    }
    static void updateToLocalGC(int from = -1, int to = -1);
    static void readProbabilities(int parIndex);
    static void readAllParameters();
    static void updateParameters(int idx);
    static void storeGCPars(int idx);
    static Integer getD() {return d;}
    static double getMAL() {return mal.doubleValue();}
    static void resetModelCount(){introncount = 0;};
    static double getGeoProb(){return geoProb;}
    static Double getAssMotifProbThreshold(double q) {return assMotif->getProbThreshold(q);}
    static double getMeanIntrLen();
private:
    static void initSnippetProbs();
    void processSequence( const char* start, const char* end);
    /**
     *
     */
    void buildProbabilities ( const AnnoSequence* annoseq );
    void buildLenDist       ( const AnnoSequence* annoseq );
    void processASS         ( const char* dna, int pos, Boolean withMotif=true);
    void processDSS         ( const char* dna, int pos );
    void makeDSSProbs       ( );
    void storeIntronLengths ( const AnnoSequence* annoseq);
    void printLengthQuantiles();
    void initCountVars      ( );
    void readSpliceSites    ();

public:    
    static Integer         k;
    static PatMMGroup      emiprobs;
    static PatMMGroup      *GCemiprobs;
    static BinnedMMGroup   dssBinProbs;
    static BinnedMMGroup   assBinProbs;
    static vector<Double>  lenDist;
private:
    StateType              itype;
    Integer                gweight;
    char*                  codon;
    static int beginOfBioIntron, endOfBioIntron;
    static vector<Integer> emicount;
    static vector<Integer> intlencount;
    static Integer         introns;         // number of introns
    static Integer         introns_d;       // number of introns of length <= d
    static Integer         d;
    static Integer         minwindowcount;
    static double          slope_of_bandwidth;
    static Double          patpseudo;
    static Double          probShortIntron; // probability that intron is at most d long
    static Double          *GCprobShortIntron; // array for each GC content class
    static Double          mal;             // mean additional length (if longer than d)
    static Double          *GCmal;             // array for each GC content class
    static Integer         introncount;
    static Integer         gesbasen;

    // variables related to splicing
    static Integer         c_ass;
    static Integer         c_dss;
    static vector<Integer> asscount;
    static vector<Integer> dsscount;
    static Integer         ass_upwindow_size;
    static vector<Double>  assprobs;
    static vector<Double>  dssprobs;
    static Motif           *assMotif;    // basecounts of the window before the ass
    static Motif           *GCassMotif;  // array of Motifs, one for each GC content class
    static Boolean         hasSpliceSites;
    static Double          asspseudo;         // pseudocount for patterns in acceptor splice sites
    static Double          dsspseudo;         // pseudocount for patterns in donor splice sites 
    static Double          dssneighborfactor; // taken from the prob of the neighbor patterns
    static Integer         ass_motif_memory;  // order of the markov chain in the ass upstream motif
    static Integer         ass_motif_radius;  // radius of the pooling window in the ass upstream motif
    static double          non_gt_dss_prob;
    static double          non_ag_ass_prob;
    static SnippetProbs    *snippetProbs, *rSnippetProbs;  
    static bool            initAlgorithmsCalled, haveSnippetProbs;
    static int             lastParIndex; // GC-index of current parameter set
    static Integer         verbosity;
    static double          geoProb;
    static int             ass_outside; // Constant::ass_upwindow_size + Constant::ass_start + ASS_MIDDLE;
};

class IntronModelError : public ProjectError {
public:
    IntronModelError(string msg) : ProjectError(msg) {}
};

#endif    //  _INTRONMODEL_HH

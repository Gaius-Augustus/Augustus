/**********************************************************************
 * file:    namgene.hh
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  the class NAMGene is AUGUSTUS' entry point for the algorithms
 * authors: Mario Stanke (mario@gobics.de), Stafilarakis
 *
 *********************************************************************/

#ifndef _NAM_GENE_HH
#define _NAM_GENE_HH
#include <climits>

// project includes
#include "extrinsicinfo.hh"  // SequenceFeatureCollection
#include "vitmatrix.hh"
#include "matrix.hh"
#include "pp_scoring.hh"
#include "statemodel.hh"


class NAMGeneError : public ProjectError {
public:
    NAMGeneError( string msg ) : ProjectError( msg ) { }
};


class NAMGene {
public:
    NAMGene();
    ~NAMGene() {}

    void readModelProbabilities( int number = 1);
    StatePath* getSampledPath(const char *dna, const char* seqname = NULL);
    StatePath* getViterbiPath(const char *dna, const char* seqname = NULL);
    StatePath* getTrainViterbiPath(const char *dna, SequenceFeatureCollection *sfc);
    /**
     * Predict the path through the HMM and compute the predicted genes.
     * strand: plusstrand, minusstrand or bothstrands
     */
    Transcript* doViterbiPiecewise(SequenceFeatureCollection& sfc, AnnoSequence *annoseq, Strand strand);

    Double getEmissionProbability();

    /**
     * Get the matrix with the variables calculated by the forward algorithm.
     * @return  The forward variables.
     * @see forwardAlgorithm
     */
    const ViterbiMatrixType& getForwardVariables( ) {
	return forward;
    } 
    
    /**
     * Get the matrix with the variables calculated by the viterbi algorithm.
     * @return  The viterbi variables.
     * @see viterbiAlgorithm
     */
    const ViterbiMatrixType& getViterbiVariables( ) {
	return viterbi;
    }

    /*
     * get the type corresponding to a state number
     * @param   i the thumber of a state, as in the properties
     * @return  the enumeration type
     */
    StateType getStateType(int i);
    
    int getStateIndex(StateType type);
    Double getPathEmiProb(StatePath *path, const char *dna, SequenceFeatureCollection& sfc, int countStart=-1, int countEnd=-1);
    void setNeedForwardTable(bool b){needForwardTable = b;}
    // set the path and emiProbs for all annotations in annoseq (for option scoreTx)
    void setPathAndProb(AnnoSequence *annoseq, FeatureCollection &extrinsicFeatures);
    void setAllTranscripts(list<Transcript*> *tl) {sampledTxs = tl;}
    list<Transcript*> *getAllTranscripts() {return sampledTxs;}
    void getPrepareModels(const char *dna, int len) {prepareModels(dna, len);}
private:
    /**
     * Start the viterbi algorithm with the given DNA sequence.
     *
     * @param   dna The DNA sequence to be used in the algorithm.
     * @return   the viterbi path
     */
    void viterbiAndForward(const char* dna, bool useProfile=false);
    
    /*
     * repeatedly do the viterbi algorithm on pieces of dna, no strands
     */
    list<AltGene>* getStepGenes(AnnoSequence *annoseq, SequenceFeatureCollection& sfc, Strand strand, bool onlyViterbi=true);

    list<AltGene>* findGenes(const char *dna, Strand strand, bool onlyViterbi=true);
    int getNextCutEndPoint(const char *dna, int beginPos, int maxstep, SequenceFeatureCollection& sfc);
    void readTransAndInitProbs( );
    void readOvlpLenDist( );
    void checkProbsConsistency( );
    void computeReachableStates( );
    void createStateModels( );
    void setStatesInitialProbs();
    void initAlgorithms(); // called for new DNA
    void updateToLocalGCEach(int idx, int from = -1, int to = -1);
    void printDPMatrix(); // print Viterbi matrix
    /**
     * Read in the parameters for all the models from files
     *
     */
    void prepareModels(const char*dna, int len);

    int tryFindCutEndPoint(StatePath *condensedExamPath, int examIntervalStart, int examIntervalEnd, list<Feature> *groupGaps, bool onlyInternalIR);
private:
    /// @doc The forward variables matrix.
    ViterbiMatrixType      forward;
    /// @doc The viterbi variables matrix.
    ViterbiMatrixType      viterbi; 
    /// @doc The transitions matrix.
    Matrix<Double>      transitions;
    /// @doc The initial probabilities
    vector<Double>      initProbs;
    /// @doc The initial probabilities
    vector<Double>      termProbs;
    // Whether the state is reachable at all with prob > 0
    vector<Boolean> stateReachable;
    /// @doc The array with the model states.
    vector<StateModel*> states;  
    /// @doc The count of the model states
    Integer             statecount;
    /// @doc mapping from state numbers to state types
    vector<StateType>   stateMap;
    /// @doc index of the parameters read into the models (based on sequence composition)
    int                 lastParIndex;
    /// overlap length distributions for bacterial model
    Boolean proteinOutput;
    Boolean codSeqOutput;
    Boolean noInFrameStop;
    double minmeanexonintronprob;
    double minexonintronprob;  // lower bound on probabilities of all exons and introns in the coding region
    int maxtracks;       // maximum reported number of transcripts at the same sequence position
    int sampleiterations;
    bool alternatives_from_sampling;
    bool alternatives_from_evidence;
    bool mea_prediction;
    bool mea_eval;
    bool needForwardTable;
    bool show_progress;
    PP::SubstateModel* profileModel;
    ContentStairs cs; // holds the local GC content class per position in the currently examined DNA
    int curGCIdx; // current index of GC content class
    list<Transcript*> *sampledTxs; // = alltranscripts stored for MultSpeciesMode
};


#endif  // _NAM_GENE_HH

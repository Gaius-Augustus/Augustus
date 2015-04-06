/***********************************************************************
 * file:    statemodel.hh
 * licence: Artistic Licence, see file LICENCE.TXT or
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  base interface class for the state model classes
 * authors: Emmanouil Stafilarakis, Mario Stanke (mario@gobics.de), 
 *          Oliver Keller
 *
 **********************************************************************/

#ifndef _STATEMODEL_HH
#define _STATEMODEL_HH

// project includes
#include "matrix.hh"
#include "vitmatrix.hh"
#include "merkmal.hh"



/*
 * constants used in splitting into state/substate pairs
 * we use the lower byte for the state, and the upper 3 for the substate
 * SHIFT_LEFT is used to build the full id from the substate id
 * SHIFT_RIGHT is used to retrieve back the substate id
 */
#define SHIFT_SIZE (sizeof (signed char)*8)
#define SHIFT_LEFT(x) ((x) << SHIFT_SIZE)
#define SHIFT_RIGHT(x) ((x) >> SHIFT_SIZE)
#define MAX_STATECOUNT SHIFT_LEFT(1)


/*
 * getFullStateId:  merge state/substate pairs into single integers (used in calls of viterbi)
 *                  precondition: state >= -1, substate >= -1
 *                  if substate is -1 (none), state remains unchanged
 * getStatePair:    split integers back into state/substate pairs
 *
 */
inline int getFullStateId(int state, SubstateId substate = SubstateId()) { 
    return SHIFT_LEFT (substate.fullId()) + state;
}
inline void getStatePair(int fullState, int& state, SubstateId& substate) {
    state = (signed char)(fullState);
    substate.set(SHIFT_RIGHT(fullState -state));
}

/*
 * Predecessor in the state transition Graph
 */
struct Ancestor{
    Ancestor(int newpos=0, Double newval=0.0) : pos(newpos), val(newval) {} 
    int     pos;
    Double  val;
};


/*
 * This is the base interface common to all state model classes
 * (ExonModel, IntronModel, IGenicModel, UTRModel)
 *
 */
class StateModel {
protected:
    StateModel() {}  // do not create StateModel object
public:
    void initPredecessors(Matrix<Double>&, int self);

    // virtual methods to be implemented by the specialised classes
    virtual void registerPars(Parameters* parameters) {}
    virtual void buildModel(const AnnoSequence* annoseq, int parIndex) =0;
    virtual void printProbabilities(int zusNumber=1, BaseCount *bc = NULL, const char* suffix = NULL) =0;
    virtual StateType getStateType() const = 0;
    virtual void viterbiForwardAndSampling(ViterbiMatrixType&, ViterbiMatrixType&, int, int,
					   AlgorithmVariant, OptionListItem&) = 0;
    virtual Double emiProbUnderModel(int , int) const = 0;
    virtual void initAlgorithms(Matrix<Double>&, int) = 0;
    virtual void updateToLocalGCEach(Matrix<Double>& trans, int cur) { };
    virtual ~StateModel() {} // nothing to do here since class is purely abstract

    // class functions
    static void init();
    static StateModel* newStateModelPtr (const char* path);
    static void determineShortPatterns(const vector<Integer>& patcounts, int k, int minCount);
    static void makeProbsFromCounts(vector<Double > &patprobs , const vector<Integer > &patcounts, 
				    int k, Double pseudocount, Boolean shorten = false);
    static void computeEmiFromPat(const vector<Double>& patprobs, vector<Double>& emiprobs, Integer k);
    static void prepareViterbi(const char* dna, int len, const vector<StateType> &stateMap);
    static void readProbabilities(int);
    static void resetPars();
    // this is done when gcIdx changes once per state class
    static void updateToLocalGC(int idx, int from = -1, int to = -1);
    static void readAllParameters();
    static void storeGCPars(int);
    static void resetModelCounts();
    static bool isPossibleDSS(int pos) {
	return pos >= 1 && pos <= dnalen-2 &&
	    (onGenDSS(sequence + pos) || 
	     seqFeatColl->isHintedDSS(pos, plusstrand));
    } 
    static bool isPossibleRDSS(int pos) {
	return pos >= 1 && pos <= dnalen-2 &&
	    (onGenRDSS(sequence + pos -1) ||
	     seqFeatColl->isHintedDSS(pos, minusstrand));
    } 
    static bool isPossibleASS(int pos) {
	return pos >= 1 && pos <= dnalen-2 &&
	    (onASS(sequence + pos -1) ||
	     seqFeatColl->isHintedASS(pos, plusstrand));
    } 
    static bool isPossibleRASS(int pos) {
	return pos >= 1 && pos <= dnalen-2 &&
	    (onRASS(sequence + pos) ||
	     seqFeatColl->isHintedASS(pos, minusstrand));
    } 
    static void setSFC(SequenceFeatureCollection *psfc) {
	seqFeatColl = psfc;
    }
    static void setPP(PP::SubstateModel* mdl) {
	profileModel = mdl;
    }
    static void setCountRegion(int from, int to){countStart = from; countEnd = to;}
    static int getActiveWindowStart(int);
    static void setGCIdx(int idx) {gcIdx = idx;}
    static void setContentStairs(ContentStairs *stairs) {cs = stairs;}
    static int getGCIdx(int at){if (cs) return cs->idx[at]; else return -1;}
protected:
    // variable unique to each model
    vector<Ancestor>  ancestor;    // predecessor in the state transition graph

    // class variables shared by all models
    static const vector<StateType>* stateMap;  // needed in exonmodel 
    static const char* sequence;   // the sequence currently examined
    static int dnalen;
    static SequenceFeatureCollection* seqFeatColl;
    static vector<Boolean>* shortpattern;
    static PP::SubstateModel* profileModel;
    static int countStart, countEnd; // this is the range of positions over which features are counted in CRF training
    static int activeWinLen;         // states ending before the active window will be deleted if they are not yet used
    static int gcIdx; // GC content class index for all states
    static ContentStairs *cs;
}; // class StateModel


/*
 * classes Snippet*
 * intelligently store and retreive the sequence emission probabilities of the sequence from a to b
 * for common pairs of a and b
 */
class SnippetListItem {
public:
    SnippetListItem (Double p, int length) { prob = p; len = length; next = NULL;}
    ~SnippetListItem() {
	if (next)
	    delete next;
    }
    Double prob;
    int len;
    SnippetListItem *next;
};

class SnippetList {
public:
    Double getProb(int base, int &partlen);
    SnippetList() {first = last = NULL;}
    ~SnippetList(){}
    SnippetListItem *first, *last;
};

class SnippetProbs {
public:
    SnippetProbs(const char* dna, int k, bool forwardStrand=true){
	this->dna = dna;
	this->k = k;
	this->forwardStrand = forwardStrand;
	if (dna)
	  n = strlen(dna);
	else 
	  n = 0;
	snippetlist = new SnippetList*[n];
	for (int i=0; i<n; i++)
	    snippetlist[i] = NULL;
    }
    Double getSeqProb(int base, int len);
    ~SnippetProbs () {
	if (snippetlist) {
	    for (int i=0; i<n; i++) {
		if (snippetlist[i]){
		    delete snippetlist[i]->first;
		    delete snippetlist[i];
		}
	    }
	    delete [] snippetlist;
	}
    }
    void setEmiProbs(vector<Double> *emiprobs) {
	this->emiprobs = emiprobs;
    }

    void addProb(int base, int len, Double p);

protected:
    const char *dna;
    int n;
    vector<Double> *emiprobs;
    int k;
    SnippetList **snippetlist;
    bool forwardStrand;
private:
    Double getElemSeqProb(int base, int len);
};

/*
 * SegProbs - another class for caching probabilities of sequence segments
 * idea: cumulative product => get segment prob in constant time via a ratio 
 */
class SegProbs {
public:
    SegProbs(const char* dna, int k, bool forwardStrand=true) : s2i(k+1) {
	this->dna = dna;
	this->k = k;
	this->forwardStrand = forwardStrand;
	if (dna)
	    n = strlen(dna);
	else 
	    n = 0;
	cumProds = new Double[n+1];
    }
    ~SegProbs () { delete [] cumProds; }
    void setEmiProbs(vector<Double> *emiprobs, int from = -1, int to = -1);
    Double getSeqProb(int from, int to);
    int getN() { return n; }
protected:
    const char *dna;
    int n;
    int k;
    Seq2Int s2i;
    vector<Double> *emiprobs;
    Double *cumProds;
    bool forwardStrand;
};

/*
 * data structure to store possible endOfPred positions to iterate directly
 * only over those endOfPred positions, that are possible in the inner loop of
 * the Viterbi algorithm
 */
class EOPList {
public:
   void clear() { possibleEndOfPreds.clear(); }
   void init() {
      eopit = possibleEndOfPreds.begin();
      inCache = false;
   }
   void decrement(int &endOfPred);
   void update(int endOfPred);
   list<int> possibleEndOfPreds;
   list<int>::iterator eopit;
   bool inCache;
};


#endif  //  _STATEMODEL_HH

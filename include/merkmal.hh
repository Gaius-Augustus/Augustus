/*****************************************************************************\
 * Filename : merkmal.cc
 * Author   : Mario Stanke, mstanke@gwdg.de
 *
 * Description: Features (German: Merkmale) for the Conditional Random Field
 *
 * Date       |   Author        |  Changes
 *------------|-----------------|------------------------------------------
 * 27.07.08   | Mario Stanke    | creation of the classes
 \****************************************************************************/

#ifndef _MERKMAL_HH
#define _MERKMAL_HH

// project includes
#include "types.hh"
#include "extrinsicinfo.hh"
#include "gene.hh"

// standard C/C++ includes
#include <iostream>
#include <vector>
#include <iomanip>  // for setw

using namespace std;


// Forward declarations
class Parameters;
class StateModel;


/*
 * MMGroup (MM = Merkmal = Feature)
 * The derived classes contain a group of parameters, e.g. a motif or an emission probability matrix.
 * As an example take the emission probabilities of exons. Here we have one parameter for every reading
 * frame f and every pattern number p.
 * The derived class FramedEmiProb contains the matrix of parameters and implements a function
 * int getIndex(f,p) for intuitive access to the linear vector of parameters.
 */
class MMGroup {
 public:
    MMGroup(string name){parname = name; parameters = NULL;};
    MMGroup(){parameters = NULL;};
    virtual ~MMGroup() {};
    // derived classes implement:
    // Double getFactor(...); e.g. getFactor(int frame, int pattern)
    // int getIndex(...); e.g. getIndex(int frame, int pattern)
    void registerPars(Parameters* parameters, double alphavalue = 1.0); // adds this group of features to the Parameters
    virtual string verbalDescription(int index) = 0;
    virtual Double* getFactor(int index) = 0;
    virtual void smooth(){} // make parameters of this feature group smooth, consistant or monotonic, etc.
    void addCount(int index, int count=1);
    int getOffset() {return offset;}
    void setOffset(int offset) { this->offset = offset; }
    void setName(string name) {parname = name;}
    virtual int getNumPars() {return numPars;}
protected:
    string parname;                        // name for the parameters group, e.g. "exon emiprobs" or "ass motif"
    int numPars;                           // number of parameters
private:
    Parameters* parameters;                // Parameters object that this group belongs to
    int offset;                            // offset in the Parameters vector where the features of this group start
};

/*
 * PatMMGroup
 * Contains a vector of parameters. Is used in particular for intron emiprobs.
 */
class PatMMGroup : public MMGroup{
public:
    int getNumPars() {numPars = probs.size(); return numPars;}
    PatMMGroup(string name) : MMGroup(name){};
    PatMMGroup() : MMGroup(){};
    Double* getFactor(int index);
    string verbalDescription(int index);
    Double getMinProb(float qthresh);
    // data members
    int order;
    vector<Double> probs;
};


/*
 * FramedPatMMGroup
 * Contains a vector of parameters for each frame. Is used in particular for exon emiprobs.
 */
class FramedPatMMGroup : public MMGroup{
public:
    ~FramedPatMMGroup(){}
    int getNumPars() {numPars = probs[0].size() + probs[1].size() + probs[2].size(); return numPars;}
    FramedPatMMGroup() : MMGroup(){};
    FramedPatMMGroup(string name) : MMGroup(name){};
    void getFramePat(int index, int &frame, int &pattern);
    inline int getIndex(int frame, int pattern){ return frame * probs[0].size() + pattern;}
    Double* getFactor(int index);
    string verbalDescription(int index);
    // data members
    int order;
    vector<Double> probs[3];
};

/*
 * BinnedMMGroup
 * Contains features for bins of a probability (or a score). E.g. when a probability is mapped to a piecewise constant function.
 * Used to make features for signal models where a pattern probability previously existed. Avoids overfitting.
 */
class BinnedMMGroup : public MMGroup{
public:
    BinnedMMGroup(string name, int monotonic = 0) : MMGroup(name){this->monotonic = monotonic;};
    BinnedMMGroup() : MMGroup(){reset();};
    void reset(){nbins=0; origprobs.clear(); bb.clear(); avprobs.clear();}
    void removeOrigs(){origprobs.clear();}
    void addProb(Double p){origprobs.push_back(p);};
    void trainBins(int maxNumBins);
    void monotonify(); // change weights so function is monotonic (again)
    void printBoundaries();
    int getIndex(Double p);
    Double factor(Double p){if (nbins==0) return p; return avprobs[getIndex(p)];}
    void write(ostream &out);
    void read(ifstream &in);
    Double getMinProb(float qthresh);
    // inherited
    Double* getFactor(int index);
    string verbalDescription(int index);
    void smooth(){ monotonify(); }
    int getNumPars() {numPars = avprobs.size(); return numPars;}
    // data members
    int monotonic;            // +1: monotonic increasing, -1: monotonic decreasing, 0 not monotonic
    int nbins;                // number of bins
    vector<Double> origprobs; // original probability set used for "training" of bb and avprobs
    vector<Double> bb;        // bin boundaries,
    vector<Double> avprobs;   // probs for bins [-infty, bb[0]), [bb[0],bb[1]), ... , [bb[nbins-1], infty]
};

/*
 * Parameters
 * holds all parameters that are trained through the CRF
 */
class Parameters {
 public:
    Parameters(){
	weights.reserve(5000);             // adding parameters up to this number won't force a resize of the vector
	counts.reserve(5000);
	alpha.reserve(5000);
	parptrs = NULL;
	mmgroups = NULL;
    }
    ~Parameters(){}
    int addMMGroup(MMGroup* mmgroup, double alphavalue = 1.0);      // returns the offset
    void addMerkmal(Double* parptr, double alphavalue = 1.0);
    void addCount(int index, int count=1) { counts[index] += count; }
    void resetCounts(){ for (int i=0; i<counts.size(); i++) counts[i]=0; };
    void addWeights(vector<double> &h);
    void setWeights(vector<double> &v);
    void smoothFeatures();
    void updatePars();       // set pars to the exp of weights
    void updateWeights();    // set weights to the log of pars
    double getWeight(int i){ return weights[i]; }
    double getAlpha(int i){return alpha[i]; }
    string verbalDescription(int index);
    void print(int numprint, bool countsOnly) { print<int>(counts, numprint, countsOnly); }
    template <class T>
    void print(vector<T> &counts, int numprint, bool countsOnly); // numprint < 0 => print all

//    void print(vector<int> &counts, bool verbose, bool countsOnly);
    vector<int> getCounts(){ return counts; }
    int size() { return weights.size(); }
    vector<double> getWeights(){ return weights; }
 private:
    vector<double> weights;               // the weightsparameters of the CRF (weight vector), logarithms of parameters
    vector<int> counts;                   // holds counts of how often each feature was used in a parse
    vector<double> alpha;                 // regularization weight: alpha[j] large => weight[j] is allowed to vary more
    vector<Double*> *parptrs;             // pointer to vector of pointers, that point to the actual Double variables used in the prediction
    vector<MMGroup*> *mmgroups;           // ordered list of all contained MMGroups
};

template <class T>
void Parameters::print(vector<T> &counts, int numprint, bool countsOnly){
    T countsum = 0;
    double weightsum = 0.0;

    cout << "[CRF parameters]: ";
    if (!mmgroups){
	cout << " empty " << endl;
	return;
    }

    cout << weights.size() << " parameters" << endl;
    int grpidx=0;
    MMGroup* mmgroup = (*mmgroups)[grpidx];
    for (int index=0; index < weights.size(); index++){
	if (numprint > 0 && index == numprint/2 && weights.size() > numprint)
	    cout << "..." << endl;
	if (numprint < 0 || index < numprint/2 || weights.size()-index <= numprint/2) {
	    cout << index << "\t";
	    while(mmgroup->getOffset() + mmgroup->getNumPars() <= index && grpidx < mmgroups->size()-1){
		grpidx++;
		mmgroup = (*mmgroups)[grpidx];
	    }
	    if (grpidx >= mmgroups->size() || mmgroup->getOffset() > index){
		throw ProjectError("Parameters::print: inconsistency");
	    }
	    cout << left << setw(30) << mmgroup->verbalDescription(index - mmgroup->getOffset()) << "\t";
	    cout << "alpha=" << alpha[index];
	    cout << "\tc=" << counts[index];
	    if (!countsOnly)
		cout << "\tw=" << setprecision(6) << weights[index] << "\tp=" << *(*parptrs)[index];
	    cout << endl;
	}
	countsum += counts[index];
	weightsum += weights[index];
    }
    cout << "***\t" << setw(30) << "average over all features" << "\t";
    cout << "c=" << countsum;
    if (!countsOnly)
	cout << "\tw=" << weightsum/weights.size() << "\tp=" << LLDouble::exp(weightsum/weights.size());
    cout << endl;
}


/*
 * CRF
 * implements functions for training Conditional Random Fields
 */
class CRF {
public:
    static void train(Parameters* startPars, vector<AnnoSequence*> trainGenes, FeatureCollection &extrinsicFeatures);
    static void onlineLargeMarginTraining(Parameters* startPars, vector<AnnoSequence*> trainGenes, FeatureCollection &extrinsicFeatures);
    static void improvedIterativeScaling(Parameters* startPars, vector<AnnoSequence*> trainGenes, FeatureCollection &extrinsicFeatures);
    template <class U, class V>
    static void compareWeights(Parameters* parameters, vector<U> &startWeights, vector<V> &endWeights, int numPrint=100);
    static double lossFct(Transcript* correct, Transcript* predicted);
    static void setEvalAnnoSeqs(vector<AnnoSequence*> seqs){evalAnnoSeqs = seqs;}
    // this function is only temporary for saving intermediate parameter files
    static void setPrintPars(int ostatecount, StateModel **ostates){
	statecount = ostatecount;
	states = ostates;
    }
private:
    static vector<double> capOutliers(vector<double> bs);
    static vector<AnnoSequence*> evalAnnoSeqs;
    static FeatureCollection *evalExtrinsicFeatures;
    static int statecount;
    static StateModel **states;
};

struct LessDereference {
    template <class T>
    bool operator()(const T * lhs, const T * rhs) const {
        return *lhs < *rhs;
    }
};

template <class U, class V>
void CRF::compareWeights(Parameters* parameters, vector<U> &startWeights, vector<V> &endWeights, int numPrint){
    int n = numPrint;
    if (n > startWeights.size())
        n = startWeights.size();
    if (n > 0)
	cout << "sorted list of the " << n << " most significant changes." << endl;
    vector<double> absDiffs;
    absDiffs.reserve(startWeights.size());
    for (int i=0; i < startWeights.size(); i++)
        absDiffs.push_back(abs(startWeights[i] - endWeights[i]));

    std::vector<const double *> pointer;
    pointer.reserve(absDiffs.size());
    const double *const start = &absDiffs[0];
    const double *const end = start + absDiffs.size();

    for (const double * iter = start; iter != end; ++iter)
        pointer.push_back(iter);


    std::sort(pointer.begin(), pointer.end(), LessDereference());
    for (int i = pointer.size()-1; i >= pointer.size()-n && i >= 0; i--) {
        const double * p = pointer[i];
        int idx = p-start;
	if (endWeights[idx] != startWeights[idx]){
	    cout << parameters->verbalDescription(idx) << "\t" << startWeights[idx] << " --> " << endWeights[idx]
		 << " diff= " << endWeights[idx] - startWeights[idx] << "\t" << setprecision(4) << LLDouble::exp(endWeights[idx] - startWeights[idx]) << endl;
	}
    }
}

#endif  //  _MERKMAL_HH

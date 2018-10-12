/*
 * igenicmodel.hh
 *
 * License: Artistic License, see file LICENSE.TXT or 
 *          https://opensource.org/licenses/artistic-license-1.0
 */
 
#ifndef __IGENICMODEL_HH
#define __IGENICMODEL_HH

#include "statemodel.hh"

/**
 * @author Mario Stanke
 * 
 */
class IGenicModelError : public ProjectError {
public:
    IGenicModelError(string msg) : ProjectError(msg) {}
};

/**
 * @brief intergenic region
 * 
 * @author Mario Stanke
 * 
 */
class IGenicModel : public StateModel {
public:
    IGenicModel() : gweight(1) {
	init();
    }
    ~IGenicModel() {}
    
    StateType getStateType() const {
	return igenic;
    }

    void buildModel         ( const AnnoSequence* annoseq, int parIndex );
    void registerPars       ( Parameters* parameters);
    void printProbabilities ( int zusNumber, BaseCount *bc, const char* suffix = NULL );
    void initAlgorithms     ( Matrix<Double>&, int);
    void viterbiForwardAndSampling(ViterbiMatrixType&, ViterbiMatrixType&, int, int,
				   AlgorithmVariant, OptionListItem&);
    Double emiProbUnderModel(int begin, int end) const;
    static void init();
    static void setSFC(SequenceFeatureCollection *sfc) {
	seqFeatColl = sfc;
    }
    static void resetPars() {}
    static void updateToLocalGC(int from, int to);
    static void readProbabilities(int zusNumber);
    static void readAllParameters();
    static void storeGCPars(int idx);
    static double getGeoProb(){return geoProb;}
    static vector<double> getNucleotideProbs(){return nucProbs[0];}

private:
    void processSequence( const char* start, const char* end );
public:    
  static Integer         k;           // UTR intron and intron
  static PatMMGroup      emiprobs;    // can use this data.
  static PatMMGroup      *GCemiprobs; // array for each GC content class
private:
  Integer                gweight;
  static Double          patpseudocount;
  static Integer         gesbasen;
  static vector<Integer> emicount;
  static vector<vector<Double> > Pls;
  static vector<vector<Double> >* GCPls;
  static vector<vector<double> > nucProbs;
  static int             lastParIndex; // GC-index of current parameter set
  static int             verbosity;
  static double          geoProb;
};

#endif    //  __IGENICMODEL_HH

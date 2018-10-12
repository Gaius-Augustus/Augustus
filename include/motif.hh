/*
 * motif.hh
 *
 * License: Artistic License, see file LICENSE.TXT or 
 *          https://opensource.org/licenses/artistic-license-1.0
 */

#ifndef _MOTIF_HH
#define _MOTIF_HH

// project includes
#include "matrix.hh"
#include "geneticcode.hh"  // for Seq2Int

/**
 * @details The sequences are weighed according to the relative frequencies of the 4 nucleotides in them 
 * and the relative frequencies of the nucleotides in the input sequence <br>
 * equalWeights     : all are equally important <br>
 * gcContentClass   : consider only sequences in the same class of gc-Content (isochore) as the
 *                    input sequence (like Genscan does) <br>
 * multiNormalKernel: use an arbitrary weighing matrix as inverse of the covariance matrix
 *                    of a multvariate normal distribution <br>
 * 
 * @see BaseCount
 * 
 * @author Mario Stanke
 */
enum WeighingType { equalWeights, gcContentClasses, multiNormalKernel};

/**
 * @author Mario Stanke
 */
class BaseCount{
public:
    int a;
    int c;
    int g;
    int t;
    static WeighingType weithType;
    static Matrix<double> weighingMatrix;

    double ra, rc, rg, rt;
    BaseCount();
    BaseCount(int a, int c, int g, int t);
    BaseCount(const char *sequence, int len=-1);
    ~BaseCount() {};
    static void init();
    void normalize();
    void addSequence(const char *sequence, int len);
    void addCharacter(char nucleotide, bool subtract=false);
    void reverse();
    static int weight( BaseCount bc1, BaseCount bc2);
    static double doubleWeight( BaseCount bc1, BaseCount bc2);
    static int gcContentWeight(BaseCount bc1, BaseCount bc2);
    static int gcContentClass (double gcContent);
    static int gcContentClassWeight(BaseCount bc1, BaseCount bc2);
    static double multiNormalKernelWeight(BaseCount bc1, BaseCount bc2);
    static void setWeightMatrix(string matrixFileName);
    static double phi(double x, double sigma);
  
};

ostream& operator<<( ostream& out, const BaseCount& bc );

/**
 * @author Mario Stanke
 */
struct Composition{
    double a;
    double c;
    double g;
    double t;
    
    Composition() {
	a = 0.0;
	c = 0.0;
	g = 0.0;
	t = 0.0;
    }
    Composition(BaseCount bc) {
	double sum = bc.a + bc.c + bc.g + bc.t;
	a = bc.a / sum;
	c = bc.c / sum;
	g = bc.g / sum;
	t = bc.t / sum;
    }
};

/**
 * @author Mario Stanke
 */
class Motif {
public:
    int n;
    int k;
    int numSeqs;    // number of unweighed sequences 
    int neighbors;
    int pseudocount;
    
    Motif() :
	n(0), k(0), neighbors(0),
	pseudocount(1),
	windowProbs(NULL),
	windowCounts(NULL),
	s2i(0) {}
    Motif & operator = (const Motif & other);

    Motif(int length, int memory=0, int pseudocount = 1, int neighbors = 0);
    ~Motif();
    /*
     * add one sequence to the training set of the motif
     * seq is the beginning of the motiv, but
     * seq[-k] ... seq[n-1] or seq[0] ... seq[n+k-1] (reverse case) must be accessible!
     */
    void addSequence(const char* seq, int weight = 1, bool reverse=false);
    void makeProbs();
    void printProbs();
    Double seqProb(const char* seq, bool reverse=false, bool complement=false);
    void write(ofstream &out);
    void read(ifstream &in);
    void clearCounts();
    char* getSampleDNA();
    Double getProbThreshold (double q, int numSamples = 10000);
private:
    vector<Double> *windowProbs;
    vector<int>    *windowCounts;
    Seq2Int s2i;
};

/**
 * @author Mario Stanke
 */
class ContentDecomposition {
public:
    int n;
    BaseCount *zus;
    ContentDecomposition() : n(0), zus(0) { 
	setProperties();
    }
    BaseCount getBaseCount(int i);
    int getNearestBaseCountIndex(BaseCount bc);
private:
    void setProperties();
    void makeDecomposition();
};

/**
 * @brief holds the stepwise constant function of GC content class indices
 * 
 * @author Mario Stanke
 */
class ContentStairs {
public:
  ContentStairs();
  ~ContentStairs();
  void computeStairs(const char* dna);
  int getNextStep(int from);
  int *idx; // GC content class index for each position of dna
  const char *dna; // just a pointer, to check whether update is necessary
  int n; // dna length
  map<int,int> nextStep;
private:
  int GCwinsize;
};

#endif

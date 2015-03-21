/*****************************************************************************\
 * Filename : motif.hh
 * Author   : Mario Stanke
 * Project  : HMM
 *
 * Copyright: ©Stanke
 *
 * Description: 
 *
 *
 * Date       |   Author              |  Changes
 *------------|-----------------------|----------------------------------------
 * 08.08.2002 | Mario Stanke          | Creation of the file.
\*****************************************************************************/

#ifndef _MOTIF_HH
#define _MOTIF_HH

// project includes
#include "matrix.hh"
#include "geneticcode.hh"  // for Seq2Int


/*
 * The sequences are weighed according to the relative frequencies of the 4 nucleotides in them 
 * and the relative frequencies of the nukleotides in the input sequence
 * equalWeights     : all are equally important
 * gcContentClass   : consider only sequences in the same class of gc-Content (isochore) as the
 *                    input sequence (like Genscan does)
 * multiNormalKernel: use an arbitrary weighing matrix as inverse of the covariance matrix
 *                    of a multvariate normal distribution
 * see class BaseCount
 */
enum WeighingType { equalWeights, gcContentClasses, multiNormalKernel};

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

/*
 * class Motif
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

/*
 * ContentStairs
 * holds the stepwize constant function of GC content class indices
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

/*****************************************************************************\
 * Filename : contTimeMC.hh
 * Authors  : Mario Stanke
 *
 * Description: continuous-time Markov chains for sequence evolution
 *              codon models
 *
 * Date       |   Author              |  Changes
 *------------|-----------------------|------------------------------------------
 * 17.02.2013 | Mario Stanke          | creation of the class
 * 27.05.2013 | Stefanie König        | creation of the class ExonEvo
 *            |                       | and a base class Evo from which ExonEvo
 *            |                       | and CodonEvo inherit
\******************************************************************************/
#ifndef _CONTTIMEMC_HH
#define _CONTTIMEMC_HH

// project includes
#include "geneticcode.hh"
#include "matrix.hh"

// standard C/C++ includes
#include <iostream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>

// forward declarations
class PhyloTree;

using namespace std;

/*
 * abstract base class for codon and exon evolution
 */
class Evo {
public:
    Evo(int s) : states(s), m(0) {};
    virtual ~Evo();
    int getNumStates(){return states;}

    /* Determine the branch lengths for which matrices P should be stored, at most m.
     * For trees with few species, say <10, this could be all branch lengths b occuring
     * in the tree. A value of m=-1 means to store all lengths in b. Memory may be saved
     * by building m clusters roughly representing all lengths in b.
     */
    void setBranchLengths(vector<double> b, int m=-1);
    void printBranchLengths();

    virtual void computeLogPmatrices()=0; // precomputes and stores the array of matrices
    virtual void addBranchLength(double b)=0; //add new branch length b and matrix P(b)

    double getPi(int i) const {return pi[i];}  //equilibrium frequency of state i
    double getLogPi(int i) const {return log(getPi(i));}

    /*
     * Returns a pointer to the probability matrix for which t comes closest.
     * u is the index to a second parameter for which the probability matrices
     * are stored for different values.
     * in the case of codon evolution u is the index to omega (dN/ds, selection)
     * in the case of exon evolution probability matrices currently only
     * depend on t and thus u is set to 0. It is however conceivable
     * to allow different rates of exon loss in future.
     */
    gsl_matrix *getSubMatrixP(int u, double t); 
    gsl_matrix *getSubMatrixLogP(int u, double t); 

protected:
    int findClosestIndex(vector<double> &v, double val);

protected:
    int states; //number of states (64 in codon model and (currently) 2 in exon model)
    int m; // number of branch lengths (times) for which P's are stored
    double *pi;
    vector<double> times; // sorted vector of branch lengths
    Matrix<gsl_matrix *> allPs; // parametrized probability matrices
    Matrix<gsl_matrix *> allLogPs; // parametrized log probability matrices
};

/*
 * class CodonEvo
 * Computes and stores 64x64 codon substitution probability matrices
 * P(t,omega,kappa,pi) = exp(Q(omega,kappa,pi)*t)
 * for a sample of values for t and omega (dN/dS, selection)
 * and a fixed given value of kappa (transition/transversion ratio)
 * and vector pi (codon usage).
 * Matrices are precomputed and retrieval of P(t,omega) is then a constant time lookup
 * as this is needed very often during estimation of omega for a given tree.
 */
class CodonEvo : public Evo {
public:
    CodonEvo() : Evo(64), k(0), kappa(4.0) {};
    ~CodonEvo();
    void setKappa(double kappa){ this->kappa = kappa;} // transition/transversion ratio
    void setPi(double *pi); // codon usage
   
    /*
     * Chooses k values for omega around 1, e.g. 0.8, 1, 1.2 for k=3.
     * Later, one can see which of these values of omega gives the maximum likelihood.
     */
    void setOmegas(int k);
    void setPrior(double sigma = 0.2);
    double getPrior(int u){return omegaPrior[u];}
    vector<double> *getOmegas(){return &omegas;}
    double getOmega(int u){return omegas[u];}
    int getK(){ return k;}
    void printOmegas();

    void computeLogPmatrices(); // precomputes and stores the array of matrices

    /*
     * Returns a pointer to the 64x64 substitution probability matrix
     * for which t and omega come closest to scored values.
     */
    gsl_matrix *getSubMatrixLogP(double omega, double t);
    
    /*
     * Computes the substitution probability 
     * P( X(t)=to | X(0)=from, omega), where 1 <= from,to <= 64 are codons
     * If you want to call this for many values of 'from' and 'to', use rather getSubMatrixP above.
     */
    double subLogProb(int from, int to, double omega, double t){
	gsl_matrix *P = getSubMatrixLogP(omega, t);
	return gsl_matrix_get(P, from, to);
    }
    /*
     * log likelihood of exon candidate pair e1, e1 (required: equal number of codons, gapless)
     * Used for testing. The pruning algorithm requires tuples of codons instead.
     */
    double logLik(const char *e1, const char *e2, int u, double t); // output variables
    /*
     * Estimate selection omega on pair of codon sequences.
     */
    double estOmegaOnSeqPair(const char *e1, const char *e2, double t, // input variables
			     int& numAliCodons, int &numSynSubst, int &numNonSynSubst); // output variables
    /*
     * Estimate omega on a sequence of codon tuples.
     * only for testing, may need adjustment
     */
    double estOmegaOnSeqTuple(vector<string> &seqtuple, PhyloTree *tree,
			      int &subst); // output variables
    // store log likeliehoods for all omegas for one codon tuple
    vector<double> loglikForCodonTuple(vector<string> &seqtuple, PhyloTree *tree);
    // for GRK proposal
    double graphOmegaOnCodonAli(vector<string> &seqtuple, PhyloTree *tree);
    /* 
     * add new branch length b
     * this function is currently not needed, since the pruning algorithm does change the phylogenetic tree
     */
    void addBranchLength(double b){} 

private:
    int k; // number of different omega values for which P's are stored
    double kappa;
    vector<double> omegas; // sorted vector of omegas (contains values below, around and above 1)
    vector<double> omegaPrior; // prior distribution on omega, centered at 1
};

/*
 * 64x64 rate matrix Q for codon substitutions as in Yang, "Computational Molecular Evolution", (2.7)
 */

gsl_matrix *getCodonRateMatrix(double *pi,    // codon usage, normalized vector with 64 elements
			       double omega,  // dN/dS, nonsynonymous/synonymous ratio
			       double kappa); // transition/transversion ratio, usually >1

/*
 * perform a decompososition of the rate matrix as Q = U * diag(lambda) * U^{-1}
 * returns 0 if sucessful
 */
int eigendecompose(gsl_matrix *Q,       // input
		   double *pi,          // input, must be same as in construction of Q
		   gsl_vector *&lambda, // output, memory for lambda, U and Uinv is allocated within the function
		   gsl_matrix *&U,      // output
		   gsl_matrix *&Uinv);  // output

/*
 * compute the matrix exponential P(t) = exp(Q*t),
 * where Q is given by U, diag(lambda) and U^{-1} as computed by eigendecompose
 */

gsl_matrix *expQt(double t, gsl_vector *lambda, gsl_matrix *U, gsl_matrix *Uinv);

void printCodonMatrix(gsl_matrix *M);

/*
 * computes the element-wise natrual logarithm of a matrix P
 */
gsl_matrix *log(gsl_matrix *P, int states);

/*
 * class ExonEvo
 * Computes and stores 2x2 probability matrices for exon gain/loss events
 * P(t) = exp(Q(lambda, mu)*t) for a sample of values for t
 * and a fixed given value for lambda (rate of exon gain) and mu (rate of exon loss)
 * Matrices are precomputed and retrieval of P(t) is then a constant time lookup
 */
class ExonEvo : public Evo{

public:
    ExonEvo() : Evo(2)  {
	setLambda();
	setMu();
	setPhyloFactor();
	setPi();
    };
    void setPi();
    void setPhyloFactor();
    void setPhyloFactor(double k){phylo_factor=k;}
    double getPhyloFactor() const{return phylo_factor;}
    void setMu();
    double getMu() const{return mu;}
    void setLambda();
    double getLambda() const{return lambda;}
    gsl_matrix *computeP(double t); //computes transition matrix P for time t
    void computeLogPmatrices();
    void addBranchLength(double b); //add new branch length b and matrix P(b)
    double minBranchLength() const {return times.front();}
   
private:
    double mu;            // rate for exon loss
    double lambda;        // rate for exon gain
    double phylo_factor;  // tuning parameter to fit phylogenetic score to scores in the graph

};

/*
 * class Parsimony
 * transitions from state i to state j are independ of the branch length
 * and cost -1 if i!=j and 0 if i==j.
 * can be used to calculate the minimum nunber of codon substitutions
 * for a given tree topology
 */
class Parsimony : public Evo{

public:
    Parsimony() : Evo(64) {
	this->m = 1;
	this->times.push_back(1);
	this->pi = new double[64];
	for(int i=0; i<64; i++){
	    pi[i]=1;
	}
    };
    void computeLogPmatrices();
    void addBranchLength(double b){} 
};
#endif    // _CONTTIMEMC_HH

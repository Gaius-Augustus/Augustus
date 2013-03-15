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

using namespace std;

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
class CodonEvo {
public:
    CodonEvo() : k(0), m(0), kappa(2) {};
    ~CodonEvo();
    void setKappa(double kappa){ this->kappa = kappa;} // transition/transversion ratio
    void setPi(double *pi); // codon usage

    /* Determine the branch lengths for which matrices P should be stored, at most m.
     * For trees with few species, say <10, this could be all branch lengths b occuring
     * in the tree. A value of m=-1 means to store all lengths in b. Memory may be saved
     * by building m clusters roughly representing all lengths in b.
     */
    void setBranchLengths(vector<double> &b, int m=-1);
    void printBranchLengths();

    /*
     * Chooses k values for omega around 1, e.g. 0.8, 1, 1.2 for k=3.
     * Later, one can see which of these values of omega gives the maximum likelihood.
     */
    void setOmegas(int k);
    vector<double> *getOmegas(){return &omegas;}
    double getOmega(int u){return omegas[u];}
    int getK(){ return k;}
    void printOmegas();

    void computeLogPmatrices(); // precomputes and stores the array of matrices

    /*
     * Returns a pointer to the 64x64 substitution probability matrix
     * for which t and omega come closest to scored values.
     * The second version is for the (typical) case when omega is already given by an index u to omegas;
     */
    gsl_matrix *getSubMatrixLogP(double omega, double t);
    gsl_matrix *getSubMatrixLogP(int u, double t);
    
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
private:
    int findClosestIndex(vector<double> &v, double val);

private:
    int k; // number of different omega values for which P's are stored
    int m; // number of branch lengths (times) for which P's are stored
    double kappa;
    double *pi;
    vector<double> times; // sorted vector of branch lengths
    vector<double> omegas; // sorted vector of omegas (contains values below, around and above 1)
    Matrix<gsl_matrix *> allPs; // parametrized log probability matrices
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
#endif    // _CONTTIMEMC_HH

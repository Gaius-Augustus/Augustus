/*
 * contTimeMC.hh
 * @author Mario Stanke
 * @author Stefanie König
 *
 * License: Artistic License, see file LICENSE.TXT or 
 *          https://opensource.org/licenses/artistic-license-1.0
 */

#ifndef _CONTTIMEMC_HH
#define _CONTTIMEMC_HH

#include "matrix.hh"

// standard C/C++ includes
#include <cmath>
#include <iostream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>


using namespace std;

/**
 * @brief abstract base class for codon and exon evolution
 * 
 * @author Mario Stanke
 * @author Stefanie König
 */
class Evo {
public:
  Evo(int s) : states(s), m(0), pi(NULL) {};
    virtual ~Evo();
    int getNumStates(){return states;}

    /* Determine the branch lengths for which matrices P should be stored, at most m.
     * For trees with few species, say <10, this could be all branch lengths b occurring
     * in the tree. A value of m=-1 means to store all lengths in b. Memory may be saved
     * by building m clusters roughly representing all lengths in b.
     */
    void setBranchLengths(vector<double> b, int m=-1);
    void printBranchLengths();

    virtual void getRateMatrices() {}; // precomputes or reads the rate matrices
    virtual void computeLogPmatrices()=0; // precomputes and stores the array of matrices
    virtual void addBranchLength(double b)=0; //add new branch length b and matrix P(b)

    double getPi(int i) const {return pi[i];}  //equilibrium frequency of state i
    double getLogPi(int i) const {return log(getPi(i));}

    gsl_matrix *getSubMatrixQ(int u){ return allQs[u]; }
    
    /*
     * Returns a pointer to the probability matrix for which t comes closest.
     * u is the index to a second parameter for which the probability matrices
     * are stored for different values.
     * in the case of codon evolution u is the index to omega (dN/dS, selection)
     * in the case of exon evolution probability matrices currently only
     * depend on t and thus u is set to 0. It is however conceivable
     * to allow different rates of exon loss in future.
     */
    gsl_matrix *getSubMatrixP(int u, double t);
    gsl_matrix *getSubMatrixLogP(int u, double t); 

    /*
     * compute the matrix exponential P(t) = exp(Q*t),
     * where Q is given by U, diag(lambda) and U^{-1} as computed by eigendecompose
     */
    gsl_matrix *expQt(double t, gsl_vector *lambda, gsl_matrix *U, gsl_matrix *Uinv);

protected:
    int findClosestIndex(vector<double> &v, double val);

protected:
    const int states; //number of states (64 in codon model and (currently) 2 in exon model)
    int m; // number of branch lengths (times) for which P's are stored
    double *pi;
    vector<double> times; // sorted vector of branch lengths
    vector<gsl_matrix *> allQs; // rate matrices
    Matrix<gsl_matrix *> allPs; // Parameterized probability matrices
    Matrix<gsl_matrix *> allLogPs; // Parameterized log probability matrices

};

/*
 * perform a decompososition of the rate matrix as Q = U * diag(lambda) * U^{-1}
 * returns 0 if sucessful
 */
int eigendecompose(gsl_matrix *Q,       // input
		   double *pi,          // input, must be same as in construction of Q
		   gsl_vector *&lambda, // output, memory for lambda, U and Uinv is allocated within the function
		   gsl_matrix *&U,      // output
		   gsl_matrix *&Uinv);  // output

// eigen decomposition was successful, if Q - U * diag(lambda) * Uinv is close to the 0-matrix
void validateEigenDecomp(gsl_matrix *Q, gsl_vector *lambda, gsl_matrix *U, gsl_matrix *Uinv, int states);

void printCodonMatrix(gsl_matrix *M);

/*
 * computes the element-wise natrual logarithm of a matrix P
 */
gsl_matrix *log(gsl_matrix *P, int states);

/*
 * @brief class ExonEvo
 * @details Computes and stores 2x2 probability matrices for exon gain/loss events
 * P(t) = exp(Q(lambda, mu)*t) for a sample of values for t
 * and a fixed given value for lambda (rate of exon gain) and mu (rate of exon loss)
 * Matrices are precomputed and retrieval of P(t) is then a constant time lookup
 * 
 * @author Mario Stanke
 * @author Stefanie König
 */
class ExonEvo : public Evo{

public:
    ExonEvo(int num_states=2) : Evo(num_states), U(NULL), Uinv(NULL), l(NULL)  {
	setLambda();
	setMu();
	setAliErr();
	// mu, lambda and ali_error have to be set before calculating equilibrium frequencies Pi
	setPi();
    };
    ~ExonEvo(){
        if (U)
            gsl_matrix_free(U);
        if (Uinv)
            gsl_matrix_free(Uinv);
        if (l)
            gsl_vector_free(l);
    }
    void setPi();
    void setLambda(double lambda = 0.0001);
    void setMu(double mu = 0.0001);
    void setAliErr(double ali_error = 0.1);
    double getLambda() const{return lambda;}
    double getMu() const{return mu;}
    double getAliErr() const{return ali_error;}

    void getRateMatrices(); // precomputes or reads the rate matrices
    void computeLogPmatrices();
    void addBranchLength(double b); //add new branch length b and matrix P(b)
    double minBranchLength() const {return times.front();}
    int eigendecompose(gsl_matrix *Q);
    gsl_matrix *expQt(double t) {return Evo::expQt(t,l,U,Uinv);}
    gsl_matrix *getExonRateMatrix();
    void validateEigenDecomp(gsl_matrix *Q){::validateEigenDecomp(Q,l,U,Uinv,states);}
   
private:
    double mu;            // rate for exon loss
    double lambda;        // rate for exon gain
    double ali_error;     // rate for alignment error

    // eigen decomposition of exon rate matrix Q = U * diag(l_1,...,l_n) * Uinv
    gsl_matrix *U;
    gsl_matrix *Uinv;
    gsl_vector *l;
};

/*
 * @brief class Parsimony
 * @details transitions from state i to state j are independent of the branch length
 * and cost -1 if i!=j and 0 if i==j.
 * can be used to calculate the minimum number of codon substitutions
 * for a given tree topology
 * 
 * @author Mario Stanke
 * @author Stefanie König
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

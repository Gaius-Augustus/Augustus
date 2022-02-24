/*
 * codonevo.hh
 *
 * License: Artistic License, see file LICENSE.TXT or 
 *          https://opensource.org/licenses/artistic-license-1.0
 */

#ifndef _CODONEVODISCR_HH
#define _CODONEVODISCR_HH

#include "contTimeMC.hh"
#include "geneticcode.hh"

enum MatrixFormat {fmt_GSL, fmt_TXT};

//forward declarations
class PhyloTree;

/**
 * @brief class CodonEvoDiscr
 * @details Computes and stores 64x64 codon substitution probability matrices
 * P(t,omega,kappa,pi) = exp(Q(omega,kappa,pi)*t)
 * for a sample of values for t and omega (dN/dS, selection)
 * and a fixed given value of kappa (transition/transversion ratio)
 * and vector pi (codon usage).
 * Matrices are precomputed and retrieval of P(t,omega) is then a constant time lookup
 * as this is needed very often during estimation of omega for a given tree.
 * 
 * @author Mario Stanke
 * @author Stefanie König
 * adapted by Giovanna Migliorelli
 * */

// clamsa related code
class CodonEvoDiscr : public Evo {
public:
    CodonEvoDiscr() : Evo(64), M(0), Theta(NULL), bias(NULL), w(NULL) {};
    ~CodonEvoDiscr();

    /*
     * in the discriminative settings, k stands for the number of models matrices were computed on the base of
     */
    void setM();
    int getM(){ return M;}
    void writeRateMatrices(string filename); // serialize to file
    void readMatrices(MatrixFormat fmt = fmt_TXT);

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
  
    // store log likeliehoods for all omegas for one codon tuple
    vector<double> loglikForCodonTuple(vector<string> &seqtuple, PhyloTree *ctree){
        int dummysubs = 0;
        return loglikForCodonTuple(seqtuple, ctree, NULL, dummysubs);
    }
    
    /* loglikForCodonTuple
    * @param[in] seqtuple input codon alignment with as many rows as ctree has leaves
    * @param[in] ctree codon scaled tree
    * @param[in] tree optionally a second tree for running the Fitch alg
    * @return the a vector of k log likelihoods 
    */
    vector<double> loglikForCodonTuple(vector<string> &seqtuple, PhyloTree *ctree, PhyloTree *tree, int &subs);
    /* 
     * add new branch length b
     * this function is currently not needed, since the pruning algorithm does change the phylogenetic tree
     */
    void addBranchLength(double b){} 

    void produceComb(vector<vector<int> >& config){}
    double getProb(const vector<double> &lls);
    
private:
    int M; // number of rate matrices for which P's are stored
    vector<double*> piArr;  // differently from the case of omega here each model M has its own pi
    gsl_matrix *Theta; // M x 2 matrix from ClaMSA logreg layer
    double *bias;  // vector of 2 biases from ClaMSA logreg layer
    double *w; // weight vector summarizing equivalently Theta and bias
};

#endif    // _CODONEVODISCR_HH

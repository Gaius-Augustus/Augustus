/*
 * codonevo.hh
 *
 * License: Artistic License, see file LICENSE.TXT or 
 *          https://opensource.org/licenses/artistic-license-1.0
 */

#ifndef _CODONEVO_HH
#define _CODONEVO_HH

#include "contTimeMC.hh"
#include "geneticcode.hh"

//forward declarations
class PhyloTree;

/**
 * @brief class CodonEvo
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
 */
class CodonEvo : public Evo {
public:
    CodonEvo() : Evo(64), k(0), kappa(4.0) {
        // initialize pi with the uniform distribution on all sense codons
        pi = new double[64];
        double normsum = 0.0;
        for (unsigned c=0; c<64; c++){
            if (GeneticCode::translate(c) == '*')
                pi[c] = 0.0;
            else
                pi[c] = 1.0;
            normsum += pi[c];
        }
        
        for (int c=0; c<64; c++)
            pi[c] /= normsum; 
    };
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
    void setAAPostProbs();
    void getRateMatrices(); // precomputes or reads the rate matrices
    void writeRateMatrices(string filename); // serialize to file
    void readRateMatrices(string filename);  // deserialize from file
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

  vector<double> loglikForCodonTuple(vector<string> &seqtuple, PhyloTree *ctree){
    int i = 0;
    return loglikForCodonTuple(seqtuple, ctree, NULL, i);
  }
  vector<double> loglikForCodonTuple(vector<string> &seqtuple, PhyloTree *ctree, PhyloTree *tree, int &subs);
    // for GRK proposal
  double graphOmegaOnCodonAli(vector<string> &seqtuple, PhyloTree *tree, int refSpeciesIdx = 0);
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
    vector<double> aaUsage; // amino acid usage for incorporation into rate matrix
    vector<vector<double> > aaPostProb; // retreived from BLOSUM (amino acid substitution rate matrix)
};

/*
 * 64x64 rate matrix Q for codon substitutions as in Yang, "Computational Molecular Evolution", (2.7)
 */

gsl_matrix *getCodonRateMatrix(double *pi,    // codon usage, normalized vector with 64 elements
			       double omega,  // dN/dS, nonsynonymous/synonymous ratio
			       double kappa, // transition/transversion ratio, usually >1
			       vector<vector<double> > *aaPostProb = NULL); // posterior probs for AA substitutions

gsl_matrix *getNonCodingRateMatrix(vector<double> *pi_nuc, double kappa); // rate matrix for non-coding model

#endif    // _CODONEVO_HH

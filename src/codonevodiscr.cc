/*
 * codonevo.cc
 *
 * License: Artistic License, see file LICENSE.TXT or 
 *          https://opensource.org/licenses/artistic-license-1.0
 */

// project includes
#include "codonevodiscr.hh"
#include "geneMSA.hh"
#include "phylotree.hh"
#include "igenicmodel.hh"
#include "properties.hh"

#include <stdio.h>


void CodonEvoDiscr::setM(){
    int M = 3; // number of clamsa models
    try {
        M = Properties::getIntProperty("/CompPred/clamsa_numModels");
    } catch(...){}
    
    cout << "CodonEvoDiscr::setM = " << M << endl;
    for (int u=0; u<piArr.size(); ++u)
        delete[] piArr[u];
    piArr.clear();

    if (M>0){
        piArr.resize(M, NULL);
        for (int u=0; u<piArr.size(); ++u)
            piArr[u] = new double[64];
    }
    this->M = M;
}

void CodonEvoDiscr::writeRateMatrices(string filename){
    gsl_matrix *Q;
    int status;
    FILE *Qfile;
    Qfile = fopen(filename.c_str(), "w");
    if (!Qfile)
        throw ProjectError("Could not write to rate matrix file " + filename);
    
    for (int u=0; u<M; u++){
        Q = allQs[u];
        status = gsl_matrix_fprintf(Qfile, Q, "%g");
        if (status)
            throw ProjectError("Could not write " + itoa(u+1) + "-th rate matrix to " + filename);
    }
    fclose(Qfile);
}

void CodonEvoDiscr::readMatrices(MatrixFormat fmt){ 
    string clamsa_Q, clamsa_pi, clamsa_Theta, clamsa_biases;
    try {
        clamsa_Q = Properties::getProperty("/CompPred/clamsa_Q");
    } catch (...) {
        clamsa_Q = "rates-Q.txt";
    }
    try {
        clamsa_pi = Properties::getProperty("/CompPred/clamsa_pi");
    } catch (...) {
        clamsa_pi = "rates-pi.txt";
    }
    try {
        clamsa_Theta = Properties::getProperty("/CompPred/clamsa_Theta");
    } catch (...) {
        clamsa_Theta = "Theta.txt";
    }
    try {
        clamsa_biases = Properties::getProperty("/CompPred/clamsa_biases");
    } catch (...) {
        clamsa_biases = "biases.txt";
    }

    
    if (allQs.size() != M) // empty upon initialization
        allQs.assign(M, NULL);
    
    if (fmt == fmt_GSL){
        int status;
        FILE *Qfile;
        Qfile = fopen(clamsa_Q.c_str(), "r");
        if (!Qfile)
            throw ProjectError("Could not read rate matrix file " + clamsa_Q);
        
        for (int u=0; u<M; u++){
            if (allQs[u] == NULL){
                allQs[u] = gsl_matrix_calloc (64, 64); // initialize Q as
                                                    // the zero-matrix
            }
            status = gsl_matrix_fscanf(Qfile, allQs[u]);
            if (status)
                throw ProjectError("Could not read " + itoa(u+1) +
                                "-th rate matrix from " + clamsa_Q);
        }
        fclose(Qfile);
    } else { // fmt == fmt_TXT
        ifstream Qfile(clamsa_Q);
        
        if (!Qfile.is_open())
            throw ProjectError("Could not read rate matrix file " + clamsa_Q);
        
        double Qij;
        for (int i, j, u=0; u<M; u++){
            if (allQs[u] == NULL)
                allQs[u] = gsl_matrix_calloc (64, 64); // initialize Q as
                                                    // the zero-matrix
            for (i=0; i<64;i++){
                for (j=0; j<64;j++){
                    Qfile >> Qij;
                    gsl_matrix_set(allQs[u], i, j, Qij);
                }
            }
        }
        Qfile.close();
    }

    // read stationary distributions
    ifstream Pfile(clamsa_pi);
    if (!Pfile.is_open())
        throw ProjectError("Could not read initial distribution file " + clamsa_pi);
    
    for (int i,u=0; u<M; u++)
        for (i=0;i<64;i++)
            Pfile >> piArr[u][i];
    Pfile.close();


    // read Theta
    ifstream thetaFile(clamsa_Theta);
        
    if (!thetaFile.is_open())
        throw ProjectError("Could not read logreg matrix file " + clamsa_Theta);
        
    double d;
    if (Theta == NULL)
        Theta = gsl_matrix_calloc (M, 2); // initialize Theta
    for (int u=0; u<M; u++) {
        for (int j=0; j<2; j++){
            thetaFile >> d;
            gsl_matrix_set(Theta, u, j, d);
        }
    }
    thetaFile.close();
    
    // read biases
    ifstream biasesFile(clamsa_biases);
        
    if (!biasesFile.is_open())
        throw ProjectError("Could not read logreg biases file " + clamsa_biases);
        
    if (bias == NULL)
        bias = new double[2]; // 0: negative, 1: positive class
    for (int j=0; j<2; j++)
        biasesFile >> bias[j];
    biasesFile.close();

    // compute weights for logistic regression that output a single probability for being positive (class 1)
    // This is converting from multiclass logistic regression with 2 classes to scalar predictions as described
    // in external-use-of-rate-matrices.pdf of the ClaMSA package.
    // Note that ClaMSA uses the _negative_ mean log likelihood as input to the last Dense layer.
    if (w == NULL)
        w = new double[M+1]; // index 0 is bias, 1 is class 0, ...
    cout << "w: ";
    w[0] = bias[1] - bias[0];
    cout << w[0];
    for (int u=0; u<M; u++){
        // weights of log-likelihoods are multiplied by -1 as -LogP is ised in ClaMSA
        w[1+u] = gsl_matrix_get(Theta, u, 0) - gsl_matrix_get(Theta, u, 1); // intentionally the other way around as for
                                                                            // the bias
        cout << "\t" << w[1+u];
    }
    cout << endl;

    // print a 2x3 sample so it can be checked manually
    cout << "Q[0, 0:1, 0:2] = " << endl << setw(5)
         << gsl_matrix_get(allQs[0], 0, 0) << "\t" << gsl_matrix_get(allQs[0], 0, 1) << "\t" << gsl_matrix_get(allQs[0], 0, 2) << endl
         << gsl_matrix_get(allQs[0], 1, 0) << "\t" << gsl_matrix_get(allQs[0], 1, 1) << "\t" << gsl_matrix_get(allQs[0], 1, 2) << endl;

    /*
     * plausibility tests of parameters
     */
    double tolerance = 1e-4;
    cout << "Running plausibility test for discriminative rate matrix parameters ..." << endl;
    // test row(Q) = 0
    double sum;
    bool foundError = false;
    
    for (int u=0; u<M; u++){
        for (int i=0; i<64; i++){
            sum = 0;
            for (int j=0; j<64; j++)
                sum += gsl_matrix_get(allQs[u], i, j);
            if (std::abs(sum) > tolerance){
                cerr << "Test: The row Q[" << u << "," << i << ":] sums up to " << sum << endl;
                foundError = true;
            }
        }
    }
    if (foundError){
        cerr << "CodonEvoDiscr::readMatrices At least one row of a discriminatively trained codon rate matrix does not sum up to 1.0"
             << endl << "A possible reason is a wrong matrix format." << endl;
    } else {
        cout << "Rate matrix rows sum up to 0. OK" << endl;
    }
    
    
    // test sum(piArr) = 1
    foundError = false;
    for (int u=0; u<M; u++){
        sum = 0;
        for (int i=0; i<64; i++)
            sum += piArr[u][i];
        if (std::abs(sum-1.0) > tolerance){
            cerr << "Test pi[" << u << ",:] sums up to " << sum << endl;
            foundError = true;
        }
    }
    if (foundError){
        cerr << "CodonEvoDiscr::readMatrices At least one stationary distribution (pi) does not sum up to 1.0" << endl;
    } else {
        cout << "Stationary distributions sum up to 1. OK" << endl;
    }
    
    // test pi Q = 0
    foundError = false;
    for (int u=0; u<M; u++){
        for (int j=0; j<64; j++){
            sum = 0.0;
            for (int i=0; i<64; i++){
                sum += piArr[u][i] * gsl_matrix_get(allQs[u], i, j);
            }
            if (std::abs(sum) > tolerance){
                cerr << "pi*Q[" << u << "," << j << "] sums up to " << sum << endl;
                foundError = true;
            }
        }
    }
    if (foundError){
        cerr << "CodonEvoDiscr::readMatrices pi*Q is not the zero vector."
             << "Therefore pi is not the stationary distribution implicitly given by the rate matrix Q." << endl;
    } else {
        cout << "pi is the stationary distribution for rate matrix Q. OK" << endl;
    }
    
    cout << "expected number of mutations per time unit:" << endl;
    for (int u=0; u<M; u++){
        double scale = 0.0; // factor by which matrix must be scaled to achieve one expected mutation per time unit
        for (int i=0; i<64; i++)
            scale += - gsl_matrix_get(allQs[u], i, i) * piArr[u][i];
        cout << u << "\t" << scale << endl;
    }
}

/* 
 * Precompute and store the array of matrices.
 * Takes approximate k * m * 0.015 seconds on greif1: Time for eigendecompose is small compared to time for expQt.
 */
/* 
 * Precompute and store the array of matrices.
 * Takes approximate k * m * 0.015 seconds on greif1: Time for eigendecompose is small compared to time for expQt.
 */
void CodonEvoDiscr::computeLogPmatrices(){
    double t;
    gsl_matrix *Q, *U, *Uinv, *P;
    gsl_vector *lambda;
    int status;
    allPs.assign(M, m, NULL); // the model indexes the rows (M), times indexes the columns (m)
    allLogPs.assign(M, m, NULL);
    
    // loop over the models
    for (int u=0; u<M; u++){
        Q = allQs[u];
        pi = piArr[u];

        status = eigendecompose(Q, pi, lambda, U, Uinv);
        if (status) {
            stringstream s;
            s << "Spectral decomposition of rate matrix for model=" << u << " failed.";
            throw ProjectError(s.str());
        }
        for (int v=0; v<m; v++){
            t = times[v]; // time

            P = expQt(t, lambda, U, Uinv);
            // store P

            allPs[u][v] = P;
            allLogPs[u][v] = log(P, states);
            //#ifdef DEBUG
            cout << "codon rate matrix log P(t=" << t << ", model=" << u << ")" << endl;
            //printCodonMatrix(allLogPs[u][v]);
            //#endif

            /*for(int i=0; i<64; i++){
              cout<<"vvv\t"<<Seq2Int(3).inv(i)<<"\t"<<omega<<"\t"<<t<<"\t"<<gsl_matrix_get(allPs[u][v],i,i)<<endl;
              }*/
        }
        gsl_matrix_free(U);
        gsl_matrix_free(Uinv);
        gsl_vector_free(lambda);
    }
}

/*
 * Returns a pointer to the 64x64 substitution probability matrix
 * for which omega and t come closest to scored values.
 */
gsl_matrix *CodonEvoDiscr::getSubMatrixLogP(double omega, double t){
    // determine index u into vector of omegas such that omegas[u] is close to omega
    //int u = findClosestIndex(omegas, omega);
    //return Evo::getSubMatrixLogP(u, t);
    return NULL;
}

/*
 * log likelihood of exon candidate pair e1, e2 (required: equal number of codons, gapless)
 * Used for testing. The pruning algorithm requires tuples of codons instead.
 */
double CodonEvoDiscr::logLik(const char *e1, const char *e2, int u, double t) {
    gsl_matrix *logP = Evo::getSubMatrixLogP(u, t);
    int n = strlen(e1)/3; // number of nucleotide triples
    int from, to; // from: codon in e1, to: codon in e2
    Seq2Int s2i(3);
    double loglik = 0.0;
    for (int i=0; i<n; i++){
	try {
	    from = s2i(e1 + 3*i);
	    to   = s2i(e2 + 3*i);
	    loglik += gsl_matrix_get(logP, from, to);
	} catch (...){} // gap or n characters
    }
    return loglik;
}

// destructor
CodonEvoDiscr::~CodonEvoDiscr(){
    for (int u=0; u<M;++u)
        if (piArr[u] != NULL)
            delete piArr[u];
    if (Theta){
        gsl_matrix_free(Theta);
        Theta = NULL;
    }
    if (bias){
        delete [] bias;
        bias = NULL;
    }
    if (w){
        delete [] w;
        w = NULL;
    }
    pi = NULL;
}

// get vector of log likelihoods of omegas for one single codon tuple
vector<double> CodonEvoDiscr::loglikForCodonTuple(vector<string> &seqtuple, PhyloTree *ctree, PhyloTree *tree, int &subs){
    if (seqtuple.size() != ctree->numSpecies())
        throw ProjectError("CodonEvoDiscr::logLikForCodonTuple: inconsistent number of species.");

    // make sure that all codons are present (strings of length 3, "---" is allowed)
    for (const string &codon : seqtuple) {
        if (codon.length() != 3)
            throw ProjectError("CodonEvoDiscr::loglikForCodonTuple: codon tuple has not length 3");
    }
    
    vector<double> logliks(M, 0.0), test_logliks(M, 0.0);
    Seq2Int s2i(3);
    vector<int> codontuple(ctree->numSpecies(), 64); // 64 = missing codon                                                     
    bit_vector bv(ctree->numSpecies(), 0);
    int numCodons = 0;
    // cout << "ctree size = " << ctree->numSpecies() << "\tcodontuple(int): ";
    for (size_t s=0; s < ctree->numSpecies(); s++){
        if (seqtuple[s].size() > 0)
            try {
                codontuple[s] = s2i(seqtuple[s].c_str());
                numCodons++;
                bv[s] = 1;
            } catch(...){} // gap or n character                                                                              
        // cout << codontuple[s] << "|";
    }

    if (numCodons >= 2){ // no computation if nothing can be compared
        if (tree){
            // cout << "### original species tree = ";
            tree->printTree();
            // calculate number of substitutions by fitch algorithm
            PhyloTree *tr;       
            auto topit = GeneMSA::topologies.find(bv);
            if (topit == GeneMSA::topologies.end()) // insert new topology
                tr = new PhyloTree(*tree);
            else
                tr  = new PhyloTree(*topit->second);
            tr->prune(bv, NULL);
            GeneMSA::topologies.insert(pair<bit_vector,PhyloTree*>(bv, tr));
            subs = tr->fitch(codontuple);
        }
        // calculate log likelihood for each of ClaMSAs rate matrix
        for (int u=0; u < M; u++){ // loop over rate matrices
            pi = piArr[u];
	    logliks[u] = ctree->pruningAlgor(codontuple, this, u);
            // cout << "loglik " << u << " = " << setprecision(6) << logliks[u] << endl;
        }
        //for (const string &codon : seqtuple)
        //    cout << "codon " << codon << endl;
        //double p = getProb(logliks);
        // cout << "clamsaProb = " << p << endl;
    }
    return logliks;
}


double CodonEvoDiscr::getProb(const vector<double> &lls){
    assert(lls.size() == M);
    double z = w[0];
    for (int i=0; i<M; ++i){
       z += w[i+1] * lls[i];
       //cout << "ll[" << i << "]=" << setprecision(5) << lls[i] << "\t";
    }
    //cout << " logit=" << z << endl;
    // apply sigmoid function
    z = 1.0/(1.0 + exp(-z));
    return z;
}

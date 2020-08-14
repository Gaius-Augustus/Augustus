/*
 * contTimeMC.cc
 *
 * License: Artistic License, see file LICENSE.TXT or 
 *          https://opensource.org/licenses/artistic-license-1.0
 */

// project includes
#include "contTimeMC.hh"
#include "geneticcode.hh"
#include "properties.hh"

// standard C/C++ includes
#include <iostream>
#include <iomanip>
#include <vector> 
#include <algorithm>

// destructor
Evo::~Evo(){
    times.clear();
    
    for (int i=0; i < allPs.getColSize(); i++)
	for (int j=0; j< allPs.getRowSize(); j++){
	    gsl_matrix_free(allPs[i][j]);
	}
    allPs.assign(0, 0, NULL);

    for (int i=0; i < allLogPs.getColSize(); i++)
    	for (int j=0; j< allLogPs.getRowSize(); j++){
    	    gsl_matrix_free(allLogPs[i][j]);
    	}
    allLogPs.assign(0, 0, NULL);

    for (size_t u=0; u < allQs.size(); u++)
        gsl_matrix_free(allQs[u]);
    allQs.assign(0, NULL);
    
    if (pi != NULL)
        delete[] pi;
}

/* Determine the branch lengths for which matrices P should be stored, at most m.
 * For trees with few species, say <10, this could be all branch lengths b occurring
 * in the tree. A value of m=-1 means to store all lengths in b. For large trees, memory 
 * is saved by storing only m times roughly representing all lengths in b.
 */
void Evo::setBranchLengths(vector<double> b, int m){
    this->m = m;
    // make a sorted copy of the unique times in b
    times = b;
    sort(times.begin(), times.end());
    vector<double>::iterator it = unique(times.begin(), times.end());
    times.resize(distance(times.begin(), it));

    if (this->m < 0 || this->m > times.size())
	this->m = times.size(); // as default take all values of b

    if (this->m < times.size()){
	// Store fewer values than are in b. 
	if (this->m == 1){
	    // take mean branch length of b if there shall be only one branch length in 'times'
	    double sum = 0.0;
	    for (vector<double>::iterator it = b.begin(); it != b.end(); it++)
		sum += *it;
	    times.clear();
	    times.push_back(sum / b.size());
	} else {
	    // For now, very simply take m equidistant branch lengths between min and max
	    double min = times.front();
	    double step = (times.back() - min) / (this->m-1);
	    times.clear();
	    for (int i=0; i < this->m; i++){
		times.push_back(min + i*step);
	    }
	}
    }
}

void Evo::printBranchLengths(){
  cout << "branch lengths: " << endl;
  for (int i=0; i<times.size(); i++)
	cout << i << " " << times[i] << endl;
}

// determines index j into vector v such that v[j] is closest to val
int Evo::findClosestIndex(vector<double> &v, double val){
    // do linear search, logarithmic search is probably not worth it with the vector sizes we expect
    int j = 0;
    double dist = 1000; // infinity
    for (int i=0; i < v.size(); i++)
	if (abs(v[i]-val) < dist){
	    dist = abs(v[i]-val);
	    j = i;
	}
    return j;
}


gsl_matrix *Evo::getSubMatrixP(int u, double t){
    // determine index v into vector of branch lengths such that times[v] is close to t
    int v = findClosestIndex(times, t);
    return allPs[u][v];
}

gsl_matrix *Evo::getSubMatrixLogP(int u, double t){
    // determine index v into vector of branch lengths such that times[v] is close to t
    int v = findClosestIndex(times, t);
    return allLogPs[u][v];
}


/*
 * perform a decompososition of the rate matrix as Q = U * diag(lambda) * U^{-1}
 */
int eigendecompose(gsl_matrix *Q,        // input
		    double *pi,          // input, must be same as in construction of Q
		    gsl_vector *&lambda, // output, memory for lambda, U and Uinv is allocated within the function
		    gsl_matrix *&U,      // output
		    gsl_matrix *&Uinv){  // output
    
    // Codon usages pi[j] that vanish are a problem because of division by zero.
    // Therefore, for the purpose of this decomposition, treat any such pi[j] as if 1 instead. This is 
    // not a problem as the j-th column and row of Q is a zero column anyways.
    
    // Compute B = diag(pi^{1/2}) * Q * diag(pi^{-1/2}) which is symmetrical, as Q is time-reversible
    gsl_matrix* B = gsl_matrix_calloc (64, 64);
    for (int i=0; i<64; i++){
	for (int j=0; j<64; j++){
	    if (pi[i] > 0.0 && pi[j] > 0.0){
		double bij = gsl_matrix_get(Q, i, j) * sqrt(pi[i]) / sqrt(pi[j]); // Q vanishes otherwise anyways
		gsl_matrix_set(B, i, j, bij);
	    }
	}
    }

    gsl_eigen_symmv_workspace *w = gsl_eigen_symmv_alloc (64); // work space
    lambda = gsl_vector_calloc (64);
    U = gsl_matrix_calloc (64, 64);
    Uinv = gsl_matrix_calloc (64, 64);
    int status = gsl_eigen_symmv (B, lambda, U, w); // now: B = U * diag(lambda) * U^t, for orthogonal U, i.e. U^t = U^{-1}
    if (status)
	return status;
    gsl_eigen_symmv_free(w);
    gsl_matrix_free(B);

    // reverse transformation to Q = diag(pi^{-1/2}) * B * diag(pi^{1/2}) = [diag(pi^{-1/2}) * U] * diag(lambda) * [U^t * diag(pi^{1/2})]
    // reuse space of U
    // set Uinv <- U^t * diag(pi^{1/2})
    for (int j=0; j<64; j++){
	double factor = (pi[j]>0.0)? sqrt(pi[j]) : 1;
	for (int i=0; i<64; i++){
	    gsl_matrix_set(Uinv, i, j, factor * gsl_matrix_get(U, j, i)); // (this U is orthogonal)
	}
    }
    // and U <- diag(pi^{-1/2}) * U
    for (int i=0; i<64; i++){
	double factor = (pi[i] > 0.0)? 1/sqrt(pi[i]) : 1;
	for (int j=0; j<64; j++){
	    gsl_matrix_set(U, i, j, factor * gsl_matrix_get(U, i, j));
	}
    }
    //validateEigenDecomp(Q,lambda,U,Uinv,64);
    return 0;
}

// to test whether the decomposition is correct, look at whether Q - U * diag(lambda) * Uinv is close to the 0-matrix
void validateEigenDecomp(gsl_matrix *Q, gsl_vector *lambda, gsl_matrix *U, gsl_matrix *Uinv, int states){
    int N = states;
    cout << "residue Q - U * diag(lambda) * Uinv" << endl;
    for (int i=0; i<N; i++){
	for (int j=0; j<N; j++){
	    double residue = gsl_matrix_get(Q, i, j);
	    for (int k=0; k<N; k++){
		residue -= gsl_matrix_get(U, i, k) * gsl_vector_get(lambda, k) * gsl_matrix_get(Uinv, k, j);
	    }
	    cout << setw(6) << fixed << setprecision(3) << residue << " ";
	}
	cout << endl;
    }
}

/*
 * compute the matrix exponential P(t) = exp(Q*t),
 * where Q is given by U, diag(lambda) and U^{-1} as computed by eigendecompose
 * Values of P are stored in log-space, i.e. ln P is stored (elementwise natural logarithm).
 */
gsl_matrix *Evo::expQt(double t, gsl_vector *lambda, gsl_matrix *U, gsl_matrix *Uinv){

    int N = states;

    gsl_matrix* P = gsl_matrix_alloc (N, N);    
    double Pij;
    // P(t) = exp(Q*t) = U * diag(exp(lambda_1 * t), ..., exp(lambda_n * t)) * Uinv
    for (int i=0; i<N; i++){
	for (int j=0; j<N; j++){
	    Pij = 0;
	    for (int k=0; k<N; k++)
		Pij += gsl_matrix_get(U, i, k) * exp(gsl_vector_get(lambda, k)*t) * gsl_matrix_get(Uinv, k, j);
	    if (Pij>0.0)
		gsl_matrix_set(P, i, j, Pij);
	    else if (Pij < -0.001){
		cerr << Pij << endl;
		throw ProjectError("negative probability in substitution matrix");
	    } else {
		gsl_matrix_set(P, i, j, 0.0);
	    }
	}
    }
    if(states == 3 || states == 4){
	gsl_matrix_set(P, 0, 0, gsl_matrix_get(P, 0, 0));
	gsl_matrix_set(P, 0, 2, gsl_matrix_get(P, 0, 0));     
    }
    if(states == 4){
	gsl_matrix_set(P, 0, 3, gsl_matrix_get(P, 0, 2));     
	gsl_matrix_set(P, 1, 3, gsl_matrix_get(P, 1, 2));
	gsl_matrix_set(P, 2, 3, 1);
    }
    return P;
}

/*
 * print a 64x64 codon matrix with neat row and column labels by codon and amino acid
 */
void printCodonMatrix(gsl_matrix *M) {    // M is usually Q or P = exp(Qt)
    cout << "      ";
    for (int j=0; j<64; j++)
	cout << "   " << Seq2Int(3).inv(j) << "  "; // head line codon
    cout << endl << "      ";;
    for (int j=0; j<64; j++)
	cout << "   " << GeneticCode::translate(j) << "    "; // head line amino acid
    cout << endl;
    for (int i=0; i<64; i++){
	cout << Seq2Int(3).inv(i) << " " << GeneticCode::translate(i) << " ";
	for (int j=0; j<64; j++){
	    double m = gsl_matrix_get(M,i,j);
	    if ( 0 <= m && m < 1e-6)
		cout << setw(7) << " "; // probability 0
	    else if (m == -numeric_limits<double>::max())
		cout << setw(8) << " -inf "; // probability 0 in log-space
	    else
		cout << setw(7) << fixed << setprecision(3) << gsl_matrix_get(M,i,j) << " ";
	}
	cout << endl;
    }
}

/*
 * computes the element-wise natrual logarithm of a matrix P
 */
gsl_matrix *log(gsl_matrix *P, int states){

    gsl_matrix* LogP = gsl_matrix_calloc(states, states); 
    for (int i=0; i<states; i++){
     	for (int j=0; j<states; j++){
     	    double Pij = gsl_matrix_get(P,i,j);
    	    if (Pij>0.0)
    		gsl_matrix_set(LogP, i, j, log(Pij));
    	    else if (Pij < -0.001){
    		cerr << Pij << endl;
    		throw ProjectError("negative probability in substitution matrix");
    	    } else {
    		gsl_matrix_set(LogP, i, j, -numeric_limits<double>::max()); // ln 0 = -infty       
    	    }
     	}	
    }
    return LogP;
}

void ExonEvo::setLambda(double lambda){
    try {
	this->lambda = Properties::getdoubleProperty("/CompPred/exon_gain");
    } catch (...) {
        this->lambda = lambda;
    }
    
    if (this->lambda <= 0.0)
	throw ProjectError("the rates for exon loss/gain have to be positive");
}

void ExonEvo::setMu(double mu){
    try {
	this->mu  = Properties::getdoubleProperty("/CompPred/exon_loss");
    } catch (...) {
        this->mu = mu;
    }

    if (this->mu <= 0.0)
	throw ProjectError("the rates for exon loss/gain have to be positive");
}

void ExonEvo::setAliErr(double ali_error){
    try {
	this->ali_error = Properties::getdoubleProperty("/CompPred/ali_error");
    } catch (...) {
        this->ali_error = ali_error;
    }

    if (this->ali_error <= 0.0)
	throw ProjectError("the rate for alignment errors has to be positive");
}


void ExonEvo::getRateMatrices(){
    allQs.assign(1, getExonRateMatrix());
}


void ExonEvo::computeLogPmatrices(){
    gsl_matrix *Q = allQs[0];

    // decompose rate matrix Q = U * diag(l_1,...,l_n) * Uinv
    int status = eigendecompose(Q);
    if (status) {
        stringstream s;
        s << "Spectral decomposition of exon rate matrix for failed.";
        throw ProjectError(s.str());
    }

    // for all branch lengths t, compute transition matrices P(t) and logP(t)
    allPs.assign(1, m, NULL);
    allLogPs.assign(1, m, NULL);
    for (int v=0; v<m; v++){
        double t = times[v]; // time
        gsl_matrix *P = expQt(t);

        /*
          cout << "printing P for t=" << t << endl;
          for (int i = 0; i < states; ++i){
          for(int j = 0 ; j < states; j++){
          cout << gsl_matrix_get (P, i, j) << " ";
          }
          cout << endl;
          }
        */
        allPs[0][v] = P;
        allLogPs[0][v] = log(P, states);
    }
}

void ExonEvo::setPi(){
    if (states == 2){
	this->pi = new double[2];
	this->pi[0] = (mu / (lambda + mu));
	this->pi[1] = (lambda / (lambda + mu));
    }
    else if (states == 3){
	this->pi = new double[3];
	this->pi[0] = (mu / (lambda + mu));
        this->pi[1] = (lambda / (lambda + mu));
	this->pi[2] = 0.0;
	/*double denominator = 3.0 * (lambda + mu + ali_error);
        this->pi[0] = (ali_error + (2*mu)) / denominator;
        this->pi[1] = (ali_error + (2*lambda)) / denominator;
	this->pi[2] = 1.0 / 3.0;*/
    }
    else if (states == 4){
	this->pi = new double[4];
	this->pi[0] = (mu / (lambda + mu));
        this->pi[1] = (lambda / (lambda + mu));
	this->pi[2] = 0.0;
	this->pi[3] = 0.0;
	/*double denominator = 3.0 * (lambda + mu + ali_error);
        this->pi[0] = (ali_error + (2*mu)) / denominator;
        this->pi[1] = (ali_error + (2*lambda)) / denominator;
	this->pi[2] = 1.0 / 3.0;*/
    }
    else{
	throw ProjectError("internal error: wrong number of states in ExonEvo model: choose either 2,3 or 4.");
    }
}  

void ExonEvo::addBranchLength(double b){

    bool isIncluded = false;

    for (int i=0; i < times.size(); i++){
	if(times[i] == b){
	    isIncluded = true;
	    break;
	}
    }
    if(!isIncluded){ // add branch length
	times.push_back(b);
	gsl_matrix *P= expQt(b);
	vector<gsl_matrix*> row = allPs.getRow(0);
	row.push_back(P);
	allPs.assign(1, row.size(), NULL);
	for (int v=0; v<row.size(); v++){
	    allPs[0][v]=row[v];	
	}
	gsl_matrix *LogP=log(P,states);
	row = allLogPs.getRow(0);
	row.push_back(LogP);
	allLogPs.assign(1, row.size(), NULL);
	for (int v=0; v<row.size(); v++){
	    allLogPs[0][v]=row[v];	
	}
    }
}

gsl_matrix *ExonEvo::getExonRateMatrix(){
    const int N = states;

    gsl_matrix *Q = gsl_matrix_alloc (N, N);

    gsl_matrix_set (Q, 0, 1, lambda); // exon gain with rate lambda
    gsl_matrix_set (Q, 1, 0, mu);     // exon log with rate mu
    if(N > 2){                        // rate for alignment errors: ali_error >> mu,lambda 
	gsl_matrix_set (Q, 0, 2, ali_error);
	gsl_matrix_set (Q, 1, 2, ali_error);
	gsl_matrix_set (Q, 2, 0, 0);
	gsl_matrix_set (Q, 2, 1, 0);
    }
    if(N > 3){
	gsl_matrix_set (Q, 0, 3, 0);
	gsl_matrix_set (Q, 1, 3, 0);
	gsl_matrix_set (Q, 2, 3, 0);
	gsl_matrix_set (Q, 3, 0, 0);
	gsl_matrix_set (Q, 3, 1, 0);
	gsl_matrix_set (Q, 3, 2, 0);
    }
    // set diagonal elements to the negative sum of the other rates in the row
    for (int i = 0; i < N; ++i){
	double Pii = 0;
	for(int j = 0 ; j < N; j++){
	    if(i != j)
		Pii -= gsl_matrix_get (Q, i, j);
	}
	gsl_matrix_set (Q, i, i, Pii);
    }
    /*cout << "printing Q" << endl;
    for (int i = 0; i < N; ++i){
	for(int j = 0 ; j < N; j++){
	    cout << gsl_matrix_get (Q, i, j) << " ";
	}
	cout << endl;
	}*/
    return Q;
}
    

int ExonEvo::eigendecompose(gsl_matrix *Q){

    const int N = states; // dimension of matrices: N x N

    // allocate memory for temporary matrices/vectors ...
    gsl_vector_complex *l_complex = gsl_vector_complex_alloc (N);    // complex vector of eigenvalues
    gsl_matrix_complex *U_complex = gsl_matrix_complex_alloc (N, N); // complex matrix of eigenvectors 
    gsl_matrix *U_LU_decomp = gsl_matrix_alloc (N, N);               // LU decomposition of matrix U (required for calcuating Uinv)
    gsl_eigen_nonsymmv_workspace * w = gsl_eigen_nonsymmv_alloc (N);

    // ... and permanent stuff: Q = U * diag(l_1,...,l_n) * Uinv 
    l = gsl_vector_alloc (N);                                        // real vector of eigenvalues
    U = gsl_matrix_alloc (N, N);                                     // real matrix of eigenvectors
    Uinv = gsl_matrix_alloc (N, N);                                  // inverse of U

    int s; // signum for LU decomposition
    gsl_permutation * perm = gsl_permutation_alloc (N); // permutation matrix of LU decomposition
    
    int status = gsl_eigen_nonsymmv (Q, l_complex, U_complex, w); // eigen decomposition
    if (status)
	return status;
    
    // check if all eigenvalues and eigenvectors have real values
    for (int i = 0; i < N; ++i){
	gsl_complex z = gsl_vector_complex_get (l_complex, i);
	if(GSL_IMAG(z) != 0)
	    throw ProjectError("Matrix diagonalization: complex eigen values.");
	gsl_vector_set (l, i, GSL_REAL(z));
        for(int j = 0 ; j < N; j++){
	    z = gsl_matrix_complex_get(U_complex, i, j);
	    if(GSL_IMAG(z) != 0)
		throw ProjectError("Matrix diagonalization: complex eigen vector.");
	    gsl_matrix_set (U, i, j, GSL_REAL(z)); 
	    gsl_matrix_set (U_LU_decomp, i, j, GSL_REAL(z));
	}
    }

    status = gsl_linalg_LU_decomp (U_LU_decomp, perm, &s); // LU decomposition (required for calculating Uinv)
    if (status) 
        return status;

    // Invert matrix U_LU_decomp
    status = gsl_linalg_LU_invert (U_LU_decomp, perm, Uinv);
    if (status)
        return status;

    // free memory of temporary stuff
    gsl_permutation_free (perm);
    gsl_vector_complex_free(l_complex);
    gsl_matrix_complex_free(U_complex);
    gsl_matrix_free(U_LU_decomp);
    gsl_eigen_nonsymmv_free (w);

    // check if decomposition is correct
    //validateEigenDecomp(Q);

    return 0;
}

void Parsimony::computeLogPmatrices(){

    allLogPs.assign(1, m, NULL);
    gsl_matrix* P = gsl_matrix_calloc(states, states);
    for (int i=0; i<states; i++){
	for (int j=0; j<states; j++){
	    if(i != j){
		gsl_matrix_set(P, i, j, -1);
	    }
	    else{
		gsl_matrix_set(P, i, j, 0);
	    }
	}
    }
    allLogPs[0][0]=P;
}

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
#include "phylotree.hh"
#include "geneMSA.hh"
#include "igenicmodel.hh"

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
    if(pi != NULL)
      delete[] pi;
}

/* Determine the branch lengths for which matrices P should be stored, at most m.
 * For trees with few species, say <10, this could be all branch lengths b occuring
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

void CodonEvo::setPi(double *pi){
    this->pi = new double[64];
    for (int i=0;i<64; i++)
	this->pi[i] = pi[i];
}

/*
 * Chooses k values for omega around 1, for which matrices P will be stored.
 * Examples:
 * k   set of omegas
 * 1   1 (does not make sense for selection detection)
 * 2   0.67, 1
 * 3   0.67, 1, 1.5
 * 4   0.5, 0.75, 1, 1.33
 * 5   0.5, 0.75, 1, 1.333, 2
 * 6   0.4, 0.6, 0.8, 1, 1.25, 1.67
 * 10  0.29, 0.43, 0.57, 0.71, 0.86, 1, 1.17, 1.4, 1.75, 2.33
 * Later, one can see which of these values of omega gives the maximum likelihood.
 */
void CodonEvo::setOmegas(int k){
    omegas.clear();
    int c = 2; // must be >0, the larger this constant the closer the omegas are to 1
    if (k<1)
	return;
    int r = (int) (0.5 + (double) (k-1)/2);
    int rr = k-1-r;
    for (int i=r; i>=1; i--)
	omegas.push_back(1 - (double) i/(c+r));
    omegas.push_back(1);
    for (int i=1; i<=rr; i++)
	omegas.push_back(1.0 / (1 - (double) i/(c+r)));
    if(Constant::useNonCodingModel){
      this->k = k+1;
      omegas.push_back(0);
    }
    else
      this->k = k;
    /*cout<<"omega: ";
    for(int j=0; j<omegas.size(); j++){
      cout<<omegas[j]<<"\t";
    }
    cout<<endl;
    */
}

/*
 * use a normal distribution with mean 1 and standard deviation sigma as prior for omega
 */
void CodonEvo::setPrior(double sigma ){
    omegaPrior.resize(k, 0);
    double sum = 0;
    for (int i=0; i<k; i++){
	double omega = omegas[i];
	if (omega > 1.0){
	  //omega = 1.0 / omega;
	  omega = sqrt(-2*log(1/omega) * pow(sigma,2) + pow(1/omega - 1,2)) + 1;
	}
	double t = (omega-1.0)/sigma;
	sum += omegaPrior[i] = exp(-t*t/2);
    }
    // normalize
    for (int i=0; i<k; i++)
	omegaPrior[i] /= sum;
}

void CodonEvo::printOmegas(){
  cout << "index\tomega\tprior" << endl;
    for (int i=0; i < omegas.size(); i++)
      cout << i << "\t" << omegas[i] << "\t" << omegaPrior[i] << endl;
}

/*
 * Read BLOSUM from file. BLOSUM file contains joint probabilities of pairs.
 * Use them to calculate posterior probabilities of AA substitutions.
 */

void CodonEvo::setAAPostProbs(){
  aaUsage.resize(20,0);
  vector<vector<double> > b(20, vector<double>(20));
  ifstream blosumfile;
  string filename = Constant::configPath + "profile/default.qij";
  string buffer = "";
  //char aaNames[20] = {""};
  int aaIndex[20] = {-1};
  char aa = 'X';
  blosumfile.open(filename.c_str(), ifstream::in);
  if(!blosumfile){
    string errmsg = "Could not open the BLOSUM file " + filename + ".";
    throw PropertiesError(errmsg);
  }
  while(blosumfile >> buffer && buffer.substr(0,1) == "#"){
    blosumfile.ignore( numeric_limits<streamsize>::max(), '\n');
  }
  //aaNames[GeneticCode::get_aa_from_symbol(buffer[0])] = buffer[0];
  aaIndex[0] = GeneticCode::get_aa_from_symbol(buffer[0]);
  for(int i=1; i<20; i++){
    blosumfile >> aa;
    //aaNames[GeneticCode::get_aa_from_symbol(aa)] = aa;
    aaIndex[i] = GeneticCode::get_aa_from_symbol(aa);
  }
  
  for(int i = 0; i < 20; i++){
    for(int j = 0; j <= i; j++){
      blosumfile >> b[aaIndex[i]][aaIndex[j]];
    }
  }
  for(int i = 0; i < 20; i++){
    for(int j = i+1; j < 20; j++){
      b[aaIndex[i]][aaIndex[j]] = b[aaIndex[j]][aaIndex[i]];
    }
  }

  for(int i = 0; i < 20; i++){
    for(int j = 0; j < 20; j++){
      aaUsage[aaIndex[i]] += b[aaIndex[j]][aaIndex[i]];  // set aaUsage 
    }
  }
   
  /***
  // print Matrix
  cout << "blosum: " << endl << "\t";
  for (int i = 0; i < 20; i++)
    cout << aaNames[i] << "\t";
  cout << endl;
  for(int i = 0; i < 20; i++) {
    cout << aaNames[i] << "\t";                                                                                                   
    for (int j = 0; j < 20; j++) {
      cout << b[i][j] << "\t";
    }
    cout << endl;
  }
  ***/
  
  // compute posterior probabilities: P(i|j) = P(i,j)/P(j)
  for(int i = 0; i < 20; i++){
    for(int j = 0; j < 20; j++){
      b[i][j] = (b[i][j] / aaUsage[j]); 
    }
  }

  aaPostProb = b;
}



/* 
 * Precompute and store the array of matrices.
 * Takes approximate k * m * 0.015 seconds on greif1: Time for eigendecompose is small compared to time for expQt.
 */
void CodonEvo::computeLogPmatrices(){
    double omega;
    double t;
    gsl_matrix *Q, *U, *Uinv, *P;
    gsl_vector *lambda;
    int status;
    allPs.assign(k, m, NULL); // omegas index the rows, times index the columns
    allLogPs.assign(k, m, NULL);
    for (int u=0; u<k; u++){
	omega = omegas[u];
	// compute decomposition of Q, which does not require t yet
	if(u < k){
	  if(!aaPostProb.empty()) // use amino acid score in codon rate matrix 
	    Q = getCodonRateMatrix(pi, omega, kappa, &aaPostProb);
	  else
	    Q = getCodonRateMatrix(pi, omega, kappa);
	}else{
	  vector<double> pi_nuc = IGenicModel::getNucleotideProbs();
	  Q = getNonCodingRateMatrix(&pi_nuc, kappa);
	}
	status = eigendecompose(Q, pi, lambda, U, Uinv);
	if (status) {
	    stringstream s;
	    s << "Spectral decomposition of rate matrix for omega=" << omega << " failed.";
	    throw ProjectError(s.str());
	}
	for (int v=0; v<m; v++){
	    t = times[v]; // time
	    P = expQt(t, lambda, U, Uinv);
	    // store P
	    allPs[u][v] = P;
	    allLogPs[u][v] = log(P,states);
	    //#ifdef DEBUG
	    //cout << "codon rate matrix log P(t=" << t << ", omega=" << omega << ")" << endl;
	    //printCodonMatrix(allLogPs[u][v]);
	    //#endif

	    /*for(int i=0; i<64; i++){
	       cout<<"vvv\t"<<Seq2Int(3).inv(i)<<"\t"<<omega<<"\t"<<t<<"\t"<<gsl_matrix_get(allPs[u][v],i,i)<<endl;
	       }*/

	}
	gsl_matrix_free(U);
	gsl_matrix_free(Uinv);
	gsl_matrix_free(Q);
	gsl_vector_free(lambda);
    }
}

/*
 * Returns a pointer to the 64x64 substitution probability matrix
 * for which omega and t come closest to scored values.
 */
gsl_matrix *CodonEvo::getSubMatrixLogP(double omega, double t){
    // determine index u into vector of omegas such that omegas[u] is close to omega
    int u = findClosestIndex(omegas, omega);
    return Evo::getSubMatrixLogP(u, t);
}

/*
 * log likelihood of exon candidate pair e1, e2 (required: equal number of codons, gapless)
 * Used for testing. The pruning algorithm requires tuples of codons instead.
 */
double CodonEvo::logLik(const char *e1, const char *e2, int u, double t) {
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

double CodonEvo::estOmegaOnSeqPair(const char *e1, const char *e2, double t, // input variables
		int& numAliCodons, int &numSynSubst, int &numNonSynSubst){ // output variables
    if (!e1 || !e2 || strlen(e1) != strlen(e2) || strlen(e1)%3 != 0)
	throw ProjectError("CodonEvo::logLik: wrong exon lengths");
    int maxU = 0; // index to omegas
    double loglik, ML = -numeric_limits<double>::max();
    for (int u=0; u < k; u++){
	loglik = logLik(e1, e2, u, t);
	if (loglik > ML){
	    ML = loglik;
	    maxU = u;
	}
    }
    
    numAliCodons = 0; // number of aligned codons (gapless in both exons)
    numSynSubst = 0;
    numNonSynSubst = 0;
    int n = strlen(e1)/3; // number of nucleotide triples
    int from, to; // from: codon in e1, to: codon in e2
    Seq2Int s2i(3);
    for (int i=0; i<n; i++){
	try {
	    from = s2i(e1 + 3*i);
	    to   = s2i(e2 + 3*i);
	    numAliCodons++;
	    if (from != to){
		if (GeneticCode::translate(from) == GeneticCode::translate(to))
		    numSynSubst++;
		else 
		    numNonSynSubst++;
	    }
	} catch (...){} // gap or n characters
    }
    return omegas[maxU];
}

// destructor
CodonEvo::~CodonEvo(){
    omegas.clear();
}

/*
 * 64x64 rate matrix Q for codon substitutions as in Yang, "Computational Molecular Evolution", (2.7)
 */

gsl_matrix *getCodonRateMatrix(double *pi,    // codon usage, normalized vector with 64 elements
			       double omega,  // dN/dS, nonsynonymous/synonymous ratio
			       double kappa, // transition/transversion ratio, usually >1
			       vector<vector<double> > *aaPostProb){ // posterior probabilities of aa substitutions
    gsl_matrix* Q = gsl_matrix_calloc (64, 64); // initialize Q as the zero-matrix
    Seq2Int s2i(3);

    // store and print number of synonymous and nonsynonymous substitutions for each codon                                           
    vector<int> num_syn(64,0);
    vector<int> num_nonsyn(64,0);
    vector<int> num_wk(64,0);

    // numerator and denominator of lambda, the normalizing factor for using amino acid scores in codon rate matrix
    vector<double> lambda1(64,0);
    vector<double> lambda2(64,0);

    // initialize off-diagonal elements
    // codon i (row)    abc
    // codon j (col)    dbc
    // f                012
    int codoni[3], codonj[3];
    for (int a=0; a<4; a++){
	codoni[0] = a;
	for (int b=0; b<4; b++){
	    codoni[1] = b;
	    for (int c=0; c<4; c++){
		codoni[2] = c;
		int i = 16*codoni[0] + 4*codoni[1] + codoni[2]; // in {0,...63}
		// go through the 3x3 codons that differ from codon i by a single mutation
		for (int f=0;f<3; f++){ // position of mutation
		    for (int d=0; d<4; d++){ // d = substituted base
			if (d != codoni[f]){
			    codonj[0] = codoni[0];
			    codonj[1] = codoni[1];
			    codonj[2] = codoni[2];
			    codonj[f] = d;
			    int j = 16*codonj[0] + 4*codonj[1] + codonj[2]; // in {0,...63}
			    double qij = pi[j];
			    if (!(pi[i] > 0.0))
				qij = 0.0; // rows corresponding to stop codons are also set to 0
			    if (GeneticCode::is_purine(codonj[f]) == GeneticCode::is_purine(codoni[f])){
				qij *= kappa; // transition
			    }
			    if (GeneticCode::translate(i) != GeneticCode::translate(j)){
			      qij *= omega; // non-synonymous subst.
			      if(aaPostProb != NULL && GeneticCode::map[i] != -1 && GeneticCode::map[j] != -1){
				//qij *= (*aaPostProb)[GeneticCode::map[i]][GeneticCode::map[j]];
				lambda1[i] += pi[j];  
				lambda2[i] += pi[j] * (*aaPostProb)[GeneticCode::map[i]][GeneticCode::map[j]];
			      }
			      num_nonsyn[i]++;
			    }else{
			      num_syn[i]++;
			    }
			    if(GeneticCode::is_purine(codonj[f]) == GeneticCode::is_purine(codoni[f]) && GeneticCode::translate(i) != GeneticCode::translate(j)){
			      num_wk[i]++;
			    }
			    gsl_matrix_set (Q, i, j, qij);
			}
		    }
		}
	    }
	}
    }
    /*
    cout<<"------------ number of synonymous and nonsynonymous substitutions for each codon ----------------"<<endl;
    cout<<"codon\tnum_syn\tnum_nonsyn\tnum_wk\tcodon_usage"<<endl;
    for(int i=0; i<64; i++){
      cout<<Seq2Int(3).inv(i)<<"\t"<<num_syn[i]<<"\t"<<num_nonsyn[i]<<"\t"<<num_wk[i]<<"\t"<<pi[i]<<endl;
    }
    */

    // add posterior prob. of AAs and lambda to non-syn. substitutions if animo acid score is used
    if(aaPostProb != NULL){
      double qij = 0;
      for (int i=0; i<64; i++){
	for (int j=0; j<64; j++){
          if (GeneticCode::translate(i) != GeneticCode::translate(j) && GeneticCode::map[i] != -1 && GeneticCode::map[j] != -1){
            qij = gsl_matrix_get(Q, i, j);
    	    qij *= (*aaPostProb)[GeneticCode::map[i]][GeneticCode::map[j]] * lambda1[i]/lambda2[i];
    	    gsl_matrix_set(Q, i, j, qij);
          }
        }
      }
    }

    // initialize diagonal elements (each row must sum to 0)
    double rowsum;
    double scale = 0.0; // factor by which matrix must be scaled to achieve one expected mutation per time unit
    for (int i=0; i<64; i++){
	rowsum = 0.0;
	for (int j=0; j<64; j++)
	    rowsum += gsl_matrix_get(Q, i, j);
	gsl_matrix_set(Q, i, i, -rowsum); // q_{i,i} = - sum_{j!=i} q_{i,j}
	scale += rowsum * pi[i];
    }
    if (scale != 0.0)
	scale = 1.0/scale;
    else
	scale = 1.0; // should not happen
    // normalize Q so that the exptected number of codon mutations per time unit is 1
    gsl_matrix_scale(Q, scale);
    return Q;
}

/*
 * compute rate matrix for non-codong model
 */

gsl_matrix *getNonCodingRateMatrix(vector<double> *pi_nuc, double kappa){

  gsl_matrix* Q = gsl_matrix_calloc (64, 64); // initialize Q as the zero-matrix
  Seq2Int s2i(3);
  vector<double> pi(64,0);

  // initialize off-diagonal elements                                                                                               
  // codon i (row)    abc                                                                                                           
  // codon j (col)    dbc                                                                                                           
  // f                012                                                                                                           
  int codoni[3], codonj[3];
  for (int a=0; a<4; a++){
    codoni[0] = a;
    for (int b=0; b<4; b++){
      codoni[1] = b;
      for (int c=0; c<4; c++){
	codoni[2] = c;
	int i = 16*codoni[0] + 4*codoni[1] + codoni[2]; // in {0,...63}                                                       
	// go through the 3x3 codons that differ from codon i by a single mutation
	for (int f=0;f<3; f++){ // position of mutation
	  for (int d=0; d<4; d++){ // d = substituted base
	    if (d != codoni[f]){
	      codonj[0] = codoni[0];
	      codonj[1] = codoni[1];
	      codonj[2] = codoni[2];
	      codonj[f] = d;
	      int j = 16*codonj[0] + 4*codonj[1] + codonj[2]; // in {0,...63}
	      double qij = (*pi_nuc)[codonj[f]];
	      pi[i] = (*pi_nuc)[codonj[f]];
	      if (GeneticCode::is_purine(codonj[f]) == GeneticCode::is_purine(codoni[f])){
		qij *= kappa; // transition
	      }
	      gsl_matrix_set (Q, i, j, qij);
	    }
	  }
	}
      }
    }
  }

  // initialize diagonal elements (each row must sum to 0)
  double rowsum;
  double scale = 0.0; // factor by which matrix must be scaled to achieve one expected mutation per time unit                       
  for (int i=0; i<64; i++){
    rowsum = 0.0;
    for (int j=0; j<64; j++)
      rowsum += gsl_matrix_get(Q, i, j);
    gsl_matrix_set(Q, i, i, -rowsum); // q_{i,i} = - sum_{j!=i} q_{i,j}                                                           
    scale += rowsum * pi[i];
  }
  if (scale != 0.0)
    scale = 1.0/scale;
  else
    scale = 1.0; // should not happen                                                                                             
  // normalize Q so that the exptected number of codon mutations per time unit is 1                                                 
  gsl_matrix_scale(Q, scale);

  return Q;

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

void ExonEvo::setLambda(){

    try {
	lambda = Properties::getdoubleProperty("/CompPred/exon_gain");
    } catch (...) {
	lambda = 0.0001;
    }
    if(lambda <= 0.0){
	throw ProjectError("the rates for exon loss/gain have to be positive");
    }
}

void ExonEvo::setMu(){
  
    try {
	mu  = Properties::getdoubleProperty("/CompPred/exon_loss");
    } catch (...) {
	mu  = 0.0001;
    }
    if(mu <= 0.0){
	throw ProjectError("the rates for exon loss/gain have to be positive");
    }
}

void ExonEvo::setAliErr(){
  
    try {
	ali_error = Properties::getdoubleProperty("/CompPred/ali_error");
    } catch (...) {
	ali_error  = 0.1;
    }
    if(ali_error <= 0.0){
	throw ProjectError("the rate for alignment errors has to be positive");
    }
}

void ExonEvo::computeLogPmatrices(){

    gsl_matrix *Q = getExonRateMatrix();

    // decompose rate matrix Q = U * diag(l_1,...,l_n) * Uinv
    int status = eigendecompose(Q);
    if (status) {
	stringstream s;
	s << "Spectral decomposition of exon rate matrix for failed.";
	throw ProjectError(s.str());
    }
    gsl_matrix_free(Q);

    // for all branch lengths t, compute transition matrices P(t) and logP(t)
    allPs.assign(1, m, NULL);
    allLogPs.assign(1, m, NULL);
    for (int v=0; v<m; v++){
	double t = times[v]; // time
	gsl_matrix *P = expQt(t);
	/*cout << "printing P for t=" << t << endl;
	for (int i = 0; i < states; ++i){
	    for(int j = 0 ; j < states; j++){
		cout << gsl_matrix_get (P, i, j) << " ";
	    }
	    cout << endl;
	    }*/
        allPs[0][v]=P;
	allLogPs[0][v]=log(P,states); 
    }
}

void ExonEvo::setPi(){

    if(states == 2){
	this->pi = new double[2];
	this->pi[0] = (mu / (lambda + mu));
	this->pi[1] = (lambda / (lambda + mu));
    }
    else if(states == 3){
	this->pi = new double[3];
	this->pi[0] = (mu / (lambda + mu));
        this->pi[1] = (lambda / (lambda + mu));
	this->pi[2] = 0.0;
	/*double denominator = 3.0 * (lambda + mu + ali_error);
        this->pi[0] = (ali_error + (2*mu)) / denominator;
        this->pi[1] = (ali_error + (2*lambda)) / denominator;
	this->pi[2] = 1.0 / 3.0;*/
    }
    else if(states == 4){
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

/* 
 * Estimate omega on a sequence of codon tuples.
 * only for testing, may need adjustments
 * TODO:
 * - change data structure for codon tuples to vectors of integers, more efficient
 * - save log-likelihood of codon tuples in a cache, so that calculation has to be done only once.
 * - scale branch lengths
 */
double CodonEvo::estOmegaOnSeqTuple(vector<string> &seqtuple, PhyloTree *tree,
				    int &subst){ //output variables
    if (seqtuple.size() != tree->numSpecies())
	throw ProjectError("CodonEvo::estOmegaOnSeqTuple: inconsistent number of species.");
    for(int i=1; i<seqtuple.size();i++){
	if(seqtuple[0].length() != seqtuple[i].length()){
	    throw ProjectError("CodonEvo::estOmegaOnSeqTuple: wrong exon lengths");
	}
    }

    int n = seqtuple[0].length()/3; // number of nucleotide triples
    int maxU = 0; // index to omegas
    double ML = -numeric_limits<double>::max();
    Evo* evo = this;
    for (int u=0; u < k; u++){ // loop over omegas
	Seq2Int s2i(3);
	double loglik = 0.0;
	for (int i=0; i<n; i++){
	    vector<int> codontuple(tree->numSpecies(), 64); // 64 = missing codon
	    int numCodons = 0;
	    for(size_t s=0; s < tree->numSpecies(); s++){
		if (seqtuple[s].size()>0)
		    try {
			codontuple[s] = s2i(seqtuple[s].c_str() + 3*i);
			numCodons++;
		    } catch(...){} // gap or n character
	    }
	    if (numCodons >= 2){
		loglik += tree->pruningAlgor(codontuple, evo, u);
	    }
	}
	//	cout << "loglikelihood(omega=" << omegas[u] << ")= " << setPrecision(4) << loglik << endl;
	if (loglik > ML){
	    ML = loglik;
	    maxU = u;
	}
    }
    // count number of substitutions
    subst = 0;
    Seq2Int s2i(3);

    // settings to reduce MAP algorithm to Fitch Algorithm
    vector<double> weights(tree->numSpecies(),0);
    Parsimony parsi;
    parsi.computeLogPmatrices();
    Evo *parsi_base = &parsi;
    
    for (int i=0; i<n; i++){
	vector<int> codontuple(tree->numSpecies(),64);
	int numCodons=0;
	for(size_t s=0; s < tree->numSpecies(); s++){
	    if (seqtuple[s].size() > 0)
		try {
		    codontuple[s] = s2i(seqtuple[s].c_str() + 3*i);
		    numCodons++;
		} catch(...){} // gap or n character
	}
	if(numCodons >= 2){
	    PhyloTree temp(*tree); // only use a copy of the tree !!!
	    subst += -temp.MAP(codontuple, weights, parsi_base, true); // Fitch Algorithm 
	}
    }
    return omegas[maxU];
}

// get vector of log likelihoods of omegas for one single codon tuple

vector<double> CodonEvo::loglikForCodonTuple(vector<string> &seqtuple, PhyloTree *ctree, PhyloTree *tree, int &subs){
    if (seqtuple.size() != ctree->numSpecies())
        throw ProjectError("CodonEvo::logLikForCodonTuple: inconsistent number of species.");
    for(int i=1; i<seqtuple.size();i++){
        if(seqtuple[0].length() != 3 || seqtuple[i].length() != 3){
            throw ProjectError("CodonEvo::loglikForCodonTuple: codon tuple has not length 3");
        }
    }
    int numCodons;
    vector<double> logliks(k, 0.0);
    Seq2Int s2i(3);
    vector<int> codontuple(ctree->numSpecies(), 64); // 64 = missing codon                                                     
    bit_vector bv(ctree->numSpecies(), 0);
    numCodons = 0;
    //cout << "codontuple(int): ";
    for(size_t s=0; s < ctree->numSpecies(); s++){
      if (seqtuple[s].size()>0)
	try {
	  codontuple[s] = s2i(seqtuple[s].c_str());
	  numCodons++;
	  bv[s] = 1;
	} catch(...){} // gap or n character                                                                              
      //cout << codontuple[s] << "|";
    }
    //cout << endl;
    if (numCodons >= 2){
      if(tree){
	//cout << "### original species tree = ";
	//tree->printTree();
	// calculate number of substitutions by fitch algorithm
	PhyloTree *tr;       
	unordered_map<bit_vector,PhyloTree*, boost::hash<bit_vector>>::iterator topit = GeneMSA::topologies.find(bv);
	if(topit == GeneMSA::topologies.end()) // insert new topology
	  tr = new PhyloTree(*tree);
	else
	  tr  = new PhyloTree(*topit->second);
	tr->prune(bv,NULL);
	GeneMSA::topologies.insert(pair<bit_vector,PhyloTree*>(bv,tr));
	subs = tr->fitch(codontuple);
      }
      // calculate log likelihood of omegas
      for (int u=0; u < k; u++){ // loop over omegas	
	logliks[u] = ctree->pruningAlgor(codontuple, this, u);
      }
    }
    return logliks;
}

double CodonEvo::graphOmegaOnCodonAli(vector<string> &seqtuple, PhyloTree *tree, int refSpeciesIdx){
  if (seqtuple.size() != tree->numSpecies())
    throw ProjectError("CodonEvo::graphOmegaOnCodonAli: inconsistent number of species. Tree and alignment file must have the same (number of) species.");
  
  vector<string> speciesnames;
  tree->getSpeciesNames(speciesnames);
  int n = 0;  // number of nucleotide triples

  for(int i=0; i<seqtuple.size();i++){
    if(seqtuple[i].length() != 0){
      n = seqtuple[i].length()/3;
      break;
    }
  }
  int numCodons;
  double Eomega, loglik, sum;
  Seq2Int s2i(3);
  vector<double> logliks(k, 0.0);
  vector<double> postprobs(k, 0.0);
  // output: stats of codon alignment sites 
  vector<char> aminoAcidsRef(n,'\0');
  vector<int> refPos(n,-1);
  vector<int> numSubst(n,-1);
  vector<double> postmeanSite(n,0.0);
  vector<double> stdPostmeanSite(n,0.0);
  vector<double> postProb_gt1(n,0.0);
  int codonIdx = 0;
  int refCodonIdx = 0;  

  for (int i=0; i<n; i++){
    vector<double> logliksSite(k,0.0);
    Eomega = sum = 0.0;      
    vector<int> codontuple(tree->numSpecies(), 64); // 64 = missing codon
    numCodons = 0;
    PhyloTree pruned_tr(*tree);
    for(size_t s=0; s < tree->numSpecies(); s++){
      if (seqtuple[s].size()>0){
	try {
	  codontuple[s] = s2i(seqtuple[s].c_str() + 3*i);
	  numCodons++;
	  if(s == refSpeciesIdx){
	    refPos[codonIdx] = refCodonIdx;
	    refCodonIdx++;
	  }
	} catch(...){ // gap or n character
	  pruned_tr.drop(speciesnames[s]);
	} 
      }else{ // species missing in alignfile
	pruned_tr.drop(speciesnames[s]);
      }
    }
    if(codontuple[refSpeciesIdx] == 64)
      aminoAcidsRef[codonIdx] = '-';
    else
      aminoAcidsRef[codonIdx] = GeneticCode::translate(codontuple[refSpeciesIdx]);
    
    if (numCodons >= 2){
      numSubst[codonIdx] = pruned_tr.fitch(codontuple);
      for (int u=0; u < k; u++){ // loop over omegas
	loglik = pruned_tr.pruningAlgor(codontuple, this, u);
  	logliks[u] += loglik;
	logliksSite[u] += loglik;
      }    
      double max_loglik = *max_element(logliksSite.begin(), logliksSite.end());
      for (int u=0; u < k; u++){
	sum += postprobs[u] = exp(logliksSite[u] - max_loglik) * omegaPrior[u];
      }
      for (int u=0; u < k; u++){
	postprobs[u] /= sum;
      }
      for (int u=0; u < k; u++){
	Eomega += postprobs[u] * omegas[u];
      }
      postmeanSite[codonIdx] = Eomega;
      for (int u=0; u < k; u++)
	stdPostmeanSite[codonIdx] += postprobs[u] * pow(omegas[u] - Eomega, 2);
      
      /*** 
      for(size_t s=0; s < tree->numSpecies(); s++){
	if (seqtuple[s].size()>0)
	  cout << codontuple[s] << "\t" << seqtuple[s].substr(i*3,3) << endl;
      }
      cout << "index\tomega\tloglik\tmaxloglik\tpostProb\tsum\tpostMean" << endl;
      for(int u=0; u < k; u++){
	cout << u << "\t" << omegas[u] << "\t" << logliksSite[u] << "\t" << max_loglik << "\t" << postprobs[u] << "\t" << sum << "\t" << Eomega << endl;
      }
      ***/

      for(int u=0; u<k; u++){
	if(omegas[u] > 1){
	  postProb_gt1[codonIdx] += postprobs[u];
	}
      }    
    }
    codonIdx++;
  }
  
  cout << setw(10) << "ali_pos" << setw(10) << "ref_pos" << setw(10) << "AS_ref" << setw(10) << "Pr(w>1)" << setw(10) << "post_mean" << setw(4) << "+-" << setw(10) << "SE_for_w" << setw(10) << "num_subst" << endl;
  for(int i = 0; i < codonIdx; i++){
    cout << setw(10)<< i << setw(10) << refPos[i] << setw(10) << aminoAcidsRef[i];
    if(postProb_gt1[i] > 0){
      cout << setw(13) << postProb_gt1[i];
      if(postProb_gt1[i] > 0.95)
	cout << "*";
      if(postProb_gt1[i] > 0.99)
	cout << "*";
      cout << setw(10) << postmeanSite[i] << setw(4) << "+-" << setw(13) << stdPostmeanSite[i] << setw(11) << numSubst[i] << endl;  
    }
    else{
      cout << setw(13) << -1 << setw(10) << -1 << setw(4) << "+-" << setw(13) << -1 << setw(11) << -1 << endl;
    }  
  }
  
    // posterior mean estimate of omega
    sum = 0.0;
    double max_loglik = *max_element(logliks.begin(), logliks.end());
    for (int u=0; u < k; u++){
      sum += postprobs[u] = exp(logliks[u] - max_loglik) * omegaPrior[u];
    }
    //        cout << "posterior distribution of omega" << endl;
    for (int u=0; u < k; u++){
      postprobs[u] /= sum;
      //	cout << omegas[u] << " " << postprobs[u] << endl;
    }
    Eomega = 0.0;
    for (int u=0; u < k; u++){
      Eomega += postprobs[u] * omegas[u];
    }
    /***
        cout << "whole alignment: " << endl << "index\tomega\tloglik\tmaxloglik\tpostProb\tsum\tpostMean" << endl;
    for(int u=0; u < k; u++){
      cout << u << "\t" << omegas[u] << "\t" << logliks[u] << "\t" << max_loglik << "\t" << postprobs[u] << "\t" << sum << "\t" << Eomega << endl;
    }
    ***/
 
    cout << endl << "# reference species: " << speciesnames[refSpeciesIdx] << endl << "# posterior mean estimate of omega for whole alignment : " << Eomega << endl;
    return Eomega;
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

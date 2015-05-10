/**********************************************************************
 * file:    motif.cc
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  Weight matrices and other motifs
 * authors: Mario Stanke, mario@gobics.de
 *
 * date    |   author      |  changes 
 * --------|---------------|------------------------------------------ 
 * 08.08.02| Mario Stanke  | creation of the file
 * 10.05.05| Mario Stanke  | making gc_range_min, gc_range_max variable
 **********************************************************************/

#include "motif.hh"

// project includes
#include "properties.hh"
#include "projectio.hh" // for comment
#include "types.hh"

// standard C/C++ includes
#include <iostream>
#include <fstream>
#include <iomanip>  // for setw, setprecision


WeighingType BaseCount::weithType = equalWeights;
Matrix<double> BaseCount::weighingMatrix = 0;

BaseCount::BaseCount(){
    a = c = g = t = 0;
    ra = rc = rg = rt = 0.25;
}
BaseCount::BaseCount(int a, int c, int g, int t) {
    this->a = a;
    this->c = c;
    this->g = g;
    this->t = t;
    normalize();
}
BaseCount::BaseCount(const char *sequence, int len) {
    a = c = g = t = 0;
    ra = rc = rg = rt = 0.25;
    if (len<0) 
	len = strlen(sequence);
    addSequence(sequence, len);
    normalize();
}


void BaseCount::init(){
    
    int wtNumber = 1;
    WeighingType wt = equalWeights;
    try {
	wtNumber = Properties::getIntProperty( "/BaseCount/weighingType" );
    } catch (ProjectError) {
	cerr << "couldn't read the '/BaseCount/weighingType' property" << endl;
    }
    switch (wtNumber) {
	case 1: wt = equalWeights; break;
	case 2: wt = gcContentClasses; break;
	case 3: wt = multiNormalKernel; break;
    }
    weithType = wt;
    if (wt == multiNormalKernel) { 
	string weightMatrixFileName;
	try {
	    weightMatrixFileName = Properties::getProperty( "/BaseCount/weightMatrixFile" );
	} catch (ProjectError) {
	    cerr << "couldn't read the '/BaseCount/weightMatrixFile' property" << endl;
	}
	setWeightMatrix(weightMatrixFileName.c_str());
    }
}
void BaseCount::addSequence(const char *sequence, int len){
    for (int i=0; (i < len) && (sequence[i] != '\0'); i++ )
	addCharacter(sequence[i]);
}

void BaseCount::addCharacter(char nucleotide, bool subtract){
  switch( nucleotide ){
  case 'a': case 'A':
    if (subtract) a--; else a++; break;
  case 'c': case 'C':
    if (subtract) c--; else c++; break;
  case 'g': case 'G':
    if (subtract) g--; else g++; break;
  case 't': case 'T':
    if (subtract) t--; else t++; break;
  }
}

// double BaseCount::relFreq(int i){
//     switch(i) {
// 	case 0: return ra;
// 	case 1: return rc;
// 	case 2: return rg;
// 	case 3: return rt;
// 	default: 
// 	    throw InvalidNukleotideError("");
//     }
// }

int BaseCount::weight( BaseCount bc1, BaseCount bc2) {
    return (int) doubleWeight(bc1, bc2);
}

double BaseCount::doubleWeight( BaseCount bc1, BaseCount bc2) {
    bc1.normalize();
    bc2.normalize();
    switch (weithType) {
	case equalWeights:
	    return 1;
	case gcContentClasses:
	    return gcContentClassWeight(bc1, bc2);
	case multiNormalKernel: 
	    if (weighingMatrix.getColSize() != 4 || weighingMatrix.getRowSize() != 4){
		throw ProjectError("multinormal weighing, but have no 4x4 weight matrix!");
	    }
	    return multiNormalKernelWeight(bc1, bc2);
	default:
	    return 1;
    }
}
    
int BaseCount::gcContentWeight(BaseCount bc1, BaseCount bc2) {
    double gcContent1, gcContent2, gcContentDiff;
    gcContent1 = bc1.rg + bc1.rc;
    gcContent2 = bc2.rg + bc2.rc;
    gcContentDiff = gcContent1 - gcContent2;
	
    return (int) (1 * phi(gcContentDiff, 0.05));
}
/*
 * Use only the genes in the same isochore (gc content class)
 */

int BaseCount::gcContentClassWeight(BaseCount bc1, BaseCount bc2) {
    double gcContent1, gcContent2;
    gcContent1 = bc1.rg + bc1.rc;
    gcContent2 = bc2.rg + bc2.rc;
	
    if (gcContentClass(gcContent1) == gcContentClass(gcContent2)) {
	return 1;
    } else
	return 0;
}

int BaseCount::gcContentClass (double gcContent) {
    if (gcContent < .43) 
	return 0;
    if (gcContent < .51)
	return 1;
    if (gcContent < .57) 
	return 2;
    return 3;
}

double BaseCount::multiNormalKernelWeight(BaseCount bc1, BaseCount bc2){
    double z[4]; // a vector of differences of the nucleotide frequencies
    z[0]= bc1.ra - bc2.ra;
    z[1]= bc1.rc - bc2.rc;
    z[2]= bc1.rg - bc2.rg;
    z[3]= bc1.rt - bc2.rt;

    /* 
     * compute z * weighingMatrix * z^t manually
     */
    double temp[4] = {0,0,0,0};     // = z * weighingMatrix
    for (int j=0; j<4; j++) {
	for(int i=0; i<4; i++) {
	    temp[j] += z[i] * weighingMatrix[i][j];
	}
    }
    double t = 0;
    for (int i=0; i<4; i++) {
	t += temp[i]*z[i]; 
    }
    return (1 + 9 * exp(-t));
}

void BaseCount::setWeightMatrix(string matrixFileName){
    ifstream istrm((Constant::fullSpeciesPath() + matrixFileName).c_str());
    if( istrm ){    
	weighingMatrix.assign(4, 4);
	istrm >> comment;
	for (int z=0; z<4; z++) {
	    for (int s=0; s<4; s++) {
		istrm >> weighingMatrix[z][s];
	    }
	}
	// cout <<"Gewichtungsmatrix: \n" << weighingMatrix << endl;
	istrm.close();
    } else {
	string errorMess("Could't open the file with the weight matrix: ");
	errorMess.append(matrixFileName);
	throw ProjectError(errorMess.c_str());
    }

}

// density of centered normal distribution with variance sigma^2
double BaseCount::phi(double x, double sigma) {
    return 0.39894228 / sigma * exp(-x*x/2/sigma/sigma);
}

void BaseCount::normalize(){
    double sum = a+c+g+t;
    if (sum > 0.0){
	ra = a / sum;
	rc = c / sum;
	rg = g / sum;
	rt = t / sum;
    }
}

void BaseCount::reverse() {
    BaseCount temp = *this;
    a = temp.t;
    c = temp.g;
    g = temp.c;
    t = temp.a;
    ra = temp.rt;
    rc = temp.rg;
    rg = temp.rc;
    rt = temp.ra;
}

ostream& operator<<( ostream& out, const BaseCount& bc ){
    out << "(" << setprecision(3) << bc.ra << ", " 
	<< setprecision(3) << bc.rc << ", " 
	<< setprecision(3) << bc.rg << ", " 
	<< setprecision(3) << bc.rt << ")" ;
    return out;
}

Motif& Motif::operator= (const Motif & other){
  if (this != &other) { // protect against invalid self-assignment
    n = other.n;
    k = other.k;
    numSeqs = other.numSeqs;
    neighbors = other.neighbors;
    pseudocount = other.pseudocount;
    s2i = other.s2i;
    if (windowProbs)
      delete [] windowProbs;
    if (windowCounts)
      delete [] windowCounts;
    windowProbs = new vector<Double>[n];
    windowCounts = new vector<int>[n];
    for (int i=0; i<n; i++) {
      windowProbs[i] =  other.windowProbs[i];
      windowCounts[i] = other.windowCounts[i];
    }
  }
  // by convention, return *this
  return *this;
}

Motif::Motif(int length, int memory, int pseudocount, int neighbors) :
    n(length),
    k(memory),
    numSeqs(0),
    s2i(memory+1)
{
    this->pseudocount = pseudocount;
    this->neighbors = neighbors;
    windowProbs = new vector<Double>[n];
    windowCounts = new vector<int>[n];
    for (int i=0; i<n; i++) {
	windowProbs[i].assign(POWER4TOTHE(k+1), 0);
	windowCounts[i].assign(POWER4TOTHE(k+1), pseudocount);
    }
}

Motif::~Motif(){
    if (windowProbs)
	delete [] windowProbs;
    if (windowCounts)
	delete [] windowCounts;    
}

/*
 * add one sequence to the training set of the motif
 * seq is the beginning of the motiv, but
 * seq[-k] ... seq[n-1] or seq[0] ... seq[n+k-1] (reverse case) must be accessible!
 */
void Motif::addSequence(const char* seq, int weight, bool reverse){
  int pn;
  for (int i=0; i<n; i++) {
    try {
      pn = (reverse)? s2i.rev(seq+i) : s2i(seq+i-k);
      windowCounts[i][pn] += weight;
      // also add the pattern to the neighbors
      for (int j=1; j <= neighbors; j++) {
	if (i-j >= 0)
	  windowCounts[i-j][pn] += weight;
	if (i+j < n)
	  windowCounts[i+j][pn] += weight;
      }
    } catch (InvalidNucleotideError e) {}
  } 
  numSeqs += 1; // do not weight here
}

/*
 * evaluate a sequence 
 * seq is the beginning of the motif, but
 * for restrictions see addSequence
 */
Double Motif::seqProb(const char* seq, bool reverse, bool complement){
    Double prob(1);
    int pn;
    for (int i=0; i<n; i++) {
	try {
	    if (!reverse) {
		if (!complement) // normal case
		    pn = s2i(seq+i-k);
		else
		    throw ProjectError("Motif::seqProb: unknown alternative");
		prob *= windowProbs[i][pn];		
	    } else {
		if (complement)  // reverse complement
		    pn = s2i.rc(seq+i);
		else   // non-complement reverse
		    pn = s2i.rev(seq+i);
		prob *= windowProbs[n-1-i][pn];
	    }
	} catch (InvalidNucleotideError e) {
	    prob *= (float) 0.25;
	}
    }
    return prob;
}

void Motif::makeProbs() { 
    int sum;
    for (int i=0; i<n; i++) {
	for(int j=0; j<windowProbs[0].size(); j += 4) {
	    sum = windowCounts[i][j] + windowCounts[i][j+1]
		+ windowCounts[i][j+2] + windowCounts[i][j+3];
	    if (sum > 0) {
		windowProbs[i][j]   = (double)windowCounts[i][j] / sum;
		windowProbs[i][j+1] = (double)windowCounts[i][j+1] / sum;
		windowProbs[i][j+2] = (double)windowCounts[i][j+2] / sum;
		windowProbs[i][j+3] = (double)windowCounts[i][j+3] / sum;
	    } else {
		windowProbs[i][j] = windowProbs[i][j+1] = windowProbs[i][j+2] = windowProbs[i][j+3] =  0.25;
	    }
	}
    }
}

void Motif::clearCounts(){
    for (int i=0; i<n; i++) {
	windowCounts[i].assign(POWER4TOTHE(k+1), pseudocount);
    }
    numSeqs = 0;
}

void Motif::printProbs(){
    if (windowProbs) {
	for (int i=0; i<n; i++) {  
	    cout << i << "   \t";
	}
	cout << endl;
	for (int j=0; j<windowProbs[0].size(); j++) {
	    for (int i=0; i<n; i++) {  
		cout << windowProbs[i][j] <<"\t";
	    } 
	    cout << endl;
	}
    }
}

void Motif::write(ofstream &out){
  out << "# width of motif, n=\n" << n << endl;
  out << "# order of markov model, k=\n" << k << endl;
  out << "# markov chain emission probabilities" << endl;
  for (int i=0; i< n; i++) {
    out << setw(2) << i << "  ";
    for (int j=0; j < windowProbs[0].size(); j++) {
      out << windowProbs[i][j];
      if (j+1<windowProbs[0].size())
	out << "\t";
    }
    out << endl;
  }
}

void Motif::read(ifstream &in){
  if (!in)
    throw ProjectError("Error when reading in motif.");
  in >> comment >> n;
  in >> comment >> k;
    
  if (windowProbs)
    delete [] windowProbs;
  if (windowCounts)
    delete [] windowCounts;
  s2i = Seq2Int(k+1);
  windowProbs = new vector<Double>[n];
  for (int i=0; i<n; i++) {
    windowProbs[i].assign(POWER4TOTHE(k+1), 0);
  }
  
  int dummy;
  for (int i=0; i<n; i++) {
    in >> comment >> dummy;
    for (int j=0; j < windowProbs[i].size(); j++) {
      in >> windowProbs[i][j];
    }
  }
}

/*
 *
 * xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 * | k bases | n bases (actual motif)    |
 * The first k bases come from a uniform distribution as they are not parametrized.
 * The n bases from the motif are samples from the non-homogeneous markov chain of the motif. P_i(X_i | X_{i-1}, ..., X_{i-k})
 */
char* Motif::getSampleDNA() {
    char bases[5] = "acgt"; // 0->a, 1->c, 2->g, 3->t
    int dnaLength = n + k;
    char *sequence;

    sequence = new char[dnaLength+1];
    // first k bases from uniform distribution
    for (int i=0; i<k; i++) {
        sequence[i] = bases[(int) (4.0 * rand() / (1.0 + RAND_MAX))];
    }
    // sample 1 position at a time
    for (int i=k; i<dnaLength; i++) {
        double r=rand() / (1.0 + RAND_MAX);
        Double cumProb = 0.0;
        for (int base = 0; base < 4 && cumProb <= r; base++) {// try new base at position i, keep that base when cumulative prob. exceeds r
            sequence[i] = bases[base];
            int pn = s2i(sequence+i-k);
            cumProb = cumProb + windowProbs[i-k][pn];
        }
    }
    sequence[dnaLength] = '\0';
    return sequence;
}

Double Motif::getProbThreshold (double q, int numSamples){
    char* seq;
    vector<Double> all_probs;
    srand(1);
    for (int i=0; i<numSamples; i++) {
        seq = getSampleDNA();
        //cout << i << "   " << seq << endl;
        Double prob = Motif::seqProb(seq + k);
        all_probs.push_back(prob);
        delete [] seq;
    }
    return quantile(all_probs, q);
}

void ContentDecomposition::setProperties(){
    n = Constant::decomp_num_at * Constant::decomp_num_gc * Constant::decomp_num_steps;
    zus = new BaseCount[n];
    makeDecomposition();
}

void ContentDecomposition::makeDecomposition(){
    // gc content between default values gc_range_min=0.32 and gc_range_max=0.73    
    double quot = 1.25;
    int a = Constant::decomp_num_at, b = Constant::decomp_num_gc;
    
    for (int i=0; i<Constant::decomp_num_steps; i++) {
	double gc = Constant::gc_range_min + (Constant::gc_range_max - Constant::gc_range_min) * (i+1)/(Constant::decomp_num_steps+1);
	double at = 1-gc;
	double quot_at, quot_cg;
	
	for (int e=0; e<a; e++) {
	    for (int f=0; f<b; f++) {
		quot_at = (2-quot) + (2*(quot-1)) * (e+1)/(a+1);
		quot_cg = (2-quot) + (2*(quot-1)) * (f+1)/(b+1);	
		
		zus[a*b*i + e*b+f].ra = at/(1+quot_at);
		zus[a*b*i + e*b+f].rt = at/(1+1/quot_at);
		zus[a*b*i + e*b+f].rc = gc/(1+quot_cg);
		zus[a*b*i + e*b+f].rg = gc/(1+1/quot_cg);

	    }
	}
    }
}    	

BaseCount ContentDecomposition::getBaseCount(int i){
    return zus[i];
}

int ContentDecomposition::getNearestBaseCountIndex(BaseCount bc) {
    double maxWeight = -1;
    double weight;
    int ret = -1;
    for (int i=0; i<n; i++) {
	weight = BaseCount::doubleWeight(bc, zus[i]);
	if (weight > maxWeight) {
	    maxWeight = weight;
	    ret = i;
	}
    }
    return ret;
}

/*
 * ContentStairs constructor
 */
ContentStairs::ContentStairs(){
  idx = NULL;
  dna = NULL;
  try {
    GCwinsize = Properties::getIntProperty( "GCwinsize" );
  } catch (ProjectError) {
    GCwinsize = 10000; // default. 5000 is also good 
    // GCwinsize = -1 means the same as infinity, take whole sequence.
  }
}

/*
 * ContentStairs destructor
 */
ContentStairs::~ContentStairs(){
    if (idx)
	delete [] idx;
}

/*
 * ContentStairs constructor
 *
 * dna xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 *                     (        |         )
 *                         GCwinsize
 *     (|             )
 *        GCwinsize
 *   
 * The local GC index is computated by taking the average of the window and then taking the best match.
 * The averaging window is at least min(dnalen, GCwinsize), i.e. in particular the first and last win/2 bases
 * are assigned the same index.
 */

void ContentStairs::computeStairs(const char* dna){
  ContentDecomposition cd;
  n = strlen(dna);
  this->dna = dna; // keep only the pointer
  nextStep.clear();
  // to smooth out when GC content teeter-totters between two classes
  int totterywin = 1000; // a little more accurate than 200 and also than 2000
  if (idx)
    delete [] idx;
  idx = new int[n];
  for (int i=0; i < n; i++)
    idx[i] = -1;
  int x;
  int win = GCwinsize;
  if (win > n || win < 1)
    win = n;
  
  // first window goes from 0 to win-1
  BaseCount bc(dna, win);
  x = cd.getNearestBaseCountIndex(bc);
  for (int i=0; i <= int(win/2); i++)
    idx[i] = x;
  
  // mittle part for which the window is a proper substring of dna
  // window goes from i - int(win/2) to i + int((win+1)/2) - 1
  for (int i = int(win/2) + 1; i <= n - int((win+1)/2); i++){ // shift window by 1 to the right
    //    cout << "i=" << i << " adding " << dna[i + int((win+1)/2) - 1] << " removing " << dna[i - int(win/2) - 1] << endl;
    bc.addCharacter(dna[i + int((win+1)/2) - 1]);       // add new base to the right
    bc.addCharacter(dna[i - int(win/2) - 1], true); // subtract old base to the left
    //    bc.normalize();
    idx[i] = x = cd.getNearestBaseCountIndex(bc);
    //    cout << i << " " << bc << endl;
  }
  // last window goes from n-win to n-1
  for (int i= n - int((win+1)/2) + 1; i < n; i++)
    idx[i] = x;

  // do another smoothing step to prevent too frequent back-and forth switching,
  // when a region is on the boundary between two classes
  //                     ------
  // ---  ---  -- -- ----
  //    --   --  -  -
  //    ^ do not switch here
  //
  x = -2;
  int lastStep = 0; // position where the class changed last time
  for (int i=0; i < n; i++)
    if (idx[i] != x) {
      if (i - lastStep < totterywin && lastStep > 0 && idx[lastStep-1] == idx[i]){
	// have a short stretch and the same (different) class on both sides
	for (int j = lastStep; j<i; j++)
	  idx[j] = idx[i];
      }
      lastStep = i;
      x = idx[i];
    }
  // store nextSteps
  x = 0;
  for (int i=1; i < n; i++)
      if (idx[i] != idx[x]) {
	  nextStep[x] = i;
	  x = i;
      }
  nextStep[x] = n;
  if (nextStep[0] == 0)
      nextStep[0] = n;
  if (nextStep[1] == 0)
      nextStep[1] = nextStep[0]; // this is an exception so the Viterbi algorithm can start at 1 instead of 0

  //  for (map<int,int>::iterator it = nextStep.begin(); it != nextStep.end(); ++it)
  //    cout << "# " << it->first << " => " << it->second << " : " << idx[it->first] << '\n';
}

int ContentStairs::getNextStep(int from) {
    map<int, int>::iterator it = nextStep.find(from);
    if (it != nextStep.end())
	return it->second;
    else 
	return n;
}

/* 	Function definitions
 
	Created: 20-September-2011	
	Last modified: 8-February-2012
*/

#ifndef _FILTERBAM_HH
#define _FILTERBAM_HH

#include <api/BamReader.h>  
#include <api/BamAlignment.h>
#include <iostream> 
#include <vector>
#include <string>
#include <stdio.h>
#include <stdexcept>
#include <sstream>
#include <stdlib.h>
#include <time.h>
#include <iomanip>
#include <getopt.h> 
#include<list> 

using namespace std;  
using namespace BamTools;    
 
string getReferenceName(const RefVector &refData, int RefID);

// For computing Coverage and Insertions
uint32_t sumMandIOperations(const vector<CigarOp> &cigar);
uint32_t sumDandIOperations(const vector<CigarOp> &cigar);
void printElapsedTime(int tEnd, int tStart); 

// For option initialization
struct globalOptions_t {
	bool best;
	bool help;
	bool noIntrons;
	bool paired;
	bool uniq;
	bool verbose;  
  	bool pairwiseAlignments;	
	int insertLimit;
	int maxIntronLen;
	int maxSortesTest; 
	int minCover;
	int minId;
	float uniqThresh;
	const char* commonGeneFile;
  	const char* inputFile;
  	const char* outputFile;
	const char* pairBedFile;
};

globalOptions_t initOptions(int argc, char *argv[]);


// Class definition 
class MatePairs
{
   public:
  	  int alIt;
  	  int alJit;
  	  float score;
      MatePairs();
      MatePairs(const MatePairs &);
      MatePairs(int it, int jit, float score);
      ~MatePairs(){};
 	  void setValues(int it, int jit, float score);
      MatePairs &operator=(const MatePairs &rhs);
      int operator==(const MatePairs &rhs) const;
      int operator<(const MatePairs &rhs) const;
};

void printMatePairs(vector<MatePairs> matepairs, vector<BamAlignment> &qali);


// Class definition 
class PairednessCoverage
{
   public:
	  int coord; // pStart or pEnd coordinates
  	  int label; // {-1} for pEnd and {1} for pStart
  	  string chr; // chromosome name
      PairednessCoverage();
      PairednessCoverage(const PairednessCoverage &);
  PairednessCoverage(int coord, int label, string chr);
      ~PairednessCoverage(){};
   	  void setValues(int coord, int label, string chr);
      PairednessCoverage &operator=(const PairednessCoverage &rhs);
      int operator==(const PairednessCoverage &rhs) const;
      int operator<(const PairednessCoverage &rhs) const;
};
#endif

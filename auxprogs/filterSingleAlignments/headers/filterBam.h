/* 	Function definitions
 
	Created: 20-September-2011	
	Last modified: 11-April-2012
*/

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
 
// For debugging
string printCigarOperations(vector<CigarOp> cigar);
string getReferenceName(vector<RefData> refData, int RefID);
string getReferenceLength(vector<RefData> refData, int RefID);

// For computing Coverage and Insertions
uint32_t sumMandIOperations(vector<CigarOp> cigar, string printFlag);
uint32_t sumDandIOperations(vector<CigarOp> cigar, string printFlag);
int printElapsedTime(int tEnd, int tStart); 

// For option initialisation
struct globalOptions_t {
	bool best;
	bool help;
	bool noIntrons;
	bool paired;
	bool uniq;
	bool verbose;  	
	int insertLimit;
	int maxIntronLen;
	int maxSortesTest; 
	int minCover;
	int minId;
	int minIntronLen;
	float uniqThresh;
	const char* commonGeneFile;
  	const char* inputFile;
  	const char* outputFile;
	const char* pairBedFile;
};

globalOptions_t initOptions(int argc, char *argv[]);


#ifndef MATEPAIRS_H
#define MATEPAIRS_H
// Class definition 
class MatePairs
{
   friend ostream &operator<<(ostream &, const MatePairs &);

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
  	  bool operator() (const MatePairs &lhs, const MatePairs &rhs) const;
};
#endif

void printMatePairs(vector<MatePairs> someList, vector<BamAlignment> &qali);
void printList(list<MatePairs> someList);

#ifndef SINGLEALIGNMENT_H
#define SINGLEALIGNMENT_H
// Class definition 
class SingleAlignment
{
   friend ostream &operator<<(ostream &, const SingleAlignment &);

   public:
  		// Attributes
  	  	BamAlignment *al;
	  	const RefVector *refData;
		float coverage;
		float percId;
		float score;

		// Methods
		SingleAlignment();
		SingleAlignment(const SingleAlignment &);
		SingleAlignment(BamAlignment *al, float coverage, float percId, const RefVector *refData);
		~SingleAlignment(){};
		void setValues(BamAlignment *al, float coverage, float percId, const RefVector *refData);

		// Operators
		SingleAlignment &operator=(const SingleAlignment &rhs);
		bool operator()(const SingleAlignment &lhs, const SingleAlignment &rhs) const;
};
#endif


// This function is required for built-in STL list functions like sort
// We just compare equality of mates in terms of their coverage
bool byScore(SingleAlignment *i, SingleAlignment *j);

// Sends to "cout" the contents of a vector of singleAlignments
void printQali(vector<SingleAlignment*> &vecSingleAlignments);
void printQaliOld(vector<SingleAlignment*> vecSingleAlignments);

/* 	Function definitions
 
	Created: 20-September-2011	
	Last modified: 3-January-2012
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
#ifndef QALISIZE
#define QALISIZE 8
#endif
#ifndef SIZE
#define SIZE 3
#endif

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
	const char* pairbedFile;
};

globalOptions_t initOptions(int argc, char *argv[]);



///////////////////////////////////////////////////////
// QALI vector class definition
///////////////////////////////////////////////////////

#ifndef QALIv_H
#define QALIv_H
// Class definition 
class Qaliv
{
 
	public:
   	  string name[QALISIZE];
	  Qaliv();
  	  Qaliv(string *message);
  	  friend bool compare_group(const Qaliv& a, const Qaliv& b);
  	  friend bool compare_value(const Qaliv& a, const Qaliv& b);
	  friend ostream& operator<<(ostream& os, const Qaliv& x);

};
#endif

void printVector(vector<Qaliv> someVector);
void printFormattedVector(vector<Qaliv> someVector);
bool compare_criteria_1(const Qaliv& a, const Qaliv& b); 
bool compare_criteria_2(const Qaliv& a, const Qaliv& b);

///////////////////////////////////////////////////////


#ifndef MATEPAIRS_H
#define MATEPAIRS_H
// Class definition 
class MatePairs
{
   friend ostream &operator<<(ostream &, const MatePairs &);

   public:
      float values[SIZE];
      MatePairs();
      MatePairs(const MatePairs &);
      ~MatePairs(){};
 	  void push(float values[SIZE]);
      MatePairs &operator=(const MatePairs &rhs);
      int operator==(const MatePairs &rhs) const;
      int operator<(const MatePairs &rhs) const;
};
#endif

void printList(list<MatePairs> someList);
void printFormattedList(list<MatePairs> someList);




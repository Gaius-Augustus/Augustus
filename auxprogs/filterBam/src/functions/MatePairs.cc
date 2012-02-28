/*	Implementation of the MATEPAIRS class

	Created: 7-November-2011	
	Last modified: 10-February-2012
*/

#include <iostream>
#include <list>
#include <string>
#include <vector>
#include <api/BamAlignment.h>
#include <cassert>
#include <algorithm>
#include <functional>
#include <map>

using namespace std;
using namespace BamTools;



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

// Constructor
MatePairs::MatePairs()   // Constructor
{
  alIt = 0;
  alJit = 0;
  score = 0;
}

// Constructor through values
MatePairs::MatePairs(const MatePairs &copyin)   // Copy constructor to handle pass by value.
{
	this->alIt = copyin.alIt; 	               
	this->alJit = copyin.alJit; 	               
	this->score = copyin.score; 	               
}

// Initialisation method
MatePairs::MatePairs(int it, int jit, float score)
{
	this->alIt = it;
	this->alJit = jit;
	this->score = score;
}

// Initialisation method
void MatePairs::setValues(int it, int jit, float score)
{
	this->alIt = it;
	this->alJit = jit;
	this->score = score;
}

// Ostream operator
ostream &operator<<(ostream &output, const MatePairs &aaa)
{
  output << aaa.alIt << ' ' << aaa.alJit << aaa.score;
  // output << endl;
  return output;
}

// Assignment operator
MatePairs& MatePairs::operator=(const MatePairs &rhs)
{  
  this->alIt = rhs.alIt;
  this->alJit = rhs.alJit;
  this->score = rhs.score;
  return *this;
}

// Equality comparison operator
int MatePairs::operator==(const MatePairs &rhs) const
{
	if( this->score != rhs.score) return 0;
	return 1;
}

// This function is required for built-in STL list functions like sort
// We just compare equality of mates in terms of their score
int MatePairs::operator<(const MatePairs &rhs) const
{
   if( this->score >= rhs.score ) return 1;
   return 0;
}


void printMatePairs(vector<MatePairs> matepairs, vector<BamAlignment> &qali)
{
  int it, jit;
  float score;

  for (int iter=0; iter < matepairs.size(); iter++)
	{
	  it = matepairs.at(iter).alIt;
	  jit = matepairs.at(iter).alJit;
	  score = matepairs.at(iter).score;
	  cout << "(" << it << "," << jit << ") = (" 
		   << qali.at(it).Name << "," << qali.at(jit).Name << "),"			 	
		   << " scoreMate=" << score << endl;
	}
}



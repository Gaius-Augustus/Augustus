/*	Implementation of the MATEPAIRS class

	Created: 7-November-2011	
	Last modified: 10-February-2012
*/

#include "filterBam.h"
#include "bamaccess.hh"
#include <iostream>
#include <list>
#include <string>
#include <vector>
#include <cassert>
#include <algorithm>
#include <functional>
#include <map>

using namespace std;

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

// Initialization method
MatePairs::MatePairs(int it, int jit, float score)
{
	this->alIt = it;
	this->alJit = jit;
	this->score = score;
}

// Initialization method
void MatePairs::setValues(int it, int jit, float score)
{
	this->alIt = it;
	this->alJit = jit;
	this->score = score;
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

void printMatePairs(vector<MatePairs> &matepairs, vector<BamAlignmentRecord_> &qali)
{
  int it, jit;
  float score;

  for (unsigned int iter=0; iter < matepairs.size(); iter++)
	{
	  it = matepairs.at(iter).alIt;
	  jit = matepairs.at(iter).alJit;
	  score = matepairs.at(iter).score;
	  cout << "(" << it << "," << jit << ") = (" 
		   << qali.at(it)->getQueryName()  << "," << qali.at(jit)->getQueryName()  << "),"			 	
		   << " scoreMate=" << score << endl;
	}
}



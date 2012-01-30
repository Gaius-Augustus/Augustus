/*	Implementation of a SINGLEREAD class to store information for this type of reads

	Created: 30-January-2012	
	Last modified: 30-January-2012
*/

#include <iostream>
#include <list>
#include <string>
#include <vector>
#include <api/BamAlignment.h>
#include <cassert>
#include <algorithm>
#include <functional>

using namespace std;
using namespace BamTools;


// Class definition 
class SingleRead
{
   friend ostream &operator<<(ostream &, const SingleRead &);

   public:
  	  BamAlignment al;
  	  float coverage;
  	  float percId;
  	  float score;

      SingleRead();
      SingleRead(const SingleRead &);
  	  SingleRead(BamAlignment al, float coverage, float percId);
      ~SingleRead(){};
 	  void setValues(BamAlignment al, float coverage, float percId);
      SingleRead &operator=(const SingleRead &rhs);
      int operator==(const SingleRead &rhs) const;
      int operator<(const SingleRead &rhs) const;
  	  bool operator() (const SingleRead &lhs, const SingleRead &rhs) const;
};

// Constructor
SingleRead::SingleRead()   // Constructor
{
  	  al();
  	  coverage = 0;
  	  percId = 0;
  	  score = coverage+percId;
}

// Constructor through values
SingleRead::SingleRead(const SingleRead &copyin)   // Copy constructor to handle pass by value.
{
	this->al = copyin.al; 	               
	this->coverage = copyin.coverage; 	               
	this->percId = copyin.percId; 	               
	this->score = copyin.score; 	               
}

// Initialisation method
SingleRead::SingleRead(BamAlignment al, float coverage, float percId)
{
	this->al = al; 	               
	this->coverage = coverage; 	               
	this->percId = percId; 	               
	this->score = score; 	               
}


// Initialisation method
void SingleRead::setValues(BamAlignment al, float coverage, float percId)
{
	this->al = al; 	               
	this->coverage = coverage; 	               
	this->percId = percId; 	               
	this->score = score; 	               
}

// Ostream operator
ostream &operator<<(ostream &output, const SingleRead &aaa)
{
  output << aaa.al.Name << ' ' << aaa.coverage << ' ' << aaa.percId << ' ' << aaa.score;
  // output << endl;
  return output;
}

// Assignment operator
SingleRead& SingleRead::operator=(const SingleRead &rhs)
{  
  this->al = rhs.al; 	               
  this->coverage = rhs.coverage; 	               
  this->percId = rhs.percId; 	               
  this->score = rhs.score; 	               
  return *this;
}

// Equality comparison operator
int SingleRead::operator==(const SingleRead &rhs) const
{
	if( this->score != rhs.score) return 0;
	return 1;
}

// This function is required for built-in STL list functions like sort
// We just compare equality of mates in terms of their score
int SingleRead::operator<(const SingleRead &rhs) const
{
  if( this->score < rhs.score ) return 1;
   return 0;
}


void toString(SingleRead read)
{
  BamAlignment al = read.al;
  float coverage, percId, score;
  cout << al.Name << "," << al.RefId << "," << al.coverage << "," << al.percId << "," << al.score << endl; 
}

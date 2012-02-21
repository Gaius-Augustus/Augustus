/*	Implementation of a SINGLEREAD class to store information for this type of reads

	Created: 30-January-2012	
	Last modified: 11-February-2012
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
  	  BamAlignment *al;
  	  float coverage;
  	  float percId;
  	  float score;

      SingleRead();
      SingleRead(const SingleRead &);
  	  SingleRead(BamAlignment *al, float coverage, float percId);
      ~SingleRead(){};
 	  void setValues(BamAlignment *al, float coverage, float percId);
  	  void toString();
  	  void size();
      // SingleRead &operator=(const SingleRead &rhs);
      // int operator==(const SingleRead &rhs) const;
      // int operator<(const SingleRead &rhs) const;
  	  // bool operator() (const SingleRead &lhs, const SingleRead &rhs) const;
};

// Constructor
SingleRead::SingleRead()   // Constructor
{
  	  al;
  	  coverage;
  	  percId;
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
SingleRead::SingleRead(BamAlignment *al, float coverage, float percId)
{
	this->al = al; 	               
	this->coverage = coverage; 	               
	this->percId = percId; 	               
	this->score = coverage+percId; 	               
}


// Initialisation method
void SingleRead::setValues(BamAlignment *al, float coverage, float percId)
{
	this->al = al; 	               
	this->coverage = coverage; 	               
	this->percId = percId; 	               
	this->score = coverage+percId; 	               
}

// Printing method
void SingleRead::toString()
{
  BamAlignment *al = this->al;
  float coverage = this->coverage; 
  float percId = this->percId;
  float score = this->score;

  if (this->al->IsMapped())
	{
  	  cout << al->Name << ' ' << al->RefID << ' ' << coverage << ' ' << percId << ' ' << score << endl;
	} else {
		cout << "Unmapped alignment" << endl;
  	}
}

// Size of elements
void SingleRead::size()
{
  BamAlignment *al = this->al;
  float coverage = this->coverage; 
  float percId = this->percId;
  float score = this->score;

  cout << "size(SingleRead)=" << sizeof(al) + sizeof(coverage) + sizeof(percId) + sizeof(score) << "bytes." << endl;

}

// Ostream operator
ostream &operator<<(ostream &output, const SingleRead &sr)
{
  if (sr.al->IsMapped())
	{
	  output << sr.al->Name << ' ' << sr.al->RefID << ' ' << sr.coverage << ' ' << sr.percId << ' ' << sr.score;
	} else {
	  output << "Unmapped alignment";
    }
  return output;
}

// // Assignment operator
// SingleRead& SingleRead::operator=(const SingleRead &rhs)
// {  
//   this->al = rhs.al; 	               
//   this->coverage = rhs.coverage; 	               
//   this->percId = rhs.percId; 	               
//   this->score = rhs.score; 	               
//   return *this;
// }

// // Equality comparison operator
// int SingleRead::operator==(const SingleRead &rhs) const
// {
// 	if( this->score != rhs.score) return 0;
// 	return 1;
// }

// // This function is required for built-in STL list functions like sort
// // We just compare equality of mates in terms of their score
// int SingleRead::operator<(const SingleRead &rhs) const
// {
//   if( this->score < rhs.score ) return 1;
//    return 0;
// }



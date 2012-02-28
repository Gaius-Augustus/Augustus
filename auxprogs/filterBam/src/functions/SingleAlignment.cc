/*	Implementation of a SINGLEALIGNMENT class to store information for this type of alignments

	Created: 30-January-2012	
	Last modified: 22-February-2012
*/

#include <iostream>
#include <list>
#include <string>
#include <vector>
#include <api/BamAlignment.h>
#include <api/BamReader.h>
#include <cassert>
#include <algorithm>
#include <functional>

using namespace std;
using namespace BamTools;

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
      	int operator==(const SingleAlignment &rhs) const;
      	int operator<(const SingleAlignment &rhs) const;
  	  	bool operator() (const SingleAlignment &lhs, const SingleAlignment &rhs) const;
};

// Constructor
SingleAlignment::SingleAlignment()   // Constructor
{
  	  al;
	  refData;
  	  coverage = 0;
  	  percId = 0;
	  score = coverage+percId;
}

// Constructor through another SingleAlignment
SingleAlignment::SingleAlignment(const SingleAlignment &copyin)   // Copy constructor to handle pass by value.
{
	this->al = copyin.al; 	               
	this->coverage = copyin.coverage; 	               
	this->percId = copyin.percId; 	               
	this->score = copyin.score; 
	this->refData = copyin.refData;	               
}

// Constructor through value assignment
SingleAlignment::SingleAlignment(BamAlignment *al, float coverage, float percId, const RefVector *refData)
{
	this->al = al; 	               
	this->coverage = coverage; 	               
	this->percId = percId; 	               
	this->score = coverage+percId; 	               
	this->refData = refData;
}


// Initialisation method
void SingleAlignment::setValues(BamAlignment *al, float coverage, float percId, const RefVector *refData)
{
	this->al = al; 	               
	this->coverage = coverage; 	               
	this->percId = percId; 	               
	this->score = coverage+percId; 	               
	this->refData = refData; 
}


// Ostream operator
ostream &operator<<(ostream &output, const SingleAlignment &sr)
{
  int size = sizeof(sr.al) + sizeof(sr.coverage) + sizeof(sr.percId) + sizeof(sr.score) + sizeof(sr.refData);

  if (sr.al->IsMapped())
	{
	  output << sr.al->Name << ' ' << (*sr.refData).at(sr.al->RefID).RefName << ' ' << sr.coverage << ' ' << sr.percId << ' ' << sr.score << ' ' << size << " bytes";
	} else {
	  output << "Unmapped alignment";
    }
  return output;
}

// Assignment operator
SingleAlignment& SingleAlignment::operator=(const SingleAlignment &rhs)
{  
  this->al = rhs.al; 	               
  this->coverage = rhs.coverage; 	               
  this->percId = rhs.percId; 	               
  this->score = rhs.score; 	               
  return *this;
}

// Equality comparison operator
int SingleAlignment::operator==(const SingleAlignment &rhs) const
{
	if( this->score != rhs.score) return 0;
	return 1;
}

// This function is required for built-in STL list functions like sort
// We just compare equality of mates in terms of their score
int SingleAlignment::operator<(const SingleAlignment &rhs) const
{
  if( this->score < rhs.score ) return 1;
   return 0;
}



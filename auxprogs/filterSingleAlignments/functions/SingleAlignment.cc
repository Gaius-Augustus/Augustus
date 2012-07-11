/*	Implementation of a SINGLEALIGNMENT class to store information for this type of alignments

	Created: 30-January-2012	
	Last modified: 11-April-2012
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
		bool operator()(const SingleAlignment &lhs, const SingleAlignment &rhs) const;

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
	this->refData = copyin.refData;	               
	this->coverage = copyin.coverage; 	               
	this->percId = copyin.percId; 	               
	this->score = copyin.score; 
}

// Constructor through value assignment
SingleAlignment::SingleAlignment(BamAlignment *al, float coverage, float percId, const RefVector *refData)
{
    this->al = al; 	               
	this->refData = refData;
	this->coverage = coverage; 	               
	this->percId = percId; 	               
	this->score = coverage+percId; 	               
}


// Initialisation method
void SingleAlignment::setValues(BamAlignment *al, float coverage, float percId, const RefVector *refData)
{
    this->al = al; 	               
	this->refData = refData; 
	this->coverage = coverage; 	               
	this->percId = percId; 	               
	this->score = coverage+percId; 	               
}


// Ostream operator
ostream &operator<<(ostream &output, const SingleAlignment &sr)
{
  // sizeof(sr.al) = 8 bytes // sizeof(sr.refData) = 8 bytes // sizeof(sr.coverage) = 4 bytes
  // sizeof(sr.percId) = 4 bytes // sizeof(sr.score) = 4 bytes
  int size = sizeof(sr.al) + sizeof(sr.refData) + sizeof(sr.coverage) + sizeof(sr.percId) + sizeof(sr.score);

  if (sr.al->IsMapped())
	{
	  output << (*sr.al).Name << " " << (*sr.refData).at(sr.al->RefID).RefName << " " << sr.coverage << " " << sr.percId << ", score= " << sr.score << " " << size << " bytes" << ", RefID=" << sr.al->RefID; //RefLength=" << (*sr.refData).at(sr.al->RefID).RefLength << " Size pointed by *singleAlignment: " << size << " bytes";
	} else {
	  output << "Unmapped alignment";
    }
  return output;
}

// Assignment operator
SingleAlignment& SingleAlignment::operator=(const SingleAlignment &rhs)
{  
  // this->al = rhs.al; 	               
  // this->refData = rhs.refData; 
  this->coverage = rhs.coverage; 	               
  this->percId = rhs.percId; 	               
  this->score = rhs.score; 	               
  return *this;
}


// Assignment operator
bool SingleAlignment::operator() (const SingleAlignment &lhs, const SingleAlignment &rhs) const
{  
  return lhs.score < rhs.score; 	               
}


// We just compare equality of mates in terms of their score
bool byScore(SingleAlignment *i, SingleAlignment *j) 
{ 
	return ((*i).score > (*j).score); 
}


// Sends to "cout" the contents of a vector of singleAlignments
void printQaliOld(vector<SingleAlignment*> vecSingleAlignments)
{
  // This is the size of the vector container
  cout << "Within printQaliOld, vecSingleAlignments is located at " << &vecSingleAlignments << endl;
  cout << "The size of the input argument is " << sizeof(vecSingleAlignments) << " bytes" << endl;
  for (int it=0; it < vecSingleAlignments.size(); it++)
	{
	  cout << (*vecSingleAlignments.at(it)) << endl;
	}
}

// Sends to "cout" the contents of a vector of singleAlignments
void printQali(vector<SingleAlignment*> &vecSingleAlignments)
{
  // This is the size of the vector container
  cout << "Within printQali, vecSingleAlignments is located at " << &vecSingleAlignments << endl;
  cout << "The size of the input argument is " << sizeof(vecSingleAlignments) << " bytes" << endl;
  for (int it=0; it < vecSingleAlignments.size(); it++)
	{
	  cout << (*vecSingleAlignments.at(it)) << endl;
	}
}

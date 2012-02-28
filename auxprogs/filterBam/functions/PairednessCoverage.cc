/*	Implementation of a bespoke class that contains pairedness coverage information

	Created: 29-January-2012	
	Last modified: 29-January-2012
*/

#include <iostream>
#include <list>
#include <string>
#include <vector>
#include <cassert>
#include <algorithm>
#include <functional>

using namespace std;


// Class definition 
class PairednessCoverage
{
   friend ostream &operator<<(ostream &, const PairednessCoverage &);

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
  	  bool operator() (const PairednessCoverage &lhs, const PairednessCoverage &rhs) const;
};

// Constructor
PairednessCoverage::PairednessCoverage()   // Constructor
{
  coord = 0;
  label = 0;
  chr = "";
}

// Constructor through values
PairednessCoverage::PairednessCoverage(const PairednessCoverage &copyin)   // Copy constructor to handle pass by value.
{
	this->coord = copyin.coord; 	               
	this->label = copyin.label; 	               
	this->chr = copyin.chr; 	               
}

// Initialisation method
PairednessCoverage::PairednessCoverage(int coord, int label, string chr)
{
	this->coord = coord; 	               
	this->label = label; 	               
	this->chr = chr; 	               
}


// Initialisation method
void PairednessCoverage::setValues(int coord, int label, string chr)
{
	this->coord = coord;
	this->label = label;
	this->chr = chr;
}

// Ostream operator
ostream &operator<<(ostream &output, const PairednessCoverage &aaa)
{
  output << aaa.chr << ' '<< aaa.coord << ' ' << aaa.label;
  // output << endl;
  return output;
}

// Assignment operator
PairednessCoverage& PairednessCoverage::operator=(const PairednessCoverage &rhs)
{  
	this->coord = rhs.coord; 	               
	this->label = rhs.label; 	               
	this->chr = rhs.chr;
  	return *this;
}

// Equality comparison operator
int PairednessCoverage::operator==(const PairednessCoverage &rhs) const
{
	if( this->coord != rhs.coord) return 0;
	return 1;
}

// This function is required for built-in STL list functions like sort
// We just compare equality of mates in terms of their score
int PairednessCoverage::operator<(const PairednessCoverage &rhs) const
{
   if( this->coord < rhs.coord ) return 1;
   return 0;
}




/*	Implementation of a bespoke class that contains pairedness coverage information

	Created: 29-January-2012	
	Last modified: 29-January-2012
*/

#include "filterBam.h"
#include <iostream>
#include <list>
#include <string>
#include <vector>
#include <cassert>
#include <algorithm>
#include <functional>

using namespace std;

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

// Initialization method
PairednessCoverage::PairednessCoverage(int coord, int label, string chr)
{
	this->coord = coord; 	               
	this->label = label; 	               
	this->chr = chr; 	               
}


// Initialization method
void PairednessCoverage::setValues(int coord, int label, string chr)
{
	this->coord = coord;
	this->label = label;
	this->chr = chr;
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




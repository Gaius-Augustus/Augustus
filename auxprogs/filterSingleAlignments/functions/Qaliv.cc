/*	Implementation of the QALI class

	Created: 22-October-2011	
	Last modified: 2-January-2012
*/

#include <api/BamReader.h>
#include <api/BamAlignment.h>
#include <iostream>
#include <vector> 
#include <string>
#include <cctype>
#include <stdio.h>
#include <algorithm>
#include <iterator>
#include <cstdlib>
#define QALISIZE 8

using namespace std;


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

// Constructor
Qaliv::Qaliv()   
{
  for (int it=0; it<QALISIZE; it++)
	{
	  name[it] = "";
	}
}

// Constructor through string array
Qaliv::Qaliv(string *message)
{
  for (int it=0; it<QALISIZE; it++)
	{
	  name[it] = *(message+it);
	}
 
}

// Sorting by $tName
bool compare_criteria_1(const Qaliv& a, const Qaliv& b) 
{ 
	return a.name[1] < b.name[1]; 
}

// Sorting by $tStart
bool compare_criteria_2(const Qaliv& a, const Qaliv& b) 
{ 
	return a.name[4] < b.name[4]; 
}

// Ostream operator
ostream &operator<<(ostream &output, const Qaliv &x)
{
  for (int it=0; it<QALISIZE; it++)
	{
	  output << x.name[it] << ' ';
	}
  output << endl;	
   return output;
}


void printVector(vector<Qaliv> someVector)
{
  vector<Qaliv>::iterator it = someVector.begin();
  cout << " ";
  for (it; it != someVector.end(); it++)
	{
	  cout << *it << " ";
	}
}

void printFormattedVector(vector<Qaliv> someVector)
{
  vector<Qaliv>::iterator it = someVector.begin();
  int jit;
  for (it; it != someVector.end(); it++)
	{
	  cout << "\t[ ";
	  for (jit=0; jit<QALISIZE; jit++)
		{
		  cout << (*it).name[jit] << " ";
		}
	  cout << " ]," << endl;
	}
}

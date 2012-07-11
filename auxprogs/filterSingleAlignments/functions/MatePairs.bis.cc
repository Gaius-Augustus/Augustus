/*	Implementation of the MATEPAIRS class

	Created: 7-November-2011	
	Last modified: 16-November-2011
*/

#include <iostream>
#include <list>
#include <string>
#define SIZE 3	

using namespace std;


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
      ~MatePairs(){};
 	  void push(float values[SIZE]);
      MatePairs &operator=(const MatePairs &rhs);
      int operator==(const MatePairs &rhs) const;
      int operator<(const MatePairs &rhs) const;
};

// Constructor
MatePairs::MatePairs()   // Constructor
{
  for (int it=0; it<SIZE; it++)
	{
	  values[it] = 0;
	}
}

// Constructor through values
MatePairs::MatePairs(const MatePairs &copyin)   // Copy constructor to handle pass by value.
{                             
	for (int it=0; it<SIZE; it++)
		{
	  	this->values[it] = copyin.values[it]; 	               
		}  
}

// push method
void MatePairs::push(float values[SIZE])
{
  for (int it=0; it<SIZE; it++)
	{
  	this->values[it] = values[it];
	}
}

// Ostream operator
ostream &operator<<(ostream &output, const MatePairs &aaa)
{
	for (int it=0; it<SIZE; it++)
		{
		  output << aaa.values[it] << ' ';
		}
  	// output << endl;	
   	return output;
}

// Assignment operator
MatePairs& MatePairs::operator=(const MatePairs &rhs)
{  
	for (int it=0; it<SIZE; it++)
	  {
		this->values[it] = rhs.values[it];
	  }
   	return *this;
}

// Equality comparison operator
int MatePairs::operator==(const MatePairs &rhs) const
{
  int it = 0;
  while (it < SIZE)
	{
   	if( this->values[it] != rhs.values[it]) 
	  {return 0;}
	it++;
	}	  
}

// This function is required for built-in STL list functions like sort
int MatePairs::operator<(const MatePairs &rhs) const
{
   //if( this->values[1] == rhs.values[1] && this->values[4] <  rhs.values[4] ) return 1;
   if( this->values[2] <  rhs.values[2] ) return 1;
   return 0;
}


void printList(list<MatePairs> someList)
{
  list<MatePairs>::iterator it = someList.begin();
  cout << " ";
  for (it; it != someList.end(); it++)
	{
	  cout << *it << " " << endl;
	}
}

void printFormattedList(list<MatePairs> someList)
{
  list<MatePairs>::iterator it = someList.begin();
  int jit;
  for (it; it != someList.end(); it++)
	{
	  cout << "\t[ ";
	  for (jit=0; jit<SIZE; jit++)
		{
		  cout << (*it).values[jit] << " ";
		}
	  cout << " ]," << endl;
	}
}

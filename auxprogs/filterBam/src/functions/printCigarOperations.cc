/*	Prints and returns a string with the operations 
	and lengths of an alignment's CIGAR string 

	Created: Friday, September 16, 2011
	Last modified:  7-October-2011
*/
 
#include <api/BamReader.h> 
#include <api/BamAlignment.h>
#include <iostream>
#include <vector>
#include <string>   
#include <stdio.h>   
#include <boost/lexical_cast.hpp>
 
using namespace std; 
using namespace BamTools; 
 
string printCigarOperations(vector<CigarOp> cigar)
{
	int cigarSize = cigar.size();
 	uint32_t cigarLength;
  	char cigarType;
	string cigarOperations = "";;

	// Scanning through all CIGAR operations
	for (int it=0; it<cigarSize; it++)
		{
		// Length is number of bases
		cigarLength = cigar.at(it).Length;	
		// Type means operations, i.e. MINDSHP 
		cigarType = cigar.at(it).Type;  
		// Building output string
		cigarOperations = boost::lexical_cast<string>(cigarLength);
		cigarOperations.append(";");
		cigarOperations.append(1, cigarType);
		printf("%s", cigarOperations.c_str());
		} // end for

	return cigarOperations;
}

/*	Computes Insertions wrt Query and Reference through the summation 
	of operations D and I in the CIGAR string

	Created: 26-September-2011
	Last modified: 11-November-2011
*/

#include <api/BamReader.h>
#include <api/BamAlignment.h>
#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>

using namespace std;
using namespace BamTools;

// In SAM, sumDandIOperations = $qBaseInsert+$tBaseInsert
uint32_t sumDandIOperations(vector<CigarOp> cigar, string printFlag)
{
	int cigarSize = cigar.size();
	uint32_t sumDandI = 0;
 	uint32_t cigarLength;
  	char cigarType;

	// Scanning through all CIGAR operations
	for (int it=0; it<cigarSize; it++)
		{
		// Length is number of bases
		cigarLength = cigar.at(it).Length;	
		// Type means operations, i.e. MINDSHP 
		cigarType = cigar.at(it).Type;  

		if (cigarType == 'D' || cigarType == 'I')
			{
		  	// Sum and print operations
			sumDandI = cigarLength + sumDandI;	
			if (!printFlag.compare("print"))
			  {cout << cigarLength << ";" << cigarType << ",";}
			}	
		} // end for

	return sumDandI;
}

/*	Sums and (if desired), prints the total length of the M and I
	cigar operations

	Created: Friday, September 16, 2011
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

// In SAM, sumMandIOperations = $qEnd-$qStart in PSL
uint32_t sumMandIOperations(vector<CigarOp> cigar, string printFlag)
{
	int cigarSize = cigar.size();
	uint32_t sumMandI = 0;
 	uint32_t cigarLength;
  	char cigarType;

	// Scanning through all CIGAR operations
	for (int it=0; it<cigarSize; it++)
		{
		// Length is number of bases
		cigarLength = cigar.at(it).Length;	
		// Type means operations, i.e. MINDSHP 
		cigarType = cigar.at(it).Type;  

		if (cigarType == 'M' || cigarType == 'I')
			{
		  	// Sum and print operations
			sumMandI = cigarLength + sumMandI;	
			if (!printFlag.compare("print"))
			  {cout << cigarLength << ";" << cigarType << ",";}
			}	
		} // end for

	return sumMandI;
}

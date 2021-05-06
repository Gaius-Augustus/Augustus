/*	Computes Insertions wrt Query and Reference through the summation 
	of operations D and I in the CIGAR string

	Created: 26-September-2011
	Last modified: 11-November-2011
*/

#include <SeqLib/BamReader.h>
#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>

using namespace std;
using namespace SeqLib;

// In SAM, sumDandIOperations = $qBaseInsert+$tBaseInsert
uint32_t sumDandIOperations(const Cigar &cigar)
{
	uint32_t sumDandI = 0;

	// Scanning through all CIGAR operations
	for (auto it : cigar)
		{
		if (it.Type() == 'D' || it.Type() == 'I')
			sumDandI = it.Length() + sumDandI;	
		} // end for

	return sumDandI;
}

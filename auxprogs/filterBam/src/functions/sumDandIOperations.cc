/*	Computes Insertions wrt Query and Reference through the summation 
	of operations D and I in the CIGAR string

	Created: 26-September-2011
	Last modified: 11-November-2011
*/

#include <api/BamAlignment.h>
#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>

using namespace std;
using namespace BamTools;

// In SAM, sumDandIOperations = $qBaseInsert+$tBaseInsert
uint32_t sumDandIOperations(const vector<CigarOp> &cigar)
{
    uint32_t sumDandI = 0;

    // Scanning through all CIGAR operations
    for (auto it : cigar)
    {   // Length is number of bases
        // Type means operations, i.e. MINDSHP 
        if (it.Type == 'D' || it.Type == 'I')
        {
            sumDandI += it.Length;
        }
    } // end for

    return sumDandI;
}

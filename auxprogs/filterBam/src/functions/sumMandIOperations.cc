/** Sums the total length of the M and I cigar operations.
*/

#include <api/BamAlignment.h>
#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>

using namespace std;
using namespace BamTools;

// In SAM, sumMandIOperations = $qEnd-$qStart in PSL
uint32_t sumMandIOperations(const vector<CigarOp> &cigar)
{
    uint32_t sumMandI = 0;

    // Scanning through all CIGAR operations
    for (auto it : cigar)
    {   // Length is number of bases
        // Type means operations, i.e. MINDSHP 
        if (it.Type == 'M' || it.Type == 'I')
        {
            sumMandI += it.Length;
        }
    } // end for
    
    return sumMandI;
}

/*	Returns string with length of the reference of an alignment sequence

	Created: 22-September-2011
	Last modified: 6-October-2011
*/

#include <api/BamReader.h>
#include <api/BamAlignment.h>
#include <vector>
#include <string>
#include <stdio.h>
#include <stdexcept>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace BamTools;

string getReferenceLength(vector<RefData> refData, int RefID)
{
	string rLength;
	try {
	  	rLength = boost::lexical_cast<string>(refData.at(RefID).RefLength);
		} catch (out_of_range& oor) {
	  	// BAM spec: refID=-1 for a read without a mapping position
	  	rLength.append("[printReferenceLength: ");
		rLength.append(oor.what());
		rLength.append("]"); 
	}
	// printf("%s", rLength.c_str());

	return rLength;
}

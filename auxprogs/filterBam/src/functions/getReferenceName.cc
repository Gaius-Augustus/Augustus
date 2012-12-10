/*	Returns string with name of the reference of an alignment sequence

	Created: Friday, September 16, 2011
	Last modified: 7-October-2011
*/

#include <api/BamReader.h>
#include <api/BamAlignment.h>
#include <vector>
#include <string>
#include <stdio.h>
#include <stdexcept>
#include <sstream>

using namespace std;
using namespace BamTools;

// string getReferenceName(vector<RefData> refData, int RefID)
string getReferenceName(const RefVector &refData, int RefID)
{
	string refName = "";
	try {
	  	refName = refData.at(RefID).RefName;
		} catch (out_of_range& oor) {
	  	// BAM spec: refID=-1 for a read without a mapping position
	  	refName.append("[printReferenceName: ");
		refName.append(oor.what());
		refName.append("]"); 
	}
	// printf("%s", rName.c_str());

	return refName;
}

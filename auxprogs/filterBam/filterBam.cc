/*	Filters the coverage of a BAM alignment file

	Created: 28-September-2011
	Last modified: 20-December-2011
*/ 
 
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <api/BamAlignment.h> 
#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>
#include <stdexcept>
#include <time.h>
#include <getopt.h> 
#include <sstream>
#include <stdlib.h>
#include <iomanip> 
#include <boost/lexical_cast.hpp>
#include <unordered_map>

// Default options
int help = 0;
int maxintronlen = 500000;
int minintronlen = 35;
int minId = 92;
int minCover = 80;
int insertLimit = 10;
bool nointrons = false;

using namespace BamTools;
using namespace std;

string bool_cast(const bool b) {
    ostringstream ss;
    ss << boolalpha << b;
    return ss.str();
}
uint32_t sumMandIOperations(vector<CigarOp> cigar, string printFlag);
uint32_t sumDandIOperations(vector<CigarOp> cigar, string printFlag);
string getReferenceName(vector<RefData> refData, int RefID);
int printElapsedTime(int tEnd, int tStart);


// Definition of global variables
static const char *optString = "abcdefghijklmnv:h?";
extern int opterr; // Display error if opterr=0

struct globalOptions_t {
	int help;
	int maxintronlen;
	int minintronlen;
	int maxSortesTest; 
	int minId;
	int minCover;
	int insertLimit;
	float uniqthresh;
	bool paired;
	bool uniq;
	bool nointrons;
	bool verbose;  	
};

// longOpts array
static const struct option longOpts[] = {
    { "help", no_argument, NULL, 'h' },
    { "maxintronlen", required_argument, NULL, 'a' },
	{ "minintronlen", required_argument, NULL, 'b' },
  	{ "minId", required_argument, NULL, 'd' },
    { "minCover", required_argument, NULL, 'e' },
	{ "insertLimit", required_argument, NULL, 'f' },
    { "nointrons", no_argument, NULL, 'm' },
    { NULL, no_argument, NULL, 0 }
};


// Display usage when --help
void displayUsage(int argc, char *argv[])
{
 	cout <<  "-----------------------------------------------------------------------------------" << endl;
	cout <<  " USAGE: " << argv[0] << " in.bam out.bam OPTIONS [option args] \n";
  	cout <<  "-----------------------------------------------------------------------------------" << endl;
	cout <<  " Available options are                                                             " << endl;
  	cout <<  "-----------------------------------------------------------------------------------" << endl;
	cout <<  "  NOTE: 			   								                                 " << endl;
	cout <<  "  The present version of this filter works with sorted and unsorted BAM files\n"     << endl;
  	cout <<  "-----------------------------------------------------------------------------------" << endl;
	cout <<  "  --help             display this menu" << endl;
	cout <<  "  --maxintronlen=n   maximal separation of paired reads (default " << maxintronlen << ")" << endl;
	cout <<  "  --minintronlen=n   minimal     ''     ''   ''    ''   (default " << minintronlen << ")" << endl;
	cout <<  "  --minId=n          minimal percentage of identity (default " << minId << ")" << endl;
	cout <<  "  --minCover=n       minimal percentage of coverage of the query read (default " << 
			 	minCover << ")" << endl;
	cout <<  "  --insertLimit=n    maximum assumed size of inserts (default " << insertLimit << ")" << endl;
	cout <<  "  --nointrons        do not allow longer gaps -for RNA-RNA alignments- (default " << 
				nointrons << endl; 
	cout << endl;
}


// Options initialisation
globalOptions_t initOptions(int argc, char *argv[])
{
  	int opt = 0, longIndex;
	struct globalOptions_t globalOptions;

	// Display usage if only ./program is provided
	if (argc==1)
	  {
		displayUsage(argc, argv);    
		exit(EXIT_FAILURE);
	  }

    // Initialize globalOptions with default values
	globalOptions.help = help;
	globalOptions.maxintronlen = maxintronlen;
	globalOptions.minintronlen = minintronlen;
	globalOptions.minId = minId;
	globalOptions.minCover = minCover;
	globalOptions.insertLimit = insertLimit;
	globalOptions.nointrons = nointrons;

	// Capturing options
	opt = getopt_long_only(argc, argv, optString, longOpts, &longIndex);
	// Scanning through options
    while( opt != -1 ) {
        switch( opt ) {                
            case 'a'	:	globalOptions.maxintronlen = atoi(optarg);	break;
            case 'b'	:	globalOptions.minintronlen = atoi(optarg);	break;
  			case 'd'	:	globalOptions.minId = atoi(optarg);			break;
			case 'e'	:	globalOptions.minCover = atoi(optarg);		break;
            case 'f'	:	globalOptions.insertLimit = atoi(optarg);	break;
			case 'm'	:	globalOptions.nointrons = 1;			    break;
			case 'h'	:   // fall-through is intentional
            case '?'	:	displayUsage(argc, argv);					exit(EXIT_FAILURE);	// break;          
			default		:	/* N.B. You won't actually get here*/		break;
        } // end switch
        
        opt = getopt_long(argc, argv, optString, longOpts, &longIndex);
    } // end while    

	return globalOptions;
}

 
int main(int argc, char *argv[])
{
  // Variable definition
  const char *fileName = argv[1];
  const char *outputFileName = argv[2];
  BamReader reader;
  BamWriter writer;
  BamAlignment al;
  // Time stamps
  time_t tStart, tEnd, tElapsed; 
  // Alignment attributes and variables
  string qName;
  uint32_t qLength;
  uint32_t RefID;
  uint32_t editDistance;
  vector<CigarOp> cigar;
  // Reference data
  vector<RefData> refData;
  string rName;
  // Cigar data
  char cigarType;
  uint32_t cigarLength;
  uint32_t sumMandI;
  int cigarSize; 
  // Filter variables
  float coverage;
  float percid;
  int32_t sumDandI;  // sumDandI = qBaseInsert+tBaseInsert, in some cases
  // Counters
  int outIntrons = 0;
  int outMinId = 0;
  int outMinCover = 0;
  int line = 0;
  // Initialising options
  struct globalOptions_t globalOptions;
  globalOptions = initOptions(argc, argv);
  int minCover = globalOptions.minCover;
  int minId = globalOptions.minId;
  int insertLimit = globalOptions.insertLimit;
  bool nointrons = globalOptions.nointrons;

  // Starting timer
  tStart = time(NULL);    

  // Opening BAM file
  if (!reader.Open(fileName)) 
	{
	cerr << "Could not open input BAM files. " << endl;
	return 0;
  	}

  // Retrieve 'metadata' from BAM files, these are required by BamWriter
  refData = reader.GetReferenceData();
  const SamHeader header = reader.GetHeader();

  // Open BamWriter, HEADER and REFERENCE will be written into the file
  if (!writer.Open(outputFileName, header, refData)) 
	{
	cerr << "Could not open output BAM file" << endl;
	return 0;
	}

  // Scanning through alignments
  nextAlignment:
  	while (reader.GetNextAlignment(al)) 
		{

		  // Displays line progress
		  line++;
		  if (line%100000==1)
			{
			  printf("\rProcessed line: %d", line);
			}

		  // Fetching alignment information
		  cigar = al.CigarData;
		  cigarSize = cigar.size();
		  qName = al.Name; 			  // query name
		  qLength = al.Length; 		  // query length
		  RefID = al.RefID;			  // ID of reference seq. (later used)
		  sumMandI = 0;				  // Equiv to $qEnd-$qStart in PSL
		  sumDandI = 0;
		  al.GetTag("NM", editDistance); // edit distance (see SAM spec)


		  // Working only on reads whose Reference seq ID has a defined mapping
		  rName = getReferenceName(refData, RefID);
		  if (rName.find("printReferenceName")!=-1)
		  	{	
		  	  goto nextAlignment;
		  	}

		  // Computing the pctge of Identity
		  percid = 100*(qLength-editDistance)/qLength;
		  RefID = al.RefID;
		  // Skip if not complying with percid
		  if (percid < minId)
		  	{
		  	  outMinId++;
		  	  goto nextAlignment;
		  	}

		  // Computing the coverage
		  sumMandI = sumMandIOperations(cigar, "no");
		  coverage = (float)100*sumMandI/qLength;  	
		  // Skip if not complying with minCover
		  if (coverage < minCover)
		  	{	
		  	  outMinCover++;
		  	  goto nextAlignment;
		  	}
	
		  // Number of inserts
		  sumDandI = sumDandIOperations(cigar, "no");
		  // Save if complying with insertLimit
		  if (nointrons && sumDandI > insertLimit)
		  	{
		  	  outIntrons++;
		  	  goto nextAlignment;
		  	}

		  // Saving cases that do comply
		  writer.SaveAlignment(al);
	
		} // end while	
  
	// Displaying results to STDOUT
	cout <<  "\n        filtered:" << endl;
	cout <<  "----------------:" << endl;
	cout <<  "percent identity: " << outMinId << endl;
	cout <<  "coverage        : " << outMinCover << endl;
	if (nointrons)
	  {
		cout <<  "nointrons       : " << outIntrons << endl;
	  }

	// Displaying command line
	for (int it=0; it<argc;it++)
	  {
		cout << argv[it] << " ";
	  }
	cout << endl;

	// Closing file handles
	reader.Close();
	writer.Close();

	// Ending timer and printing elapsed time
	tEnd = time(NULL);    
	tElapsed = printElapsedTime(tEnd, tStart);     

} // end main



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

string getReferenceName(vector<RefData> refData, int RefID)
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

int printElapsedTime(int tEnd, int tStart)
{
  int tElapsed;
  tElapsed = difftime(tEnd, tStart);     
  
  cout << "Elapsed time: " << setw(3) << tElapsed << " seconds.\n"; 

  return tElapsed;
}

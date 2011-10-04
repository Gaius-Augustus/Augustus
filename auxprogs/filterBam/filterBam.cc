/*	Filters the coverage of a BAM alignment file

	Created: 28-September-2011
	Last modified: 4-October-2011
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

// Default options
#define MAXINTRONLEN  500000
#define MININTRONLEN  35
#define MAXSORTESTEST  100000 
#define MINID  92
#define MINCOVER  80
#define UNIQTHRESH  0.96 
#define UNIQ  0
#define NOINTRONS  0
#define BEST  0
#define COMMONGENEFILE
#define PAIRBEDFILE
#define PAIRED  0
#define VERBOSE  0
#define INSERTLIMIT 10
#define HELP  0

using namespace BamTools;
using namespace std;

// Definition of global variables
static const char *optString = "xictunbcdpvl:h?";
extern int opterr; // Display error if opterr=0

struct globalOptions_t {
	int maxintronlen;
	int minintronlen;
	int maxSortesTest; 
	int minId;
	int minCover;
	float uniqthresh;
	int uniq;
	bool nointrons;
	int best;
	string commongenefile;
	string pairbedfile;
	int paired;
	int verbose;
	int help;
  	int insertLimit;
	string qnamestems; 
	string cmdline;
};

globalOptions_t initOptions(int argc, char *argv[]);

// longOpts array
static const struct option longOpts[] = {
    { "help", no_argument, NULL, 'h' },
    { "maxintronlen", required_argument, NULL, 'x' },
    { "minId", required_argument, NULL, 'i' },
    { "minCover", required_argument, NULL, 'c' },
    { "uniqthresh", required_argument, NULL, 't' },
    { "uniq", required_argument, NULL, 'u' },
    { "nointrons", no_argument, NULL, 'n' },
    { "best", required_argument, NULL, 'b' },
    { "commongenefile", required_argument, NULL, 'g' },
    { "pairbedfile", required_argument, NULL, 'd' },
    { "paired", required_argument, NULL, 'p' },
    { "verbose", required_argument, NULL, 'v' },
	{ "insertLimit", required_argument, NULL, 'l'},
    { NULL, no_argument, NULL, 0 }
};


// Display usage when --help
void displayUsage(int argc, char *argv[])
{
	cout << "USAGE: available options are" << endl;
  	cout << "---------------------" << endl;
	cout <<  "Usage: " << argv[0] << " in.bam out.bam\n";
	cout <<  "  PREREQUISITE: Does BAM need to be sorted??\n";
	cout <<  "                	   If so, sort with samtools\n";
	cout <<  "  if option 'paired' is used then it expects .f,.r or /1,/2 suffixes of mate pairs\n";
	cout <<  "  \n";
	cout <<  "  options:\n";
	cout <<  "  --pairbed=s        file name of pairedness coverage:\n";
	cout <<  "                     a .bed format file in which for each position the number of filtered\n";
	cout <<  "                     read pairs is reported that contain the position in or between the reads\n";
	cout <<  "  --minId=n          minimal percentage of identity (default 92)\n";
	cout <<  "  --minCover=n       minimal percentage of coverage of the query read (default 80)\n";
	cout <<  "  --uniq             take only best match and only, when second best is much worse (default false)\n";
	cout <<  "  --uniqthresh       threshold % for uniq, second best must be at most this fraction of best (default .96)\n";
	cout <<  "  --best             output all best matches that satisfy minId and minCover\n";
	cout <<  "  --commongenefile=s file name in which to write cases where one read maps to several different genes\n";
	cout <<  "  --nointrons        do not allow longer gaps (for RNA-RNA alignments)\n";
	cout <<  "  --paired           require that paired reads are on opposite strands of same target(default false)\n";
	cout <<  "  --maxintronlen=n   maximal separation of paired reads (default 500000\n";
	cout <<  "  --verbose          output debugging info (default false)\n";
	cout << "help            : --h or --help" << endl;
    exit(EXIT_FAILURE);
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
	  }

    // Initialize globalOptions 
	globalOptions.maxintronlen = MAXINTRONLEN;
	globalOptions.minintronlen = MININTRONLEN;
	globalOptions.maxSortesTest = MAXSORTESTEST; 
	globalOptions.minId = MINID;
	globalOptions.minCover = MINCOVER;
	globalOptions.uniqthresh = UNIQTHRESH; 
	globalOptions.uniq = UNIQ;
	globalOptions.nointrons = NOINTRONS;
	globalOptions.best = BEST;
	globalOptions.commongenefile = "";
	globalOptions.pairbedfile = "";
	globalOptions.paired = PAIRED;
	globalOptions.verbose = VERBOSE;
	globalOptions.help = HELP;
	globalOptions.insertLimit = INSERTLIMIT;
	globalOptions.qnamestems = ""; 
	globalOptions.cmdline = " ";

	// Capturing options
	opt = getopt_long_only(argc, argv, optString, longOpts, &longIndex);
	// Scanning through options
    while( opt != -1 ) {
        switch( opt ) {                
            case 'x':
			  	globalOptions.maxintronlen = atoi(optarg);
                break;
			case 'i':
			  	globalOptions.minId = atoi(optarg);
				break;
			case 'c':
			  	globalOptions.minCover = atoi(optarg);
			    break;
			case 't':
			  	globalOptions.uniqthresh = atof(optarg);
			    break;
			case 'u':
			  	globalOptions.uniq = atoi(optarg);
			    break;
			case 'n':
			  	globalOptions.nointrons = 1;
			    break;
			case 'b':
			  	globalOptions.uniqthresh = atoi(optarg);
			    break;
			case 'g':
			  	globalOptions.commongenefile = optarg;
			    break;
			case 'd':
			  	globalOptions.pairbedfile = optarg;
			    break;
			case 'p':
			  	globalOptions.paired = atoi(optarg);
			    break;
			case 'v':
			  	globalOptions.verbose = atoi(optarg);
			    break;
			case 'l':
			  	globalOptions.insertLimit = atoi(optarg);
			    break;

			case 'h':   // fall-through is intentional
            case '?':
			  	displayUsage(argc, argv);
                break;
               
            default:
			  	// N.B. You won't actually get here
                break;
        } // end switch
        
        opt = getopt_long(argc, argv, optString, longOpts, &longIndex);
    } // end while    

	return globalOptions;
}


// Definition of other auxiliary functions
int printElapsedTime(int tEnd, int tStart)
{
  int tElapsed;
  tElapsed = difftime(tEnd, tStart);     
  
  cout << "Elapsed time: " << setw(3) << tElapsed << " seconds.\n"; 

  return tElapsed;
}


// In SAM, sumMandIOperations = $qEnd-$qStart in PSL
int sumMandIOperations(vector<CigarOp> cigar, string printFlag)
{
	int cigarSize = cigar.size();
	int sumMandI = 0;
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
int sumDandIOperations(vector<CigarOp> cigar, string printFlag)
{
	int cigarSize = cigar.size();
	int sumDandI = 0;
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




int main(int argc, char *argv[])
{
  // Variable definition
  const char *inputFileName = argv[argc-2];
  const char *outputFileName = argv[argc-1];
  BamReader reader;
  BamWriter writer;
  BamAlignment al;
  time_t tStart, tEnd, tElapsed; 
  vector<CigarOp> cigar;
  char cigarType;
  uint32_t cigarLength, qLength, RefID, editDistance;
  string qName, rName;
  int cigarSize, sumMandI;
  float coverage, percid;
  int32_t InsertSize;  // InsertSize = qBaseInsert+tBaseInsert, in some cases
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
  if (!reader.Open(inputFileName)) 
	{
	cerr << "Could not open input BAM files. " << endl;
	return 0;
  	}

  // Retrieve 'metadata' from BAM files, these are required by BamWriter
  const RefVector refData = reader.GetReferenceData();
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
		qName = al.Name.substr(0,14); // query name
		qLength = al.Length; 		  // query length
		RefID = al.RefID;			  // ID of reference seq. (later used)
		sumMandI = 0;				  // Equiv to $qEnd-$qStart in PSL
		InsertSize = 0;
		al.GetTag("NM", editDistance); // edit distance (see SAM spec)

		// Pctge of Identity
		percid = 100*(qLength-editDistance)/qLength;
		RefID = al.RefID;
		// Save if coplying with percid
		if (percid < minId)
	  	{
			outMinId++;
			goto nextAlignment;
	  	}

		// Computing the coverage
		sumMandI = sumMandIOperations(cigar, "no");
		coverage = (float)100*sumMandI/qLength;  	
		// Save if complying with minCover
		if (coverage < minCover)
		  {	
			outMinCover++;
			goto nextAlignment;
		  }	

		// Number of inserts
		InsertSize = sumDandIOperations(cigar, "no");
		// Save if complying with insertLimit
		if (nointrons && InsertSize > insertLimit)
		  {
			outIntrons++;
			goto nextAlignment;
		  }
	
		// Save in case all criteria was satisfied
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


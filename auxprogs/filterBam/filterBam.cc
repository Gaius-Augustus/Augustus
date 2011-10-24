/*	Filters the coverage of a BAM alignment file

	Created: 28-September-2011
	Last modified: 24-October-2011
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
#define QALISIZE 8

// Default options
#define HELP  0
#define MAXINTRONLEN  500000
#define MININTRONLEN  35
#define MAXSORTESTEST  100000 
#define MINID  92
#define MINCOVER  80
#define INSERTLIMIT 10
#define UNIQTHRESH  0.96 
#define COMMONGENEFILE
#define PAIRBEDFILE
#define PAIRED  0
#define UNIQ  0
#define NOINTRONS  0
#define BEST  0
#define VERBOSE  0

using namespace BamTools;
using namespace std;

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
	string commongenefile;
	string pairbedfile;
	bool paired;
	bool uniq;
	bool nointrons;
	bool best;
	bool verbose;  	
};

globalOptions_t initOptions(int argc, char *argv[]);

// longOpts array
static const struct option longOpts[] = {
    { "help", no_argument, NULL, 'h' },
    { "maxintronlen", required_argument, NULL, 'a' },
	{ "minintronlen", required_argument, NULL, 'b' },
	{ "maxSortesTest", required_argument, NULL, 'c' },
  	{ "minId", required_argument, NULL, 'd' },
    { "minCover", required_argument, NULL, 'e' },
	{ "insertLimit", required_argument, NULL, 'f' },
    { "uniqthresh", required_argument, NULL, 'g' }, 	// 'h' option is for "help" 
    { "commongenefile", required_argument, NULL, 'i' }, // 'h' option is for "help"
    { "pairbedfile", required_argument, NULL, 'j' },
    { "paired", no_argument, NULL, 'k' },
    { "uniq", no_argument, NULL, 'l' },
    { "nointrons", no_argument, NULL, 'm' },
    { "best", no_argument, NULL, 'n' },
    { "verbose", no_argument, NULL, 'v' },
    { NULL, no_argument, NULL, 0 }
};


// Display usage when --help
void displayUsage(int argc, char *argv[])
{
	cout <<  " Usage: " << argv[0] << " in.bam out.bam\n";
  	cout <<  "--------------------------------------------------" << endl;
	cout <<  " Available options are                            " << endl;
  	cout <<  "--------------------------------------------------" << endl;
	cout <<  "  PREREQUISITE: Does BAM need to be sorted??\n";
	cout <<  "                	   If so, sort with samtools\n";
  	cout <<  "--------------------------------------------------" << endl;
	cout <<  "  --help             display this menu" << endl;
	cout <<  "  --maxintronlen=n   maximal separation of paired reads (default " << MAXINTRONLEN << ")" << endl;
	cout <<  "  --minintronlen=n   minimal     ''     ''   ''    ''   (default " << MININTRONLEN << ")" << endl;
	cout <<  "  --maxSortesTest=n  maximal sortedness (default " << MAXSORTESTEST << ")" << endl;
	cout <<  "  --minId=n          minimal percentage of identity (default " << MINID << ")" << endl;
	cout <<  "  --minCover=n       minimal percentage of coverage of the query read (default " << 
			 	MINCOVER << ")" << endl;
	cout <<  "  --insertlimit=n    maximum assumed size of inserts (default " << INSERTLIMIT << ")" << endl;
	cout <<  "  --uniqthresh=n     threshold % for uniq, second best must be at most this" << endl;
	cout <<  "                     fraction of best (default " << UNIQTHRESH << ") " << endl;
	cout <<  "  --commongenefile=s file name in which to write cases where one read maps to \n" <<
	         "                     several different genes" << endl;
	cout <<  "  --pairbedfile=s    file name of pairedness coverage:" << endl;
	cout <<  "                     options:"  << endl;
	cout <<  "	\t\t a .bed format file in which for each position the number of" << endl;
	cout <<  "	\t\t filtered read pairs is reported that contain the position in" << endl; 
	cout <<  "  \t\t\t or between the reads" << endl;
	cout <<  "  --paired           require that paired reads are on opposite strands of same target" << endl;
	cout <<  "  \t\t     (default " << PAIRED << ")" << endl;
	cout <<  "                     NOTE: " << endl;
	cout <<  "	\t\t if option 'paired' is used then it expects .f,.r or /1,/2" << endl; 
	cout <<  "	\t\t suffixes of mate pairs" << endl;
	cout <<  "  --uniq             take only best match and only, when second best is much worse (default " << 
			 	UNIQ << ")" << endl;
	cout <<  "  --nointrons        do not allow longer gaps -for RNA-RNA alignments- (default " << 
				NOINTRONS << endl; 
	cout <<  "  --best             output all best matches that satisfy minId and minCover (default " << 
				BEST << ")" << endl;
	cout <<  "  --verbose          output debugging info (default " << VERBOSE << ")" << endl;
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
	globalOptions.help = HELP;
	globalOptions.maxintronlen = MAXINTRONLEN;
	globalOptions.minintronlen = MININTRONLEN;
	globalOptions.maxSortesTest = MAXSORTESTEST; 
	globalOptions.minId = MINID;
	globalOptions.minCover = MINCOVER;
	globalOptions.insertLimit = INSERTLIMIT;
	globalOptions.uniqthresh = UNIQTHRESH; 
	globalOptions.commongenefile = "";
	globalOptions.pairbedfile = "";
	globalOptions.paired = PAIRED;
	globalOptions.uniq = UNIQ;
	globalOptions.nointrons = NOINTRONS;
	globalOptions.best = BEST;
	globalOptions.verbose = VERBOSE;

	// Capturing options
	opt = getopt_long_only(argc, argv, optString, longOpts, &longIndex);
	// Scanning through options
    while( opt != -1 ) {
        switch( opt ) {                
            case 'a':
			  	globalOptions.maxintronlen = atoi(optarg);
                break;
            case 'b':
			  	globalOptions.minintronlen = atoi(optarg);
                break;
            case 'c':
			  	globalOptions.maxSortesTest = atoi(optarg);
                break;
			case 'd':
			  	globalOptions.minId = atoi(optarg);
				break;
			case 'e':
			  	globalOptions.minCover = atoi(optarg);
			    break;
            case 'f':
			  	globalOptions.insertLimit = atoi(optarg);
                break;
			case 'g':
			  	globalOptions.uniqthresh = atof(optarg);
			    break;
			case 'i':
			  	globalOptions.commongenefile = optarg;
			    break;
			case 'j':
			  	globalOptions.pairbedfile = optarg;
			    break;
			case 'k':
			  	globalOptions.paired = 1;
			    break;
			case 'l':
			  	globalOptions.uniq = 1;
			    break;
			case 'm':
			  	globalOptions.nointrons = 1;
			    break;
			case 'n':
			  	globalOptions.best = 1;
			    break;
			case 'v':
			  	globalOptions.verbose = 1;
			    break;

			case 'h':   // fall-through is intentional
            case '?':
			  	displayUsage(argc, argv);
    			exit(EXIT_FAILURE);
                // break;              
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
  unordered_map<string, int> qNameStems;
  unordered_map<string, int>::iterator it;
  time_t tStart;
  time_t tEnd;
  time_t tElapsed; 
  vector<CigarOp> cigar;
  char cigarType;
  uint32_t cigarLength;
  uint32_t qLength;
  uint32_t RefID;
  uint32_t editDistance;
  string rName;
  string qName;
  string qNameStem = ""; 
  string oldqNameStem = "";
  // Qali variables
  string qSuffix;
  bool strand;
  int rstart, rend;
  std::stringstream ss_rstart, ss_rend, ss_percid, ss_coverage;
  int cigarSize;
  int sumMandI;
  float coverage;
  float percid;
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
  bool paired = globalOptions.paired;
  int maxSortesTest = globalOptions.maxSortesTest;

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

		/////////////////////////////////////////////////////////////////
		//     Checking of file pairedness STARTS here 
 		/////////////////////////////////////////////////////////////////
	  	qName = al.Name;
	  	qNameStem = qName; 
		
		if (paired)
		  {	
			qNameStem = qName.substr(0, qName.find("/"));
			qSuffix = qName.substr(qName.find("/")+1, qName.length());	
		  } 

		// Printing info
		printf("%s, %s, %s, %s\n", qName.c_str(), qNameStem.c_str(), qSuffix.c_str(), oldqNameStem.c_str());
 		
		// Verifying file is sorted in ascending order
		if (oldqNameStem.compare("") && oldqNameStem.compare(qNameStem))   
		  { 
			// Checking whether 10th field is sorted in ascending order
			if (line <= maxSortesTest && qNameStems[qNameStem])   
			  {
				cerr << "Input file not sorted by query name!" << qNameStem.c_str() << 
				  " occurred previously. Set LC_ALL=C and sort -k 10,10" << endl; 	
				exit(1);
			  }		  	
			/* 	Calling of: "processQuery() if (@qali);" goes here */
		  }  // end outer if
		/////////////////////////////////////////////////////////////////
		//     Checking of file pairedness continues in next section 
 		/////////////////////////////////////////////////////////////////


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


     	/* push @qali, [$_, $targetname, $qsuffix, $strand, $tstart, $tend, $percid, $coverage]; 
		   Generating array of strings to do equivalent to: */
		ss_rstart << rstart; ss_rend << rend; ss_percid << percid; ss_coverage << coverage;
		string qali[QALISIZE] = {rName, qNameStem, qSuffix, (char*)strand, ss_rstart.str(), 
								 ss_rend.str(), ss_percid.str(), ss_coverage.str()};

		/////////////////////////////////////////////////////////////////
		//     Checking of file pairedness continues here 
		/////////////////////////////////////////////////////////////////
		oldqNameStem = qNameStem;
		if (line <= maxSortesTest)
		  {
			qNameStems[qNameStem] = 1;
		  }
		/////////////////////////////////////////////////////////////////
		//     Checking of file pairedness ENDS here 
		/////////////////////////////////////////////////////////////////

	
		// Save in case all criteria was satisfied
		writer.SaveAlignment(al);

	} // end while	
  
	/* 	Calling of: "processQuery() if ($qnamestem ne "");" goes here */


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
	cout << "Cmd line: " << endl;
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


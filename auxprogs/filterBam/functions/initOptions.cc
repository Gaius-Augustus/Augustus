/*	Auxiliary function for the example of how to use getopt_long() 


	Created: 27-September-2011
	Last modified:  9-January-2012
*/

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h> 
#include <getopt.h> 

// Default options
bool best = false;
bool help = false;
bool noIntrons = false;
bool paired = false;
bool uniq = false;
bool verbose = false;
int insertLimit = 10;
int maxIntronLen = 500000;
int maxSortesTest = 100000;
int minCover = 80;
int minId = 92;
int minIntronLen = 35;
float uniqThresh = 0.96;
const char* commonGeneFile;
const char* inputFile;
const char* outputFile;
const char* pairbedFile;


using namespace std;

// Definition of global variables
static const char *optString = "a:b:c:d:e:f:g:i:j:k:l:m:n:o:p:q:h?";
extern int opterr; // Display error if opterr=0

#ifndef GLOBALOPTIONS_T
#define GLOBALOPTIONS_T
struct globalOptions_t {
	bool best;
	bool help;
	bool noIntrons;
	bool paired;
	bool uniq;
	bool verbose;  	
	int insertLimit;
	int maxIntronLen;
	int maxSortesTest; 
	int minCover;
	int minId;
	int minIntronLen;
	float uniqThresh;
	const char* commonGeneFile;
  	const char* inputFile;
  	const char* outputFile;
	const char* pairbedFile;
};
#endif



// longOpts array
static const struct option longOpts[] = {
    { "best", no_argument, NULL, 'a' },
    { "help", no_argument, NULL, 'h' }, 
    { "noIntrons", no_argument, NULL, 'b' },
    { "paired", no_argument, NULL, 'c' },
    { "uniq", no_argument, NULL, 'd' },
    { "verbose", no_argument, NULL, 'e' },
	{ "insertLimit", required_argument, NULL, 'f' },
    { "maxIntronLen", required_argument, NULL, 'g' },
	{ "maxSortesTest", required_argument, NULL, 'i' },
    { "minCover", required_argument, NULL, 'j' },
  	{ "minId", required_argument, NULL, 'k' },
	{ "minIntronLen", required_argument, NULL, 'l' },
    { "uniqThresh", required_argument, NULL, 'm' }, 	 
    { "commonGeneFile", required_argument, NULL, 'n' }, 
    { "in", required_argument, NULL, 'o' }, 
    { "out", required_argument, NULL, 'p' }, 
    { "pairbedFile", required_argument, NULL, 'q' },
    { NULL, no_argument, NULL, 0 }
};


// Display usage when --help
void displayUsage(int argc, char *argv[])
{
	cout <<  " Usage: " << argv[0] << " --in in.bam --out out.bam [options]\n";
  	cout <<  "--------------------------------------------------" << endl;
	cout <<  " Available options are                            " << endl;
  	cout <<  "--------------------------------------------------" << endl;
	cout <<  "  PREREQUISITE: file must be sorted lexicographically by 'queryname', with e.g.\n" <<
		   	 "                sort -k 1,1 but be aware: LC_ALL may have to be set to C\n" << 
			 " 		because sort ignores characters like ':' \n";
  	cout <<  "--------------------------------------------------" << endl;
	cout <<  "  --best             output all best matches that satisfy minId and minCover (default " << 
					best << ")" << endl;
	cout <<  "  --help             display this menu" << endl;
	cout <<  "  --noIntrons        do not allow longer gaps -for RNA-RNA alignments- (default " << 
					noIntrons << endl;
	cout <<  "  --paired           require that paired reads are on opposite strands of same target" << endl;
	cout <<  "  \t\t     (default " << paired << ")" << endl;
	cout <<  "                     NOTE: " << endl;
	cout <<  "	\t\t if option 'paired' is used then it expects .f,.r or /1,/2" << endl; 
	cout <<  "	\t\t suffixes of mate pairs" << endl;
	cout <<  "  --uniq             take only best match and only, when second best is much worse " << 
					" (default " << uniq << ")" << endl;
	cout <<  "  --verbose          output debugging info (default " << verbose << ")" << endl;
	cout <<  "  --insertLimit=n    maximum assumed size of inserts (default " << insertLimit << ")" << endl;
	cout <<  "  --maxIntronLen=n   maximal separation of paired reads (default " << maxIntronLen << 
					")" << endl;
	cout <<  "  --maxSortesTest=n  maximal sortedness (default " << maxSortesTest << ")" << endl;
	cout <<  "  --minCover=n       minimal percentage of coverage of the query read (default " << 
			 		minCover << ")" << endl;
	cout <<  "  --minId=n          minimal percentage of identity (default " << minId << ")" << endl;
	cout <<  "  --minIntronLen=n   minimal     ''     ''   ''    ''   (default " << minIntronLen << 
					")" << endl;
	cout <<  "  --uniqThresh=n     threshold % for uniq, second best must be at most this" << endl;
	cout <<  "                     fraction of best (default " << uniqThresh << ") " << endl;
	cout <<  "  --commonGeneFile=s file name in which to write cases where one read maps to \n" <<
	         "                     several different genes" << endl;
	cout <<  "  --pairbedFile=s    file name of pairedness coverage:" << endl;
	cout <<  "                     options:"  << endl;
	cout <<  "	\t\t a .bed format file in which for each position the number of" << endl;
	cout <<  "	\t\t filtered read pairs is reported that contain the position in" << endl; 
	cout <<  "  \t\t\t or between the reads" << endl;
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
	globalOptions.best = best;
	globalOptions.help = help;
	globalOptions.noIntrons = noIntrons;
	globalOptions.paired = paired;
	globalOptions.uniq = uniq;
	globalOptions.verbose = verbose;
	globalOptions.insertLimit = insertLimit;
	globalOptions.maxIntronLen = maxIntronLen;
	globalOptions.maxSortesTest = maxSortesTest; 
	globalOptions.minCover = minCover;
	globalOptions.minId = minId;
	globalOptions.minIntronLen = minIntronLen;
	globalOptions.uniqThresh = uniqThresh; 
	globalOptions.commonGeneFile = "";
	globalOptions.inputFile = "";
	globalOptions.outputFile = "";
	globalOptions.pairbedFile = "";

	// Capturing options
	opt = getopt_long_only(argc, argv, optString, longOpts, &longIndex);
	// Scanning through options
    while( opt != -1 ) {
        switch( opt ) {                
			case 'a'	:	globalOptions.best = true;					break;
			case 'b'	:	globalOptions.noIntrons = true;				break;
			case 'c'	:	globalOptions.paired = true;			    break;
			case 'd'	:	globalOptions.uniq = true;					break;
			case 'e'	:	globalOptions.verbose = true;				break;
            case 'f'	:	globalOptions.insertLimit = atoi(optarg);	break;
            case 'g'	:	globalOptions.maxIntronLen = atoi(optarg);	break;
            case 'i'	:	globalOptions.maxSortesTest = atoi(optarg);	break;
			case 'j'	:	globalOptions.minCover = atoi(optarg);		break;
			case 'k'	:	globalOptions.minId = atoi(optarg);			break;
            case 'l'	:	globalOptions.minIntronLen = atoi(optarg);	break;
			case 'm'	:	globalOptions.uniqThresh = atof(optarg);	break;
			case 'n'	:	globalOptions.commonGeneFile = optarg;		break;
			case 'o'	:	globalOptions.inputFile = optarg;			break;
			case 'p'	:	globalOptions.outputFile = optarg;			break;
			case 'q'	:	globalOptions.pairbedFile = optarg;			break;
			case 'h':   // HELP: fall-through is intentional
            case '?':
				globalOptions.help = true;	
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

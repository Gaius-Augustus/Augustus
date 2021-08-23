/*	Auxiliary function for the example of how to use getopt_long() 


	Created: 27-September-2011
	Last modified: 27-January-20121
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
bool pairwiseAlignments = false;
int insertLimit = 10;
int maxIntronLen = 500000;
int maxSortesTest = 100000;
int minCover = 80;
int minId = 92;
float uniqThresh = 0.96;
const char* commonGeneFile;
const char* inputFile;
const char* outputFile;
const char* pairBedFile;


using namespace std;

// Definition of global variables
static const char *optString = "b:c:d:e:g:h?:i:l:n:o:p:q:s:u:v:w:x";

#ifndef GLOBALOPTIONS_T
#define GLOBALOPTIONS_T
struct globalOptions_t {
	bool best;
	bool help;
	bool noIntrons;
	bool paired;
	bool uniq;
	bool verbose;  
  	bool pairwiseAlignments;	
	int insertLimit;
	int maxIntronLen;
	int maxSortesTest; 
	int minCover;
	int minId;
	float uniqThresh;
	const char* commonGeneFile;
  	const char* inputFile;
  	const char* outputFile;
	const char* pairBedFile;
};
#endif



// longOpts array
static const struct option longOpts[] = {
    { "best", no_argument, NULL, 'b' },
    { "help", no_argument, NULL, 'h' }, 
    { "noIntrons", no_argument, NULL, 'n' },
    { "paired", no_argument, NULL, 'p' },
    { "uniq", no_argument, NULL, 'u' },
    { "verbose", no_argument, NULL, 'v' },
    { "insertLimit", required_argument, NULL, 'l' },
    { "maxIntronLen", required_argument, NULL, 'x' },
    { "maxSortesTest", required_argument, NULL, 's' },
    { "minCover", required_argument, NULL, 'c' },
    { "minId", required_argument, NULL, 'e' },
    { "uniqThresh", required_argument, NULL, 'q' }, 	 
    { "commonGeneFile", required_argument, NULL, 'g' }, 
    { "in", required_argument, NULL, 'i' }, 
    { "out", required_argument, NULL, 'o' }, 
    { "pairBedFile", required_argument, NULL, 'd' },
    { "pairwiseAlignments", no_argument, NULL, 'w' },
    { NULL, no_argument, NULL, 0 }
};


// Display usage when --help
void displayUsage(char *argv[])
{
    cout <<  " Usage: " << argv[0] << " --in in.bam --out out.bam [options]\n";
    cout <<  "--------------------------------------------------" << endl;
    cout <<  " PREREQUISITE: (File sorted by query name)       " << endl; 
    cout <<  " File must be sorted lexicographically by 'queryname', with e.g.\n" << 
             "\n  1) sort -k 1,1 [be aware: 'export LC_ALL=C' might be used " << 
             " because sort ignores characters like ':' \n" << 
             "      Note: bear in mind that this will require converting your BAM file into SAM.\n" << 
             "\n  2) samtools and bamtools provide facilities to do the sorting,\n" <<
             "      but they are not guaranteed to work because of the problem mentioned above.\n" << 
             "\n  3) In the case of samtools, the command is: 'samtools sort [-n] file.bam'\n" <<  
             "      [-n] should sort by query name, just as 'sort -k 10,10' would do in a PSL file.\n" << 
             "      Without options, the sorting will be done by reference name and target coordinate,\n" << 
             "      just as a 'sort -n -k 16,16 | sort -k 14,14' would do with PSL.\n" << 
             "      For more information check the man page included in samtools distribution.\n" <<   
             "\n  4) bamtools can also sort bam files 'bamtools sort -queryname -in file.bam'," << 
             "\n      but only provides the option to do it by queryname." << endl;
    cout <<  "--------------------------------------------------" << endl;
    cout <<  " Available options are                            " << endl;
    cout <<  "--------------------------------------------------" << endl;
    cout <<  "  --uniq               keep only the best match, remove all matches, if the second best is not much worse" << endl;
    cout <<  "  --uniqThresh f       threshold % for uniq, second best must be lower than this" << endl;
    cout <<  "                       fraction of best to keep the best match (default " << uniqThresh << ")" << endl;
    cout <<  "  --best               output all best matches that satisfy minId and minCover" << endl;
    cout <<  "  --minId n            minimal percentage of identity (default " << minId << ")" << endl;
    cout <<  "  --minCover n         minimal percentage of coverage of the query read (default " << minCover << ")" << endl;
    cout <<  "  --noIntrons          do not allow longer gaps -for RNA-RNA alignments-" << endl;
    cout <<  "  --insertLimit n      maximum assumed size of inserts (default " << insertLimit << ")" << endl;
    cout <<  "  --maxSortesTest n    test if input file is sorted by query name for this number of alignments (default " << maxSortesTest << ")" << endl;
    cout <<  "  --paired             require that paired reads are on opposite strands of same target" << endl;
    cout <<  "                       Requires alignment names to contain the suffixes /1,/2 or /f,/r" << endl;
    cout <<  "  --pairwiseAlignments use in case alignments were done in pairwise fashion" << endl;
    cout <<  "  --maxIntronLen n     maximal separation of paired reads (default " << maxIntronLen << ")" << endl;
    cout <<  "  --pairBedFile s      file name of pairedness coverage:" << endl;
    cout <<  "                         a .bed format file in which for each position the number of" << endl;
    cout <<  "                         filtered read pairs is reported that contain the position in" << endl; 
    cout <<  "                         or between the reads" << endl;
    cout <<  "  --commonGeneFile s   file name in which to write cases where one read maps to" << endl;
    cout <<  "                       several different genes" << endl;
    cout <<  "  --verbose            output debugging info" << endl;
    cout <<  "  --help               display this menu" << endl;
}


// Options initialization
globalOptions_t initOptions(int argc, char *argv[])
{
  	int opt = 0, longIndex;
	struct globalOptions_t globalOptions;

	// Display usage if only ./program is provided
	if (argc==1)
	  {
		displayUsage(argv);
		exit(EXIT_FAILURE);
	  }

    // Initialize globalOptions with default values
	globalOptions.best = best;
	globalOptions.help = help;
	globalOptions.noIntrons = noIntrons;
	globalOptions.paired = paired;
	globalOptions.uniq = uniq;
	globalOptions.verbose = verbose;
	globalOptions.pairwiseAlignments = verbose;
	globalOptions.insertLimit = insertLimit;
	globalOptions.maxIntronLen = maxIntronLen;
	globalOptions.maxSortesTest = maxSortesTest; 
	globalOptions.minCover = minCover;
	globalOptions.minId = minId;
	globalOptions.uniqThresh = uniqThresh; 
	globalOptions.commonGeneFile = "";
	globalOptions.inputFile = "";
	globalOptions.outputFile = "";
	globalOptions.pairBedFile = "";

	// Capturing options
	opt = getopt_long_only(argc, argv, optString, longOpts, &longIndex);
	// Scanning through options
    while( opt != -1 ) {
        switch( opt ) {                
            case 'b'	:	globalOptions.best = true;                  break;
            case 'n'	:	globalOptions.noIntrons = true;             break;
            case 'p'	:	globalOptions.paired = true;                break;
            case 'u'	:	globalOptions.uniq = true;                  break;
            case 'v'	:	globalOptions.verbose = true;               break;
            case 'w'	:	globalOptions.pairwiseAlignments = true;    break;
            case 'l'	:	globalOptions.insertLimit = atoi(optarg);   break;
            case 'x'	:	globalOptions.maxIntronLen = atoi(optarg);  break;
            case 's'	:	globalOptions.maxSortesTest = atoi(optarg); break;
            case 'c'	:	globalOptions.minCover = atoi(optarg);      break;
            case 'e'	:	globalOptions.minId = atoi(optarg);         break;
            case 'q'	:	globalOptions.uniqThresh = atof(optarg);    break;
            case 'g'	:	globalOptions.commonGeneFile = optarg;      break;
            case 'i'	:	globalOptions.inputFile = optarg;           break;
            case 'o'	:	globalOptions.outputFile = optarg;          break;
            case 'd'	:	globalOptions.pairBedFile = optarg;         break;
            case 'h':   // HELP: fall-through is intentional
            case '?':
				globalOptions.help = true;	
			  	displayUsage(argv);
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

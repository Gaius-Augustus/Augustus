/**********************************************************************
 * Load the sequences from a flat file in fasta format into a database.
 * file:    faload2db.cc
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 *
 * authors: Mario Stanke, mario.stanke@uni-greifswald.de
 *
 * date    |   author      |  changes
 * --------|---------------|------------------------------------------
 * 30.08.12| Mario Stanke  | creation of the file
 **********************************************************************/

// project includes
//#include "types.hh"

// standard C/C++ includes
#include <string>
#include <iostream>     /* for printf */
#include <getopt.h>     /* for getopt_long; standard getopt is in unistd.h */
#include <stdlib.h>     /* for exit() */
using namespace std;

void printUsage();

int chunksize = 50000;

/*
 * main
 */
int main( int argc, char* argv[] ){
    int c;
    int help = 0;
    string species;
    string dbaccess;
    string genomefname;
    static struct option long_options[] = {
        {"species", 1, 0, 's'},
        {"dbaccess", 1, 0, 'd'},
        {"help", 0, 0, 'h'},
        {NULL, 0, NULL, 0}
    };
    int option_index = 0;
    while ((c = getopt_long(argc, argv, "s:d:h", long_options, &option_index)) != -1) {
        switch (c) {
        case 's':
	    species = optarg;
            break;
        case 'd':
	    dbaccess = optarg;
            break;
        case 'h':
            help = 1;
            break;
        default:
            break;
        }
    }
    if (optind < argc-1) {
        cerr << "More than one option without name: ";
        while (optind < argc)
            cerr << " " << argv[optind++];
	cerr << endl << "Only the sequence file name does not require a parameter name." << endl;
	if (!help)
	    printUsage();
	exit(1);
    }
    if (optind == argc-1)
	genomefname = argv[optind];
    else {
	cerr << "Missing sequence file name." << endl;
	printUsage();
	exit(1);
    }
    if (help) {
	printUsage();
	exit(1);
    }
    if (species.empty()){
	cerr << "Missing species name. Required parameter." << endl;
	printUsage();
	exit(1);
    }
    if (dbaccess.empty()){
	cerr << "Missing database access info. dbaccess is a required parameter." << endl;
	printUsage();
	exit(1);
    }
    
    cout << "species=" << species << endl
	 << "dbaccess=" << dbaccess << endl
	 << "genomefname=" << genomefname << endl;
}


void printUsage(){
    cout << "usage:\n\
faload2db [parameters] --species=SPECIES --dbaccess=dbname,host,user,passwd  inputfilename\n\
\n\
inputfilename refers to a genome file in FASTA format\n\
SPECIES is the same identifier as is used in the treefile and alnfile\n\
\n\
dbname,host,user,passwd are the name of the SQL database, the host name or IP, the database user and password\n\
When storing genomes of multiple organisms call this program repeatedly for each one.\n\
\n\
parameters:\n\
--help        print this usage info\n\
--chunksize   the sequences in the input genome are split into chunks of this size so\n\
              that subsequent retrievals of small sequence ranges do not require to read\n\
              the complete - potentially much longer - chromosome. (default " << chunksize << ")\n";
}

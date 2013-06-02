// ***************************************************************************
// bamtools.cpp (c) 2010 Derek Barnett, Erik Garrison
// Marth Lab, Department of Biology, Boston College
// ---------------------------------------------------------------------------
// Last modified: 12 October 2012 (DB)
// ---------------------------------------------------------------------------
// Integrates a number of BamTools functionalities into a single executable.
// ***************************************************************************

#include "bamtools_illumina.h"
#include "bamtools_version.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
using namespace BamTools;
using namespace std;

// bamtools subtool names
static const string ILLUMINA  = "illumina";

// bamtools help/version constants
static const string HELP          = "help";
static const string LONG_HELP     = "--help";
static const string SHORT_HELP    = "-h";
static const string VERSION       = "version";
static const string LONG_VERSION  = "--version";
static const string SHORT_VERSION = "-v";

// determine if string is a help constant
static bool IsHelp(char* str) {
    return ( str == HELP ||
             str == LONG_HELP ||
             str == SHORT_HELP );
}

// determine if string is a version constant
static bool IsVersion(char* str) {
    return ( str == VERSION ||
             str == LONG_VERSION ||
             str == SHORT_VERSION );
}

// subtool factory method
AbstractTool* CreateTool(const string& arg) {
  
    // determine tool type based on arg
    if ( arg == ILLUMINA )  return new IlluminaTool;

    // unknown arg
    return 0;
}

// print help info
int Help(int argc, char* argv[]) {
  
    // check for 'bamtools help COMMAND' to print tool-specific help message
    if (argc > 2) {
        
	// determine desired sub-tool
        AbstractTool* tool = CreateTool( argv[2] );

        // if tool known, print its help screen
        if ( tool ) return tool->Help();
    }

    // print general BamTools help message
    cerr << endl;
    cerr << "usage: bamtools [--help] COMMAND [ARGS]" << endl;
    cerr << endl;
    cerr << "Available bamtools commands:" << endl;
    cerr << "\tillumina         filters alignments coming from illumina paired-end reads" << endl;
    cerr << "See Augustus documentation for more information." << endl;
    cerr << endl;
    return EXIT_SUCCESS;
}

// print version info
int Version(void) {

    stringstream versionStream("");
    versionStream << BAMTOOLS_VERSION_MAJOR << "."
                  << BAMTOOLS_VERSION_MINOR << "."
                  << BAMTOOLS_VERSION_BUILD;

    cout << endl;
    cout << "illuminaFilter" << versionStream.str() << endl;
    cout << "Written by Tonatiuh Pena Centeno" << endl;
    cout << endl;
    return EXIT_SUCCESS;
}

// toolkit entry point
int main(int argc, char* argv[]) {

    // just 'bamtools'
    if ( (argc == 1) ) return Help(argc, argv);
    
    // 'bamtools help', 'bamtools --help', or 'bamtools -h'
    if ( IsHelp(argv[1]) ) return Help(argc, argv); 
    
    // 'bamtools version', 'bamtools --version', or 'bamtools -v'
    if ( IsVersion(argv[1]) ) return Version(); 
        
    // determine desired sub-tool, run if found
    AbstractTool* tool = CreateTool( argv[1] );
    if ( tool ) return tool->Run(argc, argv);

    // no tool matched, show help
    return Help(argc, argv);
}

/**********************************************************************
 * Bulk loading of sequences and hints into an SQLite database.
 * file:    load2sqlitedb.cc
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 *
 * authors: Mario Stanke, mario.stanke@uni-greifswald.de
 *
 * date    |   author      |  changes
 * --------|---------------|------------------------------------------
 * 24.06.14|Stefanie Koenig| creation of the file
 **********************************************************************/

// Project includes
#include "fasta.hh"
#include "hints.hh"
#include "sqliteDB.hh"
#include "projectio.hh"

// standard C/C++ includes
#include <string>
#include <iostream>
#include <fstream>
#include <getopt.h>     /* for getopt_long; standard getopt is in unistd.h */
#include <stdlib.h>     /* for exit() */
#include <exception>

#include <stdio.h>

using namespace std;

int chunksize = 50000;

void printUsage();

/*
 * main
 */
int main( int argc, char* argv[] ){
    int c;
    int help = 0;
    string species;
    string dbfile;
    string fastafile;

    static struct option long_options[] = {
        {"species", 1, 0, 's'},
        {"dbaccess", 1, 0, 'd'},
        {"help", 0, 0, 'h'},
	{"chunksize", 1, 0, 'c'},
        {NULL, 0, NULL, 0}
    };
    int option_index = 0;
    while ((c = getopt_long(argc, argv, "s:d:hc:", long_options, &option_index)) != -1) {
        switch (c) {
        case 's':
	    species = optarg;
            break;
        case 'd':
	    dbfile = optarg;
            break;
        case 'h':
            help = 1;
            break;
        case 'c':
	    chunksize = atoi(optarg);
            break;
        default:
            break;
        }
    }
    if (optind < argc-2) {
        cerr << "More than two options without name: ";
        while (optind < argc)
            cerr << " " << argv[optind++];
	cerr << endl << "Only the sequence/hints file name does not require a parameter name." << endl;
	if (!help)
	    printUsage();
	exit(1);
    }
    if (optind == argc-1)
	fastafile = argv[optind];
    else {
	cerr << "Missing sequence/hints file name." << endl;
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
    if (dbfile.empty()){
	cerr << "Missing database file. dbaccess is a required parameter." << endl;
	printUsage();
	exit(1);
    }
    if (chunksize < 2){
	cerr << "Chunksize too small (" << chunksize << "). Should be roughly in the oder of a gene's length." << endl;
	printUsage();
	exit(1);
    }
    if (chunksize > 1000000){
	cerr << "Chunksize too big (" << chunksize << "). " << endl;
	printUsage();
	exit(1);
    }

    // open database
    SQLiteDB db(dbfile.c_str(),crw);
    try {
	db.createTableSpeciesnames();
	db.createTableSeqnames();

	ifstream ifstrm;
	ifstrm.open(fastafile.c_str());
	if( !ifstrm ){
	    cerr << "Could not open input file " << fastafile << endl;
	    exit(1);
	}
	// if input file is in Fasta format, try to load sequences
	if(isFasta(ifstrm)){
	    cout << "Looks like " << fastafile << " is in fasta format." << endl;
	    db.createTableGenomes();

	    // drop index on genomes table for faster insertion
	    db.exec("DROP INDEX IF EXISTS gidx;");

	    int seqCount = 0, chunkCount = 0;
	    unsigned int lenCount = 0;
	    int speciesid = db.getSpeciesID(species);
	
	    db.beginTransaction();

	    Statement stmt1(&db); 
	    stmt1.prepare("INSERT INTO genomes VALUES(?1,?2,?3,?4,?5,?6,?7);"); // prepared statement for inserting sequence chunks

	    Statement stmt2(&db);
	    stmt2.prepare("INSERT INTO seqnames (speciesid,seqname) VALUES(?1,?2);"); // prepared statement for inserting fasta headers
	    	    
	    while (ifstrm){
		char *name = NULL;
		readFastaHeader(ifstrm, name);
		try{
		    stmt2.bindInt(1,speciesid);
		    stmt2.bindText(2,name);
		    stmt2.step();
		    stmt2.reset();
		}catch(const char* err){
		    cerr << "failed inserting sequence "<< name << " for species "<< species <<" (ID="<<speciesid<<")"<< endl;
		    cerr << "Is it possible that the sequence is already in the database?" << endl;
		    cerr << "If so, please delete it, before reloading it" << endl;
		    throw err;
		}
		int seqnr = db.lastInsertID();
		int start = 0;
		string line;
		
  
		while (ifstrm && ifstrm.peek( ) != '>'){
		    streampos file_start = ifstrm.tellg(), file_end = ifstrm.tellg();
		    int length = 0;
		    while(length < chunksize && ifstrm && ifstrm.peek( ) != '>'){
			if(getline(ifstrm, line)){
			    length+=line.size();
			    file_end = ifstrm.tellg();
			}
		    }
		    stmt1.bindInt(2,seqnr);
		    stmt1.bindInt(3,speciesid);
		    stmt1.bindInt(4,start);
		    stmt1.bindInt(5,start+length-1);
		    stmt1.bindInt64(6,(uint64_t)file_start);
		    stmt1.bindInt(7,file_end-file_start);
		    stmt1.step();
		    stmt1.reset();
  
		    chunkCount++;
		    start += chunksize;
		    lenCount +=length;
		}
		delete name;
		seqCount++;
	    }
	    db.endTransaction();
		
	    // rebuild index on genomes table
	    db.exec("CREATE INDEX gidx ON genomes(speciesid,seqnr,start,end);");

	    if (seqCount > 0)
		cout << "Inserted " << chunkCount << " chunks of " << seqCount << " sequences (total length "
		     << lenCount << " bp)." << endl;
	    else
		cout << "No sequences found. Nothing inserted into database." << endl;
	        
	}
	else if(isGFF(ifstrm)){ // if input file is in GFF format, try to load hints
	    cout << "Looks like " << fastafile << " is in gff format." << endl;
	    db.createTableHints();
	    db.createTableFeatureTypes();
	    int hintCount=0; // number of hints inserted into the database
	    
	    // drop index on hints table (for faster insertion)
	    db.exec("DROP INDEX IF EXISTS hidx");

	    // bulk insert of all hints
	    db.beginTransaction();

	    Statement stmt(&db);
	    stmt.prepare("INSERT INTO hints (speciesid,seqnr,source,start,end,score,type,strand,frame,priority,grp,mult,esource) \
		       select speciesid,seqnr,?3,?4,?5,?6,?7,?8,?9,?10,?11 \
                       ,?12,?13 from speciesnames natural join seqnames \
                       where speciesname=?1 and seqname=?2;");
	    
	    while(ifstrm){
		Feature f;
		ifstrm >> f >> comment >> ws;
		try{
		    stmt.bindText(1,species.c_str());
		    stmt.bindText(2,f.seqname.c_str());
		    stmt.bindText(3,f.source.c_str());
		    stmt.bindInt(4,f.start);
		    stmt.bindInt(5,f.end);
		    stmt.bindDouble(6,f.score);
		    stmt.bindInt(7,f.type);
		    switch (f.strand) {
		    case plusstrand : stmt.bindText(8,"+"); break;
		    case minusstrand : stmt.bindText(8,"-"); break;
		    default : stmt.bindText(8,".");
		    }
		    switch (f.frame) {
		    case 0 :  stmt.bindText(9,"0"); break;
		    case 1 :  stmt.bindText(9,"1"); break;
		    case 2 :  stmt.bindText(9,"2"); break;
		    default :  stmt.bindText(9,".");
		    }
		    stmt.bindInt(10,f.priority);
		    stmt.bindText(11,f.groupname.c_str());
		    stmt.bindInt(12,f.mult);
		    stmt.bindText(13,f.esource.c_str());
		    stmt.step();
		    stmt.reset();
		}catch(const char* error){
		    cerr << "insert failed on\n" << f << endl;
		    cerr << error << endl;
		    exit(1);
		}
		if(db.numChanges() != 1){
		    cerr << "insert failed on\n" << f << endl;
		    cerr << "hints can only be inserted for sequences in the database." << endl;
			exit(1);
		}
		hintCount++;
	    }	
	    db.endTransaction();
	    
	    // rebuild index on hints table
	    db.exec("CREATE INDEX hidx ON hints(speciesid,seqnr,start,end);");

	    if(hintCount > 0)
		cout << "Inserted " << hintCount<< " hints for species " << species << endl; 
	    else
		cout << "No hints found. Nothing inserted into database." << endl;
	}
	else{
	    cout << fastafile << " is neither in gff nor in fasta format." << endl;
	}
	ifstrm.close();   
    } catch(const char* err){
        cerr << "\n" <<  argv[0] << ": ERROR\n\t" << err << "\n\n";
	exit(1);
    }
}

void printUsage(){
    cout << "usage:\n\
load2sqlitedb [parameters] --species=SPECIES --dbaccess=database.db fastafile\n\
\n\
fastafile refers to a genome file in FASTA format or a hints file in GFF format\n\
SPECIES is the same identifier as is used in the treefile and alnfile parameters to augustus.\n\
\n\
database.db is the name of the database that will be opened or created if it does not exist already.\n\
When storing genomes/hints of multiple organisms call this program repeatedly for each one.\n\
A single table with the structure is created.\n\
\n\
parameters:\n\
--help        print this usage info\n\
--chunksize   this option is only relevant when loading a sequence file\n\
              the sequences in the input genome are split into chunks of this size so\n\
              that subsequent retrievals of small sequence ranges do not require to read\n\
              the complete - potentially much longer - chromosome. (<= 1000000, default " << chunksize << ")\n\
\n\
example:\n\
     load2sqlitedb --species=chicken --dbaccess=chicken.db chickengenome.fa\n\
     load2sqlitedb --species=chicken --dbaccess=chicken.db chickenhints.gff\n";
}

/*
 * load2sqlitedb.cc
 *
 * License: Artistic License, see file LICENSE.TXT or 
 *          https://opensource.org/licenses/artistic-license-1.0
 * 
 * Description: Bulk loading of sequences and hints into an SQLite database.
 */

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
    int noIdx = 0;
    int makeIdx = 0;
    bool clean = false;
    string species;
    string dbfile;
    string fastafile;

    static struct option long_options[] = {
        {"species", 1, 0, 's'},
        {"dbaccess", 1, 0, 'd'},
        {"help", 0, 0, 'h'},
	{"chunksize", 1, 0, 'c'},
        {"noIdx", 0, 0, 'i'},
        {"makeIdx", 0, 0, 'm'},
	{"clean", 0, 0, 'r'},
        {NULL, 0, NULL, 0}
    };
    int option_index = 0;
    while ((c = getopt_long(argc, argv, "s:d:hc:imr", long_options, &option_index)) != -1) {
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
	case 'i':
            noIdx = 1;
	    break;
	case 'm':
            makeIdx = 1;
	    break;
	case 'r':
	    clean = true;
        default:
            break;
        }
    }
    if (dbfile.empty()){
	cerr << "Missing database file. dbaccess is a required parameter." << endl;
	printUsage();
	exit(1);
    }
    if (noIdx && makeIdx){
	cerr << "You can only use one of the options --noIdx or --makeIdx at a time." << endl;
	exit(1);
    }
    if(makeIdx){ // only build indices
	SQLiteDB db(dbfile.c_str(),crw);
	try{
	    cout << "Creating indices on genomes tables."<< endl;	    
	    db.exec("CREATE INDEX IF NOT EXISTS gidx ON genomes(speciesid,seqnr,start,end);");
	    if(db.tableExists("hints")){
		cout << "Creating indices on hints tables."<< endl;	    		
		db.exec("CREATE INDEX IF NOT EXISTS hidx ON hints(speciesid,seqnr,start,end);");
	    }
	}catch(const char* err){
	    cerr << "Failed creating indices." << endl;
	    cerr << err << endl;
	    exit(1);
	}
	exit(0);
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
    if (optind == argc-1){
	fastafile = argv[optind];
	// expand the ~ to the $Home directory
	if (fastafile.length()>0 && fastafile[0]=='~')
	    fastafile.replace(0,1,getenv("HOME"));
    }
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
    if (chunksize < 2){
	cerr << "Chunksize too small (" << chunksize << "). Should be roughly in the order of a gene's length." << endl;
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
	db.createTableGenomes();
	db.createTableHints();
	db.createTableFeatureTypes();

	ifstream ifstrm;
	ifstrm.open(fastafile.c_str());
	if( !ifstrm ){
	    cerr << "Could not open input file " << fastafile << endl;
	    exit(1);
	}
	// if input file is in Fasta format, try to load sequences
	if(isFasta(ifstrm)){
	    cout << "Looks like " << fastafile << " is in fasta format." << endl;

	    // drop index on genomes table for faster insertion
	    db.exec("DROP INDEX IF EXISTS gidx;");

	    int seqCount = 0, chunkCount = 0;
	    unsigned int lenCount = 0;
	    int speciesid = db.getSpeciesID(species,clean);
	
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
		    cerr << "If you want to replace the existing genome, please make a clean load with option --clean" << endl;
		    throw err;
		}
		int seqnr = db.lastInsertID();
		int start = 0;
		int length = 0;
		string line;
		char c;
		
		streampos file_start = ifstrm.tellg();
		streampos file_end = file_start;
		streampos pos = file_end;
	
		while (ifstrm && ifstrm.peek() != '>' && ifstrm.peek() != EOF){
		    ifstrm.get(c);
		    pos += streamoff(1);
		    if(isalpha(c)){
			length++;
			file_end = pos;
		    }
		    if(length >= chunksize){
			stmt1.bindInt(2,seqnr);
			stmt1.bindInt(3,speciesid);
			stmt1.bindInt(4,start);
			stmt1.bindInt(5,start+length-1);
			stmt1.bindInt64(6,(uint64_t)file_start);
			stmt1.bindInt(7,file_end-file_start);
			stmt1.step();
			stmt1.reset();
  
			chunkCount++;
			lenCount +=length;

			start += length;
			// next file_start position is the character preceeding
			// the next non-whitespace character
			while (ifstrm && !isalpha(ifstrm.peek()) && ifstrm.peek() != '>' && ifstrm.peek() != EOF){
			    ifstrm.get(c);
			    pos += streamoff(1);
			}
			file_start = pos;
			length = 0;
		    }
		}			
		// last chunk
		if(length > 0){
		    stmt1.bindInt(2,seqnr);
		    stmt1.bindInt(3,speciesid);
		    stmt1.bindInt(4,start);
		    stmt1.bindInt(5,start+length-1);
		    stmt1.bindInt64(6,(uint64_t)file_start);
		    stmt1.bindInt(7,file_end-file_start);
		    stmt1.step();
		    stmt1.reset();
		
		    chunkCount++;
		    lenCount +=length;
		}
		delete name;
		seqCount++;
		
	    }
	    db.endTransaction();
		
	    // rebuild index on genomes table
	    if(!noIdx)
		db.exec("CREATE INDEX gidx ON genomes(speciesid,seqnr,start,end);");

	    if (seqCount > 0)
		cout << "Inserted " << chunkCount << " chunks of " << seqCount << " sequences (total length "
		     << lenCount << " bp)." << endl;
	    else
		cout << "No sequences found. Nothing inserted into database." << endl;
	        
	}
	else if(isGFF(ifstrm)){ // if input file is in GFF format, try to load hints
	    cout << "Looks like " << fastafile << " is in gff format." << endl;

	    int hintCount=0; // number of hints inserted into the database
	    
	    // drop index on hints table (for faster insertion)
	    db.exec("DROP INDEX IF EXISTS hidx");

	    if(clean){
	      int id = db.getSpeciesID(species,false,true);
	      if(id >= 0){
		db.deleteHints(id); // delete existing hints from DB
	        cout << "Deleted existing hints for " << species << " from database." << endl;
	      }
	    }
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
	    if(!noIdx)
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
    exit(0);
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
--noIdx       use this flag to suppress the building of indices on the database tables.\n\
              If you are going to load several genomes and/or hint files in a row, this option\n\
              is recommended to speed up the loading. But make sure to build indices with\n\
              --makeIdx after all genomes/hints are loaded. Otherwise, data retrieval operations\n\
              can be very slow.\n\
--makeIdx     use this flag to build the indices on the database tables after loading several\n\
              genomes and/or hint files with --noIdx. Only call this once for all species, e.g.\n\
              load2sqlitedb --makeIdx --dbaccess=database.db\n\
--clean       makes a clean load deleting existing hints/genome for the species from the database.\n\
              When called with a gff file, only the hints for the species are delete, but not the genome.\n\
              When called with a fasta file, both hints and genome for the species are deleted.\n\
\n\
examples:\n\
     load2sqlitedb --species=chicken --dbaccess=chicken.db chickengenome.fa\n\
     load2sqlitedb --species=chicken --dbaccess=chicken.db chickenhints.gff\n";
}

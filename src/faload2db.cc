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
#include "fasta.hh"

// standard C/C++ includes
#include <string>
#include <iostream>
#include <fstream>
#include <getopt.h>     /* for getopt_long; standard getopt is in unistd.h */
#include <stdlib.h>     /* for exit() */
#include <mysql++.h>
#include <exception>
#include <ssqls.h>

using namespace std;

//sql_create_5(genomes,1,5,int,seqid,string,sequence,int,start,int,end,string,species)

int chunksize = 50000;
mysqlpp::Connection con;

void printUsage();
void connectDB(string dbaccess);
void createTable();
int insertSeq(char *sequence, char* name, int length, string species);

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
	    dbaccess = optarg;
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

    try {
	connectDB(dbaccess);
    } catch (const char *m) {
	cerr << "Database connection error:" << endl << m << endl;
	exit(1);
    }
    try {
	createTable();

	ifstream ifstrm;
	ifstrm.open(genomefname.c_str());
	if( !ifstrm )
	    throw string("Could not open input file \"") + genomefname + "\"!";

	char *sequence = NULL, *name = NULL;
	int length, seqCount = 0, chunkCount = 0, lenCount = 0;
	readOneFastaSeq(ifstrm, sequence, name, length);
	while (sequence){
	    chunkCount += insertSeq(sequence, name, length, species);
	    seqCount++;
	    lenCount += length;
	    delete sequence;
	    delete name;
	    sequence = name = NULL;
	    readOneFastaSeq(ifstrm, sequence, name, length);
	}
	if (seqCount > 0)
	    cout << "Inserted " << chunkCount << " chunks of " << seqCount << " sequences (total length "
		 << lenCount << " bp)." << endl;
	else
	    cout << "No sequences found. Nothing inserted into database." << endl;
    } catch( string err ){
        cerr << "\n" <<  argv[0] << ": ERROR\n\t" << err << "\n\n";
	exit(1);
    } 
}


void printUsage(){
    cout << "usage:\n\
faload2db [parameters] --species=SPECIES --dbaccess=dbname,host,user,passwd  inputfilename\n\
\n\
inputfilename refers to a genome file in FASTA format\n\
SPECIES is the same identifier as is used in the treefile and alnfile parameters to augustus.\n\
\n\
dbname,host,user,passwd are the name of the SQL database, the host name or IP, the database user and password\n\
When storing genomes of multiple organisms call this program repeatedly for each one.\n\
A single table with the structure\n\
\
is created.\n\
\n\
parameters:\n\
--help        print this usage info\n\
--chunksize   the sequences in the input genome are split into chunks of this size so\n\
              that subsequent retrievals of small sequence ranges do not require to read\n\
              the complete - potentially much longer - chromosome. (<= 1000000, default " << chunksize << ")\n\
\n\
example:\n\
     faload2db --species=chicken --dbaccess=birds,localhost,mario,dF$n.E chickengenome.fa\n\
\n\
Example code for database creation before calling faload2db:\n\
mysql -u root -p\n\
create database birds;\n\
select password('dF$n.E');\n\
create user `mario`@`%` identified by password '*72E7F76393830492EB3C58AD730188708BD72DE1'; /* or whatever the password code is*/\n\
grant all privileges on birds.* to mario@'%';\n";
}

void connectDB(string dbaccess){
    string db_name, host, user, passwd;

    string::size_type start=0, end;
    end = dbaccess.find(','); // string 'dbaccess' is delimited by ','
    if (end == string::npos)
	throw ("Database name missing in dbaccess.");
    db_name = dbaccess.substr(start, end-start);
    start = end + 1;
    end = dbaccess.find(',', start);
    if (end == string::npos)
	throw ("Host name missing in dbaccess.");
    host = dbaccess.substr(start, end-start);
    start = end + 1;
    end = dbaccess.find(',', start);
    if (end == string::npos)
	throw ("User name missing in dbaccess.");
    user = dbaccess.substr(start, end-start);
    start = end + 1;
    end = dbaccess.find(',', start);
    if (end == string::npos)
	passwd = dbaccess.substr(start);
    else 
	passwd = dbaccess.substr(start, end-start); // in case dbaccess also includes the port number
    
    try {
	cout << "Trying to connect to database " << db_name << " on server "
	     << host << " as user " << user << " using password " << passwd << " ..." << endl;
	con.connect(db_name.c_str(), host.c_str(), user.c_str(), passwd.c_str());
    } catch(const mysqlpp::BadQuery& e){
	cout << "Connection error: " << e.what() << endl;
    }
}

/*
 * create table `genomes` if it does not already exist
 */
void createTable(){
    mysqlpp::Query query = con.query("CREATE TABLE IF NOT EXISTS genomes (\
        seqid int(10) unsigned NOT NULL AUTO_INCREMENT,			\
        dnaseq longtext NOT NULL,					\
        start int(9) unsigned NOT NULL,\
        end int(9) unsigned NOT NULL,\
        species varchar(50),\
        PRIMARY KEY (seqid),\
        KEY region (species,start,end)\
      ) ENGINE=MyISAM DEFAULT CHARSET=latin1 MAX_ROWS=1000000 AVG_ROW_LENGTH=50000;");
   query.execute();
}

/*
 * insert one sequence into database
 * return the number of chunks that sequence was cut into
 */
int insertSeq(char *sequence, char *name, int length, string species){
    int chunks = 0;
    
    cout << "inserting " << name << " len=" << length << "..." << endl;

    mysqlpp::Query query = con.query();

    //genomes row(mysqlpp::null, sequence, 0, length-1, species);
    //query.insert(row);

    query << "INSERT INTO genomes (dnaseq,start,end,species) VALUES(\""
	  << sequence << "\",0," << length-1 << ",\"" << species << "\")";

    cout << "Executing" << endl << query.str() << endl;
    query.execute();

    return chunks;
}

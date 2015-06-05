/**********************************************************************
 * Load the sequences from a flat file in fasta format into a database.
 * file:    load2db.cc
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 *
 * authors: Mario Stanke, mario.stanke@uni-greifswald.de
 *
 * date    |   author      |  changes
 * --------|---------------|------------------------------------------
 * 30.08.12| Mario Stanke  | creation of the file
 * 20.11.13|Stefanie Koenig| loading extrinsic evidence into database
 **********************************************************************/

// Project includes
#include "fasta.hh"
#include "projectio.hh"
#include "hints.hh"

// standard C/C++ includes
#include <string>
#include <iostream>
#include <fstream>
#include <getopt.h>     /* for getopt_long; standard getopt is in unistd.h */
#include <stdlib.h>     /* for exit() */
#include <mysql++.h>
#include <exception>
#include <ssqls.h>

//#include <boost/iostreams/filtering_stream.hpp>
//#include <boost/iostreams/filter/gzip.hpp>
//#include <boost/iostreams/copy.hpp>

using namespace std;
//using boost::iostreams::filtering_istream;
//using boost::iostreams::gzip_decompressor;

//sql_create_6(genomes,1,6,int,seqid,string,sequence,string,seqname,int,start,int,end,string,species)

int chunksize = 50000;
mysqlpp::Connection con;

void printUsage();
void connectDB(string dbaccess);
void createTableGenomes();
void createTableSpeciesnames();
void createTableSeqnames();
void createTableHints();
void createTableFeatureTypes();
int getSpeciesID(string species);
int getSeqNr(char* seqname,int speciesid);
int insertSeqName(char* seqname,int speciesid);
int insertSeq(string sequence, char* name, int length, int speciesid);
void insertHint(Feature &f, string species);
void checkConsistency(); // checks consistency of the 'hints' data with the 'genomes' data
void removeDuplicates();

/*
 * main
 */
int main( int argc, char* argv[] ){
    int c;
    int help = 0;
    string species;
    string dbaccess;
    string filename;
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
	cerr << endl << "Only the sequence/hints file name does not require a parameter name." << endl;
	if (!help)
	    printUsage();
	exit(1);
    }
    if (optind == argc-1)
	filename = argv[optind];
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
	createTableSpeciesnames();
	createTableSeqnames();

	ifstream ifstrm;
	ifstrm.open(filename.c_str());
	if( !ifstrm )
	    throw string("Could not open input file \"") + filename + "\"!";

	/*filtering_istream zin;
	try {  
	    zin.push(gzip_decompressor());
	    zin.push(ifstrm);
	    zin.peek();
	    if (!zin)
		throw("Could not read first character assuming gzip format.");
	    cout << "Looks like " << filename << " is a gzip file. Deflating..." << endl;
        } catch (...) { // boost::iostreams::gzip_error& 
	    // not a gzip file or ill-formatted
	    zin.reset();
	    ifstrm.seekg(0);
	    zin.push(ifstrm);
	    }*/
	
	//if(isFasta(zin)){
	if(isFasta(ifstrm)){
	    cout << "Looks like " << filename << " is in fasta format." << endl;
	    createTableGenomes();
	    char *sequence = NULL, *name = NULL;
	    int length = 0, seqCount = 0, chunkCount = 0;
	    unsigned int lenCount = 0;
	    //readOneFastaSeq(zin, sequence, name, length);
	    readOneFastaSeq(ifstrm, sequence, name, length);
	    int speciesid = getSpeciesID(species);
	    while (sequence){
		chunkCount += insertSeq(sequence, name, length, speciesid);
		seqCount++;
		lenCount += length;
		delete sequence;
		delete name;
		sequence = name = NULL;
		//readOneFastaSeq(zin, sequence, name, length);
		readOneFastaSeq(ifstrm, sequence, name, length);
	    }
	    if (seqCount > 0)
		cout << "Inserted " << chunkCount << " chunks of " << seqCount << " sequences (total length "
		     << lenCount << " bp)." << endl;
	    else
		cout << "No sequences found. Nothing inserted into database." << endl;
	    
	}
	else if(isGFF(ifstrm)){
	    cout << "Looks like " << filename << " is in gff format." << endl;
	    createTableHints();
	    createTableFeatureTypes();
	    map<string,int> hintCount; // stores for each sequence the number of hints inserted into the database
	    while(ifstrm){
		Feature f;
		try{
		    ifstrm >> f >> comment >> ws;
		    if(f.type != -1){
			insertHint(f,species);
			hintCount[f.seqname]++;
		    }
		} catch (ProjectError e){}
	    }
	    if(!hintCount.empty()){
		cout << "inserted" << endl; 
		for(map<string,int>::iterator it= hintCount.begin(); it != hintCount.end(); it++){
		    cout << it->second << " hints for " << it->first << endl;
		}
		//removeDuplicates(); dublicastes can only be removed when multiplicity is updated
	    }
	    else{
		cout << "No hints found. Nothing inserted into database." << endl;
	    }
	    checkConsistency();
	}
	else{
	    cout << filename << " is neither in gff nor in fasta format." << endl;
	}
	ifstrm.close();   
    } catch( string err ){
        cerr << "\n" <<  argv[0] << ": ERROR\n\t" << err << "\n\n";
	exit(1);
    } 
}

void printUsage(){
    cout << "usage:\n\
load2db [parameters] --species=SPECIES --dbaccess=dbname,host,user,passwd  inputfilename\n\
\n\
inputfilename refers to a genome file in FASTA format or a hints file in GFF format\n\
SPECIES is the same identifier as is used in the treefile and alnfile parameters to augustus.\n\
\n\
dbname,host,user,passwd are the name of the SQL database, the host name or IP, the database user and password\n\
When storing genomes/hints of multiple organisms call this program repeatedly for each one.\n\
A single table with the structure\n\
\
is created.\n\
\n\
parameters:\n\
--help        print this usage info\n\
--chunksize   this option is only relevant when loading a sequence file\n\
              the sequences in the input genome are split into chunks of this size so\n\
              that subsequent retrievals of small sequence ranges do not require to read\n\
              the complete - potentially much longer - chromosome. (<= 1000000, default " << chunksize << ")\n\
\n\
example:\n\
     load2db --species=chicken --dbaccess=birds,localhost,mario,dF$n.E chickengenome.fa\n\
     load2db --species=chicken --dbaccess=birds,localhost,mario,dF$n.E chickenhints.gff\n\
\n\
Example code for database creation before calling load2db:\n\
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
void createTableGenomes(){
    mysqlpp::Query query = con.query("CREATE TABLE IF NOT EXISTS genomes (\
        seqid int(10) unsigned NOT NULL AUTO_INCREMENT,			\
        dnaseq longtext NOT NULL,					\
        seqnr int(10) unsigned,					\
        start int(9) unsigned NOT NULL,\
        end int(9) unsigned NOT NULL,\
        speciesid int(10) unsigned,\
        PRIMARY KEY (seqid),\
        KEY region (speciesid,seqnr,start,end),\
        FOREIGN KEY (speciesid,seqnr) REFERENCES seqnames(speciesid,seqnr)\
      ) ENGINE=MyISAM DEFAULT CHARSET=latin1 MAX_ROWS=10000000 AVG_ROW_LENGTH=50000;");
   query.execute();
}
/*
 * create table 'speciesnames' if it does not already exist
 */
void createTableSpeciesnames(){
    mysqlpp::Query query = con.query("CREATE TABLE IF NOT EXISTS speciesnames (\
        speciesid int(10) unsigned NOT NULL AUTO_INCREMENT PRIMARY KEY,\
        speciesname varchar(50) UNIQUE KEY\
      ) ENGINE=MyISAM DEFAULT CHARSET=latin1 MAX_ROWS=1000000 AVG_ROW_LENGTH=50000;");
    query.execute();
}
/*                                                                                                                                                                                                                
 * create table 'seqnames' if it does not already exist                                                                                                                                                        
 */
void createTableSeqnames(){
    mysqlpp::Query query = con.query("CREATE TABLE IF NOT EXISTS seqnames (\
        seqnr int(10) unsigned AUTO_INCREMENT,\
        speciesid int(10) unsigned REFERENCES speciesnames(speciesid),\
        seqname varchar(50),\
        PRIMARY KEY(speciesid,seqnr),\
        UNIQUE KEY(speciesid,seqname)\
      ) ENGINE=MyISAM DEFAULT CHARSET=latin1 MAX_ROWS=1000000 AVG_ROW_LENGTH=50000;");
    query.execute();
}
/*
 * create table 'hints' if it does not already exist
 */
void createTableHints(){
    mysqlpp::Query query = con.query("CREATE TABLE IF NOT EXISTS hints (\
        hintid int(10) unsigned AUTO_INCREMENT PRIMARY KEY,\
        speciesid int(10) unsigned,\
        seqnr int(10) unsigned,\
        source varchar(50),\
        start int(9) unsigned NOT NULL,\
        end int(9) unsigned NOT NULL,\
        score float NOT NULL DEFAULT 0.0,\
        type tinyint unsigned NOT NULL,\
        strand enum('+','-','.') NOT NULL DEFAULT '.',\
        frame enum ('0','1','2','.') NOT NULL DEFAULT '.',\
        priority smallint NOT NULL DEFAULT -1,\
        grp varchar(100) DEFAULT '',\
        mult smallint unsigned DEFAULT 1,\
        esource varchar(10) NOT NULL,\
        KEY region (speciesid,seqnr,start,end),\
        FOREIGN KEY (speciesid,seqnr) REFERENCES seqnames(speciesid,seqnr)\
 ) ENGINE=MyISAM DEFAULT CHARSET=latin1 MAX_ROWS=1000000000 AVG_ROW_LENGTH=50000;");
    query.execute();
}

void createTableFeatureTypes(){
    mysqlpp::Query query = con.query("CREATE TABLE IF NOT EXISTS featuretypes (\
        typeid tinyint unsigned PRIMARY KEY,\
        typename enum('start','stop','ass','dss','tss','tts','exonpart','exon','intronpart','intron',\
        'irpart','CDS','CDSpart','UTR','UTRpart','nonexonpart','genicpart') NOT NULL\
 ) ENGINE=MyISAM DEFAULT CHARSET=latin1 MAX_ROWS=1000000 AVG_ROW_LENGTH=50000;");
    query.execute();
    for(int i=0; i < 17; i++){
	query << "INSERT IGNORE INTO featuretypes VALUES(" << i << "," << (i+1) << ")";
	//cout << "Executing" << endl << query.str() << endl;
	query.execute();
    }
}

/*
 * insert one sequence into database
 * return the number of chunks that sequence was cut into
 * start and end are 0-based and inclusive.
 */
int insertSeq(string sequence, char *name, int length, int speciesid){
    int chunks = 0;
    int start = 0, end;
    int seqnr = insertSeqName(name,speciesid);
    while (start < length) {
	end = (start + chunksize < length)? (start + chunksize - 1) : length - 1;
	mysqlpp::Query query = con.query();
	query << "INSERT INTO genomes (dnaseq,seqnr,speciesid,start,end) VALUES(\""
	      << sequence.substr(start, end-start+1) << "\"," << seqnr << "," 
	      << speciesid << "," << start << "," << end << ")";
	//cout << "Executing" << endl << query.str() << endl;
	query.execute();
	chunks++;
	start += chunksize;
    }
    return chunks;
}
/*
 * insert one hint into database
 * start and end are 0-based and inclusive
 */
void insertHint(Feature& f, string species){
    mysqlpp::Query query = con.query();
    int speciesid = getSpeciesID(species);
    int seqnr = getSeqNr(&(f.seqname[0]),speciesid);
    string attributes = "";
    string values = "";
    switch (f.strand) {
        case plusstrand: attributes+="strand,"; values+=("\"+\","); break;
        case minusstrand: attributes+="strand,"; values+=("\"-\","); break;
        default : break;
    }
    if (f.frame != -1){
	attributes+="frame,";
	values+=("\"" + itoa(f.frame) + "\",");
    }
    if (f.groupname != ""){
	attributes+="grp,";
	values+=("\"" + f.groupname + "\",");
    }
    if (f.mult > 1){
	attributes+="mult,";
	values+=(itoa(f.mult) + ",");
    }
    if (f.priority >= 0){
	attributes+="priority,";
	values+=(itoa(f.priority) + ",");
    }
    query << "INSERT INTO hints (speciesid,seqnr,source,start,end,score,type," << attributes << "esource) VALUES ("
	  << speciesid << "," << seqnr <<  ",\"" << f.source << "\"," << f.start << "," << f.end << "," << f.score << ","
          << f.type << "," << values << "\""<< f.esource << "\")";
    //cout << "Executing" << endl << query.str() << endl;
    query.execute();
  
}
/*
 * returns the speciesid for a given species name
 * if the species name is not in the database it is inserted
 */
int getSpeciesID(string species){
    mysqlpp::Query query = con.query();
    query << "SELECT speciesid FROM speciesnames WHERE speciesname='" <<species<<"'";;
    //cout << "Executing" << endl << query.str() << endl;
    mysqlpp::StoreQueryResult res = query.store();
    if(!res.empty()){
	return res[0]["speciesid"];
    }
    else{
	// insert species name
	query << "INSERT INTO speciesnames (speciesname) VALUES (\"" << species << "\")";
	//cout << "Executing" << endl << query.str() << endl;
	query.execute();
	return query.insert_id();
    }
}
/*
 * returns the seqnr for a given sequence name and species id
 * it the sequence name is not in the database, it will be inserted
 */
int getSeqNr(char* seqname,int speciesid){
    mysqlpp::Query query = con.query();
    query << "SELECT seqnr FROM seqnames WHERE seqname='" << seqname << "' AND speciesid=" << speciesid;
    // cout << "Executing" << endl << query.str() << endl;
    mysqlpp::StoreQueryResult res = query.store();
    if (!res.empty()){
        return res[0]["seqnr"];
    } else{
        query << "INSERT INTO seqnames (speciesid,seqname) VALUES (" << speciesid << ",\"" << seqname << "\")";
	// cout << "Executing" << endl << query.str() << endl;
        query.execute();
        return query.insert_id();
    }
}
int insertSeqName(char* seqname,int speciesid){
    try{
	mysqlpp::Query query = con.query();
	query << "INSERT INTO seqnames (speciesid,seqname) VALUES (" << speciesid << ",\"" << seqname << "\")";
	//cout << "Executing" << endl << query.str() << endl;
	query.execute();
	return query.insert_id();
    } catch(const mysqlpp::BadQuery& e){
        cerr <<"sequence "<< seqname <<" already exists in database" << endl;
	cerr <<"Please delete sequence, before reloading"<<endl;
	exit(1);
    }
}

void checkConsistency(){
    mysqlpp::Query query = con.query("SELECT DISTINCT speciesname,seqname FROM hints as h,\
        speciesnames as s,seqnames as n where h.speciesid=s.speciesid AND\
        s.speciesid=n.speciesid AND h.seqnr=n.seqnr AND NOT EXISTS\
        (SELECT DISTINCT speciesid,seqnr FROM genomes as g WHERE g.speciesid=h.speciesid AND g.seqnr=h.seqnr);");
    mysqlpp::StoreQueryResult res = query.store();
    if(!res.empty()){
	cout << "WARNING: the following sequences are not in the database although we have hints for them" << endl;
	for(size_t i=0;i<res.num_rows();i++){
	    cerr<<res[i]["seqname"]<<" (species "<<res[i]["speciesname"]<<")"<<endl;
	}
    }
}
void removeDuplicates(){
    cout << "removing duplicate hints"<< endl;
    mysqlpp::Query query = con.query(" ALTER IGNORE TABLE hints ADD UNIQUE\
        (speciesid,seqnr,source,start,end,score,type,strand,frame, priority,grp,mult,esource);");
    mysqlpp::SimpleResult res = query.execute();
    cout << res.info() << endl;
}

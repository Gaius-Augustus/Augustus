/**********************************************************************
 * retrieve genomic sequences from a mysql database, 
 * file:    getSeq.cc
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 *
 * authors: Stefanie Koenig
 *
 * date    |   author      |  changes
 * --------|---------------|------------------------------------------
 * 20.11.13|Stefanie Koenig| creation of the file
 **********************************************************************/

// Project includes
#include "randseqaccess.hh"

// standard C/C++ includes
#include <getopt.h>     /* for getopt_long; standard getopt is in unistd.h */
#include <stdlib.h>     /* for exit() */

using namespace std;

int fold=60; // line width of output sequence

void printUsage();
void printSeq(string sequence,int length);
/*
 * main
 */
int main( int argc, char* argv[] ){

    int intMax = std::numeric_limits<int>::max();
   
    int c;
    string species, seqname;
    int start = 1;
    int end = intMax;
    Strand strand=plusstrand;
    int help=0;

    static struct option long_options[] = {
	{"dbaccess",1, 0, 'c'},
        {"species",1, 0, 's'},
	{"speciesfilenames", 1, 0,'f'},
        {"seq",1, 0, 'q'},
        {"start",1, 0, 'a'},
        {"end",1, 0, 'b'},
	{"rc",0,0,'r'},
	{"help",0,0,'h'},
        {0,0,0,0}
    };
    int option_index = 0;
    while ((c = getopt_long(argc, argv, "c:s:f:q:a:b:rh", long_options, &option_index)) != -1) {
        switch(c)
            {
	    case 'c':
		Constant::dbaccess = optarg;
                break;
            case 's':
                species = optarg;
                break;
	    case 'f':
		Constant::speciesfilenames = optarg;
		break;
            case 'q':
                seqname = optarg;
                break;
            case 'a':
                start = atoi(optarg);
                break;
            case 'b':
                end = atoi(optarg);
                break;
            case 'r':
                strand=minusstrand;
                break;    
	    case 'h':
		help=1;
		break;
	    }
    }

    if (help){
	printUsage();
	exit(1);
    }
    if (species.empty()){
	cerr << "Missing species name. Required parameter." << endl;
	printUsage();
	exit(1);
    }
    if (seqname.empty()){
	cerr << "Missing sequence name. Required parameter." << endl;
	printUsage();
	exit(1);
    }
    
    if (Constant::dbaccess.empty()){
	cerr << "Missing database access info. dbaccess is a required parameter." << endl;
	printUsage();
	exit(1);
    }

    if(end < start || start < 1 || end < 1){
	cerr << "Not a genomic interval. Typo in the start or end coordinate?" << endl;
	exit(1);
    }
    
    DbSeqAccess *rsa=NULL;

    if (Constant::dbaccess.find(',') != string::npos){ // assuming mysql access
	cerr << "assuming a MySQL database" << endl;
#ifdef AMYSQL
	rsa = new MysqlAccess;
#else
	cerr << "Database access not possible with this compiled version. Please recompile with flag MYSQL." << endl;
	exit(1);
#endif
	
    }
    else if(Constant::dbaccess.find('.') != string::npos){ // assuming sqlite access
	cerr << "assuming an SQLite database" << endl;
	if(Constant::speciesfilenames.empty()){
	    cerr << "Missing parameter speciesfilenames." << endl;
	    exit(1);
	}
#ifdef SQLITE
	rsa = new SQLiteAccess(Constant::dbaccess.c_str());
#else
	cerr <<"Database access not possible with this compiled version. Please recompile with flag SQLITE." << endl;
	exit(1);
#endif
    }

    try{
	AnnoSequence *annoseq = rsa->getSeq(species, seqname, start-1, end-1, strand);
	if(annoseq){
	    if(annoseq->offset + annoseq->length < end ){
		if(end < intMax){
		    cerr <<"Warning: end position "<< end <<" is past the end of the sequence." << endl;
		}
		end = annoseq->offset + annoseq->length;
		cerr <<"Retrieving "<< seqname << ":" << start << "-" << end << endl;
	    }
	    // print fasta header
	    cout << ">"<< seqname << " " << start << " " << end;
	    if(strand == plusstrand)
		cout << " +" <<endl;
	    else
		cout << " -" <<endl;
	    // print sequence
	    printSeq(annoseq->sequence, end-start+1);
	}
    }
    catch(ProjectError e){
	cerr << "random sequence access failed on " << species << ", " << seqname << ":" 
	     << start << "-" << end << endl;;
	cerr << e.getMessage() << endl;
	exit(1);
    }
    exit(1);
}

void printUsage(){
    cerr << "usage:\n\
getSeq [parameters] --species=SPECIES --seq=SEQUENCE --dbaccess=dbname,host,user,passwd \n\
\n\
SPECIES is the species identifier used when loading the sequence into the database\n\
SEQUENCE is the ID of the sequence to retrieve\n\
dbname,host,user,passwd are the name of the SQL database, the host name or IP, the database user and password\n\
\n\
parameters:\n\
--help        print this usage info\n\
--rc          output the reverse complement of the sequence\n\
--start=N     retrieve subsequence starting at position N (coordinates are 1-based)\n\
--end=N       retrieve subsequence ending at position N (coordinates are 1-based)\n\
\n\
example:\n\
     getSeq --species=hg19 --seq=chr21 --dbaccess=saeuger,localhost,cgp,AVglssd8 \n\
     getSeq --species=hg19 --seq=chr21 --start=47870612  --end=48086047 --rc --dbaccess=saeuger,localhost,cgp,AVglssd8 \n\
\n";
}

void printSeq(string sequence,int length){
    int start=0, end;
    while(start < length){
	end = (start + fold < length)? (start + fold - 1) : length - 1;
	cout << sequence.substr(start,end-start+1) << endl;
	start+=fold;
    }
}

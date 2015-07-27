/**********************************************************************
 * file:    main.cc
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 *
 * authors: Stefanie Koenig
 *
 * date    |   author      |  changes
 * --------|---------------|------------------------------------------
 * 22.01.15|Stefanie Koenig| creation of the file
 **********************************************************************/

// standard C/C++ includes
#include <getopt.h>     /* for getopt_long; standard getopt is in unistd.h */
#include <iostream>
#include <thread>
#include <functional>   /* for std::ref */

// project includes
#include "genome.hh"
#include "projectio.hh"


using namespace std;


void printUsage();
map<string, pair<string,string> > getFileNames(string listfile);

/*
 * main
 */
int main( int argc, char* argv[] ){

    int c;
    int help=0;
    string srcGenome;
    string halfile;
    string tmpdir;
    string outdir;
    string gtfs;
    string halParam = "";
    string halLiftover_exec = "";
    string homGeneFile = ""; 
    size_t maxCpus=1;

    static struct option long_options[] = {
	{"gtfs", 1, 0,'g'},
        {"srcGenome",1, 0, 's'},
	{"halfile",1, 0, 'a'},
 	{"tmpdir",1, 0, 't'},
 	{"outdir",1, 0, 'o'},
 	{"halLiftover_exec_dir",1, 0, 'e'},
	{"cpus",1, 0, 'n'},
	{"noDupes",0,0,'d'},
	{"printHomologs",1,0,'m'},
	{"help",0,0,'h'},
        {0,0,0,0}
    };
    int option_index = 0;
    while ((c = getopt_long(argc, argv, "g:s:a:t:o:e:n:dm:h", long_options, &option_index)) != -1) {
        switch(c)
            {
	    case 'g':
		gtfs = optarg;
                break;
            case 's':
                srcGenome = optarg;
                break;
	    case 'a':
		halfile = optarg;
		break;   
	    case 't':
		tmpdir = optarg;
		break;   
	    case 'o':
		outdir = optarg;
		break;   
	    case 'e':
		 halLiftover_exec = optarg;
		break;   
	    case 'n':
		maxCpus = atoi(optarg);
		break;   
	    case 'd':
		halParam+="--noDupes ";
		break;
	    case 'm':
		homGeneFile = optarg;
		break;   
	    case 'h':
		help=1;
		break;
	    }
    }
    
    try{
	if (help){
	    printUsage();
	    exit(1);
	}
	if (halfile.empty()){
	    cerr << "Missing input hal file. Required parameter." << endl;
	    printUsage();
	    exit(1);
	}
	if (tmpdir.empty())
	    tmpdir="tmp";
	tmpdir = expandDir(tmpdir);
	createDir(tmpdir);
	
	halLiftover_exec = expandDir(halLiftover_exec);
	halLiftover_exec +="halLiftover";
	
	if(!outdir.empty()){
	    outdir = expandDir(outdir);
	    createDir(outdir);
	}
	
	if (gtfs.empty()){
	    cerr << "Missing input gtf file names. Required parameter." << endl;
	    printUsage();
	    exit(1);
	}
	
	if(maxCpus < 1){
	    maxCpus = 1;
	    cerr << "number of cpus must be at least 1. Proceeding with --cpus=1" << endl; 
	}
	// check if the boost library is installed (only required for option --printHomologs)
	if(!homGeneFile.empty()){
#ifndef BOOST
	    throw ProjectError("The option --printHomologs requires the boost C++ libary.\n"
                               "Please install the boost library, e.g. using the APT package manager\n\n"
                               "sudo apt-get install libboost-graph-dev\n\n"
                               "Then edit the Makefile by setting the flag BOOST = true and recompile homGeneMapping.\n");
#endif
	}
	/*
	 * check if external program halLiftOver is executable
	 */
	string cmd = "which " + halLiftover_exec;
	string path = exec(cmd.c_str());
	if(path.empty())
	    throw ProjectError("halLiftover is not executable.\n" 
			      "Please add the directory which contains the executable halLiftover to the\n" 
			       "PATH environment variable or specify the path with --halLiftover_exec_dir.\n");
	
	/*
	 * parsing of input gene files
	 */    
	map<string, pair<string, string> > filenames = getFileNames (gtfs);
	Genome::setNumGenomes(filenames.size());
	
	vector<Genome> genomes;
	for(map<string, pair<string, string> >::iterator it = filenames.begin(); it != filenames.end(); it++){
	    Genome genome(it->first, genomes.size());
	    genome.parseGTF(it->second.first);
	    if(!it->second.second.empty()) // read hints file if specified
		genome.parseExtrinsicGFF(it->second.second);
	    genome.setTmpDir(tmpdir);
	    genome.printBed(); // print sequence coordinates, that need to be mapped to the other genomes, to file
	    genomes.push_back(genome);
	}
	
	/*
	 * haLiftover from each genome to each other genome (quadratic to the number of genomes)
	 * TODO:
	 * check if liftover is symmetric (i.e. i -> j is the same as j -> i )
	 * in this case, runtime can be reduced by running the second loop only for j > i 
	 */
	cout << "halLiftover starts. Processing" << endl;
        if(maxCpus > 1){ // in parallel
	    vector<thread> th;
	    for(int i = 0; i < genomes.size(); i++){
		cout << genomes[i].getName() << endl;
		for(int j = 0; j < genomes.size(); j++){
		    if(i != j){
			if(th.size() == maxCpus){ // wait for all running threads to finish
			    // synchronize threads
			    for(int i = 0; i < th.size(); i++){
				th[i].join();
			    }
			    th.clear();
			}
			// halLiftover
			th.push_back(thread(&Genome::liftOverTo, &genomes[i], ref(genomes[j]), halfile, halLiftover_exec, halParam));
		    }
		}
	    }
	    // synchronize threads
	    for(int i = 0; i < th.size(); i++){
		th[i].join();
	    }
	}
	else{ // sequentially
            for(int i = 0; i < genomes.size(); i++){
                cout << genomes[i].getName() << endl;           
                for(int j = 0; j < genomes.size(); j++){
                    if(i != j){
                        // halLiftover
                        genomes[i].liftOverTo(genomes[j], halfile, halLiftover_exec, halParam);
                    }
                }
            }
        }
	for(int i = 0; i < genomes.size(); i++){
	    for(int j = 0; j < genomes.size(); j++){
		if(i != j){
		    genomes[i].readBed(genomes[j]); // reading in bed files with mapped coordinates
		}
	    }
	    // identify homologous gene features (exons/introns)
	    genomes[i].mapGeneFeatures(genomes);
	    // print extended gene files with homology information
	    genomes[i].printGFF(outdir,genomes);
	}
	// print a list with homologous transcript IDs, e.g.
	// # 0     dana
	// # 1     dere
	// # 2     dgri
	// # 3     dmel
	// # 4     dmoj
	// # 5     dper
	// (0, jg4139.t1), (0, jg4140.t1), (1, jg7797.t1), (2, jg3247.t1), (4, jg6720.t1), (5, jg313.t1)
	// (1, jg14269.t1), (3, jg89.t1) (5, jg290.t1)
	// ...
#ifdef BOOST
	if(!homGeneFile.empty()){
	    printHomGeneList(homGeneFile,genomes);
	}
#endif
	for(int i = 0; i < genomes.size(); i++){
	    genomes[i].destroyGeneList();
	    genomes[i].destroyHintList();
	}

    } catch( ProjectError& err ){
	cerr << "\n" <<  argv[0] << ": ERROR\n" << err.getMessage( ) << "\n\n";
	exit(1);
    }
    return 0;
}

void printUsage(){
    cerr << "homGeneMapping takes a set of gene predictions of different genomes and a hal\n\
alignment of the genomes and prints a summary for each gene, e.g.\n\
- how many of its exons/introns are in agreement with genes of other genomes\n\
- how many of its exons/introns are supported by extrinsic evidence from any of the genomes\n\
- a list of geneids of homologous genes\n\n\
usage:\n\
homGeneMapping [Options] --gtfs=gffilenames.tbl --halfile=aln.hal\n\
\n\
ARGUMENTS:\n\
--halfile=aln.hal             input hal file\n\
--gtfs=gtffilenames.tbl       a text file containing the locations of the input gene files\n\
                              and optionally the hints files (both in GTF format).\n\
                              The file is formatted as follows:\n\n\
                              name_of_genome_1  path/to/genefile/of/genome_1  path/to/hintsfile/of/genome_1\n\
                              name_of_genome_2  path/to/genefile/of/genome_2  path/to/hintsfile/of/genome_2\n\
                              ...\n\
                              name_of_genome_N  path/to/genefile/of/genome_N  path/to/hintsfile/of/genome_N\n\
\n\
OPTIONS:\n\
--help                        print this usage info\n\
--cpus=N                      N is the number of CPUs to use (default: 1)\n\
--noDupes                     do not map between duplications in hal graph. (default: off)\n\
--halLiftover_exec_dir=DIR    Directory that contains the executable halLiftover\n\
                              If not specified it must be in $PATH environment variable.\n\
--tmpdir=DIR                  a temporary file directory that stores lifted over files. (default 'tmp/' in current directory)\n\
--outdir=DIR                  file direcory that stores output gene files. (default: current directory)\n\
--printHomologs=FILE          prints disjunct sets of homologous transcripts to FILE, e.g.\n\
                              # 0     dana\n\
                              # 1     dere\n\
                              # 2     dgri\n\
                              # 3     dmel\n\
                              # 4     dmoj\n\
                              # 5     dper\n\
                              (0, jg4139.t1) (0, jg4140.t1) (1, jg7797.t1) (2, jg3247.t1) (4, jg6720.t1) (5, jg313.t1)\n\
                              (1, jg14269.t1) (3, jg89.t1) (5, jg290.t1)\n\
                              ...\n\
                              Two transcripts are in the same set, if all their exons/introns are homologs and their are\n\
                              no additional exons/introns.\n\
                              This option requires the Boost C++ Library\n\
\n\
example:\n\
homGeneMapping --noDupes --halLiftover_exec_dir=~/tools/progressiveCactus/submodules/hal/bin --gtfs=gtffilenames.tbl --halfile=msca.hal\n\
\n";
}

/*
 * read an input file of format:
 * human <TAB> human.genepred.gtf <TAB> human.hints.gff
 * mouse <TAB> mouse.genepred.gtf <TAB> mouse.hints.gff
 * cow   <TAB> cow.genepred.gtf   <TAB> cow.hints.gff
 * ...
 */
map<string,pair<string,string> > getFileNames (string listfile){
    map<string,pair<string, string> > filenames;
    ifstream ifstrm(listfile.c_str());
    if (ifstrm.is_open()){
        char buf[256];
	while(ifstrm.getline(buf,255)){
	    stringstream stm(buf);
            string species, genefile;
            if(stm >> species >> genefile){
		// expand home
		genefile=expandHome(genefile);
		string hintsfile;
		if(stm >> hintsfile)
		    hintsfile=expandHome(hintsfile);
		filenames[species] = make_pair(genefile,hintsfile);
	    }
	    else
		throw ProjectError(listfile + " has wrong format in line\n" + buf + "\nCorrect format:\n\n"
				   "name_of_genome_1  path/to/gtf/of/genome_1  path/to/hintsfile/of/genome_1\n"
				   "name_of_genome_2  path/to/gtf/of/genome_2  path/to/hintsfile/of/genome_2\n...\n"
				   "name_of_genome_N  path/to/gtf/of/genome_N  path/to/hintsfile/of/genome_N\n\n"
				   "the last column is optional.\n");
	}
        ifstrm.close();
    }
    else
	throw ProjectError("Could not open input file " + listfile + ".\n");
    return filenames;
}

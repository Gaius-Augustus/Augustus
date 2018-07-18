/**
 * \file main.cpp
 */

#include "Genomic_Data.hpp"
#include "Compute_UTRs.hpp"
#include "Supporting_Methods.hpp"
#include "Splice_Sites.hpp"
#include "UTRs.hpp"
#include "Test.hpp"
#include "Global.hpp"

#include <string>
#include <cstring>
#include <vector>

#include <iostream>
#include <sstream>
#include <iostream>
#include <getopt.h>

#include <boost/assign/std/vector.hpp>

using namespace std;
using namespace boost::assign;

stringstream g_output;

#ifndef TEST_MODE

int main(int argc, char* argv[]) {
	//argv[0] is the path and name of the program itself

	//required options
	string scaffold_fname;
	string coding_region_fname;
	string introns_fname;
	string repeat_fname;
	string wiggle_fname;

	bool zero_cov = false;

	int opt;
	int option_index = -1;

	static option long_options[] = {
		{"in-scaffold-file", 1, 0, 'S'},
		{"in-coding-region-file", 1, 0, 'C'},
		{"in-intron-file", 1, 0, 'I'},
		{"in-repeat-file", 1, 0, 'R'},
		{"in-wiggle-file", 1, 0, 'W'},
		{"out-file-name", 1, 0, 'o'},
		{"smoothing-window-size", 1, 0, 'w'},
		{"read-length", 1, 0, 'r'},
		{"limit", 1, 0, 'l'},
		{"drop-window-size", 1, 0, 'v'},
		{"minimum-length", 1, 0, 'n'},
		{"minimum-average-coverage", 1, 0, 'c'},
		{"percent-window", 1, 0, 'p'},
		{"percent-intron", 1, 0, 'i'},
		{"percent-multiplicity", 1, 0, 'm'},
		{"splice-sites", 1, 0, 's'},
		{"zero-coverage", 0, 0, 'z'},
		{"format-examples", 0, 0, 'f'},
		{"help", 0, 0, 'h'},
		{0, 0, 0, 0}
	};

	//run getopt to get the parameter for ZeroCov
	while ( (opt = getopt_long(argc, argv, "S:C:I:R:W:o:w:r:l:v:n:c:p:i:m:s:zfh", long_options, &option_index)) != -1) {
		if (opt == 'z') {
			zero_cov = true;
		}
	}

	//optional options and default values
	string output_fname = "utrs.gff";
	unsigned read_length = 150;
	unsigned smoothing_window_size = !zero_cov ? 150 : 1;
	unsigned limit = 5000;
	unsigned drop_window_size = !zero_cov ? 50 : 10;
	unsigned min_length = 2;
	unsigned min_average_cov = !zero_cov ? 10 : 1;
	double p_win = !zero_cov ? 0.6 : 0.1;
	double p_int = 0.5;
	double p_mult = 0.1;
	vector< pair <string, string> > splice_sites;
	splice_sites += make_pair("gt","ag");
	bool help = false;
	bool format_examples = false;

	string raw_splice_sites; //only shortly used for parsing of splice sites

	reset_getopt(); //reset getopt, to call getopt() a second time

	while ( (opt = getopt_long(argc, argv, "S:C:I:R:W:o:w:r:l:v:n:c:p:i:m:s:zfh", long_options, &option_index)) != -1) {
		switch(opt) {
			case 'S':
				scaffold_fname = optarg;
				break;
			case 'C':
				coding_region_fname = optarg;
				break;
			case 'I':
				introns_fname = optarg;
				break;
			case 'R':
				repeat_fname = optarg;
				break;
			case 'W':
				wiggle_fname = optarg;
				break;
			case 'o':
				output_fname = optarg;
				break;
			case 'w':
				smoothing_window_size = atoi(optarg);
				break;
			case 'r':
				read_length = atoi(optarg);
				break;
			case 'l':
				limit = atoi(optarg);
				break;
			case 'v':
				drop_window_size = atoi(optarg);
				break;
			case 'n':
				min_length = atoi(optarg);
				break;
			case 'c':
				min_average_cov = atoi(optarg);
				break;
			case 'p':
				p_win = atof(optarg);
				break;
			case 'i':
				p_int = atof(optarg);
				break;
			case 'm':
				p_mult = atof(optarg);
				break;
			case 's':
				raw_splice_sites = optarg;
				//parse splice sites
				splice_sites = parse_sp_sites(raw_splice_sites);
				break;
			case 'z':
				zero_cov = true;
				break;
			case 'f':
		    	format_examples = true;
		    	break;
			case 'h':
				help = true;
				break;
		}
	}

	bool no_options = (argc == 1) ? true : false;
	bool required_option_missing = false;

	//testing if obligatory options are set
	if ( scaffold_fname.empty() && !no_options && !help && !format_examples  ) {
		ERR_STRM << "Error (in main.cpp): option --in-scaffold-file/-S is an obligatory argument and was not set." << endl;
		cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;
		required_option_missing = true;
	}

	if ( coding_region_fname.empty() && !no_options && !help && !format_examples ) {
		ERR_STRM << "Error (in main.cpp): option --in-coding-region-file/-C is an obligatory argument and was not set." << endl;
		cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;
		required_option_missing = true;
	}

	if ( introns_fname.empty() && !no_options && !help && !format_examples ) {
		ERR_STRM << "Error (in main.cpp): option --in-intron--file/-I is an obligatory argument and was not set." << endl;
		cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;
		required_option_missing = true;
	}

	if ( wiggle_fname.empty() && !no_options && !help && !format_examples ) {
		ERR_STRM << "Error (in main.cpp): option --in-wiggle-file/-W is an obligatory argument and was not set." << endl;
		cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;
		required_option_missing = true;
	}

	//help printed if no option was set or a required option was missing or the option for help was set
	if ( ( help || required_option_missing || no_options )) {
		cout << endl;
		cout << "utrrnaseq - A tool for identifying UTRs of annotated genes on the basis of RNA-Seq data." << endl;
		cout << "            Starting from known start/stop codons of genes, coverage with RNA-Seq data in" << endl;
		cout << "            potential UTRs is monitored for a drastic drop, and such points are defined as" << endl;
		cout << "	     ends of UTRs. Spliced alignments are incorported as introns into UTRs. Drops" << endl;
		cout << "            in coverage due to repeat masking are not reported as UTR endpoints. Only UTR" << endl;
		cout << "            models with evidence from RNA-Seq coverage are reported." << endl;
		cout << endl;
		cout << "Obligatory arguments are" << endl;
		cout << endl;
		cout << "--long option  			-short option   data type       description" << endl;
		cout << "--in-scaffold-file		-G		s      		scaffold file in FASTA-format. Required." << endl;
		cout << "--in-coding-region-file		-C		s      		file with start and stop features in GTF/GFF format. Required." << endl;
		cout << "--in-intron-file		-I		s      		intron file in GTF/GFF format. Required." << endl;
		cout << "--in-wiggle-file		-W		s      		wiggle file in WIG-format. Required." << endl;
		cout << endl;
		cout << "Optional arguments are:" << endl;
		cout << endl;
		cout << "--long option			-short option	data type	description" << endl;
		cout << "--in-repeat-file		-R		s      		repeat file in GTF/GFF format." << endl;
		cout << "--out-file-name	 		-o		s		output filename. Optional. Default Value: 'UTRs.gff'" << endl;
		cout << "--smoothing-window-size		-w		i		smoothing window size. Optional. Default Value: 150" << endl;
		cout << "--read-length	 		-r		i		read length of RNA-Seq data. Optional. Default Value: 150" << endl;
		cout << "--limit	 			-l		i		maximal distance from computation start. Optional. Default Value: 5000" << endl;
		cout << "--drop-window-size		-v		i		window size after UTR end. Optional. Default Value: 50" << endl;
		cout << "--minimum-length		-n		i		minimal UTR length. Optional. Default Value: 2" << endl;
		cout << "--minimum-average-coverage	-c		i		minimal average UTR coverage. Optional. Default Value: 10" << endl;
		cout << "--percent-window		-p		d		percentage of window coverage after UTR. Optional. Default Value: 0.6" << endl;
		cout << "--percent-intron		-i		d		percentage of coverage in introns. Optional. Default Value: 0.5" << endl;
		cout << "--percent-multiplicity		-m		d		percentage of multiplicity of introns. Optional. Default Value: 0.1" << endl;
		cout << "--splices-sites			-s		s 		accepted splice sites. If 'all' is chosen, no splice site filtering is done."
			 << "Optional. Default Value: GT_AG" << endl;
		cout << "--zero-coverage	 		-z		none		Determination of UTRs based on zero coverage. Optional. Default Value: false" << endl;
		cout << "--format-examples		-f		none		Only print format examples of input files. Optional. Default Value: false" << endl;
		cout << "--help	 			-h		none		Produce help message. Optional. Default Value: false" << endl;

		if( required_option_missing && !no_options ){
		  ERR_STRM << "Error (in main.cpp): One or more required option missing!" << endl;
		  cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << endl;
		  abort();
		}

		if (!format_examples) {
			return 0;
		}
	}

	//format examples printed if option for format exampled is set
	if( format_examples ) {
	  cout << endl;
	  cout << "Format examples for input files" << endl;
	  cout << endl;
	  cout << "--in-scaffold-file/-S contains the scaffold sequence of the target species. It must be in (multiple) FASTA format, e.g.:" << endl;
	  cout << endl;
	  cout << ">chr19" << endl;
	  cout << "CCGTCCTGGTGCCCAACATGGAGAAGTACCTGACTTCTGATGTGGCGAAGAAGCTCAACG" << endl;
	  cout << "GCGCGAGGAGCGAGGAGTCGGGCTCCGAGGAGAAGGTTGAGAAGGTGTAATGGAATGTTG" << endl;
	  cout << "GAAGCGATGTCGGTCTTGTTTGCTCATTTGTTCGTCTGATTGATGGTTTCCAACTAGCTG" << endl;
	  cout << ">chr1" << endl;
	  cout << "CCGTCCTGGTGCCCAACATGGAGAAGTACCTGACTTCTGATGTGGCGAAGAAGCTCAACG" << endl;
	  cout << "GCGCGAGGAGCGAGGAGTCGGGCTCCGAGGAGAAGGTTGAGAAGGTGTAATGGAATGTTG" << endl;
	  cout << endl;
	  cout << "--in-coding-region-file/-C contains the start and end positions of annotated genes/CDS (without UTR). It must be in tabular-separated 9-column gff format. ";
	  cout << "It must contain the features start_codon and stop_codon, and the last column must contain gene and transcript IDs, e.g.:" << endl;
	  cout << endl;
	  cout << "chr19\thg19_refGene\tstart_codon\t110679\t110681\t0\t+\t.\tgene_id \"NM_001005240\"; transcript_id \"NM_001005240\";" << endl;
	  cout << "chr19\thg19_refGene\tstop_codon\t111594\t111596\t0\t+\t.\tgene_id \"NM_001005240\"; transcript_id \"NM_001005240\";" << endl;
	  cout << endl;
	  cout << "--in-intron-file/-I gives information about possible introns that was derived from aligning RNA-Seq data against the genome. ";
	  cout << "It must be in tabular-separated 9-column gff format. ";

	  cout << "It must contain the feature intron. The last column must contain a tag \"mult=\" or \"grp=\" or \"group=\" ";
	  cout << "where the integer gives the coverage of this particular intron by RNA-Seq data and the four character ";
	  cout << "string defines the splice site pattern. Currently, only the splice site pattern GTAG is supported. E.g.: " << endl;
	  cout << endl;

	  cout << "chr19\tb2h_lib1\tintron\t66500\t70927\t0\t-\t.\tmult=6_GTAG" << endl;
	  cout << "chr19\tb2h_lib1\tintron\t70815\t72170\t0\t+\t.\tmult=1_GTAG" << endl;
	  cout << endl;
	  cout << "--in-repeat-file/-R contains sequence segments that belong to repetitive genomic regions. ";
	  cout << "It must be in tabular-separated 9-column gff format. ";
	  cout << "It must contain the feature nonexonpart.:" << endl;
	  cout << "chr19\trepmask\tnonexonpart\t10001189\t10001346\t0\t.\t.\tsrc=RM" << endl;
	  cout << "chr19\trepmask\tnonexonpart\t10001347\t10001650\t0\t.\t.\tsrc=RM" << endl;
	  cout << endl;
	  cout << "--in-wiggle-file/-W contains RNA-Seq coverage information in variable step or fixed step wiggle format, e.g.:" << endl;
	  cout << "track name=trimmed.sf.psl type=wiggle_0" << endl;
	  cout << "variableStep chrom=chr19" << endl;
	  cout << "60001 200" << endl;
	  cout << "60004 150" << endl;
	  cout << "60005 152" << endl;     
	  cout << endl;

	  return 0;
	}

	UTRs::parameters().output_fname = output_fname;
	UTRs::parameters().window_size = smoothing_window_size;
	UTRs::parameters().read_length = read_length;
	UTRs::parameters().limit = limit;
	UTRs::parameters().drop_window_size = drop_window_size;
	UTRs::parameters().min_length = min_length;
	UTRs::parameters().min_average_cov = min_average_cov;
	UTRs::parameters().p_win = p_win;
	UTRs::parameters().p_int = p_int;
	UTRs::parameters().p_mult = p_mult;
	UTRs::parameters().zero_cov = zero_cov; //for optional computation


	char* em_mem = new char[16384];   // Reserve emergency memory that can be deleted just in case we run out of memory
	try {
		if (repeat_fname.empty())
			Genomic_Data::initialize(scaffold_fname, coding_region_fname, introns_fname, "", splice_sites, false);
		else
			Genomic_Data::initialize(scaffold_fname, coding_region_fname, introns_fname, repeat_fname, splice_sites, true);

		compute_UTRs(wiggle_fname);
		cout << "Computation of UTRs finished!" << endl;
		UTRs::output(); //output printed
	}
	catch (bad_alloc& ex) {
		delete[] em_mem;
        cerr << "utrrnaseq ran out of memory";
        exit(1);
	}
	catch (exception& ex) {
		delete[] em_mem;
        cerr << "utrrnaseq crashed for unknown reasons";
        exit(1);
	}

	cout << "Finished!" << endl;
	return 0;
}

#endif


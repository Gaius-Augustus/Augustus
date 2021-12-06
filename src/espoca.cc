/*
 * espoca.cc
 *
 * License: Artistic License, see file LICENSE.TXT or 
 *          https://opensource.org/licenses/artistic-license-1.0
 * 
 */

// project includes
#include "types.hh"
#include "gene.hh"
#include "genbank.hh"
#include "namgene.hh"
#include "evaluation.hh"
#include "statemodel.hh"
#include "codonMSA.hh"
#ifdef COMPGENEPRED
#include "compgenepred.hh"
#endif
// standard C/C++ includes
#include <fstream>
#include <sys/stat.h>
#include <getopt.h>     /* for getopt_long; standard getopt is in unistd.h */
#include <stdlib.h>     /* for exit() */


/*
 * Set the global variables according to command line options or configuration file
 */

void printUsage();
void printIntro();

/*
 * main
 */
int main( int argc, char* argv[] ){
    string     configfile;
    string     commandline;
    int        errorcode = 0;
    string     species;
    int        help = 0;
    double     branchlength = 0.2;

    LLDouble::setOutputPrecision(3);


    for (int i=0; i<argc; i++){
	commandline += argv[i];
	if (i<argc-1) 
	    commandline += " ";
    }

    static struct option long_options[] = {
      {"species",1, 0, 's'},
      {"alnfile",1, 0, 'a'},
      {"treefile",1, 0, 't'},
      {"help",0,0,'h'},
      {"useAminoAcidRates",1,0,'r'},
      {"branchlength",1,0,'b'},
      {0,0,0,0}
    };

    int option_index = 0;
    int c; 
    while ((c = getopt_long(argc, argv, "s:a:t:rb:h", long_options, &option_index)) != -1) {
      switch(c)
	{
	case 's':
	  species = optarg;
	  break;
	case 'a':
	  Constant::codonalnfile = optarg;
	  break;
	case 't':
	  Constant::treefile = optarg;
	  break;
	case 'r':
	  Constant::useAArates = true;
	  break;
	case 'b':
	  branchlength = atof(optarg);
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
      cerr << "Missing species name." << endl;
      printUsage();
      exit(1);
    }
    if (Constant::codonalnfile.empty()){
      cerr << "Missing codon alignment file." << endl;
      printUsage();
      exit(1);
    }
    if(Constant::treefile.empty()){
      cerr << "Warning: No treefile specified. Using startree with branchlength of 0.2." << endl;
    }

    try{
      Properties::init( argc, argv );
      Constant::init();
      //      Gene::init();
      GeneticCode::init();
      if (Properties::hasProperty("translation_table")){
	GeneticCode::chooseTranslationTable(Properties::getIntProperty("translation_table"));
      }

      // setParameters(); // NOTE: need Constant and GeneticCode to be initialized first
      StateModel::init();   // set global parameters of state models	  
		
	// calculate omega on input codon alignment
       
	clock_t start;
	start = clock();
	//  CompGenePred cgp;
	//cgp.start();

	CodonMSA cAli(Constant::codonalnfile, branchlength);
	printIntro();
	cAli.printOmegaStats();
	cout << "# total time: " << (double) (clock()-start) / CLOCKS_PER_SEC << "s" << endl;
	
	
	//	if (verbosity>2)
	cout << "# command line:" << endl << "# " << commandline << endl;

    } catch( ProjectError& err ){
        cerr << "\n" <<  argv[0] << ": ERROR\n\t" << err.getMessage( ) << "\n\n";
        errorcode=1;
    } catch ( HelpException help ) {
    	cerr << help.message << endl;
    }
    //    if (outputfile.is_open())
    //	outputfile.close();
    //    if (errorfile.is_open())
    //	errorfile.close();
    return errorcode;
}


void printUsage(){

  cerr << "\nESPOCA - Estimate Selective Pressure on Codon Alignments\n\n";
  cerr << "USAGE:\nespoca [options] --species=SPECIES --alnfile=ALNFILE --treefile=TREEFILE > outfile\n\n";
  cerr << "DESCRIPTON:\n\
 SPECIES   species parameter for calculation of the codon usage. type 'augustus --species=help' to see what species are available\n\
 ALNFILE   codon alignment file in multi fasta format\n\
 TREEFILE  phylogenetic tree with branchlength in newick format (startree is used if not specified)\n\n";
  cerr << "OPTIONS:\n\
 --help         print this usage\n";
}

void printIntro(){

  cout << "#\n# ESPOCA - Estimate Selective Pressure on Codon Alignments.\n#\n";
  cout << "# Description of the table columns:\n# ali_pos   ref_pos   AS_ref    Pr(w>1)   post_mean +-  SE_for_w  num_subst\n";
  cout << "# 1. ali_pos     position of codon site in the alignment\n";
  cout << "# 2. ref_pos     position of codon in reference species (first species in the alignment file), -1 if gap in reference\n";
  cout << "# 3. AS_ref      amino acid of reference sequence at ref_pos\n";
  cout << "# 4. Pr(w>1)     probability of omega > 1 at alipos (*: Pr(w>1) > 0.90, **: Pr(w>1) > 0.95)\n";
  cout << "# 5. post_mean   posterior mean estimate of omega at ali_pos\n";
  cout << "# 6. SE_for_w    standard deviation of omega at ali_pos\n";
  cout << "# 7. num_subst   number of substitution calculated by the Fitch algorithm\n\n";
}

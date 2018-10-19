/*
 * evaluate.cc
 *
 * License: Artistic License, see file LICENSE.TXT or 
 *          https://opensource.org/licenses/artistic-license-1.0
 */

// project includes
#include "gene.hh"
#include "genbank.hh"
#include "evaluation.hh"

// standard C/C++ includes
#include <iostream>


/*
 * main
 */

int main( int argc, char* argv[] ){
    try{

	if( argc <= 2 ){
	    cerr << "Usage:\nevaluate annofile predictfile" << endl;
	    exit(1);
	}
	Evaluation eval;
	GBProcessor gbankanno(argv[1]), gbankpred(argv[2]);
	if (gbankanno.fileType() != genbank || gbankpred.fileType() != genbank)
	    throw ProjectError("Input files not in genbank format.");

	AnnoSequence *annoseq = gbankanno.getAnnoSequenceList();
	AnnoSequence *predseq = gbankpred.getAnnoSequenceList();
	while( annoseq ){
	    if (!predseq || strcmp(predseq->seqname, annoseq->seqname) != 0){
		throw ProjectError("The two input files didn't contain the "
				   "same number of sequences in the same order.");
	    }
	    //cout << "anno:" << endl;
	    //annoseq->anno->printGFF();
	    //cout << "pred:" << endl;
	    //predseq->anno->printGFF();
	    eval.addToEvaluation(predseq->anno->genes, annoseq->anno->genes, bothstrands);
	    predseq = predseq->next;
	    annoseq = annoseq->next;
	}
	eval.print();
    } catch( ProjectError& err ){
        cerr << "\n" <<  argv[0] << ": ERROR\n\t" << err.getMessage( ) << "\n\n";
        return 1;
    }
    return 0;
}

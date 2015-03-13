/**********************************************************************
 * file:    igenictrain.cc
 * licence: 
 *          
 * descr.:  training for intergenic region
 * authors: Mario Stanke, mario@gobics.de
 *
 * date    |   author      |  changes 
 * --------|---------------|------------------------------------------
 * 27.03.06| Mario Stanke  | took introns in cds out of the training
 **********************************************************************/

#include "igenicmodel.hh"

// project includes
#include "gene.hh"
#include "properties.hh"
#include "merkmal.hh"
#include "intronmodel.hh" // so splice sites models can be reused here

// standard C/C++ includes
#include <fstream>
#include <sstream>

/*
 * ===[ IGenicModel::buildModel ]=========================================
 */
void IGenicModel::buildModel( const AnnoSequence* annoseq, int parIndex ){
    if (verbosity)
	cout << " ** building model for intergenic region *IR*" << endl;
    if (!annoseq->anno) 
	throw ProjectError("Tried to train IGenicModel with 0 sequences.");
    
    State* exon, *lastexon;
    const Gene* curgene;
    gesbasen = 0;
    emicount.assign( POWER4TOTHE( k+1 ), 0 );
    int igenic_winlen = 1000000;  // neutral effect
    int igenic_start, igenic_end;
    const char *seqfile = NULL;
    try {
	seqfile = Properties::getProperty("/IGenicModel/seqfile");
    } catch (...) {
    }
	
    if (!seqfile) {
	const AnnoSequence *as = annoseq;
	while (as){
	    const char* sequence = as->sequence;
	    curgene = dynamic_cast<Gene*> (as->anno->genes); // assume here that sequences with multiple genes have already been split
	    if (curgene){
		gweight = curgene->weight;
	    
		if (curgene->exons) { 
		    //cout << "gene " << curgene->id << ", weight " << gweight << " " << curgene->geneBegin() << ".." << curgene->geneEnd() << " seqlen=" << strlen(sequence) << endl;
		    // upstream UTR introns and igenic region
		    igenic_start = curgene->geneBegin() - igenic_winlen;
		    if (igenic_start<0)
			igenic_start=0;
		    processSequence( sequence + igenic_start + k, sequence + curgene->geneBegin()-1-20); // could be more accurate here 
		    if (curgene->utr5exons) {
			lastexon = curgene->utr5exons;
			exon = lastexon->next;
			while (exon) {
			    processSequence (sequence + lastexon->end + 1 + k + Constant::dss_end + DSS_MIDDLE, 
					     sequence + exon->begin - 1 - Constant::ass_start - ASS_MIDDLE 
					     - Constant::ass_upwindow_size);
			    lastexon = exon;
			    exon = exon->next;
			}
		    }
	      
		    // for unknown reasons, things get significantly worse when the intron sequences are used for training
		    /* when using the human training set
		       State *in = curgene->introns;
		       while( in ){
		       processSequence(sequence + in->begin + Constant::dss_end + DSS_MIDDLE + k,
		       sequence + in->end - ASS_MIDDLE - Constant::ass_start -
		       Constant::ass_upwindow_size);
		       in = in->next;
		       }*/
	    
		    // downstream UTR introns and igenic region
		    if (curgene->utr3exons) {
			lastexon = curgene->utr3exons;
			exon = lastexon->next;
			while (exon) {
			    processSequence (sequence + lastexon->end + 1 + k + Constant::dss_end + DSS_MIDDLE, 
					     sequence + exon->begin - 1 - Constant::ass_start - ASS_MIDDLE 
					     - Constant::ass_upwindow_size);
			    lastexon = exon;
			    exon = exon->next;
			}
		    }
		    igenic_end = curgene->geneEnd() + igenic_winlen;
		    if (igenic_end > (int) strlen(sequence)-1)
			igenic_end = (int) strlen(sequence)-1;
		    processSequence( sequence + curgene->geneEnd() +  1 + k,
				 sequence + igenic_end );
		}
	    }
	    as = as->next;
	}
    } else { // read in the intergenic sequences from a separate training file just with ir sequences
	string outfile = Constant::fullSpeciesPath() + seqfile;
	if (verbosity)
	    cout << "Building content model only from sequences in " << outfile << endl;
	ifstream istrm(outfile.c_str());
	char buff[100000];
	if (istrm) {
	    do {
		istrm.getline( buff, 99999 );
		if (strncmp (buff, ">", 1) != 0) {
		    cout << "reading in sequence of length " << strlen(buff) << endl;
		    processSequence(buff + k, buff + strlen(buff)-1);
		}
	    } while( istrm ); 
	    istrm.clear( ios::failbit );
	} else {
	    cerr << "IGenicModel: Couldn't open file " << seqfile << endl;
	}
    }
  
    Pls.resize( k+1 );
    Pls[k].assign( emicount.size(), 0.0);
    makeProbsFromCounts(Pls[k], emicount, k, patpseudocount);
 
    for( int i = k-1; i >= 0; i-- ){
	int size = POWER4TOTHE( i+1 );	
	Pls[i].resize(size);
	for( int j = 0; j < size; j++ ) {
            Double pi1 = Pls[i+1][j];
            Double pi2 = Pls[i+1][j+1*size];
            Double pi3 = Pls[i+1][j+2*size];
            Double pi4 = Pls[i+1][j+3*size];
            Pls[i][j] = pi1+pi2+pi3+pi4;
        }
    }
    emiprobs.probs.resize(Pls[k].size());
    emiprobs.order = k;
    computeEmiFromPat(Pls[k], emiprobs.probs, k);

    lastParIndex = parIndex;
    if (verbosity) {
	cout << "number of bases used for training: " << gesbasen << endl;
	cout << "end *IR*" << endl;
    }
}

/*
 * ===[ IGenicModel::registerPars  ]===================================
 */
void IGenicModel::registerPars (Parameters* parameters){
  try {
    Constant::CRFtrainIgenic = Properties::getBoolProperty("CRFtrainIgenic");
  } catch (...) {}

  if (!parameters || !Constant::CRFtrainIgenic || Constant::tieIgenicIntron)
    return;
  for (int idx=0; idx < Constant::decomp_num_steps; idx++)
    parameters->addMMGroup(&GCemiprobs[idx]);
}


/*
 * ===[ IGenicModel::processSequence ]====================================
 */
void IGenicModel::processSequence( const char* start, const char* end ){
    Seq2Int s2i( k+1 );
    for( ; start <= end; start++){
	try {
	    emicount[ s2i( start-k ) ] += gweight ;
	    emicount[ s2i.rc( start-k ) ] += gweight ;
	    gesbasen += 2*gweight;
	} catch (InvalidNucleotideError e) {}
    }
}

/*
 * storeGCPars
 * store GC dependent model probabilities in the respective array variables
 */
void IGenicModel::storeGCPars(int idx){
  GCemiprobs[idx] = emiprobs;
  std::stringstream out;
  out << "igenic emiprob gc" << (idx+1);
  GCemiprobs[idx].setName(out.str());
  GCPls[idx] = Pls;
}

/*
 * ===[ IGenicModel::printProbabilities ]=================================
 */
void IGenicModel::printProbabilities( int idx, BaseCount *bc, const char* suffix){
    try{  
	const char* fname = Properties::getProperty("/IGenicModel/outfile");
	string filename = Constant::fullSpeciesPath() + fname;
	if (suffix)
	    filename += suffix;

	if (idx < 0) {       // do nothing but deleting the outfile
	    ofstream ostrm( filename.c_str(), ios::trunc);
	    ostrm.close();
	    return;
	}
	ofstream ostrm( filename.c_str(), ios::app);

        if( ostrm ){
            Seq2Int s2i(k+1);
	    if (verbosity)
		cout << "Writing intergenic region model parameters [" << idx+1 << "] to file " 
		     << filename << "." << endl;
	    
	    ostrm << "[" << idx+1 << "]" << endl;
	    if (bc)
		ostrm << "# (a,c,g,t)= " << *bc << endl;
	
            ostrm << "#\n# Probabilities file for the intergenic region model\n#\n" << endl;
            ostrm << "# k =\n" << k << endl;
            ostrm << "\n# The P_l's\n[P_ls]" << endl;
            for( int i = 0; i <= k; i++ ){
                ostrm << "# l=\n" <<  i << endl;
                Seq2Int s2i(i+1);
                ostrm << "# Values" << endl;
                for( int j = 0; j < GCPls[idx][i].size(); j++ )
                    ostrm << s2i.INV(j) << "\t" << GCPls[idx][i][j] << endl;
            }
	    
            // the emission probabilities that are trained through CRF-training
	    if (Constant::tieIgenicIntron && IntronModel::GCemiprobs[gcIdx].probs.size() != GCemiprobs[gcIdx].probs.size())
	      throw ProjectError("Cannot have different order k for igenic and intron model when option tieIgenicIntron is used.");
	    int precision = LLDouble::getOutputPrecision();
	    LLDouble::setOutputPrecision(4);
            ostrm << "\n[EMISSION]" << endl;
            ostrm << "\n# Vector size (4^(k+1))\n" << GCemiprobs[idx].probs.size()
                  << "\n# Probabilities\n";
            for( int i = 0; i < GCemiprobs[idx].probs.size(); i++ )
	      ostrm << s2i.INV(i) << "\t" << (Constant::tieIgenicIntron? IntronModel::GCemiprobs[gcIdx].probs[i] : GCemiprobs[idx].probs[i]) << endl;
	    LLDouble::setOutputPrecision(precision);
	    
	    /************* frequencies of the (k+1)-patterns
	    
	    ostrm << "# patterns" << endl;
	    for( int j = 0; j < emicount.size(); j++ )
		ostrm << s2i.INV(j) << "\t" << emicount[j] << endl;
	    
		//  *************/

            ostrm.close();
        }
    }
    catch( ProjectError ) { 
	if (verbosity)
	    cout << "Intergenic Region model parameters not saved." << endl;
    }
}

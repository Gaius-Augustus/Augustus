/**********************************************************************
 * file:    exontrain.cc
 * licence: 
 *          
 * descr.:  training for exon model parameters
 * authors: Mario Stanke, mario@gobics.de
 *
 * date    |   author      |  changes
 * --------|---------------|------------------------------------------
 * 2.1.06  | Mario Stanke  | creating of file
 **********************************************************************/

#include "exonmodel.hh"

// project includes
#include "commontrain.hh"
#include "geneticcode.hh"
#include "gene.hh"
#include "projectio.hh"
#include "properties.hh"

// standard C/C++ includes
#include <fstream>
#include <sstream>


/*
 * ===[ ExonModel::buildModel ]===========================================
 */
void ExonModel::buildModel( const AnnoSequence* annoseq, int parIndex){
    if (!annoseq->anno) 
	throw ProjectError("Tried to train Exonmodel with 0 sequences.");
    
    if (verbosity)
	cout << " ** building model for exons *EXON*" << endl;
    tiswins = new list<string>;
    gesbasen[0] = gesbasen[1] = gesbasen[2] = 0;
    for (int f=0; f<3; f++) {
	patterncount[f].assign(POWER4TOTHE( k+1 ), 0 );
	initpatterncount[f].assign(POWER4TOTHE( k+1 ), 0 );
	etpatterncount[f].assign(POWER4TOTHE( k+1 ), 0 );    
    }
    transInitMotif = new Motif(trans_init_window, tis_motif_memory, 1, tis_motif_radius);
    if (etMotif){
	delete etMotif[0];
	delete etMotif[1];
	delete etMotif[2];
    }
    etMotif = new Motif*[3];
    etMotif[0] = new Motif(Constant::et_coding_len, etorder, etpseudocount);
    etMotif[1] = new Motif(Constant::et_coding_len, etorder, etpseudocount);
    etMotif[2] = new Motif(Constant::et_coding_len, etorder, etpseudocount);
//    delete baumwelch;
//    baumwelch = new BaumWelch(k, 3, numModels);

    ochrecount = ambercount = opalcount = 0;
    if (!hasLenDist) {
	lenCountSingle.assign(exonLenD+1,0);
	lenCountInitial.assign(exonLenD+1,0);
	lenCountInternal.assign(exonLenD+1,0);
	lenCountTerminal.assign(exonLenD+1,0);
	numExonsOfType.assign(NUM_TYPES, 0);    
	numHugeExonsOfType.assign(NUM_TYPES, 0);    
    }
    const Gene * curgene;
    const AnnoSequence *as = annoseq;
    while (as){
	sequence = as->sequence;
	curgene = dynamic_cast<Gene*> (as->anno->genes); // assume here that sequences with multiple genes have already been split
	if (curgene){ // only training with coding genes
	    gweight = curgene->weight;
	    if (curgene->clength % 3 == 0) {
		try{
		    if (curgene->exons)
			processExons( curgene );
		}
		catch( ExonModelError& exerr ){
		    cerr << "ExonModel::buildModel( " << curgene->id << " ):\t"
			 << exerr.getMessage() << endl << flush;
		}
	    } else {
		if (verbosity)
		    cerr << "gene " << curgene->geneid << " transcr. " << curgene->id << " in sequence " << curgene->seqname << ": " 
			 << "coding length not a multiple of 3. Skipping..." << endl;
	    }
	}
      as = as->next;
    } 

/*
    // AA patterns
    if (!hasAAdep) {
	AADependency aasmall(1);
	aadep.setK(2);
	aadep.makeTransProbs(gene);
	aasmall.makeTransProbs(gene);
	 *
	 *  determine, whether a higher order markov chain brings different emission probabilities
	 *  by computing the Quotients P(A2|A1, A0) / P(A2|A1)
	 *
	for (int pn = 0; pn < aadep.trans.size(); pn++ ){
	    int shortpn;
	    shortpn = pn % ((int) (pow(20, aasmall.k+1) + 0.5));
	    aadep.trans[pn] = aadep.trans[pn] / aasmall.trans[shortpn];
	}
    }
*/

    buildProbabilities( );
    delete tiswins; // free memory
    lastParIndex = parIndex;
}

/*
 * ===[ ExonModel::buildProbabilities ]===================================
 */ 
void ExonModel::buildProbabilities( ){
    //cerr << "entering buildProbabilities" << endl << flush;
 
    int numpatterns = patterncount[0].size();

    /*
     *  motif before translation initiation "ATG"
     */
    transInitMotif->makeProbs();
    //cout << "TransInit Motif" << endl;
    if (Constant::CRFtrainTIS){
      // make TIS CRF features
      transInitBinProbs.reset();
      // fill in orig probs for initial adjustment
      // first the probs of true TIS windows
      for (list<string>::iterator sit = tiswins->begin(); sit != tiswins->end(); sit++){
	transInitBinProbs.addProb(transInitMotif->seqProb(sit->c_str() + tis_motif_memory));
      }
      // then some probs of false TIS windows (10% random sequences)
      for (int i=0; i < 1 + tiswins->size()/10; i++){
	char *r = getRandomDNA(tis_motif_memory + trans_init_window);
	transInitBinProbs.addProb(transInitMotif->seqProb(r + tis_motif_memory));
	delete r;
      }
      transInitBinProbs.trainBins(Constant::tis_maxbinsize);
      transInitBinProbs.printBoundaries();
      //transInitBinProbs.removeOrigs();
    }
    // start codons
    GeneticCode::trainStartCodonProbs(startcounts);
    /*
     * etMotif = motif before the donor splice site.
     * the index is the reading frame at the beginning (left) of the motif
     */
    etMotif[0]->makeProbs();
    etMotif[1]->makeProbs();
    etMotif[2]->makeProbs();    

   /*
    * Emissions probabilities and probabilities for shorter patterns are computed from frequencies 
    * of patterns of length k+1
    * Pls[i][f][p] is the probability of pattern number p of length (i+1) in frame f
    * The frame is the position of the last nucleotide. 
    */
    if (verbosity)
	cout << " number of bases in the reading frames: " 
	     << gesbasen[0]<< " " << gesbasen[1]<< " " << gesbasen[2]<< endl;
  
    Pls.assign( k+1, 3);
 
    for (int f=0; f<3; f++) {
	Pls[k][f].resize(numpatterns);
	if (verbosity > 1)
	    cout << "--- frame = " << f << " ---    minPatSum = " << minPatSum << endl;
	if (minPatSum > 0) { 
	    determineShortPatterns(patterncount[f], k, minPatSum);
	    makeProbsFromCounts(Pls[k][f], patterncount[f], k, patpseudo, true);
	} else 
	    makeProbsFromCounts(Pls[k][f], patterncount[f], k, patpseudo, false);
    }
    
    // initial patterns, only emission probabilities
    for (int f=0; f<3; f++) {
	vector<Double> tempinitpatprobs(numpatterns);
	initemiprobs[f].resize(numpatterns);
	if (verbosity > 1)
	    cout << "--- initial frame = " << f << " ---    minPatSum = " << minPatSum << endl;
	if (minPatSum > 0) { 
	    determineShortPatterns(initpatterncount[f], k, minPatSum);
	    makeProbsFromCounts(tempinitpatprobs, initpatterncount[f], k, patpseudo, true);
	} else 
	    makeProbsFromCounts(tempinitpatprobs, initpatterncount[f], k, patpseudo, false);
	    
	    // make emission probabilities of the pattern probabilities and discard the latters
	computeEmiFromPat(tempinitpatprobs, initemiprobs[f], k);	
    }
    // end of initial patterns

    // internal exon terminal patterns, only emission probabilities
    for (int f=0; f<3; f++) {
	vector<Double> tempetpatprobs(numpatterns);
	etemiprobs[f].resize(numpatterns);
	if (verbosity > 1)
		cout << "--- internal exon terminal frame = " << f << " ---    minPatSum = " << minPatSum << endl;
	    if (minPatSum > 0) { 
		determineShortPatterns(etpatterncount[f], k, minPatSum);
		makeProbsFromCounts(tempetpatprobs, etpatterncount[f], k, patpseudo, true);
	    } else 
		makeProbsFromCounts(tempetpatprobs, etpatterncount[f], k, patpseudo, false);
	    // make emission probabilities of the pattern probabilities and discard the latters
	    computeEmiFromPat(tempetpatprobs, etemiprobs[f], k);	
    }
    // end of internal exon terminal patterns
    

    /* 
     * Compute lower order pattern probs from the ones above,
     * eg AAAA from AAAAA AAAAC AAAAG AAAAT
     */
    computeLowerOrderPats();

    emiprobs.probs[0].resize(Pls[k][0].size());
    emiprobs.probs[1].resize(Pls[k][1].size());
    emiprobs.probs[2].resize(Pls[k][2].size());
    emiprobs.order = k;
    computeEmiFromPat(Pls[k][0], emiprobs.probs[0], k);
    computeEmiFromPat(Pls[k][1], emiprobs.probs[1], k);
    computeEmiFromPat(Pls[k][2], emiprobs.probs[2], k);


    if (!hasLenDist) 
	computeLengthDistributions();

    /*
     * Print out the numer of exons of each of the 8 types.
     * Only needed for making the transition matrix by hand.
     */
    if (verbosity)
	for (int type=singleG; type <= terminal ; type++) {
	    cout << stateTypeNames[type] << " : " << numExonsOfType[type] << endl;
	}
    hasLenDist = true;
    if (verbosity>1){
	int totalcount = ochrecount + ambercount + opalcount;
	cout << "Frequency of stop codons:" << endl;
	if (totalcount>0) {
	    cout << "tag:" << setw(5) << ambercount << " (" << setprecision(3) << (float) ambercount/totalcount << ")" << endl;
	    cout << "taa:" << setw(5) << ochrecount << " (" << setprecision(3) << (float) ochrecount/totalcount << ")" << endl;
	    cout << "tga:" << setw(5) << opalcount << " (" << setprecision(3) << (float) opalcount/totalcount << ")" << endl;
	}
    }

    if (verbosity)
	cout << "end *EXON*" << endl;
}

/*
 * ===[ ExonModel::registerPars  ]===================================
 */
void ExonModel::registerPars (Parameters* parameters){
  int tisalpha = 30;
  try {
    tisalpha = Properties::getIntProperty("/ExonModel/tisalpha");
  } catch (...){}
  if (!parameters)
    return;
  if (Constant::CRFtrainTIS)
    for (int idx=0; idx < Constant::decomp_num_steps; idx++)
      parameters->addMMGroup(&GCtransInitBinProbs[idx], tisalpha);
  if (Constant::CRFtrainCDS)
    for (int idx=0; idx < Constant::decomp_num_steps; idx++)
      parameters->addMMGroup(&GCemiprobs[idx]);
}

/*
 * ===[ ExonModel::computeLowerOrderPats ]===================================
 */
void ExonModel::computeLowerOrderPats(){
   /* 
     * Compute lower order pattern probs Pls[i] (i=0..k-1) from the Pls[k]
     * eg AAAA from AAAAA CAAAA GAAAA TAAAA
     */
    int size;
    for( int i = k-1; i >= 0; i-- ){
	size = (int) (POWER4TOTHE( i+1 )+0.1);
        Pls[i][0].resize(size);
	Pls[i][1].resize(size);
	Pls[i][2].resize(size);
	for( int j = 0; j < size; j++ ){
            Double pi1 = Pls[i+1][0][j];
            Double pi2 = Pls[i+1][0][j + 1*size];
            Double pi3 = Pls[i+1][0][j + 2*size];
            Double pi4 = Pls[i+1][0][j + 3*size];
            Pls[i][0][j] = pi1+pi2+pi3+pi4;
            pi1 = Pls[i+1][1][j];
            pi2 = Pls[i+1][1][j + 1*size];
            pi3 = Pls[i+1][1][j + 2*size];
            pi4 = Pls[i+1][1][j + 3*size];
            Pls[i][1][j] = pi1+pi2+pi3+pi4;
	    pi1 = Pls[i+1][2][j];
            pi2 = Pls[i+1][2][j + 1*size];
            pi3 = Pls[i+1][2][j + 2*size];
            pi4 = Pls[i+1][2][j + 3*size];
            Pls[i][2][j] = pi1+pi2+pi3+pi4;
        }
    }
}


/*
 * ===[ ExonModel::storeGCPars ]=======================================
 * store GC dependent model probabilities in the respective array variables
 */
void ExonModel::storeGCPars(int idx){
  GCPls[idx] = Pls;
  GCemiprobs[idx] = emiprobs;
  std::stringstream out;
  out << "exon emiprob gc" << (idx+1);
  GCemiprobs[idx].setName(out.str());
  GCinitemiprobs[idx][0] = initemiprobs[0];
  GCinitemiprobs[idx][1] =  initemiprobs[1];
  GCinitemiprobs[idx][2] =  initemiprobs[2];
  GCetemiprobs[idx][0] =  etemiprobs[0];
  GCetemiprobs[idx][1] =  etemiprobs[1];
  GCetemiprobs[idx][2] =  etemiprobs[2];
  if (transInitMotif)
    GCtransInitMotif[idx] =  *transInitMotif;
  if (transInitBinProbs.nbins>0){
    GCtransInitBinProbs[idx] = transInitBinProbs;
    std::stringstream id;
    id << "tis bin gc" << (idx+1);
    GCtransInitBinProbs[idx].setName(id.str());
					    }
  if (etMotif){
    *GCetMotif[idx][0] = *etMotif[0];
    *GCetMotif[idx][1] = *etMotif[1];
    *GCetMotif[idx][2] = *etMotif[2];
  }
}

/*
 * printProbabilities
 * idx - the index of the GC content class
 * suffix - an identifyer to append to the parameter file name in case there are multiple variants generated during training
 */
void ExonModel::printProbabilities(int idx, BaseCount *bc, const char* suffix){
    Seq2Int s2i_e(k+1);

        try{ 
	    const char* fname = Properties::getProperty( "/ExonModel/outfile" );
	    string filename = Constant::fullSpeciesPath() + fname;
	    if (suffix)
		filename += suffix;
	    if (idx < 0) { // do nothing but deleting the outfile
		ofstream ofstrm(filename.c_str(), ios::trunc);
		ofstrm.close();
		return;
	    }
            ofstream ofstrm(filename.c_str(), ios::app);
	    if (verbosity)
  	        cout << "Writing exon model parameters [" << idx+1 << "] to file " << filename << "." << endl;
//	    LLDouble::setOutputPrecision(3);
	    if (idx == 0 ) { 
		ofstrm << "#exon model parameters\n# begin of content independent part" << endl;
		
		ofstrm << "# start codon probabilities\n[STARTCODONS]" << endl;
		GeneticCode::writeStart(ofstrm);

		ofstrm << "\n# Length distributions\n[LENGTH]" << endl;
		ofstrm << "# maximal individually stored length probability =\n" << exonLenD << endl;
		ofstrm << "# slope of smoothing bandwidth =\n" << slope_of_bandwidth << endl;
		ofstrm << "# smoothing minwindowcount =\n" << minwindowcount << endl;
		ofstrm << "# length single  initial  internal  terminal" << endl;
		ofstrm << "# total number of exons of above types" << endl;
		ofstrm << "       " << numSingle << setw(15) << numInitial << setw(15) 
		       << numInternal <<setw(15) << numTerminal << endl;
		ofstrm << "# number of exons exceeding length d" << endl;
		ofstrm << "       " << numHugeSingle << setw(15) << numHugeInitial << setw(15) 
		       << numHugeInternal << setw(15) << numHugeTerminal << endl;
		
		ofstrm << "# 1000 P(len=k), k=0,1,..., " << exonLenD << endl;
		for( int i = 0; i <= exonLenD; i++ ){
		    ofstrm << i 
			   << "\t" << 1000 * lenDistSingle[i] << "\t" << 1000 * lenDistInitial[i] 
			   << "\t" << 1000 * lenDistInternal[i] << "\t" << 1000 * lenDistTerminal[i] 
			   << endl;
		}

/*
		ofstrm << "\n# actual lengths " << endl;
		for( int i = 0; i <= exonLenD; i++ ){
		    ofstrm << i 
			   << "\t" << lenCountSingle[i] << "\t" << lenCountInitial[i] 
			   << "\t" << lenCountInternal[i] << "\t" << lenCountTerminal[i] 
			   << endl;
		}

*/
		ofstrm << "# end of content independent part" << endl;
	    }

	    ofstrm << "\n# data set number\n[" << idx+1 << "]" << endl;
	    if (bc)
		ofstrm << "# (a,c,g,t)= " << *bc << endl;
            ofstrm << "#\n# Probabilities file for the exon model\n#\n"
                   << endl;


            ofstrm << "\n# Die P_l's\n[P_ls]" << endl;
            ofstrm << "# k = " << k << endl;
            for( int i = 0; i <= k; i++ ){
                ofstrm << "# l=\n" <<  i << endl;
                Seq2Int s2i(i+1);
                ofstrm << "# Values" << endl;
                for( int j = 0; j < GCPls[idx][i][0].size(); j++ )
                    ofstrm << s2i.inv(j) << "\t" << GCPls[idx][i][0][j]
                           << "\t     " << GCPls[idx][i][1][j] << "\t     "
                           << GCPls[idx][i][2][j] << endl;
            }


	    /************* frequencies of the (k+1)-patterns

	    Seq2Int s2i(k+1);
	    ofstrm << "\n# normal patterns" << endl;
	    for( int j = 0; j < patterncount[0].size(); j++ )
		ofstrm << s2i.inv(j) << "\t" << patterncount[0][j]
		       << "\t     " << patterncount[1][j] << "\t     "
		       << patterncount[2][j] << endl;
	    
	    ofstrm << "\n# initial patterns" << endl;
	    for( int j = 0; j < initpatterncount[0].size(); j++ )
		ofstrm << s2i.inv(j) << "\t" << initpatterncount[0][j]
		       << "\t     " << initpatterncount[1][j] << "\t     "
		       << initpatterncount[2][j] << endl;
	    
	    ofstrm << "\n# internal exon terminal patterns" << endl;
	    for( int j = 0; j < etpatterncount[0].size(); j++ )
		ofstrm << s2i.inv(j) << "\t" << etpatterncount[0][j]
		       << "\t     " << etpatterncount[1][j] << "\t     "
		       << etpatterncount[2][j] << endl;

		       //  ************/


	    ofstrm << "\n# translation initiation motif\n[TRANSINIT]" << endl;
	    GCtransInitMotif[idx].write(ofstrm);
	    if (GCtransInitBinProbs[idx].nbins > 0){
	      ofstrm << "\n#\n# TIS binned probabilities\n[TRANSINITBIN]\n";
	      GCtransInitBinProbs[idx].write(ofstrm);
	    }
	    ofstrm << "\n# dss upstream motif, reading frame 0(reverse)\n[ETMOTIF0]" << endl;
	    GCetMotif[idx][0]->write(ofstrm);
	    ofstrm << "\n# dss upstream motif, reading frame 1(reverse)\n[ETMOTIF1]" << endl;
	    GCetMotif[idx][1]->write(ofstrm);
	    ofstrm << "\n# dss upstream motif, reading frame 2(reverse)\n[ETMOTIF2]" << endl;
	    GCetMotif[idx][2]->write(ofstrm);
	    
/*
	    if (verbosity) 
		cout << "Writing amino acid patterns to file " << fname << "." << endl;
	    aadep.printTrans(ofstrm);
*/
	    // for human readability only: emission probabilities 

	    int precision = LLDouble::getOutputPrecision();
	    LLDouble::setOutputPrecision(4);
	    ofstrm << "\n\n#\n# Emission probabilities\n#\n[EMISSION]\n";
	    ofstrm << "# Size of vector\n" << GCemiprobs[idx].probs[0].size() << endl;
	    ofstrm << "# k : order of the markov model\n" << k << endl;   
	    ofstrm << "# patpseudocount (pseudocount of sequence patterns)\n" << patpseudo << endl;
	    ofstrm << "# Probabilities\n# Format: pattern win0 win1 win2\n";
	    for( int i = 0; i < GCemiprobs[idx].probs[0].size(); i++ ){
		ofstrm << s2i_e.inv(i);
		for (int f=0; f<3; f++)
		    ofstrm << "\t" << GCemiprobs[idx].probs[f][i];
		ofstrm << endl;
	    }
	    LLDouble::setOutputPrecision(precision);
	    
	    ofstrm << "\n\n#\n# Initial emission probabilities\n#\n[INITEMISSION]\n";
	    ofstrm << "# Size of vector\n" << GCinitemiprobs[idx][0].size() << endl;
	    ofstrm << "# k : order of the markov model\n" << k << endl;   
	    ofstrm << "# patpseudocount (pseudocount of sequence patterns)\n" << patpseudo << endl;
	    ofstrm << "# Probabilities\n# Format: pattern win0 win1 win2\n";
	    for( int i = 0; i < GCinitemiprobs[idx][0].size(); i++ ){
		ofstrm << s2i_e.inv(i);
		for (int f=0; f<3; f++)
		    ofstrm << "     \t" << GCinitemiprobs[idx][f][i];
		ofstrm << endl;
	    }
	    
	    ofstrm << "\n\n#\n# Internal exon terminal emission probabilities\n#\n[ETEMISSION]\n";
	    ofstrm << "# Size of vector\n" << GCetemiprobs[idx][0].size() << endl;
	    ofstrm << "# k : order of the markov model\n" << k << endl;   
	    ofstrm << "# patpseudocount (pseudocount of sequence patterns)\n" << patpseudo << endl;
	    ofstrm << "# Probabilities\n# Format: pattern win0 win1 win2\n";
	    for( int i = 0; i < GCetemiprobs[idx][0].size(); i++ ){
		ofstrm << s2i_e.inv(i);
		for (int f=0; f<3; f++)
		    ofstrm << "     \t" << GCetemiprobs[idx][f][i];
		ofstrm << endl;
	    }



	    /*
	    ofstrm << "# number of reading frames\n" << 3 << endl;
	    ofstrm << "# number of types\n" << numModels << endl;
	    ofstrm << "# probabilities of the emissions" << endl;
	    for (int i=0; i < emission[0].getRowSize(); i++) {
		ofstrm << s2i_e.inv(i) << "  ";
		for (int t=0; t < numModels; t++) {
		    for (int f=0; f < 3; f++) {
			ofstrm << setw(12) << emission[t][f][i];
		    }
		    if (t < numModels) {
			ofstrm << "  ";
		    }
		}
		ofstrm << endl;
	    }
	    */
	    ofstrm.close();
	}
	catch( ProjectError ) {
	    if (verbosity)    
		cout << "Exon model parameters not saved." << endl;
	}
}



/*
 * ===[ ExonModel::computeLengthDistributions]====================================
 */
void ExonModel::computeLengthDistributions( ){

    Smooth sm(minwindowcount, slope_of_bandwidth);
  
    /*
     * compute length distributions
     *
     */
    lenDistSingle.assign(Constant::Constant::max_exon_len + 1, 0.0);
    lenDistInitial.assign(Constant::max_exon_len + 1, 0.0);
    lenDistInternal.assign(Constant::max_exon_len + 1, 0.0);
    lenDistTerminal.assign(Constant::max_exon_len + 1, 0.0);
    

    numHugeSingle   = numHugeExonsOfType[singleG];
    numHugeInitial  = numHugeExonsOfType[initial0] +  numHugeExonsOfType[initial1] 
	+ numHugeExonsOfType[initial2];
    numHugeInternal = numHugeExonsOfType[internal0] +  numHugeExonsOfType[internal1] 
	+ numHugeExonsOfType[internal2];
    numHugeTerminal = numHugeExonsOfType[terminal];

    numSingle   = numExonsOfType[singleG];
    numInitial  = numExonsOfType[initial0] +  numExonsOfType[initial1] + numExonsOfType[initial2];
    numInternal = numExonsOfType[internal0] +  numExonsOfType[internal1] + numExonsOfType[internal2];
    numTerminal = numExonsOfType[terminal];
   
    // compute mean exon lengths
    if (verbosity) {
	long int sum;
	int num, i;
	cout << "single, initial, internal, terminal mean exon lengths :" << endl;
	sum = 0; num = 0;
	for (i=0; i<= exonLenD; i++) {
	    sum += i*lenCountSingle[i]; num += lenCountSingle[i];
	}
	if (num>0)
	    cout << sum/num << "\t";
	else 
	    cout << "n.a." << "\t";
	sum = 0; num = 0;
	for (i=0; i<= exonLenD; i++) {
	    sum += i*lenCountInitial[i]; num += lenCountInitial[i];
	}
	if (num>0)
	    cout << sum/num << "\t";
	else 
	    cout << "n.a." << "\t";
	sum = 0; num = 0;
	for (i=0; i<= exonLenD; i++) {
	    sum += i*lenCountInternal[i]; num += lenCountInternal[i];
	}
	if (num>0)
	    cout << sum/num << "\t";
	else 
	    cout << "n.a." << "\t";
	sum = 0; num = 0;
	for (i=0; i<= exonLenD; i++) {
	    sum += i*lenCountTerminal[i]; num += lenCountTerminal[i];
	}
	if (num>0)
	    cout << sum/num << endl;
	else 
	    cout << "n.a." << endl;
    }


/*
    if (numHugeSingle>0) {
	cerr << numHugeSingle << " single exons were "
	     << " longer than " << exonLenD << endl;
    }
    if (numHugeInitial>0) {
	cerr << numHugeInitial << " initial exons were "
	     << " longer than " << exonLenD << endl;
    }
    if (numHugeInternal>0) {
	cerr << numHugeInternal << " internal exons were "
	     << " longer than " << exonLenD << endl;
    } 
    if (numHugeTerminal>0) {
	cerr << numHugeTerminal << " terminal exons were "
	     << " longer than " << exonLenD << endl;
    }
*/
    sm.smoothCounts(lenCountSingle, lenDistSingle);
    sm.smoothCounts(lenCountInitial, lenDistInitial);
    sm.smoothCounts(lenCountInternal, lenDistInternal);
    sm.smoothCounts(lenCountTerminal, lenDistTerminal);

    // make single exons shorter than the minimal coding length impossible in the first place
    for (int i=0; i < Constant::min_coding_len && i < lenDistSingle.size(); i++)
	lenDistSingle[i] = 0.0;

    scaleDblVector(lenDistSingle, Double(numSingle-numHugeSingle) / numSingle);
    scaleDblVector(lenDistInitial, Double(numInitial-numHugeInitial) / numInitial);
    scaleDblVector(lenDistInternal, Double(numInternal-numHugeInternal) / numInternal);
    scaleDblVector(lenDistTerminal, Double(numTerminal-numHugeTerminal) / numTerminal);

    fillTailsOfLengthDistributions();
}



/*
 * ===[ ExonModel::processExons ]=========================================
 */
void ExonModel::processExons( const Gene* gene ){
    State   *exon = gene->exons;
    curwin = 0;
    if( !exon->next ){
	try{
	    processSingleExon( exon );
	} catch (ExonModelError e) {
	    if (verbosity)
	      cerr << "gene " << gene->geneid << " transcr. " << gene->id << " in sequence " << gene->seqname << ": " << e.getMessage() << endl;
	}
    }
    else{
	try{
	    processInitialExon( exon );
	} catch( ExonModelError e ){
	    if (verbosity)
	      cerr << "gene " << gene->geneid << " transcr. " << gene->id << " in sequence " << gene->seqname << ": " << e.getMessage() << endl;
	}
	exon = exon->next;
	//prewin=-1;
	while( exon->next ){
	    try{
		processInternalExon( exon );
	    } catch( ExonModelError e ){
	      if (verbosity)
		    cerr << "gene " << gene->geneid << " transcr. " << gene->id << " in sequence " << gene->seqname << ": " << e.getMessage() << endl;
	    }
	    exon = exon->next;
	} 
	try{
	    processTerminalExon( exon );
	} catch( ExonModelError e ){
	    if (verbosity)
		cerr << "gene " << gene->geneid << " transcr. " << gene->id << " in sequence " << gene->seqname << ": " << e.getMessage() << endl;
	}
    }
}


/*
 * ===[ ExonModel::processSingleExon ]====================================
 */
void ExonModel::processSingleExon( const State* exon) {
    int length = exon->end-exon->begin+1; 
    if (length < STARTCODON_LEN + STOPCODON_LEN) {
	throw ExonModelError("Single training exon too short.");
    } 
    if (!GeneticCode::isStartcodon(sequence+exon->begin)) {
      string msg("Single exon gene does not begin with start codon but with ");
      msg += string(sequence+exon->begin, 3);
	throw ExonModelError(msg);
    }
    startcounts[Seq2Int(3)(sequence+exon->begin)]++; // count start codon
    
    if (length > trans_init_window && exon->begin >= trans_init_window + tis_motif_memory){
	transInitMotif->addSequence(sequence + exon->begin - trans_init_window, gweight);
	tiswins->push_back(string(sequence + exon->begin - trans_init_window - tis_motif_memory, 
				  tis_motif_memory + trans_init_window));
    }
    curwin = STARTCODON_LEN + k;

    const char *beginOfInnerSeq, *endOfInnerSeq, *endOfInit;
    beginOfInnerSeq =  sequence + exon->begin + STARTCODON_LEN + k;
    endOfInnerSeq = sequence + exon->end - STOPCODON_LEN;
    endOfInit = beginOfInnerSeq + Constant::init_coding_len - 1;

    processInnerSequence(beginOfInnerSeq, endOfInit, 2); // init part
    processInnerSequence(endOfInit+1, endOfInnerSeq, 0); // normal 

    int stppos = exon->end - STOPCODON_LEN + 1;
    if (strncmp(sequence + stppos, "taa", 3)==0)
	ochrecount++;
    else if (strncmp(sequence + stppos, "tag", 3)==0)
	ambercount++;
    else if (strncmp(sequence + stppos, "tga", 3)==0)
	opalcount++;
    else 
	throw ExonModelError("Single exon doesn't end in stop codon. Variable stopCodonExcludedFromCDS set right?");
    
    if (!hasLenDist) {
      if (length <= exonLenD) {
	lenCountSingle[length] += 1; // gweight;
      } else {
	numHugeExonsOfType[singleG] += 1; // gweight;
      }
      numExonsOfType[singleG] += 1; // gweight ;
    }
}


/*
 * ===[ ExonModel::processInitialExon ]=====================================
 */
void ExonModel::processInitialExon( const State* exon){
    int oldwin = curwin;
    int length = exon->end-exon->begin+1; 
    

    if( exon->end - exon->begin + 1 < STARTCODON_LEN )
    {
	curwin = curwin + length;
	throw ExonModelError( "Initial exon has length < 3!" );
    }
    if (!GeneticCode::isStartcodon(sequence+exon->begin)) {
	curwin = curwin + length;
	string msg("Initial exon does not begin with start codon but with ");
	msg += string(sequence+exon->begin, 3);
	throw ExonModelError(msg);
    } 
    startcounts[Seq2Int(3)(sequence+exon->begin)]++; // count start codon

    if (length > trans_init_window && exon->begin >= trans_init_window + tis_motif_memory){
	transInitMotif->addSequence(sequence + exon->begin - trans_init_window, gweight);
	tiswins->push_back(string(sequence + exon->begin - trans_init_window - tis_motif_memory, 
				  tis_motif_memory + trans_init_window));
		
    }
    
    if (exon->end - Constant::dss_start - Constant::et_coding_len + 1 >=0 )
	etMotif[mod3(oldwin + length - Constant::dss_start - Constant::et_coding_len)]->addSequence(
	    sequence + exon->end - Constant::dss_start - Constant::et_coding_len + 1, gweight, true);

    /* test output for seqlogo/pictogram
    static char *fenster = new char[Constant::et_coding_len+1];
    strncpy( fenster, sequence + exon->end - Constant::dss_start - Constant::et_coding_len + 1, Constant::et_coding_len+1 );
    fenster[Constant::et_coding_len]='\0';
    cout << "initial  etMotif " << mod3(oldwin + length - Constant::dss_start - Constant::et_coding_len) << " " <<fenster << endl;
    */


    curwin = curwin + STARTCODON_LEN + k;
    try {
	const char *beginOfInnerSeq, *endOfInit, *endOfInnerSeq;
	beginOfInnerSeq = sequence + exon->begin + STARTCODON_LEN + k;
	endOfInnerSeq = sequence + exon->end - Constant::dss_start;
	endOfInit = beginOfInnerSeq + Constant::init_coding_len - 1;
	if (endOfInit > endOfInnerSeq) // at most up to the end of the sequence
	    endOfInit = endOfInnerSeq;
	
	processInnerSequence(beginOfInnerSeq, endOfInit, 2); // init part
	processInnerSequence(endOfInit+1, endOfInnerSeq);    // normal

    } catch( ExonModelError& e ){
	cerr << "ExonModel::processInitialExon: " << e.getMessage() << endl;
	throw e;
    }
    curwin = oldwin + length;

    if (!hasLenDist) {
	if (length <= exonLenD) {
	    lenCountInitial[length] += 1; // gweight;
	} else {
	    numHugeExonsOfType[initialExon(mod3(curwin))] += 1; //gweight;
	}  
	numExonsOfType[initialExon(mod3(curwin))] += 1; //gweight;
    }
}


/*
 * ===[ ExonModel::processInternalExon ]=====================================
 *
 * cumcodinglength is the coding length of all exons before this one and is needed 
 * for training the init content models
 */
void ExonModel::processInternalExon( const State* exon){
    int oldwin = curwin;
    int length = exon->end-exon->begin+1;
    
    curwin += Constant::ass_end + k;  
    try{

	const char *beginOfInnerSeq, *beginOfET, *endOfInnerSeq;
	beginOfInnerSeq = sequence + exon->begin + Constant::ass_end + k;
	endOfInnerSeq = sequence + exon->end - Constant::dss_start;
	beginOfET = endOfInnerSeq - Constant::et_coding_len + 1;
	if (beginOfET < beginOfInnerSeq)
	    beginOfET = beginOfInnerSeq;
	
	processInnerSequence(beginOfInnerSeq, beginOfET - 1);  // normal
	processInnerSequence(beginOfET, endOfInnerSeq, 3);     // internal exon term part

	etMotif[mod3(oldwin + length - Constant::dss_start - Constant::et_coding_len)]->addSequence(
	    sequence + exon->end - Constant::dss_start - Constant::et_coding_len + 1, gweight, true);
	/* test output for seqlogo/pictogram
	static char *fenster = new char[Constant::et_coding_len+1];
	strncpy( fenster, sequence + exon->end - Constant::dss_start - Constant::et_coding_len + 1, Constant::et_coding_len+1 );
	fenster[Constant::et_coding_len]='\0';
	cout << "internal etMotif " << mod3(oldwin + length - Constant::dss_start - Constant::et_coding_len) << " " <<fenster << endl;
	*/
	
    } catch( ExonModelError e ){
	cerr << "ExonModel::processInternalExon: " << e.getMessage() << endl;
	throw e;
    }
    
    curwin = oldwin + length;
    //transition of phases: cout << "internal " << mod3(oldwin) << "->" << mod3(curwin) << endl;
    
/* for determining transition probabilities
    if (prewin>=0) {
	cout << prewin << curwin%3 << endl;
    } 
    prewin=curwin%3;*/

    if (!hasLenDist) {
	if (length <= exonLenD) {
	    lenCountInternal[length] += 1; // gweight;
	} else {
	    numHugeExonsOfType[internalExon(mod3(curwin))] += 1; //gweight;
	}
	numExonsOfType[internalExon(mod3(curwin))] += 1; //gweight;
    }
}


/*
 * ===[ ExonModel::processTerminalExon ]=======================================
 */
void ExonModel::processTerminalExon( const State* exon){
    int length = exon->end-exon->begin+1; 
    const char *beginOfInnerSeq, *endOfInnerSeq;

    curwin += Constant::ass_end + k;
    
    beginOfInnerSeq = sequence + exon->begin + Constant::ass_end + k;
    endOfInnerSeq = sequence + exon->end - STOPCODON_LEN;

    processInnerSequence( beginOfInnerSeq, endOfInnerSeq);

    int stppos = exon->end - STOPCODON_LEN + 1;
    if (strncmp(sequence + stppos, "taa", 3)==0)
	ochrecount++;
    else if (strncmp(sequence + stppos, "tag", 3)==0)
	ambercount++;
    else if (strncmp(sequence + stppos, "tga", 3)==0)
	opalcount++;
    else 
	throw ExonModelError("Terminal exon doesn't end in stop codon. Variable stopCodonExcludedFromCDS set right?");

    if (!hasLenDist) {
	if (length <= exonLenD) {
	    lenCountTerminal[length] += 1; // gweight;
	} else {
	    numHugeExonsOfType[terminal] += 1; // gweight;
	}
	numExonsOfType[terminal] += 1; // gweight;
    }
}


/*
 * ===[ ExonModel::processInnerSequence ]=================================
 * 'begin' is the first nucleotide, of which the emission is counted
 * 'end' is the last
 */
void ExonModel::processInnerSequence(const char *begin, const char* end, int modeltype){
    if (begin > end)
	return;
    Seq2Int s2i(k+1);
    for( ; begin <= end; begin++){
	curwin %= 3;
	if (curwin==0 && end-begin>=STOPCODON_LEN && GeneticCode::isStopcodon(begin)){
	    throw ExonModelError("in-frame stop codon");
	}
	try {
	    int pn = s2i(begin-k);
	    if (modeltype == 0)
		patterncount[ curwin ][pn] += gweight;
	    else if (modeltype == 2 ) // initial part of the coding region
		initpatterncount[ curwin ][pn] += gweight;
	    else if (modeltype == 3 ) // terminal part of an internal exon
		etpatterncount[ curwin ][pn] += gweight;
	    else 
		throw ProjectError("ExonModel::processInnerSequence: invalid model type");
	    gesbasen[curwin] += gweight;
	} catch (InvalidNucleotideError e) {}
	curwin++;
    }
}



/*
 * ===[ ExonModel::computeCMFromBW ]===================================
 *

void ExonModel::computeCMFromBW(){
     *
     * Emission Probabilities for order k emissions only
     * They are computed using the content models from the baumwelch object
     *
    Matrix< vector<Double> > emission(numModels, 3);

    for (int t=0; t<numModels; t++) {
    	ContentModel cm = baumwelch->getContentModel(t);
	for (int f=0; f<3; f++) {
	    emission[t][f].resize(cm.getNumPatterns());
	    computeEmiFromPat(cm.copyPatProb(f), emission[t][f], k);
	}
    }

     *
     * Initialize modelStartProbs - the vector of probabilities of the
     * different content models.
     *
    if (modelStartProbs)
	delete [] modelStartProbs;
    modelStartProbs = new Double[numModels];
    for (int t=0; t<numModels; t++) {
	modelStartProbs[t] = baumwelch->getModelTypeProb(t);
    }
}
*/

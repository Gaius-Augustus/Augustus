/**********************************************************************
 * file:    intronmodel.cc
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  
 * authors: Mario Stanke, mario@gobics.de
 *
 * date    |   author      |  changes 
 * --------|---------------|------------------------------------------
 * 06.05.02| Mario Stanke  | debugging and testing
 * 23.08.02| Mario Stanke  | performance improvement in seqProb 
 * 01.11.02| Mario Stanke  | introducing emiProbUnderModel
 * 07.11.02| Mario Stanke  | fixing length distr. bug in equalD intron
 * 02.01.03| Mario Stanke  | changing intron to model with 4 states (overlapping short and long introns)
 * ??.01.03| Mario Stanke  | changing intron to model with 5 states 
 * 26.03.03| Mario Stanke  | don't distinguish 'long' and 'short' splice sites anymore
 * 26.03.03| Mario Stanke  | introducing shadow states on reverse strand
 * 05.08.03| Mario Stanke  | check for in-frame stop codons in dss
 * 18.03.04| Mario Stanke  | bugfix in aSSProb: set emiProb before returning 0
 * 19.05.05| Mario Stanke  | bugfix in aSSProb: set oldpos before returning 0
 * 29.06.06| Mario Stanke  | bugfix for extremely short intron hints
 * 09.05.07| Mario Stanke  | bufix in aSSProb: missing re-initialization lead to crashes
 * 21.05.07| Mario Stanke  | bufix in emiProbUnderModel, allow right-truncated short introns on the reverse strand
 **********************************************************************/

#include "intronmodel.hh"

// project includes
#include "properties.hh"
#include "projectio.hh"  // for comment, goto_line_after
#include "motif.hh"
#include "extrinsicinfo.hh"
#include "pp_scoring.hh"

// standard C/C++ includes
#include <fstream>
#include <climits>

/*
 * Initialisation of static data members
 */
vector<Integer> IntronModel::emicount;
vector<Integer> IntronModel::intlencount;
Integer         IntronModel::introns = 0;    // number of introns
Integer         IntronModel::introns_d = 0;  // number of introns of length <= d
Integer         IntronModel::k = 4;
double          IntronModel::slope_of_bandwidth = 0.1;
Integer         IntronModel::minwindowcount = 1;
Integer         IntronModel::d = 500;
Double          IntronModel::patpseudo = 0;
Double          IntronModel::probShortIntron = 0;
Double*         IntronModel::GCprobShortIntron = NULL;
Double          IntronModel::mal = 1210;
Double*         IntronModel::GCmal = NULL;
Integer         IntronModel::introncount = 0;
PatMMGroup      IntronModel::emiprobs("intron emiprob");
PatMMGroup*     IntronModel::GCemiprobs = NULL;
vector<Double>  IntronModel::lenDist;
Integer         IntronModel::gesbasen = 0;
Integer         IntronModel::ass_motif_memory = 3;
Integer         IntronModel::ass_motif_radius = 3;
double          IntronModel::non_gt_dss_prob = 0.001;
double          IntronModel::non_ag_ass_prob = 0.001;
double          IntronModel::geoProb = 0.9997447;


Integer IntronModel::c_ass = 0;
Integer IntronModel::c_dss = 0;
vector<Integer> IntronModel::asscount;
vector<Integer> IntronModel::dsscount;
Integer         IntronModel::ass_upwindow_size = 20;
vector<Double> IntronModel::assprobs;
vector<Double> IntronModel::dssprobs;
Motif*         IntronModel::assMotif = NULL;            // motif window before the assy
Motif*         IntronModel::GCassMotif = NULL;          // array of Motifs, one for each GC content class
Boolean        IntronModel::hasSpliceSites = false;
Double         IntronModel::asspseudo = .1;                  // pseudocount for patterns in acceptor splice sites
Double         IntronModel::dsspseudo = .1;                  // pseudocount for patterns in donor splice sites 
Double         IntronModel::dssneighborfactor = (float) 0.01; 
SnippetProbs*  IntronModel::snippetProbs;
SnippetProbs*  IntronModel::rSnippetProbs;
bool           IntronModel::initAlgorithmsCalled = false;
bool           IntronModel::haveSnippetProbs = false;
int            IntronModel::lastParIndex = -1;
Integer        IntronModel::verbosity;
BinnedMMGroup  IntronModel::dssBinProbs("dss", 1); // features from binned dss probs (for CRF training), monotonic increasing
BinnedMMGroup  IntronModel::assBinProbs("ass", 1); // features from binned ass probs (for CRF training), monotonic increasing
int            IntronModel::beginOfBioIntron = 0;
int            IntronModel::endOfBioIntron = 0;
int            IntronModel::ass_outside = 0;
/*
 * IntronModel constructor
 */
IntronModel::IntronModel() : gweight (1) {
    itype = toStateType( Properties::getProperty("/IntronModel/type", introncount++) );
    codon = new char[4];
    codon[3] = '\0';
}

/*
 * IntronModel destructor
 */
IntronModel::~IntronModel( ){
    if( --introncount == 0 ) {
      lastParIndex = -1;
      //delete assMotif;
      if (GCassMotif)
	delete [] GCassMotif;
      if (snippetProbs)
	delete snippetProbs;
      if (rSnippetProbs)
	delete rSnippetProbs;
      assMotif = GCassMotif = NULL;
      snippetProbs = rSnippetProbs = NULL;
      delete [] codon;
    }
}

/*
 * IntronModel initialisation of class variables
 */
void IntronModel::init() {
    // Initialize class variables from properties values
    Properties::assignProperty("/IntronModel/verbosity", verbosity);
    Properties::assignProperty("/IntronModel/patpseudocount", patpseudo);
    Properties::assignProperty("/IntronModel/k", k);
    try{
	slope_of_bandwidth = Properties::getdoubleProperty( "/IntronModel/slope_of_bandwidth");
    }catch( ProjectError e) { 
	cerr << e.getMessage();
    }
    try{
	minwindowcount = Properties::getIntProperty( "/IntronModel/minwindowcount");
    }catch( ProjectError e) { 
	cerr << e.getMessage();
    }
    try{
	asspseudo = Properties::getDoubleProperty( "/IntronModel/asspseudocount" );
    }catch( ProjectError e) { 
	cerr << e.getMessage();
    }
    try{
	dsspseudo = Properties::getDoubleProperty( "/IntronModel/dsspseudocount" );
    }catch( ProjectError e) { 
	cerr << e.getMessage();
    }
    try{
	dssneighborfactor = Properties::getDoubleProperty( "/IntronModel/dssneighborfactor" );
    }catch( ProjectError e) { 
	cerr << e.getMessage();
    }
    try{
	d = Properties::getIntProperty( "/IntronModel/d" );
    }catch( ProjectError e) { 
	cerr << e.getMessage();
    }
    try{
	ass_motif_memory = Properties::getIntProperty( "/IntronModel/ass_motif_memory" );
    }catch( ProjectError e) { 
	cerr << e.getMessage();
    }
    try{
	ass_motif_radius = Properties::getIntProperty( "/IntronModel/ass_motif_radius" );
    }catch( ProjectError e) { 
	cerr << e.getMessage();
    }
    Properties::assignProperty("/IntronModel/non_gt_dss_prob", non_gt_dss_prob);
    Properties::assignProperty("/IntronModel/non_ag_ass_prob", non_ag_ass_prob);

    ass_upwindow_size = Constant::ass_upwindow_size;
    ass_outside = ass_upwindow_size + Constant::ass_start + ASS_MIDDLE;

    if (d < DSS_MIDDLE + Constant::dss_end + ass_upwindow_size + Constant::ass_start + ASS_MIDDLE) {
	throw IntronModelError("Inconsistent intron length parameters. Please increase /IntronModel/d or decrease /IntronModel/ass_motif_memory.");
    }
    // reserve space for GC content dependent arrays
    if (!GCprobShortIntron)
      GCprobShortIntron = new Double[Constant::decomp_num_steps];
    if (!GCmal)
      GCmal = new Double[Constant::decomp_num_steps];
    if (!GCassMotif)
      GCassMotif = new Motif[Constant::decomp_num_steps];
    if (!GCemiprobs)
      GCemiprobs = new PatMMGroup[Constant::decomp_num_steps];
}



/*
 * IntronModel::readProbabilities
 */
void IntronModel::readProbabilities( int parIndex ){

    if (introncount == 0 || parIndex == lastParIndex)
	return;
    
    string filename = 
	Constant::fullSpeciesPath() +
	Properties::getProperty("/IntronModel/infile");
    ifstream istrm(filename.c_str());
    if( istrm ){
	int size;
	Seq2Int s2i_a( Constant::ass_size() );
	Seq2Int s2i_d( Constant::dss_size() );
	Double dbl;

	if (!hasSpliceSites) {
	    istrm >> goto_line_after( "[ASS]" ); 
	    istrm >> comment >> size >> comment >> c_ass >> comment >> asspseudo; 
	    
	    assprobs.assign( size, Double(asspseudo)/Double(c_ass + asspseudo * size) );
	    //    cerr << "\t c_ass = " << Double(emipseudo)/Double(c_ass+emipseudo*size) << endl;
	    while( istrm >> comment >> ws, istrm && istrm.peek() != '[' ){
		int pn = s2i_a.read(istrm);
		istrm >> dbl >> comment;
		assprobs[pn] = dbl / 1000;
	    }
	    streampos spos = istrm.tellg();
            istrm >> goto_line_after( "[ASSBIN]" );
            if (!istrm) {
                istrm.clear();
                istrm.seekg(spos); // go back to where you were
		Constant::ass_maxbinsize = 0; // no binning at all
            } else {
                assBinProbs.read(istrm);
                //cout << "number of bins: " << assBinProbs.nbins << endl;
                //assBinProbs.printBoundaries();
            }

	    istrm >> goto_line_after( "[DSS]" );
	    istrm >> comment >> size >> comment >> c_dss >> comment >> dsspseudo; 	    
	    dssprobs.assign( size, 0);
	    istrm >> comment >> ws;
	    for (int pn=0; pn<size; pn++){
		int dummy_pn = s2i_d.read(istrm);
		if (pn != dummy_pn) 
		    throw ProjectError("IntronModel::readProbabilities:  Error reading file " + filename);
		istrm >> dbl >> comment;
		dbl = dbl / 1000;
		dssprobs[pn] = dbl;
	    }
	    spos = istrm.tellg();
	    istrm >> goto_line_after( "[DSSBIN]" );
	    if (!istrm) {
		istrm.clear();
		istrm.seekg(spos); // go back to where you were
		Constant::dss_maxbinsize = 0; // no binning at all
	    } else {
		dssBinProbs.read(istrm);
		//cout << "number of bins: " << dssBinProbs.nbins << endl;
		//dssBinProbs.printBoundaries();
	    }

	    hasSpliceSites = true;
	    
	    // Die Längenwahrscheinlichkeiten einlesen
	    //------------------------------------------
	    istrm >> goto_line_after( "[LENGTH]" );
	    istrm >> comment >> d;
	    lenDist.assign(d+1, 0);
	    for( int i = 0; i < lenDist.size(); i++ ){
		istrm >> comment >> dbl;
		lenDist[i] = dbl / 1000;
	    }
	}

	char zusString[6];
	sprintf(zusString, "[%d]", parIndex);
	istrm >> goto_line_after(zusString);
	
	// Read in transition probabilities 
	// ---------------------------------------
	istrm >> goto_line_after( "[TRANSITION]" );
	istrm >> comment >> probShortIntron;
	istrm >> comment >> mal;


        // read in the emission probabilities
        //--------------------------------------------
        istrm >> goto_line_after( "[EMISSION]" );
        istrm >> comment >> size; 	
	istrm >> comment >> k;
	istrm >> comment >> patpseudo >> comment;
	Seq2Int s2i(k+1);
	emiprobs.probs.resize(size);
	emiprobs.order = k;
	for (int i=0; i< size; i++) {
            istrm >> comment; // comment is needed here
	    int pn = s2i.read(istrm);
            istrm >> emiprobs.probs[pn];
        }


	//delete assMotif; // remove this because the real memory is reserved with the array now
	
	assMotif = new Motif();
	istrm >> goto_line_after( "[ASSMOTIF]" );
	assMotif->read(istrm);

        istrm.close();
	lastParIndex = parIndex;
    } else 
	throw ProjectError("IntronModel::readProbabilites: Couldn't open file " + filename);
    // Nachstehendes unsinniges statement bewirkt, dass eine exception gefangen wird: ?!?
    if (false)
	cerr << "exiting IntronModel::readProbabilities... " << endl << flush;
}

/*
 * readAllParameters
 * read in species_intron_probs.pbl parameter files and store parameters for all gc content classes
 */
void IntronModel::readAllParameters(){
    if (introncount == 0)
	return;
    
    string filename = Constant::fullSpeciesPath() + Properties::getProperty("/IntronModel/infile");
    ifstream istrm(filename.c_str());
    if( istrm ){
	int size;
	Seq2Int s2i_a( Constant::ass_size() );
	Seq2Int s2i_d( Constant::dss_size() );
	Double dbl;

	if (!hasSpliceSites) {
	    istrm >> goto_line_after( "[ASS]" ); 
	    istrm >> comment >> size >> comment >> c_ass >> comment >> asspseudo; 
	    
	    assprobs.assign( size, Double(asspseudo)/Double(c_ass + asspseudo * size) );
	    //    cerr << "\t c_ass = " << Double(emipseudo)/Double(c_ass+emipseudo*size) << endl;
	    while( istrm >> comment >> ws, istrm && istrm.peek() != '[' ){
		int pn = s2i_a.read(istrm);
		istrm >> dbl >> comment;
		assprobs[pn] = dbl / 1000;
	    }
	    streampos spos = istrm.tellg();
            istrm >> goto_line_after( "[ASSBIN]" );
            if (!istrm) {
                istrm.clear();
                istrm.seekg(spos); // go back to where you were
		Constant::ass_maxbinsize = 0; // no binning at all
            } else {
                assBinProbs.read(istrm);
                //cout << "number of bins: " << assBinProbs.nbins << endl;
                //assBinProbs.printBoundaries();
            }

	    istrm >> goto_line_after( "[DSS]" );
	    istrm >> comment >> size >> comment >> c_dss >> comment >> dsspseudo; 	    
	    dssprobs.assign( size, 0);
	    istrm >> comment >> ws;
	    for (int pn=0; pn<size; pn++){
		int dummy_pn = s2i_d.read(istrm);
		if (pn != dummy_pn) 
		    throw ProjectError("IntronModel::readProbabilities:  Error reading file " + filename);
		istrm >> dbl >> comment;
		dbl = dbl / 1000;
		dssprobs[pn] = dbl;
	    }
	    spos = istrm.tellg();
	    istrm >> goto_line_after( "[DSSBIN]" );
	    if (!istrm) {
		istrm.clear();
		istrm.seekg(spos); // go back to where you were
		Constant::dss_maxbinsize = 0; // no binning at all
	    } else {
		dssBinProbs.read(istrm);
		//cout << "number of bins: " << dssBinProbs.nbins << endl;
		//dssBinProbs.printBoundaries();
	    }

	    hasSpliceSites = true;
	    
	    // read in length distribution
	    //------------------------------------------
	    istrm >> goto_line_after( "[LENGTH]" );
	    istrm >> comment >> d;
	    lenDist.assign(d+1, 0);
	    for( int i = 0; i < lenDist.size(); i++ ){
		istrm >> comment >> dbl;
		lenDist[i] = dbl / 1000;
	    }
	}
	
	// Start of GC content dependent parameters section
	// ------------------------------------------------
	char zusString[6];
	/*	if (GCprobShortIntron) 
	  delete [] GCprobShortIntron;
	if (GCmal)
	  delete [] GCmal;
	if (GCassMotif)
	  delete [] GCassMotif;
	if (GCemiprobs)
	delete [] GCemiprobs;*/

	for (int idx = 0; idx < Constant::decomp_num_steps; idx++) {
	  sprintf(zusString, "[%d]", idx+1);
	  istrm >> goto_line_after(zusString);
	  
	  // read in transition probability parameters
	  istrm >> goto_line_after( "[TRANSITION]" );
	  istrm >> comment >> GCprobShortIntron[idx];
	  istrm >> comment >> GCmal[idx];
	  
	  
	  // read in the emission probabilities
	  GCemiprobs[idx].setName(string("intron emiprob gc") + itoa(idx+1));
	  istrm >> goto_line_after( "[EMISSION]" );
	  istrm >> comment >> size; 	
	  istrm >> comment >> k;
	  istrm >> comment >> patpseudo >> comment;
	  Seq2Int s2i(k+1);
	  GCemiprobs[idx].probs.resize(size);
	  GCemiprobs[idx].order = k;
	  for (int i=0; i< size; i++) {
            istrm >> comment; // comment is needed here
	    int pn = s2i.read(istrm);
            istrm >> GCemiprobs[idx].probs[pn];
	    if (!Constant::contentmodels)
		GCemiprobs[idx].probs[pn] = 0.25; // use uniform distribution
	  }
	  
	  istrm >> goto_line_after( "[ASSMOTIF]" );
	  GCassMotif[idx].read(istrm);
	}
	// this default setting is only relevant for exon candidate filtering when doing comparative gene
	// finding
	assMotif = &GCassMotif[(int) (Constant::decomp_num_steps/2)];
	istrm.close();
    } else 
	throw ProjectError("IntronModel::readAllParameters: Couldn't open file " + filename);
}

/*
 * IntronModel::initSnippetProbs
 */
void IntronModel::initSnippetProbs() {
    if (snippetProbs)
	delete snippetProbs;
    snippetProbs = new SnippetProbs(sequence, k);
    if (rSnippetProbs)
	delete rSnippetProbs;
    rSnippetProbs = new SnippetProbs(sequence, k, false);
    haveSnippetProbs = true;
}

/*
 * IntronModel::initAlgorithms
 *
 * makes a correction on the transition matrix "trans" and on the vector of ancestors
 */
void IntronModel::initAlgorithms( Matrix<Double>& trans, int cur){
}

// call state-dependent (cur) stuff after gcIdx was set correctly
void IntronModel::updateToLocalGCEach( Matrix<Double>& trans, int cur ) {
    if (!initAlgorithmsCalled) {
	seqProb(-1, 0);
	aSSProb(-1, false); // initialize ass probs
	dSSProb(-1, false); // initialize dss probs
	initAlgorithmsCalled = true;
    }
    /*
     * correction of the transition probabilities into lessDx and equalx (x=0,1,2)
     * which are considered parameters of the intron model
     * böser HACK: we assume that the transition matrix has for each frame
     * and each state i only nonzero transitions in both lessD and equalD or in none of them
     *
     */
    Double factor = 0;
    if (itype == lessD0 || itype == lessD1 || itype == lessD2 || 
	itype == rlessD0 || itype == rlessD1 || itype == rlessD2)
	factor = probShortIntron;
    else if (itype == equalD0  || itype == equalD1 || itype == equalD2 || 
	     itype == requalD0  || itype == requalD1 || itype == requalD2)
	factor = (1 - probShortIntron);
    if (factor > 0) 
	for (int i = 0; i < trans.getColSize(); i++) 
	    if (trans[i][cur] > 0)
		trans[i][cur] = factor;
    
    /*
     * the probability of returning to a "geometric" intron state is reset
     * according to the mal (mean additional length) value
     */
    if (itype == geometric0 || itype == geometric1 || itype == geometric2 || 
	itype == rgeometric0 || itype == rgeometric1 || itype == rgeometric2) {
        if (mal > 0.0){
	    trans[cur][cur] = 1 - 1 / mal;
  	    geoProb = trans[cur][cur].doubleValue();
	}
	// normalize the rest of the line
	Double sum = 0;
	for (int i = 0; i < trans.getRowSize(); i++) {
	    if (i != cur)
		sum += trans[cur][i];
	}
	if (sum > 0) {
	    for (int i = 0; i < trans.getRowSize(); i++) {
		if (i != cur)
		    trans[cur][i] /= mal * sum;
	    }
	}
    }
}

/*
 * IntronModel::initAlgorithms
 *
 * makes a correction on the transition matrix "trans" and on the vector of ancestors
 */
void IntronModel::updateToLocalGC(int from, int to){
    // set these parameters to the one of the GC content index
    emiprobs = GCemiprobs[gcIdx];
    assMotif = &GCassMotif[gcIdx];
    probShortIntron = GCprobShortIntron[gcIdx];
    mal = GCmal[gcIdx];
    snippetProbs->setEmiProbs(&emiprobs.probs);
    rSnippetProbs->setEmiProbs(&emiprobs.probs);
}


/*
 * IntronModel::viterbiForwardAndSampling
 */
void IntronModel::viterbiForwardAndSampling(ViterbiMatrixType& viterbi,
					    ViterbiMatrixType& forward,
					    int state,
					    int base,
					    AlgorithmVariant algovar,
					    OptionListItem& oli) {
    /*
     * maximal length of intron state with specific length distribution
     * see also seqProb!!
     */
    int dStateLen = d - DSS_MIDDLE - Constant::dss_end - Constant::ass_start 
	- ASS_MIDDLE - ass_upwindow_size;
    Double predProb;
    Double fwdsum(0); 
    Double emiProb, extrinsicQuot(1);
    Double transEmiProb;
    vector<Ancestor>::const_iterator it;
    int endOfPred;
    Double maxPredProb = 0;
    SubstateId substate;
    ViterbiSubmapType substates(state); // do not activate yet
    bool checkSubstates = false;
    if (profileModel) {
	checkSubstates = doViterbi(algovar);
	if (algovar == doBacktracking)
	    getStatePair(state, state, substate);
    }
    OptionsList *optionslist = NULL;
    if (algovar == doSampling)
	optionslist = new OptionsList();

    if (itype == lessD0 || itype == lessD1 || itype == lessD2 ||
	itype == rlessD0 || itype == rlessD1 || itype == rlessD2) {
	/*
	 * State with variable length
	 */
	if (isOnFStrand(itype))
	    endOfBioIntron = base + ass_upwindow_size + Constant::ass_start + ASS_MIDDLE;
	else
	    endOfBioIntron = base + Constant::dss_end + DSS_MIDDLE;

	/*
	 * For performance: check whether an acceptor/donor splice site could follow
	 */
	if ((isOnFStrand(itype) && (endOfBioIntron - ASS_MIDDLE + 1 < dnalen - 1 && !isPossibleASS(endOfBioIntron))) ||
	    (!isOnFStrand(itype) && (endOfBioIntron - DSS_MIDDLE + 1 < dnalen - 1 && !isPossibleRDSS(endOfBioIntron)))){
	    return;
	}
	if (itype != lessD0 && itype != rlessD2 && endOfBioIntron < dnalen - 2){ // a codon is spliced, fill bases in second exon
	   switch (itype)
	      {
	      case lessD1: // *|....|**
		 codon[1] = *(sequence + endOfBioIntron + 1);
		 codon[2] = *(sequence + endOfBioIntron + 2);
		 break;
	      case lessD2: // **|....|*
		 codon[2] = *(sequence + endOfBioIntron + 1);
		 break;
	      case rlessD0: // **|....|*
		 codon[0] = wcComplement(*(sequence + endOfBioIntron + 1));
		 break;
	      case rlessD1: // *|....|**
		 codon[0] = wcComplement(*(sequence + endOfBioIntron + 2));
		 codon[1] = wcComplement(*(sequence + endOfBioIntron + 1));
		 break;
	      default:;
	      }
	}
	/*
	 * maximal length of intron state with specific length distribution
	 * see also seqProb!!
	 */	    
	int leftMostEndOfPred = base - dStateLen;
	if (leftMostEndOfPred < 0)
	    leftMostEndOfPred = 0;
	   
	for (endOfPred = base-1; endOfPred >= leftMostEndOfPred; endOfPred--){
	    ViterbiColumnType& predVit = algovar == doSampling ? forward[endOfPred] : viterbi[endOfPred];
	    // compute the maximum over the predecessor states probs times transition probability
	    /*
	     * check whether the starting position has a positive entry at all
	     */
	    for(it = ancestor.begin(); 
		it != ancestor.end() && !predVit.has(it->pos);
		++it);
	    if (it == ancestor.end()) 
		continue;
	    emiProb = emiProbUnderModel(endOfPred+1, base);
	    // every intron gets a malus, the short introns here, the long introns in the emission of the equalD state
	    // when the intron is supported, which is checked below, then only the bonus applies.
	    if (seqFeatColl)
		extrinsicQuot = seqFeatColl->collection->malus(intronF);
	    do {
		if (!predVit.has(it->pos)) continue;
		transEmiProb = it->val * emiProb * extrinsicQuot;
		predProb = predVit.get(it->pos, substate) * transEmiProb;
		if (needForwardTable(algovar)) {
		    Double fwdsummand = forward[endOfPred].get(it->pos) * transEmiProb.heated();
		    if (algovar == doSampling) {
			if (fwdsummand > 0)
			    optionslist->add(it->pos, endOfPred, fwdsummand);
			continue;
		    }
		    fwdsum += fwdsummand;
		}
		if (checkSubstates && predVit.hasSubstates(it->pos)) {
		    if (substates.inactive())
			substates.activate();
		    substates.merge(predVit.getSubstates(it->pos), transEmiProb);
		}
		if (predProb > maxPredProb) {
		    maxPredProb = predProb;
		    oli.state = it->pos;
		    oli.base = endOfPred;
		}
	    } while (++it != ancestor.end());
	}

	/*
	 * Check whether there are intron hints ending at the end of the intron.
	 * If yes, check for each hint whether that supported option is more likely 
	 * than the current most likely predecessor option.
	 */
	if (seqFeatColl){
	    Feature *intronList = seqFeatColl->getFeatureListAt(intronF, endOfBioIntron, 
								isOnFStrand(itype)? plusstrand : minusstrand);
	   
	    int oldEndOfPred = -INT_MAX;
 	    for (Feature *ihint = intronList; ihint != NULL; ihint = ihint->next) {
		if (isOnFStrand(itype))
		    endOfPred = ihint->start - 1 + DSS_MIDDLE + Constant::dss_end;
		else 
		    endOfPred = ihint->start - 1 + ass_upwindow_size + Constant::ass_start + ASS_MIDDLE;
		if (endOfPred >= 0 && ihint->end - ihint->start + 1 >= ass_upwindow_size + Constant::ass_start + ASS_MIDDLE + DSS_MIDDLE + Constant::dss_end) {
		    ViterbiColumnType& predVit = viterbi[endOfPred];
		    emiProb = emiProbUnderModel(endOfPred+1, base);
		    if (endOfPred == oldEndOfPred)
			extrinsicQuot *= ihint->bonus;
		    else 
			extrinsicQuot = ihint->bonus;
		    //cout << "IntronModel extrQ=" << extrinsicQuot << " " << ihint->start << ".." << ihint->end << endl;
		    emiProb *= extrinsicQuot; // option gets a bonus because it complies with the intron hint
		    for( it = ancestor.begin(); it != ancestor.end(); ++it ){
			transEmiProb = it->val * emiProb;
			predProb = predVit.get(it->pos, substate) * transEmiProb;
			if (needForwardTable(algovar)) {
			    // ACHTUNG: this isn't correct in the model
			    // but the effect should be very small
			    Double fwdsummand = forward[endOfPred].get(it->pos) * transEmiProb.heated();
			    if (algovar == doSampling) {
				if (fwdsummand > 0)
				    optionslist->add(it->pos, endOfPred, fwdsummand);
				continue;
			    }
			    fwdsum += fwdsummand;
			}
			if (checkSubstates && predVit.hasSubstates(it->pos)) {
			    if (substates.inactive())
				substates.activate();
			    substates.merge(predVit.getSubstates(it->pos), transEmiProb);
			} 
			//cout << "checking whether long intron shortcut is better..." << predProb << " " << maxPredProb << "\n";
			if (predProb > maxPredProb) {
			    maxPredProb = predProb;
			    oli.state = it->pos;
			    oli.base = endOfPred;
			}
		    }
		    oldEndOfPred = endOfPred;
		}
	    }
	}
    } else {
	/*
	 * State with fixed length
	 */
	bool abort = false;
	switch( itype ){
	    case longdss0: case longdss1: case longdss2:
		endOfPred = base - Constant::dss_whole_size();
		abort = (endOfPred < 0) || !isPossibleDSS(base - Constant::dss_end - DSS_MIDDLE + 1);
		break;
	    case equalD0: case equalD1: case equalD2: case requalD0: case requalD1: case requalD2:
		endOfPred = base - dStateLen;
		abort = (endOfPred < 0);
		break;
	    case geometric0: case geometric1: case geometric2: case rgeometric0: case rgeometric1: case rgeometric2:
		endOfPred = base - 1;
		// compute the maximum over the (two) predecessor states probs times transition probability
		break;
	    case longass0: case longass1: case longass2:
		endOfPred = base - Constant::ass_whole_size() - ass_upwindow_size;
		abort = (endOfPred < 0 || !isPossibleASS(base - Constant::ass_end));
		break;
	    case rlongdss0: case rlongdss1: case rlongdss2:
		endOfPred = base - Constant::dss_whole_size();
		abort = ((endOfPred < 0) || !isPossibleRDSS(base - Constant::dss_start));
		break;
	    case rlongass0: case rlongass1: case rlongass2:
		endOfPred = base - Constant::ass_whole_size() - ass_upwindow_size;
		abort = (endOfPred < 0 || !isPossibleRASS(base - ass_upwindow_size - Constant::ass_start - ASS_MIDDLE + 1));
		break;
	    default:
		throw IntronModelError("Unknown intron type.");
	}

	/*
	 * return immediately if coordinates are out of range
	 */
	if (abort) {
// 	    viterbi[base].erase(state);
// 	    if (needForwardTable(algovar)) 
// 		forward[base].erase(state);
	    return;
	}

	/*
	 * check whether the starting position has a positive entry at all
	 */
	ViterbiColumnType& predVar = algovar==doSampling? 
	    forward[endOfPred] : viterbi[endOfPred];
	for(it = ancestor.begin(); 
	    it != ancestor.end() && !predVar.has(it->pos);
	    ++it);

	/*
	 * return immediately if no valid ancestor exists
	 */
	if (it == ancestor.end()) {
// 	    viterbi[base].erase(state);
// 	    if (needForwardTable(algovar)) 
// 		forward[base].erase(state);
	    return;
	}
	
	emiProb = emiProbUnderModel(endOfPred+1, base);
	if (emiProb == 0) 
	    return;
	oli.base = endOfPred;

	// compute the maximum over the predecessor states probs times transition probability
// 	if (substate >= 0 && (itype == longass1 || itype == longass2))
// 	    substate = profileModel->advanceId(substate, -1);
	ViterbiSubmapType* firstPredSubstates = 0;
	do {
	    transEmiProb = it->val * emiProb; 
	    predProb = predVar.get(it->pos, substate) * transEmiProb;
	    if (needForwardTable(algovar)) {
		Double fwdsummand = forward[endOfPred].get(it->pos) * transEmiProb.heated();
		if (algovar == doSampling) {
		    if (fwdsummand > 0) 
			optionslist->add(it->pos, endOfPred, fwdsummand);
		    continue;
		}
		fwdsum += fwdsummand;
	    }
	    if (checkSubstates) 
		try {
		    if (substates.empty()) {
			// link it to predecessor map
			firstPredSubstates = &predVar.getSubstates(it->pos);
			substates.link(*firstPredSubstates, transEmiProb);
			substates.cloneMapIfLinkcountExceeded(firstPredSubstates);
		    } else {
			ViterbiSubmapType& predSubstates = predVar.getSubstates(it->pos);
			substates.cloneMap(firstPredSubstates);
			substates.merge(predSubstates, transEmiProb);
		    }
		} catch (NoSubmapFoundError) {}
	    if (predProb > maxPredProb) {
		maxPredProb = predProb;
		oli.state = it->pos;
	    }
	} while (++it != ancestor.end());

	if (checkSubstates && !substates.empty() && (isLongAssIntron(itype) || isRLongDssIntron(itype))) {
	    if (substates.is_linked())
		substates.cloneMap(firstPredSubstates);
	    int firstBase;
	    switch (itype) {
		case longass0:
		    firstBase = base - Constant::ass_end + 1; break;
		case longass1: 
		    firstBase = base - Constant::ass_end + 3; break;
		case longass2:
		    firstBase = base - Constant::ass_end + 2; break;
		case rlongdss0:
		    firstBase = base - Constant::dss_start +2 ; break;
		case rlongdss1:
		    firstBase = base - Constant::dss_start +3 ; break;
		case rlongdss2:
		    firstBase = base - Constant::dss_start +1 ; break;
		default:
		    throw ProjectError("Internal Error (bad state) in IntronModel::viterbi");
	    }
	    profileModel->scoreAssOrRDss(substates, !isOnFStrand(itype), firstBase, maxPredProb);
	}
    } // state with fixed length


    switch (algovar) {
	case doViterbiAndForward:
	    if (fwdsum > 0)
		forward[base][state] = fwdsum;
	case doViterbiOnly:
	    if (maxPredProb > 0)
		viterbi[base][state] = maxPredProb;
	    if (!substates.empty()) {
		if (!substates.is_linked())
		    for (ViterbiSubmapType::iterator it = substates.begin();
			 it != substates.end();
			 ++it) {
			ViterbiSubmapType* pred = it->second.predMap;
			if (pred) {
#ifdef DEBUG
			    if (it->first.slot >= MAX_BLOCKCOUNT)
				throw ProjectError("Intronmodel::viterbi: bad range for substate");
			    if (pred->get(it->first) == 0)
				throw ProjectError("Intronmodel::viterbi: substate with bad predecessor");
#endif
			    pred->succCount(it->first)++;
			}
#ifdef DEBUG		
			else if (it->first.slot < MAX_BLOCKCOUNT) 
			    throw ProjectError("Intronmodel::viterbi: substate w/o predecessor");
#endif

		    }
		viterbi[base].addSubstates(substates);
	    } 
	    return;
	case doSampling:
	    optionslist->prepareSampling();
	    try {
		oli = optionslist->sample();
	    } catch (ProjectError e) {
		cerr << "Sampling error in intron model. state=" << state << " base=" << base << endl;
		throw e;
	    }
	    delete optionslist;
	    return;
	case doBacktracking:
	    oli.state = getFullStateId(oli.state, substate);
	    return; // backtracking: do nothing
    } // end switch
} // end IntronModel::viterbiForwardAndSampling


Double IntronModel::emiProbUnderModel (int begin, int end) const {
    Seq2Int s2i(k+1);
    static Double returnProb, extrinsicQuot, restSeqProb, lenPartProb;
    returnProb = extrinsicQuot = 1;
    if (inCRFTraining){
	seqProb(-1, -1); // forget all saved information in static variables
	if (itype == longass0 || itype == longass1 || itype == longass2 || itype == rlongass0 || itype == rlongass1 || itype == rlongass2){
	    aSSProb(-1, false);
	    dSSProb(-1, false);
	}
    }
    // intron range for extrinsic bonus: intronic positions with (begin,end) boundaries
    int intronBegin=begin, intronEnd=end;
    switch( itype ){
	case longdss0: case longdss1: case longdss2: 
	    if (!(end-begin+1 == Constant::dss_whole_size())) {
		cerr << begin << ".." << end << " " << Constant::dss_whole_size() << endl;
		throw ProjectError("Assertion failed in IntronModel::emiProbUnderModel. DSS");
	    }
	    returnProb = dSSProb(end - Constant::dss_whole_size() + 1, true);
	    intronBegin = end - DSS_MIDDLE - Constant::dss_end + 1;
	    break;
	case rlongdss0: case rlongdss1: case rlongdss2:
	    if (!(end-begin+1 == Constant::dss_whole_size())) 
		throw ProjectError("Assertion failed in IntronModel::emiProbUnderModel. rDSS");
	    returnProb = dSSProb(end - Constant::dss_whole_size() + 1, false);
	    intronEnd = end - Constant::dss_start;
	    break;
	case equalD0: case equalD1: case equalD2: case requalD0: case requalD1: case requalD2:
	    returnProb = seqProb(begin, end);
	    // long unsupported introns must be punished
	    if (seqFeatColl)
		extrinsicQuot *= seqFeatColl->collection->malus(intronF);
	    break;
	case geometric0: case geometric1: case geometric2: case rgeometric0: case rgeometric1: case rgeometric2:
	    /*
	     * compute the pure emission probabilities of the intron model 
	     * of the sequence from begin to end. in viterbi: begin=end
	     */
	    returnProb = 1;
	    for(; begin<=end; begin++) {
		if (begin >= k) {
		    try {
			int pn = s2i(sequence + begin - k);
			returnProb *= emiprobs.probs[pn];
			if (inCRFTraining && (countEnd < 0 || (begin >= countStart && begin <= countEnd)))
			    GCemiprobs[gcIdx].addCount(pn);
		    } catch (InvalidNucleotideError e) {
		       returnProb *= (float) 0.25;
		    }
		} else { // at the very beginning of the sequence
		   returnProb *= (float) 0.25;
		}
	    }
	    break;
	case longass0: case longass1: case longass2:
	    returnProb = aSSProb(end - Constant::ass_whole_size() - ass_upwindow_size + 1, true);
	    intronEnd = end - Constant::ass_end;
	    break;
	case rlongass0: case rlongass1: case rlongass2:
	    returnProb = aSSProb(end - Constant::ass_whole_size() - ass_upwindow_size + 1, false);
	    intronBegin = begin + Constant::ass_end;
	    break;
	case lessD0: case lessD1: case lessD2: case rlessD0: case rlessD1: case rlessD2:
	   { 
	      if (isOnFStrand(itype)) {
		 beginOfBioIntron = begin - Constant::dss_end - DSS_MIDDLE;
		 if (beginOfBioIntron >= 0 && !isPossibleDSS(beginOfBioIntron))
		    return 0;
	      } else {
		 beginOfBioIntron = begin - ass_outside; // ass_outside = ass_upwindow_size + Constant::ass_start + ASS_MIDDLE;
		 if (beginOfBioIntron >=0 && !isPossibleRASS(beginOfBioIntron))
		    return 0;
	      }
	      if (itype != lessD0 && itype != rlessD2 && beginOfBioIntron > 1){ // a codon is spliced, check whether it is a stop codon
		 switch (itype)
		    {
		    case lessD1: // *|....|**
		       codon[0] = *(sequence + beginOfBioIntron - 1);
		       break;
		    case lessD2: // **|....|*
		       codon[0] = *(sequence + beginOfBioIntron - 2);
		       codon[1] = *(sequence + beginOfBioIntron - 1);
		       break;
		    case rlessD0: // **|....|*
		       codon[1] = wcComplement(*(sequence + beginOfBioIntron - 1));
		       codon[2] = wcComplement(*(sequence + beginOfBioIntron - 2));
		       break;
		    case rlessD1: // *|....|**
		       codon[2] = wcComplement(*(sequence + beginOfBioIntron - 1));
		       break;
		    default:;
		    }
		 if (GeneticCode::isStopcodon(codon)) // prevent in-frame stop codons spliced by SHORT or HINTED introns (not other long ones)
		    return 0;
	      }
	    int intronLength = endOfBioIntron - beginOfBioIntron + 1;
	    
	    if (intronLength <= d){
	       restSeqProb = seqProb(begin, end);
	       lenPartProb = lenDist[intronLength];
	    } else {
	       // this case applies only for long introns supported by hints as other long introns
	       // are handled by a geometric state
	       // Use gc index as in the normal case
	       // |---- d ------|-----------------geo---------------------
	       //    use this:  ^    use local gc index here
	       lenPartProb = lenDist[d] * pow(1.0 - 1.0/mal.doubleValue(), (int) (intronLength-d));
	       restSeqProb = 1;
	       Seq2Int s2i(k+1);
	       int idx = getGCIdx((begin+d < dnalen)? begin+d : dnalen-1);
	       for (int a=begin; a < begin+d; a++)
		  try {
		     restSeqProb *= (a >= k && a < dnalen - 1)? GCemiprobs[idx].probs[s2i(sequence + a - k)] : 0.25;
		  } catch(...) {
		     restSeqProb *= .25;
		  }
	       for (int a=begin+d; a <=end; a++){
		  if (idx != getGCIdx(a))
		     idx = getGCIdx(a);
		  try {
		     restSeqProb *= (a >= k && a < dnalen - 1)? GCemiprobs[idx].probs[s2i(sequence + a - k)] : 0.25;
		  } catch(...) {
		     restSeqProb *= (float) .25;
		  }
	       }
	    }
	    returnProb = lenPartProb * restSeqProb;
	    break;
	}
	default: 
	    throw IntronModelError("Unknown intron type.");
    }
    /*
     * compute extrinsicQuot (bonus for hints)
     * intronpart hint:
     * an intron gets the bonus for each position covered by an intronpart hint 
     * (counted multiply for overlapping hints)
     */
    if (seqFeatColl){
	if (intronBegin < 0) 
	    intronBegin=0;
	int coveredBegin = intronEnd+1;
	int coveredEnd = intronBegin-1;
	Feature *part;
	Feature *intronList = seqFeatColl->getFeatureListOvlpingRange(power2id[intronpartF] | power2id[nonexonpartF], intronBegin, 
								      intronEnd, 
								      isOnFStrand(itype)? plusstrand : minusstrand);
	for (part = intronList; part!= NULL; part = part->next){
	    if (part->start< coveredBegin)
		coveredBegin = part->start;
	    if (part->end > coveredEnd)
		coveredEnd = part->end;
	    for (int i=intronBegin; i<=intronEnd; i++) {
		if (part->start<=i && part->end>=i)
		    extrinsicQuot *= part->bonus;
	    }
	}
	if (coveredEnd > intronEnd)
	    coveredEnd = intronEnd;
	if (coveredBegin < intronBegin)
	    coveredBegin = intronBegin;
	if (seqFeatColl->collection->hasHintsFile)
	  extrinsicQuot *= pow (seqFeatColl->collection->malus(intronpartF), intronEnd-coveredEnd+coveredBegin-intronBegin);
    }
    return returnProb * extrinsicQuot;
}


/*
 * computes the probability of the emission of the sequence from left to right
 * left and right included
 */

Double IntronModel::seqProb(int left, int right) const {
    static Double seqProb = 1;
    static int oldleft=-1, oldright=-1;
    static vector<Double> seqProbs;
    
    int curpos;
    Seq2Int s2i(k+1);
    if (left < 0) {   // new initialisation
	seqProb = 1;
	oldleft = oldright = -1;
	seqProbs.assign(d - DSS_MIDDLE - Constant::dss_end - Constant::ass_start 
			- ASS_MIDDLE - ass_upwindow_size + 1, 0);
	return 1;
    }
    if (left > right)
	return 1;
    
    if (!inCRFTraining){
	if (itype == lessD0 || itype == lessD1 || itype == lessD2 /*|| itype == equalD0 || itype == equalD1 || itype == equalD2*/)
	    return snippetProbs->getSeqProb(right, right-left+1);
	if (itype == rlessD0 || itype == rlessD1 || itype == rlessD2 /*|| itype == requalD0 || itype == requalD1 || itype == requalD2*/)
	    return rSnippetProbs->getSeqProb(right, right-left+1);
	
	if (right == oldright && left <= oldleft &&  right-left < seqProbs.size()) {  
	    for (curpos = oldleft-1; curpos >= left; curpos--){
		try {
		    seqProb *= emiprobs.probs[s2i(sequence + curpos - k)];
		} catch (InvalidNucleotideError e) {
		    seqProb *= (float) 0.25;
		}
		seqProbs[right-curpos] = seqProb;
	    }
	    oldleft = left;
	    return seqProb;
	}
	if (right == oldright && left > oldleft && right-left < seqProbs.size()) {
	    // retrieve from the stored table
	    return seqProbs[right - left];
	}
    }
    
    // compute everything new
    //cerr << "IntronModel:: seqProb: compute everything new: " << right - left << endl;
    seqProb = 1;
    for (curpos = right; curpos >= left; curpos--){
	try {
	    if (curpos - k >= 0){
		int pn = s2i(sequence + curpos - k);
		seqProb *= emiprobs.probs[pn];
		if (inCRFTraining && (countEnd < 0 || (curpos >= countStart && curpos <= countEnd)))
		    GCemiprobs[gcIdx].addCount(pn);
	    } else
		seqProb *= (float) 0.25;
	} catch (InvalidNucleotideError e) {
	    seqProb *= (float) 0.25;
	}
	if ( right-curpos < seqProbs.size())
	    seqProbs[right-curpos] = seqProb;
    }
    oldleft = left;
    oldright = right;
    return seqProb;
}

/*
 * base is the first position of the motif window (forward case)
 * intron ...mmmmmmmmmmmmmmmmmmmmm **AG|*** ... exon
 *            motif window         pattern
 * frequencies of (a,c,g,t) in the 30 positions before ass: (0.27, 0.211, 0.107, 0.412)
 */
Double IntronModel::aSSProb(int base, bool forwardStrand){
    static char *astr = new char[Constant::ass_size()+1]; // forget the delete
    static Seq2Int s2i(Constant::ass_size());
    static Double patternProb(1), motifProb(1), emiProb(0);
    static map<int, Double> memoF, memoR; // red-black trees to memoize values that have been computed before
    static map<int, Double>::const_iterator got;

    if (base < 0 || memoF.size() > 1000) {
        // reset static values
	// otherwise it may reuse values from the previous sequence
	memoF.clear();
        memoR.clear();
	if (base < 0)
	    return 0;
    }

    // return memoized values, if available
    if (forwardStrand && (got = memoF.find(base), got != memoF.end()))
	return got->second;
    if (!forwardStrand && (got = memoR.find(base), got != memoR.end()))
    	return got->second;
    
    bool nonAG;
    if (forwardStrand) {
	int asspos = base + ass_upwindow_size + Constant::ass_start;
	if (!isPossibleASS(asspos +1)) {
	    //memoF[base] = 0;
	    return 0;
	}
	nonAG = !onASS(sequence + asspos);
	strncpy(astr, sequence + base + ass_upwindow_size, Constant::ass_start);
	strncpy(astr + Constant::ass_start, sequence + asspos + ASS_MIDDLE, Constant::ass_end);
	// determine motifProb, the probability of the motif
	motifProb = (base >= assMotif->k)  ? 
	    assMotif->seqProb(sequence + base) : 0;
    } else {
	int asspos = base + Constant::ass_end;
	if (!isPossibleRASS(asspos)) {
	    //memoR[base] = 0;
	    return 0;
	}
	nonAG = !onRASS(sequence + asspos);
	putReverseComplement(astr, sequence + asspos + ASS_MIDDLE, Constant::ass_start);
	putReverseComplement(astr + Constant::ass_start, sequence + base, Constant::ass_end);
	int motifstart =  base + Constant::ass_whole_size();
	int motifend =  motifstart + ass_upwindow_size;
	motifProb = motifend + assMotif->k < dnalen ? 
	    assMotif->seqProb(sequence + motifstart, true, true) :
	    pow(.25, (int) ass_upwindow_size);
    }
    astr[Constant::ass_size()] = '\0';
    try {
	Double assprob = assprobs[ s2i(astr) ];
	if (nonAG) assprob *= non_ag_ass_prob;
        if (assBinProbs.nbins < 1) {
            patternProb = assprob; // standard HMM probabilities
	} else {
	    int idx = assBinProbs.getIndex(assprob);
	    if (inCRFTraining && (countEnd < 0 || (base >= countStart && base <= countEnd)))
		assBinProbs.addCount(idx);
	    patternProb = assBinProbs.avprobs[idx];
	}
    } catch (InvalidNucleotideError e) {
	patternProb = 0.001 * pow(.25, (int) Constant::ass_size());
    }

    emiProb = motifProb * patternProb;
    if (forwardStrand)
	memoF[base] = emiProb;	
    else
	memoR[base] = emiProb;
    return emiProb;
}


/*
 * base is the first position of the pattern
 * exon ...  ***|GT**** ...intron
 */
Double IntronModel::dSSProb(int base, bool forwardStrand){
    static char *astr = new char[Constant::dss_size()+1]; // forget the delete
    static Seq2Int s2i(Constant::dss_size());
    static map<int, Double> memoF, memoR; // red-black trees to memoize values that have been computed before
    static map<int, Double>::const_iterator got;
    
    if (base < 0 || memoF.size() > 500){
	memoF.clear(); // reset static values
	memoR.clear();
	if (base < 0)
	    return 0;
    }
    
    // return memoized values, if available
    if (forwardStrand && (got = memoF.find(base), got != memoF.end()))
	return got->second;
    if (!forwardStrand && (got = memoR.find(base), got != memoR.end()))
    	return got->second;

    bool nonGT;
    if (forwardStrand) { // forward strand
	int dsspos = base + Constant::dss_start;
	if (!isPossibleDSS(dsspos))
	    return 0;
	nonGT = !onDSS(sequence + dsspos);
	strncpy(astr, sequence + base, Constant::dss_start);
	strncpy(astr + Constant::dss_start, sequence + dsspos + DSS_MIDDLE, Constant::dss_end);
    } else { // reverse complement
	int dsspos = base + Constant::dss_end;
	if (!isPossibleRDSS(dsspos + 1))
	    return 0;
	nonGT = !onRDSS(sequence + dsspos);
	putReverseComplement(astr, sequence + dsspos + DSS_MIDDLE, Constant::dss_start);
	putReverseComplement(astr + Constant::dss_start, sequence + base, Constant::dss_end);
    }
    astr[Constant::dss_size()] = '\0';
    try {
	Double dssprob = dssprobs[ s2i(astr) ];
	if (nonGT) dssprob *= non_gt_dss_prob;
	if (dssBinProbs.nbins >= 1) {
	    int idx = dssBinProbs.getIndex(dssprob);
	    if (inCRFTraining && (countEnd < 0 || (base >= countStart && base <= countEnd)))
		dssBinProbs.addCount(idx);
	    dssprob = dssBinProbs.avprobs[idx];
	}
	if (forwardStrand)
	    memoF[base] = dssprob;
	else
	    memoR[base] = dssprob;
	return dssprob; // standard HMM probabilities
    } catch (InvalidNucleotideError e) {
	return 0; // don't predict splice site when there is an unknown nucleotide
    }
}

double IntronModel::getMeanIntrLen(){

    Double mil_long_intr = (mal + d);
    Double mil_short_intr = 0.0;
    for( int i = 0; i < lenDist.size(); i++ ){
	mil_short_intr += (lenDist[i]*i);
    }
    Double mil = (mil_long_intr * (1.0 - probShortIntron)) + (mil_short_intr * probShortIntron);
    return mil.doubleValue();
}

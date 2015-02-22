/**********************************************************************
 * file:    ncmodel.cc
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  noncoding model (all transcribed genes that are not translated to proteins)
 * authors: Mario Stanke (mario.stanke@uni-greisfwald.de)
 * 
 * date    |   author      |  changes 
 * --------|---------------|------------------------------------------ 
 * 21.02.15| Mario Stanke  | creation of the file
 **********************************************************************/

#include "ncmodel.hh"

// project includes
#include "gene.hh"
#include "projectio.hh"
#include "motif.hh"
#include "intronmodel.hh" // so splice sites models can be reused here
#include "utrmodel.hh"    // so TSS/TTS signal models can be reused here
#include "extrinsicinfo.hh"
#include "merkmal.hh"

#include <climits>

/*
 * Initialisation of static data members
 */
Integer         NcModel::nccount = 0;
vector<Double>  NcModel::lenDistInternal;     // length distribution of non-single exons
vector<Double>  NcModel::lenDistSingle;       // length distribution of single exons
Boolean         NcModel::hasLenDist = false;
Boolean         NcModel::hasLenDist = false;
Boolean         NcModel::initAlgorithmsCalled = false;
int             NcModel::boundSpacing = 10; // without hints 5' and 3' transcript end only every ttsSpacing bases, for speed

/*
 * NcModel constructor
 */
NcModel::NcModel() {
    nctype = toStateType( Properties::getProperty("/NcModel/type", utrcount++) );
}

/*
 * destructor
 */
NcModel::~NcModel( ){
}

/*
 * NcModel initialisation of class variables
 */
void NcModel::init() {
    try {
	if (!Properties::getBoolProperty(NONCODING_KEY)) // not run with --nc=On
	    return;
    } catch(ProjectError e) {
	return;
    }
}


/*
 * NcModel::readProbabilities (of one specific gc content class)
 */
void NcModel::readProbabilities( int parIndex ){
    // there are no specific noncoding parameters (all pars are reused from intergenic and coding)
}

/*
 * readAllParameters
 */
void NcModel::readAllParameters(){
}

/*
 * NcModel::initSnippetProbs
 */
void NcModel::initSnippetProbs() {
    if (initSnippetProbs5)
	delete initSnippetProbs5;
    if (snippetProbs5) 
	delete snippetProbs5;
    if (rInitSnippetProbs5)
	delete rInitSnippetProbs5;
    if (rSnippetProbs5)
	delete rSnippetProbs5;
    if (rSnippetProbs3)
	delete rSnippetProbs3;
    if (intronSnippetProbs)
	delete intronSnippetProbs;

    initSnippetProbs5 = new SnippetProbs(sequence, k);
    snippetProbs5 = new SnippetProbs(sequence, k);
    rInitSnippetProbs5 = new SnippetProbs(sequence, k, false);
    rSnippetProbs5 = new SnippetProbs(sequence, k, false);
    rSnippetProbs3 = new SnippetProbs(sequence, k, false);
    intronSnippetProbs = new SnippetProbs(sequence, IntronModel::k);
    haveSnippetProbs = true;
}

/*
 * NcModel::initAlgorithms
 *
 * makes a correction on the transition matrix "trans" and the vector of ancestors
 * this is called after initViterbiAlgorithms
 */
void NcModel::initAlgorithms( Matrix<Double>& trans, int cur){
    if (utype == utr5intron)
	pUtr5Intron = trans[cur][cur].doubleValue();
    if (utype == utr3intron)
	pUtr3Intron = trans[cur][cur].doubleValue();
    if (utype == rutr5intron)
	prUtr5Intron = trans[cur][cur].doubleValue();
    if (utype == rutr3intron)
	prUtr3Intron = trans[cur][cur].doubleValue();

    if (!initAlgorithmsCalled) {
      seqProb(-1,-1, false, -1);
      if (tssProbsPlus.size() != dnalen+1){
	tssProbsPlus.assign(dnalen+1, -1.0);
	tssProbsMinus.assign(dnalen+1, -1.0);
      } 
 
      initSnippetProbs5->setEmiProbs(&utr5init_emiprobs.probs);
      snippetProbs5->setEmiProbs(&utr5_emiprobs.probs);
      rInitSnippetProbs5->setEmiProbs(&utr5init_emiprobs.probs);
      rSnippetProbs5->setEmiProbs(&utr5_emiprobs.probs);
      rSnippetProbs3->setEmiProbs(&utr3_emiprobs.probs);
      intronSnippetProbs->setEmiProbs(&IntronModel::emiprobs.probs);
    }
    initAlgorithmsCalled = true;
    haveSnippetProbs = false;
}

/*
 * NcModel::viterbiForwardAndSampling
 */
void Nc::viterbiForwardAndSampling( ViterbiMatrixType& viterbi,
					  ViterbiMatrixType& forward,
					  int state,
					  int base,
					  AlgorithmVariant algovar,
					  OptionListItem& oli) const {
  /* 
   *             | [begin signal]|                        | [end signal] |
   *  pred.State |TSS            | markov chain           | TTS          |
   *              ASS                                       DSS
   *              rTTS                                      rTSS
   *              rDSS                                      rASS
   *  --------------------------------------------------------------------------
   *  predProb   |         notEndPartProb                 | endPartProb  |
   *              <--------------------- lenPartProb -------------------->
   */ 
    OptionsList *optionslist = NULL;
    Feature* extrinsicexons = NULL;
    int endOfPred, leftMostEndOfPred, rightMostEndOfPred, beginOfEndPart, endOfBioExon;
    Double maxPredProb, predProb, endPartProb, notEndPartProb;
    Double fwdsum, fwdsummand;
    Double emiProb, extrinsicQuot, transEmiProb;
    vector<Ancestor>::const_iterator it;
    maxPredProb = fwdsum = 0.0;
    extrinsicQuot = 1.0;
    if (algovar == doSampling)
	optionslist = new OptionsList();
    
    getEndPositions(base, beginOfEndPart, endOfBioExon);
    switch (nctype) {
	case utr5single:
	    leftMostEndOfPred = base - (max_exon_length - Constant::trans_init_window + Constant::tss_upwindow_size);
	    rightMostEndOfPred = base - Constant::tss_upwindow_size - tss_end - 1 // normal offset so that signals do not overlap
	      + Constant::trans_init_window + tss_end; // these terms to allow 5'UTR down to a min. len. of 1 (overlapping tss and tis signals), (e.g. C.elegans has very short 5'UTRs)
	    if (rightMostEndOfPred > base-1)
	      rightMostEndOfPred = base-1; // prevent infinite loop in Viterbi
	    break;
	case rutr5single:
	    leftMostEndOfPred = base - (max_exon_length - Constant::trans_init_window + Constant::tss_upwindow_size);
	    rightMostEndOfPred = base - Constant::tss_upwindow_size - 1 + Constant::trans_init_window; // see comment above
            if (rightMostEndOfPred > base-1)
              rightMostEndOfPred = base-1; 
	    break;
	case utr5init:
	    leftMostEndOfPred = base - (max_exon_length + DSS_MIDDLE + Constant::dss_end + Constant::tss_upwindow_size);
	    rightMostEndOfPred = base - Constant::tss_upwindow_size - tss_end - Constant::dss_whole_size();
	    break;
	case rutr5init:
	    leftMostEndOfPred = base - (max_exon_length + DSS_MIDDLE + Constant::dss_end + Constant::tss_upwindow_size);
	    rightMostEndOfPred = base - Constant::tss_upwindow_size - tss_end - Constant::dss_whole_size();
	    break;
	case utr5internal:
	    leftMostEndOfPred = base - (max_exon_length + DSS_MIDDLE + Constant::dss_end + Constant::ass_upwindow_size + Constant::ass_start + ASS_MIDDLE);
	    rightMostEndOfPred = base - Constant::dss_whole_size() - Constant::ass_upwindow_size - Constant::ass_whole_size();
	    break;
	case rutr5internal:
	    leftMostEndOfPred = base - (max_exon_length + DSS_MIDDLE + Constant::dss_end + Constant::ass_upwindow_size + Constant::ass_start + ASS_MIDDLE);
	    rightMostEndOfPred = base - Constant::dss_whole_size() - Constant::ass_upwindow_size - Constant::ass_whole_size();
	    break;
	case utr5term:
	    leftMostEndOfPred = base - (max_exon_length - Constant::trans_init_window + Constant::ass_upwindow_size + Constant::ass_start + ASS_MIDDLE);
	    rightMostEndOfPred = base - Constant::ass_upwindow_size - Constant::ass_whole_size();
	    if (- Constant::ass_upwindow_size - Constant::ass_whole_size() + Constant::trans_init_window + Constant::ass_end < 0) 
		// this check is neccessary, otherwise endOfPred could be larger than base
		rightMostEndOfPred = base - Constant::ass_upwindow_size - Constant::ass_whole_size()
		    + Constant::trans_init_window + Constant::ass_end; // this is so 5' terminal exons can be shorter than trans_init_window + ass_end, negative length of state!
	                                                               // in the most extreme case, the UTR exon has length 0
	    break;
	case rutr5term:
	    leftMostEndOfPred = base - (max_exon_length - Constant::trans_init_window + Constant::ass_upwindow_size + Constant::ass_start + ASS_MIDDLE);
	    rightMostEndOfPred = base - Constant::ass_upwindow_size - Constant::ass_whole_size();
	    if (- Constant::ass_upwindow_size - Constant::ass_whole_size() + Constant::trans_init_window + Constant::ass_end < 0) 
		// this check is neccessary, otherwise endOfPred could be larger than base
		rightMostEndOfPred = base - Constant::ass_upwindow_size - Constant::ass_whole_size()
		    + Constant::trans_init_window + Constant::ass_end; // this is so 5' terminal exons can be shorter than trans_init_window, negative length of state!
	    break;
	case utr5intron: case rutr5intron: case utr3intron: case rutr3intron:
	    leftMostEndOfPred = base - 1;
	    rightMostEndOfPred = base -1;
	    break;
	case utr3single:
	    leftMostEndOfPred = base - max3singlelength;
	    if (base != dnalen-1)
		rightMostEndOfPred = base - Constant::d_polyasig_cleavage - aataaa_boxlen;
	    else
		rightMostEndOfPred = base - 1; // right-truncated, can have any lenght>0
	    break;
	case rutr3single:
	    leftMostEndOfPred = base - max3singlelength;
	    rightMostEndOfPred = base - Constant::d_polyasig_cleavage - aataaa_boxlen;
	    break;
	case utr3init:
	    leftMostEndOfPred = base - (max_exon_length + DSS_MIDDLE + Constant::dss_end);
	    rightMostEndOfPred = base - Constant::dss_end - DSS_MIDDLE; // was: base - Constant::dss_whole_size();
	    // allow overlap of models, so utr3init can be shorter than dss_whole_size
	    break;
	case rutr3init:
	    leftMostEndOfPred = base - (max_exon_length + DSS_MIDDLE + Constant::dss_end);
	    //rightMostEndOfPred = base - Constant::dss_whole_size();
	    rightMostEndOfPred = base - Constant::dss_end - DSS_MIDDLE; // allow overlap of models, so rutr3init can be shorter than dss_whole_size
	    break;
	case utr3internal:
	    leftMostEndOfPred = base - (max_exon_length + DSS_MIDDLE + Constant::dss_end + ASS_MIDDLE + Constant::ass_start + Constant::ass_upwindow_size );
	    rightMostEndOfPred = base - Constant::dss_whole_size() - Constant::ass_upwindow_size - Constant::ass_whole_size();
	    break;
	case rutr3internal:
	    leftMostEndOfPred = base - (max_exon_length + DSS_MIDDLE + Constant::dss_end + Constant::ass_upwindow_size + Constant::ass_start + ASS_MIDDLE);
	    rightMostEndOfPred = base - Constant::dss_whole_size() - Constant::ass_upwindow_size - Constant::ass_whole_size();
	    break;
	case utr3term:
	    leftMostEndOfPred = base - (max3termlength + ASS_MIDDLE + Constant::ass_start + Constant::ass_upwindow_size );
	    if (base != dnalen-1)
		rightMostEndOfPred = base - Constant::d_polyasig_cleavage - aataaa_boxlen - Constant::ass_whole_size() - Constant::ass_upwindow_size;
	    else 
		rightMostEndOfPred = base - Constant::ass_whole_size() - Constant::ass_upwindow_size;
	    break;
	case rutr3term:
	    leftMostEndOfPred = base - (max3termlength + ASS_MIDDLE + Constant::ass_start + Constant::ass_upwindow_size );
	    rightMostEndOfPred = base - Constant::d_polyasig_cleavage - aataaa_boxlen - Constant::ass_whole_size() - Constant::ass_upwindow_size;
	    break;
	default:
	    leftMostEndOfPred = base - 1;
	    rightMostEndOfPred = base -1;
    }
    
    if (beginOfEndPart >= 0)
	endPartProb = endPartEmiProb(beginOfEndPart, base, endOfBioExon);
    else 
	endPartProb = 0;

    if (endPartProb == 0){
	// viterbi[base].erase(state);
	// if (needForwardTable(algovar))
	//     forward[base].erase(state);
	return;
    }
    
    if (!isLongUTRIntron(utype)){
    
	if (utype == utr5single || utype == utr5init){ // the tss window is allowed to be before the dna start
	    if (leftMostEndOfPred < - Constant::tss_upwindow_size)
		leftMostEndOfPred = - Constant::tss_upwindow_size;
	} else if (utype == rutr3single || utype == rutr3term){// the reverse polyA signal is allowed to be before the dna start
	    if (leftMostEndOfPred < - UtrModel::aataaa_boxlen - Constant::d_polyasig_cleavage)
		leftMostEndOfPred = - UtrModel::aataaa_boxlen - Constant::d_polyasig_cleavage;
	} else 
	    if (leftMostEndOfPred < 0)
		leftMostEndOfPred = 0;

	/*
	 * get the extrinsic exonpart information about parts falling in this range
	 */
	if (!(isIntron(utype)))
	    extrinsicexons = seqFeatColl->getExonListOvlpingRange(leftMostEndOfPred + 1 - Constant::trans_init_window, // to be safe for all cases
								  endOfBioExon,
								  isOnFStrand(utype)? plusstrand : minusstrand);

	
	seqProb(-1,-1,false,-1);      // initialize the static variables
	
	for (endOfPred = rightMostEndOfPred; endOfPred >= leftMostEndOfPred; endOfPred--){
	    const ViterbiColumnType& predVit = algovar == doSampling ? 
		(endOfPred > 0 ? forward[endOfPred] : forward[0]) :
		(endOfPred > 0 ? viterbi[endOfPred] : viterbi[0]);
	    // compute the maximum over the predecessor states probs times transition probability
	    /*
	     * check whether the starting position has a positive entry at all
	     */
	    for (it = ancestor.begin(); it != ancestor.end() && predVit[it->pos]==0; ++it);
	    if (it == ancestor.end()) continue;
	    notEndPartProb = notEndPartEmiProb(endOfPred+1, beginOfEndPart-1, endOfBioExon, extrinsicexons);
	    if (notEndPartProb <= 0.0) continue;
	    emiProb = notEndPartProb * endPartProb;
	    do {
		transEmiProb = it->val * emiProb;
		if ((utype == utr5intron || utype == rutr5intron || utype == utr3intron || utype == rutr3intron)
		    && (it->pos != state || endOfPred == 0)) // transitions into an intron are punished by malus
		    transEmiProb *= seqFeatColl->collection->malus(intronF);
		predProb = predVit[it->pos] * transEmiProb;
		if (needForwardTable(algovar)) { 
                    // endOfPred < 0 appears in left truncated state
		    fwdsummand = forward[endOfPred>=0? endOfPred:0].get(it->pos) * transEmiProb.heated();
		    fwdsum += fwdsummand;
		    if (algovar == doSampling && fwdsummand > 0)
			optionslist->add(it->pos, endOfPred, fwdsummand);
		} 
		if (predProb > maxPredProb) {
		    maxPredProb = predProb;
		    oli.state = it->pos;
		    oli.base = endOfPred;
		}
	    } while (++it != ancestor.end());
	}
    } else if (seqFeatColl){
	/*
	 * UTR intron state with variable length. Only used for UTR introns exactly matching an intron hint.
	 */
	int endOfBioIntron;
	if (utype == utr5intronvar || utype == utr3intronvar)
	    endOfBioIntron = base + Constant::ass_upwindow_size + Constant::ass_start + ASS_MIDDLE;
	else
	    endOfBioIntron = base + Constant::dss_end + DSS_MIDDLE;

	Feature *intronList = seqFeatColl->getFeatureListAt(intronF, endOfBioIntron, 
							    isOnFStrand(utype)? plusstrand : minusstrand);
	int oldEndOfPred = -INT_MAX;
	for (Feature *ihint = intronList; ihint != NULL; ihint = ihint->next) {
	    //cout << "UTRModel. Checking intron from hint: " << ihint->start << ".." << ihint->end << endl;
	    if (utype == utr5intronvar || utype == utr3intronvar)
		endOfPred = ihint->start - 1 + DSS_MIDDLE + Constant::dss_end;
	    else 
		endOfPred = ihint->start - 1 + Constant::ass_upwindow_size + Constant::ass_start + ASS_MIDDLE;
	    if (endOfPred < 0 || 
		ihint->end - ihint->start + 1 < Constant::ass_upwindow_size 
		+ Constant::ass_start + ASS_MIDDLE + DSS_MIDDLE + Constant::dss_end) 
		continue;
	    const ViterbiColumnType& predVit = viterbi[endOfPred];
	    emiProb = emiProbUnderModel(endOfPred+1, base);
	    //cout << "biolen= " << ihint->end - ihint->start + 1 << " emiProb = " << emiProb << endl;
	    if (endOfPred == oldEndOfPred)
		extrinsicQuot *= ihint->bonus;
	    else 
		extrinsicQuot = ihint->bonus;
	    emiProb *= extrinsicQuot; // option gets a bonus because it complies with the intron hint
	    for( it = ancestor.begin(); it != ancestor.end(); ++it ){
		transEmiProb = it->val * emiProb;
		predProb = predVit[it->pos] * transEmiProb;
		if (needForwardTable(algovar)) {
		    // ACHTUNG: this isn't correct in the model
		    // but the effect should be very small
		    fwdsummand = forward[endOfPred].get(it->pos) * transEmiProb.heated();
		    fwdsum += fwdsummand;
		    if (algovar==doSampling && fwdsummand > 0)
			optionslist->add(it->pos, endOfPred, fwdsummand);
		}
		if (predProb > maxPredProb) {
		    maxPredProb = predProb;
		    oli.state = it->pos;
		    oli.base = endOfPred;
		}
	    }
	    oldEndOfPred = endOfPred;
	}
    }

    switch (algovar) {
	case doSampling:
	    optionslist->prepareSampling();
	    try {
		oli = optionslist->sample();
	    } catch (ProjectError e) {
		cerr << "Sampling error in UTR model. state=" << state << " base=" << base << endl;
		throw e;
	    }
	    delete optionslist;
	    return;
	case doViterbiAndForward:
	    if (fwdsum > 0)
		forward[base][state] = fwdsum;
	case doViterbiOnly:
	    if (maxPredProb > 0)
		viterbi[base][state] = maxPredProb;
	    return;
	default:
	    // backtracking: do nothing here
	    return;
    }
}

/*
 * ===[ UtrModel::endPartEmiProb ]=====================================
 *
 * Compute the probability that the end part ends at "end" in the test sequence.
 * Return 0 if it is impossible.   
 */
Double UtrModel::endPartEmiProb(int begin, int end, int endOfBioExon) const {
    Double endPartProb = 1, extrinsicQuot = 1;
    switch (utype) {
	case utr5single: case utr5term:
	    if ((endOfBioExon + 3 <= dnalen - 1) && !GeneticCode::isStartcodon(sequence+endOfBioExon+1))
		endPartProb = 0.0;
	    break;
	case utr5internal: case utr5init: case utr3internal: case utr3init:
	    endPartProb = IntronModel::dSSProb(end - Constant::dss_whole_size() + 1, true);
	    break;
	case rutr5internal: case rutr5term: case rutr3internal: case rutr3term:
	    endPartProb = IntronModel::aSSProb(end - Constant::ass_upwindow_size - Constant::ass_whole_size() + 1, false);
	    break;
	case rutr5single: case rutr5init:
	    endPartProb = tssProb(begin);
	    break;
	case utr3single: case utr3term:
	    if (end == dnalen-1)
		return 1.0;
	    if (begin < 0 || begin + aataaa_boxlen -1 >= dnalen)
		return 0.0;
	    endPartProb = ttsProbPlus[begin];
	    break;
	case rutr3single: case rutr3init:
	    if ((end + 3 > dnalen - 1) || !GeneticCode::isRCStopcodon(sequence + end + 1))
		endPartProb = 0.0;
	    break;
	case utr5intronvar: case utr3intronvar: {
	    int asspos = end + Constant::ass_upwindow_size + Constant::ass_start + ASS_MIDDLE;
	    if (asspos >= dnalen || !isPossibleASS(asspos))
		endPartProb = 0.0;
	    break;
	}
	case rutr5intronvar: case rutr3intronvar:
	    if (!isPossibleRDSS(end + Constant::dss_end + DSS_MIDDLE))
		endPartProb = 0.0;
	    break;
	default:;
    }
    
    if (endPartProb > 0.0) {
	/*
	 * dss hints
	 */
	if (utype == utr5internal || utype == utr5init || utype == utr3internal || utype == utr3init){
	    Feature *feature = seqFeatColl->getFeatureListContaining(A_SET_FLAG(dssF), endOfBioExon+1, plusstrand);
	    if (feature)
		while (feature) {
		    extrinsicQuot *= feature->bonus;
		    feature = feature->next;
		}
	    else if (seqFeatColl->collection->hasHintsFile)
		extrinsicQuot *= seqFeatColl->collection->malus(dssF);
	}
	/*
	 * ass hints
	 */
	if (utype == rutr5term || utype == rutr5internal || utype == rutr3internal || utype == rutr3term){
	    Feature *feature = seqFeatColl->getFeatureListContaining(A_SET_FLAG(assF), endOfBioExon+1, minusstrand);
	    if (feature)
		while (feature) {
		    extrinsicQuot *= feature->bonus;
		    feature = feature->next;
		}
	    else if (seqFeatColl->collection->hasHintsFile)
		extrinsicQuot *= seqFeatColl->collection->malus(assF);
	}

	/*
	 * intronpart bonus for the part of the intron that is handled in the exon states	 
	 */
	if (utype != utr3single && utype != utr3term && utype != rutr5single && utype != rutr5init
	    && endOfBioExon < end && !isIntron(utype)) {
	    /*
	     * an intron gets the bonus for each position covered by an intronpart hint 
	     * (counted multiply for overlapping hints)
	     */
	    Feature *part, *intronList = seqFeatColl->getFeatureListOvlpingRange(A_SET_FLAG(intronpartF) | A_SET_FLAG(nonexonpartF), endOfBioExon+1, end,
										 isOnFStrand(utype)? plusstrand : minusstrand);
	    for (int i=endOfBioExon+1; i <= end; i++) {
		for (part = intronList; part!= NULL; part = part->next){
		    if (part->start<=i && part->end>=i){
			extrinsicQuot *= part->bonus;
		    }
		}
	    }
	}
    }
    return endPartProb * extrinsicQuot;
}

/*
 * notEndPartEmiProb(int begin, int end)
 * endOfMiddle is the position right before the downstream signal
 */
Double UtrModel::notEndPartEmiProb(int begin, int endOfMiddle, int endOfBioExon, Feature *exonparts) const {
    Double beginPartProb = 1, middlePartProb = 1, lenProb = 1;
    Double extrinsicQuot = 1;
    int beginOfMiddle, beginOfBioExon=-1;
    Seq2Int s2i_intron(IntronModel::k+1);

    switch( utype ){
	case utr5single:
	    beginOfMiddle = begin + Constant::tss_upwindow_size + tss_end;
	    if (endOfMiddle - beginOfMiddle + 1 >= 0)
	      middlePartProb = seqProb(beginOfMiddle, endOfMiddle, false, 0);
	    else {
	      middlePartProb = pow (2.0,  -(endOfMiddle - beginOfMiddle + 1)); // overlap between 5term and first coding exon
	    }
	    beginOfBioExon = begin + Constant::tss_upwindow_size;
	    lenProb = lenDist5Single[endOfBioExon - beginOfBioExon + 1];
	    if (begin >= 0) {
		beginPartProb = tssProb(begin);
	    } else {
		beginPartProb = pow (.25, beginOfMiddle-1); // part of tss model is before start of dna
		if (begin + Constant::tss_upwindow_size == 0) // tail probability
		    lenProb = tailLenDist5Single[endOfMiddle - begin + 1 + Constant::trans_init_window - Constant::tss_upwindow_size];
	    }
	    break;
	case utr5init:
	    beginOfMiddle = begin + Constant::tss_upwindow_size + tss_end;
	    middlePartProb = seqProb(beginOfMiddle, endOfMiddle, false, 0);
	    beginOfBioExon = begin + Constant::tss_upwindow_size;	    
	    lenProb = lenDist5Initial[endOfBioExon - beginOfBioExon + 1];
	    if (begin >= 0) {
		beginPartProb = tssProb(begin);  
	    } else {
		beginPartProb = pow (.25, beginOfMiddle-1); // part of tss model before start of dna
		if (begin + Constant::tss_upwindow_size == 0) // tail probability
		    lenProb = tailLenDist5Single[endOfBioExon - beginOfBioExon + 1];
	    }
	    break;
	case utr5internal:
	    beginPartProb = IntronModel::aSSProb(begin, true);
	    beginOfBioExon = begin + Constant::ass_upwindow_size + Constant::ass_start + ASS_MIDDLE;
	    if (beginPartProb>0.0) {
		beginOfMiddle = begin + Constant::ass_upwindow_size + Constant::ass_whole_size();
		middlePartProb = seqProb(beginOfMiddle, endOfMiddle, false, 1);
		lenProb = lenDist5Internal[endOfBioExon - beginOfBioExon + 1];
	    }
	    break;
	case rutr5internal:
	    beginPartProb = IntronModel::dSSProb(begin, false);
	    beginOfBioExon = begin + Constant::dss_end + DSS_MIDDLE;
	    if (beginPartProb>0.0) {
		beginOfMiddle = begin + Constant::dss_whole_size();
		middlePartProb = seqProb(beginOfMiddle, endOfMiddle, true, 1);
		lenProb = lenDist5Internal[endOfBioExon - beginOfBioExon + 1];
	    }
	    break;
	case utr5term:
	    beginOfBioExon = begin + Constant::ass_upwindow_size + Constant::ass_start + ASS_MIDDLE;	    
	    if (beginOfBioExon >= dnalen) 
		beginPartProb = 0.0;
	    else
		beginPartProb = IntronModel::aSSProb(begin, true);
	    if (beginPartProb>0.0) {
		beginOfMiddle = begin + Constant::ass_upwindow_size + Constant::ass_whole_size();
		if (endOfMiddle - beginOfMiddle + 1 >= 0)
		    middlePartProb = seqProb(beginOfMiddle, endOfMiddle, false, 1);
		else {
		    middlePartProb = pow (4.0,  -(endOfMiddle - beginOfMiddle + 1)); // overlap between 5term and first coding exon 
		}
		lenProb = lenDist5Terminal[endOfBioExon - beginOfBioExon + 1];
	    }
	    break;
	case rutr5term:
	    beginOfMiddle = begin;
	    beginOfBioExon = begin - Constant::trans_init_window;
	    if (endOfMiddle - beginOfMiddle + 1 >= 0)
		middlePartProb = seqProb(beginOfMiddle, endOfMiddle, true, 1);
	    else {
		middlePartProb = pow (4.0,  -(endOfMiddle - beginOfMiddle + 1)); // overlap between 5term and first coding exon 
	    }
	    lenProb = lenDist5Terminal[endOfBioExon - beginOfBioExon + 1];
	    break;
	case utr5intron: case rutr5intron:
	    beginPartProb = 1.0;
	    for (int pos = begin; pos <= endOfMiddle; pos++)
		if (pos-k >= 0)
		    try {
			middlePartProb *= IntronModel::emiprobs.probs[s2i_intron(sequence + pos - k)]; // strand does not matter!
		    } catch (InvalidNucleotideError e) {
			middlePartProb *= 0.25;
		    }
		else
		   middlePartProb *= 0.25; 
	    break;
	case rutr5single:
	    beginOfMiddle = begin;
	    beginOfBioExon = begin - Constant::trans_init_window;
	    if (endOfMiddle - beginOfMiddle + 1 >= 0)
	      middlePartProb = seqProb(beginOfMiddle, endOfMiddle, true, 0);
            else {
              middlePartProb = pow (2.0,  -(endOfMiddle - beginOfMiddle + 1)); // overlap between 5term and first coding exon
            }
	    lenProb = lenDist5Single[endOfBioExon - beginOfBioExon + 1];
	    break;
	case rutr5init:
	    beginPartProb = IntronModel::dSSProb(begin, false);
	    beginOfBioExon = begin + Constant::dss_end + DSS_MIDDLE;
	    if (beginPartProb>0.0){
		beginOfMiddle = begin + Constant::dss_whole_size();
		middlePartProb = seqProb(beginOfMiddle, endOfMiddle, true, 0);
		lenProb = lenDist5Initial[endOfBioExon - beginOfBioExon + 1];
	    }
	    break;
	case utr3single:
	    beginOfMiddle = begin;
	    beginOfBioExon = begin;
	    middlePartProb = seqProb(beginOfMiddle, endOfMiddle, false, 2);
	    if (endOfBioExon != dnalen-1)
		lenProb = lenDist3Single[endOfBioExon - beginOfBioExon + 1];
	    else {
		lenProb = tailLenDist3Single[endOfBioExon - beginOfBioExon + 1];
	    }
	    break;
	case rutr3single:
	    beginOfBioExon = begin;
	    beginOfMiddle = begin + aataaa_boxlen + Constant::d_polyasig_cleavage;
	    if (begin > 0) {
		beginPartProb = ttsProbMinus[begin + Constant::d_polyasig_cleavage];
		lenProb = lenDist3Single[endOfBioExon - beginOfBioExon + 1];
	    } else {
		if (beginOfMiddle > 0 )
		    beginPartProb = pow (.25, beginOfMiddle-1);// part of reverse tts model is before start of dna
		else 
		    beginPartProb = 1.0;
		lenProb = tailLenDist3Single[endOfBioExon - beginOfBioExon + 1];
	    }
	    if (beginPartProb > 0.0)
		middlePartProb = seqProb(beginOfMiddle, endOfMiddle, true, 2);
	    break;
	case utr3init:
	    beginOfMiddle = begin;
	    beginOfBioExon = begin;
	    if (endOfMiddle - beginOfMiddle + 1 >= 0)
		middlePartProb = seqProb(beginOfMiddle, endOfMiddle, false, 2);
	    else
		middlePartProb = pow (4.0,  -(endOfMiddle - beginOfMiddle + 1)); // overlap between utr3init and last coding exon
	    lenProb = lenDist3Initial[endOfBioExon - beginOfBioExon + 1];
	    break;
	case rutr3init:
	    beginPartProb = IntronModel::dSSProb(begin, false);
	    beginOfBioExon = begin + Constant::dss_end + DSS_MIDDLE;
	    if (beginPartProb>0.0) {
		beginOfMiddle = begin + Constant::dss_whole_size();
		if (endOfMiddle - beginOfMiddle + 1 >= 0) {
		    middlePartProb = seqProb(beginOfMiddle, endOfMiddle, true, 2);
		} else {
		    middlePartProb = pow (4.0,  -(endOfMiddle - beginOfMiddle + 1)); // overlap between rutr3init and last coding exon
		}
		lenProb = lenDist3Initial[endOfBioExon - beginOfBioExon + 1];
	    }
	    break;
	case utr3internal:
	    beginPartProb = IntronModel::aSSProb(begin, true);
	    beginOfBioExon = begin + Constant::ass_upwindow_size + Constant::ass_start + ASS_MIDDLE;
	    if (beginPartProb>0.0) {
		beginOfMiddle = begin + Constant::ass_upwindow_size + Constant::ass_whole_size();
		middlePartProb = seqProb(beginOfMiddle, endOfMiddle, false, 2);
		lenProb = lenDist3Internal[endOfBioExon - beginOfBioExon + 1];
	    }
	    break;
	case rutr3internal:
	    beginPartProb = IntronModel::dSSProb(begin, false);
	    beginOfBioExon = begin + Constant::dss_end + DSS_MIDDLE ;
	    if (beginPartProb>0.0) {
		beginOfMiddle = begin + Constant::dss_whole_size();
		middlePartProb = seqProb(beginOfMiddle, endOfMiddle, true, 2);
		lenProb = lenDist3Internal[endOfBioExon - beginOfBioExon + 1];
	    }
	    break;
	case utr3term:
	    beginPartProb = IntronModel::aSSProb(begin, true);
	    beginOfBioExon = begin + Constant::ass_upwindow_size + Constant::ass_start + ASS_MIDDLE;
	    if (beginPartProb>0.0) {
		beginOfMiddle = begin + Constant::ass_upwindow_size + Constant::ass_whole_size();
		middlePartProb = seqProb(beginOfMiddle, endOfMiddle, false, 2);
		if (endOfBioExon != dnalen-1)
		    lenProb = lenDist3Terminal[endOfBioExon - beginOfBioExon + 1];
		else 
		    lenProb = tailLenDist3Single[endOfBioExon - beginOfBioExon + 1];
	    }
	    break;
	case rutr3term:
	    beginOfBioExon = begin;
	    beginOfMiddle = begin + aataaa_boxlen + Constant::d_polyasig_cleavage;
	    if (begin > 0) {
		beginPartProb = ttsProbMinus[begin + Constant::d_polyasig_cleavage];
	    } else {
		beginPartProb = pow (.25, beginOfMiddle-1);// part of reverse tts model is before start of dna
	    }
	    if (beginPartProb > 0.0){
		middlePartProb = seqProb(beginOfMiddle, endOfMiddle, true, 2);
		lenProb = lenDist3Terminal[endOfBioExon - beginOfBioExon + 1];
	    }
	    break;
	case utr3intron: case rutr3intron:
	    beginPartProb = 1.0;
	    // begin == endOfMiddle
	    for (int pos = begin; pos <= endOfMiddle; pos++)
		if (pos-k >= 0)
		    try {
			middlePartProb = IntronModel::emiprobs.probs[s2i_intron(sequence + pos - k)]; // strand does not matter!
		    } catch (InvalidNucleotideError e) {
			middlePartProb = 0.25;
		    }
		else 
		    middlePartProb = 0.25;
	    break;
	case utr5intronvar: case utr3intronvar: case rutr5intronvar: case rutr3intronvar:
	    beginPartProb = longIntronProb(begin, endOfMiddle); // includes length prob
	    break;
	default:;
    }
    Double sequenceProb = beginPartProb * middlePartProb * lenProb;
    if (!(sequenceProb > 0.0))
	return 0.0;
    /*
     *                           extrinsicQuot
     */
    Strand strand = isOnFStrand(utype)? plusstrand : minusstrand;
    
    /*      EXON
     *
     * Multiply a bonus/malus to extrinsicQuot for every exonpart hint that is
     * covered by this biological exon
     */
    if (isExon(utype)) {
        int numEPendingInExon=0, numUPendingInExon=0, nep=0; // just used for malus
	bool UTRFSupported = false, exonFSupported = false;
	Double partBonus = 1.0;
	Strand strand = isOnFStrand(utype)? plusstrand : minusstrand;
	for (Feature *part = exonparts; part!= NULL; part = part->next){
	    if (part->type == exonpartF || part->type == UTRpartF){
		if (part->type == exonpartF && part->end >= beginOfBioExon && part->end <= endOfBioExon)
		    numEPendingInExon++;
		if (part->type == UTRpartF && part->end >= beginOfBioExon && part->end <= endOfBioExon)
		    numUPendingInExon++;
		if (strand == part->strand || part->strand == STRAND_UNKNOWN){
		    if (part->start >= beginOfBioExon && part->end <= endOfBioExon) {
			partBonus *= part->bonus;
			nep += 1;// TODO add multiplicity of hint here
		    } else if (part->type == exonpartF && 
			       (((utype == utr5single || utype == utr5term || utype == rutr3single || utype == rutr3init) &&
				 part->start >= beginOfBioExon && part->start <= endOfBioExon) ||                         // overlaps at left end of CDS,
				((utype == rutr5single || utype == rutr5term || utype == utr3single || utype == utr3init) &&
				 part->end <= endOfBioExon && part->end >= beginOfBioExon))){                               // overlaps at right end of CDS,
		      partBonus *= sqrt(part->bonus);  // give half the bonus factor
		      nep += 1;
		    }
		}
	    }
	    if (part->type == UTRF && part->start == beginOfBioExon && part->end == endOfBioExon && strand == part->strand){
		    extrinsicQuot *= part->bonus;
		    UTRFSupported = true;
	    }
	    if (part->type == exonF && strand == part->strand) {
		if ((utype == utr5init || utype == utr5internal || utype == utr3internal || utype == utr3term
		     || utype == rutr5init || utype == rutr5internal || utype == rutr3internal || utype == rutr3term) &&
		    part->start == beginOfBioExon && part->end == endOfBioExon) {
		    extrinsicQuot *= part->bonus;
		    exonFSupported = true;
		}
		if ((utype == utr3init || utr3single || utype == rutr5term || rutr5single) && (part->start < beginOfBioExon && part->end == endOfBioExon) ){
		    extrinsicQuot *= sqrt(part->bonus);
		    exonFSupported = true;
		}
		if ((utype == rutr3init || rutr3single || utype == utr5term || utr5single) && (part->start == beginOfBioExon && part->end > endOfBioExon) ){
		    extrinsicQuot *= sqrt(part->bonus);
		    exonFSupported = true;
		}
	    }
	}
	extrinsicQuot *= partBonus;
	/*
	 * Malus computation
	 */
	if (seqFeatColl && nep >=5) {
	  int zeroCov = seqFeatColl->numZeroCov(beginOfBioExon, endOfBioExon, UTRpartF, strand);
	  Double localPartMalus = seqFeatColl->collection->localPartMalus(UTRpartF, zeroCov, partBonus, nep);
	  if (localPartMalus < 1.0/partBonus) // at least have ab initio probabilities
	    localPartMalus = 1.0/partBonus;
	  extrinsicQuot *= localPartMalus;
	}
	if (seqFeatColl && seqFeatColl->collection->hasHintsFile) {
	    /* We have searched for extrinsic features.
	     * Then multiply the malus for each position of the exon.
	     * We should exclude those (few) positions where an exonpart ends.
	     * Exons longer than their exonpart hint have an incentive to become shorter.
	     */
	    if (endOfBioExon-beginOfBioExon + 1 - numEPendingInExon > 0)
		extrinsicQuot *= seqFeatColl->collection->partMalus(exonpartF, endOfBioExon-beginOfBioExon + 1 - numEPendingInExon);
	    if (endOfBioExon-beginOfBioExon + 1 - numUPendingInExon > 0)
		extrinsicQuot *= seqFeatColl->collection->partMalus(UTRpartF, endOfBioExon-beginOfBioExon + 1 - numUPendingInExon);
	    if (!exonFSupported)
		extrinsicQuot *= seqFeatColl->collection->malus(exonF);
	    if (!UTRFSupported)
		extrinsicQuot *= seqFeatColl->collection->malus(UTRF);
	}
	/*
	 * dss hints
	 */
	if (utype == rutr5internal || utype == rutr5init || utype == rutr3init || utype == rutr3internal){
	    Feature *feature = seqFeatColl->getFeatureListContaining(A_SET_FLAG(dssF), beginOfBioExon-1, minusstrand);
	    if (feature)
		while (feature) {
		    extrinsicQuot *= feature->bonus;
		    feature = feature->next;
		}
	    else if (seqFeatColl->collection->hasHintsFile)
		extrinsicQuot *= seqFeatColl->collection->malus(dssF);
	}
	/*
	 * ass hints
	 */
	if (utype == utr5internal || utype == utr5term || utype == utr3internal || utype == utr3term){
	    Feature *feature = seqFeatColl->getFeatureListContaining(A_SET_FLAG(assF), beginOfBioExon-1, plusstrand);
	    if (feature)
		while (feature) {
		    extrinsicQuot *= feature->bonus;
		    feature = feature->next;
		}
	    else if (seqFeatColl->collection->hasHintsFile)
		extrinsicQuot *= seqFeatColl->collection->malus(assF);
	}

    }

    /*
     *       INTRON
     */
    if (isIntron(utype)) {
	/*
	 * an intron gets the bonus for each position covered by an intronpart hint 
	 * (counted multiply for overlapping hints)
	 */
	Feature *intronList = seqFeatColl->getFeatureListOvlpingRange(A_SET_FLAG(intronpartF) | A_SET_FLAG(nonexonpartF), begin, endOfMiddle, strand);
	for (int i=begin; i <= endOfMiddle; i++)
	    for (Feature *part = intronList; part!= NULL; part = part->next){
		if (part->start <= i && part->end >= i)
		    extrinsicQuot *= part->bonus;
	}
    } else if (utype == utr5internal || utype == utr5term || utype == utr3internal || utype == utr3term ||
	       utype == rutr5internal || utype == rutr5init || utype == rutr3internal || utype == rutr3init)  {
	/*
	 * intronpart bonus for the part of the intron that is handled in the exon states
	 */
	Feature *part, *intronList = seqFeatColl->getFeatureListOvlpingRange(A_SET_FLAG(intronpartF) | A_SET_FLAG(nonexonpartF), begin, beginOfBioExon-1, strand);
	for (int i=begin; i <= beginOfBioExon-1; i++) {
	    for (part = intronList; part!= NULL; part = part->next){
		if (part->start<=i && part->end>=i){
		    extrinsicQuot *= part->bonus;
		}
	    }
	} 
    }

    return sequenceProb * extrinsicQuot;
}


Double UtrModel::emiProbUnderModel (int begin, int end) const {
    int beginOfEndPart, endOfBioExon;
    Feature *extrinsicexons = NULL;
    Double endProb, notEndProb, emiProb;

    getEndPositions(end, beginOfEndPart, endOfBioExon);

    if (!(isIntron(utype)))
	extrinsicexons = seqFeatColl->getExonListInRange(begin - Constant::trans_init_window, // to be safe for all cases
							 endOfBioExon,
							 isOnFStrand(utype)? plusstrand : minusstrand);
    if (inCRFTraining)
      seqProb(-1, -1, false, -1); // forget all saved information in static variables


    endProb = endPartEmiProb(beginOfEndPart, end, endOfBioExon);
    notEndProb = notEndPartEmiProb(begin, beginOfEndPart-1, endOfBioExon, extrinsicexons);
    emiProb = notEndProb * endProb;
    return emiProb;
}

// 'end' is the column in the DP table
void NcModel::getEndPositions (int end, int &beginOfEndPart, int &endOfBioExon) const {
    switch (nctype) {
    case ncsingle: case ncterm: // transcription end
	beginOfEndPart = end + 1; // no end signal
	endOfBioExon = end;
	break;
    case rncsingle: case rncinit: // transcription start
	beginOfEndPart = end + 1;
	endOfBioExon = end;
	break;
    case ncinternal: case ncinit: // donor splice site after exon
	beginOfEndPart = end - Constant::dss_whole_size() + 1;
	endOfBioExon = end - Constant::dss_end - DSS_MIDDLE;
	break;

    case rutr5init:
	    beginOfEndPart = end - Constant::tss_upwindow_size - tss_end + 1;
	    endOfBioExon = end - Constant::tss_upwindow_size;
	    break;
	case utr5internal:
	    beginOfEndPart = end - Constant::dss_whole_size() + 1;
	    endOfBioExon = end - Constant::dss_end - DSS_MIDDLE;
	    break;
	case rutr5internal:
	    beginOfEndPart = end - Constant::ass_whole_size() - Constant::ass_upwindow_size + 1;
	    endOfBioExon = end - Constant::ass_upwindow_size - Constant::ass_start - ASS_MIDDLE;
	    break;
	case utr5term:
	    beginOfEndPart = end + 1; // no end signal
	    endOfBioExon = end + Constant::trans_init_window;	    
	    break;
	case rutr5term:
	    beginOfEndPart = end - Constant::ass_whole_size() - Constant::ass_upwindow_size + 1;
	    endOfBioExon = end - Constant::ass_upwindow_size - Constant::ass_start - ASS_MIDDLE;
	    break;
	case utr5intron: case rutr5intron:
	    beginOfEndPart = end + 1;
	    endOfBioExon = end;
	    break;
	case rutr3single: case rutr3init:
	    beginOfEndPart = end + 1;
	    endOfBioExon = end;
	    break;
	case utr3single: case utr3term:
	    if (end != dnalen-1) {
		beginOfEndPart = end - Constant::d_polyasig_cleavage - aataaa_boxlen + 1;
		endOfBioExon = end;
	    } else {// end part outside sequence, right-truncated 3'UTR exon
		beginOfEndPart = dnalen;
		endOfBioExon = dnalen-1 ;
	    }
	    break;	
	case utr3init:
	    beginOfEndPart = end - Constant::dss_whole_size() + 1;
	    endOfBioExon = end - Constant::dss_end - DSS_MIDDLE;
	    break;
	case utr3internal:
	    beginOfEndPart = end - Constant::dss_whole_size() + 1;
	    endOfBioExon = end - Constant::dss_end - DSS_MIDDLE;	    
	    break;
	case rutr3internal:
	    beginOfEndPart = end - Constant::ass_whole_size() - Constant::ass_upwindow_size + 1;
	    endOfBioExon = end - Constant::ass_upwindow_size - Constant::ass_start - ASS_MIDDLE;
	    break;
	case rutr3term:
	    beginOfEndPart = end - Constant::ass_whole_size() - Constant::ass_upwindow_size + 1;
	    endOfBioExon = end - Constant::ass_upwindow_size - Constant::ass_start - ASS_MIDDLE;
	    break;
	default:
	    beginOfEndPart = end + 1;
	    endOfBioExon = end;
    }
}

/*
 * computes the probability of the emission of the sequence from left to right
 * left and right included
 */

Double UtrModel::seqProb(int left, int right, bool reverse, int type) const {
    static Double seqProb = 1.0; 
    static int oldleft = -1;
    static int oldright = -1;
    static int oldtype = -1; //type 0=5' initial/single, 1= 5', 2=3'
    static bool oldReverse = false;
    int curpos, pn;
    Seq2Int s2i(k+1);
    if (left == -1 && right == -1) {   // new initialization
	seqProb = 1.0;
	oldleft= -1;
	oldright= -1;
	oldtype=-1;
	return 1.0;
    }
    if (left > right)
	return 1.0;
  
    if (utype == rutr5single || utype == rutr5init) {
	seqProb = rInitSnippetProbs5->getSeqProb(right, right-left+1);
	return seqProb;
    }
    if (utype == rutr3term) {
	seqProb = rSnippetProbs3->getSeqProb(right, right-left+1);
	return seqProb;
    }

    //  if (utype == rutr5internal)
//	return rSnippetProbs->getSeqProb(right, right-left+1);
    if (right == oldright && left <= oldleft && reverse == oldReverse && type == oldtype) {
	for (curpos = oldleft-1; curpos >= left; curpos--){
	    try {
		if (curpos < 0 || (!reverse && curpos-k < 0))
		    seqProb *= .25;
		else {
		    pn = reverse? s2i.rc(sequence+curpos) : s2i(sequence+curpos-k);
		    if (type == 0){
		      seqProb *= utr5init_emiprobs.probs[pn];
		    } else if (type == 1){
			seqProb *= utr5_emiprobs.probs[pn];
		    } else if (type == 2) {
			seqProb *= utr3_emiprobs.probs[pn];
		    } else {
			seqProb *= IntronModel::emiprobs.probs[pn]; //for testing purposes
		    }
		}
	    } catch (InvalidNucleotideError e) {
		seqProb *= .25; //  0.25, 1/4
	    }
	}
	oldleft = left;
	return seqProb;
    }
    
    // compute everything new
    seqProb = 1.0;
    for (curpos = right; curpos >= left; curpos--) {
	try {
	    if (curpos < 0 || (!reverse && curpos-k < 0))
		seqProb *= 0.25;
	    else {
		pn = reverse? s2i.rc(sequence+curpos) : s2i(sequence+curpos-k);
		if (type == 0)
		    seqProb *= utr5init_emiprobs.probs[pn];
		else if (type == 1)
		    seqProb *= utr5_emiprobs.probs[pn];
		else if (type == 2)
		    seqProb *= utr3_emiprobs.probs[pn];
		else
		    seqProb *= IntronModel::emiprobs.probs[pn]; //for testing purposes
		if (inCRFTraining && (countEnd < 0 || (curpos >= countStart && curpos <= countEnd))){
		  if (type == 0)
		    GCutr5init_emiprobs[gcIdx].addCount(pn);
		  else if (type == 1)
		    GCutr5_emiprobs[gcIdx].addCount(pn);
		  else if (type == 2)
		    GCutr3_emiprobs[gcIdx].addCount(pn);
		  else
		    IntronModel::GCemiprobs[gcIdx].addCount(pn);
		}
	    }
	} catch (InvalidNucleotideError e) {
	    seqProb *= 0.25; // 0.25 1/4
	}
    }
    oldleft = left;
    oldright = right;
    oldReverse = reverse;
    oldtype = type;
    return seqProb;
}

Double UtrModel::tssupSeqProb (int left, int right, bool reverse) const {
    static Double seqProb;
    static int curpos;
    static Seq2Int s2i(tssup_k+1);
    seqProb = 1.0;

    for (curpos = right; curpos >= left; curpos--) {
	try {
	    if (!reverse && curpos-tssup_k >= 0)
		seqProb *= tssup_emiprobs[s2i(sequence+curpos-tssup_k)];
	    else if (reverse && curpos >= 0 && curpos + tssup_k < dnalen)
		seqProb *= tssup_emiprobs[s2i.rc(sequence+curpos)];
	    else 
		seqProb *= 0.25;
	} catch (InvalidNucleotideError e) {
	    seqProb *= 0.25;
	}
    }
    return seqProb;
}


/*
 * computes the probability of the transcription start site:
 * tss_upwindow, tata box, tss site motifs
 *
 * This Model is an exception the way it is programmed now.
 * It is not a probability distribution over the set of sequences, because
 * it distinguished the two cases: with and without TATA box
 */
Double UtrModel::tssProb(int left) const { // TODO: store results for later to be more efficient
    static Double tssMotifProb, tataMotifProb, tssupwinProb, prob, extrinsicProb;
    static bool hasTATA;
    static int reltatapos;
    static int tatapos;
    static int transstart;
    int right = left + Constant::tss_upwindow_size + tss_end - 1;
    transstart = isOnFStrand(utype)? right - tss_end + 1 : left + tss_end - 1;
    extrinsicProb = 1.0;

    Feature *tsshints = seqFeatColl->getFeatureListContaining(A_SET_FLAG(tssF), transstart, isOnFStrand(utype)? plusstrand : minusstrand);
    if (tsshints) {
	while (tsshints) {
	    extrinsicProb  *= tsshints->distance_faded_bonus(transstart);
	    tsshints = tsshints->next;
	}
    } else if (seqFeatColl->collection->hasHintsFile){
	extrinsicProb = seqFeatColl->collection->malus(tssF);
    } else
	extrinsicProb = 1.0;

    // TEMP, for speed: let transcription start be possible only every ttsSpacing-th base
    if (left % ttsSpacing != 0 && !(extrinsicProb > 1.0))
	return 0.0;
    
    if (isOnFStrand(utype)){
	if (tssProbsPlus[left] > - 0.5) // have stored value
	    return tssProbsPlus[left];
	reltatapos = findTATA(sequence + right - tss_end - d_tss_tata_max + 1, d_tss_tata_max - d_tss_tata_min - 1);
	hasTATA = (reltatapos >= 0);
	if (hasTATA){
	    tatapos = right - tss_end - d_tss_tata_max + 1 + reltatapos;
	    tssMotifProb = tssMotifTATA->seqProb(sequence + right - tss_end - tss_start + 1);
	    tataMotifProb = tataMotif->seqProb(sequence + tatapos - tata_start);
	    tssupwinProb = tssupSeqProb(left, tatapos - tata_start - 1, false)*
		tssupSeqProb(tatapos + tata_end, right - tss_end - tss_start, false);    
	} else {
	    tssMotifProb = tssMotif->seqProb(sequence + right - tss_end - tss_start + 1);
	    tataMotifProb = 1.0;
	    tssupwinProb = tssupSeqProb(left, right - tss_end - tss_start, false);
	}
	prob = tssMotifProb * tataMotifProb * tssupwinProb;
	tssProbsPlus[left] = prob;
    } else { // reverse strand
	if (tssProbsMinus[left] > - 0.5) // have stored value
	    return tssProbsMinus[left];
	reltatapos = findTATA(sequence + left + tss_end + d_tss_tata_max - 1, d_tss_tata_max - d_tss_tata_min - 1, true);
	hasTATA = (reltatapos <= 0);
	if (hasTATA){
	    tatapos = left + tss_end + d_tss_tata_max - 1 + reltatapos;
	    tssMotifProb = tssMotifTATA->seqProb(sequence + left, true, true);
	    tataMotifProb = tataMotif->seqProb(sequence + tatapos - tata_end + 1, true, true);
	    tssupwinProb = tssupSeqProb(left + tata_end + tata_start - 1, tatapos - tata_end, true)*
		tssupSeqProb(tatapos + tata_start + 1, right, true);
	} else {
	    tssMotifProb = tssMotif->seqProb(sequence + left, true, true);
	    tataMotifProb = 1.0;
	    tssupwinProb = tssupSeqProb(left + tss_end + tss_start, right, true);
	}
	prob = tssMotifProb * tataMotifProb * tssupwinProb;
	tssProbsMinus[left] = prob;
    }

    prob *= extrinsicProb;
    if (isOnFStrand(utype)){
	tssProbsPlus[left] = prob;
    } else {
	tssProbsMinus[left] = prob;
    }
    return prob;
}

/*
 * Precomputes the probability of the transcription termination site:
 * the aataaa (polyA signal) which is shortly upstream of the actual tts.
 * plus the d_polyasig_cleavage bases downstream of the polyA signal
 */
void UtrModel::computeTtsProbs(){
    Double extrinsicProb;
    Double prob;
    Double randProb = 1.0/POWER4TOTHE(aataaa_boxlen);
    Feature *ttshints;
    int ttspos, aataaa_box_begin;
    Seq2Int s2i_aataaa(aataaa_boxlen);

    // plus strand
    for (aataaa_box_begin = 0; aataaa_box_begin <= dnalen; aataaa_box_begin++) {
	// plus strand
	ttspos = aataaa_box_begin + aataaa_boxlen + Constant::d_polyasig_cleavage - 1;
	if (ttspos >= dnalen) {
	    ttsProbPlus[aataaa_box_begin] = 0;
	} else {
	    ttspos = aataaa_box_begin + aataaa_boxlen + Constant::d_polyasig_cleavage - 1;
	    ttshints = seqFeatColl->getFeatureListContaining(A_SET_FLAG(ttsF), ttspos, plusstrand);
	    extrinsicProb = 1.0;
	    if (ttshints) {
		while (ttshints) {
		    extrinsicProb  *= ttshints->distance_faded_bonus(ttspos);
		    ttshints = ttshints->next;
		}
	    } else if (seqFeatColl->collection->hasHintsFile)
		extrinsicProb = seqFeatColl->collection->malus(ttsF);
	    try {
		prob = aataaa_probs[s2i_aataaa(sequence + aataaa_box_begin)];
		prob *= prob_polya;
	    } catch (InvalidNucleotideError e) {
		prob = 0.0;
	    }
	    if ((extrinsicProb > 1.0 || aataaa_box_begin % ttsSpacing == 0) && prob == 0)//if no aataaa like pattern: allow tts every ttsSpacing-th base
		prob = (1.0-prob_polya) * randProb; // randprob = 1/4^6
	    if (prob > 0.0) { // compute prob of downstream window up to 'tts'
		prob *= ttsMotif->seqProb(sequence + aataaa_box_begin + aataaa_boxlen);
	    }
	    ttsProbPlus[aataaa_box_begin] = prob * extrinsicProb;
	}
	// minus strand
	ttspos = aataaa_box_begin - Constant::d_polyasig_cleavage;
	if (ttspos < 0 || aataaa_box_begin + aataaa_boxlen - 1 >= dnalen) {
	    ttsProbPlus[aataaa_box_begin] = 0;
	} else {
	    ttshints = seqFeatColl->getFeatureListContaining(A_SET_FLAG(ttsF), ttspos, minusstrand);
	    extrinsicProb = 1.0;
	    if (ttshints) {
		while (ttshints) {
		    extrinsicProb  *= ttshints->distance_faded_bonus(ttspos);
		    ttshints = ttshints->next;
		}
	} else if (seqFeatColl->collection->hasHintsFile)
	    extrinsicProb = seqFeatColl->collection->malus(ttsF);
        try {
	    prob = aataaa_probs[s2i_aataaa.rc(sequence + aataaa_box_begin)];
	    prob *= prob_polya;
	} catch (InvalidNucleotideError e) {
	    prob = 0.0;
	}
	if ((extrinsicProb > 1.0 || aataaa_box_begin % ttsSpacing == 0) && prob == 0)
	    prob = (1.0-prob_polya) * randProb;
	if (prob > 0.0) { // compute prob of downstream window up to 'tts'
	    prob *= ttsMotif->seqProb(sequence + ttspos, true, true);
	}
	ttsProbMinus[aataaa_box_begin] = prob * extrinsicProb;
	}
    }
}

/*
 * longIntronProb
 * Probability of a UTR intron (that is supported by a hint).
 */
Double UtrModel::longIntronProb(int internalBegin, int internalEnd) const {
    Double seqProb(1.0);
    Double lenProb;
    int internalIntronLen = internalEnd - internalBegin + 1;
    double p;
    switch (utype) {
	case utr5intronvar: p = pUtr5Intron; break;
	case utr3intronvar: p = pUtr3Intron; break;
	case rutr5intronvar: p = prUtr5Intron; break;
	case rutr3intronvar: p = prUtr3Intron; break;
	default: throw ProjectError("UtrModel::longIntronProb: Unknown alternative.");
    }
    lenProb = pow(p, internalIntronLen - 1) * (1.0-p);
    seqProb = intronSnippetProbs->getSeqProb(internalEnd, internalIntronLen);
/*    
    // malus for intron bases not covered by intronparts
    int coveredBegin = internalEnd+1;
    int coveredEnd = internalBegin-1;
    Feature *part, *intronList = seqFeatColl->getFeatureListOvlpingRange(A_SET_FLAG(intronpartF) | A_SET_FLAG(nonexonpartF), internalBegin, 
								      internalEnd, isOnFStrand(utype)? plusstrand : minusstrand);
    for (part = intronList; part!= NULL; part = part->next){
	if (part->start< coveredBegin)
	    coveredBegin = part->start;
	if (part->end > coveredEnd)
	    coveredEnd = part->end;
	for (int i=internalBegin; i<=internalEnd; i++) {
	    if (part->start<=i && part->end>=i){
		extrinsicQuot *= part->bonus;
	    }
	}
    }
    if (coveredEnd > internalEnd)
	coveredEnd = internalEnd;
    if (coveredBegin < internalBegin)
	coveredBegin = internalBegin;
    if (seqFeatColl->collection->hasHintsFile)
	extrinsicQuot *= pow (seqFeatColl->collection->malus(intronpartF), internalEnd-coveredEnd + coveredBegin-internalBegin);
    //cout << "iternallen " << internalIntronLen << " lenProb = " << lenProb << endl;
    */
    return seqProb * lenProb;
}

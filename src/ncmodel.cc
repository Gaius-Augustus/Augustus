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
Boolean         NcModel::haveSnippetProbs = false;
SnippetProbs*   NcModel::snippetProbs = NULL;
int             NcModel::boundSpacing = 10; // without hints 5' and 3' transcript end only every ttsSpacing bases, for speed
vector<Double>  NcModel::ttsProbPlus;
vector<Double>  NcModel::ttsProbMinus;
vector<Double>  NcModel::tssProbPlus;
vector<Double>  NcModel::tssProbMinus;

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
    if (--nccount == 0){
	if (snippetProbs) 
	    delete snippetProbs;
    }
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
    if (snippetProbs)
	delete snippetProbs;

    snippetProbs = new SnippetProbs(sequence, IntronModel::k);
    haveSnippetProbs = true;
}

/*
 * NcModel::initAlgorithms
 *
 * makes a correction on the transition matrix "trans" and the vector of ancestors
 * this is called after initViterbiAlgorithms
 */
void NcModel::initAlgorithms( Matrix<Double>& trans, int cur){
    /*
    if (utype == utr5intron)
	pUtr5Intron = trans[cur][cur].doubleValue();
    if (utype == utr3intron)
	pUtr3Intron = trans[cur][cur].doubleValue();
    if (utype == rutr5intron)
	prUtr5Intron = trans[cur][cur].doubleValue();
    if (utype == rutr3intron)
    prUtr3Intron = trans[cur][cur].doubleValue();*/

    if (!initAlgorithmsCalled) {
	precomputeTxEndProbs();
    }
    initAlgorithmsCalled = true;
    haveSnippetProbs = false;
}

/*
 * computeLengthDistributions
 * use a parametric distribution (negative binomial)
 */
void NcModel::computeLengthDistributions(){
    Double mean = 200;
    Double disp = 0.5; // > 0, the larger the dispersion, the larger the variance
    double r = 1/disp;
    Double p = mean / (mean + r);
    // number of successes until r failures occur, success prob is p
    // variance = mean + dispersion * mean^2
    // use recursive formular for NB-distribution to fill in the table
    lenDistSingle.resize(Constant::max_exon_len+1, 0.0);
    lenDistSingle[0] = pow(1-p, r);
    for (int k=1; i <= Constant::max_exon_len; k++){
	lenDistSingle[k] = lenDistSingle[k-1] * p * (k+r) / (k+1);
    }
    cout << "lenDistSingle[100]  = " << lenDistSingle[100]  << " should be 0.003660507" << endl;
    cout << "lenDistSingle[1000] = " << lenDistSingle[1000] << " should be 4.681851e-06" << endl;
    cout << "lenDistSingle[5000] = " << lenDistSingle[5000] << " should be 1.212119e-22" << endl;
    lenDistInternal = lenDistSingle;
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
    switch (nctype)
	{
	case ncsingle: case rncsingle:
	    leftMostEndOfPred = base - Constant::max_exon_len;
	    rightMostEndOfPred = base - 1;
	    break;
	case ncinit: case rncinit:
	    leftMostEndOfPred = base - (Constant::max_exon_len + DSS_MIDDLE + Constant::dss_end);
	    rightMostEndOfPred = base - Constant::dss_whole_size();
	    break;
	case ncinternal: case rncinternal:
	    leftMostEndOfPred = base - (Constant::max_exon_len + DSS_MIDDLE + Constant::dss_end + Constant::ass_upwindow_size + Constant::ass_start + ASS_MIDDLE);
	    rightMostEndOfPred = base - DSS_MIDDLE - Constant::dss_end - Constant::ass_upwindow_size - Constant::ass_start - ASS_MIDDLE - 1; // => minimum exon len = 1
	    break;
	case ncterm: case rncterm:
	    leftMostEndOfPred = base - (Constant::max_exon_len + Constant::ass_upwindow_size + Constant::ass_start + ASS_MIDDLE);
	    rightMostEndOfPred = base - Constant::ass_upwindow_size - Constant::ass_whole_size();
	    break;
	default: # 4 intron states
	    leftMostEndOfPred = base - 1;
	    rightMostEndOfPred = base - 1;
	}
    
    if (beginOfEndPart >= 0)
	endPartProb = endPartEmiProb(beginOfEndPart, base, endOfBioExon);
    else 
	endPartProb = 0;

    if (endPartProb == 0)
	return;
    
    if (nctype != ncintronvar && nctype != rncintronvar){
	if (leftMostEndOfPred < 0)
	    leftMostEndOfPred = 0;

	/*
	 * get the extrinsic exonpart information about parts falling in this range
	 */
	if (!(isNcIntron(nctype)))
	    extrinsicexons = seqFeatColl->getExonListOvlpingRange(leftMostEndOfPred,
								  endOfBioExon,
								  isOnFStrand(utype)? plusstrand : minusstrand);

	
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
		if ((nctype == ncintron || nctype == rncintron)
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
	 * intron state with variable length. Only used for introns exactly matching an intron hint.
	 */
	int endOfBioIntron;
	if (nctype == ncintronvar)
	    endOfBioIntron = base + Constant::ass_upwindow_size + Constant::ass_start + ASS_MIDDLE;
	else
	    endOfBioIntron = base + Constant::dss_end + DSS_MIDDLE;

	Feature *intronList = seqFeatColl->getFeatureListAt(intronF, endOfBioIntron, 
							    isOnFStrand(utype)? plusstrand : minusstrand);
	int oldEndOfPred = -INT_MAX;
	for (Feature *ihint = intronList; ihint != NULL; ihint = ihint->next) {
	    if (nctype == ncintronvar)
		endOfPred = ihint->start - 1 + DSS_MIDDLE + Constant::dss_end;
	    else
		endOfPred = ihint->start - 1 + Constant::ass_upwindow_size + Constant::ass_start + ASS_MIDDLE;
	    if (endOfPred < 0 ||
		ihint->end - ihint->start + 1 < Constant::ass_upwindow_size 
		+ Constant::ass_start + ASS_MIDDLE + DSS_MIDDLE + Constant::dss_end) 
		continue;
	    const ViterbiColumnType& predVit = viterbi[endOfPred];
	    emiProb = emiProbUnderModel(endOfPred+1, base);
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
		cerr << "Sampling error in noncoding model. state=" << state << " base=" << base << endl;
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
 * ===[ NcModel::endPartEmiProb ]=====================================
 *
 * Compute the probability that the end part ends at "end" in the test sequence.
 * Return 0 if it is impossible.   
 */
Double NcModel::endPartEmiProb(int begin, int end, int endOfBioExon) const {
    Double endPartProb = 1, extrinsicQuot = 1;
    switch (nctype) 
	{
	case ncsingle: case ncterm: case rncsingle: case rncinit: 
	    if (end % boundSpacing != 0) // allow transcript starts only on a grid for efficiency
		endPartProb = 0.0;
	    break;
	case ncinit: case ncinternal:
	    endPartProb = IntronModel::dSSProb(begin, true);
	    break;
	case rncterm: case ncinternal:
	    endPartProb = IntronModel::aSSProb(begin, false);
	    break;
	case ncintronvar:
	    if (!isPossibleASS(end + Constant::ass_upwindow_size + Constant::ass_start + ASS_MIDDLE))
		endPartProb = 0.0;
	    break;
	case rncintronvar:
	    if (!isPossibleRDSS(end + Constant::dss_end + DSS_MIDDLE))
		endPartProb = 0.0;
	    break;
	default:;
    }
    
    if (endPartProb > 0.0) {
	/*
	 * dss hints
	 */
	if (nctype == ncinternal || nctype == ncinit){
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
	if (nctype == rncinternal || nctype == rncterm){
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
	if ((nctype == ncinit || nctype == ncinternal || nctype == rncterm || nctype == ncinternal) // splice site follows
	    && endOfBioExon < end) {
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
Double NcModel::notEndPartEmiProb(int begin, int endOfMiddle, int endOfBioExon, Feature *exonparts) const {
    Double beginPartProb = 1, middlePartProb = 1, lenProb = 1;
    Double extrinsicQuot = 1;
    int beginOfMiddle, beginOfBioExon=-1;
    Seq2Int s2i_intron(IntronModel::k+1);

    switch( nctype ){
	case ncsingle:
	    middlePartProb = seqProb(begin, endOfMiddle);
	    beginOfBioExon = begin;
	    lenProb = lenDistSingle[endOfBioExon - beginOfBioExon + 1];
	    beginPartProb = tssProbPlus[begin];
	    break;
	case ncinit:
	    middlePartProb = seqProb(begin, endOfMiddle, false);
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

// begin = endOfPred + 1
Double NcModel::emiProbUnderModel (int begin, int end) const {
    int beginOfEndPart, endOfBioExon;
    Feature *extrinsicexons = NULL;
    Double endProb, notEndProb, emiProb;

    getEndPositions(end, beginOfEndPart, endOfBioExon);

    if (!(isNcIntron(nctype)))
	extrinsicexons = seqFeatColl->getExonListInRange(begin, // to be safe for all cases
							 endOfBioExon,
							 isOnFStrand(nctype)? plusstrand : minusstrand);
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
    case rncinternal: case rncterm: // reverse acceptor splice site after exon
	beginOfEndPart = end - Constant::ass_whole_size() - Constant::ass_upwindow_size + 1;
	endOfBioExon = end - Constant::ass_upwindow_size - Constant::ass_start - ASS_MIDDLE;
	break;
    default: // an intron state: ncintron, ncintronvar, rncintron, rncintronvar
	    beginOfEndPart = end + 1;
	    endOfBioExon = end;
}

/*
 * computes the probability of the emission of the sequence from left to right
 * left and right included
 */

Double NcModel::seqProb(int left, int right) const {
    if (left > right)
	return 1.0;
    
    return snippetProbs(right, right-left+1);
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
 * here only from hints
 */
 Double NcModel::tssProb(int tsspos) const { // TODO: store results for later to be more efficient
     Double extrinsicProb(1.0);

     Feature *tsshints = seqFeatColl->getFeatureListContaining(A_SET_FLAG(tssF), tsspos, isOnFStrand(utype)? plusstrand : minusstrand);
     if (tsshints) {
	 while (tsshints) {
	     extrinsicProb  *= tsshints->distance_faded_bonus(tsspos);
	     tsshints = tsshints->next;
	 }
     } else if (seqFeatColl->collection->hasHintsFile){
	 extrinsicProb = seqFeatColl->collection->malus(tssF);
     } else
	 extrinsicProb = 1.0;

     // for speed: let transcription start be possible only every ttsSpacing-th base
     if (tsspos % boundSpacing != 0 && !(extrinsicProb > 1.0))
	 return 0.0;
    
     return extrinsicProb;
 }

/*
 * Precomputes the probability of the transcription termination and initiation sites
 * on the whole sequence. Uses the tss and tts hints, spacing (% 10) for efficiency.
 * TODO: Disallow transcript ends that are nowhere near a hint: 
 * Exploit that nc genes unsupported by hints are not considered.
 */
void NcModel::precomputeTxEndProbs(){
    Feature *f;
    int pos;

    tssProbsPlus.assign(dnalen+1, 0.0);
    tssProbsMinus.assign(dnalen+1, 0.0);
    ttsProbsPlus.assign(dnalen+1, 0.0);
    ttsProbsMinus.assign(dnalen+1, 0.0);

    // in addition allow and incentivize transcript boundaries when hinted to
    f = seqFeatColl->getAllActiveFeatures(tssF);
    while (f) {
	for (pos = f->start; pos <= f->end; pos++){
	    if (f->strand == plusstrand || f->strand == bothstrands){
		if (tssProbPlus[pos] == 0)
		    tssProbPlus[pos] = f->distance_faded_bonus(pos);
		else 
		    tssProbPlus[pos] *= f->distance_faded_bonus(pos);
	    }
	    if (f->strand == minusstrand || f->strand == bothstrands){
		if (tssProbMinus[pos] == 0)
		    tssProbMinus[pos] = f->distance_faded_bonus(pos);
		else 
		    tssProbMinus[pos] *= f->distance_faded_bonus(pos);
	    }	    
	}
	f = f->next;
    }
    f = seqFeatColl->getAllActiveFeatures(ttsF);
    while (f) {
	for (pos = f->start; pos <= f->end; pos++){
	    if (f->strand == plusstrand || f->strand == bothstrands){
		if (ttsProbPlus[pos] == 0)
		    ttsProbPlus[pos] = f->distance_faded_bonus(pos);
		else 
		    ttsProbPlus[pos] *= f->distance_faded_bonus(pos);
	    }
	    if (f->strand == minusstrand || f->strand == bothstrands){
		if (ttsProbMinus[pos] == 0)
		    ttsProbMinus[pos] = f->distance_faded_bonus(pos);
		else 
		    ttsProbMinus[pos] *= f->distance_faded_bonus(pos);
	    }	    
	}
	f = f->next;
    }
    // allow transcript boundaries every 10th base (boundSpacing) for efficiency
    Double tssMalus = seqFeatColl->collection->malus(tssF);
    Double ttsMalus = seqFeatColl->collection->malus(ttsF);
    for (pos = 0; pos <= dnalen; pos += boundSpacing) {
	if (tssProbPlus[pos] == 0)
	    tssProbPlus[pos] = tssMalus;
	if (tssProbMinus[pos] == 0)
	    tssProbMinus[pos] = tssMalus;
	if (ttsProbPlus[pos] == 0)
	    ttsProbPlus[pos] = ttsMalus;
	if (ttsProbMinus[pos] == 0)
	    ttsProbMinus[pos] = ttsMalus;
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

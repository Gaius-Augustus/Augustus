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
#include "exonmodel.hh"   // so coding exon length distributions can be reused here
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
Boolean         NcModel::initAlgorithmsCalled = false;
Boolean         NcModel::haveSnippetProbs = false;
SegProbs*       NcModel::segProbs = NULL;
int             NcModel::boundSpacing = 10; // without hints 5' and 3' transcript end only every ttsSpacing bases, for speed
vector<Double>  NcModel::ttsProbPlus;
vector<Double>  NcModel::ttsProbMinus;
vector<Double>  NcModel::tssProbPlus;
vector<Double>  NcModel::tssProbMinus;

/*
 * NcModel constructor
 */
NcModel::NcModel() {
    nctype = toStateType( Properties::getProperty("/NcModel/type", nccount++) );
    strand = isOnFStrand(nctype)? plusstrand : minusstrand;
}

/*
 * destructor
 */
NcModel::~NcModel( ){
    if (--nccount == 0){
	if (segProbs)
	    delete segProbs;
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
 * registerPars
 */
void NcModel::registerPars (Parameters* parameters){
}

/*
 * readAllParameters
 */
void NcModel::readAllParameters(){
}

/*
 * buildModel
 */
void NcModel::buildModel( const AnnoSequence* annoseq, int parIndex){
}

/*
 * printProbabilities
 */
void NcModel::printProbabilities( int idx, BaseCount *bc, const char* suffix ){
}

/*
 * NcModel::initSnippetProbs
 */
void NcModel::initSnippetProbs() {
    if (segProbs)
	delete segProbs;
    segProbs = new SegProbs(sequence, IntronModel::k);
    haveSnippetProbs = true;
}

/*
 * NcModel::initAlgorithms
 *
 * makes a correction on the transition matrix "trans" and the vector of ancestors
 * this is called after initViterbiAlgorithms
 */
void NcModel::initAlgorithms( Matrix<Double>& trans, int cur){
    if (nctype == ncintron || nctype == rncintron)
	pIntron = trans[cur][cur].doubleValue();
    else 
	pIntron = 0.999;

    if (!initAlgorithmsCalled) {
	precomputeTxEndProbs();
	computeLengthDistributions();
    }
    eop.clear();
    initAlgorithmsCalled = true;
}

/*
 * NcModel::updateToLocalGC
 */
void NcModel::updateToLocalGC(int from, int to){
   segProbs->setEmiProbs(&IntronModel::emiprobs.probs, from, to);
}

/*
 * computeLengthDistributions
 * use a parametric distribution (negative binomial)
 */
void NcModel::computeLengthDistributions(){
    double mean = 200;
    double disp = 0.5; // > 0, the larger the dispersion, the larger the variance
    double r = 1/disp;
    double p = mean / (mean + r);
    int minSingleExonLen = 11; // a more or less random lower bound
    // number of successes until r failures occur, success prob is p
    // variance = mean + dispersion * mean^2
    // use recursive formular for NB-distribution to fill in the table
    lenDistSingle.resize(Constant::max_exon_len+1, 0);
    lenDistSingle[0] = pow(1-p, r);
    for (int k=1; k <= Constant::max_exon_len; k++)
	lenDistSingle[k] = lenDistSingle[k-1] * p * (k+r-1) / k;
    lenDistInternal = ExonModel::lenDistInternal; // was: = lenDistSingle; causes problems in distinguishing coding from nc
    for (int k=0; k < minSingleExonLen; k++)
	lenDistSingle[k] = 0;
}

/*
 * NcModel::viterbiForwardAndSampling
 */
void NcModel::viterbiForwardAndSampling( ViterbiMatrixType& viterbi,
					 ViterbiMatrixType& forward,
					 int state,
					 int base,
					 AlgorithmVariant algovar,
					 OptionListItem& oli) {
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
	default: // 4 intron states
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
	if (!(isNcIntron(nctype))){
	    extrinsicexons = seqFeatColl->getFeatureListOvlpingRange(A_SET_FLAG(exonF) | A_SET_FLAG(exonpartF),
								     leftMostEndOfPred, endOfBioExon, strand);
	    if (true) { // TODO option: allowOnlyExonHintedNCExons
	       // determine first hinted exon
	       int minE = rightMostEndOfPred + 1;
	       for (Feature *f = extrinsicexons; f != NULL; f = f->next)
		  if (f->start < minE)
		     minE = f->start;
	       if (minE > leftMostEndOfPred){
		  leftMostEndOfPred = minE;
		  if (leftMostEndOfPred > rightMostEndOfPred - 200){
		     leftMostEndOfPred = rightMostEndOfPred - 200;
		     if (leftMostEndOfPred < 0)
			leftMostEndOfPred = 0;
		  }
	       }
	    }
	}
	if (algovar == doBacktracking || algovar == doSampling)
	   eop.clear();
	eop.init(); // initialize iterator on list with possible endOfPreds
	
	for (endOfPred = rightMostEndOfPred; endOfPred >= leftMostEndOfPred; eop.decrement(endOfPred)){
	    const ViterbiColumnType& predVit = algovar == doSampling ? 
		(endOfPred > 0 ? forward[endOfPred] : forward[0]) :
		(endOfPred > 0 ? viterbi[endOfPred] : viterbi[0]);
	    // compute the maximum over the predecessor states probs times transition probability

	    // check whether the starting position has a positive entry at all
	    for (it = ancestor.begin(); it != ancestor.end() && predVit[it->pos]==0; ++it);
	    if (it == ancestor.end()) continue;

	    notEndPartProb = notEndPartEmiProb(endOfPred+1, beginOfEndPart-1, endOfBioExon, extrinsicexons);
	    if (notEndPartProb == 0)
	       continue;
	    try {
	       eop.update(endOfPred);
	    } catch (ProjectError e) {
	       cerr << "Error in EOPList::update.\tbase=" << base << "\t" << stateTypeNames[nctype]
		    << "\tnotEndPartProb=" << notEndPartProb << endl;
	       throw e;
	    }

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
							    isOnFStrand(nctype)? plusstrand : minusstrand);
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
	case ncsingle: case ncterm: 
	    endPartProb = ttsProbPlus[end]; break;
	case rncsingle: case rncinit: 
	    endPartProb = tssProbMinus[end]; break;
	case ncinit: case ncinternal:
	    endPartProb = IntronModel::dSSProb(begin, true);
	    break;
	case rncterm: case rncinternal:
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
	default:
	    ;
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
										 isOnFStrand(nctype)? plusstrand : minusstrand);
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

    switch( nctype )
	{
	case ncsingle:
	    beginOfBioExon = begin;
	    beginPartProb = tssProbPlus[begin];
	    if (beginPartProb > 0) {
	       middlePartProb = segProbs->getSeqProb(begin, endOfMiddle);
	       lenProb = lenDistSingle[endOfBioExon - beginOfBioExon + 1];
	    }
	    break;
	case ncinit:
	    beginOfBioExon = begin;
	    beginPartProb = tssProbPlus[beginOfBioExon];
	    if (beginPartProb > 0) {
	       middlePartProb = segProbs->getSeqProb(begin, endOfMiddle);
	       lenProb = lenDistInternal[endOfBioExon - beginOfBioExon + 1];
	    }
	    break;
	case ncinternal:
	    beginPartProb = IntronModel::aSSProb(begin, true);
	    beginOfBioExon = begin + Constant::ass_upwindow_size + Constant::ass_start + ASS_MIDDLE;
	    if (beginPartProb > 0) {
		beginOfMiddle = begin + Constant::ass_upwindow_size + Constant::ass_whole_size();
		middlePartProb = segProbs->getSeqProb(beginOfMiddle, endOfMiddle);
		lenProb = lenDistInternal[endOfBioExon - beginOfBioExon + 1];
	    }
	    break;
	case ncterm:
	    beginOfBioExon = begin + Constant::ass_upwindow_size + Constant::ass_start + ASS_MIDDLE;	    
	    if (beginOfBioExon >= dnalen) 
		beginPartProb = 0;
	    else
		beginPartProb = IntronModel::aSSProb(begin, true);
	    if (beginPartProb > 0) {
		beginOfMiddle = begin + Constant::ass_upwindow_size + Constant::ass_whole_size();
		if (endOfMiddle - beginOfMiddle + 1 >= 0)
		    middlePartProb = segProbs->getSeqProb(beginOfMiddle, endOfMiddle);
		else {
		    middlePartProb = pow (4.0,  -(endOfMiddle - beginOfMiddle + 1));
		}
		lenProb = lenDistInternal[endOfBioExon - beginOfBioExon + 1];
	    }
	    break;
	case rncsingle:
	    beginOfBioExon = begin;
	    beginPartProb = ttsProbMinus[beginOfBioExon];
	    if (beginPartProb > 0) {
	       middlePartProb = segProbs->getSeqProb(begin, endOfMiddle);
	       lenProb = lenDistSingle[endOfBioExon - beginOfBioExon + 1];
	    }
	    break;
	case rncinternal:
	    beginPartProb = IntronModel::dSSProb(begin, false);
	    beginOfBioExon = begin + Constant::dss_end + DSS_MIDDLE;
	    if (beginPartProb > 0) {
		beginOfMiddle = begin + Constant::dss_whole_size();
		middlePartProb = segProbs->getSeqProb(beginOfMiddle, endOfMiddle);
		lenProb = lenDistInternal[endOfBioExon - beginOfBioExon + 1];
	    }
	    break;
	case rncterm:
	    beginOfBioExon = begin;
	    beginPartProb = ttsProbMinus[beginOfBioExon];
	    if (beginPartProb > 0) {
	       if (endOfMiddle - begin + 1 >= 0)
		  middlePartProb = segProbs->getSeqProb(begin, endOfMiddle);
	       else {
		  middlePartProb = pow (4.0,  -(endOfMiddle - begin + 1));
	       }
	       lenProb = lenDistInternal[endOfBioExon - beginOfBioExon + 1];
	    }
	    break;
	case rncinit:
	    beginPartProb = IntronModel::dSSProb(begin, false);
	    beginOfBioExon = begin + Constant::dss_end + DSS_MIDDLE;
	    if (beginPartProb > 0){
		beginOfMiddle = begin + Constant::dss_whole_size();
		middlePartProb = segProbs->getSeqProb(beginOfMiddle, endOfMiddle);
		lenProb = lenDistInternal[endOfBioExon - beginOfBioExon + 1];
	    }
	    break;
	case ncintron: case rncintron: // strand does not matter
	    for (int pos = begin; pos <= endOfMiddle; pos++)
		if (pos - IntronModel::k >= 0)
		    try {
			middlePartProb *= IntronModel::emiprobs.probs[s2i_intron(sequence + pos - IntronModel::k)]; 
		    } catch (InvalidNucleotideError e) {
		       middlePartProb *= (float) 0.25;
		    }
		else
		   middlePartProb *= (float) 0.25;
	    break;
	case ncintronvar: case rncintronvar: // introns that are supported by a hint
	    {
   	        int intronLen = endOfBioExon - begin + 1;
	        if (nctype == ncintronvar)
		  intronLen += Constant::dss_end + DSS_MIDDLE;
		else 
		    intronLen += Constant::ass_upwindow_size + Constant::ass_start + ASS_MIDDLE;
		lenProb = IntronModel::lenDist[IntronModel::getD()] * pow(1.0 - 1.0/IntronModel::getMAL(), (int) (intronLen-IntronModel::getD()));
		// was before an independent distribution for nc introns:
		// int internalIntronLen = endOfMiddle - begin + 1;
		// lenProb = pow(pIntron, internalIntronLen - 1) * (1.0-pIntron);
		// it appears to be better for distinguishing coding from nc to use the same distribution
		middlePartProb = segProbs->getSeqProb(begin, endOfMiddle);
		break;
	    }
	default: ;
	}
    
    Double sequenceProb = beginPartProb * middlePartProb * lenProb;
    if (!(sequenceProb > 0))
       return sequenceProb;
    /*
     *                           extrinsicQuot
     */
    
    /*      EXON
     *
     * Multiply a bonus/malus to extrinsicQuot for every exonpart hint that is
     * covered by this biological exon
     */
    if (isExon(nctype)) {
        int numEPendingInExon=0, nep=0; // just used for malus
	bool exonFSupported = false;
	Double partBonus = 1;
	for (Feature *part = exonparts; part != NULL; part = part->next){
	    if (part->type == exonpartF){
		if (part->type == exonpartF && part->end >= beginOfBioExon && part->end <= endOfBioExon)
		    numEPendingInExon++;
		if (strand == part->strand || part->strand == STRAND_UNKNOWN){
		    if (part->start >= beginOfBioExon && part->end <= endOfBioExon) {
			partBonus *= part->bonus;
			nep += 1;// TODO add multiplicity of hint here
		    }
		}
	    }
	    if (part->type == exonF && strand == part->strand && part->start == beginOfBioExon && part->end == endOfBioExon) {
		extrinsicQuot *= part->bonus;
		exonFSupported = true;
	    }
	}
	extrinsicQuot *= partBonus;
	/*
	 * Malus computation
	 */
	if (seqFeatColl && nep >=1) { // was 5
	    int zeroCov = seqFeatColl->numZeroCov(beginOfBioExon, endOfBioExon, exonpartF, strand);
	    Double localPartMalus = seqFeatColl->collection->localPartMalus(exonpartF, zeroCov, partBonus, nep);
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
	    if (!exonFSupported)
		extrinsicQuot *= seqFeatColl->collection->malus(exonF);
	}
	/*
	 * dss hints
	 */
	if (nctype == rncinternal || nctype == rncinit){
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
	if (nctype == ncinternal || nctype == ncterm){
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
    if (isIntron(nctype)) {
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
    } else if (nctype == ncinternal || nctype == ncterm || nctype == rncinternal || nctype == rncinit) {
	/*
	 * intronpart bonus for the part of the intron that is handled in the exon states
	 */
	Feature *part, *intronList = seqFeatColl->getFeatureListOvlpingRange(A_SET_FLAG(intronpartF) | A_SET_FLAG(nonexonpartF), begin, beginOfBioExon-1, strand);
	if (intronList){
	   for (int i=begin; i <= beginOfBioExon-1; i++) {
	      for (part = intronList; part != NULL; part = part->next){
		 if (part->start <= i && part->end >= i){
		    extrinsicQuot *= part->bonus;
		 }
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
    if (inCRFTraining){
	// TODO; make sure there are no old values reused
	// seqProb(-1, -1, false, -1); // forget all saved information in static variables
    }
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
}

/*
 * computes the probability of the emission of the sequence from left to right
 * left and right included
 *

Double NcModel::seqProb(int left, int right) const {
    if (left > right)
	return 1;
    
    return segProbs->getSeqProb(left, right);
    }*/

/*
 * Precomputes the probability of the transcription termination and initiation sites
 * on the whole sequence. Uses the tss and tts hints, spacing (% 10) for efficiency.
 * TODO: Disallow transcript ends that are nowhere near a hint: 
 * Exploit that nc genes unsupported by hints are not considered.
 */
void NcModel::precomputeTxEndProbs(){
    Feature *f;
    int pos;

    tssProbPlus.assign(dnalen+1, 0);
    tssProbMinus.assign(dnalen+1, 0);
    ttsProbPlus.assign(dnalen+1, 0);
    ttsProbMinus.assign(dnalen+1, 0);

    // allow and incentivize transcript boundaries when hinted to
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

    // allow transcript boundaries every boundSpacing-th base for efficiency
    // but only near exonpart hint boundaries
    Double tssMalus = seqFeatColl->collection->malus(tssF);
    Double ttsMalus = seqFeatColl->collection->malus(ttsF);
    f = seqFeatColl->getAllActiveFeatures(exonpartF);
    vector<int> v;
    v.reserve(3);
    while (f) {
	// allow tss and tts near every exonpart boundary
	v.clear();
	v.push_back(f->start);
	v.push_back(boundSpacing * (f->start/boundSpacing));
	v.push_back(boundSpacing * (1 + f->start/boundSpacing));
	for (vector<int>::iterator p = v.begin(); p!= v.end(); ++p) {
	    if (*p >= 0 && *p <= dnalen){ // p a possible left boundary
		if (tssProbPlus[*p] == 0)
		    tssProbPlus[*p] = tssMalus;
		if (ttsProbMinus[*p] == 0)
		    ttsProbMinus[*p] = ttsMalus;
	    }
	}
	v.clear();
	v.push_back(f->end);
	v.push_back(boundSpacing * (f->end/boundSpacing));
	v.push_back(boundSpacing * (1 + f->end/boundSpacing));
	for (vector<int>::iterator p = v.begin(); p!= v.end(); ++p){
	    if (*p >= 0 && *p <= dnalen){ // p a possible right boundary
		if (tssProbMinus[*p] == 0)
		    tssProbMinus[*p] = tssMalus;
		if (ttsProbPlus[*p] == 0)
		    ttsProbPlus[*p] = ttsMalus;
	    }
	}
	f = f->next;
    }
}

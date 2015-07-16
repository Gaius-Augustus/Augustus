/**********************************************************************
 * file:    statemodel.cc
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  base class for exonmodel, intronmodel, igenicmodel and utrmodel
 * authors: Mario Stanke (mario@gobics.de), Stafilarakis, Oliver Keller
 * 
 * date    |   author      |  changes 
 * --------|---------------|------------------------------------------ 
 * 17.11.01| Stafilarakis  | creation of the file
 *         | Mario Stanke  | computeEmiFromPat, makeProbsFromCounts, determineShortPatterns
 * 14.04.03| Mario Stanke  | threw out the shared object stuff
 * 21.09.05| Mario Stanke  | untranslated region
 **********************************************************************/

#include "statemodel.hh"

// project includes
#include "exonmodel.hh"
#include "intronmodel.hh"
#include "igenicmodel.hh"
#include "utrmodel.hh"
#include "ncmodel.hh"
#include "pp_scoring.hh"
#include "extrinsicinfo.hh"

// standard C/C++ includes
#include <iostream>


// definition of class variables
vector<Boolean>*           StateModel::shortpattern = NULL;
int                        StateModel::dnalen = -1;
int                        StateModel::countStart = -1;
int                        StateModel::countEnd = -1;
const char*                StateModel::sequence = NULL;
SequenceFeatureCollection* StateModel::seqFeatColl = NULL;
PP::SubstateModel*         StateModel::profileModel = NULL;
const vector<StateType>*   StateModel::stateMap = NULL;
int                        StateModel::activeWinLen = 1;
int                        StateModel::gcIdx = 0;
ContentStairs*             StateModel::cs = NULL;

/* --- StateModel methods ------------------------------------------ */

/*
 * initPredecessors
 */
void StateModel::initPredecessors(Matrix<Double>& trans, int self) {
    ancestor.clear();
    for( int i = 0; i < trans.getColSize(); i++ )
        if( trans[i][self] != 0 ) 
	    ancestor.push_back(Ancestor(i, trans[i][self]));
}

/*
 * init
 */
void StateModel::init() {
    ExonModel::init();
    IntronModel::init();
    UtrModel::init();
    IGenicModel::init();
    NcModel::init();

    // by linking submaps we reach at most as far as the length of 
    // equalD state + MAX_LINKCOUNT times geometric states + splicesite states
    // add this to the maximum exon length
    activeWinLen = ExonModel::getMaxStateLen() + IntronModel::getD() + MAX_LINKCOUNT;
}

/**
 * newStateModelPtr: Get a pointer to a StateModel object.
 *
 * creates a new StateModel descendant; the type is determined by
 * by the string _path_, currently one of:
 *  intronmodel, exonmodel, igenicmodel, utrmodel, ncmodel
 *
 */
StateModel* StateModel::newStateModelPtr(const char* path) {
    string intronMdlStr = "intronmodel", exonMdlStr = "exonmodel", igenicMdlStr = "igenicmodel", utrMdlStr = "utrmodel", ncMdlStr = "ncmodel";
    string pathString(path);
    try {
	if (pathString == intronMdlStr)
	    return new IntronModel();
	else if (pathString == exonMdlStr)
	    return new ExonModel();
	else if (pathString == igenicMdlStr)
	    return new IGenicModel();
	else if (pathString == utrMdlStr)
	    return new UtrModel();
	else if (pathString == ncMdlStr)
	    return new NcModel();
    } catch (ProjectError e) {
	cerr << e.getMessage() << endl;
    }
    throw ProjectError("StateModelPtr: error creating model \"" + string(path) + "\".");
    return (StateModel*) 0;
}

	    
/*
 * StateModel::determineShortPatterns
 * for i=0, .. , k-1 is patcounts[i] the number of patterns with base4 representation i
 * A ** midpn ** A
 * C             C
 * G             G
 * T             T
 *
 * If for some x the number sum_y #(midpn|y) is less than minCount, then 
 * the by 1 shorter pattern is used later: shortpattern[midpn] = true;
 */
void StateModel::determineShortPatterns(const vector<Integer>& patcounts, int k, int minCount){
    if (k<2) {
	cerr << "StateModel::determineShortPatterns: "
	     << "Ignoring the possibility of shortening the patterns." << endl;
	return;
    }
    int n;
    int shortsize = POWER4TOTHE(k);
    if (!shortpattern) {
	shortpattern = new vector<Boolean>(shortsize, false);
    }
    for (int pn = 0; pn < shortsize; pn++) {
	n = patcounts[4*pn + 0] + patcounts[4*pn + 1] + patcounts[4*pn + 2] + patcounts[4*pn + 3];
	(*shortpattern)[pn] = (n < minCount);
    }
    // cout << " #shortpatterns = " << numshort << " / " << shortsize << endl;
}

void StateModel::makeProbsFromCounts(vector<Double> &patprobs , const vector<Integer> &patcounts, 
				     int k, Double pseudocount, Boolean shorten){
    int size = POWER4TOTHE(k+1);
    if (size != patcounts.size() || size != patprobs.size() )
	throw ProjectError("StateModel::makeProbsFromCounts: size doesn't match k");
    if (shorten && !shortpattern) {
	cerr << "StateModel::makeProbsFromCounts: Couldn't shorten the patterns." <<endl;
	cerr << "Don't know which to shorten yet." << endl;
	shorten = false;
    }

    Double sum, normsum = 0.0;
    if (!shorten || k<2) { // normal estimation 
	for(int pn=0; pn < size; pn+=4) {
	    sum = patcounts[pn+0]+patcounts[pn+1]+patcounts[pn+2]+patcounts[pn+3];
	    for (int i=0; i<4; i++) {
		normsum += patprobs[pn+i] = patcounts[pn+i] + pseudocount;
	    }
	}
    } else { // shorten the pattern by 1 if specified in the shortpattern vector
	int shortsize = POWER4TOTHE(k-1);
	for(int pn=0; pn < shortpattern->size(); pn++) {
	    if (!((*shortpattern)[pn])) { // do it normally here also
		sum = patcounts[4*pn+0]+patcounts[4*pn+1]+patcounts[4*pn+2]+patcounts[4*pn+3];
		for (int i=0; i<4; i++) {
		    normsum += patprobs[4*pn+i] = patcounts[4*pn+i] + pseudocount;
		}
	    } else { // pool the data for 4 different nucs at the beginning
		sum = 0.0;
		int midpn = pn % shortsize;
		for (int a=0; a<4; a++) {// build the sum = *|midpn|*
		    for(int b=0; b<4; b++) {
			sum += patcounts[4*(a*shortsize+midpn) + b];
		    }  
		}
		for(int b=0; b<4; b++) {
		    // p is the number of patterns "a|midpn|b" for any a
		    Double p = patcounts[4*midpn + b]+
			patcounts[4*(  shortsize+midpn) + b]+
			patcounts[4*(2*shortsize+midpn) + b]+
			patcounts[4*(3*shortsize+midpn) + b];
			normsum += patprobs[4*pn+b] = p/4 + pseudocount;
		    
		}
	    }
	}
    }
    // normalize
    for (int pn=0; pn < size; pn++) 
	patprobs[pn] /= normsum;
}

void StateModel::computeEmiFromPat(const vector<Double>& patprobs, vector<Double>& emiprobs, Integer k){
    /*
     * Emission Probabilities for order k emissions only
     * They are computed using the Pls also.
     */
    int size = POWER4TOTHE(k+1);
    for( int i = 0; i < size; i += 4 ){
	Double sum = patprobs[i] + patprobs[i+1] + patprobs[i+2] + patprobs[i+3]; 
	for( int nuk = 0; nuk < 4; nuk++ )
	    if (k>0) {
		emiprobs[i+nuk] = patprobs[i+nuk]/sum;
	    } else {
		emiprobs[i+nuk] = patprobs[i+nuk];
	    }
    }
}

// this is stuff that needs to be called once for each new sequence
void StateModel::prepareViterbi(const char* dna, int len, 
				const vector<StateType>& smap) {
    sequence = dna;
    dnalen = len;
    stateMap = &smap; // needed in exonmodel to determine predecessor type
    if (profileModel) {
	PP::DNA::initSeq(dna, len);
	profileModel->initScores();
	profileModel->initHitSequences();
    }
    ExonModel::setORF();
}

void StateModel::readProbabilities(int newParIndex) {
    ExonModel::readProbabilities(newParIndex);
    IntronModel::readProbabilities(newParIndex);
    IGenicModel::readProbabilities(newParIndex);
    UtrModel::readProbabilities(newParIndex);
}

void StateModel::resetPars(){
    ExonModel::resetPars();
    IntronModel::resetPars();
    IGenicModel::resetPars();
    UtrModel::resetPars();
    NcModel::resetPars();
}

void StateModel::updateToLocalGC(int idx, int from, int to){
    gcIdx = idx;
    ExonModel::updateToLocalGC(from, to);
    IntronModel::updateToLocalGC(from, to);
    IGenicModel::updateToLocalGC(from, to);
    if (Constant::utr_option_on)
	UtrModel::updateToLocalGC(from, to);
    if (Constant::nc_option_on)
	NcModel::updateToLocalGC(from, to);
}

void StateModel::readAllParameters(){
  ExonModel::readAllParameters();
  IntronModel::readAllParameters();
  IGenicModel::readAllParameters();
  UtrModel::readAllParameters();
  NcModel::readAllParameters();
}

void StateModel::storeGCPars(int idx){
  ExonModel::storeGCPars(idx);
  IntronModel::storeGCPars(idx);
  IGenicModel::storeGCPars(idx);
  UtrModel::storeGCPars(idx);
}

void StateModel::resetModelCounts() {
     // reset model count of the statemodels which have a count (all except igenic)
    ExonModel::resetModelCount();
    IntronModel::resetModelCount();
    UtrModel::resetModelCount();
    NcModel::resetModelCount();
}

int StateModel::getActiveWindowStart(int base) {
    int result = base - activeWinLen;
    if (seqFeatColl) {
	// check if we have hints that make lessD states longer
	// we might already be in the longass state following it
	int endOfBioIntron = base - Constant::ass_end;
	for (Feature* ihint = seqFeatColl->getFeatureListOvlpingRange(intronF, endOfBioIntron, endOfBioIntron, plusstrand); 
		 ihint != NULL; ihint = ihint->next) {
	    int endOfPred = ihint->start - 1 + DSS_MIDDLE + Constant::dss_end;
	    if (result > endOfPred)
		result = endOfPred;
	}
	// analogously for rlessD state, followed by rlongdss
	endOfBioIntron = base - Constant::dss_start;
	for (Feature* ihint = seqFeatColl->getFeatureListOvlpingRange(intronF, endOfBioIntron, endOfBioIntron, minusstrand); 
		 ihint != NULL; ihint = ihint->next) {
	    int endOfPred = ihint->start - 1 + Constant::ass_upwindow_size + Constant::ass_start + ASS_MIDDLE;
	    if (result > endOfPred)
		result = endOfPred;
	}
    }
    return result;
}



/* --- Snippet* methods -------------------------------------------- */

Double SnippetProbs::getElemSeqProb(int base, int len){
    Double seqProb = 1.0;
    Seq2Int s2i(k+1);
    if (forwardStrand) {
	for (int curpos = base; curpos >= base - len + 1; curpos--) {
	    try {
		if (curpos-k >= 0)
		    seqProb *= (*emiprobs)[s2i(dna+curpos-k)];
		else 
		    seqProb *= 0.25;
	    } catch (InvalidNucleotideError e) {
		seqProb *= 0.25;
	    }
	}
    } else {
	for (int curpos = base; curpos >= base - len + 1; curpos--) {
	    try {
		if (curpos >= 0 && curpos+k < n)
		    seqProb *= (*emiprobs)[s2i.rc(dna+curpos)];
		else 
		    seqProb *= 0.25;
	    } catch (InvalidNucleotideError e) {
		seqProb *= 0.25;
	    }
	}
    }
    return seqProb;
}

Double SnippetProbs::getSeqProb(int base, int len){
    Double p;
    if (len == 0)
	return 1.0;
    if (base<0 || base > n-1)
	throw ProjectError(string("SnippetProbs::getSeqProb: Index out of range: ") + itoa(base));
    SnippetList * sl = snippetlist[base];
    if (sl != NULL) {
	if (sl->last->len < len) {
	    p = getSeqProb(base - sl->last->len, len - sl->last->len) * sl->last->prob;
	    addProb(base, len, p);
	} else {
	    int partlen = len;
	    p = sl->getProb(base, partlen); // search for previously computed value
	    if (partlen != len) {                           // not exactly found,
		Double prest = 1.0;
		if (partlen == 0) {
		    prest = getElemSeqProb(base, len); // compute rest
		    addProb(base, len, prest); // and store
		} else {
		    prest = getSeqProb(base-partlen, len-partlen); // compute rest
		}
		p *= prest;
	    }
	}
    } else {
	p = getElemSeqProb(base, len);
	addProb(base, len, p);
    }
    return p;
}

void SnippetProbs::addProb(int base, int len, Double p) {
    SnippetListItem *sni = new SnippetListItem(p, len);
    if (snippetlist[base] == NULL) {
	snippetlist[base] = new SnippetList();
	snippetlist[base]->first = snippetlist[base]->last = sni;
    } else if (snippetlist[base]->last->len < len){ // insert at end
	snippetlist[base]->last->next = sni;
	snippetlist[base]->last = sni;
    } else { // sort into at the right position
	SnippetListItem *s = snippetlist[base]->first;
	if (s->len > len) { // insert at very beginning
	    sni->next = s;
	    snippetlist[base]->first = sni;
	    return;// bugfix: inserted 060606
	}
	while (s && s->next && s->next->len < len)
	    s = s->next;
	if (s->next && s->next->len == len) {
	    cerr << "SnippetProb: tried to add snippet of same length: prob " << p << " original prob " << s->next->prob << " len " << len << endl;
	    for (SnippetListItem *snit = snippetlist[base]->first; snit != NULL; snit=snit->next)
		cout << base <<  "\t" << snit->len << " " << snit->prob << endl;
	} else {
	    sni->next = s->next;
	    s->next = sni;
	}
    }
}


/*
 * getProb
 * if the prob for the len is not stored then partlen
 * is the largest length below len for wich the prob is stored
 * and prob is that probability.
 */
Double SnippetList::getProb(int base, int &partlen) {
    SnippetListItem *s = first;
    if (!s || s->len > partlen) {
	partlen = 0;
	return 1.0;
    }
    while (s->next && s->next->len < partlen)
	s = s->next;
    if (!s->next || s->next->len > partlen) {
	partlen = s->len;
	return s->prob;
    } else {
	return s->next->prob;
    }
}

/*
 * SegProbs::setEmiProbs
 */
void SegProbs::setEmiProbs(vector<Double> *emiprobs, int from, int to) {
    this->emiprobs = emiprobs;
    Double p;
    if (from < 0 || to < 0) { // if nothing else is requested, use the whole sequence
	from = 1;
	to = n;
    }
    if (from == 1)
	cumProds[0] = (float) .25;
    if (to > n)
	to = n;
    if (cumProds[from-1] == 0){
	cerr << "SegProbs::setEmiProbs: previous emiprobs not computed. This should not happen." << endl;
	cumProds[from-1] = 1;
    }
    for (int i = from; i <= to; i++){
	if (forwardStrand) {
	    try {
		if (i<k)
		    throw InvalidNucleotideError('\0');
		p = (*emiprobs)[s2i(dna + i - k)];
	    } catch (InvalidNucleotideError){
		p = (float) 0.25;
	    }
	} else {
	    try {
		if (i >= n-k)
		    throw InvalidNucleotideError('\0');
		p = (*emiprobs)[s2i.rc(dna + i)];
	    } catch (InvalidNucleotideError){
		p = (float) 0.25;
	    }
	}
	cumProds[i] = cumProds[i-1] * p;
    }
    //    cout << "setEmiProbs " << "cumProds[" << from << "]= " << cumProds[from] << "\t" << "cumProds[" << to << "]= " << cumProds[to] << endl;
}

/*
 * SegProbs::getSeqProb
 */
Double SegProbs::getSeqProb(int from, int to) {
    if (from == to){
	try {
	    if (forwardStrand) {
		if (to < k)
		    return (float) .25;
		return (*emiprobs)[s2i(dna + to - k)];
	    } else
		return (*emiprobs)[s2i.rc(dna + to)];
	} catch (InvalidNucleotideError) {
	    return (float) .25;
	}
    } else if (from > to)
       return 1;
    if (to > n)
	to = n;
    Double r = cumProds[to];
    if (from < 1)
	return r;
    else
	return r / cumProds[from - 1];
}

/*
 * make the Viterbi and Forward loop more efficient by memoizing the possible
 * endOfPred (left interval boundaries) for each state
 *    |     |        |  |   |           |
 *                   |  |   |           |    |  |  | (next time endPartProb > 0)
 * C                        A          D1      D2  B
 *  possibleEndOfPred is a decreasingly sorted list of endOfPred positions, for 
 *  which notEndPartProb > 0
 *  eopit points to the current list element and iterates from right to left
 */

void EOPList::decrement(int &endOfPred){
   if (inCache && eopit != possibleEndOfPreds.end() && *eopit == endOfPred
       && ++eopit != possibleEndOfPreds.end())
      endOfPred = *eopit;
   else
      endOfPred--;
}


void EOPList::update(int endOfPred){
    if (possibleEndOfPreds.empty()) {
	possibleEndOfPreds.push_front(endOfPred);
	return;
    }
    if (endOfPred == *eopit) // case A
	return;
    if (endOfPred >= possibleEndOfPreds.front()){ // case B
	if (endOfPred > possibleEndOfPreds.front())
	    possibleEndOfPreds.push_front(endOfPred);
	eopit = possibleEndOfPreds.begin();
	return;
    }
    if (endOfPred < possibleEndOfPreds.back()){ // case C
	eopit = possibleEndOfPreds.insert(possibleEndOfPreds.end(), endOfPred);
	return;
    }
    if (endOfPred < *eopit){ // case D
	list<int>::iterator nxt = eopit;
	++nxt;
	if (nxt != possibleEndOfPreds.end() && endOfPred == *nxt){ // D1
	    eopit = nxt;
	    inCache = true;
	    return;
	}
	if (nxt == possibleEndOfPreds.end() || endOfPred > *nxt){ // D2
	    eopit = possibleEndOfPreds.insert(nxt, endOfPred);
	    return;
	}
	// this line should very rarely be reached (e.g. when the support of the length distribution is not an interval)
	//cerr << "EOPList::update. Error: endOfPred = " << endOfPred << "\teopit=" << *eopit << "\tinCache=" << inCache << endl;
	while (eopit != possibleEndOfPreds.end() && *eopit > endOfPred)
	    ++eopit;
	update(endOfPred);
	return;
    }
    // this should not happen either
    throw ProjectError("Error 2 in UtrModel::updatePossibleEOPs");
}

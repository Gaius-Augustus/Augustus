/**********************************************************************
 * file:    gene.cc
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  
 * authors: Mario Stanke, mario@gobics.de
 *
 * date    |   author      |  changes 
 * --------|---------------|------------------------------------------
 * 17.06.02| Mario Stanke  | creation of the file
 * 29.11.02| Mario Stanke  | making partial genes possible
 * 21.08.05| Mario Stanke  | added Gene::printCodingSeq
 * 23.10.06| Mario Stanke  | Evidence, SrcEvidence to attach hint summary to genes
 * 29.05.07| Mario Stanke  | introduce frame modifier, so truncated exons get the right frame
 * 25.07.08| Mario Stanke  | updated getInducedStatePath to deal with UTRs
 * 11.03.12| Mario Stanke  | frame_compatible for overlapping coding regions (bacteria)
 * 25.03.12| Mario Stanke  | frame_compatible for exon and hint (corrected compatibility for peptide hints)
 * 08.09.12| Stefanie König| reverse type of reverse strand exons added in reverseGeneSequence
 **********************************************************************/

#include "gene.hh"

// project includes
#include "intronmodel.hh"  // for getD
#include "hints.hh"

// standard C/C++ includes
#include <iostream>
#include <iomanip>  // for setw, setprecision
#include <climits>

// output options
bool Gene::print_start=false;
bool Gene::print_stop=false;
bool Gene::print_introns=false;
bool Gene::print_cds=false;
bool Gene::print_exonnames=false;
bool Gene::stopCodonExcludedFromCDS=false;
bool Gene::print_utr=false;
bool Gene::print_tss=false;
bool Gene::print_tts=false;
bool Gene::gff3 = false;
bool Gene::print_blocks = false;


/* --- Evidence methods -------------------------------------------- */

bool operator<(const SrcEvidence& e1, const SrcEvidence& e2){
    return (e1.srcname < e2.srcname);
}

void Evidence::add(string source, string name){
    bool alreadyHaveSrc = false;
    list<SrcEvidence>::iterator seit = sourceEvidence.begin();
    while(!alreadyHaveSrc && seit != sourceEvidence.end()) {
	if (seit->srcname == source) {
	    alreadyHaveSrc = true;
	    seit->groupnames.push_back(name);
	    seit->freq++;
	}
	seit++;
    }
    if (!alreadyHaveSrc) {
        sourceEvidence.push_back(SrcEvidence(source));
	sourceEvidence.back().groupnames.push_back(name);
	sourceEvidence.back().freq++;
    }
    numEvidence++;
}

void Evidence::add(string source){
    bool alreadyHaveSrc = false;
    list<SrcEvidence>::iterator seit = sourceEvidence.begin();
    while(!alreadyHaveSrc && seit != sourceEvidence.end()) {
	if (seit->srcname == source) {
	    alreadyHaveSrc = true;
	    seit->freq++;
	}
	seit++;
    }
    if (!alreadyHaveSrc) {
	sourceEvidence.push_back(SrcEvidence(source));
	sourceEvidence.back().freq++;
    }
    numEvidence++;
}

void Evidence::print(){
    sourceEvidence.sort();
    for (list<SrcEvidence>::iterator seit = sourceEvidence.begin(); seit != sourceEvidence.end(); seit++) {
	cout << "# " << setw(6) << seit->srcname << ":" << setw(4) << seit->freq << " ";
	if (withNames) {
	    int grouplistlen = 0;
	    list<string>::iterator gnit = seit->groupnames.begin();
	    while( grouplistlen < 80 && gnit != seit->groupnames.end()){
		if (gnit->length()>0) {
		    if (grouplistlen == 0) {
			cout << "(";
		    } else 
			cout << ",";
		    cout << *gnit;
		    grouplistlen += gnit->length() + 1;
		}
		gnit++;
	    }
	    if (gnit != seit->groupnames.end())
		cout << ",...";
	    if (grouplistlen>0) 
		cout << ")";
	}
	cout << endl;
    }
}

/*
 * stateEvidenceSummary
 * Only count each source once for each state.
 */
Evidence *stateEvidenceSummary(State *state1, State *state2){
    Evidence *evi = new Evidence(false);
    Evidence *se;
    int numStatesWithEvidence = 0;
    while (state1) {
	se = state1->evidence;
	if (se) {
	    if (!se->sourceEvidence.empty())
		numStatesWithEvidence++;
	    for (list<SrcEvidence>::iterator srcevit = se->sourceEvidence.begin(); srcevit != se->sourceEvidence.end(); srcevit++) {
		evi->add(srcevit->srcname);
	    }
	}
	state1 = state1->next;
    }
    while (state2) {
	se = state2->evidence;
	if (se) {
	    if (!se->sourceEvidence.empty())
		numStatesWithEvidence++;
	    for (list<SrcEvidence>::iterator srcevit = se->sourceEvidence.begin(); srcevit != se->sourceEvidence.end(); srcevit++) {
		evi->add(srcevit->srcname);
	    }
	}
	state2 = state2->next;
    }

    evi->numEvidence = numStatesWithEvidence;
    return evi;
}

int lenStateList(State *head) {
    int len=0;
    while(head) {
	len++;
	head = head->next;
    }
    return len;
}


/* --- State methods ----------------------------------------------- */

// reading frame, in case it is a (coding) exon, usually depends only on the type
int State::frame(){
    return mod3(stateReadingFrames[type] + framemod);
}

bool State::frame_compatible(const Feature *hint){
  return ((hint->strand == plusstrand && strand() == plusstrand && mod3(end - hint->start + 1 - hint->frame - frame())==0)
	  || (hint->strand == minusstrand && strand() == minusstrand && mod3(end - hint->end + hint->frame + frame() + 1)==0));
}

bool frame_compatible(State *ex1, State* ex2){
  return ((isOnFStrand(ex1->type) && isOnFStrand(ex2->type) && mod3(ex2->end - ex1->end - ex2->frame() + ex1->frame()) == 0)
	  || (!isOnFStrand(ex1->type) && !isOnFStrand(ex2->type) && mod3(ex2->end - ex1->end + ex2->frame() - ex1->frame()) == 0));
}

bool State::operator< (const State &other) const {
  return (begin<other.begin || (begin==other.begin && end < other.end));
}

bool State::operator== (const State &other) const {
  return (begin==other.begin && end==other.end);
}

State *State::getBiologicalState() {
    State *bioState = new State();
    int beginShift=0, endShift=0;
    
    switch (type) {
	case igenic :
	    beginShift = 0;
	    endShift   = Constant::trans_init_window;
	    break;
	case singleG:
	    beginShift = Constant::trans_init_window;
	    endShift   = 0;
	    break;
	case rsingleG:
	    beginShift = 0;
	    endShift   = - Constant::trans_init_window;
	    break;
	case initial0: case initial1: case initial2:
	    beginShift = Constant::trans_init_window;
	    if (!(truncated & TRUNC_RIGHT))
	      endShift   = Constant::dss_start;
	    else
		bioState->framemod = mod3(-Constant::dss_start);
	    break;
	case rterminal0: case rterminal1: case rterminal2:
	    beginShift = 0;
	    if (!(truncated & TRUNC_RIGHT))
	      endShift   = Constant::ass_end;
	    else {
		bioState->framemod = mod3(Constant::ass_end);
	    }
	    break;
	case internal0: case internal1: case internal2:
	    if (!(truncated & TRUNC_LEFT))
		beginShift = - Constant::ass_end;
	    if (!(truncated & TRUNC_RIGHT))
		endShift   = Constant::dss_start;
	    else 
		bioState->framemod = mod3(-Constant::dss_start);
	    break;
	case rinternal0: case rinternal1: case rinternal2:
	    if (!(truncated & TRUNC_LEFT))
		beginShift = - Constant::dss_start;
	    if (!(truncated & TRUNC_RIGHT))
	      endShift   = Constant::ass_end;
	    else
		bioState->framemod = mod3(Constant::ass_end);
	    break;
	case terminal:
	    if (!(truncated & TRUNC_LEFT))
		beginShift = - Constant::ass_end;
	    endShift   = 0;
	    break;
	case rinitial:
	    if (!(truncated & TRUNC_LEFT))
		beginShift = - Constant::dss_start;
	    endShift   = - Constant::trans_init_window;
	    break;
	case intron_type:
	    if (!(truncated & TRUNC_LEFT))
		beginShift = Constant::dss_start;
	    else
		beginShift = - 1;
	    if (!(truncated & TRUNC_RIGHT))
		endShift   = - Constant::ass_end;
	    break;
	case rintron_type:
	    if (!(truncated & TRUNC_LEFT))
		beginShift = Constant::ass_end;
	    else
		beginShift = - 1;
	    if (!(truncated & TRUNC_RIGHT))
		endShift   = - Constant::dss_start;
	    break;
	case utr5single:
	    beginShift = Constant::tss_upwindow_size;
	    endShift   = Constant::trans_init_window;
	    break;
	case rutr5single:
	    if (!(truncated & TRUNC_LEFT))
		beginShift = - Constant::trans_init_window;
	    else
		beginShift = -begin;	
	    endShift   = - Constant::tss_upwindow_size;
	    break;
	case utr5init:
	    beginShift = Constant::tss_upwindow_size;
	    endShift   = -Constant::dss_end-DSS_MIDDLE;
	    break;
	case rutr5init:
	    beginShift = Constant::dss_end+DSS_MIDDLE;
	    endShift   = -Constant::tss_upwindow_size;
	    break;    
	case utr5internal:
	    beginShift = Constant::ass_upwindow_size + Constant::ass_start + ASS_MIDDLE;
	    endShift   = -Constant::dss_end-DSS_MIDDLE;
	    break;
	case rutr5internal:
	    beginShift = Constant::dss_end+DSS_MIDDLE;
	    endShift   = -Constant::ass_upwindow_size - Constant::ass_start - ASS_MIDDLE;
	    break;
	case utr5term:
	    beginShift = Constant::ass_upwindow_size + Constant::ass_start + ASS_MIDDLE;
	    endShift   = Constant::trans_init_window;
	    break;
	case rutr5term:
	    beginShift = -Constant::trans_init_window;
	    endShift   = -Constant::ass_upwindow_size - Constant::ass_start - ASS_MIDDLE;
	    break;
	case utr3single:
	    if ((truncated & TRUNC_LEFT) && begin == 1) // left truncated
		beginShift = -1;
	    endShift   = 0;
	    break;
	case rutr3single:
	    if (begin < 0) // left truncated
		beginShift = -begin;
	    endShift   = 0;
	    break;
	case utr3init:
	    beginShift = 0;
	    endShift   = -Constant::dss_end - DSS_MIDDLE;
	    break;
	case rutr3init:
	    beginShift = Constant::dss_end + DSS_MIDDLE;
	    endShift   = 0;
	    break;
	case utr3internal:
	    beginShift = Constant::ass_upwindow_size + Constant::ass_start + ASS_MIDDLE;
	    endShift   = - Constant::dss_end - DSS_MIDDLE;
	    break;
	case rutr3internal:
	    beginShift = Constant::dss_end + DSS_MIDDLE;
	    endShift   = - Constant::ass_upwindow_size - Constant::ass_start - ASS_MIDDLE;
	    break;
	case utr3term:
	    beginShift = Constant::ass_upwindow_size + Constant::ass_start + ASS_MIDDLE;
	    endShift   = 0;
	    break;
	case rutr3term:
	    if (begin < 0) // left truncated
		beginShift = -begin;
	    endShift   = - Constant::ass_upwindow_size - Constant::ass_start - ASS_MIDDLE;
	    break;
	default:
	    beginShift = endShift = 0;
    }
    
    bioState->begin = begin + beginShift;
    bioState->end = end + endShift;
    bioState->type = type;
    bioState->next = NULL;
    bioState->hasScore = hasScore;
    bioState->apostprob = apostprob;
    bioState->truncated = truncated;
    return bioState;
}


/*
 * setTruncFlag
 * set the left and right truncated flag if appropriate
 * This is for truncated 'interval' states.
 */
void State::setTruncFlag(int end, int predEnd, int dnalen){
    if (end == dnalen-1 && (isInitialExon(type) || isInternalExon(type) ||
			    isRTerminalExon(type) || isRInternalExon(type) ||
			    isIntron(type) || type == utr3single || type == utr3term))
	truncated |= TRUNC_RIGHT;
    if ((predEnd == -1 || predEnd == 0 ) && (isInternalExon(type) || type == terminal ||
					       isRInternalExon(type) || type == rinitial || 
					       isCodingIntron(type) || 
					     (isExon(type) && (is3UTR(type) || is5UTR(type)))))
	truncated |= TRUNC_LEFT;
}


/* --- StatePath methods ------------------------------------------- */

/**
 * @memo A path through the Hidden-Markov-Model states
 * @author Mario Stanke
 */
void StatePath::print() {
    State* cur=first;
    if (seqname != "")
	cout << seqname << endl;
    cout << setw(25)<< "type" << "\t" << "state id" << "\t" << setw(8) << "begin" << "\t" << setw(8) << "end" << "\t" << setw(8) << "meanprob" << "\t"
	 << setw(16) << "prob" << endl;
    while(cur) {
	cout << setw(25) << stateTypeNames[cur->type] << "\t" <<  setw(8) << cur->type << "\t" << setw(8)
	     << cur->begin << "\t" << setw(8) << cur->end << "\t" << setw(8);
	if (cur->end >= cur->begin)
	    cout << (cur->prob).getRoot(cur->end-cur->begin+1);
	else
	    cout << "?";
	cout << "\t" << setw(9) << cur->prob << endl;
	cur = cur->next;
    }
}

/*
 * equality operator for StatePath's
 * Two paths are considered equal if, all states are identical.
 */

bool StatePath::operator== (const StatePath &other) const {
    bool unequal = false;
    State *s1 = first, *s2 = other.first;
    
    while (!unequal && s1 && s2) {
	unequal = !((s1->begin == s2->begin) && (s1->end == s2->end) && (s1->type == s2->type));
	s1 = s1->next;
	s2 = s2->next;
    }
    if (!unequal && !s1 && !s2) 
	return true;
    else 
	return false;
}
/*
 * comparison operator for StatePath's
 * p1 < p2 iff the a posteriori probability of p1 is greater than that of p2
 */
bool StatePath::operator< (const StatePath &other) const{
    return (pathemiProb < other.pathemiProb);
}


/*
 * StatePath::projectOntoGeneSequence
 * makes a sequence of (possibly partial) genes out of this state path
 *
 */

Gene* StatePath::projectOntoGeneSequence(const char *genenames){
    State* curState = first, *anIntron, *lastIntron = NULL;
    int genenumber = 1;
    Gene *aGene = NULL, *firstGene = NULL, *lastGene=NULL;
    State *last5utrexon, *last3utrexon, *utrexon;
    list<PP::Match>::iterator protIt = proteinMatches.begin();
    // check whether we have an incomplete gene beginning with an intron in the CDS
    if (curState && isCodingIntron(curState->type)){
	anIntron = new State();
	anIntron->begin = curState->begin;
	anIntron->type = isOnFStrand(curState->type)? intron_type : rintron_type;
	anIntron->truncated |= curState->truncated;
	while(curState->next && isCodingIntron(curState->next->type))
	    curState = curState->next;
	anIntron->end = curState->end;
	anIntron->truncated |= curState->truncated; // could theoretically be both left and right truncated
	aGene = new Gene();
	aGene->introns = lastIntron = anIntron->getBiologicalState();
	aGene->transstart = lastIntron->begin; // if gene starts with intron then let this be the transcript start
    }
    while (curState != NULL) {
	while(curState && !isExon(curState->type)) {
	    curState = curState->next;
	}
	if (curState != NULL) {
	    if (!aGene)
		aGene = new Gene();
	    aGene->source = "AUGUSTUS";

	    aGene->id = genenames + itoa(genenumber);

	    aGene->seqname = seqname;
	    aGene->strand = isOnFStrand(curState->type)? plusstrand : minusstrand;
	    if (aGene->strand == minusstrand)
		aGene->frame = 2; // this is an intermediate value for the leftmost base, only default
	    genenumber++;
	    last5utrexon = last3utrexon = NULL;
	    if (is5UTR(curState->type)){
		while (curState && is5UTR(curState->type)){ // read in the 5' UTR
		    if (!last5utrexon)
			aGene->complete5utr = (curState->type == utr5single || curState->type == utr5init);
		    if (isExon(curState->type)){
			utrexon = curState->getBiologicalState();
			if (!last5utrexon)
			    aGene->utr5exons = utrexon;
			else
			    last5utrexon->next = utrexon;
			last5utrexon = utrexon;
		    }
		    curState = curState->next;
		}
	    } else if (is3UTR(curState->type)){
		while (curState && is3UTR(curState->type)){ // read in the reverse 3' UTR
		    if (!last3utrexon)
			aGene->complete3utr = (curState->type == rutr3single || curState->type == rutr3term);
		    if (isExon(curState->type)){
			utrexon = curState->getBiologicalState();
			if (!last3utrexon)
			    aGene->utr3exons = utrexon;
			else
			    last3utrexon->next = utrexon;
			last3utrexon = utrexon;
		    }
		    curState = curState->next;
		}
	    }
	    if (curState && isExon(curState->type)) { // otherwise the (incomplete) gene consists just of UTR
		if (curState->type == singleG || curState->type == rsingleG) {
		    aGene->exons = curState->getBiologicalState();
		} else { // multiple exon gene or incomplete
		    State *lastExon=0;		
		    if (!isInitialExon(curState->type) && !isRTerminalExon(curState->type)){
			aGene->complete = false; // doesn't start with initial or rterminal: incomplete
		    }
		    if (curState->type == terminal || curState->type == rinitial){
			aGene->exons = curState->getBiologicalState();
			// determine reading frame of gene
			if (aGene->strand == plusstrand){
			    aGene->frame = mod3(aGene->exons->frame() - aGene->exons->length());
			} else
			    // this is an intermediate value, the correct reading frame at the right end
			    // is determined below
			    aGene->frame = mod3(aGene->exons->frame() + aGene->exons->length());
		    } else {
			aGene->exons = lastExon = curState->getBiologicalState();
			// determine reading frame of gene
			if (aGene->strand == plusstrand){
			    aGene->frame = mod3(aGene->exons->frame() - aGene->exons->length());
			} else
			    // this is an intermediate value, the correct reading frame at the right end
			    // is determined below
			    aGene->frame = mod3(aGene->exons->frame() + aGene->exons->length());
			curState = curState->next;
			while(curState && curState->type != terminal && curState->type != rinitial) {
			    if (isIntron(curState->type)) {
				anIntron = new State();
				anIntron->begin = curState->begin;
				anIntron->type = isOnFStrand(curState->type)? intron_type : rintron_type;
				while(curState->next && isIntron(curState->next->type)){
				    curState = curState->next;
				}
				anIntron->end = curState->end;
				anIntron->truncated = curState->truncated;
				if (aGene->introns == NULL) {
				    aGene->introns = lastIntron = anIntron->getBiologicalState();
				} else {
				    lastIntron = lastIntron->next = anIntron->getBiologicalState();
				}
				if (lastIntron->end > aGene->transstart)
				  aGene->transend = lastIntron->end; // if gene ends with intron then let this be the transcript end
			    } else if (isInternalExon(curState->type) || isRInternalExon(curState->type)) {
				lastExon = lastExon->next = curState->getBiologicalState();
			    } else {
				throw ProjectError("state path doesn't constitute a valid gene");
			    }
			    curState = curState->next;
			}
			if (curState == NULL) {
			    // terminal exon or rinitial exon missing, or only terminal exon or only riniital exon - incomplete Gene
			    aGene->complete = false;
			} else {
			    // append terminal or riniital exon to list of exons
			    lastExon = lastExon->next = curState->getBiologicalState();
			} 
		    }  // type != terminal && type != rinitial
		} // multiple exon gene or incomplete
		if (curState)
		    curState = curState->next;
		last5utrexon = last3utrexon = NULL;
		if (curState) {
		    if (is5UTR(curState->type)){
			while (curState && is5UTR(curState->type)){ // read in the reverse 5' UTR
			    if (!(curState->next && is5UTR(curState->next->type)))
				aGene->complete5utr = (curState->type==rutr5single || curState->type==rutr5init);
			    if (isExon(curState->type)){
				utrexon = curState->getBiologicalState();
				if (!last5utrexon)
				    aGene->utr5exons = utrexon;
				else
				    last5utrexon->next = utrexon;
				last5utrexon = utrexon;
			    }
			    curState = curState->next;
			}
		    } else if (is3UTR(curState->type)){
			while (curState && is3UTR(curState->type)){ // read in the 3' UTR
			    if (!(curState->next && is3UTR(curState->next->type)))
				aGene->complete3utr = (curState->type == utr3single || curState->type == utr3term);
			    if (isExon(curState->type)){
				utrexon = curState->getBiologicalState();
				if (!last3utrexon)
				    aGene->utr3exons = utrexon;
				else
				    last3utrexon->next = utrexon;
				last3utrexon = utrexon;
			    }
			    curState = curState->next;
			}
		    }
		}
	    } else { //the gene consists just of UTR
		if (!Constant::reportUtrOnlyGenes) {// default, do not consider this gene
		    if (aGene)
			delete aGene;
		    aGene = NULL;
		}
	    }
	    //
	    // finish construction of gene
	    //
	    if (aGene){
		// add UTR introns as gaps between UTR exons
		State *intron;
		lastIntron = NULL;
		for (utrexon = aGene->utr5exons; utrexon != NULL && utrexon->next != NULL; utrexon = utrexon->next) {
		    intron = new State();
		    intron->begin = utrexon->end + 1;
		    intron->end = utrexon->next->begin - 1;
		    intron->type = intron_type;
		    intron->next = NULL;
		    if (lastIntron)
			lastIntron->next = intron;
		    else
			aGene->utr5introns = intron;
		    lastIntron = intron;
		}
		lastIntron = NULL;
		for (utrexon = aGene->utr3exons; utrexon != NULL && utrexon->next != NULL; utrexon = utrexon->next) {
		    intron = new State();
		    intron->begin = utrexon->end + 1;
		    intron->end = utrexon->next->begin - 1;
		    intron->type = intron_type;
		    intron->next = NULL;
		    if (lastIntron)
			lastIntron->next = intron;
		    else
			aGene->utr3introns = intron;
		    lastIntron = intron;
		}
		
		// determine coding length of aGene
		aGene->clength = 0;
		for (State* ee = aGene->exons; ee != NULL; ee = ee->next)
		    aGene->clength += ee->length();
	    
		if (aGene->strand == minusstrand)
		    aGene->frame = mod3(aGene->frame - aGene->clength + 1); // only != 0 for incomplete genes
	    
		// transcription boundaries
		if (aGene->utr5exons && (aGene->transstart < 0 || aGene->transstart > aGene->utr5exons->begin))
		    aGene->transstart = aGene->utr5exons->begin;
		if (aGene->utr3exons && (aGene->transstart < 0 || aGene->transstart > aGene->utr3exons->begin))
		    aGene->transstart = aGene->utr3exons->begin;
		if (last5utrexon && (aGene->transend < 0 || aGene->transend < last5utrexon->end))
		    aGene->transend = last5utrexon->end;
		if (last3utrexon && (aGene->transend < 0 || aGene->transend < last3utrexon->end))
		    aGene->transend = last3utrexon->end;

		// determine the length of the range of the coding region, initialize codingstart and -end
		if (aGene->exons){
		    aGene->codingstart = aGene->exons->begin;
		    aGene->codingend = aGene->lastExon()->end;
		    aGene->length = aGene->codingend - aGene->codingstart + 1;
		    while (protIt != proteinMatches.end() && protIt->firstBase <= aGene->codingend) {
			aGene->proteinMatches.push_back(*protIt);
			++protIt;
		    }
		} else{
		    aGene->length = 0;
		}

		if (aGene->codingend > aGene->transend)
		    aGene->transend = -1;
		if (aGene->codingstart >= 0 && aGene->codingstart < aGene->transstart)
		    aGene->transstart = -1;
	    
		// append this gene to list of genes
		if (firstGene == NULL) {
		    firstGene = lastGene = aGene;
		} else {
		    lastGene = lastGene->next = aGene;
		}
		aGene = NULL;
	    }
	} 
    }
    return firstGene;
}

/*
 * StatePath::getInducedStatePath
 *
 * Takes a list of genes and creates the path through the HMM states induced by these genes, i.e.
 * a path that backtranslates through projectOntoGeneSequence to the original list of genes.
 * This fails for incomplete genes.
 */

StatePath* StatePath::getInducedStatePath(Gene *genelist, int dnalen, bool printErrors){
    StatePath *path = new StatePath();
    path->intron_d = IntronModel::getD();
    int endOfPred = 0;
    Gene *curgene;
    State *exon;
    int end, frame=0;
    bool onFStrand;
    
    if (genelist) 
	path->seqname = genelist->seqname;
    else
	path->seqname = ""; 

    for (curgene = genelist; curgene != NULL; curgene = curgene->next){
	onFStrand = (curgene->strand == plusstrand);

	/*
	 * intergenic region from end of last gene until beginning of gene
	 */
	if (onFStrand) {
	    if(curgene->utr5exons) {
		end = curgene->utr5exons->begin - 1 - Constant::tss_upwindow_size;
	    } else {
		end = curgene->exons->begin - 1 - Constant::trans_init_window;
	    }
	} else {
	    if(curgene->utr3exons) {
		end = curgene->utr3exons->begin - 1; // gene starts right with tts
	    } else {
		// assuming: curgene->exons->type == rsingleG || isRTerminalExon(curgene->exons->type)
		end = curgene->exons->begin - 1; // gene starts right with stop codon
	    }
	}
	
	if (endOfPred < end || endOfPred == 0) {// have intergenic gap or first gene
	    for (int pos = endOfPred+1; pos <= end; pos++)
		path->push(new State(pos, pos, igenic));
	    endOfPred = end;

	    exon = onFStrand? curgene->utr5exons : curgene->utr3exons;
	    /*
	     * left UTR (single exon)
	     */
	    if (exon && !exon->next){
		end = exon->end - (onFStrand? Constant::trans_init_window : 0);
		path->push(new State(endOfPred+1, end, onFStrand? utr5single: rutr3single));
		endOfPred = end;
	    }
	    /*
	     * left UTR (multiple exon)
	     */
	    if (exon && exon->next){
		// leftmost UTR exon
		end = exon->end + (onFStrand? Constant::dss_end + DSS_MIDDLE : Constant::ass_upwindow_size + Constant::ass_start + ASS_MIDDLE);
		path->push(new State(endOfPred+1, end, onFStrand? utr5init: rutr3term));
		endOfPred = end;
		while (exon->next) {
		    // UTR intron
		    end = exon->next->begin - 1 - (onFStrand? Constant::ass_upwindow_size + Constant::ass_start + ASS_MIDDLE : Constant::dss_end + DSS_MIDDLE);
		    for (int pos = endOfPred+1; pos <= end; pos++)
			path->push(new State(pos, pos, onFStrand? utr5intron: rutr3intron));
		    endOfPred = end;
		    exon = exon->next;
		    if (exon->next) {
			// internal UTR exon
			end = exon->end + (onFStrand? Constant::dss_end + DSS_MIDDLE : Constant::ass_upwindow_size + Constant::ass_start + ASS_MIDDLE);
			path->push(new State(endOfPred+1, end, onFStrand? utr5internal: rutr3internal));
			endOfPred = end;
		    } else {
			// rightmost UTR exon
			end = exon->end - (onFStrand? Constant::trans_init_window : 0);
			path->push(new State(endOfPred+1, end, onFStrand? utr5term: rutr3init));
			endOfPred = end;
		    }
		}
	    }
	    /*
	     * single CDS gene ...
	     */
	    if (curgene->exons->next == NULL){
		end = curgene->exons->end + (onFStrand? 0 : Constant::trans_init_window);
		path->push(new State(endOfPred+1, end, onFStrand? singleG: rsingleG));
		endOfPred = end;
	    }
	    /*
	     * ... or multi CDS gene?
	     */
	    else {
		// initial exon or rterminal exon
		end = curgene->exons->end - (onFStrand? Constant::dss_start : Constant::ass_end);
		frame = mod3(onFStrand? curgene->exons->length() : 2 - curgene->exons->length() );
		path->push(new State(endOfPred+1, end, onFStrand? initialExon(frame) : rterminalExon(frame)));
		endOfPred = end;
		// loop over the internal exons
		for (exon = curgene->exons->next; exon->next != NULL; exon = exon->next){
		    // intron preceeding the internal exon
		    end = exon->begin + (onFStrand? Constant::ass_end : Constant::dss_start) - 1;
		    path->pushIntron(endOfPred+1, end, frame, onFStrand);
		    endOfPred = end;
		    
		    // internal exon
		    end = exon->end - (onFStrand? Constant::dss_start : Constant::ass_end);
		    frame = mod3(frame + (onFStrand? exon->length() : -exon->length()));
		    path->push(new State(endOfPred+1, end, onFStrand? internalExon(frame) : rinternalExon(frame)));
		    endOfPred = end;
		}
		// intron preceding the terminal exon
		end = exon->begin + (onFStrand? Constant::ass_end : Constant::dss_start) - 1;
		path->pushIntron(endOfPred+1, end, frame, onFStrand);
		endOfPred = end;
		// terminal or rinitial exon
		end = exon->end + (onFStrand? 0 : Constant::trans_init_window);
		frame = mod3(frame + (onFStrand? exon->length() : -exon->length()));
		if (((onFStrand && frame != 0) || (!onFStrand && frame != 2) ) && printErrors)
		       throw ProjectError(string("StatePath::getInducedStatePath reading frame error in sequence ") + path->seqname 
					  + ":" + itoa(endOfPred+1) + ".." + itoa(end));
		path->push(new State(endOfPred+1, end, onFStrand? terminal : rinitial));
		endOfPred = end;
	    }
	    exon = onFStrand? curgene->utr3exons : curgene->utr5exons;
	    /*
	     * right UTR (single exon)
	     */
	    if (exon && !exon->next){
		end = exon->end + (onFStrand? 0 : Constant::tss_upwindow_size);
		path->push(new State(endOfPred+1, end, onFStrand? utr3single: rutr5single));
		endOfPred = end;
	    }
	    /*
	     * right UTR (multiple exon)
	     */
	    if (exon && exon->next){
		// leftmost UTR exon
		end = exon->end + (onFStrand? Constant::dss_end + DSS_MIDDLE : Constant::ass_upwindow_size + Constant::ass_start + ASS_MIDDLE);
		path->push(new State(endOfPred+1, end, onFStrand? utr3init: rutr5term));
		endOfPred = end;
		while (exon->next) {
		    // UTR intron
		    end = exon->next->begin - 1 - (onFStrand? Constant::ass_upwindow_size + Constant::ass_start + ASS_MIDDLE : Constant::dss_end + DSS_MIDDLE);
		    for (int pos = endOfPred+1; pos <= end; pos++)
			path->push(new State(pos, pos, onFStrand? utr3intron: rutr5intron));
		    endOfPred = end;
		    exon = exon->next;
		    if (exon->next) {
			// internal UTR exon
			end = exon->end + (onFStrand? Constant::dss_end + DSS_MIDDLE : Constant::ass_upwindow_size + Constant::ass_start + ASS_MIDDLE);
			path->push(new State(endOfPred+1, end, onFStrand? utr3internal: rutr5internal));
			endOfPred = end;
		    } else {
			// rightmost UTR exon
			end = exon->end + (onFStrand? 0 : Constant::tss_upwindow_size);
			path->push(new State(endOfPred+1, end, onFStrand? utr3term: rutr5init));
			endOfPred = end;
		    }
		}
	    }
	} else { // overlapping genes (or immediately adjacent with 0 bases intergenic region)
	    // skip gene
	    cerr << "Gene " << curgene->geneid << " on sequence " << path->seqname << " overlaps previous gene. I will skip it." << endl;
	}
    } // next gene

    /*
     * intergenic region from endOfPred+1 until end of sequence
     */
    end = dnalen-1;
    for (int pos = endOfPred+1; pos<=end; pos++)
	path->push(new State(pos, pos, igenic));

    // reverse the order of the states
    path->reverse();
    return path;
}

/*
 *  StatePath::pushIntron
 *
 *  push all the states belonging to a biological intron onto the state path
 *  begin and end refer to the internal state boundaries
 *  states depend on the length of the intron
 */

void StatePath::pushIntron(int begin, int end, int frame, bool onFStrand){
    State st;
    int biolen = end - begin + 1 - Constant::dss_start - Constant::ass_end;
    int dStateLen = intron_d - DSS_MIDDLE - Constant::dss_end - Constant::ass_start
	- ASS_MIDDLE - Constant::ass_upwindow_size;
    int leftSignalendpos = begin + (onFStrand? Constant::dss_whole_size() : Constant::ass_whole_size() + Constant::ass_upwindow_size) - 1;
    int rightSignalbeginpos = end - (onFStrand? Constant::ass_whole_size() + Constant::ass_upwindow_size : Constant::dss_whole_size()) + 1;
    push(new State(begin, leftSignalendpos, onFStrand? longdssIntron(frame) : rlongassIntron(frame)));
    if (biolen <= intron_d) {
	// lessD intron state
	push(new State(leftSignalendpos + 1, rightSignalbeginpos - 1, onFStrand? lessDIntron(frame) : rlessDIntron(frame)));
    } else {
	push(new State(leftSignalendpos + 1, leftSignalendpos + dStateLen, onFStrand? equalDIntron(frame) : requalDIntron(frame)));
	// loop for the geometric intron states
	for (int pos = leftSignalendpos + dStateLen + 1; pos < rightSignalbeginpos; pos++) {
	    push(new State(pos, pos, onFStrand? geometricIntron(frame) : rgeometricIntron(frame)));
	}
    }
    push(new State(rightSignalbeginpos , end, onFStrand? longassIntron(frame) : rlongdssIntron(frame)));
}

void StatePath::reverse(){
    State *st, *last, *nextSt;
    StatePath second;
    if (first) {
	st = last = first;
	while(st) {
	    nextSt = st->next;
	    second.push(st);
	    st = nextSt;
	}
	last->next = NULL;
	first = second.first;
	second.first = NULL;
    }
}
 
/*
 * StatePath::condenseStatePath
 *
 * Whenever the same state occurs more than once in a row make one long state out of it.
 */

StatePath *StatePath::condenseStatePath(StatePath *oldpath){
    StatePath *path = new StatePath();
    State *oldstate, *temp;
    
    for (oldstate = oldpath->first; oldstate != NULL; oldstate = oldstate->next) {
        if (path->first != NULL && path->first->type == oldstate->type 
	    && !isCodingExon(oldstate->type)) { // no condensing when coding exons overlap
	    path->first->end = oldstate->end;
	    path->first->truncated |= oldstate->truncated;
	    path->first->prob *= oldstate->prob;
	} else {
	    temp = new State(*oldstate);
	    temp->next = NULL;
	    path->push(temp);
	}
    }
    path->reverse();
    path->pathemiProb = oldpath->pathemiProb;
    path->proteinMatches = oldpath->proteinMatches;
    if (oldpath->seqname != ""){    
	path->seqname = oldpath->seqname;
    }
    return path;
}


/* --- Gene methods ------------------------------------------------ */

/*
 * constructor
 */
Gene::Gene(const Gene& other){
    // copy all simple data members
    length = other.length;
    clength = other.clength;
    bc = other.bc;
    weight = other.weight;
    strand = other.strand;
    complete = other.complete;
    frame    = other.frame;
    transstart = other.transstart;
    transend = other.transend;
    complete5utr = other.complete5utr;
    complete3utr = other.complete3utr;
    apostprob = other.apostprob;
    hasProbs = other.hasProbs;
    throwaway = other.throwaway;
    viterbi = other.viterbi;
    codingstart = other.codingstart;
    codingend = other.codingend;
    geneid = other.geneid;
    // simply copy the pointer to the next gene
    next = other.next;
    // copy the short strings
    if (other.id != "")
	id = other.id;
    seqname = other.seqname;
    if (other.source != "") 
	source = other.source;

    // copy the exons and introns
    exons = introns = utr5exons = utr3exons = utr5introns = utr3introns = NULL;
    if (other.exons)
      exons = other.exons->cloneStateSequence();
    if (other.introns)
      introns = other.introns->cloneStateSequence();
    if (other.utr5exons)
      utr5exons = other.utr5exons->cloneStateSequence();
    if (other.utr3exons)
      utr3exons = other.utr3exons->cloneStateSequence();
    if (other.utr5introns)
	utr5introns = other.utr5introns->cloneStateSequence();
    if (other.utr3introns)
	utr3introns = other.utr3introns->cloneStateSequence();    
    // copy the hint lists
    supportingEvidence = incompatibleEvidence = CDSexonEvidence = CDSintronEvidence = UTR5stateEvidence = UTR3stateEvidence = NULL;
    if (other.supportingEvidence)
	supportingEvidence = new Evidence(*other.supportingEvidence);
    if (other.incompatibleEvidence)
	incompatibleEvidence = new Evidence(*other.incompatibleEvidence);
    if (other.CDSexonEvidence)
	CDSexonEvidence = new Evidence(*other.CDSexonEvidence);
    if (other.CDSintronEvidence)
	CDSintronEvidence = new Evidence(*other.CDSintronEvidence);
    if (other.UTR5stateEvidence)
	UTR5stateEvidence = new Evidence(*other.UTR5stateEvidence);
    if (other.UTR3stateEvidence)
	UTR3stateEvidence = new Evidence(*other.UTR3stateEvidence);
    
    proteinMatches = other.proteinMatches;
}

/* 
 * getCodingSequence
 *
 * creates a new string consisting of the concatenation of all exons
 * includes the stop codon, ends with a terminating character
 */ 
char* Gene::getCodingSequence(AnnoSequence *annoseq) const{
    char *codingseq = new char[clength+1];
    int curpos = 0;
    State *exon = exons;
    while (exon) {
	int cplength = exon->length();
	if (curpos + cplength <= clength) {
	    strncpy(codingseq + curpos, annoseq->sequence - annoseq->offset + exon->begin, cplength);
	    curpos += cplength;
	} else {
	    // shouldn't happen
	    throw ProjectError("Gene::getCodingSequence(): wrong coding length");
	}    
	exon = exon->next;
    }
    codingseq[clength]='\0';
    if (strand == plusstrand){
	return codingseq;
    }
    else {
	char* res = reverseComplement(codingseq);
	delete [] codingseq;
	return res;
    }
}

bool Gene::hasInFrameStop(AnnoSequence *annoseq) const{
    Seq2Int s2i(3);
    const char * codingSeq = getCodingSequence(annoseq);
    codingSeq += mod3(-frame); // if gene is incomplete, frame is the position of the first base
    while (strlen(codingSeq)>3){
	try {
	    if (GeneticCode::map[s2i(codingSeq)] < 0){
		//in-frame stop codon. Can happen when stop codon is parted by splice site, hopefully rare.
		return true;
	    }
	} catch (...) {}// because of masking
	codingSeq += 3;
    }
    return false;
}

// void Gene::computeBC(char *seq){
//     /*
//      * update the base count of exons and introns
//      */ 
//     State *stv[4], *st;
//     stv[0] = exons;
//     stv[1] = introns;
//     stv[2] = utr5exons;
//     stv[3] = utr3exons;
//     for (int i=0; i<4; i++){
// 	st = stv[i];
// 	while( st ){
// 	    for( int i = st->begin; i <= st->end; i++ )
// 		switch( seq[i] ){
// 		case 'a': case 'A':
// 		    st->bc.a++; break;
// 		case 'c': case 'C':
// 		    st->bc.c++; break;
// 		case 'g': case 'G':
// 		    st->bc.g++; break;
// 		case 't': case 'T':
// 		    st->bc.t++; break;
// 		}
// 	    st = st->next;
// 	}
//     }
// }

int Gene::numExons() const {
    int numExons = 0;
    for (State* ex = exons; ex; ex = ex->next)
	numExons++;
    return numExons;
}

State* Gene::lastExon() const{
    State *lastExon = exons;
    while(lastExon && lastExon->next)	    
	lastExon = lastExon->next;
    return lastExon;
}

bool Gene::identicalCDS(Gene *other){
    State *state, *ostate;
    state = exons;
    ostate = other->exons;
    while (state && ostate) {
	if (state->begin != ostate->begin  || state->end != ostate->end)
	    return false;
	state = state->next;
	ostate = ostate->next;
    }
    if (state || ostate)
	return false;
    return true;
}


/*
 * almostIdenticalTo
 * Transcript a is almost identical to transcript b iff
 * they are identical up to a small shift in the tss or tts.
 */
bool Gene::almostIdenticalTo(Gene *other){
    State *state, *ostate;
    if (strand != other->strand)
	return false;
    // check whether the CDS is identical
    if (!identicalCDS(other))
	return false;
    // check whether the 5'UTR is almost identical
    state = utr5exons;
    ostate = other->utr5exons;
    while (state && ostate) {
	if (state->begin != ostate->begin)
	    if (!(strand == plusstrand && state == utr5exons && state->begin >= ostate->begin - Constant::almost_identical_maxdiff && state->begin <= ostate->begin + Constant::almost_identical_maxdiff))
		return false;
	if (state->end != ostate->end)
	    if (!(strand == minusstrand && state->next == NULL && state->end >= ostate->end - Constant::almost_identical_maxdiff && state->end <= ostate->end + Constant::almost_identical_maxdiff))
		return false;
	state = state->next;
	ostate = ostate->next;
    }
    if (state || ostate)
	return false;   
    // check whether the 3'UTR is almost identical
    state = utr3exons;
    ostate = other->utr3exons;
    while (state && ostate) {
	if (state->begin != ostate->begin)
	    if (!(strand == minusstrand && state == utr3exons && state->begin >= ostate->begin - Constant::almost_identical_maxdiff && state->begin <= ostate->begin + Constant::almost_identical_maxdiff))
		return false;
	if (state->end != ostate->end)
	    if (!(strand == plusstrand && state->next == NULL && state->end >= ostate->end - Constant::almost_identical_maxdiff && state->end <= ostate->end + Constant::almost_identical_maxdiff))
		return false;
	state = state->next;
	ostate = ostate->next;
    }
    if (state || ostate)
	return false;
    return true;
}

void Gene::shiftCoordinates(int d){
    State *s, *sv[6];
    sv[0]=exons;
    sv[1]=introns;
    sv[2]=utr5exons;
    sv[3]=utr3exons;
    sv[4]=utr5introns;
    sv[5]=utr3introns;
  
    for (int i=0; i<6; i++){
	s = sv[i];
	while(s) {
	    s->begin += d;
	    s->end += d;
	    s = s->next;
	}
    }
    codingstart += d;
    codingend += d;
    if (transstart >=0)
	transstart +=d;
    if (transend >=0)
	transend +=d;

    for_each (proteinMatches.begin(), proteinMatches.end(), PP::Match::shift(d));
}

void Gene::addStatePostProbs(float p){
  State *st, *sv[4];
  sv[0]=exons;  sv[1]=introns;  sv[2]=utr5exons;  sv[3]=utr3exons;
  for (int i=0; i<4; i++){
    st = sv[i];
    while (st) {
      st->apostprob += p;
      st->hasScore = true;
      st = st->next;
    }
  }
}

void Gene::setStatePostProbs(float p){
  State *st, *sv[4];
  sv[0]=exons;  sv[1]=introns;  sv[2]=utr5exons;  sv[3]=utr3exons;
  for (int i=0; i<4; i++){
    st = sv[i];
    while (st) {
      st->apostprob = p;
      st->hasScore = true;
      st = st->next;
    }
  }
}
void Gene::setStateHasScore(bool has){
  State *st, *sv[4];
  sv[0]=exons;  sv[1]=introns;  sv[2]=utr5exons;  sv[3]=utr3exons;
  for (int i=0; i<4; i++){
    st = sv[i];
    while (st) {
      st->hasScore = has;
      st = st->next;
    }
  }
}

/*
 * Gene::addSampleCount
 * add p to the apostprob of the transcript as well to the apostprob of all states
 */
void Gene::addSampleCount(int k){
  State *st, *sv[4];
  sv[0]=exons;  sv[1]=introns;  sv[2]=utr5exons;  sv[3]=utr3exons;
  for (int i=0; i<4; i++){
    st = sv[i];
    while (st) {
      st->sampleCount += k;
      st = st->next;
    }
  }
}
/*
 * Gene::addSampleCount
 * add p to the apostprob of the transcript as well to the apostprob of all states
 */
void Gene::setSampleCount(int k){
  State *st, *sv[4];
  sv[0]=exons;  sv[1]=introns;  sv[2]=utr5exons;  sv[3]=utr3exons;
  for (int i=0; i<4; i++){
    st = sv[i];
    while (st) {
      st->sampleCount = k;
      st = st->next;
    }
  }
}

/*
 * Gene::normPostProb
 * divide the apostprob of the transcript and the apostprob of all states by n
 */
void Gene::normPostProb(float n){
  if (n==0.0)
    throw ProjectError("Gene::addPostProb:Tried to normalize without any obervations.");
  apostprob /= n;
  State *st, *sv[4];
  sv[0]=exons;  sv[1]=introns;  sv[2]=utr5exons;  sv[3]=utr3exons;

  for (int i=0; i<4; i++){
    st = sv[i];
    while (st) {
      st->apostprob /= n;
      st = st->next;
    }
  }
}


/*
 * updatePostProb
 *
 * Find exons and introns common to both genes and 
 * add the apostprob of other to the apostprob of this
 *
 */

void Gene::updatePostProb(Gene* other){
  
  // first check if the genes are totally disjoint
  int max, min;
  max = codingend;
  if (transend > max)
    max = transend;
  min = codingstart;
  if (transstart < min && transstart >=0)
    min = transstart;
  if (((other->codingstart > max) && (other->transstart<0 || other->transstart > max)) ||
      ((other->codingend < min) && (other->transend < min)))
    return;// transcripts do not overlap at all
  
  State *st, *otherst, *sv[4],*osv[4];
  sv[0]=exons;
  sv[1]=introns;
  sv[2]=utr5exons;
  sv[3]=utr3exons;
  osv[0]=other->exons;
  osv[1]=other->introns;
  osv[2]=other->utr5exons;
  osv[3]=other->utr3exons;
  for (int i=0; i<4; i++){
    st = sv[i];
    otherst = osv[i];
    // now mergecompare the two lists, assuming that they are sorted
    while (st && otherst){
      if (st->begin == otherst->begin && st->end==otherst->end){
	// add to the apostprob of each state the sampleCount of the other
	st->apostprob += otherst->sampleCount;
	otherst->apostprob += st->sampleCount;
	st = st->next;
	otherst = otherst->next;
      } else {
	if (st->begin < otherst->begin)
	  st = st->next;
	else
	  otherst = otherst->next;
      }
    }
  }
}

/*
 * meanStateProb
 * 
 * returns the geometric mean of the exon state posterior probabilities
 */
double Gene::meanStateProb(){
    if (!hasProbs)
	return 0.0;
    State *s, *sv[4];  
    sv[0]=exons;
    sv[1]=introns;
    sv[2]=utr5exons;
    sv[3]=utr3exons;
    double meanstateprob = 1.0;
    int numstates = 0;
    for (int i=0; i<4; i++){
	for (s = sv[i]; s!=NULL; s=s->next){
	    meanstateprob *= s->apostprob;
	    numstates++;
	}
    }
    meanstateprob = pow(meanstateprob, 1.0/numstates);
    return meanstateprob;
}

/*
 * Gene::addUTR
 * takes the mRNA from the first parameter and constructs the utrexons from it.
 * If gene has no coding exons and utr is incomplete at one end, then the incomplete end
 * is the one adjacent to the coding region, e.g.
 * this->strand == plusstrand 
 *     |-----           -----          ------|
 *  complete                             incomplete
 * is a 5'UTR.
 */

void Gene::addUTR(State *mrnaRanges, bool complete_l, bool complete_r){
  if (!complete_l && !complete_r && !exons)
    return; // no CDS and mRNA incomplete at both ends, can't use this data
  if (complete_l && complete_r && !exons)
    return; // no CDS and mRNA complete at both ends, can't use this data  
  State *last5utrexon=NULL, *last3utrexon=NULL, *st;
  transstart = INT_MAX, transend = -1;
  while (mrnaRanges){
    if (mrnaRanges->begin < transstart)
      transstart = mrnaRanges->begin;
    if (mrnaRanges->end > transend)
      transend = mrnaRanges->end;
    if (exons){ // coding region exists
      if (mrnaRanges->begin < codingstart){
	st = new State();
	st->begin = mrnaRanges->begin;
	st->end = (mrnaRanges->end < codingstart)? mrnaRanges->end : codingstart-1;
	if (strand == plusstrand){
	  if (!last5utrexon){
	    utr5exons = st;
	  } else {
	    last5utrexon->next = st;
	  }
	  last5utrexon = st;
	} else {
	  if (!last3utrexon){
	    utr3exons = st;
	  } else {
	    last3utrexon->next = st;
	  }
	  last3utrexon = st;
	}
      }
      if (mrnaRanges->end > codingend){
	st = new State();
	st->begin = (mrnaRanges->begin > codingend)? mrnaRanges->begin : codingend+1;
	st->end = mrnaRanges->end;
	if (strand == plusstrand){
	  if (!last3utrexon){
	    utr3exons = st;
	  } else {
	    last3utrexon->next = st;
	  }
	  last3utrexon = st;
	} else {
	  if (!last5utrexon){
	    utr5exons = st;
	  } else {
	    last5utrexon->next = st;
	  }
	  last5utrexon = st;
	}
      }
    } else { // no coding region, just add UTR to an 'empty' gene
      st = new State();
      st->begin = mrnaRanges->begin;
      st->end = mrnaRanges->end;
      if ((strand == plusstrand && !complete_r) || (strand == minusstrand && !complete_l)){
	if (!last5utrexon){
	  utr5exons = st;
	} else {
	  last5utrexon->next = st;
	}
	last5utrexon = st;
      } else {
	if (!last3utrexon){
	  utr3exons = st;
	} else {
	  last3utrexon->next = st;
	}
	last3utrexon = st;
      }
    }
    mrnaRanges = mrnaRanges->next;
  }
  if (exons) {
    complete5utr = (bool)((strand == plusstrand)? complete_l: complete_r);
    complete3utr = (bool)((strand == plusstrand)? complete_r: complete_l);
  } else { // assume, that it is complete
    if ((strand == plusstrand && !complete_r) || (strand == minusstrand && !complete_l)){
      complete5utr = true;
      complete3utr = false;
    } else {
      complete3utr = true;
      complete5utr = false;
    }
  }
  if (transstart == INT_MAX)
    transstart = -1;
  if (transstart > codingstart && codingstart>0)
    transstart = codingstart;
  if (transend > 0 && transend < codingend)
    transend = codingend;
}

/*
 * Gene::initializeGroupLists
 */
void Gene::compileExtrinsicEvidence(list<HintGroup>  *groupList){
    if (!groupList)
	return;    
    if (supportingEvidence)
	delete supportingEvidence;
    if (incompatibleEvidence)
	delete incompatibleEvidence;
    supportingEvidence = new Evidence(true);
    incompatibleEvidence = new Evidence(true);
    double sf;
    list<HintGroup>::iterator git;
    for (git = groupList->begin(); git != groupList->end(); git++) {
	if(!(git->getEnd() < geneBegin() || git->getBegin() > geneEnd())) {
	    sf = supportingFraction(&*git);
	    if (sf >= 1.0) {
		supportingEvidence->add(git->getSource(), git->getName());
	    } else { // TODO: in case of CDS only predictions the rest does not need to be incompatible neccessarily
		incompatibleEvidence->add(git->getSource(), git->getName());
	    }
	    addSupportedStates(&*git);
	}
    }
    CDSexonEvidence = stateEvidenceSummary(exons);
    CDSintronEvidence = stateEvidenceSummary(introns);
    UTR5stateEvidence = stateEvidenceSummary(utr5exons, utr5introns);
    UTR3stateEvidence = stateEvidenceSummary(utr3exons, utr3introns);
}


/*
 * supportingFraction
 * Compute an estimate of how well the HintGroup is compatible with the transcript.
 * The percentage of certain hints types (intron, exon, exonpart, CDS, CDSpart, ass, dss, nonirpart, nonexonpart) that is compatible with the transcript.
 */
double Gene::supportingFraction(HintGroup *group){
    list<Feature*> *hints = group->getHints();
    int supporting = 0;
    int total = 0;
    Feature *hint;
    State *state, *last5utrexon = NULL, *last3utrexon = NULL, *laststate;
    bool supports;
    if (!hints)
	return 0.0;
    for (list<Feature*>::iterator hit = hints->begin(); hit != hints->end(); hit++) {
	supports = false;
	hint = *hit;
	if ((hint->type == nonirpartF) && geneBegin() <= hint->start && geneEnd() >= hint->end)
	    supports = true;
	for (state = exons; state != NULL; state = state->next){
	    if ((hint->type == exonF || hint->type == CDSF) && hint->start == state->begin && hint->end == state->end)
		supports = true;
	    else if ((hint->type == exonpartF || hint->type == CDSpartF) && hint->start >= state->begin && hint->end <= state->end
		     && state->frame_compatible(hint))
	      supports = true;
	}
	for (state = introns; state != NULL; state = state->next){
	    if (hint->type == intronF && hint->start == state->begin && hint->end == state->end)
		supports = true;
	    else if ((hint->type == intronpartF || hint->type == nonexonpartF) && hint->start >= state->begin && hint->end <= state->end)
		supports = true;
	    else if ((hint->type == assF || hint->type == dssF) && ((hint->start <= state->begin && hint->end >= state->begin) || (hint->start <= state->end && hint->end >= state->end)))
		  supports = true;
	}
	laststate = NULL;
	for (state = utr5exons; state != NULL; state = state->next){
	    if ((hint->type == exonF || hint->type == UTRF) && hint->start == state->begin && hint->end == state->end)
		supports = true;
	    else if ((hint->type == exonpartF || hint->type == UTRpartF) && hint->start >= state->begin && hint->end <= state->end)
		supports = true;
	    else if ((hint->type == assF || hint->type == dssF) && 
		((state->next && (hint->start <= state->end+1 && hint->end >= state->end+1)) ||
		 (state != utr5exons && (hint->start <= state->begin-1 && hint->end >= state->begin-1))))
		supports = true;
	    else if (hint->type == intronF && laststate && laststate->end+1 == hint->start && state->begin-1 == hint->end)
  	        supports = true;
	    else if ((hint->type == intronpartF || hint->type == nonexonpartF) && laststate && laststate->end+1 <= hint->start && state->begin-1 >= hint->end)
	        supports = true;
	    last5utrexon = state;
	    laststate = state;
	}
	laststate = NULL;
	for (state = utr3exons; state != NULL; state = state->next){
	    if ((hint->type == exonF || hint->type == UTRF) && hint->start == state->begin && hint->end == state->end)
		supports = true;
	    else if ((hint->type == exonpartF || hint->type == UTRpartF) && hint->start >= state->begin && hint->end <= state->end)
		supports = true;
	    else if ((hint->type == assF || hint->type == dssF) && 
		((state->next && (hint->start <= state->end+1 && hint->end >= state->end+1)) ||
		 (state != utr3exons && (hint->start <= state->begin-1 && hint->end >= state->begin-1))))
		supports = true;
	    else if (hint->type == intronF && laststate && laststate->end+1 == hint->start && state->begin-1 == hint->end)
		supports = true;
	    else if ((hint->type == intronpartF || hint->type == nonexonpartF) && laststate && laststate->end+1 <= hint->start && state->begin-1 >= hint->end)
	        supports = true;
	    last3utrexon = state;
	    laststate = state;
	}
	if (hint->type == exonF || hint->type == exonpartF){
	    int exonbegin=-1, exonend=-1;
	    if (exons && exons->next == NULL) {// single CDS case
		// exon including the start codon, forward strand
		if (strand == plusstrand && last5utrexon && utr3exons){
		    exonbegin = last5utrexon->begin;
		    exonend = utr3exons->end;
		}
		// exon including the start codon, reverse strand
		if (strand == minusstrand && last3utrexon && utr5exons){
		    exonbegin = last3utrexon->begin;
		    exonend = utr5exons->end;	    
		}
	    }
	    // exon including the start codon, multi CDS exon case, forward strand
	    if (strand == plusstrand && last5utrexon && exons && exons->next) {
		exonbegin = last5utrexon->begin;
		exonend = exons->end;
	    }
	    // exon including the start codon, multi CDS exon case, reverse strand
	    if (strand == minusstrand && exons && exons->next && utr5exons){
		exonbegin = lastExon()->begin;
		exonend = utr5exons->end;
	    }
	    if (exonbegin>0 && exonend>0 && hint->type == exonF && hint->start == exonbegin && hint->end == exonend)
		supports = true;
	    if (exonbegin>0 && exonend>0 && hint->type == exonpartF && hint->start >= exonbegin && hint->end <= exonend)
		supports = true;
	    // exon including the stop codon, multi CDS exon case, forward strand
	    if (strand == plusstrand && utr3exons && exons){
		exonbegin = lastExon()->begin;
		exonend = utr3exons->end;
	    }
	    // exon including the stop codon, multi CDS exon case, reverse strand
	    if (strand == minusstrand && last3utrexon && exons){
		exonbegin = last3utrexon->begin;
		exonend = exons->end;
	    }
	    if (exonbegin>0 && exonend>0 && hint->type == exonF && hint->start == exonbegin && hint->end == exonend)
		supports = true;
	    if (exonbegin>0 && exonend>0 && hint->type == exonpartF && hint->start >= exonbegin && hint->end <= exonend)
		supports = true;	
	}
	if (hint->type == exonF || hint->type == exonpartF || hint->type == CDSF || hint->type == CDSpartF || hint->type == intronF 
	    || hint->type == intronpartF || hint->type == assF || hint->type == dssF || hint->type == UTRF || hint->type == UTRpartF || hint->type == nonirpartF || hint->type == nonexonpartF) {
	    total++;
	    if (supports)
		supporting++;
	}
    }
    if (total > 0.0)
	return (double) supporting / total;
    else
	return 0.0;
}

/*
 * Gene::addSupportedStates
 */
void Gene::addSupportedStates(HintGroup *group){
    if (!group)
	return;
    State *state;
    bool supported;
    bool contradicted;
    Feature *hint;
    string srcname = group->getSource();
    list<Feature*>::iterator hit;
    list<Feature*> *hints = group->getHints();
    if (!hints)
	return;
    if (group->getEnd() < geneBegin() || group->getBegin() > geneEnd())
	return;
    for (state = exons; state != NULL; state = state->next){
	supported = contradicted = false;
	for (hit = hints->begin(); hit != hints->end(); hit++) {
	    hint = *hit;
	    if ((hint->type == exonF || hint->type == CDSF) && hint->start == state->begin && hint->end == state->end)
		supported = true;
	    else if ((hint->type == exonpartF || hint->type == CDSpartF) && hint->start >= state->begin && hint->end <= state->end)
		supported = true;
	    else if ((hint->type == intronpartF || hint->type == intronF || hint->type == UTRF || hint->type == UTRpartF) && !(hint->start > state->end || hint->end < state->begin))
		contradicted = true;
	    // exceptions for the CDS at both ends
	    if (state == exons && hint->type == exonF && hint->end == state->end && hint->start < state->begin)
		supported = true;
	    if (state == exons && hint->type == exonpartF && hint->end <= state->end && hint->end >= state->begin)
		supported = true;
	    if (state->next == NULL && hint->type == exonF && hint->start == state->begin && hint->end >= state->end)
		supported = true;
	    if (state->next == NULL && hint->type == exonpartF && hint->start <= state->end && hint->start >= state->begin)
		supported = true;
	    if (state == exons && state->next == NULL && (hint->type == exonF || hint->type == exonpartF) && hint->start <= state->begin && hint->end >= state->end)
		supported = true;
	}
	if (supported && !contradicted)
	    state->addEvidence(srcname);
    }
    State *intronTypes[3];
    intronTypes[0] = introns;
    intronTypes[1] = utr5introns;
    intronTypes[2] = utr3introns;
    for (int i=0; i<3; i++) { // loop over the 3 intron types: CDS intron, 5'UTR intron and 3'UTR intron
	for (state = intronTypes[i]; state != NULL; state = state->next){
	    supported = contradicted = false;
	    for (hit = hints->begin(); hit != hints->end(); hit++) {
		hint = *hit;
		if (hint->type == intronF && hint->start == state->begin && hint->end == state->end)
		    supported = true;
		else if (hint->type == intronpartF && hint->start >= state->begin && hint->end <= state->end)
		    supported = true;
		else if ((hint->type == exonpartF || hint->type == exonF || hint->type == UTRF || hint->type == UTRpartF) && !(hint->start > state->end || hint->end < state->begin))
		    contradicted = true;
	    }
	    if (supported && !contradicted)
		state->addEvidence(srcname);
	}
    }
    for (state = utr5exons; state != NULL; state = state->next){
	supported = contradicted = false;
	for (hit = hints->begin(); hit != hints->end(); hit++) {
	    hint = *hit;
	    if ((hint->type == exonF || hint->type == UTRF) && hint->start == state->begin && hint->end == state->end)
		supported = true;
	    else if ((hint->type == UTRpartF || hint->type == exonpartF)&& hint->start >= state->begin && hint->end <= state->end)
		supported = true;
	    else if ((hint->type == intronpartF || hint->type == intronF || hint->type == CDSF || hint->type == CDSpartF) && !(hint->start > state->end || hint->end < state->begin))
		contradicted = true;
	    if (hint->type == exonF && (( strand == plusstrand && state->next == NULL && hint->start == state->begin && hint->end >= state->end) ||
					( strand == minusstrand && state == utr5exons && hint->end == state->end && hint->start <= state->begin)))
		supported = true;
	    if (hint->type == exonpartF && ((strand == plusstrand && state->next == NULL && hint->start >= state->begin && hint->start <= state->end) ||
					    (strand == minusstrand && state == utr5exons && hint->end >= state->begin && hint->end <= state->end)))
		supported = true;
	}
	if (supported && !contradicted)
	    state->addEvidence(srcname);
    }
    for (state = utr3exons; state != NULL; state = state->next){
	supported = contradicted = false;
	for (hit = hints->begin(); hit != hints->end(); hit++) {
	    hint = *hit;
	    if ((hint->type == exonF || hint->type == UTRF) && hint->start == state->begin && hint->end == state->end)
		supported = true;
	    else if ((hint->type == UTRpartF || hint->type == exonpartF)&& hint->start >= state->begin && hint->end <= state->end)
		supported = true;
	    else if ((hint->type == intronpartF || hint->type == intronF || hint->type == CDSF || hint->type == CDSpartF) && !(hint->start > state->end || hint->end < state->begin))
		contradicted = true;
	    if (hint->type == exonF && ((strand == plusstrand && state == utr3exons && hint->end == state->end && hint->start <= state->end)||
					(strand == minusstrand && state->next == NULL && hint->start == state->begin && hint->end >= state->end)))
		supported = true;
	    if (hint->type == exonpartF && ((strand == plusstrand && state == utr3exons && hint->end >= state->begin && hint->end <= state->end) ||
					    (strand == minusstrand && state->next == NULL && hint->start >= state->begin && hint->start <= state->end))) 
		supported = true;
	}
	if (supported && !contradicted)
	    state->addEvidence(srcname);
    }
}

/*
 * Gene::getPercentSupported
 * in [0,1]
 */
double Gene::getPercentSupported(){
    int numStates = lenStateList(exons) + lenStateList(introns) + lenStateList(utr5exons) + lenStateList(utr5introns) + lenStateList(utr3exons) + lenStateList(utr3introns);
    if (numStates <= 0)
	return 0.0;
    int numSupported = 0;
    if (CDSexonEvidence)
	numSupported += CDSexonEvidence->numEvidence;
    if (CDSintronEvidence)
	numSupported += CDSintronEvidence->numEvidence;
    if (UTR5stateEvidence)
	numSupported += UTR5stateEvidence->numEvidence;
    if (UTR3stateEvidence)
	numSupported += UTR3stateEvidence->numEvidence;
    return ((double) numSupported) / numStates;
}

/*
 * Gene::getCDSCoord (start counting with 0)
 * if comp is true, start counting from the downstream end of the CDS
 * if CDS is incomplete, count only nucleotides of the predicted CDS
 * if loc is not inside a CDS, return -1
 */
int Gene::getCDSCoord(int loc, bool comp) const {
    if (exons == NULL || exons->begin > loc)
	return -1;
    State* current = exons;
    int result = 0;
    if (comp) {
	while (loc > current->end) {
	    current = current->next;
	    if (current == NULL)
		return -1;
	}
	result = loc < current->begin ? 
	    current->end - current->begin + 1 :
	    current->end - loc;
	while ((current = current->next) != NULL) 
	    result += (current->end - current->begin + 1);
	return result;
    } else
	do {
	    if (loc < current->begin)
		return -1;
	    if (loc <= current->end)
		return result + loc - current->begin;
	    result += current->length();
	    current = current->next;
	} while (current != NULL);
    return -1;
}
    
/*
 * Gene::completeCDS
 * checks whether it has a single or initial and terminal CDS and the first and last coding exons are not truncated
 */
bool Gene::completeCDS() const {
    if (exons == NULL)
	return false;
    if (isInternalExon(exons->type) || isRInternalExon(exons->type) || 
	isInternalExon(lastExon()->type) || isRInternalExon(lastExon()->type) ||
	exons->type == terminal || isRTerminalExon(lastExon()->type) ||
	exons->type == rinitial || isInitialExon(lastExon()->type))
	return false; // at least one complete exon is missing
    if (exons->truncated & TRUNC_LEFT || lastExon()->truncated & TRUNC_RIGHT )
	return false; // an exon is truncated at the boundaries
    return true;
}

/*
 * Gene::operator<
 * Genes are sorted by transcription start.
 */
bool Gene::operator< (const Gene &other) const {
    if (geneBegin() < other.geneBegin()) 
	return true;
    return false;
}

bool Gene::operator== (const Gene &other) const {
    // check equality of coding exons, utrexons, introns
    State *st, *otherst, *sv[4],*osv[4];
    sv[0]=exons;
    sv[1]=introns;
    sv[2]=utr5exons;
    sv[3]=utr3exons;
    osv[0]=other.exons;
    osv[1]=other.introns;
    osv[2]=other.utr5exons;
    osv[3]=other.utr3exons;

    for (int i=0; i<4; i++){
	st = sv[i];
	otherst = osv[i];
	if ((!st && otherst) || (st && !otherst))
	    return false;
	while (st) {
	    if (st->begin != otherst->begin || st->end != otherst->end)
		return false;
	    st = st->next;
	    otherst = otherst->next;
	    if ((!st && otherst) || (st && !otherst))
		return false;
	}
    }
    return true;
}

void Gene::print(){
    State* curExon;
    cout << "sequence: " << seqname << " ,  gene: " << id << endl;
    for(curExon = exons; curExon != NULL; curExon = curExon->next) {
	cout << stateTypeNames[curExon->type] << ":  " << curExon->begin << " .. " << curExon->end << endl;
    }
}

void Gene::printGFF() const{
    State* curExon, *lastExon = NULL, *curIntron;
    State *first_right_utr = (strand == plusstrand)? utr3exons : utr5exons,
      *first_left_utr = (strand == plusstrand)? utr5exons : utr3exons;
    int from, to;
    string transcript_id = geneid + "." + id ;
    string parentstr = gff3 ? 
	("Parent=" + transcript_id) :
	("transcript_id \"" + transcript_id + "\"; gene_id \"" + geneid + "\";");

    // output left utr
    for(curExon = first_left_utr; curExon != NULL; curExon = curExon->next) {
	if (strand == plusstrand && curExon == utr5exons && complete5utr && print_tss) {
	    cout << seqname << "\t" << source << "\t";
	    if (gff3)
		cout << "transcription_start_site\t";
	    else
		cout << "tss\t" ;
	    cout << curExon->begin+1 << "\t" << curExon->begin+1 << "\t";
	    cout << "."; // no apostprob of signals yet
	    cout << "\t+\t.\t" << parentstr << endl;
	}
	if (strand == minusstrand && curExon == utr3exons && complete3utr && print_tts) {
	    cout << seqname << "\t" << source << "\t";
	    if (gff3)
		cout << "transcription_end_site\t";
	    else
		cout << "tts\t";
	    cout << curExon->begin+1 << "\t" << curExon->begin+1 << "\t";
	    cout << "."; // no apostprob of signals yet
	    cout << "\t-\t.\t" << parentstr << endl;
	}	
	if (print_utr) { 
	    if (curExon->end >= curExon->begin) {// UTR exon can have length 0 when start codon comes right after splice site
		cout << seqname << "\t" << source << "\t";
		if ((strand == plusstrand))
		    if (gff3)
			cout << "five_prime_utr";
		    else 
			cout << "5'-UTR";
		else 
		    if (gff3)
			cout << "three_prime_utr";
		    else 
			cout << "3'-UTR";
		cout << "\t" << curExon->begin+1 << "\t" << curExon->end+1 << "\t";
		if (curExon->hasScore)
		    cout << setprecision(3) << curExon->apostprob; // score
		else
		    cout << ".";
		cout << "\t" << ((strand == plusstrand)? '+' : '-') << "\t"; // strand
		cout << ".\t" << parentstr << endl;
	    }
	} else {
	    from = curExon->begin+1;
	    to = curExon->end+1; 
	    if (!curExon->next) {// last left utr exon
		if (exons) {
		    to = exons->end+1;
		    if (exons && !exons->next && first_right_utr)
			to = first_right_utr->end+1;
		}
	    }
	    cout << seqname << "\t" << source << "\t";
	    cout << "exon\t" << from << "\t" << to << "\t.";
	    cout << "\t" << ((strand == plusstrand)? '+' : '-') << "\t"; 
	    cout << ".\t" << parentstr << endl;
	}
   }

    // TODO: im singlestrand-modus werden bei start- und stop codons auf dem minusstrang
    // die Stränge als + angegeben. 
    if (exons) {
	// output forward strand start codon
	if (print_start && (strand == plusstrand) && (isInitialExon(exons->type) || exons->type == singleG)){
	    cout << seqname << "\t" << source << "\t" << "start_codon" << "\t"
		 << exons->begin + 1 << "\t" << exons->begin + 3
		 << "\t.\t+\t0\t" << parentstr << endl;
	}
	// output reverse strand stop codon
	if (print_stop && (strand == minusstrand) && 
	    (exons->type == terminal || exons->type == singleG || isRTerminalExon(exons->type) || exons->type == rsingleG)){
	    cout << seqname << "\t" << source << "\t" << "stop_codon" << "\t"
		 << exons->begin + 1 << "\t" << exons->begin + 3
		 << "\t.\t-\t0\t" << parentstr << endl;
	}
    }

    for(curExon = exons; curExon != NULL; curExon = curExon->next) {

	// print 'single', 'initial', 'terminal' or 'internal'
	if (print_exonnames && !gff3) {
	    cout << seqname << "\t" << source << "\t";
	    if (curExon->type == singleG || curExon->type == rsingleG)
		cout << "single";
	    else if (isInitialExon(curExon->type) || curExon->type == rinitial)
		cout << "initial";
	    else if (curExon->type == terminal || isRTerminalExon(curExon->type))
		cout << "terminal";
	    else 
		cout << "internal";
	    cout << "\t" 
		 << curExon->begin+1 << "\t" << curExon->end+1 << "\t";
	    if (curExon->hasScore)
		cout << setprecision(3) << curExon->apostprob; // score
	    else
		cout << ".";
	    cout << "\t" << ((strand == plusstrand)? '+' : '-') << "\t"; // strand
	    // reading frame, here: #nucleotides in 5' before first complete codon
	    if (strand == plusstrand)
		cout << mod3(3-(curExon->frame() - curExon->length()));
	    else
		cout << mod3(2-curExon->frame());
	    cout << "\t" << "transcript_id \"" << geneid << "." << id << "\"; gene_id " << "\"" << geneid << "\";" << endl;
	}

	lastExon = curExon;

    }
    // output introns
    if (print_introns)
	for (curIntron = introns; curIntron != NULL; curIntron = curIntron->next) {
	    cout << seqname << "\t" << source << "\t" << "intron";
	    cout << "\t" << curIntron->begin+1 << "\t" << curIntron->end+1 << "\t";
	    if (curIntron->hasScore)
		cout << setprecision(3) << curIntron->apostprob; // score
	    else
		cout << ".";
	    cout << "\t" << ((strand == plusstrand)? '+' : '-') << "\t" // strand
		 << "." << "\t"; // reading frame
	    cout << parentstr << endl;
	}

    /* 
     * CDS exons
     */

    for(curExon = exons; curExon != NULL; curExon = curExon->next) {
	/* 
	 * print CDS
	 */
	if (print_cds) {
	    int beginmod=0, endmod=0;
	    if (stopCodonExcludedFromCDS) {
		if (curExon->type == terminal || curExon->type == singleG)
		    endmod = -3;
		if (isRTerminalExon(curExon->type) || curExon->type == rsingleG)
		    beginmod = 3;
	    }
	    if ( curExon->begin+1 + beginmod <= curExon->end+1 + endmod) {
		// this check ensures that exons consisting just of one stop codon are omitted completely, when the stop codon is excluded from CDS
		cout << seqname << "\t" << source << "\t" << "CDS" << "\t"; 
		cout << curExon->begin+1 + beginmod << "\t" << curExon->end+1 + endmod << "\t";
		if (curExon->hasScore)
		    cout << setprecision(3) << curExon->apostprob; // score
		else
		    cout << ".";
		cout << "\t" << ((strand == plusstrand)? '+' : '-') << "\t"; // strand
		if (strand == plusstrand)
		    cout << mod3(3-(curExon->frame() - curExon->length())) << "\t";
		else
		    cout << mod3(2-curExon->frame()) << "\t";
		if (gff3) 
		    cout << "ID=" << geneid << "." << id << ".cds" << ";";
		cout << parentstr << endl;
	    }	
	}
	// print 'exon'
	if (Constant::utr_option_on && !print_utr) { // output 'exon' only if prediction with UTR is turned on and UTRs are not output in the other (UTR) format
	  if (curExon != exons || !first_left_utr) {// first_left_utr can be NULL in case of incomplete gene
		from = curExon->begin+1;
		to = curExon->end+1;
		if (!curExon->next && first_right_utr)
		    to = first_right_utr->end+1;
		cout << seqname << "\t" << source << "\t";
		cout << "exon\t" << from << "\t" << to << "\t.";
		cout << "\t" << ((strand == plusstrand)? '+' : '-') << "\t"; 
		cout << ".\t" << parentstr << endl;
	    }
	}

    }
    
    if (lastExon){
	// output forward strand stop codon
	if (print_stop && (strand == plusstrand) && (lastExon->type == terminal || lastExon->type == singleG)){
	    cout << seqname << "\t" << source << "\t" << "stop_codon" << "\t"
		 << lastExon->end - 1 << "\t" << lastExon->end + 1
		 << "\t.\t+\t0\t" << parentstr << endl;
	}
	// output reverse strand start codon
	if (print_start && (strand == minusstrand) && 
	    (isInitialExon(lastExon->type) || lastExon->type == singleG || lastExon->type == rinitial || lastExon->type == rsingleG)){
	    cout << seqname << "\t" << source << "\t" << "start_codon" << "\t"
		 << lastExon->end - 1 << "\t" << lastExon->end + 1
		 << "\t.\t-\t0\t" << parentstr << endl;
	}
	
    }
    // output right utr
    for(curExon = first_right_utr; curExon != NULL; curExon = curExon->next) {
	if (print_utr) {
	    if (curExon->end >= curExon->begin) {
		cout << seqname << "\t" << source << "\t";
		if ((strand == plusstrand))
		    if (gff3)
			cout << "three_prime_utr";
		    else
			cout << "3'-UTR";
		else 
		    if (gff3)
			cout << "five_prime_utr";
		    else 
			cout << "5'-UTR";
		cout << "\t" << curExon->begin+1 << "\t" << curExon->end+1 << "\t";
		if (curExon->hasScore)
		    cout << setprecision(3) << curExon->apostprob; // score
		else
		    cout << ".";
		cout << "\t" << ((strand == plusstrand)? '+' : '-') << "\t"; // strand
		cout << ".\t" << parentstr << endl;
	    }
	} else {
	    if (curExon != first_right_utr) {
		from = curExon->begin+1;
		to = curExon->end+1;
		cout << seqname << "\t" << source << "\t";
		cout << "exon\t" << from << "\t" << to << "\t.";
		cout << "\t" << ((strand == plusstrand)? '+' : '-') << "\t"; 
		cout << ".\t" << parentstr << endl;
	    }
	}
	if (!curExon->next && strand == plusstrand && complete3utr && print_tts) {
	    cout << seqname << "\t" << source << "\t";
	    if (gff3)
		cout << "transcription_end_site\t";
	    else
		cout << "tts\t";
	    cout << curExon->end+1 << "\t" << curExon->end+1 << "\t";
	    cout << "."; // no apostprob of signals yet
	    cout << "\t+\t.\t" << parentstr << endl;
	}
	if (!curExon->next && strand == minusstrand && complete5utr && print_tss) {
	    cout << seqname << "\t" << source << "\t";
	    if (gff3)
		cout << "transcription_start_site\t";
	    else
		cout << "tss\t";
	    cout << curExon->end+1 << "\t" << curExon->end+1 << "\t";
	    cout << "."; // no apostprob of signals yet
	    cout << "\t-\t.\t" << parentstr << endl;
	}	
    }

    // output protein matches
    list<PP::Match>::const_iterator it = proteinMatches.begin();
    if (it == proteinMatches.end()) 
	return;

    string match_id = "ID=pp." + transcript_id;
 
    // output blocks
    while (it != proteinMatches.end()) {
	int targetStart = ((strand == plusstrand ? 
			   getCDSCoord(it->firstBase, false) : 
			   getCDSCoord(it->lastBase, true)) - it->offset)/3;
#ifdef DEBUG
	if (targetStart == -1)
	    throw ProjectError("match outside of exon!");
#endif
	cout << seqname << "\t" << source << "\t" << "protein_match" << "\t"
	     << (it->firstBase+1) << "\t" << (it->lastBase+1) << "\t"
	     << it->score << "\t"  
	     << (it->complement? "-" : "+") << "\t" << it->phase() << "\t";
	if (gff3)
	    cout << match_id << "." << it->blockId << ";Target=" << it->blockId << " " << (it->firstCodon()+1) << " " << (it->lastCodon()+1)
		 << ";Target_start=" << targetStart << ";" << endl;   // no "Parent=" as GBrowse is bothered by this
	else
	    cout << "target \"" << it->blockId << "[" << (it->firstCodon()+1) << ".." << (it->lastCodon()+1) << "]\"; "
		 << "target_start " << targetStart << "; " << parentstr << endl;
	++it; 
   }

    // for GBrowse, output interblock regions (only if gff3 is turned on)
    if (gff3) {
	it = proteinMatches.begin();
	int blockno = (strand==plusstrand? it->blockno : it->blockno+1);
	for (curExon = exons; curExon != NULL; curExon = curExon->next) {
	    int firstBase = curExon->begin;
	    int phase = (strand == plusstrand) ? 
		mod3(3-(curExon->frame() - curExon->length())) : 0;
	    while (it != proteinMatches.end() && it->firstBase <= curExon->end) {
		if (it->firstBase > firstBase) {
#ifdef DEBUG
		    if ((strand==plusstrand && blockno != it->blockno) || (strand==minusstrand && blockno != it->blockno+1))
			throw ProjectError("gene.cc: Gene::printGFF(): error printing interblock regions!");
#endif
		    cout << seqname << "\t" << source << "\t" << "interblock_region" << "\t"
			 << (firstBase+1) << "\t" << (it->firstBase) << "\t" << "." << "\t"
			 << (strand == plusstrand? '+' : '-') << "\t" << phase << "\t" 
			 << match_id << ".iBR" << blockno  << "\n";
		}
		firstBase = it->lastBase+1; phase=0;
		blockno = (strand == plusstrand ? it->blockno + 1 : it->blockno);
		++it;
	    }
	    if (firstBase <= curExon->end)
		cout << seqname << "\t" << source << "\t" << "interblock_region" << "\t"
		     << (firstBase+1) << "\t" << (curExon->end + 1) << "\t" << "." << "\t"
		     << (strand == plusstrand? '+' : '-') << "\t" 
		     << (strand == plusstrand? phase : mod3(2-curExon->frame())) << "\t" 
		     << match_id << ".iBR" << blockno << "\n";
	}
    }
    
}

void Gene::printCodingSeq(AnnoSequence *annoseq) const {
    int linelength=100;
    int curlinelength=0;
    string codingSeq = getCodingSequence(annoseq);
    int offset = 0;    
    cout << "# coding sequence = [";
    curlinelength += 21;

    while (offset < codingSeq.length()){
      cout << codingSeq.substr(offset, linelength-curlinelength);
      offset += linelength-curlinelength;
      if (offset < codingSeq.length()) {
	cout << endl << "# ";
	curlinelength = 2;
      }

    }
    cout << "]" << endl;
}

// TODO: move this to GeneticCode
string getTranslation(const char* codingSeq) {
    string result="";
    while (codingSeq[0] && codingSeq[1] && codingSeq[2]) {
	try {
	    char aa = GeneticCode::translate(codingSeq);
	    // for backward compatibily, replace '*' by 'X'
	    if (aa != '*') 
		result.append(1,aa);
	    else if (codingSeq[3])
		result.append("X");
	} catch (...) {
	    result.append("X"); // from masking
	}
	codingSeq += 3;
    }
    return result;
}

void Gene::printProteinSeq(AnnoSequence *annoseq) const {
    int linelength=100;
    string prefix = "# protein sequence = [";
    
    // if gene is incomplete, frame is the position of the first base
    string trans = getTranslation(getCodingSequence(annoseq) + mod3(-frame)); 
    int i = linelength - prefix.length();
    cout << prefix << trans.substr(0,i);
    while (i<trans.length()) {
	cout << "\n# " << trans.substr(i, linelength - 2);
	i += (linelength -2);
    }
    cout << "]" << endl;
}

void Gene::printBlockSequences(AnnoSequence *annoseq) const {
    list<PP::Match>::const_iterator it = proteinMatches.begin(), it2;
    while (it != proteinMatches.end()) {
	int targetStart;
	string codingseq="";
	string protseq="";
	for (it2 = it; it2 != proteinMatches.end() && it2->blockno == it->blockno; it2++) {
	    codingseq.append(annoseq->sequence - annoseq->offset + it2->firstBase, 
			      it2->lastBase - it2->firstBase + 1);
	}
	--it2;
	if (strand == plusstrand) {
	    targetStart = getCDSCoord(it->firstBase, false)/3;
	    protseq = getTranslation(codingseq.c_str() + mod3(-it->offset));
	} else  {
	    targetStart = getCDSCoord(it2->lastBase, true)/3;
	    char *res = reverseComplement(codingseq.c_str());
	    protseq = getTranslation(res + mod3(-it2->offset));
	    delete[] res;
	}
	cout << "# sequence of block " << left << setw(10) << it->blockId
	     << right << setw(5) << targetStart << " [" << protseq << "] " 
	     << targetStart + protseq.length() << "\n";
	it = ++it2;
    }
}
    
void Gene::printEvidence() {
    cout << "# Evidence for and against this transcript:" << endl;
    // compatibility by transcript
    int numCDS=0, numCDSintron=0, num5UTR=0, num3UTR=0;
    if (CDSexonEvidence || CDSintronEvidence || UTR5stateEvidence || UTR3stateEvidence) {
	numCDS = lenStateList(exons);
	numCDSintron = lenStateList(introns);
	num5UTR = lenStateList(utr5exons) + lenStateList(utr5introns);
	num3UTR = lenStateList(utr3exons) + lenStateList(utr3introns);
	cout << "# % of transcript supported by hints (any source): " << setprecision(3) 
	     << 100.0 * getPercentSupported() <<  endl;
    }
    if (CDSexonEvidence && CDSintronEvidence){
	cout << "# CDS exons: " << CDSexonEvidence->numEvidence << "/" << numCDS << endl;
	CDSexonEvidence->print();
	cout << "# CDS introns: " << CDSintronEvidence->numEvidence << "/" << numCDSintron << endl;
	CDSintronEvidence->print();
    }
    if (UTR5stateEvidence && UTR3stateEvidence) {
	cout << "# 5'UTR exons and introns: " << UTR5stateEvidence->numEvidence << "/" << num5UTR << endl;
	UTR5stateEvidence->print();
	cout << "# 3'UTR exons and introns: " << UTR3stateEvidence->numEvidence << "/" << num3UTR << endl;
	UTR3stateEvidence->print();
    }   
    // compatibility by hints
    if (supportingEvidence){
	cout << "# hint groups fully obeyed: " << supportingEvidence->numEvidence << endl;
	supportingEvidence->print();
    }
    if (incompatibleEvidence){
	cout << "# incompatible hint groups: " << incompatibleEvidence->numEvidence << endl;
	incompatibleEvidence->print();
    }
}

void Gene::init() {
    Properties::assignProperty("start", print_start);
    Properties::assignProperty("stop", print_stop);
    Properties::assignProperty("tss", print_tss);
    Properties::assignProperty("tts", print_tts);
    Properties::assignProperty("introns", print_introns);
    Properties::assignProperty("gff3", gff3);
    Properties::assignProperty("stopCodonExcludedFromCDS", stopCodonExcludedFromCDS);
    Properties::assignProperty("cds", print_cds);
    Properties::assignProperty("exonnames", print_exonnames);
    Properties::assignProperty("print_utr", print_utr);
    Properties::assignProperty("print_blocks", Gene::print_blocks);
}

list<Gene> *Gene::filterGenePrediction(list<Gene> *gl, const char *seq, Strand strand, bool noInFrameStop, double minmeanexonintronprob, double minexonintronprob){
    State *s;
    list<Gene> *filteredTranscripts = new list<Gene>;
    AnnoSequence *annoseq = new AnnoSequence();
    annoseq->sequence = newstrcpy(seq);
    annoseq->offset = 0;

    for(list<Gene>::iterator git = gl->begin();git != gl->end();git++){
	bool keep = true;
	// delete gene if the combined CDS is too short, unless a CDS exon is truncated
	if (git->clength < Constant::min_coding_len && git->completeCDS())
	    keep = false;
	if (keep && !(strand == bothstrands || git->strand == strand))
	    keep = false;
	if (keep && git->throwaway)
	    keep = false;
    
	if (keep && git->hasInFrameStop(annoseq) && noInFrameStop){
	    keep = false;
	}
	if (keep && git->hasProbs) {
	    /* 
	     * filter by a posteriori probability
	     * keep gene if the geometric mean of all states posterior probabilities is at least minmeanexonintronprob
	     */
	    if (git->meanStateProb() < minmeanexonintronprob && !(Constant::keep_viterbi && git->viterbi))
		keep = false;
	}
	if (keep && git->hasProbs) {
	    /* 
	     * filter transcript by a posteriori probability of exons,
	     * keep gene if all apostprobs of coding exons are at least minexonintronprob
	     */
	    for (s = git->exons; s!=NULL; s=s->next){
		if (s->apostprob == 0.0)
		    throw ProjectError("Gene::filterGenePrediction:Detected state with posterior probability 0.");
		if (s->apostprob < minexonintronprob && !(Constant::keep_viterbi && git->viterbi))
		    keep = false;
	    } 
	    for (s = git->introns; s!=NULL; s=s->next){
		if (s->apostprob == 0.0)
		    throw ProjectError("Gene::filterGenePrediction:Detected state with posterior probability 0.");
		if (s->apostprob < minexonintronprob && !(Constant::keep_viterbi && git->viterbi))
		    keep = false;
	    }
	}

	if (keep) {
	    filteredTranscripts->push_back(*git);
	} else {
#ifdef DEBUG
	    // cerr << "Gene deleted! (clength " << git->id << ")" << endl;
#endif
	}
    }
    delete annoseq;
    return filteredTranscripts;
}

/*
 * filterTranscriptsByMaxTracks
 * 
 * With the priority given by the order of gl, delete transcripts that overlap more than maxTracks times
 * at the same position. Purpose: We then need only maxTracks for display in a genome browser.
 */
void Gene::filterTranscriptsByMaxTracks(list<Gene> *gl, int maxTracks){
    if (maxTracks < 0)
	maxTracks = INT_MAX;

    /*
     * sort from gl into intermediate gene list sorted by meanStateProb
     */
    list<Gene> *sorted = new list<Gene>;
    list<Gene>::iterator it, largestit;
    double maxMeanProb, mMP;
    while(!gl->empty()){
	// find transcript with largest meanStateProb among the rest
	// TODO: sort by a combination of extrinsic evidence and posterior probabilities
	maxMeanProb = -1.0;
	for (it = gl->begin(); it != gl->end();it++) {
	    mMP = it->meanStateProb();
	    if (mMP > maxMeanProb){
		maxMeanProb = mMP;
		largestit = it;
	    }
	    if (it->viterbi && Constant::keep_viterbi) {
		largestit = it;
		maxMeanProb = 1.0;
	    }
	}
	sorted->push_back(*largestit);
	gl->erase(largestit);
    }

    /*
     * further stable sort via the following criteria
     * for genes with the same start and end point the posterior transcript prob. defines the sorting
    
    cout << "vor Nachsortieren" << endl;
    for (it = sorted->begin(); it != sorted->end();it++) {
	cout << it->id << "\t" << it->geneBegin() << "-" << it->geneEnd() << "\t" << it->meanStateProb() << "\t" << it->apostprob << endl;
    }
    Gene *tempg;
    list<Gene>::iterator it2, tempit;
    for (it = sorted->begin(); it != sorted->end(); it++) {
	// exchange it with the most likely with the same begin and end point
	for (it2 = it, it2++; it2 != sorted->end(); it2++) {
	    if (it->geneBegin() == it2->geneBegin() && it->geneEnd() == it2->geneEnd() &&
		it->apostprob < it2->apostprob) {
		// swap it and it2
		tempg = new Gene(*it);
		tempit = it;
		tempit++;
		sorted->erase(it);
		it = sorted->insert(tempit, *it2);
		tempit = it2;
		tempit++;
		sorted->erase(it2);
		it2 = sorted->insert(tempit, *tempg);
		delete tempg;
	    }
	}
    }
    cout << "nach Nachsortieren" << endl;
    for (it = sorted->begin(); it != sorted->end();it++) {
	cout << it->id << "\t" << it->geneBegin() << "-" << it->geneEnd() << "\t" << it->meanStateProb() << "\t" << it->apostprob << endl;
    } */
    

    /*
     * Filter back into gl
     * keep a sorted list of segments with the number of overlapping transcripts in these segments
     *  -------------       ---------------
     *       ----------------------
     * 12345678901234567890123456789012345
     * is represented by the list of FreqSegments
     * (INT_MIN,0)(2,1)(7,2)(15,1)(22,2)(29,1)(INT_MAX,0)
     */
    list<FreqSegment> *fsl = new list<FreqSegment>;
    list<FreqSegment>::iterator fsit, fsit2;
    int freq = 0;
    bool tracksFull;
    fsl->push_back(FreqSegment(INT_MIN,0));
    fsl->push_back(FreqSegment(INT_MAX,0));
    for (it = sorted->begin();it != sorted->end(); it++) {
	//cout << "FreqSegmentList "; for (fsit = fsl->begin(); fsit!=fsl->end(); fsit++) cout << "(" << fsit->start << "," << fsit->freq << ")";	cout << endl;
	//cout << "inserting transcript " << it->id << "\t" << it->viterbi << "\t" << it->geneBegin() << "-" << it->geneEnd() << endl;
	// find the two segment endpoints where the left endpoint of it fits in between
	fsit = fsl->begin();
	while (fsit->start <= it->geneBegin() && fsit != fsl->end())
	    fsit++;
	fsit--;
	// check whether already maxTracks tracks are used in this area
	tracksFull = false;
	for (fsit2 = fsit; fsit2->start <= it->geneEnd();fsit2++) {
	    if (fsit2->freq >= maxTracks)
		tracksFull = true;
	}
	if (!tracksFull) { //insert new transcript only if track number is in limit
	    gl->push_back(*it);
	    if (fsit->start < it->geneBegin()){ // insert new FreqSegment at the beginning of transcript
		freq = fsit->freq;
		fsit++;
		fsit = fsl->insert(fsit, FreqSegment(it->geneBegin(), freq));
	    }
	    // add 1 to the freq of all FreqSegments up to the end of the gene
	    while (fsit->start <= it->geneEnd()){
		freq = fsit->freq;
		fsit->freq += 1;
		fsit++;
	    }
	    if (fsit->start > it->geneEnd() + 1) { // insert new FreqSegment at the end of transcript
		fsl->insert(fsit, FreqSegment(it->geneEnd()+1, freq));
	    }
	}
    }	
    //cout << "FreqSegmentList "; for (fsit = fsl->begin(); fsit!=fsl->end(); fsit++) cout << "(" << fsit->start << "," << fsit->freq << ")";	cout << endl;
    delete sorted;
}

/*
 * getGenesOnStrand
 * 
 * creates a new subsequence of the genes list with all genes on the strand 'strand'. Allocates new memory.
 */
Gene* Gene::getGenesOnStrand(Gene* genes, Strand strand){
    Gene *sublist = NULL, *lastGene = NULL;
    for (Gene *gene = genes; gene != NULL; gene = gene->next)
	if (gene->strand == strand) {
	    if (lastGene){
		lastGene->next = new Gene(*gene);
		lastGene = lastGene->next;
	    } else
		sublist = lastGene = new Gene(*gene);
	    lastGene->next = NULL;
	}
    return sublist;
}



/* --- AltGene methods --------------------------------------------- */

bool AltGene::operator< (const AltGene &other) const {
  return mincodstart < other.mincodstart;
}


void AltGene::addGene(Gene *gene) {
    if (transcripts.empty()) {
	strand = gene->strand;
	mincodstart = gene->codingstart;
	maxcodend = gene->codingend;
	transcripts.push_back(gene);
    } else {
	if (gene->strand != strand)
	    throw ProjectError("AltGene::addGene: Tried to summarize transcripts on opposite strands into gene.");
	if (gene->codingstart < mincodstart)
	    mincodstart = gene->codingstart;
	if (gene->codingend > maxcodend)
	    maxcodend = gene->codingend;
    
	bool prevExist = false;
	for(list<Gene*>::iterator git = transcripts.begin(); git != transcripts.end(); git++){
	    if (**git == *gene) {
		(*git)->addSampleCount(1);
		(*git)->addStatePostProbs(1.0);
		prevExist = true;
	    }
	}
	if (!prevExist)
	    transcripts.push_back(gene);
    }
    apostprob += gene->apostprob;
}

/*
 * A gene overlaps the transcript if they have a common coding base.
 */
bool AltGene::overlaps(Gene *gene){
  if (!gene->exons)
    return false;
  if (!((strand == gene->strand) && (gene->codingstart <= maxcodend) && (gene->codingend >= mincodstart)))
      return false;
  // check whether any coding region overlaps
  State* exon, *aexon;
  for (list<Gene*>::iterator agit = transcripts.begin(); agit != transcripts.end(); agit++){
      for (aexon = (*agit)->exons; aexon != NULL; aexon = aexon->next) {
	 for (exon = gene->exons; exon != NULL; exon = exon->next) {
	     if (!(exon->end < aexon->begin || exon->begin > aexon->end)
		 && frame_compatible(exon, aexon))
		 return true;
	 }  
      }
  }
  return false;
}

void AltGene::shiftCoordinates(int d){
  mincodstart += d;
  maxcodend += d;
  for (list<Gene*>::iterator git=transcripts.begin(); git != transcripts.end(); git++) {
    (*git)->shiftCoordinates(d);
  }
}

/*
 * sort transcripts by 
 * 1. percentage supported by hints
 * 2. a posteriori probabilities
 */
void AltGene::sortTranscripts(int numkeep){
    if (transcripts.size()<2)
	return;
    // bubblesort
    list<Gene *> sortedtr;
    list<Gene *>::iterator it, largestit;
    float maxPercentSupp, maxapostprob;
    double msp, percentSupp;
    if (numkeep<0)
	numkeep = transcripts.size(); // keep all in this case

    while(!transcripts.empty() && (int) sortedtr.size() < numkeep){
	// find best transcript within the rest
	// criteria: 1. % supported 2. mean state posterior probability
	maxPercentSupp = 0.0;
	maxapostprob = -1.0;
	largestit = transcripts.begin(); // to make sure at least one transcript is chosen anyway
	for (it=transcripts.begin(); it!=transcripts.end();it++) {
	    percentSupp = (*it)->getPercentSupported();
	    msp = (*it)->meanStateProb();
	    if (percentSupp > maxPercentSupp || ((percentSupp == maxPercentSupp) && msp > maxapostprob)){
		maxPercentSupp = percentSupp;
		maxapostprob = msp;
		largestit = it;
	    }
	}
	sortedtr.push_back(*largestit);
	transcripts.erase(largestit);
    }
    transcripts = sortedtr;
}

/*
 * deleteSuboptimalTranscripts
 * This is a bit complicated in order to deal with cases like this:
 * >-------                    --->
 * >-------        -----       --->
 *                             >-----              ----               ---->
 *                                                                     >-----        ----->
 */
void AltGene::deleteSuboptimalTranscripts(bool uniqueCDS){
    list<Gene *> suboptimalGenes;
    double percentSupp1, percentSupp2;
    list<Gene *>::iterator it1, it2;
    bool suboptimal, almostidentical;
    int lengthdiff;
    double p1, p2;
    for (it1 = transcripts.begin(); it1 != transcripts.end(); it1++) {
	for (it2 = transcripts.begin(); it2 != transcripts.end(); it2++) {
	    percentSupp1 = (*it1)->getPercentSupported();
	    percentSupp2 = (*it2)->getPercentSupported();
	    suboptimal = (!((*it1)->codingend < (*it2)->codingstart ||  (*it1)->codingstart > (*it2)->codingend) /* TODO: some exon overlaps */ &&
			  percentSupp2 < percentSupp1 * Constant::subopt_transcript_threshold);
	    lengthdiff = (*it1)->geneEnd()-(*it1)->geneBegin() - ((*it2)->geneEnd()-(*it2)->geneBegin()); 
	    p1 = (*it1)->meanStateProb();
	    p2 = (*it2)->meanStateProb();
	    almostidentical = (p1 > p2 || (p1 == p2 && lengthdiff > 0 )) && (*it1)->almostIdenticalTo(*it2);
	    if (uniqueCDS) // if user requires unique CDS
		almostidentical |= (p1 > p2 || (p1 == p2 && lengthdiff > 0 )) && (*it1)->identicalCDS(*it2);
#ifdef DEBUG
	    if (almostidentical) {
		cout << "This transcript is deleted " << endl;
		(*it2)->printGFF();
		cout << "because it is almost identical to this transcript:" << endl;
		(*it1)->printGFF();
	    }
#endif
	    if (it1 != it2 && (suboptimal || almostidentical))
		suboptimalGenes.push_back(*it2);
	}
    }
    suboptimalGenes.sort();
    suboptimalGenes.unique();
    for (it2 = suboptimalGenes.begin(); it2 != suboptimalGenes.end(); it2++) {
#ifdef DEBUG
	cout << "suboptimal or almost identical transcript " << (*it2)->getPercentSupported() << endl;
	(*it2)->printGFF();
#endif
	transcripts.remove(*it2);
    }
}

int AltGene::minTransBegin(){
  int min = INT_MAX;
  for (list<Gene *>::iterator it=transcripts.begin(); it!=transcripts.end();it++)
    if (min> (*it)->geneBegin())
      min = (*it)->geneBegin();
  return min;
}

int AltGene::maxTransEnd(){
  int max = 0;
  for (list<Gene *>::iterator it=transcripts.begin(); it!=transcripts.end();it++)
    if (max < (*it)->geneEnd())
      max = (*it)->geneEnd();
  return max;
}

/* --- Annotation methods ------------------------------------------ */

void Annotation::appendGene(Gene* gene){
    if (!gene) 
	return;
    
    if (lastGene && lastGene->next == NULL)
	lastGene->next = gene;
    else if (genes) {
	for (lastGene = genes; lastGene->next; lastGene = lastGene->next)
	    ;
	lastGene->next = gene;
    } else {
	genes = gene;
	lastGene = gene;
    }
    gene->next = NULL;
}

void Annotation::appendForwardGene(Gene* gene){
    if (!gene) 
	return;
    
    if (lastForwardGene && lastForwardGene->next == NULL)
	lastForwardGene->next = gene;
    else if (forwardGenes) {
	for (lastForwardGene = forwardGenes; lastForwardGene->next; lastForwardGene = lastForwardGene->next)
	    ;
	lastForwardGene->next = gene;
    } else {
	forwardGenes = gene;
	lastForwardGene = gene;
    }
    gene->next = NULL;
}

void Annotation::appendBackwardGene(Gene* gene){ 
    if (!gene) 
	return;
   if (lastBackwardGene && lastBackwardGene->next == NULL)
	lastBackwardGene->next = gene;
    else if (backwardGenes) {
	for (lastBackwardGene = backwardGenes; lastBackwardGene->next; lastBackwardGene = lastBackwardGene->next)
	    ;
	lastBackwardGene->next = gene;
    } else {
	backwardGenes = gene;
	lastBackwardGene = gene;
    }
    gene->next = NULL;
}

const void Annotation::printGFF() const {
    Gene *gene;
    for (gene = genes; gene != NULL; gene = gene->next) {
      gene->printGFF();
    }
}


/* --- AnnoSequence methods ---------------------------------------- */

void AnnoSequence::printGFF(){
    cout << "# Sequence " << seqname << " length=" << length << endl;
    anno->printGFF();
}

AnnoSequence* AnnoSequence::getReverseComplement(){
    AnnoSequence *resAS;
    resAS = new AnnoSequence(); 
    
    resAS->seqname = newstrcpy(seqname);
    resAS->length = length;
    resAS->bc = bc;
    resAS->next = next;
    resAS->weight = weight;
    // reverse complement the sequence
    resAS->sequence = reverseComplement(sequence);
    // reverse complement the genes
    resAS->anno = new Annotation();
    resAS->anno->genes = anno->genes->cloneGeneSequence();
    reverseGeneSequence(resAS->anno->genes, length-1);
    return resAS;
}

const AnnoSeqGeneIterator& operator++(AnnoSeqGeneIterator& gi){
    if (!gi.hasMoreElements())
	throw ProjectError("Assertion failed in AnnoSeqGeneIterator& operator++");
    gi.gene = gi.gene->next;
    bool empty = false;
    while (!gi.gene && !empty) {
	gi.annoseq = gi.annoseq->next;
	if (!gi.annoseq)
	    empty = true;
	else
	    gi.gene = gi.annoseq->anno->genes;
    }
    if (empty)
	    gi.gene = NULL;
    return gi;
}


/* --- StatePathCollection methods --------------------------------- */

/*
 * StatePathCollection::containsPath
 * Checks whether p is already in the collection. See equality operator of class StatePath
 */

bool StatePathCollection::containsPath(StatePath *p){
    bool in = false;
    StatePath *q;
    list<StatePath*>::iterator it = pathlist.begin();
    while(!in && it != pathlist.end()) {
	q = *it;
	in = (*q == *p);
	++it;
    }
    return in;
}

/*
 * StatePathCollection::positionInCollection
 * returns i if *p is the i-th largest Element in the collection, -1 if *p is not in it
 */

int StatePathCollection::positionInCollection(StatePath *p){
    int i=0;
    StatePath *q;
    bool found = false;
    list<StatePath*>::iterator it = pathlist.begin();
    while (it != pathlist.end() && !found) {
	q = *it;
	found = (*q == *p);
	i++;
	++it;
    }
    if (!found)
	return -1;
    else
	return i;
}


/* --- global functions -------------------------------------------- */

/*
 * make a list of genes (with alternative transcrips) to
 * a pure pointer list of transcripts
 */
Gene* getPtr(list<AltGene> *agl){
  if (agl->empty())
    return NULL;
  Gene *head = NULL, *prev=NULL;
  list<AltGene>::iterator ait;
  list<Gene*>::iterator it;
  
  if(agl->begin()->transcripts.empty())
    throw ProjectError("getPtr: Gene without transcripts.");
  
  for (ait = agl->begin(); ait != agl->end(); ait++){
    for (it = ait->transcripts.begin(); it != ait->transcripts.end(); it++){
      if (head == NULL) {
	head = new Gene(**it);
	prev = head;
      } else {
	prev->next = new Gene(**it);
	prev = prev->next;
      }
    }
  }
  if (prev)
    prev->next = NULL;
  return head;
}

void printGeneList(list<AltGene> *genelist, AnnoSequence *annoseq, bool withCS, bool withAA, bool withEvidence){
    string source;
    if (!genelist->empty() && !genelist->begin()->transcripts.empty()){
	source = (*(genelist->begin()->transcripts.begin()))->source;
    }
    // print gene line
    for (list<AltGene>::iterator geneit = genelist->begin(); geneit != genelist->end(); geneit++){
	cout << "# start gene " << geneit->id << endl;
	cout << geneit->seqname << "\t" << source << "\tgene\t" << (geneit->minTransBegin() + 1)<< "\t" << (geneit->maxTransEnd()+1) << "\t";
	if (geneit->hasProbs)
	    cout << setprecision(3) << geneit->apostprob;
	else 
	    cout << ".";
	cout << "\t"  << ((geneit->strand == plusstrand)? '+' : '-') << "\t.\t";
	if (Gene::gff3)
	    cout << "ID=";
	cout << geneit->id << endl;
	for (list<Gene*>::iterator trit = geneit->transcripts.begin(); trit != geneit->transcripts.end(); trit++) {
	    // print transcript lines
	    cout << geneit->seqname << "\t" << source << "\ttranscript\t" << ((*trit)->geneBegin()+1) << "\t" << ((*trit)->geneEnd()+1) << "\t";
	    if ((*trit)->hasProbs)
		cout << setprecision(3) << (*trit)->apostprob;
	    else 
		cout << ".";
	    cout << "\t" << (((*trit)->strand == plusstrand)? '+' : '-') << "\t.\t";
	    if (Gene::gff3)
		cout << "ID=" << geneit->id << "." << (*trit)->id << ";" << "Parent=" << geneit->id << endl;
	    else
		cout << geneit->id << "." << (*trit)->id << endl;

	    (*trit)->printGFF();	   
	    if (annoseq->sequence){
	      if (withCS){
		    (*trit)->printCodingSeq(annoseq);
	      }
	      if (withAA) {
		(*trit)->printProteinSeq(annoseq);
#ifndef DEBUG
		    if (Gene::print_blocks)
#endif
			(*trit)->printBlockSequences(annoseq);
	      }
	    }
	    if (withEvidence){
	      (*trit)->printEvidence();
	    }
	}
	cout << "# end gene " << geneit->id << endl << "###" << endl;
    }
}

void printGeneList(Gene* seq, AnnoSequence *annoseq, bool withCS, bool withAA){
    // make a list of AltGenes with one transcript each from seq
    list<AltGene> *genelist = new list<AltGene>;
    for (Gene* curTx = seq; curTx != NULL; curTx = curTx->next) {
	AltGene gene;
	gene.id = curTx->geneid;
	gene.seqname = string(curTx->seqname);
	gene.addGene(curTx);
	genelist->push_back(gene);
    }
    // and now call the printGeneList function above
    printGeneList(genelist, annoseq, withCS, withAA, false);
}

void printGeneSequence(Gene* seq, AnnoSequence *annoseq, bool withCS, bool withAA){
    if (!seq)
	cout << "# (none)" << endl;
    while(seq) {
	seq->printGFF();
	if (annoseq->sequence){
	  if (withCS)
	    seq->printCodingSeq(annoseq);
	  if (withAA)
	    seq->printProteinSeq(annoseq);
	}
	seq = seq->next;
    }
}
/*
list<Gene> *reverseGeneList(list<Gene> *geneList, int endpos){
    Gene *genes = getPtr(geneList);
    reverseGeneSequence(genes, endpos);
    list<Gene> *reversedList = new list<Gene>;
    while (genes){
	reversedList->push_back(*genes);
	genes = genes->next;
    }
    return reversedList;
}
*/

/*
 * Sort genes by increasing transcription start
 */

list<Gene*>* sortGenePtrList(list<Gene*> glist){
  list<Gene*> *returnlist = new list<Gene*>;
  list<Gene*>::iterator git1, git2;
  // insertion sort, quick if glist already sorted
  for (git1 = glist.begin(); git1 != glist.end(); git1++){
    git2 = returnlist->begin();
    while (git2 != returnlist->end() && (*git1)->geneBegin()>(*git2)->geneBegin()){
      git2++;
    }
    returnlist->insert(git2, *git1);
  }
  return returnlist;
}

list<AltGene> *reverseGeneList(list<AltGene> *altGeneList, int endpos){
  if (altGeneList == NULL)
    return NULL;
  list<AltGene> *reversedList = new list<AltGene>;
  for (list<AltGene>::iterator ait = altGeneList->begin(); ait != altGeneList->end(); ait++) {
    AltGene ag;
    for (list<Gene*>::iterator git = ait->transcripts.begin(); git != ait->transcripts.end(); git++){
      (*git)->next = NULL; // this changes the gene object
      reverseGeneSequence(*git, endpos);
      ag.addGene(*git);
    }
    reversedList->push_back(ag);
  }
  return reversedList;
}

/*
 * groupTranscriptsToGenes
 * summarize overlapping transcripts to genes
 */
list<AltGene>* groupTranscriptsToGenes(list<Gene> *transcripts){
    list<AltGene> *agl = new list<AltGene>;
    list<AltGene>::iterator agit;
    list<Gene>::iterator geneit;
    Gene* g;
    transcripts->sort();

    // loop over the transcripts
    for (geneit = transcripts->begin(); geneit != transcripts->end(); geneit++) {
	g = &(*geneit);
	/*
	 * Merge transcript g into the list of genes. If necessary, several previous genes have to be
	 * merged to one gene.
	 */
	agit = agl->begin();
	list<AltGene>::iterator firstOlp;
	bool overlapFound = false;
	while (agit != agl->end()) {
	    if (agit->overlaps(g)){
		if (overlapFound == false) {
		    // found first overlapping gene, just add this transcript to the existing gene
		    agit->addGene(g);
		    firstOlp = agit;
		    overlapFound = true;
		    agit++;
		} else {
		    // further overlap found for transcript
		    // transcript g overlaps ANOTHER gene, merge the genes
		    for (list<Gene*>::iterator git = agit->transcripts.begin(); git != agit->transcripts.end(); git++){
			firstOlp->addGene(*git);
		    }
		    agit = agl->erase(agit);
		}
	    } else {
		agit++;
	    }
	}
	
	if (!overlapFound){
	    // transcript does not overlap with any previously found gene, create a new gene.
	    AltGene ag;
	    ag.addGene(g);
	    ag.hasProbs = true;
	    agl->push_back(ag);
	}
    }
    return agl;
}

/*
 * reverseGeneSequence(Gene* &seq, int endpos)
 * Reverses the order of the genes and changes the positions of the exons and introns
 * relative to the sequence length 'endpos', changes the strand member variables
 * 
 */
void reverseGeneSequence(Gene* &seq, int endpos){
    Gene *head = seq, *temp;
    State *headst, *tempst;
    seq = NULL;
    int temppos;
    while (head) {
	temp = head->next;
	head->bc.reverse();	

	// reverse exons
	headst = head->exons;
	head->exons = NULL;
	while (headst) {
	    tempst = headst->next;
	    temppos = headst->begin;
	    headst->begin = endpos - headst->end;
	    headst->end = endpos - temppos;
	    switch (headst->type) {
		case initial0 : case initial1: case initial2:
		    headst->type = rinitial; break;
		case terminal : 
		    headst->type = rterminalExon(mod3(2-headst->length())); break;
		case internal0 : case internal1 : case internal2:
		    headst->type = rinternalExon(mod3(2+headst->frame() - headst->length())); break;
		case singleG : 
		    headst->type = rsingleG; break;
	        case rterminal0 : case rterminal1 : case rterminal2:
	            headst->type = terminal; break;
	        case rinitial :
		    headst->type = initialExon(mod3(headst->length())); break;
	        case rinternal0 : case rinternal1 : case rinternal2:
		    headst->type = internalExon(mod3(1+headst->frame() + headst->length())); break;
                case rsingleG : 
		    headst->type = singleG; break;
		default:;
	    }
	    // headst->bc.reverse();
	    headst->next = head->exons;
	    head->exons = headst;
	    headst = tempst;

	}

	//  reverse introns
	headst = head->introns;
	head->introns = NULL;
	while (headst) {
	    tempst = headst->next;
	    temppos = headst->begin;
	    headst->begin = endpos - headst->end;
	    headst->end = endpos - temppos;
	    // headst->bc.reverse();
	    headst->next = head->introns;
	    head->introns = headst;
	    headst = tempst;
	}
	//  reverse 5' utr
	headst = head->utr5exons;
	head->utr5exons = NULL;
	while (headst) {
	    tempst = headst->next;
	    temppos = headst->begin;
	    headst->begin = endpos - headst->end;
	    headst->end = endpos - temppos;
	    // headst->bc.reverse();
	    headst->next = head->utr5exons;
	    head->utr5exons = headst;
	    headst = tempst;
	}

	//  reverse 3' utr
	headst = head->utr3exons;
	head->utr3exons = NULL;
	while (headst) {
	    tempst = headst->next;
	    temppos = headst->begin;
	    headst->begin = endpos - headst->end;
	    headst->end = endpos - temppos;
	    // headst->bc.reverse();
	    headst->next = head->utr3exons;
	    head->utr3exons = headst;
	    headst = tempst;
	}
	// reverse boundaries: transstart, transend, codingstart, codingend
	int triangle = head->codingend;
	head->codingend = (head->codingstart >=0)? endpos-head->codingstart: -1;
	head->codingstart = (triangle>=0)? endpos-triangle: -1;
	triangle = head->transend;
	head->transend = (head->transstart >=0)? endpos-head->transstart: -1;
	head->transstart = (triangle>=0)? endpos-triangle: -1;

	// switch strand member variable
	if (head->strand == plusstrand) 
	    head->strand = minusstrand;
	else if (head->strand == minusstrand) 
	    head->strand = plusstrand;
	
	head->next = seq;
	seq = head;
	head = temp;
    }
}

/*
 * Gene::postProcessGenes
 * postprocess predicted genes
 * currently only truncate hard-masked regions at ends of UTR
 */
void postProcessGenes(list<AltGene> *genes, AnnoSequence *annoseq){
    if (!genes || !annoseq)
	return;
    try {
	if (Properties::getBoolProperty("truncateMaskedUTRs")){
	    for (list<AltGene>::iterator agit = genes->begin(); agit != genes->end(); ++agit){
		for (list<Gene*>::iterator git = agit->transcripts.begin(); git != agit->transcripts.end(); ++git){
		    (*git)->truncateMaskedUTR(annoseq);
		}
	    }
	}
    } catch(...){}
}

void Gene::truncateMaskedUTR(AnnoSequence *annoseq){
    if (utr5exons && strand == plusstrand && complete5utr){
	while (utr5exons->begin <= utr5exons->end && !isNuc(annoseq->sequence + utr5exons->begin)){
	    utr5exons->begin++;
	    transstart = utr5exons->begin;
	    complete5utr = false;
	}
    }
    if (utr3exons && strand == minusstrand && complete3utr){
	while (utr3exons->begin <= utr3exons->end && !isNuc(annoseq->sequence + utr3exons->begin)){
	    utr3exons->begin++;
	    transstart = utr3exons->begin;
	    complete3utr = false;
	}
    }
    if (utr5exons && strand == minusstrand && complete5utr){
	while (utr5exons->end >= utr5exons->begin && !isNuc(annoseq->sequence + utr5exons->end)){
	    utr5exons->end--;
	    transend = utr5exons->end;
	    complete5utr = false;
	}
    }
    if (utr3exons && strand == plusstrand && complete3utr){
	while (utr3exons->end >= utr3exons->begin && !isNuc(annoseq->sequence + utr3exons->end)){
	    utr3exons->end--;
	    transend = utr3exons->end;
	    complete3utr = false;
	}
    } 
}

void examineBaseCountOfGeneSeq(AnnoSequence  *as){
    vector<double> gccontents;
    while (as) {
	as->bc.normalize();
	gccontents.push_back(as->bc.rc + as->bc.rg);
	/*cout << setprecision(3) << (as->bc.rc + as->bc.rg) << " "
	     << "a/t=" << setprecision(3) << as->bc.ra / as->bc.rt << "\t" 
	     << "c/g=" << setprecision(3) << as->bc.rc / as->bc.rg << "\t" 
	     << "at/gc=" << setprecision(3)
	     << (as->bc.ra+ as->bc.rt) / (as->bc.rc + as->bc.rg) << "\t" 
	     << as->bc << " " << as->length << endl;*/
	as = as->next;
    }
    std::sort(gccontents.begin(), gccontents.end());
    cout << "Quantiles of the GC contents in the training set:" << endl;
    for (int i=0; i<21; i++){
	cout << i*100/20 << "%\t" << setprecision(3) << gccontents[(int) (i*(gccontents.size()-1)/20)];
	if (i%2 == 0)
	    cout << endl;
	else
	    cout << "\t";
    }
    if (Constant::decomp_num_steps > 1)
	cout << "I recommend to set /Constant/gc_range_min to " << setprecision(2) << gccontents[(int) (0.05*(gccontents.size()-1))]
	     << " and /Constant/gc_range_max to " << setprecision(2) << gccontents[(int) (0.95*(gccontents.size()-1))] << endl;
}


/* --- CodonUsage: not used right now ------------------------------ 

double CodonUsage::meanLogProb(char *seq, int len, int frame){
    double logProb=0;
    int k;
    int count=0;
    Seq2Int s2i(3);
    for (k = mod3(-frame); k<len; k+=3){
	int codon = s2i(seq+k);
	if (codonusage[codon]!=0.0){ // avoid stop codons and other problems
	    logProb += logcodonusage[s2i(seq+k)];
	    count++;
	}
    }
    return logProb/count;
}

void CodonUsage::init() {
    for (int c=0; c<64; c++) {
	codoncount[c] = 0;
	codonusage[c] = 0.0;
    }
}

void CodonUsage::addCUofSeq(char *seq, int len, int frame){
    Seq2Int s2i(3);
    int k;
    for (k = mod3(-frame); k<len; k+=3) {
	try {
	    codoncount[s2i(seq+k)]++;
	} catch (InvalidNucleotideError e) {
	    // cerr << "CodonUsage::addCUofSeq:" << e.getMessage() << endl;
	}
    }
}
 
void  CodonUsage::computeUsage(){
    int aa, i;
    int sum, numberofcodingcodons = 0;   
    GeneticCode::init();
    // compute relative codon frequencies of synonym codons
    for (aa=0; aa<20; aa++) {
	sum = 0;
	for (i=0; i<GeneticCode::codonsOfAA[aa]; i++) 
	    sum += codoncount[GeneticCode::syncodons[aa][i]];
	numberofcodingcodons += sum;
	for (i=0; i<GeneticCode::codonsOfAA[aa]; i++) {
	    int codon = GeneticCode::syncodons[aa][i];
	    codonusage[codon] = (double) codoncount[codon]/sum;
	    if (codonusage[codon]>0)
		logcodonusage[codon] = log(codonusage[codon]);
	    else {
		cerr << "Warning: codon " << codon << " has relative frequency 0." << endl;
		logcodonusage[codon] = - numeric_limits<double>::infinity();
	    }
	}
    }
    cout << "number of coding codons: " << numberofcodingcodons << endl;

    // compute amino acid frequencies
    for (aa=0; aa<20; aa++) {
	sum = 0;
	for (i=0; i<GeneticCode::codonsOfAA[aa]; i++) 
	    sum += codoncount[GeneticCode::syncodons[aa][i]];
	aaFrequencies[aa] = (double) sum/numberofcodingcodons;
    }
}

void CodonUsage::print(ostream &out){
    for (int aa=0; aa<20; aa++) {
	out << GeneticCode::aa_symbols[aa] << " freq= " << setprecision(4) << aaFrequencies[aa] << "\tcodon freqs= ";
	for (int i=0; i<GeneticCode::codonsOfAA[aa]; i++) 
	    out << "\t" << Seq2Int(3).inv(GeneticCode::syncodons[aa][i]) << ": " << setprecision(3) 
		<< codonusage[GeneticCode::syncodons[aa][i]] << " ";
	out << endl;
    }
} 
--- CodonUsage ---*/ 

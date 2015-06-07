/**********************************************************************
 * file:    extrinsicinfo.cc
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  extrinsic information/hints, user constraints
 * authors: Mario Stanke, mario@gobics.de
 *
 * date    |   author      |  changes
 * --------|---------------|------------------------------------------
 * 26.11.02| Mario Stanke  | creation of the file
 * 27.08.03| Mario Stanke  | grades of hints
 * 27.05.04| Mario Stanke  | renewed printAccuracyForSequenceSet for multi-gene sequences
 * 27.07.04| Mario Stanke  | exonpartmalus does not assume P-(pitchfork)=1 anymore
 * 24.01.05| Mario Stanke  | sort exons, so hints can be trained on nested genes
 * 11.02.05| Mario Stanke  | rename intron to intronpart
 * 23.03.05| Mario Stanke  | introduce intron hints + getFeatureListAt
 * 20.06.06| Mario Stanke  | introduce irpart, tss, tts hints, start and stop hints can be a range
 * 23.06.06| Mario Stanke  | introduce CDS, CDSpart, UTR, UTRpart hints
 * 26.07.06| Mario Stanke  | introcude nonexonpart hints
 * 28.08.06| Mario Stanke  | HintGroup
 * 10.09.06| Mario Stanke  | {Feature::,HintGroup::}compatibleWith
 * 18.10.06| Mario Stanke  | making an index to the feature lists for quick access
 * 09.02.07| Mario Stanke  | some syntax check in extrinsic.cfg file
 * 10.02.07| Mario Stanke  | enabled individual_liability for different sources
 * 07.05.07| Mario Stanke  | fixed bug in getPosFirstEndAtOrAfter. Occurred when not indexed.
 * 12.05.07| Mario Stanke  | fixed bug in computeHintedDSS that could lead to a crash when an intron hint ended exactly one base after the right sequence end
 * 22.04.08| Mario Stanke  | rewrote getExonListOvlpingRange: more efficient and now also check frame
 * 02.07.08| Mario Stanke  | 1group1gene: fill hint gaps within one group with genic hints
 **********************************************************************/

#include "extrinsicinfo.hh"

// project includes
#include "gene.hh"
#include "projectio.hh"    // for comment, goto_line_after
#include "intronmodel.hh"  // for geometric transition probability
#include "igenicmodel.hh" // for geometric transition probability


// standard C/C++ includes
#include <iomanip>  // for setw, setprecision
#include <climits>
#include <iterator> // for ostream_iterator

Bitmask     SequenceFeatureCollection::validDSS = A_SET_FLAG(Seq2Int(2)("gt")) | A_SET_FLAG(Seq2Int(2)("gc"));
Bitmask     SequenceFeatureCollection::validASS = A_SET_FLAG(Seq2Int(2)("ag"));
set<string> SequenceFeatureCollection::validHintedSites;


void SequenceFeatureCollection::initHintedSplicesites(string ssList) {
    int n = 0;
    bool ok = true;
    do {
	if (n >= ssList.length())
	    return;
	string part = ssList.substr(n,5);
	if (part.length() < 5)
	    ok = false;
	else {
	    try {
		validDSS.set(Seq2Int(2)(part.c_str()));
		validASS.set(Seq2Int(2)(part.c_str() +2));
		validHintedSites.insert(part.substr(0,4));
		ok = part[4] == ',';
 	    } catch (...) {
		ok = false;
	    }
	}
	n+=5;
    } while (ok);
    throw ProjectError("Splice site pattern has bad format (example: 'atac,ggag')");
}

void SequenceFeatureCollection::addFeature(Feature f){
    FeatureType type = f.type;
    featureLists[(int) type].push_back(f);	
    sorted = false;
}

void SequenceFeatureCollection::printFeatures(ostream& out){
    int t, size;
    list<Feature>::iterator f;
    sortFeatureLists();
    Feature sf;
    bool hasFeatures = false;
    int offset = 0;
    try {
	int predictionStart = Properties::getIntProperty( "predictionStart" );
	int predictionEnd = Properties::getIntProperty( "predictionEnd" );
	if (predictionStart == predictionEnd && predictionStart < 0)
	    offset = predictionStart;
    } catch (...) { }

    for (t=0; t < NUM_FEATURE_TYPES; t++) {
	size = featureLists[t].size();
	if (size>0) {
	    hasFeatures = true;
	    for (f=featureLists[t].begin(); f!=featureLists[t].end(); f++) {
		sf = *f;
		if (offset != 0)
		    sf.shiftCoordinates(offset, offset, false);
		out << sf << endl;
	    }
	}
    } 
    if (!hasFeatures && !Constant::MultSpeciesMode)
	out << "# (none)" << endl;
}

void SequenceFeatureCollection::sortFeatureLists(){
    if (sorted)
	return;
    for (int type = 0; type < NUM_FEATURE_TYPES; type++){
	featureLists[type].sort();
    }
    sorted = true;
}

/*
 * SequenceFeatureCollection::getFeatureAt
 * returns a pointer to the feature of type 'type' which ends at 'endPosition' and is on strand
 * 'strand' (0 = +, 1 = -) or NULL if there is none
 * This assumes that there are no two features of the same type with the same endpoint.
 */
Feature *SequenceFeatureCollection::getFeatureAt(FeatureType type, int endPosition, Strand strand) {
    list<Feature>::iterator a, e;
    a = getPosFirstEndAtOrAfter(type, endPosition);
    e = getPosFirstEndAtOrAfter(type, endPosition+K);
    while (a != e) {
	if (a->active && (a->end == endPosition) && (a->strand==strand || a->strand == STRAND_UNKNOWN || a->strand == bothstrands))
	    return &(*a);
	a++;
    }
    return NULL;
}

/*
 * SequenceFeatureCollection::getFeatureListAt
 * returns a pointer to the list of all features of type 'type' which end at 'endPosition' and are on strand
 * 'strand' (0 = +, 1 = -) or NULL if there is none
 * This should replace getFeatureAt in order to be able to deal with hints to splice variants.
 */
Feature *SequenceFeatureCollection::getFeatureListAt(FeatureType type, int endPosition, Strand strand) {
    Feature *hitlist = NULL;
    list<Feature>::iterator a, e;
    a = getPosFirstEndAtOrAfter(type, endPosition);
    e = getPosFirstEndAtOrAfter(type, endPosition+K);
    while (a != e) {
	if (a->active && (a->end==endPosition) && (a->strand==strand || a->strand == STRAND_UNKNOWN || a->strand == bothstrands)){
	    a->next = hitlist;
	    hitlist = &(*a);
	}
	a++;
    }
    return hitlist;
}

/*
 * This constructor returns a new SequenceFeatureCollection which contains only a part of this 
 * Of all Feature lists only those Features are kept which end in the interval [from, to]
 * And all coordinates are shifted by -from in the result.
 * If the rc flag is true then the reverse complement is taken, i.e. the strand is reversed and
 * the new origin is at to.
 * This is used for the piecewise prediction.
 */
SequenceFeatureCollection::SequenceFeatureCollection(SequenceFeatureCollection& other, int from, int to, bool rc){
    if (!other.sorted) 
	other.sortFeatureLists();
    featureLists = new list<Feature>[NUM_FEATURE_TYPES];
    list<Feature>::const_iterator it;
    for (int type = 0; type < NUM_FEATURE_TYPES; type++){
	for (it = other.featureLists[type].begin(); it!=other.featureLists[type].end(); it++){
	    if (it->end >= from && it->end <= to) {
		Feature ff = *it;
		ff.shiftCoordinates(from,to,rc);
		featureLists[type].push_back(ff);
	    }
	}
    }
    collection = other.collection;
    sorted = true;
    hintedSites = NULL; 
    hasLocalSSmalus = NULL;
    groupsSorted = false;
    seqlen = to-from+1;
    groupList = NULL;
    groupGaps = NULL;
    predictionScheme = NULL;
    K = BLOCKSIZE;
    firstEnd = NULL;
    lastStart = NULL;
    computeIndices();
}

/*
 * SequenceFeatureCollection::shift
 */
void SequenceFeatureCollection::shift(int offset){
    for (int type = 0; type < NUM_FEATURE_TYPES; type++)
	for (list<Feature>::iterator it = featureLists[type].begin(); it!=featureLists[type].end(); it++){
	    it->end += offset;
	    it->start += offset;
	}
    computeIndices();
}

/*
 * SequenceFeatureCollection::computeHintedSites
 */
void SequenceFeatureCollection::computeHintedSites(const char* dna) {
    if (hintedSites)
	delete [] hintedSites;
    hintedSites = new Bitmask[seqlen];
    // for (int i=0; i<seqlen; i++)
    // 	hintedSites[i] = 0;

    // go through all the dss, exon, CDS, UTR hints
    list<Feature>::const_iterator it;

    // dss
    for (it=featureLists[(int) dssF].begin(); 
	 it!=featureLists[(int) dssF].end(); it++)
    {
	if (it->strand == plusstrand || 
	    it->strand == bothstrands || it->strand == STRAND_UNKNOWN) 
	    for (int k = it->start < 0 ? 0 : it->start; k <= it->end && k < seqlen-1; k++)
		if (validDSSPattern(dna + k))
		    hintedSites[k].set(forwDSS);
	if (it->strand == minusstrand || 
	    it->strand == bothstrands || it->strand == STRAND_UNKNOWN)
	    for (int k = it->start < 1 ? 1 : it->start; k <= it->end && k < seqlen; k++)
		if (validRDSSPattern(dna + k - 1))
		    hintedSites[k].set(revDSS);
    }
    // ass
    for (it=featureLists[(int) assF].begin(); 
	 it!=featureLists[(int) assF].end(); it++)
    {
	if (it->strand == plusstrand || 
	    it->strand == bothstrands || it->strand == STRAND_UNKNOWN) 
	    for (int k = it->start < 1 ? 1 : it->start; k <= it->end && k < seqlen; k++)
		if (validASSPattern(dna + k - 1))
		    hintedSites[k].set(forwASS); 
	if (it->strand == minusstrand || 
	    it->strand == bothstrands || it->strand == STRAND_UNKNOWN)
	    for (int k = it->start < 0 ? 0 : it->start; k <= it->end && k < seqlen-1; k++)
		if (validRASSPattern(dna + k))
		    hintedSites[k].set(revASS);
    }
    // intron
    // introns are checked for full splice pattern validity
    for (it=featureLists[(int) intronF].begin(); it!=featureLists[(int) intronF].end(); it++){
	if (it->start >= 0 && it->end < seqlen) {
	    string pattern(dna + it->start,  2);
	    pattern.append(dna + it->end -1, 2);
	    if (it->strand == plusstrand || it->strand == bothstrands || it->strand == STRAND_UNKNOWN) {
		if (validSplicePattern(pattern)) {
		    hintedSites[it->start].set(forwDSS); 
		    hintedSites[it->end].set(forwASS);
		}
	    }
	    if (it->strand == minusstrand || it->strand == bothstrands || it->strand == STRAND_UNKNOWN) {
		if (validRSplicePattern(pattern)) {
		    hintedSites[it->end].set(revDSS); 
		    hintedSites[it->start].set(revASS);
		}
	    }
	}
    }
    // exon/CDS/UTR
    FeatureType exontypes[3] = { exonF, CDSF, UTRF };
    for (int i=0; i<3; i++) {
	const list<Feature>& flist = featureLists[exontypes[i]];
	for (it=flist.begin(); it!=flist.end(); it++) 
	    if (it->start > 1 && it->end < seqlen-2) {
		if (it->strand == plusstrand || it->strand == bothstrands || it->strand == STRAND_UNKNOWN) {
		    if (validDSSPattern(dna + it->end +1))
			hintedSites[it->end+1].set(forwDSS);
		    if (validASSPattern(dna + it->start -2))
			hintedSites[it->start-1].set(forwASS);
		}
		if (it->strand == minusstrand || it->strand == bothstrands || it->strand == STRAND_UNKNOWN) {
		    if (validRDSSPattern(dna + it->start -2))
			hintedSites[it->start-1].set(revDSS);
		    if (validRASSPattern(dna + it->end +1))
			hintedSites[it->end+1].set(revASS);
		}
	    }
    }
   
/*
    cout << "allowed splice sites:" << endl;
    for (int i=0; i<seqlen; i++){
	if (hintedSites[i] >0) {
	    cout << i ;
	    if (hintedSites[i][forwDSS])
		cout << " d+";
	    if (hintedSites[i][revDSS])
		cout << " d-";
	    if (hintedSites[i][forwASS])
		cout << " a+";
	    if (hintedSites[i][revASS])
		cout << " a-";

	    cout << endl;
	}
	}*/
}


/*
 * SequenceFeatureCollection::deleteFeatureAt
 * deletes the feature of type 'type' which ends at 'endPosition' and is on strand
 * 'strand' (0 = +, 1 = -) if there is none
 * I assume now that there are no two features of the same type with the same endpoint.
 */
void SequenceFeatureCollection::deleteFeatureAt(FeatureType type, int endPosition, Strand strand){
    list<Feature>::iterator f;
    bool erased=false;
    for (f=featureLists[(int) type].begin(); !erased && f!=featureLists[(int) type].end(); f++) 
	if ((f->end==endPosition) && (f->strand==strand || f->strand == STRAND_UNKNOWN || f->strand == bothstrands)){
	    featureLists[(int) type].erase(f);
	    erased = true;
	}
}

/*
 * SequenceFeatureCollection::getAllActiveFeatures
 * 
 */
Feature *SequenceFeatureCollection::getAllActiveFeatures(FeatureType type){
    Feature *head = NULL;
    list<Feature>::iterator fit;
    for (fit = featureLists[type].begin(); fit != featureLists[type].end(); fit++)
	if (fit->active){
	    fit->next = head;
	    head = &(*fit);
	}
    return head;
}

/*
 * SequenceFeatureCollection::getFeatureListInRange
 * returns a pointer to a feature of type 'type' which starts at or after 'startPosition' and ends 
 * before or at 'endPosition' and is on strand
 * seqRelFrame is the position of a first (5' -> 3') nucleotide of a codon modulo 3, -1 if frame does not matter
 * Returns a list of these features, thus changig the next pointer of the elements of this list.
 */
Feature *SequenceFeatureCollection::getFeatureListInRange(FeatureType type, int startPosition,
							  int endPosition, Strand strand, int seqRelFrame) {
    Feature *hitlist = NULL;
    list<Feature>::iterator a, e;
    a = getPosFirstEndAtOrAfter(type, startPosition);
    e = getPosStartAfter(type, endPosition);
    while (a != e) {
	if (a->active && (a->start >= startPosition) && (a->end <= endPosition) && (a->strand == strand || strand==bothstrands || a->strand == STRAND_UNKNOWN)){
	    if ((seqRelFrame == -1 ) || (a->frame == -1) ||
		((a->strand == plusstrand  || a->strand == bothstrands || a->strand == STRAND_UNKNOWN) && (mod3(a->start + a->frame - seqRelFrame) == 0)) ||
		((a->strand == minusstrand || a->strand == bothstrands || a->strand == STRAND_UNKNOWN) && (mod3(a->end - a->frame - seqRelFrame) == 0))) {
		a->next = hitlist;
		hitlist = &(*a);
	    }
	}
	a++;
    } 
    return hitlist;
}

/*
 * SequenceFeatureCollection::getFeatureListBeginningInRange
 * returns a pointer to a feature of type 'type' which starts at or after 'startPosition' and starts 
 * before or at 'endPosition' and is on strand
 * seqRelFrame is the position of a first (5' -> 3') nucleotide of a codon modulo 3, -1 if frame does not matter
 * Returns a list of these features, thus changig the next pointer of the elements of this list.
 * not (yet) needed
 */
Feature *SequenceFeatureCollection::getFeatureListBeginningInRange(FeatureType type, int startPosition, 
								   int endPosition, Strand strand, int seqRelFrame) {
    Feature *hitlist = NULL;
    list<Feature>::iterator a, e;
    a = getPosFirstEndAtOrAfter(type, startPosition);
    e = getPosStartAfter(type, endPosition);
    while (a != e) {
	if (a->active && (a->start >= startPosition) && (a->start <= endPosition) && (a->strand == strand || strand==bothstrands)){
	    if ((seqRelFrame == -1 ) || (a->frame == -1) ||
		((a->strand == plusstrand  || a->strand == bothstrands || a->strand == STRAND_UNKNOWN) && (mod3(a->start + a->frame - seqRelFrame) == 0)) ||
		((a->strand == minusstrand || a->strand == bothstrands || a->strand == STRAND_UNKNOWN) && (mod3(a->end - a->frame - seqRelFrame) == 0))) {
		a->next = hitlist;
		hitlist = &(*a);
	    }
	}
	a++;
    }
    return hitlist;
}

/*
 * SequenceFeatureCollection::getFeatureListOvlpingRange
 * returns a pointer to a feature of type 'type' whose inverval overlaps the interval
 * going from startPosition to endPosition
 * Returns a list of these features, thus changig the next pointer of the elements of this list.
 */
Feature *SequenceFeatureCollection::getFeatureListOvlpingRange(FeatureType type, int startPosition, 
							       int endPosition, Strand strand) {
    Feature *hitlist = NULL;
    list<Feature>::iterator a, e;
    a = getPosFirstEndAtOrAfter(type, startPosition);
    e = getPosStartAfter(type, endPosition);

    while (a != e) {
	if (a->active && (((a->start >= startPosition) && (a->start <= endPosition)) ||
			   ((a->start <= startPosition) && (a->end >= startPosition))) &&
	    (a->strand == strand || strand==bothstrands || a->strand == STRAND_UNKNOWN || a->strand == bothstrands)){
	    a->next = hitlist;
	    hitlist = &(*a);
	}
	a++;
    }
    return hitlist;
}

/*
 * SequenceFeatureCollection::getFeatureListOvlpingRange
 * returns a pointer to a feature of type 'type' whose inverval overlaps the interval
 * going from startPosition to endPosition
 * Returns a list of these features, thus changig the next pointer of the elements of this list.
 * This function can be used to get the lists for several feature types at the same time
 */
Feature *SequenceFeatureCollection::getFeatureListOvlpingRange(Bitmask featuretypes, int startPosition, 
							       int endPosition, Strand strand) {
   if (!collection->hasHintsFile)
      return NULL;
    Feature *hitlist = NULL;
    list<Feature>::iterator a, e;
    for (int type=0; type<NUM_FEATURE_TYPES; type++) {
	if (featuretypes[type]) {
	    a = getPosFirstEndAtOrAfter(type, startPosition);
	    e = getPosStartAfter(type, endPosition);
	    while (a != e) {
		if (a->active && (((a->start >= startPosition) && (a->start <= endPosition)) ||
		     ((a->start <= startPosition) && (a->end >= startPosition))) &&
		    (a->strand == strand || strand==bothstrands || a->strand == STRAND_UNKNOWN || a->strand == bothstrands)){
		    a->next = hitlist;
		    hitlist = &(*a);
		}
		a++;
	    }
	}
    }
    return hitlist;
}

/*
 * SequenceFeatureCollection::getFeatureListContaining
 * returns a pointer to a feature of type 'type' which starts at or before 'startPosition' and ends 
 * before or after 'endPosition' and is on strand
 * 'strand' (0 = +, 1 = -) or NULL if there is none
 */
Feature *SequenceFeatureCollection::getFeatureListContaining(Bitmask featuretypes, int position, Strand strand) {
    if (!collection->hasHintsFile)
	return NULL;
    Feature *hitlist = NULL;
    list<Feature>::iterator a, e;
    for (int type=0; type<NUM_FEATURE_TYPES; type++) {
	if (featuretypes[type]) {
	    a = getPosFirstEndAtOrAfter(type, position);
	    e = getPosStartAfter(type, position);
	    while (a != e) {
		if (a->active && (a->start <= position) && (a->end >= position) && (a->strand == strand || strand==bothstrands || a->strand == STRAND_UNKNOWN)) {
		    a->next = hitlist;
		    hitlist = &(*a);
		}
		a++;
	    }
	}
    }
    return hitlist;
}


/*
 * checkGroupConsistency
 * If one hint is inconsistent with the sequence delete it and the whole group of hints.
 * If a group is not oriented find out whether it must be on the plus or minus strand.
 */

void SequenceFeatureCollection::checkGroupConsistency(AnnoSequence *seq){
    ostringstream messages;
    int numForcedStrand = 0;
    int numDeleted = 0;
    list<HintGroup>::iterator git;
    list<Feature*>::iterator fit;
    list<Feature*> *hints;
    Feature *hint;
    int len = strlen(seq->sequence);
    bool groupOK, plusPossible, minusPossible, setOrientation;

    if (!groupList)
	return;
    for (git = groupList->begin(); git != groupList->end();) {
	groupOK = plusPossible = minusPossible = true;
	setOrientation = false;
	hints = git->getHints();
	if (hints) {
	    for (fit = hints->begin(); fit != hints->end(); fit++) {
		hint = *fit;
		if (hint->strand == plusstrand)
		    minusPossible = false;
		if (hint->strand == minusstrand)
		    plusPossible = false;
		if (hint->start >= 0 &&  hint->end < len) {
		    if (hint->type == intronF) {
		        if (hint->end < hint->start + Constant::min_intron_len-1) {// default minimal intron length is 39
			    messages << "# Error: intron hint is too short." << endl << "# " << *hint << endl;
			    groupOK = false;
			}
			string pattern(seq->sequence + hint->start,  2);
			pattern.append(seq->sequence + hint->end -1, 2);

			plusPossible = plusPossible && validSplicePattern(pattern);
			minusPossible = minusPossible && validRSplicePattern(pattern);
		    } else if (hint->type == assF) {
			bool found = false;
			for (int pos = hint->start-1; pos < hint->end && !found; pos++) 
			    found = validASSPattern(seq->sequence + pos);
			plusPossible = plusPossible && found;
			found = false;
			for (int pos = hint->start; pos <= hint->end && !found; pos++)
			    found = validRASSPattern(seq->sequence + pos);
			minusPossible = minusPossible && found;
		    }  else if (hint->type == dssF) {
			bool found = false;
			for (int pos = hint->start; pos <= hint->end && !found; pos++) 
			    found = validDSSPattern(seq->sequence + pos);
			plusPossible = plusPossible && found;
			found = false;
			for (int pos = hint->start-1; pos < hint->end && !found; pos++) 
			    found = validRDSSPattern(seq->sequence + pos);
			minusPossible = minusPossible && found;
		    } else if (hint->type == CDSF) {
			if (hint->end < hint->start) {
			    messages << "# Error: CDS hint has negative length." << endl << "# " << *hint << endl;
			    groupOK = false;
			}
			const char* cStart = seq->sequence + hint->start;
			const char* cEnd = seq->sequence + hint->end +1;
			plusPossible = plusPossible
			    && (validASSPattern(cStart-2) || GeneticCode::isStartcodon(cStart))
			    && (validDSSPattern(cEnd) || GeneticCode::isStopcodon(cEnd-3));
			minusPossible = minusPossible
			    && (validRDSSPattern(cStart-2) 
				|| GeneticCode::isRCStopcodon(cStart))
			    && (validRASSPattern(cEnd) || GeneticCode::isStartcodon(cEnd-3, true));
		    }
		}
	    }
	    /*
	     * Now make changes to strand or discard group.
	     */
	    if (!plusPossible && !minusPossible) { // no strand possible
		if (!collection->getIndividualLiability((*(hints->begin()))->esource))
		    groupOK = false; // let the whole group suffer if one hints is unsatisfyable
	    } else if (!plusPossible || !minusPossible) { // only one stand possible, set that strand for all hints
		for (fit = hints->begin(); fit != hints->end(); fit++) {
		    if (!minusPossible && (*fit)->strand != plusstrand){
			(*fit)->strand = plusstrand;
			setOrientation = true;
		    }
		    if (!plusPossible && (*fit)->strand != minusstrand){
			(*fit)->strand = minusstrand;
			setOrientation = true;
		    }
		}
		if (setOrientation) {
		    numForcedStrand++;
		    if (Constant::augustus_verbosity>3) {
			messages << "# Set orientation to " << (plusPossible? "forward" : "reverse") << " for hint group ";
			git->print(messages, false);
		    }
		}
	    } else {
		// both strands possible. do nothing.
	    }
	    if (!groupOK){
		git->setDiscardFlag(true);
		numDeleted++; 
		if (Constant::augustus_verbosity>2) {
		    messages << "# Delete group ";
		    git->print(messages, false);
		}
		git=groupList->erase(git);
	    } else {
		git++;
	    }
	} else {
	    git++;
	}
    }
    if (numForcedStrand > 0)
	messages << "# Forced unstranded hint group to the only possible strand for " << numForcedStrand << " groups." << endl;
    if (numDeleted)
	messages << "# Deleted " << numDeleted << " groups because some hint was not satisfiable." << endl;
    
    string msgstring = messages.str();
    if (msgstring.length()>0 && Constant::augustus_verbosity>0){
      	cout << msgstring;
    }
    emptyTrash();
}

#define SUGGEST_SPLICEPATTERNS_OPT \
    "# Use option --allow_hinted_splicesites to enable more splice patterns in hints." 

void SequenceFeatureCollection::warnInconsistentHints(AnnoSequence *seq){
    ostringstream messages;
    list<Feature> liste;
    /*
     * check the coordinate ranges first
     */
    int len = strlen(seq->sequence);
    bool haveOutOfBounds=false;
    list<Feature>::iterator it;
    list<Feature> badFeatures;
    for (int type = 0; type < NUM_FEATURE_TYPES; type++){
	for (it=featureLists[type].begin(); it!=featureLists[type].end(); it++){
	    if (it->end < 0 || it->start >= len) {
		if (!haveOutOfBounds) {
		    messages << "# The following hints have positions out of the boundaries of the sequence and are deleted." << endl;
		    haveOutOfBounds = true;
		}
		messages << "# " << *it << endl;
	       	if ((int)it->type != type)
		    throw ProjectError("Internal error in SequenceFeatureCollection::warnInconsistentHints.");
		badFeatures.push_back(*it);
	    }
	}
    }
    // now delete al features with coordinate out of bounds
    if (haveOutOfBounds) {
	for (it=badFeatures.begin(); it!=badFeatures.end(); it++){
	    featureLists[(int) it->type].remove(*it);
	}
    }

    badFeatures.clear();

    /*
     * check the start hints
     */
    liste = featureLists[(int) startF];
    for (it = liste.begin(); it != liste.end(); it++){
	if (it->end - it->start < 2) {
	    messages << "# Error: start hint has length < 3. Hint range should contain start codon." << endl << "# " << *it << endl;
	    badFeatures.push_back(*it);
	}
	bool hasATG = false;
	for (int i = it->start; i <= it->end-2; i++)
	    hasATG |= ((it->strand == plusstrand || it->strand == STRAND_UNKNOWN || it->strand == bothstrands) && GeneticCode::isStartcodon(seq->sequence+i)) 
		|| ((it->strand == minusstrand || it->strand == STRAND_UNKNOWN || it->strand == bothstrands) && GeneticCode::isStartcodon(seq->sequence+i, false));
	if (!hasATG){ 
	    messages << "# Error: start hint range contains no start codon" << endl << "# " << *it << endl;
	    badFeatures.push_back(*it);
	}
    }
    /*
     * check the stop hints
     */
    liste = featureLists[(int) stopF];
    for (it = liste.begin(); it != liste.end(); it++){
	if (it->end-it->start < 2) {
	    messages << "# Error: stop hint has length < 3. Hint range should contain stop codon." << endl << "# " << *it << endl;
	    badFeatures.push_back(*it);
	}	
	bool hasSTP = false;
	for (int i = it->start; i <= it->end-2; i++)
	    hasSTP |= ((it->strand == plusstrand || it->strand == STRAND_UNKNOWN || it->strand == bothstrands) && GeneticCode::isStopcodon(seq->sequence+i)) 
	      || ((it->strand == minusstrand || it->strand == STRAND_UNKNOWN || it->strand == bothstrands) && GeneticCode::isRCStopcodon(seq->sequence+i));
	if (!hasSTP){
	    badFeatures.push_back(*it);
	    messages << "# Error: stop hint range does not contain valid stopcodon" << endl << "# " << *it << endl;
	}
    }

    for (it=badFeatures.begin(); it!=badFeatures.end(); it++){
	featureLists[(int) it->type].remove(*it);
    }
    badFeatures.clear();
    /*
     * check the ASS hints
     */
    liste = featureLists[(int) assF];
    for (it = liste.begin(); it != liste.end(); it++) {
	bool fuzzy = it->end != it->start;
	if (fuzzy)
	    messages << "# Found fuzzy ass hint: " << endl << "# " << *it << endl;
	bool found = false;
	if (it->strand == plusstrand || it->strand == bothstrands || it->strand == STRAND_UNKNOWN) 
	    for (int k=it->start; k<=it->end && !found; k++)
		found = validASSPattern(seq->sequence + it->start -1);
	if (it->strand == minusstrand || it->strand == bothstrands || it->strand == STRAND_UNKNOWN) 
	    for (int k=it->start; k<=it->end && !found; k++)
		found = validRASSPattern(seq->sequence + it->start);
	if (!found) {
	    if (fuzzy)
		messages << "# Error: no valid pattern found in sequence for ass hint.";
	    else {
		messages << "# Error: invalid pattern '";
		if (it->strand == plusstrand)
		    messages << string(seq->sequence + it->start -1, 2);
		else 
		    putReverseComplement(ostream_iterator<char>(messages),
					 seq->sequence + it->start, 2);
		messages << "' found in sequence for ass hint." << endl << SUGGEST_SPLICEPATTERNS_OPT;
	    }
	    messages << endl << "# " << *it << endl;
	}
    }
    /*
     * check the DSS hints
     */
    liste = featureLists[(int) dssF];
    for (it = liste.begin(); it != liste.end(); it++) {
	bool fuzzy = it->end != it->start;
	if (fuzzy)
	    messages << "# Found fuzzy dss hint: " << endl << "# " << *it << endl;
	bool found = false;
	if (it->strand == plusstrand || it->strand == bothstrands || it->strand == STRAND_UNKNOWN) 
	    for (int k=it->start; k<=it->end && !found; k++)
		found = validDSSPattern(seq->sequence + it->start);
	if (it->strand == minusstrand || it->strand == bothstrands || it->strand == STRAND_UNKNOWN) 
	    for (int k=it->start; k<=it->end && !found; k++)
		found = validRDSSPattern(seq->sequence + it->start-1);
	if (!found) {
	    if (fuzzy)
		messages << "# Error: no valid pattern found in sequence for dss hint.";
	    else {
		messages << "# Error: invalid pattern '";
		if (it->strand == minusstrand)
		    putReverseComplement(ostream_iterator<char>(messages),
					 seq->sequence + it->start-1, 2);
		else 
		    messages << string(seq->sequence + it->start, 2);
		messages << "' found in sequence for dss hint." << endl << SUGGEST_SPLICEPATTERNS_OPT;
	    }
	    messages << endl << "# " << *it << endl;
	}
    }
    /*
     * check the exonpart hints
     */
    liste = featureLists[(int) exonpartF];
    for (it = liste.begin(); it != liste.end(); it++){
	if (it->end < it->start)
	    messages << "# Error: exonpart hint has negative length." << endl << "# " << *it << endl;

    }
    /*
     * check the exon hints
     */
    liste = featureLists[(int) exonF];
    for (it = liste.begin(); it != liste.end(); it++){
	if (it->end < it->start)
	    messages << "# Error: exon hint has negative length." << endl << "# " << *it << endl;
	if (!(it->strand != minusstrand  // could be plus
	      && (validDSSPattern(seq->sequence+it->end+1)
		  || GeneticCode::isStopcodon(seq->sequence+it->end-2))) &&
	    !(it->strand != plusstrand  // could be minus
	      && (validRDSSPattern(seq->sequence+it->start-2)
		  || GeneticCode::isRCStopcodon(seq->sequence+it->start)))) {
	    messages << "# Error: exon hint has invalid splice pattern '";
	    if (it->strand == minusstrand)
		putReverseComplement(ostream_iterator<char>(messages), 
				     seq->sequence + it->start -2, 2);
	    else
		messages << string(seq->sequence + it->end +1);
	    messages << "', and no stop codon, at 3' boundary" << endl
		     << SUGGEST_SPLICEPATTERNS_OPT << endl
		     << "# " << *it << endl;
		
	}
	if (!(it->strand != minusstrand
	      && (validASSPattern(seq->sequence+it->start-2)
		  || GeneticCode::isStartcodon(seq->sequence+it->start))) &&
	    !(it->strand !=plusstrand
	      && !(validRASSPattern(seq->sequence+it->end+1)
		   || GeneticCode::isStartcodon(seq->sequence+it->end-2, true)))) {
	    messages << "# Error: exon hint has invalid splice pattern '";
	    if (it->strand == minusstrand)
		putReverseComplement(ostream_iterator<char>(messages), 
				     seq->sequence + it->end+1, 2);
	    else
		messages << string(seq->sequence + it->start -2);
	    messages << "', and no start codon, at 5' boundary" << endl 
		     << SUGGEST_SPLICEPATTERNS_OPT << endl
		     << "# " << *it << endl;
	}
    }
    /*
     * check the intron hints
     */
    liste = featureLists[(int) intronF];
    for (it = liste.begin(); it != liste.end(); it++){
	if (it->end < it->start + 40)
	    messages << "# Error: intron hint is too short." << endl << "# " << *it << endl;
	string pattern(seq->sequence + it->start,  2);
	pattern.append(seq->sequence + it->end -2, 2);
	if (!(it->strand != minusstrand && validSplicePattern(pattern)) &&
	    !(it->strand != plusstrand && validRSplicePattern(pattern))) 
	    messages << "# Error: sequence has invalid splice pattern '" << pattern 
		     << "' at intron hint." << endl << SUGGEST_SPLICEPATTERNS_OPT << endl
		     << "# " << *it << endl;
    }

    string msgstring = messages.str();
    if (msgstring.length()>0){
	cerr << "# WARNING: Inconsistent hints for sequence " << seq->seqname 
	     << ". These hints are ignored." << endl;
	cerr << msgstring << endl;
    }
}

/*
 * cleanRedundantFeatures
 * deletes splice sites features already covered by 'exon' features
 */
void SequenceFeatureCollection::cleanRedundantFeatures() {
    //list<Feature>::iterator it;
    list<Feature> list;
    sortFeatureLists();
    list = featureLists[(int) exonF];
    // remove splice sites also supported by exonF
    //for (it = list.begin(); it != list.end(); it++){
    //	deleteFeatureAt(assF, it->start - 1, it->strand);
    //	deleteFeatureAt(dssF, it->end + 1, it->strand);
    //}
    // remove multiple Features in the lists
    for (int type = 0; type < NUM_FEATURE_TYPES; type++) 
	deleteEqualElements(featureLists[type]);
}

/*
 * delete all but the first feature of a set with the same startpoint, endpoint, frame, strand
 * requires the list to be sorted
 */
void SequenceFeatureCollection::deleteEqualElements(list<Feature> &flist){
    list<Feature>::iterator it, oldit, 
	lastF = flist.end(); //actually undefined iterator, but NULL isn't working under gcc 4.1.0
    it = flist.begin();
    while (it != flist.end()) {
	if (lastF != flist.end() && *it == *lastF) {
	    oldit = it;
	    it++;
	    // determine which one of the "equal" features should be deleted 
	    // the bonus of the remaining feature should probably be upgraded somehow
	    // but for now simply the feature with the higher bonus remains.
	    double newbonus = lastF->bonus += oldit->bonus;// TEMP HACK for testing Mario
	    if (newbonus > lastF->bonus * oldit->bonus)
	      newbonus = lastF->bonus * oldit->bonus;
	    if (lastF->bonus > oldit->bonus) {
	      lastF->bonus = newbonus;
	      lastF->mult += oldit->mult;
	      flist.erase(oldit);
	    } else {
	      oldit->bonus = newbonus;
	      oldit->mult += lastF->mult;
	      flist.erase(lastF);
	      lastF = oldit;
	    }
	} else {
	    lastF = it;
	    it++;
	}
    }
}

Feature *SequenceFeatureCollection::getExonListInRange(int startPosition, int endPosition,
						       Strand strand, int seqRelFrame) {
    Feature *e, *ep, *temp;
    // get lists for exons and exon parts 
    e  = getFeatureListInRange(exonF, startPosition, endPosition, strand, seqRelFrame);
    ep = getFeatureListInRange(exonpartF, startPosition, endPosition, strand, seqRelFrame);
    // ... and append them
    // TODO: Das hier ist uneffizient
    if (e){
	for (temp = e; temp->next != NULL; temp=temp->next) {}
	temp->next = ep;
    } else
	e=ep;
    return e;
}

/*
 * getExonListOvlpingRange
 * seqRelFrame is the position of a first (5' -> 3') nucleotide of a codon modulo 3, -1 if frame does not matter
 */
Feature *SequenceFeatureCollection::getExonListOvlpingRange(int startPosition, int endPosition,
							    Strand strand, int seqRelFrame) {
    
    Feature *first, *last, *lst, *temp;
    bool reading_frame_ok;
    lst = getFeatureListOvlpingRange(A_SET_FLAG(exonF) | A_SET_FLAG(exonpartF) | A_SET_FLAG(CDSF) | A_SET_FLAG(CDSpartF) | A_SET_FLAG(UTRF) | A_SET_FLAG(UTRpartF), 
				     startPosition, endPosition, strand);

    // now filter list by keeping only the ones with a right reading frame
    first = last = NULL;
    temp = lst;
    while (temp) {
	reading_frame_ok = true;
	if ((temp->type == CDSF || temp->type == CDSpartF || temp->type == exonF || temp->type == exonpartF) &&
	    (temp->frame >= 0 && seqRelFrame >= 0) &&
	    ((temp->strand == plusstrand && mod3(temp->frame + temp->start - seqRelFrame) != 0) ||
	     (temp->strand == minusstrand && mod3(temp->end - temp->frame - seqRelFrame) != 0)))
	    reading_frame_ok = false;
	if (reading_frame_ok) { // skip all hints with wrong reading frame
	    if (!first) {
		first = last = temp;
	    } else {
		last->next = temp;
		last = temp;
	    }
	}
	temp = temp->next;
    }
    if (last)
	last->next = NULL;
    return first;
}

/*
 * makeGroups
 * makes groups of hints accoding to the group name given in the hintsfile
 */ 
void SequenceFeatureCollection::makeGroups(){
    // HashTable<HintGroup*, string, HashString<99999>, 100000 > groupnames;
    map<string, HintGroup*> groupnames;
    Feature *hint;
    list<HintGroup>::reverse_iterator grit;
    list<Feature>::iterator it;
    int grpcount = 0, hintcount=0;
    if (groupList)
	delete groupList;
    groupList = new list<HintGroup>;
    // loop through all hints and add them to the groups
    for (int type = 0; type < NUM_FEATURE_TYPES; type++){
	for (it=featureLists[(int) type].begin(); it!=featureLists[(int) type].end(); it++){
	    hintcount++;
	    hint = &(*it);
	    if (hint->groupname == "" || groupnames.count(hint->groupname)==0 || groupList->size() == 0) {
		// create new group
		grpcount++;
		groupList->push_back(HintGroup());
		groupList->back().addFeature(hint);
		grit = groupList->rbegin();
		if (hint->groupname != ""){
		    groupnames[hint->groupname] = &(*grit);
		}
	    } else { // find the group with the same name
		groupnames[hint->groupname]->addFeature(hint);
	    }
	}
    }
    groupList->sort();
    groupsSorted = true;
}

/*
 * printGroups
 */ 
void SequenceFeatureCollection::printGroups(){
    cout << "---------------  Hint groups  -----------------" << endl;
    for (list<HintGroup>::iterator grit = groupList->begin(); grit != groupList->end(); grit++){
	cout << "name=" << grit->getName() << " members:" << grit->getSize() << " " << grit->getSource() << " range:" << grit->getBegin() << "-" << grit->getEnd();
	if (grit->getEnd()-grit->getBegin()>1000000) {
	    cout << " grosse Gruppe: " << grit->getEnd()-grit->getBegin();
	}
	cout << endl;
    }
}

/*
 * findGenicGaps
 * Find gaps within the groups and create nonirpart hints for them. The nonirpart hints are added to the group.
 * example
 * eeeeiiiiiiieeeeee                 eeeeeeiiiiiiiiiiiiiieeeee   
 *                  ggggggggggggggggg
 *                     genic hint
 * Groups with gaps occur when a 5' and 3' EST from the same clone align with a gap in between. They should be bridged to one gene.
 * Reason for not making the whole thing genicpart: EST alignments without UTR model would lead to an extension of CDS region.
 */ 
void SequenceFeatureCollection::findGenicGaps(){
   if (!groupList)
	return;
   int mingaplen = 1;
   int maxgaplen = 500000; // to prevent coincidental false grouping, e.g. across half a chromosome, to disturb the results
   double nonirpartbonus; // the inverse of this factor applies to every igenic region base in the genic gaps
   try {
     nonirpartbonus = Properties::getdoubleProperty( "bridge_genicpart_bonus" );
   } catch(...) {
     nonirpartbonus = 2.0;
     return;
   }


   bool found;
   for (list<HintGroup>::iterator grit = groupList->begin(); grit != groupList->end(); grit++){ 
//     grit->print(cout, true);
       found = false;
       list<Feature*> *featlist = grit->getHints();
       if (collection->get1group1gene((*(featlist->begin()))->esource)) {
	   // iterate through feature list from right to left (cause it is sorted by end pos.)
	   if (featlist->size() > 1) { // no gaps if there is just one hint in the group
	       list<Feature*>::reverse_iterator hit = featlist->rbegin(), prevhit = hit;
	       prevhit++;
	       while (hit != featlist->rend() && prevhit != featlist->rend()){
		   if ((*prevhit)->end >= (*hit)->start - mingaplen){// no gap between hit and prevhit
		       if ((*prevhit)->start < (*hit)->start)
			   hit = prevhit;
		   } else { // found gap between prevhit and hit
		       if ((*hit)->start - (*prevhit)->end -1 <= maxgaplen){
			   found = true;
			   Feature *genicgap = new Feature((*prevhit)->end+1, (*hit)->start-1, nonirpartF, bothstrands, -1, (*hit)->esource);
			   genicgap->seqname = (*hit)->seqname;
			   genicgap->source = (*hit)->source;
			   genicgap->bonus = nonirpartbonus;
			   genicgap->malus = 1;
			   genicgap->feature = "genicpart";
			   genicgap->groupname = (*hit)->groupname;
			   genicgap->priority = grit->getPriority();
			   grit->addFeature(genicgap);
			   addFeature(*genicgap);
		       }
		       hit = prevhit;
		   }
		   prevhit++;
	       }
	   }
	   if (found) {
	       // grit->print(cout, true);
	   }
       }
   }
}

/*
 * findGroupGaps
 * Find gaps in between different groups and store them as irpart hints.
 */ 
void SequenceFeatureCollection::findGroupGaps(){
    if (!groupList)
	return;
    if (groupGaps)
	delete groupGaps;
    groupGaps = new list<Feature>;
    groupGaps->push_back(Feature(1, seqlen, irpartF, bothstrands, -1, "groups"));
    if (!groupsSorted)
	groupList->sort();
    list<Feature>::iterator irit = groupGaps->begin();
    bool finished = false;
    for (list<HintGroup>::iterator grit = groupList->begin(); grit != groupList->end() && ! finished; grit++){
	if (grit->getGeneBegin() >=0 && grit->getGeneBegin() <= seqlen) { // exlude e.g. irpart hints at this place
	    if (grit->getGeneBegin() > irit->start && grit->getGeneEnd() < irit->end) { 
		// group proper subinterval of irpart, chop out group interval and split irpart into two irparts
		groupGaps->insert(irit, Feature(irit->start, grit->getGeneBegin()-1, irpartF, bothstrands, -1, "groups"));
		irit->start = grit->getGeneEnd() + 1;
	    } else if (irit->start >= grit->getGeneBegin() && grit->getGeneEnd() >= irit->start && grit->getGeneEnd() < irit->end ) {
		// group chops off a proper initial part of irpart
		irit->start = grit->getGeneEnd() + 1;
	    } else if (grit->getGeneBegin() <= irit->start && grit->getGeneEnd()>= irit->end) {
		// group contains irpart, delete irpart completely
		groupGaps->erase(irit);
		finished = true;
	    } else if (grit->getGeneBegin() > irit->start && grit->getGeneBegin() <= irit->end && grit->getGeneEnd() >= irit->end) {
		// group chops off a terminal part of irpart
		irit->end = grit->getGeneBegin() - 1;
		finished = true;
	    } else if (grit->getGeneEnd()<irit->start){
		// do nothing
	    } else {
		cerr << "SequenceFeatureCollection::findGroupGaps() should not happen" << endl;
		cerr << "group: " << grit->getGeneBegin() << ".." << grit->getGeneEnd() << endl;
		cerr << "irpart:" << irit->start << ".." << irit->end << endl;
	    }
	}
    }
    /*cout << "----------- group gaps ----------" << endl;
    for (irit = groupGaps->begin(); irit != groupGaps->end(); irit++) {
	cout << irit->start <<".." << irit->end << endl;
	}*/
}


/*
 *  determineInterGroupRelations
 */
void SequenceFeatureCollection::determineInterGroupRelations(){
    list<HintGroup>::iterator git1, git2;
    Feature *rascal1=NULL, *rascal2=NULL;
    if (!groupsSorted)
	groupList->sort();
    int numIncompat = 0, numWeakerThan = 0, numEqual = 0;
    bool weakerThan, compatibleWith;
    // check for absolutely equal groups, keep only one copy and store the copy number
    //cout << "groupList->size()=" << groupList->size() << endl;
    bool prevEqual;
    for (git1 = groupList->begin(); git1 != groupList->end(); git1++) {
	git2 = git1;
	git2++;
	prevEqual = true;
	/*
	 * The HintGroups are sorted such that if g1 < g2 < g3 then g1=g3 implies g1=g2.
	 * Therefore we can shorten this deletion of equal elements drastically.
	 */
	while (git2 != groupList->end() && git2->getBegin() <= git1->getEnd() && prevEqual){
	    if (*git1 == *git2) {
		//cout << "The following two groups are equal:" << endl;
		//git1->print(cout, true);
		//git2->print(cout, true);
		numEqual++;
		git1->addCopyNumber(git2->getCopyNumber());
		git2->setDiscardFlag(true);
		git2 = groupList->erase(git2);
	    } else {
		prevEqual = false;
	    }
	}
    }
#ifdef DEBUG
    cout << "done looking for equal groups. " << numEqual << " equal" << endl;
    cout << "groupList->size()=" << groupList->size() << endl;
#endif
    /*
     * Update the amount each feature conforms to the other features and delete groups with nonconformant features.
     */
    resetConformance();
    git1 = groupList->begin();
    int numTrashy=0;
    while (git1 != groupList->end()) {
	git1->updateFeatureConformance(*git1); // with itself, in case several hint groups are identical
	git2 = git1;
	git2++;
	while (git2 != groupList->end() && git2->getBegin() <= git1->getEnd()){
	    git1->updateFeatureConformance(*git2);
	    git2->updateFeatureConformance(*git1);
	    git2++;
	}
	if (git1->isTrashy()){
	  //cout << "is trashy: ";
	  //git1->print(cout, true);
	    git1 = groupList->erase(git1);
	    numTrashy++;
	} else 
	    git1++;
    }
#ifdef DEBUG
      cout << "# " << numTrashy << " hint groups are conflicting with too many others. Deleting them." << endl;
#endif
    emptyTrash();

    /*
     * determine incompatible group pairs and 'weaker than' relationships
     */
    for (git1 = groupList->begin(); git1 != groupList->end(); git1++) {
	git2 = git1;
	git2++;
	while (git2 != groupList->end() && git2->getBegin() <= git1->getEnd()){
	    // check compatibility of the two groups
	    compatibleWith = git1->compatibleWith(*git2, rascal1, rascal2, weakerThan);
	    if (!compatibleWith){
		git1->addIncompGroup(&*git2);
		git2->addIncompGroup(&*git1);
		//cout << "The following two groups are incompatible with each other:" << endl;
		//git1->print(cout);
		//git2->print(cout);
		//cout << "Reason of incompatibility: The following two hints are incompatible:" << endl;
		//cout << "# " << *rascal1 << endl;
		//cout << "# " << *rascal2 << endl;
		numIncompat++;
	    }
	    if (weakerThan) {
		numWeakerThan++;
		git1->addStrongerGroup(&*git2);
		//cout << "Group ";
		//git1->print(cout);
		//cout << "is weaker than group ";
		//git2->print(cout);
	    }
	    git2->compatibleWith(*git1, rascal1, rascal2, weakerThan); // check both orders for weakerThan
	    if (weakerThan) {
		numWeakerThan++;
		git2->addStrongerGroup(&*git1);
		//cout << "Group ";
		//git2->print(cout);
		//cout << "is weaker than group ";
		//git1->print(cout);
	    }
	    git2++;
	}
    }
#ifdef DEBUG
    cout << numIncompat   << " pairs of groups are incompatible (" << groupList->size() << " groups total)." << endl;
    cout << numWeakerThan << " ordered pairs of groups have a weakerThan relationship" << endl;
#endif
    // sort all incompGroups lists
    sortIncompGroupsOfGroups();
}

void SequenceFeatureCollection::resetConformance(){
    for (int type=0; type < NUM_FEATURE_TYPES; type++)
	for (list<Feature>::iterator it = featureLists[(int) type].begin(); it != featureLists[(int) type].end(); it++)
	    it->numContradicting = it->numSupporting = 0;
}

/*
 * emptyTrash
 * Delete features that have been marked to be discarded.
 */
void SequenceFeatureCollection::emptyTrash(){
    for (int type=0; type < NUM_FEATURE_TYPES; type++){
	list<Feature> &liste = featureLists[(int) type];
	list<Feature>::iterator it = liste.begin();
	while(it != liste.end())
	    if (it->discard)
		it = liste.erase(it);
	    else
		it++;
    }
}

/*
 * computeIndices
 * Computes indices to the feature lists so that access to the lists later is quicker.
 */
void SequenceFeatureCollection::computeIndices() {
    if (featureLists == NULL)
	return;
    if (firstEnd) {
	for (int type=0; type < NUM_FEATURE_TYPES; type++)
	    delete [] firstEnd[type];
	delete [] firstEnd;
	firstEnd = NULL;
    }
    if (lastStart) {
	for (int type=0; type < NUM_FEATURE_TYPES; type++)
	    delete [] lastStart[type];
	delete [] lastStart;
	lastStart = NULL;
    }
    if (2*K >= seqlen) { // sequence too short, don't bother making an index
	firstEnd = NULL;
	lastStart = NULL;
	return;
    }
    int numBlocks = seqlen/K+1;
    int lastChanged;
    list<Feature>::iterator fit;
    int j;
    firstEnd = new list<Feature>::iterator *[NUM_FEATURE_TYPES];
    lastStart = new list<Feature>::iterator *[NUM_FEATURE_TYPES];
    for (int type=0; type < NUM_FEATURE_TYPES; type++){
	if (featureLists[type].empty()){
	    firstEnd[type] = lastStart[type] = NULL; // don't index this hint type
	} else {
	    firstEnd[type] = new list<Feature>::iterator [numBlocks];
	    lastStart[type] = new list<Feature>::iterator [numBlocks];
	    /*
	     * create the indices based on the ends
	     */
	    lastChanged = -1;
	    for (fit = featureLists[type].begin(); fit != featureLists[type].end(); fit++) {
		if (fit->end/K > lastChanged) {// change indices from lastChanged until fit->end/K
		    for (j=lastChanged+1; j <= fit->end/K && j<numBlocks; j++){
			firstEnd[type][j] = fit;
			lastChanged = j;
		    }
		}
	    }
	    //set the rest of the indices if there is a rest
	    while(lastChanged < numBlocks-1)
		firstEnd[type][++lastChanged] = featureLists[type].end();
	    /*
	     * create the indices based on the starts
	     * This is more complicated since the feature lists are sorted by the end positions of hints.
	     */
	    int leftmostStart = seqlen+1;
	    int lastChanged = numBlocks;
	    list<Feature>::iterator lastUndershootingIt = featureLists[type].end();
	    featureLists[type].reverse(); // reverse sequence so we can iterate from the right to the left without using a reverse_iterator
	    for (fit = featureLists[type].begin(); fit != featureLists[type].end(); fit++) {
		if (fit->start < leftmostStart)
		    leftmostStart = fit->start;
		if (leftmostStart <= (lastChanged-1)*K) {// change indices from lastChanged until fit->start/K
		    for (j=lastChanged-1; j > (leftmostStart-1)/K && j>=0; j--){
			lastStart[type][j] = fit;
			lastUndershootingIt = fit;
			lastChanged = j;
		    }
		}
	    }
	    featureLists[type].reverse();
	    //set the rest of the indices if there is a rest
	    while (lastChanged > 0)
		lastStart[type][--lastChanged] = lastUndershootingIt;
 
	    /*
	    cout << " computed Indices to feature list " << type << endl;
	    for (int i=0; i<numBlocks; i++) {
		cout << " firstEnd[" << type << "][" << i*K << "]= ";
		if (firstEnd[type][i] != featureLists[type].end()) 
		    cout << *(firstEnd[type][i]) << endl;
		else 
		    cout << "Listenende" << endl;
		if (lastStart[type][i] != featureLists[type].end()) 
		    cout << "lastStart[" << type << "][" << i*K << "]= " << *(lastStart[type][i]) << endl;
		else 
		    cout << "lastStart[" << type << "][" << i*K << "]= " << "Listenende" << endl;
		    }*/
	    
	}
    }
}

/*
 * getPosFirstEndAtOrAfter
 */
list<Feature>::iterator SequenceFeatureCollection::getPosFirstEndAtOrAfter(int type, int e) {
    if (firstEnd == NULL || firstEnd[type] == NULL){
	// Don't have index. Just compute it on the fly.
	list<Feature>::iterator fit;
	for (fit = featureLists[type].begin(); fit != featureLists[type].end(); fit++) {
		if (fit->end >= e) 
		    return fit;		
	}
	return featureLists[type].end();
    }
    int i = e/K;
    if (i < 0)
	i = 0;
    if (i >= seqlen/K)
	i = seqlen/K;
    return firstEnd[type][i];
}

/*
 * getPosStartAfter
 */
list<Feature>::iterator SequenceFeatureCollection::getPosStartAfter(int type, int s) {
    if (lastStart == NULL || lastStart[type] == NULL)
	return featureLists[type].end();
    list<Feature>::iterator ret;
    int i = (s+K-1)/K;
    if (i < 0)
	i = 0;
    if (i >= seqlen/K)
	return featureLists[type].end();
    ret = lastStart[type][i];
    if (ret != featureLists[type].end())
	ret++;
    return ret;
}

/*
 * rescaleBoniByConformance
 * Adjust boni of all hints, accoding to the number of other hints that are supporting it and the number 
 * of other hints that are contradicting it.
 * A hint with a lot of contradicting other hints gets a smaller bonus.
 */
void SequenceFeatureCollection::rescaleBoniByConformance(){
    double conformance; // between 0 and 1. Close to 1 if lots of support and no contradiction.
    for (int type=0; type < NUM_FEATURE_TYPES; type++){
	list<Feature> &liste = featureLists[(int) type];
	for (list<Feature>::iterator it = liste.begin(); it != liste.end(); it++){
	    conformance = it->conformance();
	    if (it->bonus > 0.0)
		it->bonus = exp(log(it->bonus) * 2 * conformance);
	}
    }
}

/*
 * createPredictionScheme
 * This function decides which prediction runs (PredictionRuns) to make later.
 * Predictions can be made several time in the same area of the sequence if incompatible hints are found.
 * Each hint group is used at least once in a  PredictionRun such that all incompatible  hint groups are inactivated in this run.
 * The individual prediction intervals are set such that no 'genes' are separated.
 */
void SequenceFeatureCollection::createPredictionScheme(list<AltGene> *genes){
    if (predictionScheme)
	delete predictionScheme;
    predictionScheme = new PredictionScheme(seqlen);
    if (groupList == NULL || !Constant::alternatives_from_evidence)
	return; // no additional run needed
    determineInterGroupRelations();
    for (list<HintGroup>::iterator git = groupList->begin(); git != groupList->end(); git++) {
	bool overruled = false;
	bool weakerThanOther = false;
	if (git->getIncompGroups())
	    for (list<HintGroup*>::iterator incit = git->getIncompGroups()->begin(); incit != git->getIncompGroups()->end(); incit++){
		if ((*incit)->getPriority() >= 0 && git->getPriority() >= 0 && (*incit)->getPriority() > git->getPriority()) 
		    overruled = true; // incompatible with a group of higher priority, don't start a run for that
	    }
	if (git->getStrongerGroups())
	    for (list<HintGroup*>::iterator incit = git->getStrongerGroups()->begin(); incit != git->getStrongerGroups()->end(); incit++){
		if (git->getPriority() >= 0 && (*incit)->getPriority() >= git->getPriority()) 
		    weakerThanOther = true; // another group is stronger (larger alignment coverage) and has no smaller priority
	    }
	
	if (overruled) {
	  //cout << "Group overruled: ";
	  //  git->print(cout);
	}
	if (weakerThanOther) {
	  //cout << "Group weaker than other group.";
	  //  git->print(cout);
	}
	if (git->canCauseAltSplice()){
	  if (!overruled && !weakerThanOther) {
	    //cout << "Group not overruled, not weaker ";
	    // git->print(cout);
	    predictionScheme->addRun(PredictionRun(0, seqlen-1, git->getIncompGroups()));
	  }
	} else {
	  //cout << "group cannot cause altsplice:" << endl;
	  //	  git->print(cout);
	}
    }
    
    // Determine the maximum gene length and number
    int numGenes=0;
    int maxGeneLength = 1;
    for (list<AltGene>::iterator ait = genes->begin(); ait != genes->end(); ait++) {
	numGenes++;
	int len = ait->maxTransEnd() - ait->minTransBegin() + 1;
	if (len > maxGeneLength)
	    maxGeneLength = len;
    }
    //cout << "Have " << numGenes << " genes in the first run. Maximum gene length = " << maxGeneLength << endl;
    int minpadding;
    // padding is larger if we can expect longer genes
    // that is always the case except if we had lots of genes which were all short in the first run.
    minpadding = (4 * 4000 + numGenes * maxGeneLength)/(4+numGenes)/2;
    //cout << "minpadding= " << minpadding << endl;
    

    /* Now go through the list and determine the prediction interval for each PredictionRun R.
     * The prediction interval does not separate genes in 'genes' and the ranges of the groups G(R) are included in the prediction interval.
     */	
    for (list<PredictionRun>::iterator prit = predictionScheme->predictionRuns.begin(); prit != predictionScheme->predictionRuns.end(); prit++) {
#ifdef DEBUG	
	cout << "Determine prediction interval for run ";
	prit->print();
#endif
	int begin, groupsbegin = INT_MAX;
	int end, groupsend = -1;
	// First determine smallest interval that contains all groups in G(R).
	set<HintGroup*> *G = getCausingGroups(*prit);
	for(set<HintGroup*>::iterator Git = G->begin(); Git != G->end(); Git++){
	    if ((*Git)->getBegin() < groupsbegin)
		groupsbegin = (*Git)->getBegin();
	    if ((*Git)->getEnd() > groupsend)
		groupsend = (*Git)->getEnd();
	}
#ifdef DEBUG
	cout << "# smallest interval that contains all groups in G(R): " << groupsbegin << ".." << groupsend << endl;
#endif
	// Second: Increase prediction interval so no genes from the first prediction are separated.
	begin = groupsbegin;
	end = groupsend;
	for (list<AltGene>::iterator ait = genes->begin(); ait != genes->end(); ait++) {
	    if (ait->maxTransEnd() >= groupsbegin && ait->minTransBegin() < begin)
		begin = ait->minTransBegin();
	    if (ait->minTransBegin() <= groupsend && ait->maxTransEnd() > end)
		end = ait->maxTransEnd();
	}

	begin -= minpadding; 
	end += minpadding;

	// TEMP HACK here, let begin and end be a multiple of 10, so tss and tts differ also by multiples of 10
	begin = (begin/10)*10; // TODO sync with UTRmodel 10=UtrModel::ttsSpacing
	end = (end/10+1)*10;
	if (begin < 0)
	    begin = 0;
	if (end > seqlen-1)
	    end = seqlen - 1;
	if (end < begin)
	    throw ProjectError("PredictionRun has empty prediction interval");
	prit->begin = begin;
	prit->end = end;
	delete G;
    }
}

/*
 * SequenceFeatureCollection::setActiveFlag
 */
void SequenceFeatureCollection::setActiveFlag(list<HintGroup*> *groups, bool flag){
    if (!groups)
	return;
    for (list<HintGroup*>::iterator hit = groups->begin(); hit != groups->end(); hit++)
	(*hit)->setActiveFlag(flag);
}


/*
 * SequenceFeatureCollection::joinGenesFromPredRuns
 *
 * The first list of genes of genesOfRuns comes from the run with no omittedGroups.
 */
list<AltGene> *SequenceFeatureCollection::joinGenesFromPredRuns(list<list<AltGene> *> *genesOfRuns, int maxtracks, bool uniqueCDS){
    list<AltGene> *genes;
    list<Transcript*> alltranscripts, allFilteredTranscripts;
    list<Transcript*>::iterator geneit1, geneit2;
    /*
     * Some transcripts need to be deleted.
     * For a PredictionRun R let G(R) be the set of groups g of hints that have R.omittedGroups == {Groups incompatible with g}
     *
     * Delete transcript T (only) of run R iff
     * (T is supported by a group g not in G(R) or T is partial) and
     * T is not supported by a group g in G(R).
     * Partial genes are allowed at the very ends of the sequence.
     * This gets rid of bad transcripts that are likely to occurr in cases like this:
     *
     *  omitted:       <-------                                                                     ---------<
     *  G(R)                                                       <------<
     * but does not get rid of transcripts with mixed hints like
     *  omitted: <-----             ----------                                ----------<
     *           <--------          ----------                                ----------<
     *  G(R)                                                                    ---
     *
     * Transcripts from a mixture of contradictory hint groups are sometimes desired, e.g. in this case:
     *
     *   <-------      -------   --------      --------    -------
     *                                           ----      -----------    -------    ------      ------    ------<
     */
    // Add the first run with all hints to predictionScheme
    predictionScheme->predictionRuns.push_front(PredictionRun(0, seqlen, NULL, true));
    genes = new list<AltGene>;
    int runNo=0;
    int geneNo;
    set<HintGroup*> *G;
    list<PredictionRun>::iterator prit = predictionScheme->predictionRuns.begin();
    list<list<AltGene> *>::iterator generunit;
    for (generunit = genesOfRuns->begin(); generunit != genesOfRuns->end() && prit != predictionScheme->predictionRuns.end(); generunit++, prit++) {
#ifdef DEBUG
	cout << "R= ";
	prit->print();
#endif
	if (*generunit) {
	    if (Constant::alternatives_from_evidence)
		G = getCausingGroups(*prit);
	    else 
		G = NULL;
	    // now go through all transcripts and push it to alltranscripts, if it should not be excluded.
	    geneNo=1;
	    for (list<AltGene>::iterator git = (*generunit)->begin(); git != (*generunit)->end(); git++) {
		// set name according to number of run
		string str = "r" + itoa(runNo) + "g" + itoa(geneNo);
		git->id = str;
		genes->push_back(*git);
		geneNo++;
		for (list<Transcript*>::iterator trit = git->transcripts.begin(); trit != git->transcripts.end(); trit++) {		    
		    
		    if (Gene *gene = dynamic_cast<Gene*> (*trit))
			gene->compileExtrinsicEvidence(groupList);
		    if (Constant::alternatives_from_evidence) {
			// check whether transcript obeys above rule
			bool incomplete = false;
			double sf;
			double bestSfG = 0;
			double bestSfnotG = 0;
			if (!(*trit)->complete && prit->begin != 0 && prit->end != seqlen){
			    incomplete = true;
#ifdef DEBUG
			    cout << "removing incomplete gene of run " << endl;
#endif
			}
			for (list<HintGroup>::iterator grit = groupList->begin(); grit != groupList->end(); grit++){
			    if(!(grit->getEnd() < (*trit)->geneBegin() || grit->getBegin() > (*trit)->geneEnd())) {
#ifdef DEBUG				
				grit->print(cout);//AUSGABE
#endif
				sf = (*trit)->supportingFraction(&*grit);
				if (G->count(&*grit) > 0) {
#ifdef DEBUG
				    cout << "Following HintGroup in G supports above transcript:";
				    grit->print(cout);
#endif
				    if (sf > bestSfG)
					bestSfG = sf;
				} else { 
				    if (sf > bestSfnotG)
					bestSfnotG = sf;
				}
			    }
			}
#ifdef DEBUG
			cout << "bestSfG=" << bestSfG << " bestSfnotG=" << bestSfnotG << " incomplete=" << incomplete << endl;
#endif
			if ((bestSfG >= 0.8 || bestSfG >= bestSfnotG) && !incomplete) {
			    (*trit)->geneid = str;
			    alltranscripts.push_back(*trit);
			} else {
#ifdef DEBUG
			    cout << "transcript excluded: ";
			    if (incomplete)
				cout << "incomplete ";
			    cout << endl;
			    (*trit)->printGFF();
#endif
			}
		    } else {
			(*trit)->geneid = str;
			alltranscripts.push_back(*trit);
		    }
		}
	    }
	    if (G)
		delete G;
	}
	runNo++;
    }
    if (prit != predictionScheme->predictionRuns.end() || generunit != genesOfRuns->end())
	throw ProjectError("joinGenesFromPredRuns::Size of list of genes and number of prediction runs don't match!");
    
    
    
    alltranscripts.sort(ptr_comparison<Transcript>());
    // now remove multiple copies
    for (geneit1 = alltranscripts.begin(); geneit1 != alltranscripts.end(); geneit1++){
       geneit2 = geneit1;
       for (geneit2++; geneit2 != alltranscripts.end() && (*geneit2)->geneBegin() == (*geneit1)->geneBegin();)
	  if (**geneit1 == **geneit2){
	     // delete the transcript with lower meanStateProb
	     if ((*geneit2)->apostprob > (*geneit1)->apostprob) {
		// replace geneit1 with geneit2
		delete *geneit1; // delete the first transcript (deep)
		*geneit1 = *geneit2; // assignment of pointer
	     } else {
		delete *geneit2; // delete the second transcript (deep)
	     }
	     geneit2 = alltranscripts.erase(geneit2); // delete tx pointer only
	  } else
	     geneit2++;
    }
    /*
     * filter transcripts by maximum track number
     */
    Transcript::filterTranscriptsByMaxTracks(alltranscripts, maxtracks);
    genes = groupTranscriptsToGenes(alltranscripts);
    // delete imperfect transcript if there are better ones
    for (list<AltGene>::iterator ait = genes->begin(); ait != genes->end(); ait++) {
	ait->deleteSuboptimalTranscripts(uniqueCDS);
	for (list<Transcript*>::iterator it = ait->transcripts.begin(); it != ait->transcripts.end(); it++) 
	    allFilteredTranscripts.push_back(*it);
    }
    // group transcripts to genes AGAIN, as suboptimal transcripts could have chained two separate genes together
    delete genes;
    genes = groupTranscriptsToGenes(allFilteredTranscripts);
    return genes;
}

/*
 * sortIncompGroupsOfGroups
 */
void SequenceFeatureCollection::sortIncompGroupsOfGroups(){
    // sort by pointers, for '==' operator of lists
    if (!groupList)
	return;
    for (list<HintGroup>::iterator it = groupList->begin(); it != groupList->end(); it++)
	it->sortIncompGroup();
}

/*
 * SequenceFeatureCollection::getCausingGroups
 *
 * Get all hint groups whose incompGroups are identical to omittedGroups: G(R).
 * This assumes that the incompgroups and the omittedGroups are sorted.
 */
set<HintGroup*> *SequenceFeatureCollection::getCausingGroups(PredictionRun &pr){
    list<HintGroup*> *omittedGroups = pr.omittedGroups;
    set<HintGroup*> * G = new set<HintGroup*>;
#ifdef DEBUG
    cout << "G(R)= (" << endl;
#endif
    for (list<HintGroup>::iterator grit = groupList->begin(); grit != groupList->end(); grit++)
	if (pr.allHints || (grit->getIncompGroups() == NULL && omittedGroups == NULL) ||
	    (grit->getIncompGroups() != NULL && omittedGroups != NULL &&
	     *(grit->getIncompGroups()) == *(omittedGroups))){
#ifdef DEBUG
	    grit->print(cout);
#endif
	    G->insert(&*grit);
	}
#ifdef DEBUG
    cout << ")" << endl;
#endif
    return G;
}

/*
 * SequenceFeatureCollection::prepare
 * Prepare for use in predictions.
 */
void SequenceFeatureCollection::prepare(AnnoSequence *annoseq, bool print, bool withEvidence){
    if (Constant::softmasking && withEvidence){
	// check whether RM is a source key
	char *chr = annoseq->sequence;
	unsigned pos = 0;
	//annoseq->offset;
	int start, end;
	while (pos < annoseq->length){
	    while (pos < annoseq->length && !islower(chr[pos])){ pos++; }
	    if (pos < annoseq->length){
		start = end = pos;
		while (end+1 < annoseq->length && islower(chr[end+1])) { end++; }
		Feature rm(start, end, nonexonpartF, bothstrands, -1, "RM");
		rm.seqname = annoseq->seqname;
		rm.source = "softmask";
		rm.feature = "nep";
		rm.score = 0;
		rm.groupname = "";
		rm.priority = -1;
		rm.mult = 1;
		rm.gradeclass = 0;
		collection->setBonusMalus(rm);
		if (rm.bonus != 1.0)
		    addFeature(rm);
		pos = end+1;
	    }
	}
	collection->hasHintsFile = true; // this is used as an indicator to apply malus
    }
    // turn whole sequence to lowercase characters
    for (unsigned pos = 0; pos < annoseq->length; pos++)
	annoseq->sequence[pos] = tolower(annoseq->sequence[pos]);
    setSeqLen(annoseq->length);
    makeGroups();
    checkGroupConsistency(annoseq);
    findGenicGaps();
    determineInterGroupRelations();
    sortFeatureLists();
    computeIndices();
    rescaleBoniByConformance();
    findGroupGaps();
    //cleanRedundantFeatures();
    if (print) {
	cout << "# Constraints/Hints:" << endl;
	// warnInconsistentHints(annoseq);
	printFeatures(cout);
    }
}

/*
 * precompute tables for computing a local exonpart/UTRpart/CDSpart malus later
 */
void SequenceFeatureCollection::prepareLocalMalus(const char* dna){
  int sum;
  cumCovUTRpartPlus.resize(seqlen);
  cumCovUTRpartMinus.resize(seqlen);
  cumCovCDSpartPlus.resize(seqlen);
  cumCovCDSpartMinus.resize(seqlen);
  cumCovExonpartPlus.resize(seqlen);
  cumCovExonpartMinus.resize(seqlen);
  for (int i=0; i<seqlen; i++)
    cumCovUTRpartPlus[i] = cumCovUTRpartMinus[i] = cumCovExonpartMinus[i] = 0;

  sortFeatureLists();
  // UTR, plus
  vector<bool> covered(seqlen, false);
  addCumCov(covered, featureLists[UTRpartF], plusstrand);
  addCumCov(covered, featureLists[exonpartF], plusstrand);
  sum=0;
  for (int i=0; i<seqlen; i++){
    if (!covered[i]) {sum++;}
    cumCovUTRpartPlus[i] = sum;
  }
  
  // UTR, minus
  for (int i=0; i<seqlen; i++)
    covered[i] = false;
  addCumCov(covered, featureLists[UTRpartF], minusstrand);
  addCumCov(covered, featureLists[exonpartF], minusstrand);
  sum=0;
  for (int i=0; i<seqlen; i++){
    if (!covered[i]) {sum++;}
    cumCovUTRpartMinus[i] = sum;
  }

  // CDS, plus
  for (int i=0; i<seqlen; i++)
    covered[i] = false;
  addCumCov(covered, featureLists[CDSpartF], plusstrand);
  addCumCov(covered, featureLists[exonpartF], plusstrand);
  sum=0;
  for (int i=0; i<seqlen; i++){
    if (!covered[i]) {sum++;}
    cumCovCDSpartPlus[i] = sum;
  }
  // CDS, minus
  for (int i=0; i<seqlen; i++)
    covered[i] = false;
  addCumCov(covered, featureLists[CDSpartF], minusstrand);
  addCumCov(covered, featureLists[exonpartF], minusstrand);
  sum=0;
  for (int i=0; i<seqlen; i++){
    if (!covered[i]) {sum++;}
    cumCovCDSpartMinus[i] = sum;
  }
  // Exon, plus
  for (int i=0; i<seqlen; i++)
    covered[i] = false;
  addCumCov(covered, featureLists[exonpartF], plusstrand);
  sum=0;
  for (int i=0; i<seqlen; i++){
    if (!covered[i]) {sum++;}
    cumCovExonpartPlus[i] = sum;
  }
  // Exon, minus
  for (int i=0; i<seqlen; i++)
    covered[i] = false;
  addCumCov(covered, featureLists[exonpartF], minusstrand);
  sum=0;
  for (int i=0; i<seqlen; i++){
    if (!covered[i]) {sum++;}
    cumCovExonpartMinus[i] = sum;
  }
  // splice site local malus
  // mit heisser Nadel fuer rGASP gestrickt
  // uneffizient, nicht ganz genau, Strang unberuecksichtigt
  
  if (hasLocalSSmalus)
    delete [] hasLocalSSmalus;
  hasLocalSSmalus = new Bitmask[seqlen];

  int exwin = 50; // consider evidence in this window on the correct side of the ss as evidence for an exon
  int countthresh = 10; // minimum number of exonparts to count as exon evidence (multiplicities counted)
  Feature *hitlist, *f;
  int hintcount;
  for (int i = 1; i < seqlen-1; i++){
    if (validDSSPattern(dna + i) || validRDSSPattern(dna + i - 1) || validASSPattern(dna + i - 1) || validRASSPattern(dna + i)){
      hitlist =  getFeatureListOvlpingRange(A_SET_FLAG(exonpartF), i-exwin, i+exwin, bothstrands);
      if (validDSSPattern(dna + i) || validRASSPattern(dna + i)){
	f = hitlist;
	hintcount = 0;
	while (f) {
	  if (f->end >= i-exwin && f->end <= i){
	    hintcount += f->mult;
	  }
	  f = f->next;
	}
	if (hintcount >= countthresh) { // have evidence for exon on left side but no hint for SS
	  if (validDSSPattern(dna + i) && !hintedSites[i][forwDSS])
	    hasLocalSSmalus[i].set(forwDSS); // just set flag, malus is a constant for now.
	  if (validRASSPattern(dna + i) && !hintedSites[i][revASS])
	    hasLocalSSmalus[i].set(revASS);
	}
      }
      if (validRDSSPattern(dna + i - 1) || validASSPattern(dna + i - 1)){
	f = hitlist;
        hintcount = 0;
        while (f) {
          if (f->start <= i+exwin && f->start >= i){
            hintcount += f->mult;
          }
          f = f->next;
        }
        if (hintcount >= countthresh) {
          if (validRDSSPattern(dna + i - 1) && !hintedSites[i][revDSS])
            hasLocalSSmalus[i].set(revDSS);
          if (validASSPattern(dna + i - 1) && !hintedSites[i][forwASS])
            hasLocalSSmalus[i].set(forwASS);
        }
      }
      // output for testing
      /*  
      if (hasLocalSSmalus[i][forwDSS])
	cout << "forwDSS localmalus at " << i << endl;
      if (hasLocalSSmalus[i][forwASS])
        cout << "forwASS localmalus at " << i << endl;
      if (hasLocalSSmalus[i][revDSS])
        cout << "revDSS localmalus at " << i << endl;
      if (hasLocalSSmalus[i][revASS])
        cout << "revASS localmalus at " << i << endl;
      */ 
    }
  }
}

/*
 * localSSMalus
 */
double SequenceFeatureCollection::localSSMalus(FeatureType type, int pos, Strand strand){
  if (pos < 0 || pos >= seqlen ||
      (type == dssF && strand == plusstrand && !hasLocalSSmalus[pos][forwDSS]) ||
      (type == dssF && strand == minusstrand && !hasLocalSSmalus[pos][revDSS]) ||
      (type == assF && strand == plusstrand && !hasLocalSSmalus[pos][forwASS]) ||
      (type == assF && strand == minusstrand && !hasLocalSSmalus[pos][revASS]))
    return 1.0;
  return collection->localMalus(type);
}


/*
 * set flag to 1 for all positions covered by at last one feature
 */
void SequenceFeatureCollection::addCumCov(vector<bool> &cov, const list<Feature>& flist, Strand strand){
  for (list<Feature>::const_iterator it=flist.begin(); it!=flist.end(); it++) {
    if (it->strand == strand || it->strand == STRAND_UNKNOWN) {
      for (int i = it->start; i <= it->end; i++)
	if (i>=0 && i<seqlen)
	  cov[i] = true; // make this efficient later by remembering the last interval and then just flag the new regions
    }
  }
}

/*
 * set flag to 1 for all positions covered by at last one feature
 */
int SequenceFeatureCollection::numZeroCov(int start, int end, int type, Strand strand){

  if (end < 0)
    end = 0;
  if (end >= seqlen)
    end = seqlen-1;
  int numZeroCov = 0;
  if (type == UTRpartF){
    if (strand == plusstrand){
      numZeroCov = cumCovUTRpartPlus[end];
      if (start>0)
	numZeroCov -= cumCovUTRpartPlus[start-1];
    } else {
      numZeroCov = cumCovUTRpartMinus[end];
      if (start>0)
	numZeroCov -= cumCovUTRpartMinus[start-1];
    }
  } else if (type == CDSpartF){
    if (strand == plusstrand){
      numZeroCov = cumCovCDSpartPlus[end];
      if (start>0)
	numZeroCov -= cumCovCDSpartPlus[start-1];
    } else {
      numZeroCov = cumCovCDSpartMinus[end];
      if (start>0)
	numZeroCov -= cumCovCDSpartMinus[start-1];
    }
  } else if (type == exonpartF){
    if (strand == plusstrand){
      numZeroCov = cumCovExonpartPlus[end];
      if (start>0)
	numZeroCov -= cumCovExonpartPlus[start-1];
    } else {
      numZeroCov = cumCovExonpartMinus[end];
      if (start>0)
	numZeroCov -= cumCovExonpartMinus[start-1];
    }
  } else {
    // other types not supported
  }
  return numZeroCov;
}


/*
 * PredictionRun::print
 */
void PredictionRun::print(int beginPos){
    if (beginPos<0)
	cout << "PredictionRun=(" << begin << ".." << end;
    else 
	cout << "PredictionRun genomic=(" << begin+beginPos << ".." << end + beginPos;
    if (omittedGroups)
	for (list<HintGroup*>::iterator hit = omittedGroups->begin(); hit != omittedGroups->end(); hit++)
	    cout << ", " << (*hit)->getName();
    cout << ")" << endl;
}

/*
 * PredictionScheme::addRun
 * Add the run if the same run is not already included.
 * 
 */
void PredictionScheme::addRun(PredictionRun run){
    if (!run.omittedGroups)
	return;
    // check for inclusion of same set of omitted groups
    bool alreadyIn = false;
    for (list<PredictionRun>::iterator rit = predictionRuns.begin(); rit != predictionRuns.end() && !alreadyIn; rit++){
	if (*(rit->omittedGroups) == *run.omittedGroups)
	    alreadyIn = true;
    }
    if (!alreadyIn) {
	predictionRuns.push_back(run); 
	//cout << "Add run to scheme: ";
	//run.print();
    }
}

/*
 * PredictionScheme::print
 */
void PredictionScheme::print(int beginPos){
    cout << "PredictionScheme:" << endl;
    for (list<PredictionRun>::iterator rit = predictionRuns.begin(); rit != predictionRuns.end(); rit++)
	rit->print(beginPos);
}
FeatureCollection::FeatureCollection(const FeatureCollection& other){
    // copy all simple data members
    hasHintsFile=other.hasHintsFile;
    offset=other.offset;
    numSeqsWithInfo=other.numSeqsWithInfo;
    numSources=other.numSources;
    sourceKey = new string[numSources];
    individual_liability = new bool[numSources];
    oneGroupOneGene = new bool[numSources];
    for (int i=0; i<numSources; i++) {
	sourceKey[i] = other.sourceKey[i];
	individual_liability[i] = other.individual_liability[i];
	oneGroupOneGene[i] = other.oneGroupOneGene[i];
    }
    for(int i=0; i < NUM_FEATURE_TYPES; i++){
	typeInfo[i]=other.typeInfo[i];
    }
    malustable=NULL;
    localmalustable=NULL;
    // TODO deep copy of all the other stuff
}

bool FeatureCollection::skeyExists(string skey){
    int sourceNum=0;
    while (sourceNum < numSources && sourceKey[sourceNum] != skey)
	sourceNum++;
    if (sourceNum < numSources)
	return true;
    else 
	return false;
}

int FeatureCollection::esource(string skey){
    int sourceNum=0;
    while (sourceNum < numSources && sourceKey[sourceNum] != skey)
	sourceNum++;
    if (sourceNum < numSources)
	return sourceNum;
    else {
	string msg("FeatureCollection::esource: invalid source key: ");
	msg+=skey;
	throw ProjectError(msg);
    }
}

void FeatureCollection::readExtrinsicCFGFile(){ 
    string filename, skey;
    try {
	filename = Properties::getProperty( EXTRFILE_KEY ); // = "extrinsicCfgFile"
    } catch(...) {
	cerr << "Could not find parameter " EXTRFILE_KEY << endl;
	return;
    }
    
    try {
	datei.open(filename.c_str());
	if (!datei) { // second priority, look in config/extrinsic folder
	    string altfname = Properties::getConfigFilename("") + "extrinsic/" + filename.c_str();
	    datei.open(altfname.c_str());
	    if (!datei){
		throw ProjectError(string("Could neither find extrinsic config file ") + filename 
				   + " nor " + altfname +".");
	    }
	}
	datei >> comment;
	readSourceRelatedCFG(datei);
	readTypeInfo(datei);
	datei.close();
	datei.clear();
    } catch (ProjectError e) {
	cerr << e.getMessage() << endl;;
	cerr << "Could not read in file with the configuration of hints: " << filename << endl;
	datei.close();
	datei.clear();
	throw;
    }
}


void FeatureCollection::readTypeInfo(istream& datei){
    string featureName, skey;
    double bonus, malus, localMalus;
    FeatureType type;
    datei >> goto_line_after("[GENERAL]") >> comment;
    while (datei >> comment >> ws, datei && datei.peek() != '[') {
	datei >> featureName >> bonus >> malus;
	type = Feature::getFeatureType(featureName);
	if ((int)type < 0) 
	    cerr << "unknown Feature Type: " << featureName << endl;
	else {
	    // look whether there is a local malus as well after the normal malus, this is optional
	    datei >> localMalus;
	    if (datei.fail()) {// first source key
		datei.clear();
		localMalus = 1.0;
	    }
	    if (localMalus != 1.0){
		if (type == UTRpartF) {
		    cout << "# Setting UTRpart local malus: " << localMalus << endl;
		} else if (type == CDSpartF) {
		    cout << "# Setting CDSpart local malus: " << localMalus << endl;
		} else if (type == dssF) {
		    cout << "# Setting donor splice site local malus: " << localMalus << endl;
		} else if (type == assF) {
		    cout << "# Setting acceptor splice site local malus: " << localMalus << endl;
		} else if (type == exonpartF) {
		    cout << "# Setting exon local malus: " << localMalus << endl;
		} else {
			cerr << "Warning: local malus only supported for UTRpart, CDSpart, exonpart, ass and dss. Not for " << featureName << endl;
		}
		if (localMalus < 0.0 || localMalus > 1.0)
		    throw ProjectError("Local malus must be in the range 0.0-1.0.");
	    }
	    typeInfo[type] = FeatureTypeInfo(numSources, bonus, malus, localMalus);
	    for (int j=0; j < numSources; j++) {
		datei >> skey;
		/*
		 * now read the number of classes,
		 * the numclasses-1 boundaries between different classes of grades
		 * and the numclasses quotients p+(g)/p-(g) for the grades 
		 */
		try {
		    typeInfo[type].read(datei, esource(skey));
		    if (!datei) {
			cerr << "Error in syntax for type " << featureName << " source key " << skey << "." << endl;
			throw ProjectError("Syntax Error");
		    }
		} catch (ProjectError e) {
		    throw ProjectError("FeatureCollection::readExtrinsicCFGFile: " + e.getMessage());
		}
	    }
	}
    }
}
 

void FeatureCollection::readSourceRelatedCFG(istream& datei){
    char buf[256];
    int sourcenum;
    string skey;
    datei >> goto_line_after("[SOURCES]") >> comment;
    datei.getline( buf, 2048 );
    // numSources is the number of non-white space characters first
    numSources=0;
    int i=0;
    while (i< (int) strlen(buf)){
	while(!isalpha(buf[i]) && i < (int) strlen(buf)){
	    i++;
	}
	if (i < (int) strlen(buf)){
	    numSources++;
	    while(isalpha(buf[i])){
		i++;
	    }
	}
    }
    sourceKey = new string[numSources];
    individual_liability = new bool[numSources];
    oneGroupOneGene = new bool[numSources];
    for (int i=0; i<numSources; i++) {
	individual_liability[i] = false;
	oneGroupOneGene[i] = false;
    }
    numSources=0;
    // construct sourceKey array
    i=0;
    while (i< (int) strlen(buf)){
	while(!isalpha(buf[i]) && i < (int) strlen(buf)){
	    i++;
	}
	if (i < (int) strlen(buf)){
	    sourceKey[numSources] = "";
	    while(isalpha(buf[i])){
		sourceKey[numSources] += buf[i]; 
		i++;
	    }
	    numSources++;
	}
    }
    // print sourceKey array
    cout << "# Sources of extrinsic information: ";
    for (int i=0; i<numSources; i++)
 	cout << sourceKey[i] << " ";
    cout << endl;
    // check whether RM is a source key in case that softmasking is enabled
    if(Constant::softmasking && !skeyExists("RM"))
       throw ProjectError("When the softmasking option is turned on, the source RM must be specified in the extrinsic config file.\n");

    // read in other source dependent configurations
    streampos bufpos = datei.tellg();
    datei >> goto_line_after("[SOURCE-PARAMETERS]") >> comment;
    if (!datei) {
	datei.clear();
	datei.seekg(bufpos, ios::beg);
    } else {
	bool done = false;
	while (!done){
	    bufpos = datei.tellg();
	    datei.getline(buf, 255);
	    stringstream stm(buf);
	    stm >> skey;
	    if (stm) {
		if (!skeyExists(skey)) {
		    done = true;
		} else {
		    sourcenum = esource(skey);
		    string option;
		    while(stm >> option) {
			int eqpos = option.find('=');
			if (eqpos == string::npos){ // single string option
			    if (option == "individual_liability"){ // an unsatisfyable hint doesn't touch the other hints of the group
				cout << "# Setting " << option << " for " << skey << "." << endl;
				individual_liability[sourcenum] = true;
			    } else if (option == "1group1gene"){ // try not to predict igenic region in intra-group gaps
				cout << "# Setting " << option << " for " << skey << "." << endl;
				oneGroupOneGene[sourcenum] = true;
			    } else {
				throw ProjectError(string("Unknown option: ") + option + ".");
			    }
			} else { // par=value option
			    string valueStr = option.substr(eqpos+1);
			    option = option.substr(0, eqpos);
			    // example code for par=value option, not used yet
			    if (option == "localMalus"){ // extra malus in a partially covered exon
				//double value = atof(valueStr.c_str());
			    } else {
				throw ProjectError(string("Unknown option: ") + option + ".");
			    }
			}
		    }
		}
	    } else 
		done = true;
	} 
	datei.seekg(bufpos);
    }
}
	
	
void FeatureCollection::readGFFFile(const char *filename){
    /*
     * Read in the configuration file for extrinsic features.
     */    
    int predictionStart, predictionEnd, offset;
    try {
      predictionStart = Properties::getIntProperty( "predictionStart" ) - 1;
    } catch (...) {
      predictionStart = 0;
    }   
    try {
      predictionEnd = Properties::getIntProperty( "predictionEnd" ) - 1;
    } catch (...) {
      predictionEnd = INT_MAX;
    }
    
    if (predictionStart == predictionEnd && predictionStart < 0)
       offset = predictionStart + 1; // left shift input by this amount (sequence is not shifted and cut, though)
    else {
       if (predictionStart < 0)
	predictionStart = 0;
       offset = -predictionStart;
    }

    try {
	datei.open(filename);
	if( !datei ) {
	    cerr << "FeatureCollection::readGFFFile( " << filename << " ) : Could not open the file!!!" << endl;
	    throw ProjectError();
	}

	/*
	 * read in line by line
	 * 
	 */
	Feature f;
	string seqname;
	datei >> comment >> ws;
	while (datei) {
	    try{
		datei >> f >> comment >> ws;
	    } catch (ProjectError e){}
	    if ((f.end >= predictionStart && f.start <= predictionEnd) || predictionStart < 0) {
		f.start += offset;
		f.end += offset;
		setBonusMalus(f);
		if (f.bonus != 1.0){
		    seqname = f.seqname;
		    SequenceFeatureCollection*& psfc = collections[seqname];
		    if (psfc == NULL){
			psfc = new SequenceFeatureCollection(this);
			numSeqsWithInfo++;
		    }
		    psfc->addFeature(f);
		}
	    }
	}
	datei.close();
    } catch (ProjectError e) {
	cerr << e.getMessage() << endl;
	throw e;
    } catch(...) {
	cerr << "FeatureCollection::readGFFFile( " << filename << " ) : Could not read the file!!!" << endl;
	datei.close();
    }
    hasHintsFile = true;
}

/*
 * initialization of a feature with the info given 
 * in the extrinsic config table
 */
void FeatureCollection::setBonusMalus(Feature& f){
    FeatureTypeInfo& fti = typeInfo[f.type];
    int sourcenum = esource(f.esource);
    f.gradeclass = fti.gradeclass(sourcenum, f.score);
    double gradequot = fti.gradequots[sourcenum][f.gradeclass];
    // set the general values if applicable
    if (!(fti.bonus < 0)) { // general bonus/malus
	f.bonus = fti.bonus * gradequot * BONUS_FACTOR;
	f.malus = fti.malus;
	// Let the intron bonus depend on the length
	if (f.type == intronF) {
	    int laenge = f.end - f.start + 1;
	    double intronGeoProb = IntronModel::getGeoProb(); // 1-1/1730 (mal)
	    double igenicGeoProb = IGenicModel::getGeoProb(); // 0.9999;
	    f.bonus *= pow (igenicGeoProb/intronGeoProb, laenge);
	}
	if (f.mult>1){
	    // HACK: this corresponds to adding scores of identical hints, compare deleteEqualElements
	    double newbonus = pow(f.bonus, f.mult);
	    if (newbonus > f.bonus * f.mult)
		newbonus = f.bonus * f.mult;
	    f.bonus = newbonus;
	}
    } else { // individual bonus/malus
	if (f.score>0){
	    f.bonus = f.score;
	}
    }
}
/*
 * get the global malus for exonpart, etc depending on the length
 * precompute values the first time they are needed (for each type separately)
 */
double FeatureCollection::partMalus(FeatureType type, int len){
  if (malus(type) == 1.0 || len <= 0)
    return 1;
  if (!malustable) {
    malustable = new double*[NUM_FEATURE_TYPES];
    for (int t = 0; t < NUM_FEATURE_TYPES; t++)
      malustable[t] = NULL;
  }
  if (malustable[type] == NULL){
    malustable[type] = new double[maxStoreMalus];
    malustable[type][0] = 1.0;
    double cumMalus = 1.0;
    for (int i=1; i<maxStoreMalus; i++){
      cumMalus *= malus(type);
      malustable[type][i] = cumMalus;
    }
  }
  if (len < maxStoreMalus)
    return malustable[type][len];
  else
    return pow(malus(type), len);
}

/*
 * get the local malus for CDSpart and UTRpart depending on the length
 * CDSpart stands for (CDSpartF or exonpartF), UTRpart stands for (UTRpartF or exonpartF)
 * precompute values the first time they are needed (for each type separately)
 * parameter 'bonus' is not used yet. (TODO)
 */
double FeatureCollection::localPartMalus(FeatureType type, int len, Double bonus, int nindep){
  double malusfactor = typeInfo[type].localMalus;
  if (malusfactor == 1.0 || len <= 0)
    return 1.;
  if (!localmalustable) {
    localmalustable = new double*[NUM_FEATURE_TYPES];
    for (int t = 0; t < NUM_FEATURE_TYPES; t++)
      localmalustable[t] = NULL;
  }
  if (localmalustable[type] == NULL){
    localmalustable[type] = new double[maxStoreMalus];
    localmalustable[type][0] = 1.0;
    double cumMalus = 1.0;
    for (int i=1; i<maxStoreMalus; i++){
      cumMalus *= malusfactor;
      localmalustable[type][i] = cumMalus;
    }
    if (cumMalus == 0.0 && malusfactor > 0.0) // use Double instead if this happens
	cerr << "Error: underflow in FeatureCollection::localPartMalus. Please increase malus." << endl;
  }
  if (len < maxStoreMalus)
    return localmalustable[type][len];
  else
    return pow(malusfactor, len);
}

/*
 * determine the number of sequences in the set of sequences annoseq, 
 * for which we have extrinsic information
 */
int FeatureCollection::getNumCommonSeqs(AnnoSequence *annoseq){
    int ret=0;
    while( annoseq ){
	if (collections.count(annoseq->seqname)>0)
	    ret++;
	annoseq = annoseq->next;
    }
    return ret;
}


/*
 * new Version of printAccuracyForSequenceSet
 * prerequisite: genes must all be complete
 *
 */

void FeatureCollection::printAccuracyForSequenceSet(const AnnoSequence* annoseqs, bool cleanRedundancies){
/*
 * Loop over all genes in the set and count for each feature
 * the relative frequency of its occurence at positions where it is correct
 * and the relative frequency of its occurence at positions where it is not correct
 */

    int num_correct[NUM_FEATURE_TYPES];
    int num_incorrect[NUM_FEATURE_TYPES];
    int num_correct_pos[NUM_FEATURE_TYPES];
    double num_incorrect_pos[NUM_FEATURE_TYPES];
    int end_in_exon=0; //number of exonparts ending at coding base
    int numExons, len;
    int sqrExonLen = 0, codingBases = 0, allBases = 0;
   
    const AnnoSequence* curannoseq = annoseqs;
    const Transcript* curtx;
    const Gene* curG;
    list<Feature> flist;
    list<Feature>::iterator it;
    set<Feature>::iterator itms;
    bool correct;
    State* ex;
    // for assessing the meaning of grades
    Matrix< vector<int> > goodG(NUM_FEATURE_TYPES, numSources);
    Matrix< vector<int> > badG(NUM_FEATURE_TYPES, numSources);
    Matrix< vector<double> > gQuots(NUM_FEATURE_TYPES, numSources);

    cout << "Number of sources: " << numSources << endl;
    if (!cleanRedundancies)
	cerr << "#Warning: redundant hints have not been deleted." << endl;
    for (int f=0; f<NUM_FEATURE_TYPES; f++) 
	for (int s=0; s<numSources; s++) {
	    int size = typeInfo[f].gradeclassnums(s);
	    goodG[f][s].assign(size, 0);
	    badG[f][s].assign(size, 0);
	    gQuots[f][s].assign(size, 1.0);
	}

    for (int f=0; f<NUM_FEATURE_TYPES; f++) {
	num_correct[f]=num_incorrect[f]=num_correct_pos[f]=0;
	num_incorrect_pos[f]=0.0;
    }

    while (curannoseq) {
	numExons = 0; 
	int ell_ell = 0;
	if (curannoseq->anno) {
	    SequenceFeatureCollection& sfc = getSequenceFeatureCollection(curannoseq->seqname);
	    // for first assessing the bonus of the features the redundant features are not deleted
	    // but they are cleaned when running (in augustus.cc)
	    if (cleanRedundancies)
		sfc.cleanRedundantFeatures();
		
	    /*
	     * Loop over all genes in the current AnnotationSequence and make lists
	     * of starts, stops, ass, dss
	     */
	    set<Feature> startAS, stopAS, assAS, dssAS;
	    list<Feature> exonAL;

	    for (curtx = curannoseq->anno->genes; curtx != NULL; curtx = curtx->next) {
		if (!curtx->complete)
		    throw ProjectError("printAccuracyForGeneSet called for incompletely annotated gene");
		curG = dynamic_cast<const Gene*> (curtx);
		if (curG) {
		    codingBases += curG->clength;
		    for (ex = curG->exons; ex; ex = ex->next){
			ell_ell += ex->length() * (ex->length() + 1) / 2;
			exonAL.push_back(Feature(ex->begin, ex->end, exonF, curG->strand, stateReadingFrames[ex->type], "annotrain"));
		    }
		    if (curG->strand == plusstrand){
			startAS.insert(Feature(curG->exons->begin, curG->exons->begin+2, startF, plusstrand, 0, "annotrain"));
			stopAS.insert(Feature(curG->codingend-2, curG->codingend, stopF, plusstrand, 0, "annotrain"));
			for (ex = curG->exons; ex != NULL; ex = ex->next) {
			    numExons++;
			    sqrExonLen += ex->length()*(ex->length()+1)/2;
			    if (ex->next) 
				dssAS.insert(Feature(ex->end+1, ex->end+1, dssF, plusstrand, -1, "annotrain"));
			    if (ex != curG->exons)
				assAS.insert(Feature(ex->begin-1, ex->begin-1, assF, plusstrand, -1, "annotrain"));
			}
		    } else if (curG->strand == minusstrand){
			startAS.insert(Feature(curG->codingend-2, curG->codingend, startF, minusstrand, 0, "annotrain"));
			stopAS.insert(Feature(curG->exons->begin, curG->exons->begin+2, stopF, minusstrand, 0, "annotrain"));
			for (ex = curG->exons; ex != NULL; ex = ex->next) {
			    numExons++;
			    sqrExonLen += ex->length()*(ex->length()+1)/2;
			    if (ex->next)
				assAS.insert(Feature(ex->end+1, ex->end+1, assF, minusstrand, -1, "annotrain"));
			    if (ex != curG->exons)
				dssAS.insert(Feature(ex->begin-1, ex->begin-1, dssF, minusstrand, -1, "annotrain"));
			}
		    } else {
			throw ProjectError("Warning: printAccuracyForSequenceSet: have annotation with unknown strand.");
		    }
		}
	    }
	    len = curannoseq->length;
	    allBases += len;

	    num_correct_pos[startF]   += startAS.size();
	    num_incorrect_pos[startF] += 2 * howOftenOccursIt(curannoseq->sequence, STARTCODON) - startAS.size();
	    flist = sfc.getFeatureList(startF);
	    flist.sort();
	    set <Feature> startHS(flist.begin(), flist.end());	
	    set <Feature> goodstartS, badstartS; 
	    // 'start' feature

	    // bad start hints
	    set<Feature> resultS;
	    insert_iterator<set<Feature> > res_ins(resultS, resultS.begin()); 
	    badstartS  = startHS; // start with all hints bad
	    set_difference(badstartS.begin(), badstartS.end(), startAS.begin(), startAS.end(),  res_ins);
		
	    while(badstartS.size() != resultS.size()){
		badstartS = resultS;
		resultS.clear();
		res_ins = inserter (resultS, resultS.begin());
		set_difference(badstartS.begin(), badstartS.end(), startAS.begin(), startAS.end(),  res_ins);
	    } 
	    badstartS = resultS;
	    
	    for (itms = badstartS.begin(); itms != badstartS.end(); itms++) {
		num_incorrect[startF]++;
		badG[startF][esource(itms->esource)][itms->gradeclass]++;
		//cout << "start i " << it->score << endl;
	    } 


	    // good start hints
	    res_ins = inserter (goodstartS, goodstartS.begin());
	    set_difference(startHS.begin(), startHS.end(), badstartS.begin(), badstartS.end(),  res_ins);
	    
	    for (itms = goodstartS.begin(); itms != goodstartS.end(); itms++) {
		num_correct[startF]++;
		goodG[startF][esource(itms->esource)][itms->gradeclass]++;
	    } 
	    
	    // 'stop' feature
	    num_correct_pos[stopF]   += stopAS.size();
	    num_incorrect_pos[stopF] += 2*(howOftenOccursIt(curannoseq->sequence, AMBERCODON)
					   + howOftenOccursIt(curannoseq->sequence, OPALCODON)
					   + howOftenOccursIt(curannoseq->sequence, OCHRECODON)) - stopAS.size();
	    flist = sfc.getFeatureList(stopF);
	    flist.sort();
	    set <Feature> stopHS(flist.begin(), flist.end());
	    set <Feature> goodstopS, badstopS; 
	    
	    // bad stop hints
	    resultS.clear();
	    res_ins = inserter(resultS, resultS.begin()); 
	    badstopS  = stopHS; // start with all hints bad
	    set_difference(badstopS.begin(), badstopS.end(), stopAS.begin(), stopAS.end(),  res_ins);
	    
	    while(badstopS.size() != resultS.size()){
		badstopS = resultS;
		resultS.clear();
		res_ins = inserter(resultS, resultS.begin());
		set_difference(badstopS.begin(), badstopS.end(), stopAS.begin(), stopAS.end(),  res_ins);
	    } 
	    badstopS = resultS;
	    
	    for (itms = badstopS.begin(); itms != badstopS.end(); itms++) {	
		num_incorrect[stopF]++;
		badG[stopF][esource(itms->esource)][itms->gradeclass]++;
		//cout << "stop i " << it->score << endl;
	    } 
	    
	    // good stop hints
	    res_ins = inserter (goodstopS, goodstopS.begin());
	    set_difference(stopHS.begin(), stopHS.end(), badstopS.begin(), badstopS.end(),  res_ins);
	    for (itms = goodstopS.begin(); itms != goodstopS.end(); itms++) {
		num_correct[stopF]++;
		goodG[stopF][esource(itms->esource)][itms->gradeclass]++;
		//cout << "stop c " << it->score << endl;
	    } 
	
	    // 'donor splice site' feature
	    num_correct_pos[dssF]   += dssAS.size();
	    num_incorrect_pos[dssF] += 2*(howOftenOccursIt(curannoseq->sequence, DSS_SEQUENCE)) - dssAS.size();
	    
	    flist = sfc.getFeatureList(dssF);
	    flist.sort();
	    set <Feature> dssHS(flist.begin(), flist.end());
	    set <Feature> gooddssS, baddssS; 
	    
	    // bad dss hints
	    resultS.clear();
	    res_ins = inserter(resultS, resultS.begin()); 
	    baddssS  = dssHS; // start with all hints bad
	    set_difference(baddssS.begin(), baddssS.end(), dssAS.begin(), dssAS.end(),  res_ins);
		
	    while(baddssS.size() != resultS.size()){
		baddssS = resultS;
		resultS.clear();
		res_ins = inserter (resultS, resultS.begin());
		set_difference(baddssS.begin(), baddssS.end(), dssAS.begin(), dssAS.end(),  res_ins);
	    } 
	    baddssS = resultS;
		
	    for (itms = baddssS.begin(); itms != baddssS.end(); itms++) {
		num_incorrect[dssF]++;
		badG[dssF][esource(itms->esource)][itms->gradeclass]++;
		//cout << "dss i " << it->score << endl;
	    } 

	    // good dss hints
	    res_ins = inserter (gooddssS, gooddssS.begin());
	    set_difference(dssHS.begin(), dssHS.end(), baddssS.begin(), baddssS.end(),  res_ins);
		
	    for (itms = gooddssS.begin(); itms != gooddssS.end(); itms++) {
		num_correct[dssF]++;
		goodG[dssF][esource(itms->esource)][itms->gradeclass]++;
		//cout << "dss c " << it->score << endl;
		/* for ss training with predicted genes
		   if (itms->start>42 && itms->start + 42 < len) { // in bounds
		   if (itms->strand == plusstrand)
		   cout << "dss " << string(curannoseq->sequence + itms->start - 40, 0, 82) << endl;
		   else //reverse dss
		   cout << "dss " << reverseComplement(string(curannoseq->sequence + itms->start - 41, 0, 82).c_str()) << endl; 
		   }*/// end (for ss training with predicted genes)
	    }

	    // 'acceptor splice site' feature
	    num_correct_pos[assF]   += assAS.size();
	    num_incorrect_pos[assF] += 2*(howOftenOccursIt(curannoseq->sequence, ASS_SEQUENCE)) - assAS.size();

	    flist = sfc.getFeatureList(assF);
	    flist.sort();
	    set <Feature> assHS(flist.begin(), flist.end());
	    set <Feature> goodassS, badassS; 

	    // bad ass hints
	    resultS.clear();
	    res_ins = inserter(resultS, resultS.begin()); 
	    badassS  = assHS; // start with all hints bad
	    set_difference(badassS.begin(), badassS.end(), assAS.begin(), assAS.end(),  res_ins);
	    
	    while(badassS.size() != resultS.size()){
		badassS = resultS;
		resultS.clear();
		res_ins = inserter (resultS, resultS.begin());
		set_difference(badassS.begin(), badassS.end(), assAS.begin(), assAS.end(),  res_ins);
	    } 
	    badassS = resultS;
	    
	    for (itms = badassS.begin(); itms != badassS.end(); itms++) {
		num_incorrect[assF]++;
		badG[assF][esource(itms->esource)][itms->gradeclass]++;
		//cout << "ass i " << it->score << endl;
	    } 

	    // good ass hints
	    res_ins = inserter (goodassS, goodassS.begin());
	    set_difference(assHS.begin(), assHS.end(), badassS.begin(), badassS.end(),  res_ins);
		
	    for (itms = goodassS.begin(); itms != goodassS.end(); itms++) {
		num_correct[assF]++;
		goodG[assF][esource(itms->esource)][itms->gradeclass]++;
		//cout << "ass c " << it->score << endl;
		/*// for ss training with predicted genes
		  if (itms->start>42 && itms->start + 42 < len) { // in bounds
		  if (itms->strand == plusstrand)
		  cout << "ass " << string(curannoseq->sequence + itms->start - 41, 0, 82) << endl;
		  else //reverse dss
		  cout << "ass " << reverseComplement(string(curannoseq->sequence + itms->start - 40, 0, 82).c_str()) << endl; 
		  }*/// end (for ss training with predicted genes)
	    } 

	    // 'exonpart' feature
	
	    num_correct_pos[exonpartF] += ell_ell;
	    num_incorrect_pos[exonpartF] += (96.7*len*2*3 - ell_ell);  
	    // 96.7 = mean waiting time for stop codon
	    /* malus wird hier ausnahmsweise anders berechnet
	     * P(observe no exonpart in certain exon pos) = 1-#(ending in exon)/#(cod. bases)
	     * P(observe no exonpart in certain noncod. pos) = 1-#(ending not in exon)/#(noncod. bases)
	     */
	    exonAL.sort(); //eigentlich nur nötig, wenn die Exons nicht schon vorher sortiert waren, z.B. wegen verschachtelten Genen

	    flist = sfc.getFeatureList(exonpartF);
	    flist.sort();
	    list<Feature>::iterator firstPossExon = exonAL.begin();
	    int maxexonlen=10000; // assume no exon longer than this
	    for (it = flist.begin(); it != flist.end(); it++) {
		while (firstPossExon->end < it->start-maxexonlen && firstPossExon != exonAL.end())
		    firstPossExon++;
		correct = false;
		bool ending = false;
		/*
		 * check all exons in the range from firstPossExon to the first exon which ends at least 10000 to the right
		 * <--------- 10000 -------- | it | ------ 10000 -------->
		 */
		for (list<Feature>::iterator eit = firstPossExon; eit != exonAL.end() && eit->end < it->end + maxexonlen; eit++){
		    if (eit->start <= it->start && eit->end >= it->end && eit->strand == it->strand &&
			((it->frame == -1 ) || 
			 (it->strand == plusstrand && mod3(3-(eit->frame - (eit->end-eit->start+1))) == mod3(it->frame + it->start - eit->start)) ||
			 (it->strand == minusstrand && mod3(2-eit->frame) == mod3(it->frame + eit->end - it->end )))) 
			correct = true;
		    if (eit->start <= it->end && eit->end >= it->end) 
			ending = true;
		    /* for analysing the overlap of exonparts too large
		       if (ex->begin <= it->start && ex->end < it->end && ex->end >= it->start)
		       cout << "r " << it->end - ex->end << endl;
		       if (ex->begin > it->start && ex->end >= it->end && ex->begin <= it->end)
		       cout << "l " << ex->begin - it->start << endl;*/		
		}
		if (ending)
		    end_in_exon++;
		    
		if (correct) {
		    num_correct[exonpartF]++;		
		    goodG[exonpartF][esource(it->esource)][it->gradeclass]++;
		    //cout << "exonpart c " << it->score << endl;
		} else {
		    num_incorrect[exonpartF]++;
		    badG[exonpartF][esource(it->esource)][it->gradeclass]++;
		    //cout << "exonpart i " << it->score << endl;
		}
	    }

	    // 'exon' feature
	    num_correct_pos[exonF] += numExons;
	    // number of candidate exons is proportional to length
	    num_incorrect_pos[exonF] += (int) (2.811 * len); // expected number of exon candidates

	    flist = sfc.getFeatureList(exonF);
	    flist.sort();
	    firstPossExon = exonAL.begin();
	    for (it = flist.begin(); it != flist.end(); it++) {
		while (firstPossExon != exonAL.end() && firstPossExon->end < it->start)
		    firstPossExon++;
		correct = false;
		for (list<Feature>::iterator eit = firstPossExon; eit != exonAL.end() && eit->end <= it->end; eit++)
		    if (eit->start == it->start && eit->end == it->end && eit->strand == it->strand &&
			((it->frame == -1) ||
			 (it->strand == plusstrand && mod3(3-(eit->frame - (eit->end-eit->start+1))) == it->frame) ||
			 (it->strand == minusstrand && mod3(2-eit->frame) == it->frame)))
			correct = true;
		if (correct){
		    num_correct[exonF]++;
		    goodG[exonF][esource(it->esource)][it->gradeclass]++;
		    // cout << "exon c " << it->score << endl;
		} else {
		    num_incorrect[exonF]++;
		    badG[exonF][esource(it->esource)][it->gradeclass]++;
		    //cout << "exon i " << it->score << endl;
		}
	    }

	    /* 'intron' feature
	    num_correct_pos[intronpartF] += numExons - 1;
	    num_incorrect_pos[intronpartF] += len;       // 
	    if (sfc) {	
		flist = sfc.getFeatureList(intronpartF);
		for (it = flist.begin(); it != flist.end(); it++) {
		    correct = false;
		    for (State *in = curG->introns; in; in = in->next)
			if (in->begin == it->start && in->end == it->end)
			    correct = true;
		    if (correct)
			num_correct[intronpartF]++;
		    else 
			num_incorrect[intronpartF]++;
		}
	    }
	    */
	} // if (curannoseq->anno) 
	curannoseq = curannoseq->next;
    }

    // evaluate good and bad grades
    for (int f=0; f< NUM_FEATURE_TYPES; f++) {
	int numgood=0, numbad=0, numgrades=0;
	for(int s=0;s<numSources;s++) {
	    int size = typeInfo[f].gradeclassnums(s);
	    cout << sourceKey[s] << " ";
	    cout << setw(10) << featureTypeNames[f] << " good: ";
	    for (int k=0; k<size; k++) {
		numgood += goodG[f][s][k];
		numgrades++;
		cout << setw(6) << goodG[f][s][k];
	    }
	    cout << endl << setw(10) << featureTypeNames[f] << "    bad: ";	
	    for (int k=0; k<size; k++) {
		numbad  += badG[f][s][k];
		cout << setw(6) << badG[f][s][k];
	    }
	    cout << endl;
	}
	for(int s=0;s<numSources;s++) {
	    int size = typeInfo[f].gradeclassnums(s);
	    for (int k=0; k<size; k++) {
		if (s==esource("M")) { // manual anchor source type
		    gQuots[f][s][k] = BONUSSUREHINT;
		} else {
		    gQuots[f][s][k] = ((double)goodG[f][s][k]+QUOT_PSEUDOCOUNT)/(numgood+QUOT_PSEUDOCOUNT*numgrades)
			/(((double) badG[f][s][k]+QUOT_PSEUDOCOUNT)/(numbad+QUOT_PSEUDOCOUNT*numgrades));
		}
	    }
	}
    }

    cout << setw(10) << "feature" << " |" << setw(34) << "relative frequency correct" 
	 << " |" << setw(34) << "relative frequency incorrect" << " |" << endl;
    cout << setw(10) << " " << " |" << setw(34) << " " << " |" << setw(34) << " " << " |" << endl;
    for (int f=0; f< NUM_FEATURE_TYPES; f++) {
	cout << setw(10) << featureTypeNames[f] << " |"
	     << setw(8) << num_correct[f] << " /" << setw(10) << num_correct_pos[f] 
	     << setw(14) << setprecision(3) << (double) num_correct[f]/num_correct_pos[f] << " |" 
	     << setw(8) << num_incorrect[f] << " /" << setw(10) << num_incorrect_pos[f] 
	     << setw(14) << setprecision(3) << (double) num_incorrect[f]/num_incorrect_pos[f] << " |"
	     << endl; 

    }
    cout << "sum l(l+1)/2=" << sqrExonLen << endl;
    cout << "# coding bases = 2 x " << codingBases << " = " << 2*codingBases << endl;
    cout << "# bases = 2 x " << allBases << " = " << 2*allBases << endl;
    cout << "#(exonparts ending in exon)= " << end_in_exon << endl;
    double exonpartmalus = (1-((double) end_in_exon / codingBases))/
	( 1-((double) num_correct[exonpartF]+num_incorrect[exonpartF]-end_in_exon)/(allBases-codingBases)); // diese Zeile hinzugefügt 27.07.04

    cout << "--------------------configfile------------------------------" << endl;
    for (int f=0; f< NUM_FEATURE_TYPES; f++) {
	    if (f == intronpartF)
		continue;
	    double b;
	    if (num_incorrect[f] == 0) {
		if (num_correct[f] == 0)
		    b=1.0;
		else 
		    b=1e100;
	    } else
		b = (double) num_correct[f]/num_correct_pos[f] / ((double)num_incorrect[f]/num_incorrect_pos[f]);
	    cout << setw(10) << featureTypeNames[f]
		 << setw(14) << setprecision(5) << b
		 << setw(14) << setprecision(5);
	    if (f != exonpartF)
		cout << (1.0-(double) num_correct[f]/num_correct_pos[f]) / (1.0-(double) num_incorrect[f]/num_incorrect_pos[f]);
	    else 
		cout << exonpartmalus;

	    for(int s=0;s<numSources;s++) {
		int size = typeInfo[f].gradeclassnums(s);
		cout << setw(3) << sourceKey[s] << setw(5) << size;
		for (int klasse=0; klasse < size-1; klasse++)
		    cout << setw(8) << setprecision(2) << typeInfo[f].gradeclassbounds[s][klasse];
		for (int klasse=0; klasse < size; klasse++)
		    cout << setw(8) << setprecision(2) << gQuots[f][s][klasse];
	    }
	    cout << endl;
    }
}


/*
 * TODO: clean for redundancies after cheking the accuracy. Then delete the less confident type.
 * Check the configuration file syntax. E.g. when the number of grades is missing.
 */

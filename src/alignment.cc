/**********************************************************************
 * file:    alignment.cc
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  Generation of exon candidates
 * author:  Mario Stanke
 *
 * date    |   author           |  changes
 * --------|--------------------|------------------------------------------
 * 03.08.13| Mario Stanke       | creation of the file
 * 06.11.14| Lars Romoth        | overhaul of capAliSize
 *         |                    |                                           
 **********************************************************************/

#include "alignment.hh"
#include <string>
#include <iomanip>
#include <queue>
#include <climits>
#include <bitset>
#include <set>
#include <tuple>
#include "lp_lib.h"

using namespace std;

/*
 * constructur using chromosomal position of first nucleotide and row buffer as
 * read from the last column of .maf file
 */
AlignmentRow::AlignmentRow(string seqID, int chrPos, Strand strand, string rowstr) : seqID(seqID), strand(strand), cumFragLen(0) {
    int aliPos = 0;
    int len;
    int i = 0; 
    int alilen = rowstr.size();

    while(i < alilen){
	while (isGap(rowstr[i])){
	    aliPos++;
	    i++;
	}
	len = 0;
	while (i < alilen && !isGap(rowstr[i])){
	    aliPos++;
	    chrPos++;
	    len++;
	    i++;
	}
	if (len>0)
	    addFragment(chrPos - len, aliPos - len, len);
    }
    // row consists of just gaps
    if (frags.empty()){
	addFragment(chrPos, 0, 0); // add a single fragment of length 0 to start of alignment
    }
}

int AlignmentRow::chrStart() const {
    if (frags.empty())
	return -1;
    return frags[0].chrPos;
}

int AlignmentRow::chrEnd() const{
    if (frags.empty()) // this only happens temporarily during appendRow
	return -1;
    fragment last = frags[frags.size()-1];
    return last.chrPos + last.len - 1;
}

int AlignmentRow::aliEnd() const{
    if (frags.empty())
	return -1;
    fragment last = frags[frags.size()-1];
    return last.aliPos + last.len - 1;
}

void AlignmentRow::addFragment(int chrPos, int aliPos, int len){
    frags.push_back(fragment(chrPos, aliPos, len));
    cumFragLen += len;
}

void AlignmentRow::pack(){
    int n = frags.size();
    if (n <= 1)
	return;
    int i=0, // index of next fragment to extend if appropriate
	j=1; // j>i index of next fragment to be extended with if appropriate
    while (j < n){
	// check whether the i-th and j-th fragment can be joined to a single gapless one
	// neither a gap in the alignment nor in the chromosome
	if (frags[i].aliPos + frags[i].len == frags[j].aliPos && frags[i].chrPos + frags[i].len == frags[j].chrPos){
	    frags[i].len += frags[j].len;
	    j++;
	} else {
	    i++;
	    if (i < j)
		frags[i] = frags[j];
	    j++;
	}
    }
    // if successful in packing, above loop leaves a gap at the end of the 'frags' vector, close it
    frags.resize(i+1);
}

// simple left-to-right search starting from given fragment 'from'
// efficient, when many 'chrPos' are searched in left-to-right order
int AlignmentRow::getAliPos(int chrPos, vector<fragment>::const_iterator from){
    if (from == frags.end() || from->chrPos > chrPos) // chrPos the the left of alignment
	return -2;
    while (from != frags.end() && from->chrPos + from->len - 1 < chrPos)
	++from;
    if (from == frags.end())
	return -1;
    // chr          | chrPos 
    //       from->chrPos         
    //         |----------------------|
    // ali   from->aliPos
    if (chrPos < from->chrPos) // chrPos falls in an aligment gap
	return -1;
    return from->aliPos + chrPos - from->chrPos;
}

//same as above, but lets you pass the iterator, which makes multiple use possible
int AlignmentRow::getAliPos(int chrPos, vector<fragment>::const_iterator *from){
    if (*from == frags.end() || (*from)->chrPos > chrPos)
	return -2;
    while (*from != frags.end() && (*from)->chrPos + (*from)->len - 1 < chrPos)
	++*from;
    if (*from == frags.end())
	return -1;
    if (chrPos < (*from)->chrPos)
	return -1;
    return (*from)->aliPos + chrPos - (*from)->chrPos;
}

// convert from alignment to chromosomal position (inverse function of getAliPos()) 
int AlignmentRow::getChrPos(int aliPos, vector<fragment>::const_iterator from){
    if (from == frags.end() || from->aliPos > aliPos) // aliPos to the left of alignment
        return -2;
    while (from != frags.end() && from->aliPos + from->len - 1 < aliPos)
        ++from;
    if (from == frags.end())
        return -1;
    if (aliPos < from->aliPos) // aliPos falls in an aligment gap
        return -1;
    return from->chrPos + aliPos - from->aliPos;
}

/**
 * append row r2 to r1, thereby shifting all alignment coordinates of r2 by aliLen1
 * if signature string (chr and strand) is not empty, treat everything else as missing
 */
void appendRow(AlignmentRow **r1, const AlignmentRow *r2, int aliLen1, string sigstr){
    if (sigstr != "" && *r1 && (*r1)->getSignature() != sigstr)
	*r1 = NULL; // delete left alignment row if not consistent with signature
    if (*r1 == NULL && r2 != NULL && (sigstr == "" || sigstr == r2->getSignature())){
	*r1 = new AlignmentRow();
	(*r1)->seqID = r2->seqID;
	(*r1)->strand = r2->strand;
    }
    if ((*r1) && r2 &&
	((*r1)->seqID != r2->seqID || r2->chrStart() <= (*r1)->chrEnd())){
	// incomplatible sequence IDs or second row does not come after first row,
	// use the fragments with the longer chromosomal range and throw away the other alignment row
	// unless the sigstr parameter breaks the tie
	if ((*r1)->getSeqLen() >= r2->getSeqLen() || (sigstr != "" && sigstr != r2->getSignature()))
	    return; // implicitly delete the fragments of the second row r2
	else {
	    (*r1)->frags.clear(); // delete the fragments of the first row r1
	    (*r1)->seqID = r2->seqID;
	    (*r1)->strand = r2->strand;
	    (*r1)->cumFragLen = 0;
	}
    }
    if (r2 && (sigstr == "" || sigstr == r2->getSignature())) {
	// append fragments of r2, right behind the fragments of r1
	for(vector<fragment>::const_iterator it = r2->frags.begin(); it != r2->frags.end(); ++it)
	    (*r1)->addFragment(it->chrPos, it->aliPos + aliLen1, it->len);
    }
}

ostream& operator<< (ostream& strm, const AlignmentRow &row){
    strm << row.seqID << " " << row.strand << "\t" << row.chrStart() << ".." << row.chrEnd() << "\t" << row.getSeqLen() << "\t" << row.getCumFragLen();
    strm << "\t( ";
    for (vector<fragment>::const_iterator it = row.frags.begin(); it != row.frags.end(); ++it)
        strm << "(" << it->chrPos << ", " << it->aliPos << ", " << it->len << ") ";
    strm << ")";
    return strm;
}

ostream& operator<< (ostream& strm, const Alignment &a){
    strm << a.aliLen << " alignment columns" << endl;
    for (size_t i = 0; i < a.rows.size(); i++){
	strm << "row " << setw(3) << i << "\t";
	if (a.rows[i])
	    strm << *a.rows[i];
	else 
	    strm << "missing";
	strm << endl;
    }
    return strm;
}

// printing a summary of an alignment
// like above variant  "cout << alignment", but with text graphics and spaces instead of tabs
// XXXXXXXXXXXX                                                                                                                                                                                                                                                       
// XXXXXXXXXXXXXXXXXXXXXXXXX                                                                                                                                                                                                                                          
//     XXXXXXXXXXXXXXXXXXXXX 
void Alignment::printTextGraph(ostream& strm){
    const int graphWidth = 100;
    strm << aliLen << " alignment columns" << endl;
    
    // determine maximal sequence id width across rows for new tabular printing
    int maxIdLen = 2;
    for (size_t i = 0; i < rows.size(); i++){
        if (rows[i] && maxIdLen < rows[i]->seqID.length())
	    maxIdLen = rows[i]->seqID.length();
    }
    for (size_t i = 0; i < rows.size(); i++){
      strm << "row " << setw(3) << i << "\t";
      if (rows[i]) {
  	  AlignmentRow &row = *rows[i];
 	  strm << setw(maxIdLen+1) << row.seqID << " " << row.strand << setw(10) << row.chrStart() << setw(10) << row.chrEnd() << setw(10) << row.getSeqLen() << setw(8) << row.getCumFragLen();
	  bitset<graphWidth> aligned;
	  for (vector<fragment>::const_iterator it = row.frags.begin(); it != row.frags.end(); ++it) {
	      for (unsigned j = graphWidth * it->aliPos/aliLen; j < (float) graphWidth * (it->aliPos+it->len) / aliLen && j < graphWidth; j++)
	          aligned[j] = true;
	  }
	  strm << "    |";
	  for (unsigned j = 0; j < graphWidth; j++)
	      strm << (aligned[j]? "X" : " ");
	  strm << "|";
      } else  
	  strm << "missing";
      strm << endl;
    }
}

bool mergeable (Alignment *a1, Alignment *a2, int maxGapLen, float mergeableFrac, bool strong){
    int misMatches = 0;
    int maxNumMatches = 0; // weak case: number of rows, where both alignments have fragments
                           // strong case: number of rows, where at least one alignment has fragments
    // this upper bound is for shortcutting the loop, when there are (very) many species
    int veryMaxMisMatches = (1.0 - mergeableFrac) * a1->numRows();
    for (int s=0; s < a1->numRows() && misMatches <= veryMaxMisMatches; s++){
        if (a1->rows[s] && a2->rows[s]){
	    maxNumMatches++;
	    int dist = a2->rows.at(s)->chrStart() - a1->rows.at(s)->chrEnd() - 1; // 0 distance means direct adjacency
	    if (a1->rows[s]->seqID != a2->rows[s]->seqID || 
		a1->rows[s]->strand != a2->rows[s]->strand || dist < 0 || dist > maxGapLen)
		misMatches++;
	} else if (strong && (a1->rows[s] || a2->rows[s])){
	    // in the strong sense of mergeability, rows in which only 1 alignment is present count as mismatches
	    maxNumMatches++;
	    misMatches++;
	}
    }
    return (misMatches <= ((int) ((1.0 - mergeableFrac) * maxNumMatches)));
}

/*
 *
 *  |      this             |   |      other      |
 *   xxxx xxxx     xxxxx         xxxxxx  xxx xxxx
 *   xxxxxxxxxxxxxxxxxxxxxxx        xxxx   xxxxxxx
 *          NULL                  xxxxxxxxxxxxxxx
 *     xxxxxx xxxxxxxxxxxxxx           NULL
 *        xxxxxxx xxxxxxxxxx     xxxxx         xxx
 *
 *
 */
void Alignment::merge(Alignment *other, const MsaSignature *sig){
    for (size_t s=0; s<numRows(); s++)
	if (sig)
	    appendRow(&rows[s], other->rows[s], aliLen, sig->sigrows[s]);
	else
	    appendRow(&rows[s], other->rows[s], aliLen);
    aliLen += other->aliLen;
}

Alignment* mergeAliList(list<Alignment*> alis,  const MsaSignature *sig){
    if (alis.empty())
	return NULL;
    // for the species not yet fixed by sig, find the signature with the longest alignment length (majority rule)
    // this could be improved, because same-signature alignments can be unmergeable
    list<Alignment*>::iterator it = alis.begin();
    int k = (*it)->numRows();
    MsaSignature *extSig = new MsaSignature(*sig); // extended signature, not-yet-specified rows can be added
    // count the cumulative fragment length for each signature that occurs
    vector<map<string, int> > sigweights(k);
    map<string, int>::iterator sigit;
    for (; it != alis.end(); ++it){
	for (int i=0; i<k; i++){
	    if ((*it)->rows[i]){
		string sigstr = (*it)->rows[i]->getSignature();
		int weight = (*it)->rows[i]->getCumFragLen();
		sigit = sigweights[i].find(sigstr);
		sigweights[i][sigstr] += weight;
	    }
	}
    }
    // determine the extended signature
    for (int i=0; i<k; i++){
	if (sig->sigrows[i] == ""){
	    int curmax = -1;
	    for (sigit = sigweights[i].begin(); sigit != sigweights[i].end(); ++sigit){
		if (sigit->second > curmax){ // new current maximum
		    extSig->sigrows[i] = sigit->first;
		    curmax = sigit->second;
		}
	    }
	    //if (curmax >= 0)
	    //  cout << "extending signature of " << i << " by " << extSig->sigrows[i] << " weight " << curmax << endl;
	}	
    }
    // cout << "ext signature from path:\t" << extSig->sigstr() << endl;
    Alignment *ma = new Alignment(k);
    // merge all alignments with respect to extended signature into a single alignment ma
    for (it = alis.begin(); it != alis.end(); ++it){
	ma->merge(*it, extSig);
    }
    delete extSig;
    // cout << "merged alignment:" << *ma << endl;
    return ma;
}

/**
 * For each alignment in the 'alis' list, if necessary, chunk it into piece alignments
 * that have sequence ranges of size at most maxRange. The new alignments are 
 * added to the end. The maxRange limit can be violated in the untypical case that a single
 * gapless fragment exceed this length threshold.
 * The (recursively chosen) split point is before the largest gap in the interval [2/3*maxRange, maxRange].
 */
void capAliSize(list<Alignment*> &alis, int maxRange){
    typedef list<Alignment*>::iterator AliIt;
    if (alis.empty())
	return; 

    AliIt ait = alis.begin();
    int k = (*ait)->numRows(); // number of species
    for (; ait != alis.end(); ++ait){
	if ((*ait)->maxRange() > maxRange){
	    Alignment *a = *ait;
	    //cout << "splitting up alignment with maxRange= " << (*ait)->maxRange() << endl << *a << endl;
	    /* find a split "point" in a left-to-right pass wrt to alignment positions
	     *                              | candidate boundary of cut (only complete fragments are taken)
	     *                              V
	     * xxxxxxxxx                    xxxxxxxxx             xxxxx
	     * xxx           xxxxxx         xx xxx xxxxxxx        xxxxxx
	     * NULL
	     *  xxxxxx     xxxx               xxxxxxxxxxxx           xxxxx
	     *                              ^            ^        ^         
	     * fragment with leftmost begin |            |        | leftmost begin of frags following
	     *        (held in priority queue)    rightmost end of queued frags  
	     * gapsize is the sum of the genomic sequence gaps over all species
	     */

	    priority_queue<BoundaryFragment, vector<BoundaryFragment>, CompareBoundaryFragment > q;
	    vector<int> bidx(k, -1);
	    vector<int> aliPosCutScore(min(a->aliLen,maxRange + 1),0), maxIntraGapSize(k,0);
	    int maxAliPos = maxRange;
	    vector<pair<size_t,size_t>> thresholds; // <threshold,bonus>
	    thresholds.push_back(make_pair(0,0)); thresholds.push_back(make_pair(maxRange/2,maxRange/2)); thresholds.push_back(make_pair(2*maxRange/3,2*maxRange/3));
	    vector<size_t> actThresholdIndex(k,0);

	    for (size_t s=0; s<k; s++)
		if (a->rows[s] && !a->rows[s]->frags.empty()){
		    q.push(BoundaryFragment(s, a->rows[s]->frags[0].aliPos)); // add first fragment of every alignment row
		}
	    while (!q.empty()){
		// retrieve leftmost boundary fragment
		BoundaryFragment bf = q.top();
		size_t s = bf.first;
		// if the next fragment (already verified) exceeds the maxAliPos (max allowed position) cut the alignment at the position with maximum cutScore and end the while loop by emptying the queue
		if (a->rows[s]->frags[bidx[s] + 1].chrPos - a->rows[s]->frags[0].chrPos + 1 > maxAliPos){
		    //cout << a->rows[s]->seqID << ": " << a->rows[s]->frags[bidx[s] + 1].chrPos << " - " << a->rows[s]->frags[0].chrPos + 1 << " > " << maxRange << endl;

		    /*int testMax = 0;
		      for (size_t spe = 0; spe < k; spe++){
		      if (a->rows[spe] && !a->rows[spe]->frags.empty()){
		      int test = a->rows[spe]->chrEnd() - a->rows[spe]->frags[0].chrPos + 1;
		      if (test > testMax)
		      testMax = test;
		      }
		      }
		      int testLen = a->aliLen;*/

		    // update maxAliPos, if we found a smaller one
		    if (maxAliPos > a->rows[s]->frags[bidx[s] + 1].aliPos - 1){
			maxAliPos = a->rows[s]->frags[bidx[s] + 1].aliPos - 1;
		    }
		    // set every score after maxAliPos to -1 (should not be considered anymore)
		    for (size_t actAliPos = maxAliPos + 1; actAliPos < aliPosCutScore.size(); actAliPos++){
			aliPosCutScore[actAliPos] = -1;             // dont cut here! Because of this we set the score to a negative number 
		    }
		    // search for the maxCutScore and the related cutAliPos
		    int maxCutScore = 0;
		    //int lastScore = 0;
		    size_t cutAliPos = 0;
		    for (size_t actAliPos = 0; actAliPos < aliPosCutScore.size(); actAliPos++){
			//if (/*aliPosCutScore[actAliPos] > 0 &&*/ aliPosCutScore[actAliPos] != lastScore)
			//	cout << actAliPos << ": " << aliPosCutScore[actAliPos] << endl;
			//    lastScore = aliPosCutScore[actAliPos];
			if (aliPosCutScore[actAliPos] >= maxCutScore){
			    maxCutScore = aliPosCutScore[actAliPos];
			    cutAliPos = actAliPos;
			}
		    }

		    /*cout << "CutAlnPart:" << endl;
		      for (size_t sss = 0; sss < a->rows.size(); sss++){
		      if (a->rows[sss]){
		      for (int i = 0; i < a->rows[sss]->frags.size(); i++){
		      if ((i+3 >= a->rows[sss]->frags.size() || (i+3 < a->rows[sss]->frags.size() && a->rows[sss]->frags[i+3].aliPos > cutAliPos)) && (i-3 < 0 || (i-3 >= 0 && a->rows[sss]->frags[i-3].aliPos < cutAliPos)))
		      cout << "CutFragments of " << a->rows[sss]->seqID << ": ChrPos: " << a->rows[sss]->frags[i].chrPos << " bis " << a->rows[sss]->frags[i].chrEnd() << " ; AlnPos: " << a->rows[sss]->frags[i].aliPos << " bis " << a->rows[sss]->frags[i].aliPos + a->rows[sss]->frags[i].len - 1 << endl;
		      }
		      }
		      }*/
		    // with this boundary fragment the alignment becomes too long
		    // create new alignment from initial part up to cutAliPos
		    Alignment *newAli = new Alignment(k);
		    int newAliLen = 0, offset = INT_MAX;
		    for (size_t ss = 0; ss < k; ss++){
			if (a->rows[ss]){
			    // AlignmentRow *&row = newAli->rows[ss];
			    newAli->rows[ss] = new AlignmentRow();
			    newAli->rows[ss]->seqID = a->rows[ss]->seqID;
			    newAli->rows[ss]->strand = a->rows[ss]->strand;
			    vector<fragment>::iterator first = a->rows[ss]->frags.begin(), 
				last = a->rows[ss]->frags.end(), it, itCutPos;

			    // copy the initial part, save the cutPosIterator and test if a hardCut (cut throw fragment) is necessary
			    bool hardCut = false;
			    itCutPos = last;
			    for (it = first; it != last; it++){
				if ((*it).aliPos + (*it).len - 1 <= cutAliPos)
				    newAli->rows[ss]->addFragment(*it);
				else{
				    itCutPos = it;
				    if ((*it).aliPos <= cutAliPos)
					hardCut = true;
				    break;
				}
			    }
			    // do the hardCut
			    if (hardCut){
				newAli->rows[ss]->addFragment(*itCutPos);
				newAli->rows[ss]->frags.back().len = cutAliPos - (*itCutPos).aliPos + 1;
				(*itCutPos).len -= newAli->rows[ss]->frags.back().len;
				(*itCutPos).aliPos = cutAliPos + 1;
				(*itCutPos).chrPos = (*itCutPos).chrPos + newAli->rows[ss]->frags.back().len;
			    }
			    // update newAliLen as maximum of all new alignment positions
			    if (newAli->rows[ss]->aliEnd() > newAliLen)
				newAliLen = newAli->rows[ss]->aliEnd();
			    // delete the fragments just moved to the new alignment from the old alignment a
			    a->rows[ss]->frags.erase(first, itCutPos);
			    a->rows[ss]->setCumFragLen(a->rows[ss]->getCumFragLen() - newAli->rows[ss]->getCumFragLen());
			    if (a->rows[ss]->frags.empty()){
				delete a->rows[ss]; 
				a->rows[ss] = NULL;
			    }
			    // update offset as minimum of all remaining alignment positions
			    if (a->rows[ss] && a->rows[ss]->frags.size() > 0 && offset > a->rows[ss]->frags[0].aliPos)
				offset = a->rows[ss]->frags[0].aliPos;
			}
		    }
		    newAli->aliLen = newAliLen;
		    //cout << "new alignment: " << *newAli << endl;
		    // left-shift all alignment positions in the remaining alignment fragments by offset
		    a->shiftAliPositions(-offset);
		    // append remaining alignment to the end, required further shortening will be done when its turn comes in the main loop


		    /*int test1Max = 0, test2Max = 0;
		      for (size_t spe = 0; spe < k; spe++){
		      if (newAli->rows[spe] && !newAli->rows[spe]->frags.empty()){
		      int test1 = newAli->rows[spe]->chrEnd() - newAli->rows[spe]->frags[0].chrPos + 1;
		      if (test1 > test1Max)
		      test1Max = test1;
		      }
		      if (a->rows[spe] && !a->rows[spe]->frags.empty()){
		      int test2 = a->rows[spe]->chrEnd() - a->rows[spe]->frags[0].chrPos + 1;
		      if (test2 > test2Max)
		      test2Max = test2;
		      }
		      }
		      cerr << a->rows[s]->seqID << ": " << testLen << " = " << newAli->aliLen << " + " << a->aliLen << " (" << newAli->aliLen + a->aliLen << ") " << " --- " << testMax << " >= " << test1Max << " + " << test2Max << endl;*/

		    if (a->numFilledRows() > 1)
			alis.push_back(a);
		    else // remaining alignment empty or has just one row
			delete a;
		    // replace old long alignment with its initial part
		    *ait = newAli;
		    // clear the queue, so the while loop exits and the flow proceeds to the next alignment
		    while (!q.empty())
			q.pop();
		} else {
		    q.pop();
		    // actualize the fragment-index of the actual species to the actual fragment
		    bidx[s]++;
		    // if the end of the actual fragment is bigger than maxRange, then update maxAliPos
		    if (a->rows[s]->frags[bidx[s]].chrEnd() - a->rows[s]->frags[0].chrPos + 1 > maxRange){
			size_t tempMaxAliPos = a->rows[s]->frags[bidx[s]].aliPos + (maxRange - (a->rows[s]->frags[bidx[s]].chrPos - a->rows[s]->frags[0].chrPos + 1));
			if (maxAliPos > tempMaxAliPos)
			    maxAliPos = tempMaxAliPos;
		    }

		    int actFragStart = a->rows[s]->frags[bidx[s]].chrPos - a->rows[s]->frags[0].chrPos + 1;
		    int actFragEnd = a->rows[s]->frags[bidx[s]].chrEnd() - a->rows[s]->frags[0].chrPos + 1;

		    // if this is not the last fragment
		    if (bidx[s] + 1 < a->rows[s]->frags.size()){
			// calculate the size of the gap after this fragment until the next fragment
			size_t tempSize = a->rows[s]->gapLenAfterFrag(bidx[s]);
			// save the maximum gapSize of this species to add this as score for every possible cut position after the last fragment of this species
			if (maxIntraGapSize[s] < tempSize)
			    maxIntraGapSize[s] = tempSize;
			int nextFragStart = a->rows[s]->frags[bidx[s] + 1].chrPos - a->rows[s]->frags[0].chrPos + 1;
			// search for the right actual threshold (normaly this loop will only )
			while (actThresholdIndex[s]+1 < thresholds.size() && nextFragStart >= thresholds[actThresholdIndex[s]+1].first){
			    actThresholdIndex[s]++;
			}
			// from end of actual fragment to the position before the start of the next fragment
			for (size_t actAliPos = a->rows[s]->frags[bidx[s]].aliPos + a->rows[s]->frags[bidx[s]].len - 1; actAliPos < min(maxAliPos + 1, a->rows[s]->frags[bidx[s] + 1].aliPos); actAliPos++){
			    aliPosCutScore[actAliPos] += tempSize + thresholds[actThresholdIndex[s]].second;
			}

			// from start of next fragment to 1 pos before end of next fragment (or maxAliPos)
			int tempChrPos = nextFragStart - a->rows[s]->frags[bidx[s] + 1].aliPos;
			for (size_t actAliPos = a->rows[s]->frags[bidx[s] + 1].aliPos; actAliPos < min(maxAliPos + 1, a->rows[s]->frags[bidx[s] + 1].aliPos + a->rows[s]->frags[bidx[s] + 1].len - 1); actAliPos++){
			    // if the fragment lies on the next threshold or it is already the last threshold, only add bonusScore until one position before the last position in the threshold (because we cut behind this position and for cutting behind the last position we have another score)
			    if (actThresholdIndex[s] + 1 == thresholds.size() || (actThresholdIndex[s] + 1 < thresholds.size() && tempChrPos + actAliPos < thresholds[actThresholdIndex[s] + 1].first)){
				aliPosCutScore[actAliPos] += thresholds[actThresholdIndex[s]].second;
			    }else{break;}
			}
			// add the next fragment of this species and this alignment to the queue
			q.push(BoundaryFragment(s, a->rows[s]->frags[bidx[s] + 1].aliPos));
		    }else{ // if this is the last fragment add the 1 + maximum gap size of this species to score (and dont forget the threshold bonus)
			for (size_t actAliPos = a->rows[s]->frags[bidx[s]].aliPos + a->rows[s]->frags[bidx[s]].len - 1; actAliPos < aliPosCutScore.size(); actAliPos++){
			    aliPosCutScore[actAliPos] += maxIntraGapSize[s] + 1 + thresholds[actThresholdIndex[s]].second;
			}
		    }

		    // if the fragment lies exactly on the threshold we had to add the threshold score for the positions of the fragment who came behind the threshold
		    if (actFragStart < thresholds[actThresholdIndex[s]].first && actFragEnd >= thresholds[actThresholdIndex[s]].first){
			// from fragmentCutPos to one pos before end of fragment (or maxAliPos)
			for (size_t actAliPos = a->rows[s]->frags[bidx[s]].aliPos + (thresholds[actThresholdIndex[s]].first - actFragStart); actAliPos < min(maxAliPos + 1, a->rows[s]->frags[bidx[s]].aliPos + a->rows[s]->frags[bidx[s]].len - 1); actAliPos++){
			    aliPosCutScore[actAliPos] += thresholds[actThresholdIndex[s]].second;
			}
		    }
		}
	    }
	}
    } 
}

struct SegSortCriterion {
    bool operator() (const tuple<int, set<int>, set<int> > &lhs, const tuple<int, set<int>, set<int> > &rhs) {
        return get<0>(lhs) < get<0>(rhs);
    }
};

/*
 * Remove an optimally chosen subset of alignments from alis based on a geneRange overlap criterion.
 * Genomic bases in any species that are covered by some alignment in the input set but not by any in the subset
 * are penalized by uncovPen. Genomic bases that are covered by more than maxCov alignments 
 * in the subset are penalized by covPen for each additional covering alignment.
 * A minimum penalty subset of alignments is chosen and the remaining alinment are removed from alis.
 * Increasing the coverage penalty covPen punishes long overlaps more.
 * Decreasing maxCov punishes multiple overlapping geneRanges more.
 */
void reduceOvlpRanges(list<Alignment*> &alis, size_t maxCov, float covPen){
    if (alis.empty())
	return; 
    float uncovPen = 1.0; // no need to change this

    typedef list<Alignment*>::iterator AliIt;
    AliIt ait = alis.begin();
    int k = (*ait)->numRows(); // number of species
    int n = alis.size(); // number of unfiltered alignments
    /* segments will hold for each species.chr combination a sorted list of nonoverlapping 
     * parts of sequence ranges together with the set of alignment indices that cover it.
     * Before the algorithm below the triple holds 
     *    1. a genomic coordinate
     *    2. a singleton set of an index of an alignment that start at this coord
     *    3. a singleton set of an index of an alignment that ends at this coord
     * Upon finishing of the algorithm the second component is
     *    2. the set of indices of alignments that cover coord up to and excluding the next coord
     */
    map<string, list<tuple<int, set<int>, set<int> > > >  segments;
    
    int i; // alignments are identified by their index in alis list
    // insert alignment sequence ranges into lists
    for (i=0; ait != alis.end(); ++ait, i++){
	Alignment *a = *ait;
	//cout << "i=" << i << "\t";
	//a->printTextGraph(cout);
	for (size_t s=0; s<k; s++){
	    if (a->rows[s]){
		string key = string(itoa(s)) + "." + a->rows[s]->seqID;
		set<int> diffset, emptyset;
		diffset.insert(i);
		// insert the two segment boundaries
		segments[key].push_back(make_tuple(a->rows[s]->chrStart(), diffset, emptyset)); // alignment i starts
		segments[key].push_back(make_tuple(a->rows[s]->chrEnd()+1, emptyset, diffset)); // alignment i ends
	    }
	}
    }

    /* Algorithm to make segment lists overlap-free and nonredundant
     * e.g. turn --------------         into 4 nonoverlapping segments
     *                ---------
     *                   -----------
     */
    for (auto it = segments.begin(); it != segments.end(); ++it){ // loop over species.chr
	string key = it->first;
	auto &segs = it->second; // segment list for this species.chr combination, reference as list is changed
	// sort intervals by genomic coordinate
	segs.sort(SegSortCriterion());
	set<int> curSet; // set of alignment indices that cover the segment starting at current pos

	for (auto lit = segs.begin(); lit != segs.end(); ){
	    curSet.insert(get<1>(*lit).begin(), get<1>(*lit).end()); // add alignment that starts here (if any)
	    // remove (singleton) set of alignments ending here
	    if (!get<2>(*lit).empty())
		curSet.erase(*(get<2>(*lit).begin()));
	    if (next(lit) != segs.end()	&& get<0>(*lit) == get<0>(*next(lit))){ // next segment has same coord
		lit = segs.erase(lit); // remove -> always keep only last segment with this coord
	    } else {
		get<1>(*lit) = curSet;
		++lit;
	    }
	}
	/* DEBUG print alignments covering each nonoverlapping segment
	cout << key << "\t";
	for (auto lit = segs.begin(); lit != segs.end(); ++lit){
	    cout << " " << get<0>(*lit) << ":";
	    for (auto ait = get<1>(*lit).begin(); ait != get<1>(*lit).end(); ++ait)
		cout << *ait << ",";
	    cout << "\t";
	}
	cout << endl;*/
    }

    /* 
     * set up a mixed integer linear program to be solved by the lpsolve library
     */
    
    int R=0; // number of segments, two inequality constraints per segment
    for (auto it = segments.begin(); it != segments.end(); ++it) // loop over species.chr
	for (auto lit = it->second.begin(); next(lit) != it->second.end(); ++lit)	
	    if (!get<1>(*lit).empty())
		R++;
    
    lprec *lp; // linear program object
    int Ncol = n + 2*R;
    // order of variables
    // x_1, ..., x_n    binary x_i=1 iff i-th alignment is chosen in subset
    // y_1, ..., y_R    implicitly binary, y_j = 1 iff j-th segment is not covered (in optimal solution)
    // z_1, ..., z_R    implicitly int, z_j = no of coverages of j-th segment in excess of maxCov (in opt. sol.)

    lp = make_lp(0, Ncol); // model with Ncol columns and yet 0 rows
    if (lp == NULL)
	throw ProjectError("reduceOvlpRanges: Could not construct linear program");
    
    // create space for one row
    int *colno = new int[Ncol];
    REAL *row = new REAL[Ncol];

    int r = 0; // 0 <= r < R running index for segment
    int j; // j running index for constraint
    vector<int> segLen(R); // segment lengths = weights in linear program
    for (auto it = segments.begin(); it != segments.end(); ++it){ // loop over species.chr    
	auto segs = it->second; // segment list for this species.chr combination
	for (auto lit = segs.begin(); next(lit) != segs.end(); ++lit){
	    if (!get<1>(*lit).empty()){ // intermediate segments may also be uncovered by any alignment
		segLen[r] = get<0>(*next(lit)) - get<0>(*lit); // segment length, a multiplier
		// TODO: this may potentially be improved if the weighting uses the actually aligned bases in the segment
		// rather than just the length of the region spanned by the last and first aligned base.
		// cout << "segLen[" << r << "]=" << segLen[r] << endl;
		// add constraints like x2 + x5 + x7 + y1 >= 1 to penalize no coverage (y1 = 1)
		j = 0;
		for (auto ait = get<1>(*lit).begin(); ait != get<1>(*lit).end(); ++ait){
		    colno[j] = *ait + 1;
		    row[j++] = 1;
		}
		colno[j] = n + r + 1;
		row[j++] = 1;
		add_constraintex(lp, j, row, colno, GE, 1);

		// add constraints like x2 + x5 + x7 - z1 <= 3 to penalize too much coverage (z1>0)
		j = 0;
		for (auto ait = get<1>(*lit).begin(); ait != get<1>(*lit).end(); ++ait){
		    colno[j] = *ait + 1;
		    row[j++] = 1;
		}
		colno[j] = n + R + r + 1;
		row[j++] = -1;
		add_constraintex(lp, j, row, colno, LE, maxCov);
		r++;
	    }
	}
    }

    // make x_i's binary variables; x_i = 1 iff i-th alignment is chosen for subset
    for (int i=1; i <= n; i++)
	set_binary(lp, i, TRUE);

    set_add_rowmode(lp, FALSE);
    
    // objective function: sum_j noCovPen * y_j + uncovPen * z_j
    j = 0;
    for (r=0; r<R; r++){
	 colno[j] = n + R + r + 1;
	 row[j++] = covPen * segLen[r];
	 colno[j] = n + r + 1;
	 row[j++] = uncovPen * segLen[r];
    }
    set_obj_fnex(lp, j, row, colno);
    // write_LP(lp, stdout);
    set_verbose(lp, IMPORTANT);
    int ret = solve(lp);
    if (ret == 0){
	cout << "MILP objective=" << get_objective(lp) << endl;
	get_variables(lp, row);
	cout << "removing alignments ";
	for (i=0; i < n; i++)
	    if (row[i]<0.5)
		cout << i << " ";
	cout << endl;
	/*
	  for (r=0; r<R; r++){
	  if (row[n+r] + row[n+R+r] > 0)
	  cout << r+1 << "\ty=" << row[n+r] << "\tz=" << row[n+R+r] << endl; 
	  }*/

	/*
	 * remove alignments from alis that are not in the optimal subset
	 * and compute some statistics
	 */ 
	i = 0;
	int numRemovedAlis = 0;
	int cumFragLenBefore = 0, cumFragLenRemoved = 0;
	for (auto ait = alis.begin(); ait != alis.end(); i++){
	    int cfl = (*ait)->getCumFragLen();
	    cumFragLenBefore += cfl;
	    if (row[i] < 0.5){
		cumFragLenRemoved += cfl;
		ait = alis.erase(ait);
		numRemovedAlis++;
	    } else
		++ait;
	}
	int segLenBefore = 0, segLenRemoved = 0;
	for (r=0; r<R; r++){
	    segLenBefore += segLen[r];
	    if (row[n + r + 1] > 0.5)
		segLenRemoved += segLen[r];
	}

	cout << numRemovedAlis << " geneRanges=alignments (" 
	     << (0.1 * (int) (1000.0 * numRemovedAlis / (numRemovedAlis + alis.size())))
	     << "%) and " << cumFragLenRemoved << " / " << cumFragLenBefore
	     << " = " << (0.1 * (int) (1000.0 * cumFragLenRemoved / cumFragLenBefore))
	     << "% of aligned bases and " << segLenRemoved << " / " << segLenBefore << " = "
	     << (0.1 * (int) (1000.0 * segLenRemoved / segLenBefore))
	     << "% of aligment-range-covered bases were deleted because of too much redundancy." << endl;
    } else { // if this sometimes does not work it is not a big deal
	cerr << "Warning: Could not optimally solve mixed linear integer program. Return value= " << ret << endl;
    }

    // temp for debug, output filtered alignment list
    /*
    cout << "filtered alignments: " << alis.size() << endl;
    for (auto ait = alis.begin(); ait != alis.end(); ++ait){
	(*ait)->printTextGraph(cout);
	cout << endl;
	}*/

    delete [] colno;
    delete [] row;
}

int Alignment::maxRange(){
    int range, max = 0;
    for(size_t s=0; s<numRows(); s++){
	range = rows[s]? rows[s]->getSeqLen() : 0;
	if (range > max)
	    max = range;
    }
    return max;
}

int Alignment::numFilledRows() const {
    int m = 0;
    for (size_t s=0; s<rows.size(); s++)
	if (rows[s])
	    m++;
    return m;
}

int Alignment::getMaxSeqIdLen() const {
    int maxNameLen = 0;
    for(size_t s=0; s<numRows(); s++)
	if (rows[s] && rows[s]->seqID.length() > maxNameLen)
	    maxNameLen = rows[s]->seqID.length();
    
    return maxNameLen;
}

int medianChrStartEndDiff(Alignment *a1, Alignment *a2){
    if (a1->numRows() != a2->numRows())
	return 0;
    vector<int> diffs;
    for(size_t s=0; s<a1->numRows(); s++){
	AlignmentRow *r1 = a1->rows[s], *r2 = a2->rows[s];
	if (r1 && r2 && r1->strand == r2->strand && r1->seqID == r2->seqID){
	    if (r1->chrStart() - r2->chrEnd() >=0)
		diffs.push_back(r1->chrStart() - r2->chrEnd());
	}
    }
    return quantile(diffs, 0.5);
}

string Alignment::getSignature() const {
    string sig = "";
    for(size_t s=0; s<numRows(); s++)
	if (rows[s])
	    sig += itoa(s) + ":" + rows[s]->getSignature();
    return sig;
}

int Alignment::getCumFragLen() const{
    int cumFragLen = 0;
    for(size_t s=0; s<numRows(); s++)
	if (rows[s])
	    cumFragLen += rows[s]->getCumFragLen();
    return cumFragLen;
}

int Alignment::numFitting(const MsaSignature *sig) const{
    int numFitting = 0;
    for (size_t s = 0; s < numRows(); s++){
	if (sig->fits(*this, s))
	    numFitting++;
    }
    return numFitting;
}

void Alignment::shiftAliPositions(int offset){
    for(size_t s=0; s < numRows(); s++)
	if (rows[s]){
	    for (int i=0; i < rows[s]->frags.size(); i++)
		rows[s]->frags[i].aliPos += offset;
	}
    aliLen += offset;
}

// merge pairs of fragments without gap into one fragment making the alignment representation more compact
void Alignment::pack(){
    for(size_t s=0; s < numRows(); s++)
	if (rows[s])
	    rows[s]->pack();
}


int MsaSignature::maxSigStrLen = 0;


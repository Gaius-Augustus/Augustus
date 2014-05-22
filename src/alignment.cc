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
 *         |                    |                                           
 **********************************************************************/

#include "alignment.hh"
#include <string>
#include <iomanip>
#include <queue>
#include <climits>

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
    if (1){
	strm << "\t( ";
	for (vector<fragment>::const_iterator it = row.frags.begin(); it != row.frags.end(); ++it)
	    strm << "(" << it->chrPos << ", " << it->aliPos << ", " << it->len << ") ";
	strm << ")";
    }
    return strm;
}

ostream& operator<< (ostream& strm, const Alignment &a){
    strm << a.aliLen << " alignment columns" << endl;
    for (size_t i = 0; i < a.rows.size(); i++){
	strm << "row " << setw(2) << i << "\t";
	if (a.rows[i]) 
	    strm << *a.rows[i];
	else 
	    strm << "missing";
	strm << endl;
    }
    return strm;
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
	    vector<size_t> curbestidx(k, 0), bidx(k, 0); // the index into the frags vector for each row 
	    int maxgapsize = -1, maxSeqLen = -1;
	    int gapsize = 0;
	    for (size_t s=0; s<k; s++)
		if (a->rows[s] && !a->rows[s]->frags.empty()){
		    q.push(BoundaryFragment(s, a->rows[s]->frags[0].aliPos)); // add first fragment of every alignment row
		    gapsize += a->rows[s]->gapLenAfterFrag(0);
		}
	    
	    while (!q.empty()){
		// retrieve leftmost boundary fragment
		BoundaryFragment bf = q.top();
		size_t s = bf.first;
		if (a->rows[s]->frags[bidx[s]].chrEnd() - a->rows[s]->frags[0].chrPos + 1 > maxRange){
		    // with this boundary fragment the alignment becomes too long
		    // create new alignment from initial part up to curbestidx
		    Alignment *newAli = new Alignment(k);
		    int newAliLen = 0, offset = INT_MAX;
		    //cout << "boundary index: ";
		    //for (size_t ss = 0; ss < k; ss++)
		    //		cout << curbestidx[ss] << " ";
		    //cout << endl << " maxgapsize= " << maxgapsize << endl;
		    for (size_t ss = 0; ss < k; ss++){
			if (a->rows[ss]){
			    // AlignmentRow *&row = newAli->rows[ss];
			    newAli->rows[ss] = new AlignmentRow();
			    newAli->rows[ss]->seqID = a->rows[ss]->seqID;
			    newAli->rows[ss]->strand = a->rows[ss]->strand;
			    vector<fragment>::iterator first = a->rows[ss]->frags.begin(), 
				last = first + curbestidx[ss] + 1, it;
			    if (last > a->rows[ss]->frags.end())
				last = a->rows[ss]->frags.end(); // can be 1 past the end
			    // copy the initial part
			    for (it = first; it != last; it++)
				newAli->rows[ss]->addFragment(*it);
			    // update newAliLen as maximum of all new alignment positions
			    if (newAli->rows[ss]->aliEnd() > newAliLen)
				newAliLen = newAli->rows[ss]->aliEnd();
			    // delete the fragments just moved to the new alignment from the old alignment a
			    a->rows[ss]->frags.erase(first, last);
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
		    if (a->numFilledRows() > 1)
			alis.push_back(a);
		    else // remaining alignment empty of has just one row
			delete a;
		    // replace old long alignment with its initial part
		    *ait = newAli;
		    // clear the queue, so the while loop exits and the flow proceeds to the next alignment
		    while (!q.empty())
			q.pop();
		} else {
		    q.pop();
		    gapsize -= a->rows[s]->gapLenAfterFrag(bidx[s]);
		    bidx[s]++;
		    if (bidx[s] < a->rows[s]->frags.size()){
			gapsize += a->rows[s]->gapLenAfterFrag(bidx[s]);
			int seqlen = a->rows[s]->frags[bidx[s]].chrEnd() - a->rows[s]->frags[0].chrPos + 1;
			if (seqlen > maxSeqLen){
			    if (maxSeqLen < 2*maxRange/3 && seqlen >= 2*maxRange/3 && seqlen <= maxRange){
				// reached 2/3 threshold just now,
				// start determining the best cutoff point from scratch
				maxgapsize = -1; // this will be updated immediately
			    }
			    maxSeqLen = seqlen;
			}
			if (gapsize > maxgapsize){
			    // the cut is made at the widest gap in the last third
			    // unless there is no fragment in the last third, in which case
			    // it is made at the widest overall gap
			    maxgapsize = gapsize;
			    curbestidx = bidx; 
			    //cout << "new curbestidx: ";
			    //for (size_t ss = 0; ss < k; ss++)
			    //	cout << curbestidx[ss] << " ";
			    //cout << endl;
			}
			q.push(BoundaryFragment(s, a->rows[s]->frags[bidx[s]].aliPos));
		    }
		}
	    }
	}
    } 
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


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

void AlignmentRow::addFragment(int chrPos, int aliPos, int len){
    frags.push_back(fragment(chrPos, aliPos, len));
    cumFragLen += len;
}

// simple left-to-right search starting from given fragment 'from'
// effizient, when many 'chrPos' are searched in left-to-right order
int AlignmentRow::getAliPos(int chrPos, vector<fragment>::const_iterator from){
    if (from == frags.end() || from->chrPos > chrPos) // chrPos the the left of alignment
	return -2;
    while (from->chrPos + from->len - 1 < chrPos && from != frags.end())
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
    if (0){
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
    return (misMatches <= ((int) (1.0 - mergeableFrac) * maxNumMatches));
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

int MsaSignature::maxSigStrLen = 0;


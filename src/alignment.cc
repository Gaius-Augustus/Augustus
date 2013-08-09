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
AlignmentRow::AlignmentRow(string seqID, int chrPos, Strand strand, string rowstr) : seqID(seqID), strand(strand) {
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
	    frags.push_back(fragment(chrPos - len, aliPos - len, len));
    }
    // row consists of just gaps
    if (frags.empty()){
	frags.push_back(fragment(chrPos, 0, 0)); // add a single fragment of length 0 to start of alignment
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
 */
void appendRow(AlignmentRow **r1, AlignmentRow *r2, int aliLen1){
    if (*r1 == NULL && r2 != NULL){
	*r1 = new AlignmentRow();
	(*r1)->seqID = r2->seqID;
	(*r1)->strand = r2->strand;
    }
    if ((*r1) && r2 &&
	((*r1)->seqID != r2->seqID || r2->chrStart() <= (*r1)->chrEnd())){
	// incomplatible sequence IDs or second row does not come after first row,
	// use the fragments with the longer chromosomal range and throw away the other alignment row
	if ((*r1)->getSeqLen() >= r2->getSeqLen())
	    return; // implicitly delete the fragments of the second row r2
	else {
	    (*r1)->frags.clear(); // delete the fragments of the first row r1
	    (*r1)->seqID = r2->seqID;
	    (*r1)->strand = r2->strand;
	}
    }
    if (r2) {
	// append fragments of r2, right behind the fragments of r1
	for(vector<fragment>::const_iterator it = r2->frags.begin(); it != r2->frags.end(); ++it)
	    (*r1)->frags.push_back(fragment(it->chrPos, it->aliPos + aliLen1, it->len));
    }
}

ostream& operator<< (ostream& strm, const AlignmentRow &row){
    strm << row.seqID << " " << row.strand << "\t" << row.chrStart() << ".." << row.chrEnd() << "\t" << row.getSeqLen() << "\t( ";
    for (vector<fragment>::const_iterator it = row.frags.begin(); it != row.frags.end(); ++it)
	strm << "(" << it->chrPos << ", " << it->aliPos << ", " << it->len << ") ";
    strm << ")";
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

bool mergeable (Alignment *a1, Alignment *a2, int maxGapLen, float mergeableFrac){
    int remainingAllowedMismatches = (1.0 - mergeableFrac) * a1->numRows();
    for (int s=0; s < a1->numRows() && remainingAllowedMismatches >=0; s++){
        if (a1->rows[s] && a2->rows[s]){
	    int dist = a2->rows.at(s)->chrStart() - a1->rows.at(s)->chrEnd() - 1; // 0 distance means direct adjacency
	    if (dist < 0 || dist > maxGapLen)
		remainingAllowedMismatches--;
	}
    }
    return (remainingAllowedMismatches >=0);
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
void Alignment::merge(Alignment *other){
    for (size_t s=0; s<numRows(); s++)
	appendRow(&rows[s], other->rows[s], aliLen);
    aliLen += other->aliLen;
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

int Alignment::numFilledRows(){
    int m = 0;
    for (size_t s=0; s<rows.size(); s++)
	if (rows[s])
	    m++;
    return m;
}

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
	if (len>0) {
	    frags.push_back(fragment(chrPos - len, aliPos - len, len));
	}
    }
    // row consists of just gaps
    if (frags.empty()){
	frags.push_back(fragment(chrPos, 0, 0)); // add a single fragment of length 0 to start of alignment
    }
}

ostream& operator<< (ostream& strm, const AlignmentRow row){
    strm << row.seqID << " " << row.strand << " ( ";
    for (vector<fragment>::const_iterator it = row.frags.begin(); it != row.frags.end(); ++it)
	strm << "(" << it->chrPos << ", " << it->aliPos << ", " << it->len << ") ";
    strm << ")";
    return strm;
}

bool mergeable (Alignment *a1, Alignment *a2, int maxGapLen, float mergeableFrac){
    int remainingAllowedMismatches = (1.0 - mergeableFrac - 1e5) * a1->numRows;
    for (int s=0; s < a1->numRows && remainingAllowedMismatches >=0; s++){
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
    aliLen += other->aliLen;
}

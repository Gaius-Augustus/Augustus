/**********************************************************************
 * file:    alignment.hh
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  global multiple sequence alignment with efficiently stored long gaps
 * authors: Mario Stanke
 *
 *********************************************************************/

#ifndef _ALIGNMENT
#define _ALIGNMENT

#include "types.hh"

#include <vector>
#include <list>
#include <iostream>

// structure for the reading of an aligned sequence
struct block {
    int begin;
    int length;
    int previousGaps;
};

// gapless alignment fragment
class fragment {
public:
    fragment(int chrPos, int aliPos, int len) : chrPos(chrPos), aliPos(aliPos), len(len) {}
    int chrPos; // chromosomal start position of fragment, 0-based
    int aliPos; // start position of fragment in alignment, 0-based
    int len;    // fragment length
    
};

/**
 * global multiple sequence alignment
 * alignment does not contain the sequence itself, but the info, where gaps are
 */
class AlignmentRow {
public:
    AlignmentRow() {}
    AlignmentRow(string seqID, int chrPos, Strand strand, string rowbuffer);
    ~AlignmentRow(){}
    // remove all these attributes from Alexander later
    int start;    // start chromosomal position of alignment, 1-based
    int seqLen;   // length of aligned sequence fragment, not counting gaps characters
    vector<int*> cmpStarts;
    list<block> sequence;

    int end() {return start + seqLen - 1;} // last aligned position, 1-based
    int chrStart(){
	return frags.at(0).chrPos;
    }
    int chrEnd(){
	size_t n = frags.size();
	fragment last = frags.at(n-1);
	return last.chrPos + last.len - 1;
    }
    int getSeqLen(){
	return chrEnd() - chrStart() + 1;
    }
    friend ostream& operator<< (ostream& strm, const AlignmentRow);

    string seqID; // e.g. chr21
    Strand strand;
    vector<fragment> frags; // fragments are sorted by alignment positions AND by chromosomal positions (assumption for now)
};



/**
 *          1000      1010       1020      1030       1040      1050    chromosomal
 * chr21       |.........|..........|.........|..........|.........|    positions
 *               *********                 ************                 in species 1
 *                \       \               /          /
 * align chr21     *********-------------************------------
 * ment   chr7     ----*****-----------------********************
 *                     |   |                  \    a fragment    \
 *                     *****                   ********************
 * chr7           |.........|..........|.........|..........|.........| chromosomal
 *             2000      2010       2020      2030       2040      2050 positions
 *                                                                      in species 2
 * Here: Alignment a;
 *       a.aliLen = 46
 *       a.rows[0].frags = ((1002, 0, 9), (1027, 23, 12))
 *       a.rows[0].chrStart() = 1002
 *       a.rows[0].chrEnd() = 1037
 */
class Alignment {
public:
    Alignment(size_t k) : aliLen(0), rows(k, NULL), numRows(k) {} // initialize with NULL, which stand for missing AlignmentRows
    ~Alignment(){
	// Steffi: this causes a segmentation fault for more than two species. I don't know why.
	// for (int i=0; i<rows.size(); i++) 
	//    delete rows.at(i);	
    }
    friend bool mergeable (Alignment *a1, Alignment *a2, int maxGapLen, float mergeableFrac);
    void merge(Alignment *other); // append 'other' Alignment to this
public: // should rather be private
    int aliLen; // all aligned sequences are this long when gaps are included
    vector<AlignmentRow*> rows;
    size_t numRows; // = number of species
};

/*
 * b1 and b2 can be merged in that order because they are very similar and right next to each other.
 * In at least 'mergeableFrac' of the alignment rows the aligned sequenes are
 * present, refer to the same terget sequence, are on the same strand and satisfy 0 <= gaplen <= maxGapLen
 */
bool mergeable (Alignment *b1, Alignment *b2, int maxGapLen, float mergeableFrac);

inline bool isGap(char c){
    return (c == '-' || c == '.');
}

#endif  // _ALIGNMENT

/*
 * alignment.hh
 *
 * License: Artistic License, see file LICENSE.TXT or 
 *          https://opensource.org/licenses/artistic-license-1.0
 */

#ifndef _ALIGNMENT
#define _ALIGNMENT

#include "types.hh"

#include <vector>
#include <list>
#include <iostream>
#include <exception>
#include <climits>

class MsaSignature; // forward declaration

/**
 * @brief gapless alignment fragment
 * 
 * @author Mario Stanke
 */
class fragment {
public:
    fragment(int chrPos, int aliPos, int len) : chrPos(chrPos), aliPos(aliPos), len(len) {}
    fragment() {}
    int chrEnd() const { return chrPos + len - 1; }
    int chrPos; // chromosomal start position of fragment, 0-based
    int aliPos; // start position of fragment in alignment, 0-based
    int len;    // fragment length

    /*
    *   added by Giovanna Migliorelli 14.05.2020 
    *   code responsible for serialization
    *   aknowledgement : https://www.boost.org/doc/libs/1_70_0/libs/serialization/doc/tutorial.html
    */

    #ifdef TESTING
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & chrPos; 
        ar & aliPos;
        ar & len; 
    }
    #endif    
};

/**
 * @brief global multiple sequence alignment
 * 
 * @details alignment does not contain the sequence itself, but the info, where gaps are.<br>
 * Generation of exon candidates.
 * 
 * @author Mario Stanke
 */
class AlignmentRow {
public:
    AlignmentRow() : cumFragLen(0) {}
    AlignmentRow(string seqID, int chrPos, Strand strand, string rowbuffer);
    ~AlignmentRow(){}

    int chrStart() const;
    int chrEnd() const;
    int aliEnd() const;
    int getSeqLen() const { return chrEnd() - chrStart() + 1; }
    int getCumFragLen() const { return cumFragLen; }
    void setCumFragLen(int len) { cumFragLen = len;} // use with care to ensure consistency
    int gapLenAfterFrag(size_t i) const {
	if (i+1 >= frags.size())
	    return 0;
	return frags[i+1].chrPos - frags[i].chrPos - frags[i].len;
    }
    void addFragment(int chrPos, int aliPos, int len);
    void addFragment(fragment &f) { addFragment(f.chrPos, f.aliPos, f.len); }
    string getSignature() const {return seqID + ((strand == minusstrand)? "-" : "+");}
    void pack();
    friend ostream& operator<< (ostream& strm, const AlignmentRow &row);
    friend void appendRow(AlignmentRow **r1, const AlignmentRow *r2, int aliLen1, string sigstr);

    /** convert from chromosomal to alignment position
     * start search from the fragment 'from' on, i.e. assume that aliPos is not to the left of fragment *from
     * return -1, if chrPos is outside the range of these fragments
     * return -2, if position is otherwise not mappable: no fragment contains the chrPos, i.e. chrPos is in a gap
     */
    int getAliPos(int chrPos, vector<fragment>::const_iterator from);
    int getAliPos(int chrPos, vector<fragment>::const_iterator *from); // variant from Patrick Balmerth

    int getAliPos(int chrPos) { return getAliPos(chrPos, frags.begin()); }

    // convert from alignment to chromosomal position (inverse function of getAliPos())
    int getChrPos(int aliPos, vector<fragment>::const_iterator from);
    int getChrPos(int aliPos) { return getChrPos(aliPos, frags.begin()); }
    
    // data members
    string seqID; // e.g. chr21
    
    #ifdef TESTING
    int seqIDarchive;   // GM temporarily added to come around some problem with string serialization
    int chrLen;         // GM added to properly set rsa chr lengths during deserialization
    #endif

    Strand strand;
    vector<fragment> frags; // fragments are sorted by alignment positions AND by chromosomal positions (assumption for now)
private:
    int cumFragLen; // total length of all fragments, becomes a nonredundant attribute after merging of alignments

    /*
    *   added by Giovanna Migliorelli 14.05.2020 
    *   code responsible for serialization
    *   Acknowledegment : https://www.boost.org/doc/libs/1_70_0/libs/serialization/doc/tutorial.html
    */

    #ifdef TESTING
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & cumFragLen;
        ar & seqIDarchive; 
        ar & strand;
        ar & frags; 
        ar & chrLen;    
    }
    #endif
};



/**
 * @brief global multiple sequence alignment with efficiently stored long gaps.
 * 
 * @details Generation of exon candidates
 * 
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
 *
 * Coordinates are LEFT TO RIGHT, for reverse strand alignments, they refer to the REVERSE COMPLEMENT of the sequence.
 * 
 * @author Mario Stanke
 */
class Alignment {
public:
    Alignment(size_t k) : aliLen(0), rows(k, NULL) {} // initialize with NULL, which stand for missing AlignmentRows
    ~Alignment(){
	for (unsigned i=0; i<rows.size(); i++) 
	    delete rows[i];	
    }
    friend bool mergeable (Alignment *a1, Alignment *a2, int maxGapLen, float mergeableFrac, bool strong);
    friend ostream& operator<< (ostream& strm, const Alignment &a);
    void printTextGraph(ostream& strm);
    void merge(Alignment *other, const MsaSignature *sig = NULL); // append 'other' Alignment to this
    friend Alignment* mergeAliList(list<Alignment*> alis,  const MsaSignature *sig);
    friend void capAliSize(list<Alignment*> &alis, int maxRange);
    friend void reduceOvlpRanges(list<Alignment*> &alis, size_t maxCov, float covPen);
    int maxRange(); // chromosomal range, maximized over rows
    int numRows() const { return rows.size(); }
    int numFilledRows() const; // number of nonempty rows
    int getCumFragLen() const; // total length of all fragments
    int getMaxSeqIdLen() const;
    friend int medianChrStartEndDiff(Alignment *a1, Alignment *a2);
    string getSignature() const;
    int numFitting(const MsaSignature *sig) const;
    void shiftAliPositions(int offset);
    void pack(); // merge pairs of fragments without gap into one fragment making the alignment representation more compact

    #ifdef TESTING
    void convertAlignment(int r, int start, int end);
    #endif

public: // should rather be private
    int aliLen; // all aligned sequences are this long when gaps are included
    vector<AlignmentRow*> rows;

private:
    /*
    *   added by Giovanna Migliorelli 14.05.2020 
    *   code responsible for serialization and default constructor required in serialization
    *   aknowledgement : https://www.boost.org/doc/libs/1_70_0/libs/serialization/doc/tutorial.html
    */

    #ifdef TESTING
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & aliLen;
        ar & rows;
    }
    
    Alignment(){};
    #endif
};

int medianChrStartEndDiff(Alignment *a1, Alignment *a2);

// sorting operator, with respect to a given species index
struct SortCriterion {
    SortCriterion(size_t speciesIdx) : s(speciesIdx) {};
    bool operator() (Alignment* const a1, Alignment* const a2){
	// alignments, where row s is missing come last
	if (a2->rows[s] == NULL)
	    return true;
	if (a1->rows[s] == NULL)
	    return false;
	if (a1->rows[s]->seqID < a2->rows[s]->seqID)
	    return true;
	if (a1->rows[s]->seqID > a2->rows[s]->seqID)
	    return false;
	// same sequence, compare chromosomal start positions
	return (a1->rows[s]->chrStart() < a2->rows[s]->chrStart());
    }
    size_t s;
};


// sorting operator, with respect to a given species index, index version
// for sorting a vector or list of indices to another vector that holds the Alignments
struct IdxSortCriterion {
    IdxSortCriterion(vector<Alignment*> const a_, size_t speciesIdx) : s(speciesIdx), a(a_) {};
    bool operator() (int i, int j){
	// alignments, where row s is missing come last
	if (a[j]->rows[s] == NULL && a[i]->rows[s]) // this conjunction makes sorting stable where the sequence is missing
	    return true;
	if (a[i]->rows[s] == NULL)
	    return false;
	if (a[i]->rows[s]->seqID < a[j]->rows[s]->seqID)
	    return true;
	if (a[i]->rows[s]->seqID > a[j]->rows[s]->seqID)
	    return false;
	// same sequence, compare chromosomal start positions
	return (a[i]->rows[s]->chrStart() < a[j]->rows[s]->chrStart());
    }
    size_t s;
    vector<Alignment*> const &a;
};

/*
 * b1 and b2 can be merged in that order because they are very similar and right next to each other.
 * In at least 'mergeableFrac' of the alignment rows the aligned sequenes are
 * strong=false: - refer to the same terget sequence, are on the same strand and satisfy 0 <= gaplen <= maxGapLen
 *                 if both are present
 * strong=true:  - refer to the same terget sequence, are on the same strand and satisfy 0 <= gaplen <= maxGapLen
 *                 if at least one is present
 */
bool mergeable (Alignment *b1, Alignment *b2, int maxGapLen, float mergeableFrac, bool strong);

inline bool isGap(char c){
    return (c == '-' || c == '.');
}

/**
 * @brief MsaSignature is a summary of the seqId/strand combinations of the alignment
 *
 * @author Mario Stanke
 */
class MsaSignature {
public:
    string sigstr() const{
	string str;
	for (unsigned s = 0; s < sigrows.size(); ++s)
	    if (sigrows[s] != "")
		str += itoa(s) + ":" + sigrows[s];
	return str;
    }
    vector<string> sigrows; // each row contains seqID and strand, e.g. chr21+
    int numAli;
    int sumAliLen;
    int sumCumFragLen;
    int depth;
    int color;
    bool operator< (const MsaSignature &other) const {
	return (sumCumFragLen > other.sumCumFragLen);
	// before: sorting by 1) number of species and 2) sumAliLen
	// return (depth > other.depth || (depth == other.depth && sumAliLen > other.sumAliLen));
    }
    bool fits(const Alignment &a, size_t s) const {
	return a.rows[s] && (sigrows[s] == a.rows[s]->seqID + strandChar(a.rows[s]->strand));
    }
    static int maxSigStrLen;
};

bool cmpSigPtr(const MsaSignature *s1, const MsaSignature *s2);

typedef pair<size_t, int> BoundaryFragment; // (s, chrPos), where s= species index

/**
 * 
 * @author Mario Stanke
 */
class CompareBoundaryFragment {
public:
    bool operator()(BoundaryFragment& bf1, BoundaryFragment& bf2) {
	return bf1.second > bf2.second; // => sort by increasing chromosomal position
    }
};


/**
 * MsaInsertion, specifies a string to be inserted into a StringAlignment
 */
class MsaInsertion {
public:
    MsaInsertion(size_t s, size_t insertpos, string insert): s(s), insertpos(insertpos), insert(insert){ }
    MsaInsertion() { }
    /**
     * Insertion sorting criterion: Do them right-to-left.
     * Do insertions at the same position from long to short.
     */
    bool operator<(const MsaInsertion& other) const {
        if (insertpos > other.insertpos)
            return true;
        if (insertpos < other.insertpos)
            return false;
        if (insert.length() > other.insert.length())
            return true;
        if (insert.length() < other.insert.length())
            return false;
        if (s < other.s)
            return true;
        return false;
    }
    size_t s; //! index of species
    size_t insertpos; //! at which position to insert
    string insert; //! what to insert
};


/**
 * @brief global multiple sequence alignment in (standard) string representation
 * @author Mario Stanke
 */
class StringAlignment {
public:
    StringAlignment(size_t numrows) : rows(numrows, ""), k(numrows), len(0) {} // initialize with empty strings
    ~StringAlignment(){len = 0;}

    /**
     * Insert unaligned sequences into an alignment.
     * The insertion positions all refer to the MSA before insertions.
     * Example:
     *
     * a-tt-g   insert((0, 0, "ggg"), (2, 4, "gaga"),    ggga-tt-----g
     * --ctgg          (1, 4, "ttt"))                    -----ctttt-gg
     * a-ttgc   ===================================>     ---a-ttgagagc
     *
     * @param[in] insList a list of inserts at certain positions
     */
    void insert(std::list<MsaInsertion> &insList, int maxInsertLen = INT_MAX);

    /**
     * Remove all columns which in which each row that is
     * present has a gap. Example
     * 
     * ggga-tt-----g           gggatt-----g
     * -----ctttt-gg   ====>   ----ctttt-gg
     * ---a-ttgagagc           ---attgagagc
     *
     * @param[in] insList a list of inserts at certain positions
     */
    size_t removeGapOnlyCols();

    bool isGapOnlyCol(size_t col){
        for (size_t s = 0; s < k; ++s)
            if (!rows[s].empty() && rows[s].at(col) != '-')
                return false;
        return true;
    }

    /**
     * Ensure all rows have the same length and compute and set this length.
     */
    void computeLen(){
        size_t m = 0, rowlen;
        for (size_t s = 0; s < k; ++s){
            rowlen = rows[s].length();
            if (rowlen > 0){
                if (m == 0)
                    m = rowlen;
                else if (m != rowlen)
                    throw length_error("StringAlignment with rows of differing lengths");
            }
        }
        len = m;
    }

    friend ostream& operator<< (ostream& strm, const StringAlignment &msa);

// data members
    vector<string> rows;
    size_t k;
    size_t len;
};


#endif  // _ALIGNMENT

/*****************************************************************************\
 * Filename : exoncand.hh
 * Author   : Alexander Gebauer
 *
 * Description: Generation of exon candidates (=possible exons)
 *
 * An exon candidate is a sequence interval with frame and strand information
 * such that 
 * at the 5' end is either an ASS or a start codon and
 * at the 3' end is either a DSS or a stop codon.
 *
 * Date       |   Author              |  Changes
 *------------|-----------------------|------------------------------------------
 * 3.11.2011  | Mario Stanke          | creation of the file
 * 06.12.2011 | Alexander Gebauer     | definition of the stop codons
 * 17.01.2012 | Alexander Gebauer     | add class AlignSeq und struct block
 * 27.02.2012 | Alexander Gebauer     | add class ExonCandidate
 * 13.01.2011 | Mario Stanke          | revision: general set of stop codons
\******************************************************************************/

#ifndef _EXONCAND_HH
#define _EXONCAND_HH

#include <stdint.h>

#include "types.hh"
#include "exonmodel.hh" // for OpenReadingFrame
#define EXON_TYPES 17

enum ExonType {UNKNOWN_EXON = -1,
    // forward strand
    singleGene, initial_0, initial_1, initial_2, internal_0, internal_1, internal_2, terminal_exon,
    // reverse strand
    rsingleGene, rinitial_exon, rinternal_0, rinternal_1, rinternal_2, rterminal_0, rterminal_1, rterminal_2
};

bool isPlusExon(ExonType t);

extern const int exonTypeReadingFrames[EXON_TYPES-1];
extern const char* stateExonTypeIdentifiers[EXON_TYPES-1];

// converts a stateTypeIdentifier to the ExonType
ExonType toExonType(const char* str);


class ExonCandidate {
public:
    ExonCandidate(ExonType s=UNKNOWN_EXON, long int b=0, long int e=0, Double sc=0.0, Double up_sc=1.0, Double down_sc=1.0):
        type(s),
        begin(b),
        end(e),
        score(sc),
        upScore(up_sc),
        downScore(down_sc)
    {}
    ExonCandidate(ExonCandidate* other){
        begin = other->begin;
        end = other->end;
        type = other->type;
        score = other->score;
        upScore = other->upScore;
        downScore = other->downScore;
    }
    ~ExonCandidate(){}
    ExonType type;
    int begin, end;
    Double score;

    int getStart();
    int getEnd();
    Double getUpScore() const {return upScore;}
    Double getDownScore() const {return downScore;}
    void setUpScore(Double s) {upScore = s;}
    void setDownScore(Double s) {downScore = s;}
    void setScore(Double s) {score = s;}

    int getFirstCodingBase();
    int getLastCodingBase();
    int gff3Frame();
    int len() {return end-begin+1;}
    ExonType getExonType();
    Double getScore();
    int complementType();
    StateType getStateType();
    string key();
    int_fast64_t getKey(); // keys encodes all of: chrStart chrEnd type lenMod3
    bool correctType(const char* dna, int dnalen); // verify ExonType on sequence
    friend ostream& operator<<(ostream& strm, const ExonCandidate &ec);
private:
    Double upScore, downScore;
};

/*
 * assqthresh, dssqthresh are between 0 and 1 and thresholds for the inclusion of
 * acceptor/donor splice sites based on the pattern probability
 * assqthresh=0.05 means that only acceptor ss are considered
 * that have a pattern, such that 5% of true splice site patterns have lower probability.
 * The default threshold of 0 means that all splice site patterns are considered.
 */

void findExonCands(map<int_fast64_t, ExonCandidate*> &ecs, map<int_fast64_t, ExonCandidate*> &addECs, const char *dna, int minLen=1, double assmotifqthresh=0.15, double assqthresh=0.3, double dssqthresh=0.7);

//computes the score for the splice sites of an exon candidate
Double computeSpliceSiteScore(Double exonScore, Double minProb, Double maxProb); 

// create new EC from a key encoding all of: chrStart chrEnd type lenMod3
// verification of type, noInFrameStop, etc.
ExonCandidate* create(int_fast64_t key, const char* dna, int dnalen); 

#endif  //  _EXONCAND_HH

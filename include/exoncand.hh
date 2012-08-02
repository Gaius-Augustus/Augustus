/*****************************************************************************\
 * Filename : exoncand.hh
 * Author   : Alexander Gebauer
 *
 * Description: Generation of exon candidates (=possible exons)
 *
 * An exon candidate is a sequence interval with frame and strand information
 * such that 
 * at the 5' end is either an ASS or a start codon and
 * at the 3' end is either a DSS or a stop codon
 *
 * Date       |   Author              |  Changes
 *------------|-----------------------|------------------------------------------
 * 3.11.2011  | Mario Stanke          | creation of the file
 * 06.12.2011 | Alexander Gebauer     | definition of the stopcodons
 * 17.01.2012 | Alexander Gebauer     | add class AlignSeq und struct block
 * 27.02.2012 | Alexander Gebauer     | add class ExonCandidate
\******************************************************************************/

#ifndef _EXONCAND_HH
#define _EXONCAND_HH

#include "types.hh"
#include "exonmodel.hh" // for OpenReadingFrame

#define DECLARE_ON(NAME, PATTERN, COUNT)			\
        inline bool NAME(const char* dna) {				\
    return strncmp(dna, PATTERN, COUNT) == 0;	\
}

DECLARE_ON(ochre,    OCHRECODON, 3)
DECLARE_ON(amber,    AMBERCODON, 3)
DECLARE_ON(opal,     OPALCODON, 3)

#define rCAmber_SEQUENCE "cta"
#define rCOchre_SEQUENCE "tta"
#define rCOpal_SEQUENCE "tca"
DECLARE_ON(onRCOchre,     rCOchre_SEQUENCE, 3)
DECLARE_ON(onRCAmber,    rCAmber_SEQUENCE, 3)
DECLARE_ON(onRCOpal,    rCOpal_SEQUENCE, 3)

inline bool onRCStopcodon(const char* dna) {
    return
            onRCOchre(dna) || onRCOpal(dna) || onRCAmber(dna);
}

#define EXON_TYPES 17

enum ExonType{UNKNOWN_EXON = -1,
    // forward strand
    singleGene, initial_0, initial_1, initial_2, internal_0, internal_1, internal_2, terminal_exon,
    // reverse strand
    rsingleGene, rinitial_exon, rinternal_0, rinternal_1, rinternal_2, rterminal_0, rterminal_1, rterminal_2
};

extern const int exonTypeReadingFrames[EXON_TYPES];
extern const char* stateExonTypeIdentifiers[EXON_TYPES];

// converts a stateTypeIdentifier to the ExonType
ExonType toExonType(const char* str);

// structure for the reading of an aligned sequence
struct block {
    int begin;
    int length;
    int previousGaps;
};

class ExonCandidate {
public:
    ExonCandidate(ExonType s=UNKNOWN_EXON, long int b=0, long int e=0, double sc=0.0, double ass_sc=0.0, double dss_sc=0.0):
        type(s),
        begin(b),
        end(e),
        score(sc),
        assScore(ass_sc),
        dssScore(dss_sc)
    {}
    ~ExonCandidate(){}
    ExonType type;
    int begin, end;
    double score, assScore, dssScore;

    int getStart(void);
    int getEnd(void);
    ExonType getExonType(void);
    double getScore(void);
    int complementType(void);
    StateType getStateType(void);
    string createKey(void);
};

// class stores all the information about an alignment part of a single species
class AlignSeq {
public:
    AlignSeq() {}
    ~AlignSeq(){}
    string name;
    pair<string,long int> seqID; // stores the sequence ID and the length of the sequence ID
    int start, offset, seqLen, alignLen;
    Strand strand;
    vector<int*> cmpStarts;
    list<block> sequence;
};

/*
 * getExonCands: get all exon candidates
 * assqthresh, dssqthresh are between 0 and 1 and thresholds for the inclusion of
 * acceptor/donor splice sites based on the pattern probability
 * assqthresh=0.05 means that only acceptor ss are considered
 * that have a pattern, such that 5% of true splice site patterns have lower probability.
 * The default threshold of 0 means that all splice site patterns are considered.
 */

double getGC_Content(const char *dna);
// void computeIndices(list<ExonCandidate*> cand, int seqlen);

#endif  //  _EXONCAND_HH

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

/*
 * getExonCands: get all exon candidates
 * assqthresh, dssqthresh are between 0 and 1 and thresholds for the inclusion of
 * acceptor/donor splice sites based on the pattern probability
 * assqthresh=0.05 means that only acceptor ss are considered
 * that have a pattern, such that 5% of true splice site patterns have lower probability.
 * The default threshold of 0 means that all splice site patterns are considered.
 */

void getExonCands(const char* dna, float assqthresh=0, float dssqthresh=0);


#endif  //  _EXONCAND_HH

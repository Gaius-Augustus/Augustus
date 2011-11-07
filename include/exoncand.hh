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
\******************************************************************************/

#ifndef _EXONCAND_HH
#define _EXONCAND_HH

#include "types.hh"
#include "exonmodel.hh" // for OpenReadingFrame

/*
 * getExonCands: get all exon candidates
 * assqthresh, dssqthresh are thresholds for the inclusion of 
 * acceptor/donor splice sites based on the pattern probability
 * assqthresh=0.05 means that only acceptor ss are considered
 * that have a pattern, such that 5% of patterns have lower probability.
 *
 */

void getExonCands(const char* dna, float assqthresh=0, float dssqthresh=0);

#endif  //  _EXONCAND_HH

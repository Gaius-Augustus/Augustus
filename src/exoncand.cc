/**********************************************************************
 * file:    exoncand.cc
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  Generation of exon candidates
 * author:  Alexander Gebauer
 *
 * date    |   author           |  changes
 * --------|--------------------|------------------------------------------
 * 03.11.11| Mario Stanke       | creation of the file
 * 06.12.11| Alexander Gebauer  | implementation of getExonCands
 * 22.02.12| Alexander Gebauer  | implementation of readAlignments
 * 29.02.12| Alexander Gebauer  | implementation of summarizeAlignments
 * 27.03.12| Alexander Gebauer  | (implementation of computeIntervals referring to SequenceFeatureCollection::computeIndices)
 **********************************************************************/

#include "exoncand.hh"
#include "intronmodel.hh"
#include "geneticcode.hh"
#include "types.hh"
#include "motif.hh"
#include <fstream>
#include <iostream>
#include <string>
#include <sys/time.h>


const int exonTypeReadingFrames[EXON_TYPES]=
{ 0,0,1,2,0,1,2,0,          // forward exons
  2,2,0,1,2,0,1,2,          // reverse exons
};

const char* stateExonTypeIdentifiers[EXON_TYPES]=
{ "SINGLE", "INITIAL0", "INITIAL1", "INITIAL2", "INTERNAL0", "INTERNAL1", "INTERNAL2", "TERMINAL",
  "RSINGLE", "RINITIAL", "RINTERNAL0", "RINTERNAL1", "RINTERNAL2", "RTERMINAL0", "RTERMINAL1", "RTERMINAL2",
};

int ExonCandidate::getStart() {
    return begin;
}

int ExonCandidate::getEnd() {
    return end;
}

ExonType ExonCandidate::getExonType() {
    return type;
}

Double ExonCandidate::getScore() {
    return score;
}

StateType ExonCandidate::getStateType(){
    if(type >= singleGene && type <= terminal_exon)
        return StateType (type + 1);
    if(type >= rsingleGene && type <= rterminal_2)
        return StateType (type + 28);
    return TYPE_UNKNOWN;
}

// returns the complement biological exon type of a given exon type
int ExonCandidate::complementType() {
    int frontBases = 0;
    if (this->type==0 || this->type==8) {
        return (this->type + 8) % 16;
    } else if (this->type>=1 && this->type<=3) {
        return (this->type + 8 -((this->end - this->begin + 1) %3));
    } else if (this->type>=4 && this->type<=6) {
        frontBases =((this->end - this->begin + 1) - exonTypeReadingFrames[this->type]) %3;
        return (this->type + 8 - frontBases - exonTypeReadingFrames[this->type]);
    } else if (this->type==7) {
        return (this->type + 8 - ((this->end - this->begin + 1) %3));
    } else if (this->type==9) {
        return (this->type - 8 + ((this->end - this->begin + 1) %3));
    } else if (this->type>=10 && this->type<=12) {
        frontBases =(exonTypeReadingFrames[this->type] + ((this->end - this->begin + 2) %3)) %3;
        return (this->type - 6 + frontBases - exonTypeReadingFrames[this->type]);
    } else {
        return (this->type - 8  + ((this->end - this->begin + 1) %3));
    }
}

ExonType toExonType(const char* str){
    int i;
    for (i = 1; i < 9; i++)
        if (strcmp(str, stateTypeIdentifiers[i]) == 0)
            return (ExonType) (i-1);
    for(i = 36; i < 44; i++)
        if (strcmp(str, stateTypeIdentifiers[i]) == 0)
            return (ExonType)  (i-28);
    return UNKNOWN_EXON;
}

// create a key for a map function to find orthologue exons quickly
string ExonCandidate::key() {
    if(this == NULL) {
        return "no candidate";
    } else {
        return (itoa(this->begin) + ":" + itoa(this->end) + ":" + itoa(this->type));
    }
}

double getGC_Content(const char *dna) {
    int n=strlen(dna);
    double gc_content=0.0;
    for (int i=0; i<=n-1; i++) {
        if ((dna[i]=='g') || (dna[i]=='G') || (dna[i]=='C') || (dna[i]=='c'))  {
            gc_content++;
        }
    }
    gc_content=gc_content/n;
    return gc_content;
}


bool comp_test (ExonCandidate* a, ExonCandidate* b) {
    if (a->type != b->type) {
        return (a->type < b->type);
    } else if (a->end != b->end) {
        return (a->end < b->end);
    } else {
        return (a->begin < b->begin);
    }
}



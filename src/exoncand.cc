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
string ExonCandidate::createKey() {
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
    //cout<< "GC content: "<<gc_content<<endl;
    return gc_content;
}

// maybe I do not need it anymore
/*bool compType (ExonCandidate* a, ExonCandidate* b) {
	if (a->type != b->type) {
		return (a->type<b->type);
	} else{
		return (a->end>b->end);
	}
}*/

bool comp_test (ExonCandidate* a, ExonCandidate* b) {
    if (a->type != b->type) {
        return (a->type<b->type);
    } else  if (a->end != b->end) {
        return (a->end<b->end);
    } else {
        return (a->begin<b->begin);
    }
}

// I probably do not need this anymore
// computes intervals to the exon candidates vector so that access to the vector later is quicker
/*void computeIndices(list<ExonCandidate*> cand, int seqlen) {
	    //firstEnd[i] holds an iterator to cand to the first element f in the vector such that
	 	//f->end >= i*block_size, where i>=0 and i ist such that i*block_size <= seqlen
    list<ExonCandidate*>::iterator **firstEnd;
  // lastBegin[i] holds an iterator to cand to the vector such that all following elements f have
  // f->begin > i*block_size, where i>=0 and i ist such that i*block_size <= seqlen
    list<ExonCandidate*>::iterator **lastBegin;
	int block_size=20;  // block_size=20 is for testing, I will change it to 10000 or larger later
	int j;
	if (cand.empty()) {
		return;
	}
	if (firstEnd) {
		for (int type=0; type < NUM_TYPES; type++) {
			delete [] firstEnd[type];
		}
		delete [] firstEnd;
		firstEnd = NULL;
	}
	if (lastBegin) {
		for (int type=0; type < NUM_TYPES; type++) {
			delete [] lastBegin[type];
		}
		delete [] lastBegin;
		lastBegin = NULL;
	}
	if (2*block_size >= seqlen) { // sequence too short
		firstEnd = NULL;
		lastBegin = NULL;
		return;
	}
	int numBlocks = (seqlen/block_size)+1;
	int lastChanged=-1;
	firstEnd = new list<ExonCandidate*>::iterator *[NUM_TYPES];
	lastBegin = new list<ExonCandidate*>::iterator *[NUM_TYPES];
	list<ExonCandidate*>::iterator itb=cand.begin();
	list<ExonCandidate*>::iterator ite=cand.begin();

	for (int type=0; type < NUM_TYPES; type++){
		if (itb==cand.end() || (ite==cand.end())) {
			return;
		} else if (((*itb)->type!=type) || ((*ite)->type!=type)) {
		    firstEnd[type] = lastBegin[type] = NULL; // don't index this statetype because it does not exit in this sequence
		} else {
			firstEnd[type] = new list<ExonCandidate*>::iterator [numBlocks];
			lastBegin[type] = new list<ExonCandidate*>::iterator [numBlocks];

			// create the intervals based on the ends
			lastChanged = -1;
			while ((itb!=cand.end()) && ((*itb)->type == type) ) {
				if ((*itb)->end/block_size > lastChanged) {  // change indices from lastChanged until it->end/block_size
					for (j=lastChanged+1; j <= (*itb)->end/block_size && j<numBlocks; j++) {
						firstEnd[type][j] = itb;
						lastChanged = j;
					}
				}
				itb++;
			}
			//set the rest of the intervals if there is a rest
			while (lastChanged < numBlocks-1) {
				firstEnd[type][++lastChanged] = cand.end();
			}

			// create the indices based on the starts
			// This is more complicated since the vector is sorted by the end positions.
			int leftmostStart = seqlen+1;
			int lastChanged = numBlocks;
			list<ExonCandidate*>::iterator lastUndershootingIt = cand.end();

			// reverse sequence for a statetype so we can iterate from the right to the left without using a reverse_iterator
			cand.sort(compType);
			while ((ite!=cand.end()) && ((*ite)->type == type) ) {
				if ((*ite)->begin < leftmostStart) {
					leftmostStart = (*ite)->begin;
				}
				if (leftmostStart <= (lastChanged-1)*block_size) { // change indices from lastChanged until it->begin/block_size
					for (j=lastChanged-1; j > (leftmostStart-1)/block_size && j>=0; j--) {
						lastBegin[type][j] = ite;
						lastUndershootingIt = ite;
						lastChanged = j;
					}
				}
				ite++;
			}
			cand.sort(comp);
			//set the rest of the intervals if there is a rest
			while (lastChanged > 0) {
				lastBegin[type][--lastChanged] = lastUndershootingIt;
			}

		    cout << " computed indices for the intervals of the exon candidate vector  " << endl;
		    for (int i=0; i<numBlocks; i++) {
		    	cout << " firstEnd[" << type << "][" << i*block_size << "]= ";
		    	if (firstEnd[type][i] != cand.end()) {
		    		cout<< (*(firstEnd[type][i]))->begin<<"  "<<(*(firstEnd[type][i]))->end << endl;
		    	} else {
		    		cout<< "Listenende" << endl;
		    	}
		    	if (lastBegin[type][i] != cand.end()) {
		    		cout<<"lastBegin["<<type<<"]["<<i*block_size<<"]= "<< (*(lastBegin[type][i]))->begin<<"  "<<(*(lastBegin[type][i]))->end<<endl;
		    	} else {
		    		cout<< "lastBegin[" << type << "][" << i*block_size << "]= " << "Listenende" << endl;
		    	}

		    }
		}
	}
}*/

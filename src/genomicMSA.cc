/**********************************************************************
 * file:    genomicMSA.cc
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  Generation of exon candidates
 * author:  Alexander Gebauer
 *
 * date    |   author           |  changes
 * --------|--------------------|------------------------------------------
 * 04.04.12| Alexander Gebauer  | creation of the file
 * 01.07.13| Mario Stanke       | overhaul of readAlignment
 **********************************************************************/

#include "exoncand.hh"
#include "genomicMSA.hh"
#include "orthoexon.hh"
#include "phylotree.hh"
#include <fstream>
#include <iostream>
#include <string>

/*  *.maf file is read and saved into this->alignment
 *  only the species with names in the given list are considered, the rest is ignored
 */
void GenomicMSA::readAlignment(string alignFilename) {
    int index = 0;
    string rowseq, buffer;
    string completeName;
    char strandChar;
    block cur_block;
    AlignmentRow *row;
    Alignment *alignBlock;
    map<string, size_t> notExistingSpecies;
    int lenOfChr;
    string speciesName;
    string seqID;
    int chrStart;
    int seqLen; 
    Strand strand;

    numSpecies = rsa->getNumSpecies();
    ifstream Alignmentfile;
    Alignmentfile.open(alignFilename.c_str(), ifstream::in);
    if (!Alignmentfile) {
        cerr << "Could not find the alignment file " << alignFilename << "." << endl;
        throw PropertiesError( "GenomicMSA::readAlignment: Could not open this file!" );
    }

    while (!Alignmentfile.eof()) {
        Alignmentfile >> buffer;
        int numSpeciesFound = 0;
        if (buffer == "s") {
	    alignBlock = new Alignment(numSpecies); // create new empty alignment block
            // loop over the lines of the alignment block, don't search any further if enough species were found
            while (numSpeciesFound < numSpecies && buffer == "s") {
		// reads the name of the species and the seqID from the first column
		// TODO: may have to be checked/adjusted for general .maf files
		Alignmentfile >> completeName >> chrStart >> seqLen >> strandChar >> lenOfChr >> rowseq >> buffer;
		
		// split species name and sequence ID
		for (int i=0; i<completeName.length(); i++) {
		    // seperator is the point '.' for example hs19.chr21, has to be changed
		    if ((completeName[i] == '-') || (completeName[i] == '.')) { 
			speciesName = completeName.substr(0,i);
			seqID = completeName.substr(i+1, string::npos);
			// some input file have a suffix "(..)" that needs to be stripped
			string::size_type p = seqID.find_first_of("(");
			if (p != std::string::npos)
			    seqID = seqID.erase(p); 
			break;
		    }
		    if (i == completeName.length()-1) {
			speciesName = completeName;
			seqID = "unknown";
		    }
		}
		  
		if (strandChar == '+') {
		    strand = plusstrand;
		} else if (strandChar == '-') {
		    strand = minusstrand;
		} else {
		    strand = STRAND_UNKNOWN;
		}
		if (!alignBlock->aliLen)
		    alignBlock->aliLen = rowseq.length();
		else if (alignBlock->aliLen != rowseq.length()) {
		    throw ProjectError("Error in MAF in sequence " + seqID + " at position " + itoa(chrStart) 
				       + ". Alignment row does not agree in length.");
		}
		
		row = new AlignmentRow (seqID, chrStart, strand, rowseq);
		row->start = chrStart + 1; // TODO remove
		row->cmpStarts.push_back(&row->start);
		row->seqLen = seqLen;
		if (seqLen != row->getSeqLen())
		    cerr << "Inconsistenty in .maf file: Sequence length different from number of non-gap characters in row:" 
			 << endl << "speciesName" << "." << seqID << "\t" << chrStart << "\t" << seqLen << "\t" << rowseq << endl;

		// reads the aligned sequence
		cur_block.begin = row->start;
		cur_block.length = 1;
		cur_block.previousGaps = 0;
		for (int i=1; i <= rowseq.length(); i++) {
		    if ((rowseq[i-1] != '-') && ( rowseq[i] == '-')) {
			row->sequence.push_back(cur_block);
		    } else if ((rowseq[i-1]=='-')&&(rowseq[i]=='-')) {
			cur_block.previousGaps++;
		    } else if ((rowseq[i-1]=='-')&&(rowseq[i]!='-')) {
			cur_block.begin = row->start+i;
			cur_block.length=1;
			cur_block.previousGaps++;
		    } else if ((rowseq[i-1]!='-')&&(rowseq[i]!='-')) {
			if (i==rowseq.length()) {
			    row->sequence.push_back(cur_block);
			} else {
			    cur_block.length++;
			}
		    }
		}
		
		index = rsa->getIdx(speciesName);
		if (index >=0) { // species name in the white list
		    alignBlock->rows[index] = row; // place at the right position
		    // store chrLen and check whether consistent with previous chrLen
		    try{
			rsa->setLength(index, row->seqID, lenOfChr);
		    } catch (ProjectError e){
			cerr << e.getMessage() << endl << "MAF file inconsistent." << endl;
			throw e;
		    }
		    numSpeciesFound++;
		} else {
		    notExistingSpecies.insert(pair<string, size_t>(speciesName, 1));
		}
	    }
	    alignment.push_back(alignBlock);
	}
    }
    // clean up
    Alignmentfile.close();
    for (map<string,size_t>::iterator it = notExistingSpecies.begin(); it != notExistingSpecies.end(); ++it)
	cerr << "Warning: Species " << it->first << "is not included in the target list of species. These alignment lines are ingored." << endl;
}

/* printMAF
 * print alignment in .maf format, to stdout if outFname is empty string
 */
void GenomicMSA::printAlignment(string outFname){
    Alignment *aliblock;
    for (list<Alignment*>::iterator alit = alignment.begin(); alit != alignment.end(); alit++) {
	aliblock = *alit;
	cout << *aliblock << endl << endl;
	/*for(vector<AlignmentRow*>::iterator aseqit = aliblock->rows.begin(); aseqit != aliblock->rows.end(); aseqit++){
	    AlignmentRow* row = *aseqit;
	    if (row) {
		cout << setw(15) << row->seqID
		     << "\tstart=" << row->start
		     << "\tseqLen=" << row->seqLen
		     << "\t" << ((row->strand == plusstrand)? "+":"_");
		cout << "\tstarts:";
		for (vector<int *>::iterator startsit = row->cmpStarts.begin(); startsit != row->cmpStarts.end(); ++startsit){
		    cout << **startsit;
		    if (startsit+1 != row->cmpStarts.end())
			cout << ",";
		}
		cout << "\tsequence:";
		for (list<block>::iterator bit = row->sequence.begin(); bit != row->sequence.end();){
		    cout << "(" << bit->begin << "," << bit->length << "," << bit->previousGaps << ")";
		    if (++bit != row->sequence.end())
			cout << ";";
		}
	    } else {
		cout << "no alignment";
	    }
	    cout << endl;
	}
	cout << endl;*/
    }
}

/**
 * Merges pairs of alignments in trivial cases in order to reduce the alignment number without doing any pontentially false mergers.
 * The priority is on avoiding false mergers. Not all trivial mergers are guaranteed to be realized, though.
 * Trivial mergers are found by sorting with respect to a number of (random) different species.
 */
void GenomicMSA::compactify(){
    int numSortRefs = 3; // Sort wrt to at most this many species, a larger value requires more time now
                         // but may save time later by reducing alignments. For genomic alignments, where the first species is 
                         // present in ALL alignments (e.g. UCSC vertebrate), this only needs to be 1
    int numAlis = alignment.size();
    if (numSortRefs > numSpecies)
	numSortRefs = numSpecies;
    for (int i=0; i< numSortRefs; i++){
	// use the first species (probably often important to users)
	// and otherwise evenly spaced indices to the list of species
	size_t s = i * numSpecies / numSortRefs;
	alignment.sort(SortCriterion(s)); // sort by species number s
	list<Alignment*>::iterator ait, bit;
	for (ait = alignment.begin(); ait != alignment.end();) {
	    bit = ait;
	    ++bit;
	    if (bit == alignment.end())
		break;
	    // ait and bit point to neighbors in this order
	    // all rows present in both alignments must be very close neighbors (at most 3 codons distance)
	    if (mergeable(*ait, *bit, 9, 1.0)){ 
		(*ait)->merge(*bit);
		alignment.erase(bit);
	    } else {
		++ait;
	    }
	}
    }
    cout << "MAF compatification reduced the number of aligments from " << numAlis << " to " <<  alignment.size() << " (to "
	 << setprecision(3) <<  (100.0 * alignment.size() / numAlis) << "%)." << endl;
}

/**
 * 
 */
void GenomicMSA::findGeneRanges(){
    // temporary code, to be done properly later
    int maxIntronLen = 100000;
    int maxGeneLen = 1000000 - 2*GeneMSA::padding;
    int numAlis = alignment.size();
    alignment.sort(SortCriterion(0)); // sort by species number s
    list<Alignment*>::iterator ait, bit;
    for (ait = alignment.begin(); ait != alignment.end();) {
	bit = ait;
	++bit;
	if (bit == alignment.end())
	    break;
	// ait and bit point to neighbors in this order
	// all rows present in both alignments must be very close neighbors (at most 3 codons distance)
	if (((*ait)->maxRange() + (*bit)->maxRange() + maxIntronLen <= maxGeneLen)
	    && mergeable(*ait, *bit, maxIntronLen, 0.6)){ 
	    (*ait)->merge(*bit);
	    alignment.erase(bit);
	} else {
	    ++ait;
	}
    }
    cout << "findGeneRanges reduced the number of aligments from " << numAlis << " to " <<  alignment.size() << " (to "
	 << setprecision(3) <<  (100.0 * alignment.size() / numAlis) << "%)." << endl;
}


// Mario: this apparently assumes that the maf aligments are sorted so that neighboring
// alignments can be merged in this order. TODO: sorting.
// aligned sequences, which have a distance of at most maxGapLen bases will be merged
void GenomicMSA::mergeAlignment(int maxGapLen, float percentSpeciesAligned) {
    list<Alignment*>::iterator it_pos=this->alignment.begin();
    list<Alignment*>::iterator it_prev;
    list<block>::iterator it_block;

    // specieslimit is the percentage of species, which have to be in the alignmentblocks to merge them
    float specieslimit = (*it_pos)->rows.size() - (percentSpeciesAligned * (*it_pos)->rows.size());
    for (it_pos = this->alignment.begin(); it_pos != this->alignment.end(); it_pos++) {
        bool complete=true;
        int distance=-1;
        int emptyBegin=1;
        int geneStart = 0;
        int geneSeqLen = 0;
        Strand geneStrand = STRAND_UNKNOWN;
        string geneChr = "";

        if (it_pos==this->alignment.begin()) {
            it_prev=it_pos;
        } else {
            while (complete && (it_pos!=this->alignment.end())) {
                int count=0;
                bool sameDistance=true;
                for (int j=0; j<(*it_prev)->rows.size(); j++) {
                    if (((*it_prev)->rows.at(j)==NULL) || ((*it_pos)->rows.at(j)==NULL)) {
                        count++;
                        if (count>=specieslimit) {
                            complete=false;
                            break;
                        }
                    }
                    if ((*it_prev)->rows.at(j)!=NULL) {
                        geneChr = (*it_prev)->rows.at(j)->seqID;
                        geneStrand = (*it_prev)->rows.at(j)->strand;
                        geneStart = (*it_prev)->rows.at(j)->start;
                        geneSeqLen = (*it_prev)->rows.at(j)->seqLen;
                    }

                    //  max maxGapLen bases apart, on the same strand and on the same seqID, when there is an alignmentblock before
                    if (((*it_pos)->rows.at(j)!=NULL) && (geneChr != "")) {
                            if (((*it_pos)->rows.at(j)->strand != geneStrand) || ((*it_pos)->rows.at(j)->seqID != geneChr)) {
                                complete=false;
                                break;
                            }
                    }
                    if (((*it_pos)->rows.at(j)!=NULL) && (geneStart != 0)) {
                            if (((*it_pos)->rows.at(j)->start - (geneStart + geneSeqLen) > maxGapLen) || ((*it_pos)->rows.at(j)->start - (geneStart + geneSeqLen) < 0)) {
                                complete=false;
                                break;
                            }
                    }
                    // sequences of the different species have to have the same distance
                    if (((*it_prev)->rows.at(j)!=NULL) && ((*it_pos)->rows.at(j)!=NULL)) {
                       if (distance==-1) {
                            distance=(*it_pos)->rows.at(j)->start - (*it_prev)->rows.at(j)->start - (*it_prev)->rows.at(j)->seqLen;
                        } else
                            if (distance != ((*it_pos)->rows.at(j)->start - (*it_prev)->rows.at(j)->start - (*it_prev)->rows.at(j)->seqLen)) {
                            sameDistance=false;
                        }
                        // ^ operator for logical exclusive OR
                    } else if (((*it_prev)->rows.at(j)!=NULL) ^ ((*it_pos)->rows.at(j)!=NULL)) {
                        sameDistance=false;
                    }
                }

                // the information of two alignment parts become combined
                if (complete) {
                    for (int j=0; j<(*it_prev)->rows.size(); j++) {
                        if (((*it_prev)->rows.at(j)==NULL) && ((*it_pos)->rows.at(j)==NULL)) {
                            emptyBegin++;
                        } else if (((*it_prev)->rows.at(j) == NULL) && ((*it_pos)->rows.at(j) != NULL)) {
                            (*it_prev)->rows.at(j) = (*it_pos)->rows.at(j);
                            for (int k=0; k<emptyBegin; k++) {
                                (*it_prev)->rows.at(j)->cmpStarts.insert((*it_prev)->rows.at(j)->cmpStarts.begin(),NULL);
                            }
                        } else if (((*it_prev)->rows.at(j)!=NULL) && ((*it_pos)->rows.at(j)!=NULL)) {
                            if (!sameDistance) {
                                int *cmpStart_ptr;
                                cmpStart_ptr= new int;
                                *cmpStart_ptr = (*it_pos)->rows.at(j)->start + ((*it_prev)->aliLen - (*it_prev)->rows.at(j)->seqLen);
                                (*it_prev)->rows.at(j)->cmpStarts.push_back(cmpStart_ptr);
                            }
                            for (list<block>::iterator it=(*it_pos)->rows.at(j)->sequence.begin(); it!=(*it_pos)->rows.at(j)->sequence.end(); it++) {
                                it->begin=it->begin+((*it_prev)->aliLen - (*it_prev)->rows.at(j)->seqLen);
                                it->previousGaps=it->previousGaps + ((*it_prev)->aliLen - (*it_prev)->rows.at(j)->seqLen);
                            }
                            (*it_prev)->aliLen = (*it_prev)->aliLen + (*it_pos)->aliLen
                                    + (*it_pos)->rows.at(j)->start - (*it_prev)->rows.at(j)->start - (*it_prev)->rows.at(j)->seqLen;
                            (*it_prev)->rows.at(j)->seqLen = (*it_prev)->rows.at(j)->seqLen + (*it_pos)->rows.at(j)->seqLen
                                    + (*it_pos)->rows.at(j)->start - (*it_prev)->rows.at(j)->start - (*it_prev)->rows.at(j)->seqLen;

                            it_block = (*it_prev)->rows.at(j)->sequence.end();
                            (*it_prev)->rows.at(j)->sequence.splice(it_block,(*it_pos)->rows.at(j)->sequence);
                        } else {
                            (*it_prev)->rows.at(j)->cmpStarts.push_back(NULL);
                        }
                    }
		    //delete (*it_pos);
                    it_pos=this->alignment.erase(it_pos);
                }
            }
        }
        it_prev=it_pos;
    }
}


// pops the first alignment from list
GeneMSA* GenomicMSA::getNextGene() {
    if (alignment.empty())
	return NULL;
    GeneMSA *geneRange = new GeneMSA(rsa, alignment.front());
    alignment.pop_front();
    return geneRange;
}

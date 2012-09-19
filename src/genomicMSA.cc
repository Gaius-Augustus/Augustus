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
 **********************************************************************/

#include "exoncand.hh"
#include "genomicMSA.hh"
#include "orthoexon.hh"
#include <fstream>
#include <iostream>
#include <string>

// merges the merged alignment parts so that a gene is possibly in this segment
GeneMSA* GenomicMSA::getNextGene() {
    GeneMSA *ptr=new GeneMSA();
    int max_intron_length=50000;
    int max_GeneMSA_length = 1000000;
    vector<Strand >geneStrand;
    vector<string> geneChr;
    vector<int> geneStarts, geneSeqLens, seqStart;
    bool geneRange=true;
    if (this->alignment.empty()) {
        delete ptr;
        return NULL;
    }
    list<AlignmentBlock*>::iterator it_pos=this->alignment.begin();
    list<AlignmentBlock*>::iterator it_prev=it_pos;
    it_pos++;
    list<AlignmentBlock*>::iterator it_succ=it_pos;
    if (it_succ != this->alignment.end()) {
        it_succ++;
    }
    for (int i=0; i<(*it_prev)->alignSpeciesTuple.size(); i++) {
        geneStrand.push_back(STRAND_UNKNOWN);
        geneChr.push_back("");
        geneStarts.push_back(0);
        geneSeqLens.push_back(0);
        seqStart.push_back(-1);
    }
    while (geneRange && (it_pos!=this->alignment.end())) {
        for (int i=0; i<(*it_prev)->alignSpeciesTuple.size(); i++) {
            if ((*it_prev)->alignSpeciesTuple.at(i)!=NULL) {
                if ((geneChr[i] == "") && (geneStrand[i] == STRAND_UNKNOWN)) {
                    geneChr[i] = (*it_prev)->alignSpeciesTuple.at(i)->seqID.first;
                    geneStrand[i] = (*it_prev)->alignSpeciesTuple.at(i)->strand;
                }
                geneStarts[i] = (*it_prev)->alignSpeciesTuple.at(i)->start;
                geneSeqLens[i] = (*it_prev)->alignSpeciesTuple.at(i)->seqLen;
                if (seqStart[i]== -1) {
                    seqStart[i] = (*it_prev)->alignSpeciesTuple.at(i)->start;
                }
            }
            // the alignment parts have to be on the same strand, the same seqID and less than "max_intron_length" apart
            // if there is only a short possible exon part of a species different to the strand or chromosome, it will be ignored
            if (((*it_pos)->alignSpeciesTuple.at(i)!=NULL) && (geneStrand[i] != STRAND_UNKNOWN) && (geneChr[i] != "")) {
                if (((*it_pos)->alignSpeciesTuple.at(i)->strand != geneStrand[i]) || ((*it_pos)->alignSpeciesTuple.at(i)->seqID.first != geneChr[i])
                      || ((*it_pos)->alignSpeciesTuple.at(i)->start - (geneStarts[i] + geneSeqLens[i]) < 0)) {
                    if ((*it_succ)->alignSpeciesTuple.at(i) != NULL) {
                        if (((*it_succ)->alignSpeciesTuple.at(i)->strand != geneStrand[i]) || ((*it_succ)->alignSpeciesTuple.at(i)->seqID.first != geneChr[i])
                               || ((*it_succ)->alignSpeciesTuple.at(i)->start - (geneStarts[i] + geneSeqLens[i]) < 0) || (((*it_pos)->alignSpeciesTuple.at(i)->seqLen) > 1000)
                               || ((*it_succ)->alignSpeciesTuple.at(i)->start - (geneStarts[i] +  geneSeqLens[i]) > max_intron_length)
                               || ((*it_succ)->alignSpeciesTuple.at(i)->start + (*it_succ)->alignSpeciesTuple.at(i)->seqLen - seqStart[i]) > max_GeneMSA_length) {
                            geneRange=false;
                            break;
                        } else {
                            if (i != 0) {
                                (*it_pos)->alignSpeciesTuple.at(i) = NULL;
                            } else {
                                geneRange=false;
                                break;
                            }
                        }
                    } else {
                        geneRange=false;
                        break;
                    }
                }
            }
            if ((*it_pos)->alignSpeciesTuple.at(i)!=NULL && (*it_prev)->alignSpeciesTuple.at(i)!=NULL) {
                if (((*it_pos)->alignSpeciesTuple.at(i)->start - (geneStarts[i] + geneSeqLens[i]) > max_intron_length)
                    || ((*it_pos)->alignSpeciesTuple.at(i)->start + (*it_pos)->alignSpeciesTuple.at(i)->seqLen - seqStart[i]) > max_GeneMSA_length) {
                    geneRange=false;
                    break;
                }
            }
        }
        if (geneRange) {
            it_prev=it_pos;
            it_pos++;
            if (it_succ != this->alignment.end()) {
                it_succ++;
            }
        }
    }
    ptr->alignment.assign(this->alignment.begin(),it_pos);
    this->alignment.erase(this->alignment.begin(),it_pos);

    /*for (list<AlignmentBlock*>::iterator it=ptr->alignment.begin(); it!=ptr->alignment.end(); it++) {
		cout<<"ptr"<<endl;
		for (int j=0; j<(*it)->alignSpeciesTuple.size(); j++) {
			if ((*it)->alignSpeciesTuple.at(j)!=NULL) {
				//cout<<"Name: "<< (*it)->alignSpeciesTuple.at(j)->name <<endl;
				cout<<"seqID: "<< (*it)->alignSpeciesTuple.at(j)->seqID.first<<endl;
				//cout<<"seqID length: "<< (*it)->alignSpeciesTuple.at(j)->seqID.second<<endl;
				cout<<"start: "<<(*it)->alignSpeciesTuple.at(j)->start<<endl;
				cout<<"cmp_starts: ";
				for (int i=0; i<(*it)->alignSpeciesTuple.at(j)->cmpStarts.size(); i++) {
					if ((*it)->alignSpeciesTuple.at(j)->cmpStarts[i]!=NULL) {
						cout<<*((*it)->alignSpeciesTuple.at(j)->cmpStarts[i])<<"  ";
					} else {
						cout<<"NULL"<<"  ";
					}
				}
				cout<<endl;
				cout<<"sequencelength: "<<(*it)->alignSpeciesTuple.at(j)->seqLen<<endl;
				cout<<"alignmentlength: "<<(*it)->alignSpeciesTuple.at(j)->alignLen<<endl;
				cout<<"Strand: "<<(*it)->alignSpeciesTuple.at(j)->strand<<endl;
				if ((*it)->alignSpeciesTuple.at(j)->sequence.empty()) {
					cout<<"komisch"<<endl;
				}
				while (!(*it)->alignSpeciesTuple.at(j)->sequence.empty()) {
					cout<< "segment_start "<<(*it)->alignSpeciesTuple.at(j)->sequence.front().begin;
					cout<<"   length "<<(*it)->alignSpeciesTuple.at(j)->sequence.front().length;
					cout<<"   gaps "<<(*it)->alignSpeciesTuple.at(j)->sequence.front().previousGaps<<endl;
					(*it)->alignSpeciesTuple.at(j)->sequence.pop_front();
				}
			} else {
				cout<<"Nullpointer"<<endl;
			}
		}
		cout<<"ptr"<<endl;
		cout<<endl;
	}*/
    return ptr;
}

//  *.maf-file is read and will be saved as vector with pointers on a vector with pointers
void GenomicMSA::readAlignment(vector<string> speciesnames) {
    int index=0;
    string indifferent;
    string completeName;
    char str;
    string seq;
    block cur_block;
    AlignSeq *ptr;
    AlignmentBlock alignBlock;
    AlignmentBlock *complPtr;
    vector<string> notExistingSpecies;
    alignBlock.alignSpeciesTuple=vector<AlignSeq*> (speciesnames.size());

    string alignFilename = Constant::alnfile;
    ifstream Alignmentfile;
    Alignmentfile.open(alignFilename.c_str(), ifstream::in);
    if (!Alignmentfile) {
        cerr << "Could not find the alignment file " << alignFilename << "." << endl;
        throw PropertiesError( "GenomicMSA::readAlignment: Could not open this file!" ); //maybe change  PropertiesError
    }

    while (!Alignmentfile.eof()) {
        Alignmentfile>>indifferent;
        int speciesFound = 0;
        if (indifferent=="s") {
            for (int i=0; i<speciesnames.size(); i++) {
                alignBlock.alignSpeciesTuple[i]=NULL;
            }
            //for (int j=0; j<speciesnames.size(); j++) {
            while (speciesnames.size() > speciesFound) {
                if (indifferent!="s") {
                    break;
                } else {
                    indifferent="not s";
                    ptr =new AlignSeq;
                    // reads the name of the species and the seqID from the first column
                    Alignmentfile>>completeName;
                    for (int i=0; i<completeName.length(); i++) {
                        if ((completeName[i]=='-') || (completeName[i]=='.')) { //real seperator is the point '.' for example hs19.chr21, has to be changed
                            ptr->name=completeName.substr(0,i);
                            ptr->seqID.first=completeName.substr(i+1,completeName.length()- 1);
                            ptr->seqID.first=ptr->seqID.first.erase(ptr->seqID.first.find_first_of("("),ptr->seqID.first.length()-1); // version for vergl_syn testfile
                            break;
                        }
                        if (i==completeName.length()-1) {
                            ptr->name=completeName;
                            ptr->seqID.first="unknown";
                        }
                    }
                    for (int i=0; i<speciesnames.size(); i++) {
                        if (speciesnames[i] == ptr->name) {
                            index=i;
                            break;
                        } else if (i == (speciesnames.size() - 1)) {
                            index = -1;
                            if (notExistingSpecies.empty()) {
                                notExistingSpecies.push_back(ptr->name);
                            } else {
                                bool found = false;
                                for (int k=0; k<notExistingSpecies.size(); k++) {
                                    if (notExistingSpecies[k] == ptr->name) {
                                        found =true;
                                        break;
                                    }
                                }
                                if (!found) {
                                    notExistingSpecies.push_back(ptr->name);
                                    found = false;
                                }
                            }
                        }
                    }
                    Alignmentfile>>ptr->offset;
                    ptr->start=ptr->offset+1;
                    ptr->cmpStarts.push_back(&ptr->start);
                    Alignmentfile>>ptr->seqLen;
                    Alignmentfile>>str;
                    if (str == '+') {
                        ptr->strand=plusstrand;
                    } else if (str == '-') {
                        ptr->strand=minusstrand;
                    } else {
                        ptr->strand=STRAND_UNKNOWN;
                    }
                    Alignmentfile>>ptr->seqID.second;
                    Alignmentfile>>seq;
                    ptr->alignLen=seq.length();
                    Alignmentfile>>indifferent;
                    // reads the aligned sequence
                    cur_block.begin=ptr->start;
                    cur_block.length=1;
                    cur_block.previousGaps=0;
                    for (int i=1; i<=seq.length(); i++) {
                        if ((seq[i-1]!='-')&&(seq[i]=='-')) {
                            ptr->sequence.push_back(cur_block);
                        } else if ((seq[i-1]=='-')&&(seq[i]=='-')) {
                            cur_block.previousGaps++;
                        } else if ((seq[i-1]=='-')&&(seq[i]!='-')) {
                            cur_block.begin=ptr->start+i;
                            cur_block.length=1;
                            cur_block.previousGaps++;
                        } else if ((seq[i-1]!='-')&&(seq[i]!='-')) {
                            if (i==seq.length()) {
                                ptr->sequence.push_back(cur_block);

                            } else {
                                cur_block.length++;
                            }
                        }
                    }
                    if (index != -1) {
                        alignBlock.alignSpeciesTuple[index]=ptr;
                        speciesFound++;
                    }
                }
            }
            complPtr =new AlignmentBlock;
            *complPtr=alignBlock;
            this->alignment.push_back(complPtr);
        }
    }
    Alignmentfile.close();

    if (!(notExistingSpecies.empty())) {
        cout<<endl;
        for (int i=0; i<notExistingSpecies.size(); i++) {
            cout<<"species "<<notExistingSpecies[i]<<" isn't included in the tree file"<<endl;
        }
        cout<<endl;
    }

    /*for (list<AlignmentBlock*>::iterator it=this->alignment.begin(); it!=this->alignment.end(); it++) {
		cout<<"Tuplestart"<<endl;
		for (int j=0; j<(*it)->alignSpeciesTuple.size(); j++) {
			if ((*it)->alignSpeciesTuple.at(j)!=NULL) {
				cout<<"Name: "<< (*it)->alignSpeciesTuple.at(j)->name <<endl;
				cout<<"seqID: "<< (*it)->alignSpeciesTuple.at(j)->seqID.first<<endl;
				//cout<<"seqID length: "<< (*it)->alignSpeciesTuple.at(j)->seqID.second<<endl;
				cout<<"start: "<<(*it)->alignSpeciesTuple.at(j)->start<<endl;
				cout<<"cmp_starts: ";
				for (int i=0; i<(*it)->alignSpeciesTuple.at(j)->cmpStarts.size(); i++) {
					if ((*it)->alignSpeciesTuple.at(j)->cmpStarts[i]!=NULL) {
						cout<<*((*it)->alignSpeciesTuple.at(j)->cmpStarts[i])<<"  ";
					} else {
						cout<<"NULL"<<"  ";
					}
				}
				cout<<endl;
				cout<<"sequencelength: "<<(*it)->alignSpeciesTuple.at(j)->seqLen<<endl;
				cout<<"alignmentlength: "<<(*it)->alignSpeciesTuple.at(j)->alignLen<<endl;
				cout<<"Strand: "<<(*it)->alignSpeciesTuple.at(j)->strand<<endl;
				for (list<block>::iterator itb=(*it)->alignSpeciesTuple.at(j)->sequence.begin(); itb!=(*it)->alignSpeciesTuple.at(j)->sequence.end(); itb++) {
					cout<< "segment_start "<<(*itb).begin;
					cout<<"   length "<<(*itb).length;
					cout<<"   gaps "<<(*itb).previousGaps<<endl;
				}
			} else {
				cout<<"Nullpointer"<<endl;
			}
		}
		cout<<"Tupleend"<<endl;
		cout<<endl;
	}*/
}

//aligned sequences, which have a distance of at most 5 bases will be merged
void GenomicMSA::mergeAlignment(int maxGapLen, float percentSpeciesAligned) {
    list<AlignmentBlock*>::iterator it_pos=this->alignment.begin();
    list<AlignmentBlock*>::iterator it_prev;
    list<block>::iterator it_block;

    // specieslimit is the percentage of species, which have to be in the alignmentblocks to merge them
    float specieslimit = (*it_pos)->alignSpeciesTuple.size() - (percentSpeciesAligned * (*it_pos)->alignSpeciesTuple.size());
    for (it_pos=this->alignment.begin(); it_pos!=this->alignment.end(); it_pos++) {
        bool complete=true;
        //bool sameDistance=true;
        int distance=-1;
        int emptyBegin=1;
        int geneStart = 0;
        int geneSeqLen = 0;
        Strand  geneStrand = STRAND_UNKNOWN;
        string geneChr = "";

        if (it_pos==this->alignment.begin()) {
            it_prev=it_pos;
        } else {
            while (complete && (it_pos!=this->alignment.end())) {
                int count=0;
                bool sameDistance=true;
                for (int j=0; j<(*it_prev)->alignSpeciesTuple.size(); j++) {
                    if (((*it_prev)->alignSpeciesTuple.at(j)==NULL) || ((*it_pos)->alignSpeciesTuple.at(j)==NULL)) {
                        count++;
                        if (count>=specieslimit) {
                            complete=false;
                            break;
                        }
                    }
                    if ((*it_prev)->alignSpeciesTuple.at(j)!=NULL) {
                        geneChr = (*it_prev)->alignSpeciesTuple.at(j)->seqID.first;
                        geneStrand = (*it_prev)->alignSpeciesTuple.at(j)->strand;
                        geneStart = (*it_prev)->alignSpeciesTuple.at(j)->start;
                        geneSeqLen = (*it_prev)->alignSpeciesTuple.at(j)->seqLen;
                    }

                    //  max 5 bases apart, on the same strand and on the same seqID, when there is an alignmentblock before
                    if (((*it_pos)->alignSpeciesTuple.at(j)!=NULL) && (geneChr != "")) {
                            if (((*it_pos)->alignSpeciesTuple.at(j)->strand != geneStrand) || ((*it_pos)->alignSpeciesTuple.at(j)->seqID.first != geneChr)) {
                                complete=false;
                                break;
                            }
                    }
                    if (((*it_pos)->alignSpeciesTuple.at(j)!=NULL) && (geneStart != 0)) {
                            if (((*it_pos)->alignSpeciesTuple.at(j)->start - (geneStart + geneSeqLen) > maxGapLen) || ((*it_pos)->alignSpeciesTuple.at(j)->start - (geneStart + geneSeqLen) < 0)) {
                                complete=false;
                                break;
                            }
                    }
                    // sequences of the different species have to have the same distance
                    if (((*it_prev)->alignSpeciesTuple.at(j)!=NULL) && ((*it_pos)->alignSpeciesTuple.at(j)!=NULL)) {
                        if (distance==-1) {
                            distance=(*it_pos)->alignSpeciesTuple.at(j)->start - (*it_prev)->alignSpeciesTuple.at(j)->start - (*it_prev)->alignSpeciesTuple.at(j)->seqLen;
                        } else if (distance != ((*it_pos)->alignSpeciesTuple.at(j)->start - (*it_prev)->alignSpeciesTuple.at(j)->start - (*it_prev)->alignSpeciesTuple.at(j)->seqLen)) {
                            sameDistance=false;
                        }
                        // ^ operator for logical exclusive OR
                    } else if (((*it_prev)->alignSpeciesTuple.at(j)!=NULL) ^ ((*it_pos)->alignSpeciesTuple.at(j)!=NULL)) {
                        sameDistance=false;
                    }
                }

                // the informations of two alignment parts become combined
                if (complete) {
                    for (int j=0; j<(*it_prev)->alignSpeciesTuple.size(); j++) {
                        if (((*it_prev)->alignSpeciesTuple.at(j)==NULL) && ((*it_pos)->alignSpeciesTuple.at(j)==NULL)) {
                            emptyBegin++;
                        } else if (((*it_prev)->alignSpeciesTuple.at(j) == NULL) && ((*it_pos)->alignSpeciesTuple.at(j) != NULL)) {
                            (*it_prev)->alignSpeciesTuple.at(j) = (*it_pos)->alignSpeciesTuple.at(j);
                            for (int k=0; k<emptyBegin; k++) {
                                (*it_prev)->alignSpeciesTuple.at(j)->cmpStarts.insert((*it_prev)->alignSpeciesTuple.at(j)->cmpStarts.begin(),NULL);
                            }
                        } else if (((*it_prev)->alignSpeciesTuple.at(j)!=NULL) && ((*it_pos)->alignSpeciesTuple.at(j)!=NULL)) {
                            if (!sameDistance) {
                                int *cmpStart_ptr;
                                cmpStart_ptr= new int;
                                *cmpStart_ptr = (*it_pos)->alignSpeciesTuple.at(j)->start + ((*it_prev)->alignSpeciesTuple.at(j)->alignLen - (*it_prev)->alignSpeciesTuple.at(j)->seqLen);
                                (*it_prev)->alignSpeciesTuple.at(j)->cmpStarts.push_back(cmpStart_ptr);
                            }
                            /*if (((*it_prev)->alignSpeciesTuple.at(j)->start + (*it_prev)->alignSpeciesTuple.at(j)->seqLen) == (*it_pos)->alignSpeciesTuple.at(j)->start) {
                                it_block = (*it_prev)->alignSpeciesTuple.at(j)->sequence.end();
                                it_block--;
                                it_block->length =  it_block->length + (*it_pos)->alignSpeciesTuple.at(j)->sequence.begin()->length;
                                (*it_pos)->alignSpeciesTuple.at(j)->sequence.pop_front();
                            }*/
                            for (list<block>::iterator it=(*it_pos)->alignSpeciesTuple.at(j)->sequence.begin(); it!=(*it_pos)->alignSpeciesTuple.at(j)->sequence.end(); it++) {
                                it->begin=it->begin+((*it_prev)->alignSpeciesTuple.at(j)->alignLen - (*it_prev)->alignSpeciesTuple.at(j)->seqLen);
                                // maybe have to change that
                                //if (sameDistance) {
                                it->previousGaps=it->previousGaps + ((*it_prev)->alignSpeciesTuple.at(j)->alignLen - (*it_prev)->alignSpeciesTuple.at(j)->seqLen);
                                //}
                            }
                            (*it_prev)->alignSpeciesTuple.at(j)->alignLen = (*it_prev)->alignSpeciesTuple.at(j)->alignLen + (*it_pos)->alignSpeciesTuple.at(j)->alignLen
                                    + (*it_pos)->alignSpeciesTuple.at(j)->start - (*it_prev)->alignSpeciesTuple.at(j)->start - (*it_prev)->alignSpeciesTuple.at(j)->seqLen;
                            (*it_prev)->alignSpeciesTuple.at(j)->seqLen = (*it_prev)->alignSpeciesTuple.at(j)->seqLen + (*it_pos)->alignSpeciesTuple.at(j)->seqLen
                                    + (*it_pos)->alignSpeciesTuple.at(j)->start - (*it_prev)->alignSpeciesTuple.at(j)->start - (*it_prev)->alignSpeciesTuple.at(j)->seqLen;

                            it_block = (*it_prev)->alignSpeciesTuple.at(j)->sequence.end();
                            (*it_prev)->alignSpeciesTuple.at(j)->sequence.splice(it_block,(*it_pos)->alignSpeciesTuple.at(j)->sequence);
                        } else {
                            (*it_prev)->alignSpeciesTuple.at(j)->cmpStarts.push_back(NULL);
                        }
                    }
                    it_pos=this->alignment.erase(it_pos);
                }
            }
        }
        it_prev=it_pos;
    }

    /*for (list<AlignmentBlock*>::iterator it=this->alignment.begin(); it!=this->alignment.end(); it++) {
        cout<<"Tuplestart"<<endl;
        for (int j=0; j<(*it)->alignSpeciesTuple.size(); j++) {
            if ((*it)->alignSpeciesTuple.at(j)!=NULL) {
                cout<<"Name: "<< (*it)->alignSpeciesTuple.at(j)->name <<endl;
                cout<<"seqID: "<< (*it)->alignSpeciesTuple.at(j)->seqID.first<<endl;
                //cout<<"seqID length: "<< (*it)->alignSpeciesTuple.at(j)->seqID.second<<endl;
                cout<<"start: "<<(*it)->alignSpeciesTuple.at(j)->start<<endl;
                cout<<"cmp_starts: ";
                for (int i=0; i<(*it)->alignSpeciesTuple.at(j)->cmpStarts.size(); i++) {
                    if ((*it)->alignSpeciesTuple.at(j)->cmpStarts[i]!=NULL) {
                        cout<<*((*it)->alignSpeciesTuple.at(j)->cmpStarts[i])<<"  ";
                    } else {
                        cout<<"NULL"<<"  ";
                    }
                }
                cout<<endl;
                cout<<"sequencelength: "<<(*it)->alignSpeciesTuple.at(j)->seqLen<<endl;
                cout<<"alignmentlength: "<<(*it)->alignSpeciesTuple.at(j)->alignLen<<endl;
                cout<<"Strand: "<<(*it)->alignSpeciesTuple.at(j)->strand<<endl;
                for (list<block>::iterator itb=(*it)->alignSpeciesTuple.at(j)->sequence.begin(); itb!=(*it)->alignSpeciesTuple.at(j)->sequence.end(); itb++) {
                    cout<< "segment_start "<<(*itb).begin;
                    cout<<"   length "<<(*itb).length;
                    cout<<"   gaps "<<(*itb).previousGaps<<endl;
                }
            } else {
                cout<<"Nullpointer"<<endl;
            }
        }
        cout<<"Tupleend"<<endl;
        cout<<endl;
    }*/
}

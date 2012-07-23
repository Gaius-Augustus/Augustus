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
    int max_intron_length=12000;
    vector<Strand >geneStrand;
    vector<string> geneChr;
    bool geneRange=true;
    if (this->alignment.empty()) {
        delete ptr;
        return NULL;
    }
    list<AlignmentBlock*>::iterator it_pos=this->alignment.begin();
    list<AlignmentBlock*>::iterator it_prev=it_pos;
    it_pos++;
    for (int i=0; i<(*it_prev)->alignSpeciesTupel.size(); i++) {
        geneStrand.push_back(STRAND_UNKNOWN);
        geneChr.push_back("");
    }
    while (geneRange && (it_pos!=this->alignment.end())) {
        for (int i=0; i<(*it_prev)->alignSpeciesTupel.size(); i++) {
            int geneStart = 0;
            int geneSeqLen = 0;
            if ((*it_prev)->alignSpeciesTupel.at(i)!=NULL) {
                if ((geneChr[i] == "") && (geneStrand[i] == STRAND_UNKNOWN)) {
                    geneChr[i] = (*it_prev)->alignSpeciesTupel.at(i)->chromosome.first;
                    geneStrand[i] = (*it_prev)->alignSpeciesTupel.at(i)->strand;
                }
                geneStart = (*it_prev)->alignSpeciesTupel.at(i)->start;
                geneSeqLen = (*it_prev)->alignSpeciesTupel.at(i)->seqLen;
            }
            // the alignment parts have to be on the same strand, the same chromosome and less than "max_intron_length" apart
            if (((*it_pos)->alignSpeciesTupel.at(i)!=NULL)) {
                if (((*it_pos)->alignSpeciesTupel.at(i)->strand != geneStrand[i]) || ((*it_pos)->alignSpeciesTupel.at(i)->chromosome.first != geneChr[i])) {
                    geneRange=false;
                    break;
                }
            }
            if (((*it_pos)->alignSpeciesTupel.at(i)!=NULL) /*&& (geneStart != 0)*/) {
                if (((*it_pos)->alignSpeciesTupel.at(i)->start - (geneStart + geneSeqLen) > max_intron_length) || ((*it_pos)->alignSpeciesTupel.at(i)->start - (geneStart + geneSeqLen) < 0)) {
                    geneRange=false;
                    break;
                }
            }
        }
        if (geneRange) {
            it_prev=it_pos;
            it_pos++;
        }
    }
    ptr->alignment.assign(this->alignment.begin(),it_pos);
    this->alignment.erase(this->alignment.begin(),it_pos);

    /*for (list<AlignmentBlock*>::iterator it=ptr->alignment.begin(); it!=ptr->alignment.end(); it++) {
		cout<<"ptr"<<endl;
		for (int j=0; j<(*it)->alignSpeciesTupel.size(); j++) {
			if ((*it)->alignSpeciesTupel.at(j)!=NULL) {
				//cout<<"Name: "<< (*it)->alignSpeciesTupel.at(j)->name <<endl;
				cout<<"chromosome: "<< (*it)->alignSpeciesTupel.at(j)->chromosome.first<<endl;
				//cout<<"chromosome length: "<< (*it)->alignSpeciesTupel.at(j)->chromosome.second<<endl;
				cout<<"start: "<<(*it)->alignSpeciesTupel.at(j)->start<<endl;
				cout<<"cmp_starts: ";
				for (int i=0; i<(*it)->alignSpeciesTupel.at(j)->cmpStarts.size(); i++) {
					if ((*it)->alignSpeciesTupel.at(j)->cmpStarts[i]!=NULL) {
						cout<<*((*it)->alignSpeciesTupel.at(j)->cmpStarts[i])<<"  ";
					} else {
						cout<<"NULL"<<"  ";
					}
				}
				cout<<endl;
				cout<<"sequencelength: "<<(*it)->alignSpeciesTupel.at(j)->seqLen<<endl;
				cout<<"alignmentlength: "<<(*it)->alignSpeciesTupel.at(j)->alignLen<<endl;
				cout<<"Strand: "<<(*it)->alignSpeciesTupel.at(j)->strand<<endl;
				if ((*it)->alignSpeciesTupel.at(j)->sequence.empty()) {
					cout<<"komisch"<<endl;
				}
				while (!(*it)->alignSpeciesTupel.at(j)->sequence.empty()) {
					cout<< "segment_start "<<(*it)->alignSpeciesTupel.at(j)->sequence.front().begin;
					cout<<"   length "<<(*it)->alignSpeciesTupel.at(j)->sequence.front().length;
					cout<<"   gaps "<<(*it)->alignSpeciesTupel.at(j)->sequence.front().previousGaps<<endl;
					(*it)->alignSpeciesTupel.at(j)->sequence.pop_front();
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

    //read the number of species
    string speciesfile = Constant::speciesfilenames;
    ifstream SpeciesnameFile;
    SpeciesnameFile.open(speciesfile.c_str(), ifstream::in);
    if (!SpeciesnameFile) {
        cerr << "Could not find the species name file " << speciesfile<< "." << endl;
        throw PropertiesError( "GenomicMSA::readAlignment: Could not open this file!" ); //maybe change  PropertiesError
    }
    int noSpecies=0;
    string text;
    while(getline(SpeciesnameFile,text)) {
        noSpecies++;
    }
    SpeciesnameFile.close();
    alignBlock.alignSpeciesTupel=vector<AlignSeq*> (noSpecies);

    string alignFilename = Constant::alnfile;
    ifstream Alignmentfile;
    Alignmentfile.open(alignFilename.c_str(), ifstream::in);
    if (!Alignmentfile) {
        cerr << "Could not find the alignment file " << alignFilename << "." << endl;
        throw PropertiesError( "GenomicMSA::readAlignment: Could not open this file!" ); //maybe change  PropertiesError
    }

    while (!Alignmentfile.eof()) {
        Alignmentfile>>indifferent;
        if (indifferent=="s") {
            for (int i=0; i<noSpecies; i++) {
                alignBlock.alignSpeciesTupel[i]=NULL;
            }
            for (int j=0; j<noSpecies; j++) {
                if (indifferent!="s") {
                    break;
                } else {
                    indifferent="not s";
                    ptr =new AlignSeq;
                    // reads the name of the species and the chromosome from the first column
                    Alignmentfile>>completeName;
                    for (int i=0; i<completeName.length(); i++) {
                        if (completeName[i]=='-') { //real seperator is the point '.' for example hs19.chr21, has to be changed
                            ptr->name=completeName.substr(0,i);
                            ptr->chromosome.first=completeName.substr(i+1,completeName.length()- 1);
                            ptr->chromosome.first=ptr->chromosome.first.erase(ptr->chromosome.first.find_first_of("("),ptr->chromosome.first.length()-1); // version for vergl_syn testfile
                            break;
                        }
                        if (i==completeName.length()-1) {
                            ptr->name=completeName;
                            ptr->chromosome.first="unknown";
                        }
                    }
                    for (int i=0; i<speciesnames.size(); i++) {
                        if (speciesnames[i]==ptr->name) {
                            index=i;
                            break;
                        } else if (i == (speciesnames.size() - 1)) {
                            cout<<"species "<<ptr->name<<" isn't included in the dataset"<<endl;
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
                    Alignmentfile>>ptr->chromosome.second;
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
                    alignBlock.alignSpeciesTupel[index]=ptr;
                }
            }
            complPtr =new AlignmentBlock;
            *complPtr=alignBlock;
            this->alignment.push_back(complPtr);
        }
    }
    Alignmentfile.close();

    /*for (list<AlignmentBlock*>::iterator it=this->alignment.begin(); it!=this->alignment.end(); it++) {
		cout<<"tupelstart"<<endl;
		for (int j=0; j<(*it)->alignSpeciesTupel.size(); j++) {
			if ((*it)->alignSpeciesTupel.at(j)!=NULL) {
				cout<<"Name: "<< (*it)->alignSpeciesTupel.at(j)->name <<endl;
				cout<<"chromosome: "<< (*it)->alignSpeciesTupel.at(j)->chromosome.first<<endl;
				//cout<<"chromosome length: "<< (*it)->alignSpeciesTupel.at(j)->chromosome.second<<endl;
				cout<<"start: "<<(*it)->alignSpeciesTupel.at(j)->start<<endl;
				cout<<"cmp_starts: ";
				for (int i=0; i<(*it)->alignSpeciesTupel.at(j)->cmpStarts.size(); i++) {
					if ((*it)->alignSpeciesTupel.at(j)->cmpStarts[i]!=NULL) {
						cout<<*((*it)->alignSpeciesTupel.at(j)->cmpStarts[i])<<"  ";
					} else {
						cout<<"NULL"<<"  ";
					}
				}
				cout<<endl;
				cout<<"sequencelength: "<<(*it)->alignSpeciesTupel.at(j)->seqLen<<endl;
				cout<<"alignmentlength: "<<(*it)->alignSpeciesTupel.at(j)->alignLen<<endl;
				cout<<"Strand: "<<(*it)->alignSpeciesTupel.at(j)->strand<<endl;
				for (list<block>::iterator itb=(*it)->alignSpeciesTupel.at(j)->sequence.begin(); itb!=(*it)->alignSpeciesTupel.at(j)->sequence.end(); itb++) {
					cout<< "segment_start "<<(*itb).begin;
					cout<<"   length "<<(*itb).length;
					cout<<"   gaps "<<(*itb).previousGaps<<endl;
				}
			} else {
				cout<<"Nullpointer"<<endl;
			}
		}
		cout<<"tupelend"<<endl;
		cout<<endl;
	}*/
}

//aligned sequences, which have a distance of at most 5 bases will be merged
void GenomicMSA::mergeAlignment(int maxGapLen, float percentSpeciesAligned) {
    list<AlignmentBlock*>::iterator it_pos=this->alignment.begin();
    list<AlignmentBlock*>::iterator it_prev;
    //vector<bool> doesNotExistBefore(false,(*it_pos)->alignSpeciesTupel.size());
    vector<Strand >geneStrand;
    vector<string> geneChr;
    list<block>::iterator it_block;

    // specieslimit is the percentage of species, which have to be in the alignmentblocks to merge them
    float specieslimit = (*it_pos)->alignSpeciesTupel.size() - (percentSpeciesAligned * (*it_pos)->alignSpeciesTupel.size());
    for (int i=0; i<(*it_pos)->alignSpeciesTupel.size(); i++) {
        geneStrand.push_back(STRAND_UNKNOWN);
        geneChr.push_back("");
    }
    for (it_pos=this->alignment.begin(); it_pos!=this->alignment.end(); it_pos++) {
        bool complete=true;
        bool sameDistance=true;
        int distance=-1;
        int emptyBegin=1;
        int geneStart = 0;
        int geneSeqLen = 0;
        if (it_pos==this->alignment.begin()) {
            it_prev=it_pos;
        } else {
            while (complete && (it_pos!=this->alignment.end())) {
                int count=0;
                for (int j=0; j<(*it_prev)->alignSpeciesTupel.size(); j++) {
                    if (((*it_prev)->alignSpeciesTupel.at(j)==NULL) || ((*it_pos)->alignSpeciesTupel.at(j)==NULL)) {
                        count++;
                        //doesNotExistBefore[j]=true;
                        if (count>=specieslimit) {
                            complete=false;
                            break;
                        }
                    }
                    if ((*it_prev)->alignSpeciesTupel.at(j)!=NULL) {
                        if ((geneChr[j] == "") && (geneStrand[j] == STRAND_UNKNOWN)) {
                            geneChr[j] = (*it_prev)->alignSpeciesTupel.at(j)->chromosome.first;
                            geneStrand[j] = (*it_prev)->alignSpeciesTupel.at(j)->strand;
                        }
                        geneStart = (*it_prev)->alignSpeciesTupel.at(j)->start;
                        geneSeqLen = (*it_prev)->alignSpeciesTupel.at(j)->seqLen;
                    }

                    //  max 5 bases apart, on the same strand and on the same chromosome, when there is an alignmentblock before
                    if (((*it_pos)->alignSpeciesTupel.at(j)!=NULL) && (geneChr[j] != "")) {
                        //if (doesNotExistBefore[j]==false) {
                            if (((*it_pos)->alignSpeciesTupel.at(j)->strand != geneStrand[j]) || ((*it_pos)->alignSpeciesTupel.at(j)->chromosome.first != geneChr[j])) {
                                complete=false;
                                break;
                            }
                        //}
                    }
                    if (((*it_pos)->alignSpeciesTupel.at(j)!=NULL) && (geneStart != 0)) {
                        //if (doesNotExistBefore[j]==false) {
                            if (((*it_pos)->alignSpeciesTupel.at(j)->start - (geneStart + geneSeqLen) > maxGapLen) || ((*it_pos)->alignSpeciesTupel.at(j)->start - (geneStart + geneSeqLen) < 0)) {
                                complete=false;
                                break;
                            }
                        /*} else {
                            doesNotExistBefore[j]=false;
                        }*/
                    }
                    // sequences of the different species have to have the same distance
                    if (((*it_prev)->alignSpeciesTupel.at(j)!=NULL) && ((*it_pos)->alignSpeciesTupel.at(j)!=NULL)) {
                        if (distance==-1) {
                            distance=(*it_pos)->alignSpeciesTupel.at(j)->start - (*it_prev)->alignSpeciesTupel.at(j)->start - (*it_prev)->alignSpeciesTupel.at(j)->seqLen;
                        } else if (distance != ((*it_pos)->alignSpeciesTupel.at(j)->start - (*it_prev)->alignSpeciesTupel.at(j)->start - (*it_prev)->alignSpeciesTupel.at(j)->seqLen)) {
                            sameDistance=false;
                        }
                        // ^ operator for logical exclusive OR
                    } else if (((*it_prev)->alignSpeciesTupel.at(j)!=NULL) ^ ((*it_pos)->alignSpeciesTupel.at(j)!=NULL)) {
                        sameDistance=false;
                    }
                }

                // the informations of two alignment parts become combined
                if (complete) {
                    for (int j=0; j<(*it_prev)->alignSpeciesTupel.size(); j++) {
                        if (((*it_prev)->alignSpeciesTupel.at(j)==NULL) && ((*it_pos)->alignSpeciesTupel.at(j)==NULL)) {
                            emptyBegin++;
                        } else if (((*it_prev)->alignSpeciesTupel.at(j) == NULL) && ((*it_pos)->alignSpeciesTupel.at(j) != NULL)) {
                            (*it_prev)->alignSpeciesTupel.at(j) = (*it_pos)->alignSpeciesTupel.at(j);
                            for (int k=0; k<emptyBegin; k++) {
                                (*it_prev)->alignSpeciesTupel.at(j)->cmpStarts.insert((*it_prev)->alignSpeciesTupel.at(j)->cmpStarts.begin(),NULL);
                            }
                        } else if (((*it_prev)->alignSpeciesTupel.at(j)!=NULL) && ((*it_pos)->alignSpeciesTupel.at(j)!=NULL)) {
                            if (!sameDistance) {
                                int *cmpStart_ptr;
                                cmpStart_ptr= new int;
                                *cmpStart_ptr = (*it_pos)->alignSpeciesTupel.at(j)->start + ((*it_prev)->alignSpeciesTupel.at(j)->alignLen - (*it_prev)->alignSpeciesTupel.at(j)->seqLen);
                                (*it_prev)->alignSpeciesTupel.at(j)->cmpStarts.push_back(cmpStart_ptr);
                            }
                            for (list<block>::iterator it=(*it_pos)->alignSpeciesTupel.at(j)->sequence.begin(); it!=(*it_pos)->alignSpeciesTupel.at(j)->sequence.end(); it++) {
                                it->begin=it->begin+((*it_prev)->alignSpeciesTupel.at(j)->alignLen - (*it_prev)->alignSpeciesTupel.at(j)->seqLen);
                                // maybe have to change that
                                //if (sameDistance) {
                                it->previousGaps=it->previousGaps + ((*it_prev)->alignSpeciesTupel.at(j)->alignLen - (*it_prev)->alignSpeciesTupel.at(j)->seqLen);
                                //}
                            }
                            (*it_prev)->alignSpeciesTupel.at(j)->alignLen = (*it_prev)->alignSpeciesTupel.at(j)->alignLen + (*it_pos)->alignSpeciesTupel.at(j)->alignLen
                                    + (*it_pos)->alignSpeciesTupel.at(j)->start - (*it_prev)->alignSpeciesTupel.at(j)->start - (*it_prev)->alignSpeciesTupel.at(j)->seqLen;
                            (*it_prev)->alignSpeciesTupel.at(j)->seqLen = (*it_prev)->alignSpeciesTupel.at(j)->seqLen + (*it_pos)->alignSpeciesTupel.at(j)->seqLen
                                    + (*it_pos)->alignSpeciesTupel.at(j)->start - (*it_prev)->alignSpeciesTupel.at(j)->start - (*it_prev)->alignSpeciesTupel.at(j)->seqLen;
                            it_block = (*it_prev)->alignSpeciesTupel.at(j)->sequence.end();
                            (*it_prev)->alignSpeciesTupel.at(j)->sequence.splice(it_block,(*it_pos)->alignSpeciesTupel.at(j)->sequence);
                        } else {
                            (*it_prev)->alignSpeciesTupel.at(j)->cmpStarts.push_back(NULL);
                        }
                    }
                    it_pos=this->alignment.erase(it_pos);
                }
            }
        }
        it_prev=it_pos;
    }

    /*for (list<AlignmentBlock*>::iterator it=this->alignment.begin(); it!=this->alignment.end(); it++) {
        cout<<"tupelstart"<<endl;
        for (int j=0; j<(*it)->alignSpeciesTupel.size(); j++) {
            if ((*it)->alignSpeciesTupel.at(j)!=NULL) {
                cout<<"Name: "<< (*it)->alignSpeciesTupel.at(j)->name <<endl;
                cout<<"chromosome: "<< (*it)->alignSpeciesTupel.at(j)->chromosome.first<<endl;
                //cout<<"chromosome length: "<< (*it)->alignSpeciesTupel.at(j)->chromosome.second<<endl;
                cout<<"start: "<<(*it)->alignSpeciesTupel.at(j)->start<<endl;
                cout<<"cmp_starts: ";
                for (int i=0; i<(*it)->alignSpeciesTupel.at(j)->cmpStarts.size(); i++) {
                    if ((*it)->alignSpeciesTupel.at(j)->cmpStarts[i]!=NULL) {
                        cout<<*((*it)->alignSpeciesTupel.at(j)->cmpStarts[i])<<"  ";
                    } else {
                        cout<<"NULL"<<"  ";
                    }
                }
                cout<<endl;
                cout<<"sequencelength: "<<(*it)->alignSpeciesTupel.at(j)->seqLen<<endl;
                cout<<"alignmentlength: "<<(*it)->alignSpeciesTupel.at(j)->alignLen<<endl;
                cout<<"Strand: "<<(*it)->alignSpeciesTupel.at(j)->strand<<endl;
                for (list<block>::iterator itb=(*it)->alignSpeciesTupel.at(j)->sequence.begin(); itb!=(*it)->alignSpeciesTupel.at(j)->sequence.end(); itb++) {
                    cout<< "segment_start "<<(*itb).begin;
                    cout<<"   length "<<(*itb).length;
                    cout<<"   gaps "<<(*itb).previousGaps<<endl;
                }
            } else {
                cout<<"Nullpointer"<<endl;
            }
        }
        cout<<"tupelend"<<endl;
        cout<<endl;
    }*/
}

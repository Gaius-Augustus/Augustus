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
#include "geneMSA.hh"
#include "orthoexon.hh"
#include "intronmodel.hh"
#include "namgene.hh"
#include <fstream>
#include <iostream>
#include <string>

ofstream *GeneMSA::exonCands_outfile=NULL;
ofstream *GeneMSA::orthoExons_outfile=NULL;

int GeneMSA::getStart(int speciesIdx) {
    if (this->alignment.front()->alignSpeciesTupel.at(speciesIdx)!=NULL) {
        return this->alignment.front()->alignSpeciesTupel.at(speciesIdx)->offset;
    } else {
        return -1;
    }
}

int GeneMSA::getEnd(int speciesIdx){
    if (this->alignment.front()->alignSpeciesTupel.at(speciesIdx)!=NULL) {
        return this->alignment.back()->alignSpeciesTupel.at(speciesIdx)->offset + this->alignment.back()->alignSpeciesTupel.at(speciesIdx)->seqLen - 1;
    } else {
        return 0;
    }
}

string GeneMSA::getChr(int speciesIdx) {
    if (this->alignment.front()->alignSpeciesTupel.at(speciesIdx)!=NULL) {
        return this->alignment.front()->alignSpeciesTupel.at(speciesIdx)->chromosome.first;
    } else {
        return "";
    }
}

Strand GeneMSA::getStrand(int speciesIdx) {
    if (this->alignment.front()->alignSpeciesTupel.at(speciesIdx)!=NULL) {
        return this->alignment.front()->alignSpeciesTupel.at(speciesIdx)->strand;
    } else {
        return STRAND_UNKNOWN;
    }
}

list<ExonCandidate*>* GeneMSA::getExonCands(int speciesIdx) {
    if (this->exoncands.at(speciesIdx)!=NULL) {
        return exoncands.at(speciesIdx);
    } else {
        return NULL;
    }
}

list<OrthoExon> GeneMSA::getOrthoExons(){
    if (!this->orthoExonsList.empty()) {
        return this->orthoExonsList;
    } else {
        cout<<" no orthologue exons in this part of the alignment found"<<endl;
        return this->orthoExonsList;  // list is empty
    }
}

// use this compare function if you want to use ExonCandidate::computeIntervals
bool compType (ExonCandidate* a, ExonCandidate* b) {
    if (a->type != b->type) {
        return (a->type<b->type);
    } else  if (a->end != b->end) {
        return (a->end<b->end);
    } else {
        return (a->begin<b->begin);
    }
}

// compare function to find the correct reference value to compute the aligned position
bool compCmpStarts (int a, int *b) {
    if (b!=NULL) {
        return (a < *b);
    } else {
        return false;
    }
}

// compare function to find the correct aligned position of a base
bool compAlignedPos (int a, block b) {
    return (a < b.begin - b.previousGaps);
}

// compare function to find the correct real position of a base
bool compRealPos (int a, block b) {
    return (a < b.begin);
}

// compare function to sort the exon candidates by startposition
bool compBegin (ExonCandidate* a, ExonCandidate* b) {
    if (a->begin != b->begin) {
        return (a->begin<b->begin);
    } else  if (a->end != b->end) {
        return (a->end<b->end);
    } else {
        return (a->type<b->type);
    }
}


// adds the keys to the map function
map<string,ExonCandidate*>* GeneMSA::addToHash(list<ExonCandidate*> *ec) {
    map<string, ExonCandidate*> *hashCandidates =new  map<string, ExonCandidate*>;;
    for (list<ExonCandidate*>::iterator lit=ec->begin(); lit!=ec->end(); lit++) {
        (*hashCandidates)[(*lit)->createKey()] = *lit;
    }
    return hashCandidates;
}

// computes the exon candidates of a dna sequence
void GeneMSA::createExonCands(const char *dna, float assqthresh, float dssqthresh, float motifqthresh){
    int n = strlen(dna);
    int max_exon_length = 12000;
    int frame;
    Double p;
    list<int> exonStart;
    list<int> exonRCStop;
    list<int> exonASS;
    list<int> exonRDSS;
    ExonCandidate *ptr;
    list<ExonCandidate*> *candidates=new list<ExonCandidate*>;
    Double assminprob = IntronModel::assBinProbs.getMinProb(assqthresh) * IntronModel::getAssMotifProbThreshold(motifqthresh);
    Double dssminprob = IntronModel::dssBinProbs.getMinProb(dssqthresh);

    OpenReadingFrame orf(dna, max_exon_length, n);
    for (int i=0; i<=n - 1; i++) {
        // positions of all startcodons "atg"
        if (onStart(dna+i)) {
            exonStart.push_back(i + 1);
        }
        // positons of all ASSs "ag"
        if (onASS(dna+i) && (i + Constant::ass_whole_size() - Constant::ass_start < n)) {
            p = IntronModel::aSSProb(i - Constant::ass_upwindow_size - Constant::ass_start, true);
            if (p >= assminprob ) {
                exonASS.push_back(i);
            }
        }
        // positions of all reverse DSS "ac"
        if (onRDSS(dna+i) && (i + Constant::dss_whole_size() - Constant::dss_end < n)) {
            p = IntronModel::dSSProb(i - Constant::dss_end, false);
            if (p >= dssminprob) {
                exonRDSS.push_back(i);
            }
        }
        // positions of all reverse complementary stopcodons "cta, tta, tca"
        if (onRCStopcodon(dna+i)) {
            exonRCStop.push_back(i);
        }
    }
    list<int>::reverse_iterator ritStart = exonStart.rbegin();
    list<int>::reverse_iterator ritStart_cur=ritStart;
    list<int>::reverse_iterator ritASS = exonASS.rbegin();
    list<int>::reverse_iterator ritASS_cur=ritASS;
    list<int>::reverse_iterator ritRDSS = exonRDSS.rbegin();
    list<int>::reverse_iterator ritRDSS_cur=ritRDSS;
    list<int>::reverse_iterator ritRCStop = exonRCStop.rbegin();
    list<int>::reverse_iterator ritRCStop_cur=ritRCStop;

    for (int i= n - 1; i>= 2; i--) {
        // computing single genes on the forward strand with at least startcodon+codon+stopcodon
        if (ochre(dna+i) || opal(dna+i) ||amber(dna+i)) {
            ritStart=ritStart_cur;
            while ((i<*ritStart)&&(ritStart!=exonStart.rend())){
                ritStart++;
            }
            ritStart_cur=ritStart;
            int lmb = orf.leftmostExonBegin(0,i,true);
            while ((lmb<=*ritStart)&&(i-*ritStart<=max_exon_length)&&(ritStart!=exonStart.rend())) {
                if ((i-*ritStart>=4)&&((i-*ritStart+1)%3==0)) {
                    ptr = new ExonCandidate;
                    ptr->begin=*ritStart - 1;
                    ptr->end=i + 2;
                    ptr->type=singleGene;
                    candidates->push_back(ptr);
                }
                ritStart++;
            };
        }

        // computing initials on the forward strand with at least startcodon plus base
        if (onDSS(dna + i) && (i + Constant::dss_whole_size() - Constant::dss_start  < n)) {
            p = IntronModel::dSSProb(i - Constant::dss_start,true);
            for (frame=0; frame<=2; frame++) {
                ritStart=ritStart_cur;
                while ((i<*ritStart)&&(ritStart!=exonStart.rend())){
                    ritStart++;
                }
                ritStart_cur=ritStart;
                int lmb = orf.leftmostExonBegin(frame,i,true);
                while((lmb<=*ritStart)&&(i-*ritStart<=max_exon_length)&&(ritStart!=exonStart.rend())) {
                    if ((i-*ritStart>=3)&&((i-*ritStart+1)%3==frame) && (p >= dssminprob)) {
                        ptr = new ExonCandidate;
                        ptr->begin=*ritStart - 1;
                        ptr->end=i - 1;
                        if (frame==0) {
                            ptr->type=initial_0;
                        } else if (frame==1) {
                            ptr->type=initial_1;
                        } else {
                            ptr->type=initial_2;
                        }
                        candidates->push_back(ptr);
                    }
                    ritStart++;
                };
            }

            // computing internals on the forward strand with at least one codon
            for (frame=0; frame<=2; frame++) {
                ritASS=ritASS_cur;
                while ((i<*ritASS)&&(ritASS!=exonASS.rend())){
                    ritASS++;
                }
                ritASS_cur=ritASS;
                int lmb = orf.leftmostExonBegin(frame,i,true);
                while((lmb<=*ritASS)&&(i-*ritASS<=max_exon_length)&&(ritASS!=exonASS.rend())) {
                    if ((i-*ritASS>=5) && (p >= dssminprob)) {
                        ptr = new ExonCandidate;
                        ptr->begin=*ritASS + 2;
                        ptr->end=i - 1;
                        if (frame==0) {
                            ptr->type=internal_0;
                        } else if (frame==1) {
                            ptr->type=internal_1;
                        } else {
                            ptr->type=internal_2;
                        }
                        candidates->push_back(ptr);
                    }
                    ritASS++;
                };
            }
        }

        // computing terminals on the forward strand with at least one base stopcodon
        if (ochre(dna+i) || opal(dna+i) ||amber(dna+i)) {
            for (frame=0; frame<=2; frame++) {
                ritASS=ritASS_cur;
                while ((i<*ritASS)&&(ritASS!=exonASS.rend())){
                    ritASS++;
                }
                ritASS_cur=ritASS;
                while ((i-*ritASS<=max_exon_length)&&(ritASS!=exonASS.rend())) {
                    if ((i-*ritASS>=3)&&((i-*ritASS+1)%3==frame)&&(*ritASS>=orf.leftmostExonBegin(0,i,true))) {
                        ptr = new ExonCandidate;
                        ptr->begin=*ritASS + 2;
                        ptr->end=i + 2;
                        ptr->type=terminal_exon;
                        candidates->push_back(ptr);
                    }
                    ritASS++;
                };
            }
        }

        // computing single genes on the reverse strand with at least startcodon+codon+ctopcodon
        if (onRStart(dna+i)) {
            ritRCStop=ritRCStop_cur;
            while ((i<*ritRCStop)&&(ritRCStop!=exonRCStop.rend())){
                ritRCStop++;
            }
            ritRCStop_cur=ritRCStop;
            while ((i-*ritRCStop<=max_exon_length)&&(ritRCStop!=exonRCStop.rend())) {
                if ((i-*ritRCStop>=6)&&((i-*ritRCStop)%3==0)) {
                    ptr = new ExonCandidate;
                    ptr->begin=*ritRCStop;
                    ptr->end=i + 2;
                    ptr->type=rsingleGene;
                    candidates->push_back(ptr);
                    break;
                } else {
                    ritRCStop++;
                }
            };
        }

        // computing initials on the reverse strand with at least startcodon plus base
        if (onRStart(dna+i)) {
            for (frame=0; frame<=2; frame++) {
                ritRDSS=ritRDSS_cur;
                while ((i<*ritRDSS)&&(ritRDSS!=exonRDSS.rend())){
                    ritRDSS++;
                }
                ritRDSS_cur=ritRDSS;
                int lmb = orf.leftmostExonBegin(2,i,false);
                while((lmb<=*ritRDSS+2)&&(i-*ritRDSS<=max_exon_length)&&(ritRDSS!=exonRDSS.rend())) {
                    if ((i-*ritRDSS>=2)&&((i+1-*ritRDSS)%3==frame)) {
                        ptr = new ExonCandidate;
                        ptr->begin=*ritRDSS + 2;
                        ptr->end=i + 2;
                        ptr->type=rinitial_exon;
                        candidates->push_back(ptr);
                    }
                    ritRDSS++;
                };
            }
        }

        // computing internals on the reverse strand with at least a codon
        if (onRASS(dna+i) && (i + Constant::ass_upwindow_size + Constant::ass_whole_size() - Constant::ass_start < n)) {
            p = IntronModel::aSSProb(i-Constant::ass_end, false);
            for (frame=0; frame<=2; frame++) {
                ritRDSS=ritRDSS_cur;
                while ((i<*ritRDSS)&&(ritRDSS!=exonRDSS.rend())){
                    ritRDSS++;
                }
                ritRDSS_cur=ritRDSS;
                int lmb = orf.leftmostExonBegin(frame,i,false);
                while((lmb<=*ritRDSS)&&(i-*ritRDSS<=max_exon_length)&&(ritRDSS!=exonRDSS.rend())) {
                    if (i-*ritRDSS>=5 && (p >= assminprob)) {
                        ptr = new ExonCandidate;
                        ptr->begin=*ritRDSS + 2;
                        ptr->end=i - 1;
                        if (frame==0) {
                            ptr->type=rinternal_0;
                        } else if (frame==1) {
                            ptr->type=rinternal_1;
                        } else {
                            ptr->type=rinternal_2;
                        }
                        candidates->push_back(ptr);
                    }
                    ritRDSS++;
                };
            }
        }

        // computing terminals on the reverse strand with at least one base plus stopcodon
        if (onRASS(dna+i) && (i + Constant::ass_upwindow_size + Constant::ass_whole_size() - Constant::ass_start < n)) {
            p = IntronModel::aSSProb(i-Constant::ass_end, false);
            for (frame=0; frame<=2; frame++) {
                ritRCStop=ritRCStop_cur;
                while ((i<*ritRCStop)&&(ritRCStop!=exonRCStop.rend())){
                    ritRCStop++;
                }
                ritRCStop_cur=ritRCStop;
                while ((i-*ritRCStop<=max_exon_length)&&(ritRCStop!=exonRCStop.rend())) {
                    if ((i-*ritRCStop>=4)&&((i-*ritRCStop)%3==frame) && (p >= assminprob)) {
                        ptr = new ExonCandidate;
                        ptr->begin=*ritRCStop;
                        ptr->end=i - 1;
                        if (frame==0) {
                            ptr->type=rterminal_0;
                        } else if (frame==1) {
                            ptr->type=rterminal_1;
                        } else {
                            ptr->type=rterminal_2;
                        }
                        candidates->push_back(ptr);
                        break;
                    } else {
                        ritRCStop++;
                    }
                };
            }
        }
    }

    candidates->sort(compBegin);
    this->exoncands.push_back(candidates);
    existingCandidates.push_back(addToHash(candidates));

   /* if (candidates!=NULL) {
        cout<<"possible exons: "<<endl;
        for (list<ExonCandidate*>::iterator lit=candidates->begin(); lit!=candidates->end(); lit++) {
            cout<<"start: "<<(*lit)->begin;
            cout<<"     end: "<<(*lit)->end;
            cout<<"     typ "<<(*lit)->type<<endl;
        }
    }*/
    cout << "thresholds: dssminprob=" << dssminprob << " assminprob=" << assminprob << endl;
}

// computes the aligned position of a base in an alignment and the 'block' where the base is found
pair <int, int> GeneMSA::getAlignedPosition(AlignSeq *as_ptr, int pos) {
    list<block>::iterator it, cit;
    vector<int*>::iterator it_cmpStart;
    pair <int, int> alignedPos;
    it=upper_bound(as_ptr->sequence.begin(), as_ptr->sequence.end(), pos, compAlignedPos);
    it--;
    it_cmpStart=upper_bound(as_ptr->cmpStarts.begin(), as_ptr->cmpStarts.end(), (pos + it->previousGaps), compCmpStarts);
    it_cmpStart--;
    while ((*it_cmpStart)==NULL) {
        it_cmpStart--;
    }
    vector<int*>::difference_type idx = distance(as_ptr->cmpStarts.begin(), it_cmpStart);
    alignedPos.second=idx;
    alignedPos.first=pos + it->previousGaps - (*(*it_cmpStart));
    return alignedPos;
}

// computes the real position of a base dependent on its position in the alignment
int GeneMSA::getRealPosition(AlignSeq *ptr, int pos, int idx) {
    list<block>::iterator it;
    int realPos, alignedPos;
    if (ptr->cmpStarts[idx]!=NULL) {
        alignedPos = *ptr->cmpStarts[idx] + pos;
    } else {
        return -1;
    }
    it=upper_bound(ptr->sequence.begin(), ptr->sequence.end(), alignedPos, compRealPos);
    it--;
    realPos=alignedPos - it->previousGaps;
    return realPos;
}

// searches for the orthologue exons of the exon candidates of the reference species
void GeneMSA::createOrthoExons(vector<int> offsets) {
    list<AlignmentBlock*>::iterator it_ab;
    list<ExonCandidate*> cands = *getExonCands(0);
    AlignSeq *as_ptr, *ptr;
    OrthoExon oe;
    vector<ExonCandidate*> *ortho_ptr;
    string key;
    bool found, hasOrthoExon;

    it_ab=this->alignment.begin();
    for (list<ExonCandidate*>::iterator it=cands.begin(); it!=cands.end(); it++) {
        found=false;
        hasOrthoExon=false;

        while ((!found) && (!this->alignment.empty())) {
            if ((*it_ab)->alignSpeciesTupel.at(0)->start > (*it)->begin + offsets[0] + 1) {
                cout<<" exon "<<(*it)->begin + offsets[0] + 1<<".."<<(*it)->end + offsets[0] + 1<<" is outside (in front of) the aligned range"<<endl;
                found=true;
            } else if (((*it_ab)->alignSpeciesTupel.at(0)->start <= (*it)->begin + offsets[0] + 1)
                         && ((*it_ab)->alignSpeciesTupel.at(0)->start + (*it_ab)->alignSpeciesTupel.at(0)->seqLen - 1 >= (*it)->end + offsets[0] + 1)) {

                as_ptr=(*it_ab)->alignSpeciesTupel.at(0);
                if (as_ptr != NULL) {
                    ortho_ptr = new vector<ExonCandidate*>;
                    ortho_ptr->push_back(*it);
                    for (int i=1; i<(*it_ab)->alignSpeciesTupel.size(); i++) {
                        int alignedPosStart = getAlignedPosition(as_ptr, (*it)->begin + offsets[0]).first;
                        int idxStart = getAlignedPosition(as_ptr, (*it)->begin + offsets[0]).second;
                        int alignedPosEnd = getAlignedPosition(as_ptr, (*it)->end + offsets[0]).first;
                        int idxEnd = getAlignedPosition(as_ptr, (*it)->end + offsets[0]).second;

                        ptr = (*it_ab)->alignSpeciesTupel.at(i);
                        (*it)->complementType();  // ????????
                        if (ptr != NULL) {
                            int realStart = getRealPosition(ptr, alignedPosStart, idxStart);
                            int realEnd = getRealPosition(ptr, alignedPosEnd, idxEnd);
                            if (realStart!=-1 && realEnd!=-1) {
                                if (((*it)->type > -1) && ( ptr->strand == plusstrand)) {
                                    key = (itoa(realStart - offsets[i]) + ":" + itoa(realEnd - offsets[i]) + ":" + itoa((*it)->type));
                                } else if (((*it)->type > -1) && (ptr->strand == minusstrand)) {
                                    key = (itoa(realStart - offsets[i]) + ":" + itoa(realEnd - offsets[i]) + ":" + itoa((*it)->complementType()));
                                } else {
                                    key = "No key";
                                }
                                map<string, ExonCandidate*>::iterator map_it = (*existingCandidates[i]).find(key);
                                if (map_it!=(*existingCandidates[i]).end()) {
                                    ortho_ptr->push_back((*existingCandidates[i]).find(key)->second);
                                    hasOrthoExon=true;
                                } else {
                                    ortho_ptr->push_back(NULL);
                                }
                            } else {
                                ortho_ptr->push_back(NULL);
                            }
                        } else {
                            ortho_ptr->push_back(NULL);
                        }
                    }
                    oe.orthoex  = *ortho_ptr;
                    if (hasOrthoExon) {
                        this->orthoExonsList.push_back(oe);
                    }
                }
                found=true;
            } else {
                it_ab++;
                if ((*it_ab)->alignSpeciesTupel.at(0)->start <= (*it)->begin + offsets[0] + 1) {
                    this->alignment.pop_front();
                } else {
                    cout<<" exon "<<(*it)->begin + offsets[0] + 1<<".."<<(*it)->end + offsets[0] + 1<<" hasn't start and end in a coherently aligned sequence"<<endl;
                    it_ab--;
                    found=true;
                }
                if (this->alignment.empty()) {
                    cout<<" exon "<<(*it)->begin + offsets[0] + 1<<".."<<(*it)->end + offsets[0] + 1<<" is outside (behind) the aligned range"<<endl;
                }
            }
        }
    }
    /*cout<<endl;
    for (list<OrthoExon>::iterator it_oe=orthoExonsList.begin(); it_oe!=orthoExonsList.end(); it_oe++) {
        for (int j=0; j<(*it_ab)->alignSpeciesTupel.size(); j++) {
            if (it_oe->orthoex.at(j)!=NULL && (*it_ab)->alignSpeciesTupel.at(j)!=NULL) {
                cout<<"species: "<<(*it_ab)->alignSpeciesTupel.at(j)->name<<"  ";
                cout<<"start: "<<it_oe->orthoex.at(j)->begin;
                cout<<"     end: "<<it_oe->orthoex.at(j)->end;
                cout<<"     typ "<<it_oe->orthoex.at(j)->type<<endl;
            } else {
                cout<<"species ??? has no orthologue exon "<<endl;
            }
        }
        cout<<endl;
    }*/
}

 void GeneMSA::openOutputFiles(){
    string outputdirectory;  //directory for output files

    try {
        outputdirectory = Properties::getProperty("/examples/outdirectory");
    } catch (...) {
        outputdirectory = "";
    }
    outputdirectory = expandHome(outputdirectory); //replace "~" by "$HOME"

    string file_exoncand = outputdirectory + "exonCands" + ".gff";
    ofstream *os_ec = new ofstream(file_exoncand.c_str());
    if (os_ec) {
        exonCands_outfile=os_ec;
        (*os_ec) << PREAMBLE << endl;
        (*os_ec) << "#\n#-----  exon candidates  -----" << endl << "#" << endl;
    }
    string file_orthoexon = outputdirectory + "orthoExons" + ".gff";
    ofstream *os_oe = new ofstream(file_orthoexon.c_str());
    if (os_oe) {
        orthoExons_outfile=os_oe;
        (*os_oe) << PREAMBLE << endl;
        (*os_oe) << "#\n#----- orthologue exons  -----" << endl << "#" << endl;
    }
}

//writes the exon candidates of a species of a dna segment in the file 'exonCands.gff'
void GeneMSA::printExonCands(vector<int> offsets) {
    streambuf *console = cout.rdbuf();  //save old buf
    cout.rdbuf(exonCands_outfile->rdbuf());  //redirect cout to 'exonCands.gff'!

    if (!(this->exoncands.empty())) {
        for (int i=0; i<this->exoncands.size(); i++) {
            list<AlignmentBlock*>::iterator it_exists = this->alignment.end();
            for (list<AlignmentBlock*>::iterator it=this->alignment.begin(); it!=this->alignment.end(); it++) {
                if ((*it)->alignSpeciesTupel.at(i)!=NULL) {
                    cout<<" # speciesname: "<<(*it)->alignSpeciesTupel.at(i)->name <<endl; //ToDo: maybe change that, using the vector with the speciesnames instead a loop
                    it_exists = it;
                }
                break;
            }
            if ((this->exoncands.at(i)!=NULL) && (it_exists != this->alignment.end())) {
                for (list<ExonCandidate*>::iterator it_exonCands=this->exoncands.at(i)->begin(); it_exonCands!=this->exoncands.at(i)->end(); it_exonCands++) {
                    cout<<(*it_exists)->alignSpeciesTupel.at(i)->chromosome.first<< "\tAUGUSTUS\t";
                    cout<<stateExonTypeIdentifiers[(*it_exonCands)->type]<<"\t"<<(*it_exonCands)->begin + offsets[i]+1<<"\t"<<(*it_exonCands)->end + offsets[i]+1<<"\t"<<(*it_exonCands)->score<<"\t";
                    cout<<(((*it_exists)->alignSpeciesTupel.at(i)->strand == plusstrand)? '+' : '-')<<"\t"<<exonTypeReadingFrames[(*it_exonCands)->type]<<"\t";
                    cout<<"group ID"<<endl;
                }
            } else {
                cout<<"#  no exon candidates found " << endl;
            }
        }
    } else {
        cout << "#  no exon candidates found at all" << endl;
    }
    cout.rdbuf(console); //reset to standard output again
}

//writes the orthologue exons of the different species in the file 'orthoExons.gff'
void GeneMSA::printOrthoExons(vector<int> offsets) {
    streambuf *console = cout.rdbuf();
    cout.rdbuf(orthoExons_outfile->rdbuf());  //direct cout to 'orthoExons.gff'!

    if (!(this->orthoExonsList.empty())) {
        int orthoExonNo = 1;
        for (list<OrthoExon>::iterator it_oe=orthoExonsList.begin(); it_oe!=orthoExonsList.end(); it_oe++) {
            for (int j=0; j<it_oe->orthoex.size(); j++) {
                list<AlignmentBlock*>::iterator it_exists = this->alignment.end();
                if (it_oe->orthoex.at(j)!=NULL) {
                    for (list<AlignmentBlock*>::iterator it=this->alignment.begin(); it!=this->alignment.end(); it++) {
                        if ((*it)->alignSpeciesTupel.at(j)!=NULL) {
                            it_exists = it;
                        }
                        break;
                    }
                    if (it_exists != this->alignment.end()) {
                        cout<<(*it_exists)->alignSpeciesTupel.at(j)->chromosome.first<< "\tAUGUSTUS\t";
                        cout<<stateExonTypeIdentifiers[it_oe->orthoex.at(j)->type]<<"\t"<<it_oe->orthoex.at(j)->begin + offsets[j]+1<<"\t"<<it_oe->orthoex.at(j)->end + offsets[j]+1<<"\t"<<it_oe->orthoex.at(j)->score<<"\t";
                        cout<<(((*it_exists)->alignSpeciesTupel.at(j)->strand == plusstrand)? '+' : '-')<<"\t"<<exonTypeReadingFrames[it_oe->orthoex.at(j)->type]<<"\t";
                        cout<<"ID="<< orthoExonNo <<" ; "<<"species "<<"\""<<(*it_exists)->alignSpeciesTupel.at(j)->name<<"\""<<endl;
                    }
                } /*else {
                     cout<<"#  species ??? no orthologue exon found " << endl;
                }*/
            }
            cout<<"# "<<endl;
            orthoExonNo++;
        }
    } else {
        cout << "#  no orthologue exons found at all" << endl;
    }
    cout.rdbuf(console); //reset to standard output again
}

void GeneMSA::closeOutputFiles(){
    if (exonCands_outfile) {
        if(exonCands_outfile->is_open()) {
            exonCands_outfile->close();
            delete exonCands_outfile;
        }
    }
    if (orthoExons_outfile) {
        if(orthoExons_outfile->is_open()) {
            orthoExons_outfile->close();
            delete orthoExons_outfile;
        }
    }
}

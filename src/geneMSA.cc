/**********************************************************************
 * file:    geneMSA.cc
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

using namespace std;

PhyloTree *GeneMSA::tree = NULL;
int GeneMSA::utr_range = 1000;
int GeneMSA::orthoExonID = 1;
int GeneMSA::geneRangeID = 1;
//vector<int> GeneMSA::exonCandID;
vector<ofstream*> GeneMSA::exonCands_outfiles;
vector<ofstream*> GeneMSA::orthoExons_outfiles;
vector<ofstream*> GeneMSA::geneRanges_outfiles;

string GeneMSA::getName(int speciesIdx) {
    for (list<AlignmentBlock*>::iterator it=this->alignment.begin(); it!=this->alignment.end(); it++) {
        if ((*it)->alignSpeciesTuple.at(speciesIdx)!=NULL) {
            return (*it)->alignSpeciesTuple.at(speciesIdx)->name;
        }
    }
    return "";
}

long int GeneMSA::getSeqIDLength(int speciesIdx) {
    for (list<AlignmentBlock*>::iterator it=this->alignment.begin(); it!=this->alignment.end(); it++) {
        if ((*it)->alignSpeciesTuple.at(speciesIdx)!=NULL) {
            return (*it)->alignSpeciesTuple.at(speciesIdx)->seqID.second;
        }
    }
    return 0;
}

string GeneMSA::getSeqID(int speciesIdx) {
    for (list<AlignmentBlock*>::iterator it=this->alignment.begin(); it!=this->alignment.end(); it++) {
        if ((*it)->alignSpeciesTuple.at(speciesIdx)!=NULL) {
            return (*it)->alignSpeciesTuple.at(speciesIdx)->seqID.first;
        }
    }
    return "";
}

Strand GeneMSA::getStrand(int speciesIdx) {
    for (list<AlignmentBlock*>::iterator it=this->alignment.begin(); it!=this->alignment.end(); it++) {
        if ((*it)->alignSpeciesTuple.at(speciesIdx)!=NULL) {
            return (*it)->alignSpeciesTuple.at(speciesIdx)->strand;
        }
    }
    return STRAND_UNKNOWN;
}

int GeneMSA::getStart(int speciesIdx) {
    for (list<AlignmentBlock*>::iterator it=this->alignment.begin(); it!=this->alignment.end(); it++) {
        if ((*it)->alignSpeciesTuple.at(speciesIdx)!=NULL) {
            if (this->getStrand(speciesIdx)==plusstrand) {
                return (*it)->alignSpeciesTuple.at(speciesIdx)->offset - utr_range;
            } else {
                for (list<AlignmentBlock*>::reverse_iterator rit=this->alignment.rbegin(); rit!=this->alignment.rend(); rit++) {
                    if ((*rit)->alignSpeciesTuple.at(speciesIdx)!=NULL) {
                        return ((*it)->alignSpeciesTuple.at(speciesIdx)->seqID.second) - ((*rit)->alignSpeciesTuple.at(speciesIdx)->offset)
                                - ((*rit)->alignSpeciesTuple.at(speciesIdx)->seqLen) - utr_range;
                    }
                }
            }
        }
    }
    return -1;
}

int GeneMSA::getEnd(int speciesIdx){
    for (list<AlignmentBlock*>::reverse_iterator rit=this->alignment.rbegin(); rit!=this->alignment.rend(); rit++) {
        if ((*rit)->alignSpeciesTuple.at(speciesIdx)!=NULL) {
            if (this->getStrand(speciesIdx)==plusstrand) {
                return (*rit)->alignSpeciesTuple.at(speciesIdx)->offset + (*rit)->alignSpeciesTuple.at(speciesIdx)->seqLen - 1 + utr_range;
            } else {
                for (list<AlignmentBlock*>::iterator it=this->alignment.begin(); it!=this->alignment.end(); it++) {
                    if ((*it)->alignSpeciesTuple.at(speciesIdx)!=NULL) {
                        return ((*rit)->alignSpeciesTuple.at(speciesIdx)->seqID.second) - ((*it)->alignSpeciesTuple.at(speciesIdx)->offset + 1) + utr_range;
                    }
                }
            }
        }
    }
    return 0;
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
bool compRealPos (int a, block b) {
    return (a < b.begin);
}

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
    map<string, ExonCandidate*> *hashCandidates =new  map<string, ExonCandidate*>;
    if (!(ec->empty())) {
        for (list<ExonCandidate*>::iterator lit=ec->begin(); lit!=ec->end(); lit++) {
            (*hashCandidates)[(*lit)->createKey()] = *lit;
        }
        return hashCandidates;
    } else {
        return NULL;
    }
}

// computes the exon candidates of a dna sequence
void GeneMSA::createExonCands(const char *dna, double assmotifqthresh, double assqthresh, double dssqthresh){
    int n = strlen(dna);
    int max_exon_length = 12000;
    int frame;
    Double p;
    list<int> exonStart;
    list<int> exonRCStop;
    list< pair<int, Double> > exonASS;
    list< pair<int, Double> > exonRDSS;
    ExonCandidate *ec;
    list<ExonCandidate*> *candidates=new list<ExonCandidate*>;
    Properties::assignProperty("/CompPred/assmotifqthresh", assmotifqthresh);
    Properties::assignProperty("/CompPred/assqthresh", assqthresh);
    Properties::assignProperty("/CompPred/dssqthresh", dssqthresh);
    Double assminprob = IntronModel::assBinProbs.getMinProb(assqthresh) * IntronModel::getAssMotifProbThreshold(assmotifqthresh);
    Double dssminprob = IntronModel::dssBinProbs.getMinProb(dssqthresh);

    OpenReadingFrame orf(dna, max_exon_length, n);
    // preprocessing all left coordinates of an exon candidate interval
    for (int i=0; i<=n - 1; i++) {
        pair<int, Double> ss_score;
        // positions of all startcodons "atg"
        if (onStart(dna+i)) {
            exonStart.push_back(i + 1);
        }
        // positons of all ASSs "ag"
        if (onASS(dna+i) && (i + Constant::ass_whole_size() - Constant::ass_start < n)) {
            p = IntronModel::aSSProb(i - Constant::ass_upwindow_size - Constant::ass_start, true);
            if (p >= assminprob ) {
                ss_score.first = i;
                ss_score.second = p;
                exonASS.push_back(ss_score);
            }
        }
        // positions of all reverse DSS "ac"
        if (onRDSS(dna+i) && (i + Constant::dss_whole_size() - Constant::dss_end < n)) {
            p = IntronModel::dSSProb(i - Constant::dss_end, false);
            if (p >= dssminprob) {
                ss_score.first = i;
                ss_score.second = p;
                exonRDSS.push_back(ss_score);
            }
        }
        // positions of all reverse complementary stopcodons "cta, tta, tca"
        if (onRCStopcodon(dna+i)) {
            exonRCStop.push_back(i);
        }
    }
    list<int>::reverse_iterator ritStart = exonStart.rbegin();
    list<int>::reverse_iterator ritStart_cur = ritStart;
    list< pair<int, Double> >::reverse_iterator ritASS = exonASS.rbegin();
    list< pair<int, Double> >::reverse_iterator ritASS_cur = ritASS;
    list< pair<int, Double> >::reverse_iterator ritRDSS = exonRDSS.rbegin();
    list< pair<int, Double> >::reverse_iterator ritRDSS_cur = ritRDSS;
    list<int>::reverse_iterator ritRCStop = exonRCStop.rbegin();
    list<int>::reverse_iterator ritRCStop_cur = ritRCStop;

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
                if ((i-*ritStart >= 201) && ((i-*ritStart+1)%3==0)) {
                    ec = new ExonCandidate;
                    ec->begin=*ritStart - 1;
                    ec->end=i + 2;
                    ec->type=singleGene;
                    candidates->push_back(ec);
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
                        ec = new ExonCandidate;
                        ec->begin = *ritStart - 1;
                        ec->end = i - 1;
                        ec->dssScore = p;
                        if (frame == 0) {
                            ec->type = initial_0;
                        } else if (frame == 1) {
                            ec->type = initial_1;
                        } else {
                            ec->type = initial_2;
                        }
                        candidates->push_back(ec);
                    }
                    ritStart++;
                };
            }

            // computing internals on the forward strand with at least one codon
            for (frame=0; frame<=2; frame++) {
                ritASS=ritASS_cur;
                while ((i<(*ritASS).first)&&(ritASS!=exonASS.rend())){
                    ritASS++;
                }
                ritASS_cur = ritASS;
                int lmb = orf.leftmostExonBegin(frame,i,true);
                while(lmb <= (*ritASS).first && i-(*ritASS).first <= max_exon_length && ritASS != exonASS.rend()) {
                    if ((i-(*ritASS).first>=5) && (p >= dssminprob)) {
                        ec = new ExonCandidate;
                        ec->begin = (*ritASS).first + 2;
                        ec->end = i - 1;
                        ec->assScore = (*ritASS).second;
                        ec->dssScore = p;
                        if (frame == 0) {
                            ec->type = internal_0;
                        } else if (frame==1) {
                            ec->type = internal_1;
                        } else {
                            ec->type = internal_2;
                        }
                        candidates->push_back(ec);
                    }
                    ritASS++;
                };
            }
        }

        // computing terminals on the forward strand with at least one base stopcodon
        if (ochre(dna+i) || opal(dna+i) ||amber(dna+i)) {
            for (frame=0; frame<=2; frame++) {
                ritASS=ritASS_cur;
                while ((i<(*ritASS).first)&&(ritASS!=exonASS.rend())){
                    ritASS++;
                }
                ritASS_cur=ritASS;
                while ((i-(*ritASS).first <= max_exon_length)&&(ritASS!=exonASS.rend())) {
                    if ((i-(*ritASS).first>=3)&&((i-(*ritASS).first + 1)%3==frame)&&((*ritASS).first>=orf.leftmostExonBegin(0,i,true))) {
                        ec = new ExonCandidate;
                        ec->begin = (*ritASS).first + 2;
                        ec->end = i + 2;
                        ec->assScore = (*ritASS).second;
                        ec->type = terminal_exon;
                        candidates->push_back(ec);
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
            while ((i-*ritRCStop <= max_exon_length)&&(ritRCStop!=exonRCStop.rend())) {
                if ((i-*ritRCStop >= 201)&&((i-*ritRCStop)%3==0)) {
                    ec = new ExonCandidate;
                    ec->begin=*ritRCStop;
                    ec->end=i + 2;
                    ec->type=rsingleGene;
                    candidates->push_back(ec);
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
                while ((i<(*ritRDSS).first)&&(ritRDSS!=exonRDSS.rend())){
                    ritRDSS++;
                }
                ritRDSS_cur=ritRDSS;
                int lmb = orf.leftmostExonBegin(2,i,false);
                while((lmb<=(*ritRDSS).first+2)&&(i-(*ritRDSS).first<=max_exon_length)&&(ritRDSS!=exonRDSS.rend())) {
                    if ((i-(*ritRDSS).first>=2)&&((i+1-(*ritRDSS).first)%3==frame)) {
                        ec = new ExonCandidate;
                        ec->begin=(*ritRDSS).first + 2;
                        ec->end=i + 2;
                        ec->dssScore = (*ritRDSS).second;
                        ec->type=rinitial_exon;
                        candidates->push_back(ec);
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
                while ((i<(*ritRDSS).first)&&(ritRDSS!=exonRDSS.rend())){
                    ritRDSS++;
                }
                ritRDSS_cur=ritRDSS;
                int lmb = orf.leftmostExonBegin(frame,i,false);
                while((lmb<=(*ritRDSS).first)&&(i-(*ritRDSS).first<=max_exon_length)&&(ritRDSS!=exonRDSS.rend())) {
                    if (i-(*ritRDSS).first>=5 && (p >= assminprob)) {
                        ec = new ExonCandidate;
                        ec->begin=(*ritRDSS).first + 2;
                        ec->end=i - 1;
                        ec->dssScore = (*ritRDSS).second;
                        ec->assScore = p;
                        if (frame==0) {
                            ec->type=rinternal_0;
                        } else if (frame==1) {
                            ec->type=rinternal_1;
                        } else {
                            ec->type=rinternal_2;
                        }
                        candidates->push_back(ec);
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
                        ec = new ExonCandidate;
                        ec->begin=*ritRCStop;
                        ec->end=i - 1;
                        ec->assScore = p;
                        if (frame==0) {
                            ec->type=rterminal_2;
                        } else if (frame==1) {
                            ec->type=rterminal_1;
                        } else {
                            ec->type=rterminal_0;
                        }
                        candidates->push_back(ec);
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

    /*if (candidates!=NULL) {
        cout<<"possible exons: "<<endl;
        for (list<ExonCandidate*>::iterator lit=candidates->begin(); lit!=candidates->end(); lit++) {
            cout<<"start: "<<(*lit)->begin;
            cout<<"     end: "<<(*lit)->end;
            cout<<"  ass_score: "<<(*lit)->assScore;
            cout<<"  dss_score: "<<(*lit)->dssScore;
            cout<<"     typ "<<(*lit)->type<<endl;
        }
    }*/
    //cout << "thresholds: dssminprob=" << dssminprob << " assminprob=" << assminprob << endl;
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
    string key;
    bool found, hasOrthoExon;

    it_ab=this->alignment.begin();
    for (list<ExonCandidate*>::iterator ec=cands.begin(); ec!=cands.end(); ec++) {
        found=false;
        hasOrthoExon=false;
        if (!oe.orthoex.empty()) {
           oe.orthoex.clear();
           oe.orthoex.resize(this->exoncands.size());
        }
        while ((!found) && (it_ab!=this->alignment.end())) {

            if ((*it_ab)->alignSpeciesTuple.at(0)->start > (*ec)->begin + offsets[0] + 1) {
                //cout<<" exon "<<(*ec)->begin + offsets[0] + 1<<".."<<(*ec)->end + offsets[0] + 1<<" is outside (in front of) the aligned range"<<endl;
                found=true;
            } else if (((*it_ab)->alignSpeciesTuple.at(0)->start <= (*ec)->begin + offsets[0] + 1)
                         && ((*it_ab)->alignSpeciesTuple.at(0)->start + (*it_ab)->alignSpeciesTuple.at(0)->seqLen - 1 >= (*ec)->end + offsets[0] + 1)) {
                as_ptr=(*it_ab)->alignSpeciesTuple.at(0);
                if (as_ptr != NULL) {
                    oe.orthoex[0]=(*ec);
                    for (int i=1; i<(*it_ab)->alignSpeciesTuple.size(); i++) {
                        int alignedPosStart = getAlignedPosition(as_ptr, (*ec)->begin + offsets[0]).first;
                        int idxStart = getAlignedPosition(as_ptr, (*ec)->begin + offsets[0]).second;
                        int alignedPosEnd = getAlignedPosition(as_ptr, (*ec)->end + offsets[0]).first;
                        int idxEnd = getAlignedPosition(as_ptr, (*ec)->end + offsets[0]).second;

                        ptr = (*it_ab)->alignSpeciesTuple.at(i);
                        if (ptr != NULL) {
                            int realStart = getRealPosition(ptr, alignedPosStart, idxStart);
                            //cout<<"Start "<<realStart<<endl;
                            int realEnd = getRealPosition(ptr, alignedPosEnd, idxEnd);
                            if ((realStart != -1) && (realEnd != -1) && ((mod3((*ec)->end - (*ec)->begin)) == (mod3(realEnd - realStart)))) {
                                if ((*ec)->type > -1) {
                                    key = (itoa(realStart - offsets[i]) + ":" + itoa(realEnd - offsets[i]) + ":" + itoa((*ec)->type));
                                } else {
                                    key = "no key";
                                }
                                if (existingCandidates[i]!=NULL && existingCandidates[0]!=NULL) {
                                    map<string, ExonCandidate*>::iterator map_it = (*existingCandidates[i]).find(key);
                                    if (map_it!=(*existingCandidates[i]).end()) {
                                        oe.orthoex[i]=(map_it->second);
                                        hasOrthoExon=true;
                                    }
                                }
                            }
                        }
                    }
                    if (hasOrthoExon) {
                        this->orthoExonsList.push_back(oe);
                    }
                }
                found=true;
            } else {
                it_ab++;
                if ((*it_ab)->alignSpeciesTuple.at(0)->start > (*ec)->begin + offsets[0] + 1) {
                    //cout<<" exon "<<(*ec)->begin + offsets[0] + 1<<".."<<(*ec)->end + offsets[0] + 1<<" hasn't start and end in a coherently aligned sequence"<<endl;
                    it_ab--;
                    found=true;
                }
                if (this->getEnd(0) - utr_range < (*ec)->begin + offsets[0] + 1) {
                    cout<<" exon "<<(*ec)->begin + offsets[0] + 1<<".."<<(*ec)->end + offsets[0] + 1<<" is outside (behind) the aligned range"<<endl;
                }
            }
        }
    }
    //cout<<endl;
    /*for (list<OrthoExon>::iterator it_oe=orthoExonsList.begin(); it_oe!=orthoExonsList.end(); it_oe++) {
        for (int j=0; j<this->exoncands.size(); j++) {
            if (it_oe->orthoex.at(j)!=NULL && (*it_ab)->alignSpeciesTuple.at(j)!=NULL) {
                cout<<"species: "<<this->getName(j)<<"  ";
                cout<<"start: "<<it_oe->orthoex.at(j)->begin;
                cout<<"     end: "<<it_oe->orthoex.at(j)->end;
                cout<<"     typ "<<it_oe->orthoex.at(j)->type<<endl;
            } else {
                cout<<"species "<<this->getName(j)<<" has no orthologue exon "<<endl;
            }
        }
        cout<<endl;
    }*/
}

 void GeneMSA::openOutputFiles(){
    string outputdirectory;  //directory for output files
    try {
        outputdirectory = Properties::getProperty("/CompPred/outdir_orthoexons");
    } catch (...) {
        outputdirectory = "";
    }
    outputdirectory = expandHome(outputdirectory); //replace "~" by "$HOME"
    exonCands_outfiles.resize(tree->species.size());
    orthoExons_outfiles.resize(tree->species.size());
    geneRanges_outfiles.resize(tree->species.size());
    for (int i=0; i<tree->species.size(); i++) {
        string file_exoncand = outputdirectory + "exonCands." + tree->species[i] + ".gff3";
        ofstream *os_ec = new ofstream(file_exoncand.c_str());
        if (os_ec!=NULL) {
            exonCands_outfiles[i]=os_ec;
            (*os_ec) << PREAMBLE << endl;
            (*os_ec) << "#\n#-----  exon candidates  -----" << endl << "#" << endl;
        }
        string file_geneRanges = outputdirectory + "geneRanges." + tree->species[i] + ".gff3";
        ofstream *os_gr = new ofstream(file_geneRanges.c_str());
        if (os_gr!=NULL) {
            geneRanges_outfiles[i]=os_gr;
            (*os_gr) << PREAMBLE << endl;
            (*os_gr) << "#\n#-----  possible gene ranges  -----" << endl << "#" << endl;
        }

        string file_orthoexon = outputdirectory + "orthoExons." + tree->species[i] + ".gff3";
        ofstream *os_oe = new ofstream(file_orthoexon.c_str());
        if (os_oe) {
            orthoExons_outfiles[i]=os_oe;
            (*os_oe) << PREAMBLE << endl;
            (*os_oe) << "#\n#----- orthologue exons  -----" << endl << "#" << endl;
        }
    }
}

 //writes the possible gene ranges of a species in the file 'geneRanges.speciesnames.gff'
 void GeneMSA::printGeneRanges() {
     streambuf *console = cout.rdbuf();  //save old buf
     if (!(this->exoncands.empty())) {
         for (int i=0; i<this->exoncands.size(); i++) {
             cout.rdbuf(geneRanges_outfiles[i]->rdbuf());  //redirect cout to 'geneRanges.speciesname.gff'
             if (this->exoncands.at(i)!=NULL) {
                 cout << this->getSeqID(i)<<"\tGeneRange\t"<<"exon\t"<<this->getStart(i) + 1<<"\t"<< this->getEnd(i) + 1<<"\t0\t";
                 if (this->getStrand(i) == plusstrand) {
                     cout<<'+'<<"\t";
                 } else {
                     cout<<'-'<<"\t";
                 }
                 cout<<".\t"<<"Name="<<geneRangeID<<endl;
             }
         }
         geneRangeID++;
     }
     cout.rdbuf(console); //reset to standard output again
 }

//writes the exon candidates of a species of a dna segment in the file 'exonCands.species.gff'
void GeneMSA::printExonCands(vector<int> offsets) {
    //exonCandID.resize(this->exoncands.size());
    /*for (int j=0; j<exonCandID.size(); j++) {
        exonCandID[j]=1;
    }*/
    streambuf *console = cout.rdbuf();  //save old buf
    if (!(this->exoncands.empty())) {
        for (int i=0; i<this->exoncands.size(); i++) {
            cout.rdbuf(exonCands_outfiles[i]->rdbuf());  //redirect cout to 'exonCands.speciesname.gff'
            if (this->exoncands.at(i)!=NULL) {
                cout << "# sequence:\t" << this->getName(i)<<"\t"<<this->getStart(i)<< "-"<< this->getEnd(i)<<"  "<<this->getEnd(i) - this->getStart(i)<<"bp"<<endl;
                for (list<ExonCandidate*>::iterator it_exonCands=this->exoncands.at(i)->begin(); it_exonCands!=this->exoncands.at(i)->end(); it_exonCands++) {
                    cout<<this->getSeqID(i)<< "\tEC\t"<<"exon\t";
                    if (this->getStrand(i) == plusstrand) {
                        cout<<(*it_exonCands)->begin + offsets[i]+1<<"\t"<<(*it_exonCands)->end + offsets[i]+1<<"\t"<<(*it_exonCands)->score<<"\t";
                        cout<<'+'<<"\t";
                        //cout<<mod3(3-(exonTypeReadingFrames[(*it_exonCands)->type] - ((*it_exonCands)->end - (*it_exonCands)->begin + 1)))<<"\t";
                    } else {
                        cout<<this->getSeqIDLength(i) - ((*it_exonCands)->end + offsets[i])<<"\t"<<this->getSeqIDLength(i) - ((*it_exonCands)->begin+ offsets[i]) <<"\t";
                        cout<<(*it_exonCands)->score<<"\t"<<'-'<<"\t";
                        //cout<<mod3(2-(exonTypeReadingFrames[(*it_exonCands)->type]))<<"\t";
                    }
                    cout<<mod3(3-(exonTypeReadingFrames[(*it_exonCands)->type] - ((*it_exonCands)->end - (*it_exonCands)->begin + 1)))<<"\t";
                    cout<<"Name="<<stateExonTypeIdentifiers[(*it_exonCands)->type]<<endl;
                    //exonCandID[i]++;
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

//writes the orthologue exons of the different species in the files 'orthoExons.species.gff3'
void GeneMSA::printOrthoExons(vector<int> offsets) {
    streambuf *console = cout.rdbuf();
    if (!(this->orthoExonsList.empty())) {
        for (list<OrthoExon>::iterator it_oe=orthoExonsList.begin(); it_oe!=orthoExonsList.end(); it_oe++) {
            printSingleOrthoExon(*it_oe, offsets);
            orthoExonID++;
        }
    } else {
        //cout << "#  no orthologue exons found at all" << endl;
    }
    cout.rdbuf(console); //reset to standard output again
}

//writes the orthologue exons of the different species in the files 'orthoExons.species.gff3'
void GeneMSA::printSingleOrthoExon(OrthoExon &oe, vector<int> offsets) {
    for (int j=0; j<oe.orthoex.size(); j++) {
        cout.rdbuf(orthoExons_outfiles[j]->rdbuf());  //direct cout to 'orthoExons.speciesname.gff3'
        if (oe.orthoex.at(j)!=NULL) {
            cout<<this->getSeqID(j)<< "\tOE1\t"<<"exon"<<"\t";
            if (this->getStrand(j) == plusstrand) {
                cout <<oe.orthoex.at(j)->begin + offsets[j]+1<<"\t"<<oe.orthoex.at(j)->end + offsets[j]+1<<"\t"<<oe.orthoex.at(j)->score<<"\t";
                cout<<'+'<<"\t";
            } else {
                cout <<this->getSeqIDLength(j) - (oe.orthoex.at(j)->end + offsets[j])<<"\t"<<this->getSeqIDLength(j) - (oe.orthoex.at(j)->begin + offsets[j])<<"\t";
                cout<<oe.orthoex.at(j)->score<<"\t"<<'-'<<"\t";
            }
            cout<<mod3(3-(exonTypeReadingFrames[oe.orthoex.at(j)->type] - (oe.orthoex.at(j)->end - oe.orthoex.at(j)->begin + 1))) <<"\t";
            cout<<"ID="<<orthoExonID<<";Name="<<orthoExonID<<";Note="<<stateExonTypeIdentifiers[oe.orthoex.at(j)->type]<<endl;
        }
    }
}


void GeneMSA::closeOutputFiles(){
    for (int i=0; i<tree->species.size(); i++) {
        if (exonCands_outfiles[i]) {
            if(exonCands_outfiles[i]->is_open()) {
                exonCands_outfiles[i]->close();
                delete exonCands_outfiles[i];
            }
        }
        if (geneRanges_outfiles[i]) {
            if(geneRanges_outfiles[i]->is_open()) {
                geneRanges_outfiles[i]->close();
                delete geneRanges_outfiles[i];
            }
        }
        if (orthoExons_outfiles[i]) {
            if(orthoExons_outfiles[i]->is_open()) {
                orthoExons_outfiles[i]->close();
                delete orthoExons_outfiles[i];
            }
        }
    }
}

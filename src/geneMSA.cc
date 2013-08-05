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
 * 02.08.13| Mario Stanke       | rewrite of the merging of alignment blocks
 **********************************************************************/

#include "geneticcode.hh"
#include "exoncand.hh"
#include "genomicMSA.hh"
#include "geneMSA.hh"
#include "orthoexon.hh"
#include "intronmodel.hh"
#include "namgene.hh"
#include "orthograph.hh"
#include <fstream>
#include <iostream>
#include <string>
#include <sys/time.h>

using namespace std;

PhyloTree *GeneMSA::tree = NULL;
CodonEvo *GeneMSA::codonevo = NULL;
int GeneMSA::padding = 1000;
int GeneMSA::orthoExonID = 1;
int GeneMSA::geneRangeID = 1;
vector<int> GeneMSA::exonCandID;
vector<ofstream*> GeneMSA::exonCands_outfiles;
vector<ofstream*> GeneMSA::orthoExons_outfiles;
vector<ofstream*> GeneMSA::geneRanges_outfiles;
vector<ofstream*> GeneMSA::omega_outfiles;
ofstream *GeneMSA::pamlFile;

/*
 * constructor of GeneMSA
 */
GeneMSA::GeneMSA(RandSeqAccess *rsa, Alignment *a) {
    int maxDNAPieceSize = Properties::getIntProperty( "maxDNAPieceSize" );

    this->rsa = rsa;
    if (!alignment)
	throw ProjectError("Internal error in GeneMSA: Alignment missing.");
    if (alignment->numRows != rsa->getNumSpecies())
	throw ProjectError("Error in GeneMSA: Number of species in alignment is not matching the on in the tree.");
    /** construct the gene ranges, th
     * now: simple copy. TODO: extend region when apparently part of the alignment is missing8
     * human   ***********---*******************
     * mouse   *******-------******************-
     * chicken ***********----------------------
     *                        ^
     *                        | extend range here
     */
    starts.resize(alignment->numRows, -1);
    ends.resize(alignment->numRows, -1);
    for (size_t s=0; s < alignment->numRows; s++){
	if (alignment->rows[s]){
	    int chrLen = rsa->getChrLen(s, alignment->rows[s]->seqID);
	    switch (alignment->rows[s]->strand)
		{
		case plusstrand:
		    starts[s] = alignment->rows[s]->chrStart() - padding;
		    ends[s] = alignment->rows[s]->chrEnd() + padding;
		    break;
		case minusstrand:
		    starts[s] = chrLen - alignment->rows[s]->chrEnd() - padding;
		    ends[s] = chrLen - alignment->rows[s]->chrStart() + padding;
		    break;
		default:
		    throw ProjectError("GeneMSA: Unknown strand in alignment.");
		}
	    // ensure that gene predictions on each region will be done IN ONE PIECE
	    if (ends[s]-starts[s]+1 > maxDNAPieceSize){
		int tooMuch = ends[s] - starts[s] + 1 - maxDNAPieceSize;
		starts[s] += (tooMuch + 1)/2;
		ends[s] -= (tooMuch + 1)/2;
	    }
	    if (starts[s] < 0)
		starts[s] = 0;
	    if (ends[s] < 0)
		ends[s] = 0;
	    if (ends[s] >= chrLen)
		ends[s] = chrLen-1;
	    if (starts[s] >= chrLen)
		starts[s] = chrLen-1;
	}
    }
}

/*
 * destructor of GeneMSA
 */
GeneMSA::~GeneMSA(){
    if (alignment)
	delete alignment;
    for (int i=0; i<exoncands.size(); i++) {
	if (exoncands[i]!=NULL) {
	    for(list<ExonCandidate*>::iterator it = exoncands[i]->begin(); it != exoncands[i]->end(); it++){
		delete *it;
	    }
	    exoncands.at(i)->clear();
	    delete exoncands[i];
	}
    }
    for (int i=0; i<existingCandidates.size(); i++) {
	if (existingCandidates[i]!=NULL) {
	    existingCandidates.at(i)->clear();
	    delete existingCandidates[i];
	}
    }
}

int GeneMSA::getGFF3FrameForExon(ExonCandidate *ec) {
    if (isPlusExon(ec->type))
        return mod3(3-(exonTypeReadingFrames[ec->type] - (ec->end - ec->begin + 1)));
    else
        return mod3(2-exonTypeReadingFrames[ec->type]);
}



list<OrthoExon> GeneMSA::getOrthoExons(){
    if (!this->orthoExonsList.empty()) {
        return this->orthoExonsList;
    } else {
        cout<<" no orthologue exons in this part of the alignment found"<<endl;
        return this->orthoExonsList;  // list is empty
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

// compare function to sort the exon candidates by start position
bool compBegin (ExonCandidate* a, ExonCandidate* b) {
    if (a->begin != b->begin) {
        return (a->begin<b->begin);
    } else  if (a->end != b->end) {
        return (a->end<b->end);
    } else {
        return (a->type<b->type);
    }
}

string GeneMSA::reverseString(string text) {
    int i = 0;
    int n = text.length();
    while (i < (n/2)) {
        char c;
        c = (text[i]);
        text[i] = text[n-i-1];
        text[n-i-1] = c;
        i++;
    }
    return text;
}

// adds the keys to the map function
map<string,ExonCandidate*>* GeneMSA::getECHash(list<ExonCandidate*> *ec) {
    map<string, ExonCandidate*> *hashCandidates = new  map<string, ExonCandidate*>;
    if (ec->empty())
	return NULL;
	    
    for (list<ExonCandidate*>::iterator lit=ec->begin(); lit!=ec->end(); lit++)
	(*hashCandidates)[(*lit)->key()] = *lit;

    return hashCandidates;
}

// computes the score for the splice sites of an exon candidate
Double GeneMSA::computeSpliceSiteScore(Double exonScore, Double minProb, Double maxProb) {
    Double score = 0;
    score = (log(exonScore/minProb)/log(maxProb/minProb));
    return score;
}

// computes the exon candidates in a DNA sequence
// TODO: move this to exoncands.cc (Why is this here??)
void GeneMSA::createExonCands(int s, const char *dna, double assmotifqthresh, double assqthresh, double dssqthresh){
    int n = strlen(dna);
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
    Double assmaxprob = IntronModel::assBinProbs.getMinProb(0.99) * IntronModel::getAssMotifProbThreshold(0.9999);
    Double dssmaxprob = IntronModel::dssBinProbs.getMinProb(0.99);

    OpenReadingFrame orf(dna, Constant::max_exon_len, n);
    // preprocessing all left coordinates of an exon candidate interval
    for (int i=0; i<=n - 1; i++) {
        pair<int, Double> ssWithScore;
        // positions of all startcodons "atg"
        if (onStart(dna+i)) {
            exonStart.push_back(i + 1);
        }
        // positons of all ASSs "ag"
        if (onASS(dna+i) && (i + Constant::ass_whole_size() - Constant::ass_start < n)) {
            p = IntronModel::aSSProb(i - Constant::ass_upwindow_size - Constant::ass_start, true);
            if (p >= assminprob ) {
                ssWithScore.first = i;
                ssWithScore.second = computeSpliceSiteScore(p, assminprob, assmaxprob);
                exonASS.push_back(ssWithScore);
            }
        }
        // positions of all reverse DSS "ac"
        if (onRDSS(dna+i) && (i + Constant::dss_whole_size() - Constant::dss_end < n)) {
            p = IntronModel::dSSProb(i - Constant::dss_end, false);
            if (p >= dssminprob) {
                ssWithScore.first = i;
                ssWithScore.second = computeSpliceSiteScore(p, dssminprob, dssmaxprob);
                exonRDSS.push_back(ssWithScore);
            }
        }
        // positions of all reverse complementary stop codons, usually "cta, tta, tca"
        if (GeneticCode::isRCStopcodon(dna+i)) {
            exonRCStop.push_back(i);
        }
    }
    list<int>::reverse_iterator ritStart = exonStart.rbegin();
    list<int>::reverse_iterator ritStart_cur = ritStart;
    list< pair<int, Double> >::reverse_iterator ritASS = exonASS.rbegin();
    list< pair<int, Double> >::reverse_iterator ritASS_cur = ritASS;
    list< pair<int, Double> >::reverse_iterator ritRDSS = exonRDSS.rbegin();
    list< pair<int, Double> >::reverse_iterator ritRDSS_cur = ritRDSS; // remember iterator over iterations of the loop of the right end
    list<int>::reverse_iterator ritRCStop = exonRCStop.rbegin();
    list<int>::reverse_iterator ritRCStop_cur = ritRCStop;

    for (int i = n-1; i >= 2; i--) {
        // computing single genes on the forward strand with at least Constant::min_coding_len
        if (GeneticCode::isStopcodon(dna+i)) {
            ritStart = ritStart_cur;
            while ((i < *ritStart) && (ritStart != exonStart.rend())){
                ritStart++;
            }
            ritStart_cur = ritStart;
            int lmb = orf.leftmostExonBegin(0, i, true);
            while ((lmb <= *ritStart) && (i - *ritStart <= Constant::max_exon_len) && (ritStart != exonStart.rend())) {
                if ((i - *ritStart >= Constant::min_coding_len) && ((i - *ritStart+1)%3 == 0)) {
                    ec = new ExonCandidate;
                    ec->begin = *ritStart - 1;
                    ec->end = i + 2;
                    ec->type = singleGene;
                    candidates->push_back(ec);
                }
                ritStart++;
            };
        }

        // computing initial exons on the forward strand with at least startcodon plus base
        if (onDSS(dna + i) && (i + Constant::dss_whole_size() - Constant::dss_start  < n)) {
            p = IntronModel::dSSProb(i - Constant::dss_start,true);
            for (frame=0; frame<=2; frame++) {
                ritStart=ritStart_cur;
                while ((i<*ritStart)&&(ritStart!=exonStart.rend())){
                    ritStart++;
                }
                ritStart_cur = ritStart;
                int lmb = orf.leftmostExonBegin(frame,i,true);
                while((lmb <= *ritStart) && (i-*ritStart <= Constant::max_exon_len) && (ritStart!=exonStart.rend())) {
                    if ((i - *ritStart>=3) && ((i-*ritStart+1)%3==frame) && (p >= dssminprob)) {
                        ec = new ExonCandidate;
                        ec->begin = *ritStart - 1;
                        ec->end = i - 1;
                        ec->dssScore = computeSpliceSiteScore(p, dssminprob, dssmaxprob);
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
                while(lmb <= (*ritASS).first && i-(*ritASS).first <= Constant::max_exon_len && ritASS != exonASS.rend()) {
                    if ((i-(*ritASS).first>=5) && (p >= dssminprob)) {
                        ec = new ExonCandidate;
                        ec->begin = (*ritASS).first + 2;
                        ec->end = i - 1;
                        ec->assScore = (*ritASS).second;
                        ec->dssScore = computeSpliceSiteScore(p, dssminprob, dssmaxprob);
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
        if (GeneticCode::isStopcodon(dna+i)) {
            for (frame=0; frame<=2; frame++) {
                ritASS=ritASS_cur;
                while ((i<(*ritASS).first)&&(ritASS!=exonASS.rend())){
                    ritASS++;
                }
                ritASS_cur = ritASS;
                while ((i-(*ritASS).first <= Constant::max_exon_len) && (ritASS!=exonASS.rend())) {
                    if ((i-(*ritASS).first>=3) && ((i-(*ritASS).first + 1)%3==frame) && ((*ritASS).first>=orf.leftmostExonBegin(0,i,true))) {
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

        // computing single genes on the reverse strand with at least Constant::min_coding_len
        if (onRStart(dna+i)) {
            ritRCStop=ritRCStop_cur;
            while ((i < *ritRCStop) && (ritRCStop != exonRCStop.rend())){
                ritRCStop++;
            }
            ritRCStop_cur=ritRCStop;
            while ((i-*ritRCStop <=  Constant::max_exon_len) && (ritRCStop!=exonRCStop.rend())) {
                if ((i-*ritRCStop)%3 == 0) {
                    if ((i-*ritRCStop) >= Constant::min_coding_len) {
                        ec = new ExonCandidate;
                        ec->begin = *ritRCStop;
                        ec->end = i + 2;
                        ec->type=rsingleGene;
                        candidates->push_back(ec);
                        break;
                    } else {
                        break;
                    }
                } else {
                    ritRCStop++;
                }
            };
        }

        // computing initials on the reverse strand with at least start codon plus base
        if (onRStart(dna+i)) {
            for (frame=0; frame<=2; frame++) {
                ritRDSS=ritRDSS_cur;
                while ((i<(*ritRDSS).first)&&(ritRDSS!=exonRDSS.rend())){
                    ritRDSS++;
                }
                ritRDSS_cur=ritRDSS;
                int lmb = orf.leftmostExonBegin(2,i,false);
                while((lmb<=(*ritRDSS).first+2)&&(i-(*ritRDSS).first<= Constant::max_exon_len)&&(ritRDSS!=exonRDSS.rend())) {
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
                while((lmb<=(*ritRDSS).first)&&(i-(*ritRDSS).first <= Constant::max_exon_len)&&(ritRDSS!=exonRDSS.rend())) {
                    if (i-(*ritRDSS).first>=5 && (p >= assminprob)) {
                        ec = new ExonCandidate;
                        ec->begin = (*ritRDSS).first + 2;
                        ec->end = i - 1;
                        ec->dssScore = (*ritRDSS).second;
                        ec->assScore = computeSpliceSiteScore(p, assminprob, assmaxprob);
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
                while ((i-*ritRCStop <= Constant::max_exon_len) && (ritRCStop!=exonRCStop.rend())) {
                    if (i-*ritRCStop == 3) {
                        break;
                    }
                    if ((i-*ritRCStop>=4) && ((i-*ritRCStop)%3==frame) && (p >= assminprob)) {
                        ec = new ExonCandidate;
                        ec->begin = *ritRCStop;
                        ec->end = i - 1;
                        ec->assScore = computeSpliceSiteScore(p, assminprob, assmaxprob);
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
    exoncands[s] = candidates;
    existingCandidates.push_back(getECHash(candidates));
}

// computes the aligned position of a base in an alignment and the 'block' where the base is found
/*
 * steffi: the following two functions are a total mess !!!
 * I only fixed out of range iterators which occasionally caused segmentation faults
 */
pair <int, int> GeneMSA::getAlignedPosition(AlignmentRow *as_ptr, int pos) {
    list<block>::iterator it;
    vector<int*>::iterator it_cmpStart;
    pair <int, int> alignedPos;
    it = upper_bound(as_ptr->sequence.begin(), as_ptr->sequence.end(), pos, compAlignedPos);
    if(it == as_ptr->sequence.begin()){
	alignedPos.first = -1;
	alignedPos.second = 0;
	return alignedPos;
    } else {
	it--;
    }
    it_cmpStart = upper_bound(as_ptr->cmpStarts.begin(), as_ptr->cmpStarts.end(), (pos + it->previousGaps), compCmpStarts);
    if(it_cmpStart == as_ptr->cmpStarts.begin()){
	throw ProjectError("in GeneMSA::getAlignedPosition");
    } else {
	it_cmpStart--;
    }
    while ((*it_cmpStart)==NULL && it_cmpStart != as_ptr->cmpStarts.begin()){
        it_cmpStart--;
    }
    vector<int*>::difference_type idx = distance(as_ptr->cmpStarts.begin(), it_cmpStart);
    alignedPos.second = idx;
    alignedPos.first = pos + it->previousGaps - (*(*it_cmpStart));
    return alignedPos;
}

// computes the real position of a base dependent on its position in the alignment
int GeneMSA::getRealPosition(AlignmentRow *ptr, int pos, int idx) {
    list<block>::iterator it;
    int realPos, alignedPos;
    if (ptr->cmpStarts[idx] != NULL) {
        alignedPos = *ptr->cmpStarts[idx] + pos;
    } else {
        return -1;
    }
    it = upper_bound(ptr->sequence.begin(), ptr->sequence.end(), alignedPos, compRealPos);
    if(it == ptr->sequence.begin()){
	return alignedPos;
    } else{
	it--;
    }
    realPos = alignedPos - it->previousGaps;
    return realPos;
}

/**
 * createOrthoExons
 */
void GeneMSA::createOrthoExons(vector<int> offsets, float consThres) {
    cout << "creating ortho exon for alignment" << *alignment << endl;
    int k = alignment->rows.size();
    int m = 0; // the number of nonempty rows
    int aliStart, aliEnd, chrExonStart, chrExonEnd;
    string key;
    for (size_t s=0; s<k; s++)
	if (alignment->rows[s])
	    m++;
    // an ortho exon candidate must have an EC in at least this many species (any subset allowed):
    int minEC = (consThres * m > 2.0)? m * consThres + 0.9999 : 2; 
    cout << "an OE must have at least " << minEC << " ECs" << endl;

    /*
     * Store for each exon candidate in alignment space (keys of format "aliStart:aliEnd;type")
     * a list of (speciesIdx, ExonCandidate*), e.g.
     * alignedECs["100:200:1"] = {(0, ec0), (3, ec3)}
     */

    map<string, list<pair<int,ExonCandidate*> > > alignedECs;
    map<string, list<pair<int,ExonCandidate*> > >::iterator aec;
    // map all exon candidates to alignment positions, where possible
    // this search in LINEAR in the length of all exon candidates
    // + the number of all alignment fragments
    for (size_t s=0; s<k; s++){
	if (alignment->rows[s] == NULL)
	    break;
	AlignmentRow *row = alignment->rows[s];
	vector<fragment>::const_iterator from = row->frags.begin();
	for(list<ExonCandidate*>::iterator ecit = exoncands[s]->begin(); ecit != exoncands[s]->end(); ++ecit){
	    chrExonStart = (*ecit)->getStart();
	    // go the the first fragment that may contain the ec start
	    while (from != row->frags.end() && from->chrPos + from->len - 1 < chrExonStart)
		++from;
	    if (from == row->frags.end())
		break; // have searched beyond the last alignment fragment => finished
	    aliStart = row->getAliPos(chrExonStart, from);
	    if (aliStart >= 0){ // left exon boundary mappable
		chrExonEnd = (*ecit)->getEnd();
		aliEnd = row->getAliPos(chrExonEnd, from);
		if (aliEnd >= 0){
		    // both exon boundaries were mappable
		    // store the ec in the hash
		    key = itoa(aliStart) + ":" + itoa(aliEnd) + ":" + itoa((*ecit)->type);
		    cout << "could map " << rsa->getSname(s) << " " << row->seqID << ":" << chrExonStart << ".." << chrExonEnd 
			 << " to " << " alignment coordinates " << aliStart << ".." << aliEnd << " key= " << key << endl;
		    aec = alignedECs.find(key);
		    if (aec == alignedECs.end()){ // insert new list
			list<pair<int,ExonCandidate*> > e;
			e.push_back(pair<int,ExonCandidate*> (s, *ecit));
			alignedECs.insert(pair<string,list<pair<int,ExonCandidate*> > >(key, e));
		    } else {// append new entry to existing list
			aec->second.push_back(pair<int,ExonCandidate*> (s, *ecit));
		    }
		}
	    }
	}
    }
    
    /*
     * Create one ortho exon candidate for each key to which at least minEC exon candidates mapped 
     */
    for (aec = alignedECs.begin(); aec != alignedECs.end(); ++aec){
	if (aec->second.size() >= minEC){
	    OrthoExon oe;
	    oe.ID = orthoExonID;
	    orthoExonID++;
	    oe.orthoex.resize(k, NULL);
	    for (list<pair<int,ExonCandidate*> >::iterator it = aec->second.begin(); it != aec->second.end(); ++it){
		int s = it->first;
		ExonCandidate *ec = it->second;
		if (oe.orthoex[s])
		    throw ProjectError("createOrthoExons: Have two exon candidates from the same species " + rsa->getSname(s) + " with the same key.");
		oe.orthoex[s] = ec;
	    }
	    orthoExonsList.push_back(oe);
	}
    }
}

// cut off incomplete codons at both boundaries of all exon candidates
void GeneMSA::cutIncompleteCodons(vector<ExonCandidate*> &orthoexon) {
    for (int i=0; i<orthoexon.size(); i++) {
        if (orthoexon[i] != NULL) {
	    if (isPlusExon(orthoexon[i]->type)){
		// Steffi: remove Stoppcodon on forward strand
		if (orthoexon[i]->type == singleGene || orthoexon[i]->type == terminal_exon){
		    orthoexon[i]->end -= 3;
		}
		orthoexon[i]->begin += getGFF3FrameForExon(orthoexon[i]);
		orthoexon[i]->end -= exonTypeReadingFrames[orthoexon[i]->type];
	    } else {
		// Steffi: remove Stoppcodon on reverse strand
		if (orthoexon[i]->type == rsingleGene || orthoexon[i]->type >= rterminal_0){
		    orthoexon[i]->begin += 3;
		}
		orthoexon[i]->begin += mod3(orthoexon[i]->len()-2+exonTypeReadingFrames[orthoexon[i]->type]);
		orthoexon[i]->end -= 2-exonTypeReadingFrames[orthoexon[i]->type];
	    }
        }
    }
}

// designed the aligned sequence of an ortholog exon candidate
string GeneMSA::getAlignedOrthoExon(AlignmentRow *as_ptr, ExonCandidate* ec, string seq, int offset) {
    if (as_ptr != NULL) {
        int alignedPosStart = getAlignedPosition(as_ptr, ec->begin + offset + 1).first;
        int idxStart = getAlignedPosition(as_ptr, ec->begin + offset + 1).second;
        int alignedBegin = (*(as_ptr->cmpStarts[idxStart])) + alignedPosStart;
        int alignedPosEnd = getAlignedPosition(as_ptr, ec->end + offset + 1).first;
        int idxEnd = getAlignedPosition(as_ptr, ec->end + offset + 1).second;
        int alignedEnd = (*(as_ptr->cmpStarts[idxEnd])) +  alignedPosEnd;
        list<block>::iterator it = upper_bound(as_ptr->sequence.begin(), as_ptr->sequence.end(), alignedBegin, compRealPos);
        list<block>::iterator it_prev = it;
        it_prev--;
        while (((*it_prev).begin + (*it_prev).length) < alignedEnd) {
            string gap = "";
            int numberGaps = (*it).previousGaps - (*it_prev).previousGaps;
            for (int j=0; j<numberGaps; j++) {
                gap += "-";
            }
            try {
                if ((*it_prev).begin != (*it).begin - ((*it).previousGaps - (*it_prev).previousGaps)  - (*it_prev).length) {
                    break;
                }
                if (gap != "") {
                    seq.insert((*it).begin - alignedBegin - numberGaps , gap);
                }
            }
            catch (...) {
                cout << "sequence is too short, out of range error" << endl;
                break;
            }
            it++;
            it_prev++;
        }
        return seq;
    } else {
        return "";
    }
}

// reads the computed dN/dS ratio from a file and gives the exon candidate a score: (100*(1-dN/dS)
void GeneMSA::readOmega(string file){
    double omega;
    string omegaFilename = file;
    ifstream OmegaFile;
    OmegaFile.open(omegaFilename.c_str(), ifstream::in);
    if (!OmegaFile) {
        cerr << "Could not find the file with the dN/dS ratios " << omegaFilename << "." << endl;
        throw PropertiesError( "GeneMSA::readOmega: Could not open this file!" );
    } else {
            for (list<OrthoExon>::iterator it_oe=orthoExonsWithOmega.begin(); it_oe!=orthoExonsWithOmega.end(); it_oe++) {
                OmegaFile>>omega;
                for (int j=0; j<(*it_oe).orthoex.size(); j++) {
                    if (((*it_oe).orthoex.at(j)!=NULL) && (!OmegaFile.eof())) {
                        (*it_oe).orthoex.at(j)->score=100*(1-omega);
                    }
                }
            }
    }
    OmegaFile.close();
}

void GeneMSA::openOutputFiles(){
    string outputdirectory;  //directory for output files
    try {
        outputdirectory = Properties::getProperty("/CompPred/outdir_orthoexons");
    } catch (...) {
        outputdirectory = "";
    }
    outputdirectory = expandHome(outputdirectory); //replace "~" by "$HOME"
    
    exonCands_outfiles.resize(tree->numSpecies());
    orthoExons_outfiles.resize(tree->numSpecies());
    geneRanges_outfiles.resize(tree->numSpecies());
    omega_outfiles.resize(tree->numSpecies());
    vector<string> species;
    tree->getSpeciesNames(species);
    for (int i=0; i<tree->numSpecies(); i++) {
        string file_exoncand = outputdirectory + "exonCands." + species[i] + ".gff3";
        ofstream *os_ec = new ofstream(file_exoncand.c_str());
        if (os_ec!=NULL) {
            exonCands_outfiles[i]=os_ec;
            (*os_ec) << PREAMBLE << endl;
            (*os_ec) << "#\n#-----  exon candidates  -----" << endl << "#" << endl;
        }
        string file_geneRanges = outputdirectory + "geneRanges." + species[i] + ".gff3";
        ofstream *os_gr = new ofstream(file_geneRanges.c_str());
        if (os_gr!=NULL) {
            geneRanges_outfiles[i]=os_gr;
            (*os_gr) << PREAMBLE << endl;
            (*os_gr) << "#\n#-----  possible gene ranges  -----" << endl << "#" << endl;
        }
        string file_orthoexon = outputdirectory + "orthoExons." + species[i] + ".gff3";
        ofstream *os_oe = new ofstream(file_orthoexon.c_str());
        if (os_oe) {
            orthoExons_outfiles[i]=os_oe;
            (*os_oe) << PREAMBLE << endl;
            (*os_oe) << "#\n#----- ortholog exons  -----" << endl << "#" << endl;
        }
        /*string file_omega = outputdirectory + "omegaExons." + species[i] + ".gff3";
        ofstream *os_omega = new ofstream(file_omega.c_str());
        if (os_omega) {
            omega_outfiles[i]=os_omega;
            (*os_omega) << PREAMBLE << endl;
            (*os_omega) << "#\n#----- exons with a dN/dS ratio smaller than one -----" << endl << "#" << endl;
	}*/
    }
    string paml_file = outputdirectory + "pamlfile.fa";
    ofstream *os_pf = new ofstream(paml_file.c_str());
    if (os_pf!=NULL) {
        pamlFile=os_pf;
    }
}

void GeneMSA::printStats(){
    for(int s=0; s<numSpecies(); s++){
	cout << rsa->getSname(s) << "\t";
	if (alignment->rows[s])
	    cout << alignment->rows[s]->strand << "\t" << getStart(s) << "\t" << getEnd(s);
	cout << endl;
    }
}

// writes the possible gene ranges of a species in the file 'geneRanges.speciesnames.gff3'
void GeneMSA::printGeneRanges() {
    if (!(this->exoncands.empty())) {
        for (int i=0; i < this->exoncands.size(); i++) {
	    ofstream &fstrm = *geneRanges_outfiles[i]; // write to 'geneRanges.speciesname[i].gff3'
            if (this->exoncands.at(i)!=NULL) {
                fstrm << this->getSeqID(i) << "\tGeneRange\t" << "exon\t" << this->getStart(i) + 1 << "\t"
		      << this->getEnd(i) + 1 << "\t0\t";
                if (this->getStrand(i) == plusstrand)
                    fstrm <<'+'<<"\t";
		else
                    fstrm << '-' << "\t";
                fstrm << ".\t" << "Name=" << geneRangeID << endl;
            }
        }
        geneRangeID++;
    }
}

//writes the exon candidates of a species of a dna segment in the file 'exonCands.species.gff3'
void GeneMSA::printExonCands(vector<int> offsets) {
    exonCandID.resize(this->exoncands.size());
    for (int j=0; j<exonCandID.size(); j++) {
        exonCandID[j] = 1;
    }
    if (!(this->exoncands.empty())) {
        for (int i=0; i<this->exoncands.size(); i++) {
	    ofstream &fstrm = *exonCands_outfiles[i]; // write to 'exonCands.speciesname[i].gff3'
            if (this->exoncands.at(i)!=NULL) {
                fstrm << "# sequence:\t" << rsa->getSname(i) << "\t"<<this->getStart(i) + 1 << "-" 
		      << this->getEnd(i) + 1<<"  "<<this->getEnd(i) - this->getStart(i)<< "bp" << endl;
                for (list<ExonCandidate*>::iterator it_exonCands=this->exoncands.at(i)->begin(); 
		     it_exonCands!=this->exoncands.at(i)->end(); it_exonCands++) {
                    fstrm << this->getSeqID(i)<< "\tEC\t" << "exon\t";
                    if (this->getStrand(i) == plusstrand) {
                        fstrm << (*it_exonCands)->begin + offsets[i]+1 << "\t" 
			      << (*it_exonCands)->end + offsets[i]+1<<"\t"<<(*it_exonCands)->score<<"\t"
			      << '+' << "\t";
                    } else {
			int chrLen = rsa->getChrLen(i, getSeqID(i));
                        fstrm << chrLen - ((*it_exonCands)->end + offsets[i]) << "\t"
			      << chrLen - ((*it_exonCands)->begin+ offsets[i]) << "\t"
			      << (*it_exonCands)->score << "\t" << '-' << "\t";
                    }
                    fstrm << getGFF3FrameForExon(*it_exonCands) << "\t" << "ID=" << exonCandID[i] << ";"
			  << "Name=" <<stateExonTypeIdentifiers[(*it_exonCands)->type] << endl;
                    exonCandID[i]++;
                }
            } else {
                fstrm << "#  no exon candidates found " << endl;
            }
        }
    } else {
        cout << "#  no exon candidates found at all" << endl;
    }
}

//writes the orthologue exons of the different species in the files 'orthoExons.species.gff3'
void GeneMSA::printOrthoExons(RandSeqAccess *rsa, vector<int> offsets) {
    string paml;
    if (!(this->orthoExonsList.empty())) {
        try {
            Properties::assignProperty("/CompPred/paml", paml);
        }
        catch (...) {
            paml="";
        }
        for (list<OrthoExon>::iterator it_oe=orthoExonsList.begin(); it_oe!=orthoExonsList.end(); it_oe++) {
            printExonsForPamlInput(rsa, *it_oe, offsets);
	    printSingleOrthoExon(*it_oe, offsets, true, it_oe->getOmega(), it_oe->getSubst());
            if (!paml.empty()) {
                readOmega(paml);
            }
        }
        if (!paml.empty()) {
	    //printExonWithOmega(offsets);
        }
    }
}

// writes the ortholog exons of the different species in the files 'orthoExons.species.gff3'
// files: write each one to a file for its species, if false to stdout
void GeneMSA::printSingleOrthoExon(OrthoExon &oe, vector<int> offsets, bool files, double omega, int numSub) {
    streambuf *stdout = cout.rdbuf();
    for (int j=0; j<oe.orthoex.size(); j++) {
	ExonCandidate *ec = oe.orthoex.at(j);
	if (files)
	    cout.rdbuf(orthoExons_outfiles[j]->rdbuf()); // write to 'orthoExons.speciesname[j].gff3'
        if (ec != NULL) {
	    cout << getSeqID(j) << "\tOE1\t" << "exon" << "\t";
            if (this->getStrand(j) == plusstrand){ // strand of alignment
		cout << ec->begin + offsets[j]+1 << "\t" << ec->end + offsets[j]+1;
            } else {
		int chrLen = rsa->getChrLen(j, getSeqID(j));
                cout << chrLen - (ec->end + offsets[j]) << "\t"
		     << chrLen - (ec->begin + offsets[j]);
            }
	    cout << "\t" << ec->score << "\t" << (isPlusExon(ec->type)? '+' : '-') << "\t"; // strand of exon
            cout << getGFF3FrameForExon(ec) << "\t" << "ID=" << oe.ID << ";Name=" << oe.ID << ";Note=" 
		 << stateExonTypeIdentifiers[ec->type];
	    if (omega >= 0.0){
		cout << ";omega=" << omega;
		//cout << "|" << omega;  // for viewing in gBrowse use this style instead
	    }
	    if (numSub >= 0){
		cout << ";subst=" << numSub; // number of substitutions
	        //cout << "|" << numSub; // for viewing in gBrowse use this style instead
	    }
	    cout << endl;
        }
    }
   cout.rdbuf(stdout); // reset to standard output again 
}

// prints exons with an computed omega smaller than one
void GeneMSA::printExonWithOmega(vector<int> offsets) {
    streambuf *console = cout.rdbuf();
    if (!(this->orthoExonsWithOmega.empty())) {
        for (list<OrthoExon>::iterator it_oe=this->orthoExonsWithOmega.begin(); it_oe!=this->orthoExonsWithOmega.end(); it_oe++) {
            for (int j=0; j<(*it_oe).orthoex.size(); j++) {
                cout.rdbuf(omega_outfiles[j]->rdbuf()); //direct cout to 'omegaExons.speciesname.gff3'
                if (((*it_oe).orthoex.at(j)!=NULL) && ((*it_oe).orthoex.at(j)->score > 0)) {
                    cout<<this->getSeqID(j)<< "\tOO1\t"<<"exon"<<"\t";
                    if (this->getStrand(j) == plusstrand) {
                        cout <<(*it_oe).orthoex.at(j)->begin + offsets[j]+1<<"\t"<<(*it_oe).orthoex.at(j)->end + offsets[j]+1<<"\t"<<(*it_oe).orthoex.at(j)->score<<"\t";
                        cout<<'+'<<"\t";
                    } else {
			int chrLen = rsa->getChrLen(j, getSeqID(j));
                        cout << chrLen - ((*it_oe).orthoex.at(j)->end + offsets[j])<<"\t"
			     << chrLen - ((*it_oe).orthoex.at(j)->begin + offsets[j])<<"\t";
                        cout<<(*it_oe).orthoex.at(j)->score<<"\t"<<'-'<<"\t";
                    }
                    cout<<getGFF3FrameForExon((*it_oe).orthoex.at(j)) <<"\t";
                    cout<<"ID="<<(*it_oe).ID<<";Name="<<(*it_oe).ID<<";Note="<<stateExonTypeIdentifiers[(*it_oe).orthoex.at(j)->type]<<endl;
                }
            }
        }
    }
    cout.rdbuf(console); //reset to standard output again
}

// computes the sequence for the programm paml, so that every codon in the aligned sequence has no gaps and is in reading frame 0
// maybe the method looks pedestrianly, but under specific circumstances causes PAML a segmentation fault and I try to avoid it
vector <string> GeneMSA::getSeqForPaml(Alignment *it_ab,  vector<ExonCandidate*> oe, vector<string> seq, vector<int> offsets, vector<int> speciesIdx) {
    vector<string> pamlseq;
    int comparableSeq;
    vector<int> firstBaseCodon;
    vector<int> orthoExonStart;
    vector<bool> completeCodon;
    vector<int> codonCount;
    int alignedBase = 0;
    for (vector<int>::iterator it = speciesIdx.begin(); it!=speciesIdx.end(); it++) {
        orthoExonStart.push_back(getAlignedPosition(it_ab->rows.at(*it), oe[*it]->begin + offsets[*it] + 1).first);
        firstBaseCodon.push_back(oe[*it]->begin + offsets[*it] + 1);
        completeCodon.push_back(true);
        codonCount.push_back(0);
    }
    while (alignedBase < seq[0].length()) {
        comparableSeq = 0;
        for (int i = 0; i<seq.size(); i++) {
            completeCodon[i]=true;
            for (int k=alignedBase; k < alignedBase + 3; k++) {
                if (seq.at(i)[k] == '-') {
                    completeCodon[i] = false;
                }
            }
            if ((getAlignedPosition(it_ab->rows.at(speciesIdx[i]), firstBaseCodon[i]).first - orthoExonStart[i] == alignedBase) && (completeCodon[i])) {
                comparableSeq++;
            }
        }
        if (comparableSeq > 1) {
            for (int i = 0; i < completeCodon.size(); i++) {
                if ((getAlignedPosition(it_ab->rows.at(speciesIdx[i]), firstBaseCodon[i]).first - orthoExonStart[i] == alignedBase) && (completeCodon[i])) {
                    firstBaseCodon[i] = firstBaseCodon[i] + 3;
                    codonCount[i]++;
                } else {
                    seq[i].replace(alignedBase, 3, "---");
                    if ((getAlignedPosition(it_ab->rows.at(speciesIdx[i]), firstBaseCodon[i]).first) - orthoExonStart[i] <= alignedBase ) {
                        firstBaseCodon[i] = firstBaseCodon[i] + 3;
                    }
                }
            }
        } else {
            for (int i = 0; i<completeCodon.size(); i++) {
                seq[i].replace(alignedBase, 3, "---");
                if ((getAlignedPosition(it_ab->rows.at(speciesIdx[i]), firstBaseCodon[i]).first) - orthoExonStart[i] <= alignedBase ) {
                    firstBaseCodon[i] = firstBaseCodon[i] + 3;
                }
            }
        }
        alignedBase = alignedBase + 3;
    }
    for (int i=0; i<seq.size(); i++) {
        pamlseq.push_back(seq[i]);
    }
    for (int i=0; i<seq.size(); i++) {
        if (codonCount[i] < 2) {
            pamlseq.clear();
        }
    }
    return pamlseq;
}

//writes the exons for the paml input file '
void GeneMSA::printExonsForPamlInput(RandSeqAccess *rsa, OrthoExon &oe, vector<int> offsets) {
    streambuf *console = cout.rdbuf();  //save old buf
    cout.rdbuf(pamlFile->rdbuf());  //redirect cout to 'pamlFile.fa'
    vector<ExonCandidate*> orthoExon;
    if (!(orthoExonsList.empty())) {
        int noSpecies = 0;
        AnnoSequence *seqRange = NULL;
        vector<int> speciesIdx;
        int alignedOrthoExonLength=0;

        for (int j=0; j<oe.orthoex.size(); j++) {
            if (oe.orthoex[j] != NULL) {
                ExonCandidate *ec = new ExonCandidate(oe.orthoex[j]);
                orthoExon.push_back(ec);
            } else {
                orthoExon.push_back(NULL);
            }
        }
        cutIncompleteCodons(orthoExon);

        for (int i=0; i<orthoExon.size(); i++) {  
            if ((orthoExon.at(i) != NULL) && orthoExon[i]->len() >= 21 ) { //sequence has to be greater than 20 bases
                seqRange = rsa->getSeq(i, getSeqID(i), orthoExon[i]->begin + offsets[i], 
				       orthoExon[i]->end + offsets[i], getStrand(i));
                if (i == 0) {
                    alignedOrthoExonLength = getAlignedOrthoExon(alignment->rows.at(0), orthoExon[0],
								 seqRange->sequence, offsets[0]).length();
                }
                if ((((alignedOrthoExonLength) % 3) == 0) && 
		    (alignedOrthoExonLength) == (getAlignedOrthoExon(alignment->rows.at(i), orthoExon[i],
								     seqRange->sequence, offsets[i]).length())) {
                    speciesIdx.push_back(i);
                    noSpecies++;
                }
            }
        }

        if (noSpecies >= 2 /*&& speciesIdx[1]==2*/) {
            orthoExonsWithOmega.push_back(oe);
            vector<string> pamlSeq;
            for (vector<int>::iterator it = speciesIdx.begin(); it!=speciesIdx.end(); it++) {
                if (alignment->rows[*it]->strand == minusstrand) {
		    int chrLen = rsa->getChrLen(*it, getSeqID(*it));
                    seqRange = rsa->getSeq(rsa->getSname(*it), this->getSeqID(*it), chrLen - (orthoExon[*it]->end + offsets[*it] + 1),
                            chrLen - (orthoExon[*it]->begin + offsets[*it] + 1), this->getStrand(*it));
                } else {
                    seqRange = rsa->getSeq(rsa->getSname(*it), this->getSeqID(*it), orthoExon[*it]->begin + offsets[*it], orthoExon[*it]->end + offsets[*it], this->getStrand(*it));
                }
                string sequence = getAlignedOrthoExon(alignment->rows.at(*it), orthoExon[(*it)], seqRange->sequence, offsets[*it]);
                if (!isPlusExon(orthoExon[*it]->type)) {
                    int n = sequence.length();
                    for (int j=0; j<n; j++) {
                        if ( sequence[j] != '-') {
                            sequence[j] = wcComplement(sequence[j]);
                        }
                    }
                }
                pamlSeq.push_back(sequence);
            }
            pamlSeq = getSeqForPaml(alignment, orthoExon, pamlSeq, offsets, speciesIdx);
            if (! pamlSeq.empty()) {
                cout << noSpecies << "  " << alignedOrthoExonLength << endl;
                int k=0;
                for (vector<int>::iterator it = speciesIdx.begin(); it!=speciesIdx.end(); it++) {
                    cout << rsa->getSname(*it) << endl;
                    if (!isPlusExon(orthoExon[*it]->type))
                        pamlSeq[k] = reverseString(pamlSeq[k]);
                    cout << pamlSeq[k] << endl;
                    k++;
                }
                cout << endl;
		/*if (noSpecies==2){ // just for testing
		    cout << setw(8) << "codons" << setw(10) << "syn sub" << setw(12) << "nonsyn sub"
			 << setw(8) << "omega" << endl;
		    double t(0.5);
		    int numAliCodons, numSynSubst, numNonSynSubst;
		    double omega = codonevo->estOmegaOnSeqPair(pamlSeq[0].c_str(), pamlSeq[1].c_str(), t,
							       numAliCodons, numSynSubst, numNonSynSubst);
		    cout << setw(8) << numAliCodons << setw(10) << numSynSubst << setw(12) << numNonSynSubst
			 << setw(7) << omega << endl;
		    oe.setOmega(omega);
		    oe.setSubst(numSynSubst + numNonSynSubst);
		    //printSingleOrthoExon(oe, offsets, true, omega, numSynSubst + numNonSynSubst);
		    printSingleOrthoExon(oe, offsets, false, omega, numSynSubst + numNonSynSubst);
		}*/
		int subst;
		//TODO: scale branch lenghts to one substitution per codon per time unit
		double omega = codonevo->estOmegaOnSeqTuple(pamlSeq, speciesIdx, tree, subst);
		oe.setOmega(omega);
		oe.setSubst(subst);
		printSingleOrthoExon(oe, offsets, false, omega, subst);
            }
        }
    }
    for (int j=0; j<orthoExon.size(); j++) {
	if (orthoExon[j]!=NULL) {
            delete orthoExon[j];
	}
    }
    cout.rdbuf(console); //reset to standard output again
}

void GeneMSA::closeOutputFiles(){
    for (int i=0; i<tree->numSpecies(); i++) {
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
        if (omega_outfiles[i]) {
            if(omega_outfiles[i]->is_open()) {
                omega_outfiles[i]->close();
                delete omega_outfiles[i];
            }
        }
    }
    if (pamlFile) {
        if(pamlFile->is_open()) {
            pamlFile->close();
            delete pamlFile;
        }
    }
}

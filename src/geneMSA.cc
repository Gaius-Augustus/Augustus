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

/*
 * constructor of GeneMSA
 */
GeneMSA::GeneMSA(RandSeqAccess *rsa, Alignment *a) {
    int maxDNAPieceSize = Properties::getIntProperty( "maxDNAPieceSize" );

    this->rsa = rsa;
    alignment = a;
    if (!alignment)
	throw ProjectError("Internal error in GeneMSA: Alignment missing.");    
    if (alignment->numRows() != rsa->getNumSpecies())
	throw ProjectError("Error in GeneMSA: Number of species in alignment (" + itoa(alignment->numRows())
			   + ") is not matching the one in the tree (" + itoa(rsa->getNumSpecies()) + ").");
    exoncands.resize(alignment->numRows(), NULL);
    /** construct the gene ranges
     * now: simple copy. TODO: extend region when apparently part of the alignment is missing8
     * human   ***********---*******************
     * mouse   *******-------******************-
     * chicken ***********----------------------
     *                        ^
     *                        | extend range here
     */
    starts.resize(alignment->numRows(), -1);
    ends.resize(alignment->numRows(), -1);
    offsets.resize(alignment->numRows(), 0);
    for (size_t s=0; s < alignment->numRows(); s++){
	if (alignment->rows[s]){
	    int chrLen = rsa->getChrLen(s, alignment->rows[s]->seqID);
	    switch (getStrand(s))
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
	    offsets[s] = (getStrand(s) == plusstrand)? starts[s] : rsa->getChrLen(s, getSeqID(s)) - 1 - ends[s];
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
}

int GeneMSA::getGFF3FrameForExon(ExonCandidate *ec) {
    if (isPlusExon(ec->type))
        return mod3(3-(exonTypeReadingFrames[ec->type] - (ec->end - ec->begin + 1)));
    else
        return mod3(2-exonTypeReadingFrames[ec->type]);
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
    cout << "createExonCands: " << s << "\t" << assmotifqthresh 
	 << "\t" << assqthresh << "\t" << dssqthresh << endl;
    int n = strlen(dna);
    int frame;
    Double p;
    list<int> exonStart;
    list<int> exonRCStop;
    list< pair<int, Double> > exonASS;
    list< pair<int, Double> > exonRDSS;
    ExonCandidate *ec;
    list<ExonCandidate*> *candidates = new list<ExonCandidate*>;
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
    cout << "Found " << candidates->size() << " ECs on species " << rsa->getSname(s) << endl;
    exoncands[s] = candidates;
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

/**
 * createOrthoExons
 */
void GeneMSA::createOrthoExons(float consThres) {
    cout << "Creating ortho exon for alignment" << endl << *alignment << endl;
    int k = alignment->rows.size();
    int m = 0; // the number of nonempty rows
    int aliStart, aliEnd, chrExonStart, chrExonEnd;
    string key;
    for (size_t s=0; s<k; s++)
	if (alignment->rows[s])
	    m++;
    // an ortho exon candidate must have an EC in at least this many species (any subset allowed):
    int minEC = (consThres * m > 2.0)? m * consThres + 0.9999 : 2; 
    // cout << "OEs in this gene range must have at least " << minEC << " ECs" << endl;

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
	    continue;
	int offset = offsets[s];
	AlignmentRow *row = alignment->rows[s];
	vector<fragment>::const_iterator from = row->frags.begin();
	for(list<ExonCandidate*>::iterator ecit = exoncands[s]->begin(); ecit != exoncands[s]->end(); ++ecit){
	    chrExonStart = (*ecit)->getStart() + offset;
	    // go the the first fragment that may contain the ec start
	    while (from != row->frags.end() && from->chrPos + from->len - 1 < chrExonStart)
		++from;
	    if (from == row->frags.end())
		break; // have searched beyond the last alignment fragment => finished
	    aliStart = row->getAliPos(chrExonStart, from);
	    if (aliStart >= 0){ // left exon boundary mappable
		chrExonEnd = (*ecit)->getEnd() + offset;
		aliEnd = row->getAliPos(chrExonEnd, from);
		if (aliEnd >= 0){
		    // both exon boundaries were mappable
		    // store the ec in the hash
		    int lenMod3 = (chrExonEnd - chrExonStart + 1) % 3;
		    // this key is at the same time the (lexicographic) sorting criterion for printing all OrthoExons
		    stringstream key;
		    key << setfill('0') << setw(8) << itoa(aliStart) <<  " " << setfill('0') << setw(8) 
			<< itoa(aliEnd) << " " <<  setfill('0') << setw(2) << itoa((*ecit)->type) << " " << itoa(lenMod3);
		    cout << "Could map " << rsa->getSname(s) << " " << row->seqID << ":" << chrExonStart << ".." << chrExonEnd 
			 << " to " << " alignment coordinates " << aliStart << ".." << aliEnd << " key = " << key.str() << endl;
		    aec = alignedECs.find(key.str());
		    if (aec == alignedECs.end()){ // insert new list
			list<pair<int,ExonCandidate*> > e;
			e.push_back(pair<int,ExonCandidate*> (s, *ecit));
			alignedECs.insert(pair<string,list<pair<int,ExonCandidate*> > >(key.str(), e));
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
    cout << "Going through all ortho exons" << endl;
    for (aec = alignedECs.begin(); aec != alignedECs.end(); ++aec){
	if (aec->second.size() >= minEC){
	    cout << aec->first << endl;
	    OrthoExon oe;
	    oe.ID = orthoExonID;
	    orthoExonID++;
	    oe.orthoex.resize(k, NULL);
	    cout << "OE\n";
	    for (list<pair<int,ExonCandidate*> >::iterator it = aec->second.begin(); it != aec->second.end(); ++it){
		int s = it->first;
		ExonCandidate *ec = it->second;
		cout << rsa->getSname(s) << "\t" << ec->getStart() + offsets[s] << ".."
		     << ec->getEnd() + offsets[s] << "\t" << *ec << endl;
		if (oe.orthoex[s])
		    throw ProjectError("createOrthoExons: Have two exon candidates from the same species " 
				       + rsa->getSname(s) + " with the same key.");
		oe.orthoex[s] = ec;
	    }
	    orthoExonsList.push_back(oe);
	}
    }
}

// cut off incomplete codons at both boundaries of all exon candidates
// TODO: delete this
void GeneMSA::cutIncompleteCodons(OrthoExon &oe) {
    for (int i=0; i<oe.orthoex.size(); i++) {
        if (oe.orthoex[i] != NULL) {
	    if (isPlusExon(oe.orthoex[i]->type)){
		// Steffi: remove Stoppcodon on forward strand
		if (oe.orthoex[i]->type == singleGene || oe.orthoex[i]->type == terminal_exon){
		    oe.orthoex[i]->end -= 3;
		}
		oe.orthoex[i]->begin += getGFF3FrameForExon(oe.orthoex[i]);
		oe.orthoex[i]->end -= exonTypeReadingFrames[oe.orthoex[i]->type];
	    } else {
		// Steffi: remove Stoppcodon on reverse strand
		if (oe.orthoex[i]->type == rsingleGene || oe.orthoex[i]->type >= rterminal_0){
		    oe.orthoex[i]->begin += 3;
		}
		oe.orthoex[i]->begin += mod3(oe.orthoex[i]->len()-2+exonTypeReadingFrames[oe.orthoex[i]->type]);
		oe.orthoex[i]->end -= 2-exonTypeReadingFrames[oe.orthoex[i]->type];
	    }
        }
    }
}

// designed the aligned sequence of an ortholog exon candidate
// TODO: delete this
string GeneMSA::getAlignedOrthoExon(AlignmentRow *row, ExonCandidate* ec, string seq, int offset) {
    //if (row == NULL)
	return "";
    int alignedPosStart = getAlignedPosition(row, ec->begin + offset + 1).first;
    int idxStart = getAlignedPosition(row, ec->begin + offset + 1).second;
    int alignedBegin = (*(row->cmpStarts[idxStart])) + alignedPosStart;
    int alignedPosEnd = getAlignedPosition(row, ec->end + offset + 1).first;
    int idxEnd = getAlignedPosition(row, ec->end + offset + 1).second;
    int alignedEnd = (*(row->cmpStarts[idxEnd])) +  alignedPosEnd;
    list<block>::iterator it = upper_bound(row->sequence.begin(), row->sequence.end(), alignedBegin, compRealPos);
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

// writes the exon candidates of all species into the file 'exonCands.species.gff3'
void GeneMSA::printExonCands() {
    exonCandID.resize(numSpecies(), 1);

    if (!exoncands.empty()) {
        for (int s=0; s < numSpecies(); s++) {
	    ofstream &fstrm = *exonCands_outfiles[s]; // write to 'exonCands.speciesname[i].gff3'
	    list<ExonCandidate*>* sec = exoncands[s];
            if (sec) {
                fstrm << "# sequence:\t" << rsa->getSname(s) << "\t" << getStart(s) + 1 << "-" 
		      << getEnd(s) + 1 << "  " << getEnd(s) - getStart(s) << "bp" << endl;
                for (list<ExonCandidate*>::iterator ecit = sec->begin(); ecit != sec->end(); ++ecit) {
                    fstrm << getSeqID(s)<< "\tEC\t" << "exon\t";
                    if (getStrand(s) == plusstrand) {
                        fstrm << (*ecit)->begin + offsets[s] + 1 << "\t" << (*ecit)->end + offsets[s] + 1 
			      << "\t" << (*ecit)->score << "\t" << '+' << "\t";
                    } else {
			int chrLen = rsa->getChrLen(s, getSeqID(s));
                        fstrm << chrLen - ((*ecit)->end + offsets[s]) << "\t"
			      << chrLen - ((*ecit)->begin+ offsets[s]) << "\t"
			      << (*ecit)->score << "\t" << '-' << "\t";
                    }
                    fstrm << getGFF3FrameForExon(*ecit) << "\t" << "ID=" << exonCandID[s] << ";"
			  << "Name=" <<stateExonTypeIdentifiers[(*ecit)->type] << endl;
                    exonCandID[s]++;
                }
            } else {
                fstrm << "#  no exon candidates found " << endl;
            }
        }
    } else {
        cout << "#  no exon candidates found at all" << endl;
    }
}

// writes all ortholog exons of all species in the files 'orthoExons.species.gff3'
// orthoexons are sorted by alignment start coordinate
void GeneMSA::printOrthoExons(RandSeqAccess *rsa) {
    if (orthoExonsList.empty())
	return;
    for (list<OrthoExon>::iterator oeit = orthoExonsList.begin(); oeit != orthoExonsList.end(); ++oeit)
	printSingleOrthoExon(*oeit, true, oeit->getOmega(), oeit->getSubst());
}

// writes the ortholog exons on one OrthoExon into the files 'orthoExons.species.gff3'
// files: write each one to a file for its species, if false to stdout
void GeneMSA::printSingleOrthoExon(OrthoExon &oe, bool files, double omega, int numSub) {
    streambuf *stdout = cout.rdbuf();
    for (int s=0; s < numSpecies(); s++) {
	ExonCandidate *ec = oe.orthoex.at(s);
	if (files)
	    cout.rdbuf(orthoExons_outfiles[s]->rdbuf()); // write to 'orthoExons.speciesname[s].gff3'
        if (ec != NULL) {
	    cout << getSeqID(s) << "\tOE1\t" << "exon" << "\t";
            if (getStrand(s) == plusstrand){ // strand of alignment
		cout << ec->begin + offsets[s]+1 << "\t" << ec->end + offsets[s]+1;
            } else {
		int chrLen = rsa->getChrLen(s, getSeqID(s));
                cout << chrLen - (ec->end + offsets[s]) << "\t" << chrLen - (ec->begin + offsets[s]);
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

/** getCodonAlignment
 *
 * Example input (all ECs have agree in phase at both boundaries)
 *  a|c g t|t t g|a t g|t c a|a a
 *  a|c g t|t t g|a t g|t c a|a a
 *  a|c g t|t t g|a t g|t c a|a a
 *  a|c g t|t t g|a t g|t c a|a a
 *
 * Example output:
 *
 *
 *
 *
 *
 *
 */
vector<string> GeneMSA::getCodonAlignment(OrthoExon const &oe, vector<AnnoSequence> const &seqRanges) {
    vector<string> rowstrings(numSpecies(), "");
    return rowstrings;
}


// computes and sets the Omega = dN/dS attribute to all OrthoExons
void GeneMSA::computeOmegas(vector<AnnoSequence> const &seqRanges) {
    int subst;
    
    for (list<OrthoExon>::iterator oe = orthoExonsList.begin(); oe != orthoExonsList.end(); ++oe){
	// TODO: retain for each species an iterator for the first alignment fragment to start the search in
	vector<string> rowstrings = getCodonAlignment(*oe, seqRanges);

	// TODO: scale branch lenghts to one substitution per codon per time unit
	double omega = codonevo->estOmegaOnSeqTuple(rowstrings, tree, subst);
	oe->setOmega(omega);
	oe->setSubst(subst);
    }
}

void GeneMSA::closeOutputFiles(){
    for (int i=0; i<tree->numSpecies(); i++) {
        if (exonCands_outfiles[i] && exonCands_outfiles[i]->is_open()) {
	    exonCands_outfiles[i]->close();
	    delete exonCands_outfiles[i];
	}
        if (geneRanges_outfiles[i] && geneRanges_outfiles[i]->is_open()) {
	    geneRanges_outfiles[i]->close();
	    delete geneRanges_outfiles[i];
	}
        if (orthoExons_outfiles[i] && orthoExons_outfiles[i]->is_open()) {
	    orthoExons_outfiles[i]->close();
	    delete orthoExons_outfiles[i];
	}
        if (omega_outfiles[i] && omega_outfiles[i]->is_open()) {
	    omega_outfiles[i]->close();
	    delete omega_outfiles[i];
	}
    }
}

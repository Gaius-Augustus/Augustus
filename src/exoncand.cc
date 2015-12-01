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

bool isPlusExon(ExonType t){
    return t < rsingleGene;
}

const int exonTypeReadingFrames[EXON_TYPES-1] = { 0,0,1,2,0,1,2,0,         // forward exons
						  2,2,0,1,2,0,1,2};         // reverse exons

const char* stateExonTypeIdentifiers[EXON_TYPES-1] =
{ "SINGLE", "INITIAL0", "INITIAL1", "INITIAL2", "INTERNAL0", "INTERNAL1", "INTERNAL2", "TERMINAL",
  "RSINGLE", "RINITIAL", "RINTERNAL0", "RINTERNAL1", "RINTERNAL2", "RTERMINAL0", "RTERMINAL1", "RTERMINAL2"};

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

int ExonCandidate::getFirstCodingBase(){
    if (isPlusExon(type)){
	return begin + gff3Frame();
    } else {
	if (type == rsingleGene || type >= rterminal_0)
	    return begin + 3;
	return begin + mod3(len() - 2 + exonTypeReadingFrames[type]);
    }
}

int ExonCandidate::getLastCodingBase(){
    if (isPlusExon(type)){
	if (type == singleGene || type == terminal_exon)
	    return end - 3;
	return end - exonTypeReadingFrames[type];
    } else {
	return end - 2 + exonTypeReadingFrames[type];
    }
}

int ExonCandidate::gff3Frame() {
    if (isPlusExon(type))
        return mod3(3-(exonTypeReadingFrames[type] - (end - begin + 1)));
    else
        return mod3(2-exonTypeReadingFrames[type]);
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

// create a key for a map function to find ortholog exons quickly
string ExonCandidate::key() {
   return (itoa(this->begin) + ":" + itoa(this->end) + ":" + itoa(this->type));
}

// keys encodes all of: chrStart chrEnd type lenMod3
int_fast64_t ExonCandidate::getKey(){
    
    int_fast64_t start = begin;
    int_fast64_t len = end - begin;
    int lenMod3 = (end - begin + 1) % 3;
    int_fast64_t key = (start << 22 ) // 42 bits
        + ((len) << 7)                // 15 bits
        + (type << 2)                 //  5 bits
        + lenMod3;                    //  2 bits
    return key;
}


ostream& operator<<(ostream& strm, const ExonCandidate &ec){
    strm << ec.begin << "\t" << ec.end << "\ttype=" << ec.type << "\tscore=" << ec.score 
	 << "\tupScore=" << ec.getUpScore() << "\tdownScore=" << ec.getDownScore();
    return strm;
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

// computes the score for the splice sites of an exon candidate
Double computeSpliceSiteScore(Double exonScore, Double minProb, Double maxProb) {
    Double score = 0;
    score = (log(exonScore/minProb)/log(maxProb/minProb));
    return score;
}


// computes exon candidates and inserts them into the hash of ECs if they do not exist already
void findExonCands(map<int_fast64_t, ExonCandidate*> &ecs, map<int_fast64_t, ExonCandidate*> &addECs, const char *dna, int minLen, double assmotifqthresh, double assqthresh, double dssqthresh){
    int n = strlen(dna);
    int frame;
    Double p;
    list<int> exonStart;
    list<int> exonRCStop;
    list< pair<int, Double> > exonASS;
    list< pair<int, Double> > exonRDSS;
    
    map<int_fast64_t, ExonCandidate*>::iterator ecit;
    Double assminprob, dssminprob, assmaxprob, dssmaxprob;
    bool withSplicing = (IntronModel::assBinProbs.nbins > 0); // eukaryotic species?
    if (withSplicing){
	assminprob = IntronModel::assBinProbs.getMinProb(assqthresh) * IntronModel::getAssMotifProbThreshold(assmotifqthresh);
	dssminprob = IntronModel::dssBinProbs.getMinProb(dssqthresh);
	assmaxprob = IntronModel::assBinProbs.getMinProb(0.99) * IntronModel::getAssMotifProbThreshold(0.9999);
	dssmaxprob = IntronModel::dssBinProbs.getMinProb(0.99);
    } // otherwise it is a bacterium and exon candidates with splice sites don't exist
    OpenReadingFrame orf(dna, Constant::max_exon_len, n);
    // preprocessing all left coordinates of an exon candidate interval
    for (int i=0; i<=n - 1; i++) {
        pair<int, Double> ssWithScore;
        // positions of all startcodons "atg"
        if (onStart(dna+i)) {
            exonStart.push_back(i + 1);
        }
        // positons of all ASSs "ag"
        if (withSplicing && onASS(dna+i) && (i + Constant::ass_whole_size() - Constant::ass_start < n)) {
            p = IntronModel::aSSProb(i - Constant::ass_upwindow_size - Constant::ass_start, true);
            if (p >= assminprob ) {
                ssWithScore.first = i;
                ssWithScore.second = computeSpliceSiteScore(p, assminprob, assmaxprob);
                exonASS.push_back(ssWithScore);
            }
        }
        // positions of all reverse DSS "ac"
        if (withSplicing && onRDSS(dna+i) && (i + Constant::dss_whole_size() - Constant::dss_end < n)) {
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
                    ExonCandidate *ec = new ExonCandidate;
                    ec->begin = *ritStart - 1;
                    ec->end = i + 2;
                    ec->type = singleGene;
		    if(ec->len() >= minLen){
			int_fast64_t key = ec->getKey();
			ecit = ecs.find(key);
			if (ecit == ecs.end()){ // insert new EC                                                                                                                                
                            ecs.insert(pair<int_fast64_t, ExonCandidate*>(key,ec));
                            addECs.insert(pair<int_fast64_t, ExonCandidate*>(key,ec));
                        }
			else{
                            delete ec;
                        }
		    }
		    else{
			delete ec;
		    }
                }
                ritStart++;
            }
        }
	
        // computing initial exons on the forward strand with at least startcodon plus base
        if (withSplicing && onDSS(dna + i) && (i + Constant::dss_whole_size() - Constant::dss_start  < n)) {
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
                        ExonCandidate *ec = new ExonCandidate;
                        ec->begin = *ritStart - 1;
                        ec->end = i - 1;
                        ec->setDownScore(computeSpliceSiteScore(p, dssminprob, dssmaxprob));
                        if (frame == 0) {
                            ec->type = initial_0;
                        } else if (frame == 1) {
                            ec->type = initial_1;
                        } else {
                            ec->type = initial_2;
                        }
			if(ec->len() >= minLen){
			    int_fast64_t key = ec->getKey();
			    ecit = ecs.find(key);
			    if (ecit == ecs.end()){ // insert new EC                                                                                                                                
				ecs.insert(pair<int_fast64_t, ExonCandidate*>(key,ec));
				addECs.insert(pair<int_fast64_t, ExonCandidate*>(key,ec));
			    }
			    else{
				ecit->second->setDownScore(ec->getDownScore());
				delete ec;
			    }
			}
			else{
			    delete ec;
			}
		    }
                    ritStart++;
                }
            }
	    
            // computing internals on the forward strand with at least one codon
	    if (withSplicing) {
		for (frame=0; frame<=2; frame++) {
		    ritASS=ritASS_cur;
		    while (i < (*ritASS).first && ritASS!=exonASS.rend())
			ritASS++;
		    ritASS_cur = ritASS;
		    int lmb = orf.leftmostExonBegin(frame, i, true);
		    while(lmb <= (*ritASS).first && i-(*ritASS).first <= Constant::max_exon_len && ritASS != exonASS.rend()) {
			if ((i-(*ritASS).first>=5) && (p >= dssminprob)) {
			    ExonCandidate *ec = new ExonCandidate;
			    ec->begin = (*ritASS).first + 2;
			    ec->end = i - 1;
			    ec->setUpScore((*ritASS).second);
			    ec->setDownScore(computeSpliceSiteScore(p, dssminprob, dssmaxprob));
			    if (frame == 0) {
				ec->type = internal_0;
			    } else if (frame==1) {
				ec->type = internal_1;
			    } else {
				ec->type = internal_2;
			    }
			    if(ec->len() >= minLen){
				int_fast64_t key = ec->getKey();
				ecit = ecs.find(key);
				if (ecit == ecs.end()){ // insert new EC                                                                                                                                
				    ecs.insert(pair<int_fast64_t, ExonCandidate*>(key,ec));
				    addECs.insert(pair<int_fast64_t, ExonCandidate*>(key,ec));
				}
				else{
				    ecit->second->setDownScore(ec->getDownScore());
				    ecit->second->setUpScore(ec->getUpScore());
				    delete ec;
				}
			    }
			    else{
				delete ec;
			    }
			}
			ritASS++;
		    }
		}
	    }
	}
	
        // computing terminal exons on the forward strand with at least one base stopcodon
        if (withSplicing && GeneticCode::isStopcodon(dna+i)) {
            for (frame=0; frame<=2; frame++) {
                ritASS = ritASS_cur;
                while ((i<(*ritASS).first)&&(ritASS!=exonASS.rend())){
                    ritASS++;
                }
                ritASS_cur = ritASS;
                while ((i-(*ritASS).first <= Constant::max_exon_len) && (ritASS!=exonASS.rend())) {
                    if ((i-(*ritASS).first>=3) && ((i-(*ritASS).first + 1)%3==frame) && ((*ritASS).first>=orf.leftmostExonBegin(0,i,true))) {
                        ExonCandidate *ec = new ExonCandidate;
                        ec->begin = (*ritASS).first + 2;
                        ec->end = i + 2;
                        ec->setUpScore((*ritASS).second);
                        ec->type = terminal_exon;
			if(ec->len() >= minLen){
			    int_fast64_t key = ec->getKey();
			    ecit = ecs.find(key);
			    if (ecit == ecs.end()){ // insert new EC                                                                                                                                
				ecs.insert(pair<int_fast64_t, ExonCandidate*>(key,ec));
				addECs.insert(pair<int_fast64_t, ExonCandidate*>(key,ec));
			    }
			    else{
				ecit->second->setUpScore(ec->getUpScore());
				delete ec;
			    }
			}
			else{
			    delete ec;
			}
                    }
                    ritASS++;
                }
            }
        }
	
        // computing single genes on the reverse strand with at least Constant::min_coding_len
        if (onRStart(dna+i)) {
            ritRCStop=ritRCStop_cur;
            while ((i < *ritRCStop) && (ritRCStop != exonRCStop.rend())){
                ritRCStop++;
            }
            ritRCStop_cur = ritRCStop;
            while ((i-*ritRCStop <=  Constant::max_exon_len) && (ritRCStop!=exonRCStop.rend())) {
                if ((i-*ritRCStop)%3 == 0) {
                    if ((i-*ritRCStop) >= Constant::min_coding_len) {
                        ExonCandidate *ec = new ExonCandidate;
                        ec->begin = *ritRCStop;
                        ec->end = i + 2;
                        ec->type=rsingleGene;
			if(ec->len() >= minLen){
			    int_fast64_t key = ec->getKey();
			    ecit = ecs.find(key);
			    if (ecit == ecs.end()){ // insert new EC                                                                                                                                
				ecs.insert(pair<int_fast64_t, ExonCandidate*>(key,ec));
				addECs.insert(pair<int_fast64_t, ExonCandidate*>(key,ec));
			    }
			    else{
				delete ec;
			    }
			}
			else{
			    delete ec;
			}
                        break;
                    } else {
                        break;
                    }
                } else {
                    ritRCStop++;
                }
            }
        }

        // computing initials on the reverse strand with at least start codon plus base
        if (withSplicing && onRStart(dna+i)) {
            for (frame=0; frame<=2; frame++) {
                ritRDSS=ritRDSS_cur;
                while ((i<(*ritRDSS).first)&&(ritRDSS!=exonRDSS.rend())){
                    ritRDSS++;
                }
                ritRDSS_cur=ritRDSS;
                int lmb = orf.leftmostExonBegin(2,i,false);
                while((lmb<=(*ritRDSS).first+2)&&(i-(*ritRDSS).first<= Constant::max_exon_len)&&(ritRDSS!=exonRDSS.rend())) {
                    if ((i-(*ritRDSS).first>=2)&&((i+1-(*ritRDSS).first)%3==frame)) {
                        ExonCandidate *ec = new ExonCandidate;
                        ec->begin=(*ritRDSS).first + 2;
                        ec->end=i + 2;
                        ec->setDownScore((*ritRDSS).second);
                        ec->type=rinitial_exon;
			int_fast64_t key = ec->getKey();
			ecit = ecs.find(key);
			if (ecit == ecs.end()){ // insert new EC                                                                                                                                
                            ecs.insert(pair<int_fast64_t, ExonCandidate*>(key,ec));
                            addECs.insert(pair<int_fast64_t, ExonCandidate*>(key,ec));
                        }
			else{
			    ecit->second->setDownScore(ec->getDownScore());
                            delete ec;
                        }
                    }
                    ritRDSS++;
                }
            }
        }
	
        // computing internals on the reverse strand with at least a codon
        if (withSplicing && onRASS(dna+i) && (i + Constant::ass_upwindow_size + Constant::ass_whole_size() - Constant::ass_start < n)) {
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
                        ExonCandidate *ec = new ExonCandidate;
                        ec->begin = (*ritRDSS).first + 2;
                        ec->end = i - 1;
                        ec->setDownScore((*ritRDSS).second);
			ec->setUpScore(computeSpliceSiteScore(p, assminprob, assmaxprob));
                        if (frame==0) {
                            ec->type=rinternal_0;
                        } else if (frame==1) {
                            ec->type=rinternal_1;
                        } else {
                            ec->type=rinternal_2;
                        }
			if(ec->len() >= minLen){
			    int_fast64_t key = ec->getKey();
			    ecit = ecs.find(key);
			    if (ecit == ecs.end()){ // insert new EC                                                                                                                                
				ecs.insert(pair<int_fast64_t, ExonCandidate*>(key,ec));
				addECs.insert(pair<int_fast64_t, ExonCandidate*>(key,ec));
			    }
			    else{
				ecit->second->setDownScore(ec->getDownScore());
				ecit->second->setUpScore(ec->getUpScore());
				delete ec;
			    }
			}
			else{
			    delete ec;
			}
                    }
                    ritRDSS++;
                }
            }
        }

        // computing terminals on the reverse strand with at least one base plus stopcodon
        if (withSplicing && onRASS(dna+i) && (i + Constant::ass_upwindow_size + Constant::ass_whole_size() - Constant::ass_start < n)) {
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
                        ExonCandidate *ec = new ExonCandidate;
                        ec->begin = *ritRCStop;
                        ec->end = i - 1;
                        ec->setUpScore(computeSpliceSiteScore(p, assminprob, assmaxprob));
                        if (frame==0) {
                            ec->type=rterminal_2;
                        } else if (frame==1) {
                            ec->type=rterminal_1;
                        } else {
                            ec->type=rterminal_0;
                        }
			if(ec->len() >= minLen){
			    int_fast64_t key = ec->getKey();
			    ecit = ecs.find(key);
			    if (ecit == ecs.end()){ // insert new EC                                                                                                                                
				ecs.insert(pair<int_fast64_t, ExonCandidate*>(key,ec));
				addECs.insert(pair<int_fast64_t, ExonCandidate*>(key,ec));
			    }
			    else{
				ecit->second->setUpScore(ec->getUpScore());
				delete ec;
			    }
			} 
			else{
			    delete ec;
			}
                        break;
                    } else {
                        ritRCStop++;
                    }
                }
            }
        }
    }
}

ExonCandidate* create(int_fast64_t key, const char* dna, int dnalen){

    int_fast64_t chrStart = (key >> 22);
    int_fast64_t chrEnd = ((key>>7) & 0x7FFF) + chrStart; // 0x7FFF bit mask for retrieving bits 1-15 
    int len = chrEnd - chrStart + 1;
    ExonType type = (ExonType)((key>>2) & 0x1F); // 0x7FFF bit mask for retrieving bits 1-5
    int lenMod3 = (key & 0x3); // bit mask for retrieving bits 1-2
    
    // check length
    if( ((len) % 3) != lenMod3)
        return NULL;
    if(len > Constant::max_exon_len)
    	return NULL;
    if((type == singleGene || type == rsingleGene) && len < Constant::min_coding_len)
    	return NULL;
    
    ExonCandidate *ec = new ExonCandidate;
    ec->begin=chrStart;
    ec->end=chrEnd;
    ec->type=type;

    // check type (TODO: allow non-canonical, hinted SS)
    if(!ec->correctType(dna,dnalen)){
	delete ec;
	return NULL;
    }

    // check for inFrameStop
    if( GeneticCode::containsInFrameStopcodon(dna, ec->getFirstCodingBase(), ec->getLastCodingBase(), isPlusExon(type), 0) ){
        delete ec;
        return NULL;
    }
    return ec;
}

// verify ExonType on sequence (f.e. if type is initial , check whether
//  left boundary is on a start codon and right boundary is on a DSS)
bool ExonCandidate::correctType(const char* dna, int dnalen){

    bool leftTruncated = (begin == 0);
    bool rightTruncated = (end == dnalen-1);
    if(leftTruncated || rightTruncated)
	cout << "liftover of truncated exon: " << *this << endl;

    bool correctType = false;

    switch(type){
    case singleGene:
        if(onStart(dna+begin) && GeneticCode::isStopcodon(dna+end-2))
            correctType=true;
        break;
    case initial_0: case initial_1: case initial_2:
        if(onStart(dna+begin) && (rightTruncated || onDSS(dna+end+1) || onAltDSS(dna+end+1)))
            correctType=true;
        break;
    case internal_0: case internal_1: case internal_2:
        if((leftTruncated || onASS(dna+begin-2)) && (rightTruncated || onDSS(dna+end+1) || onAltDSS(dna+end+1)))
            correctType=true;
        break;
    case terminal_exon:
	if((leftTruncated || onASS(dna+begin-2)) && GeneticCode::isStopcodon(dna+end-2))
            correctType=true;
        break;
    case rsingleGene:
        if(onRStart(dna+end-2) && GeneticCode::isRCStopcodon(dna+begin))
            correctType=true;
        break;
    case rinitial_exon:
        if(onRStart(dna+end-2) && (leftTruncated || onRDSS(dna+begin-2) || onAltRDSS(dna+begin-2)))
            correctType=true;
        break;
    case rinternal_0: case rinternal_1: case rinternal_2:
        if((rightTruncated || onRASS(dna+end+1)) && (leftTruncated || onRDSS(dna+begin-2) || onAltRDSS(dna+begin-2)))
            correctType=true;
        break;
    case  rterminal_0: case rterminal_1: case rterminal_2:
        if((rightTruncated || onRASS(dna+end+1)) && GeneticCode::isRCStopcodon(dna+begin))
            correctType=true;
    default :;
    }
    return correctType;
}

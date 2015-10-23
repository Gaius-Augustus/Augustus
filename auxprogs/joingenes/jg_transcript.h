#ifndef _JG_TRANSCRIPT_HPP
#define _JG_TRANSCRIPT_HPP

#include <stdlib.h>
#include <string.h>
#include <cstring>
#include <iostream>
#include <list>
#include <unordered_map>
#include <algorithm>
#include <fstream>

using namespace std;

class Gene;
class Transcript;

#define NUM_TYPES 4
enum boundaryStatusType{TYPE_UNKNOWN = -1,		// predRange-boundaryStatusType
    NO_PROB=0, FREED=1, INDIRECT_PROBLEM=2, TOO_CLOSE=3		// INDIRECT_PROBLEM means that it is not too close, but there is an other transcript, which transcend its boundary
								// FREED means that the problem exon was erased
};

enum compareType{UNEQUAL=0, ALTERNATIVE=1, EQUAL=2
};

bool overlappingCdsWithAnything(Transcript* t1, Transcript* t2);

struct Properties{
public:
    int errordistance;
    string genemodel; 
    bool onlyCompare;
    list<string> filenames;
    list<int> priorities;
    list<int> supprList;
    string outFileName;
    bool join;
    bool selecting;
    int nrOfPrintedGenes;
    unsigned int unknownCount;
    unordered_map<string,Gene*>* geneMap;
    unordered_map<string,Transcript*>* transcriptMap;
    unsigned int minimumIntronLength;
    bool alternatives;
    bool stopincoding;
    unordered_map<string,unsigned int> warningCount;
};

class Exon{
 public:
    string chr;
    int from, to;
    string feature;
    float score;
    int frame;
    int rangeToBoundary;
    boundaryStatusType boundaryProblem;
    list <Transcript*> indirectBoundProbEnemies;

    int tooMany;
    int tooFew;
    int penalty;
    int distance;
    pair<int,int> predRange;
    bool isUtr () { return (frame == -1); }
    bool operator<(Exon const& rhs) const {
	if (from != rhs.from)
	    return (from < rhs.from);
	else
	    return (to < rhs.to);
    }
};

void output_exon_list(Transcript const* tx);							// only for semantic tests

class Gene{
 public:
    Gene(string geneID=""){
	    this->nrOfTx = 1;
	    this->nrOfPrintedTx = 0;
	    this->g_id = geneID;
	    this->printed = false;
    }
    string g_id;
    list<Transcript*> children;
    int nrOfTx;
    int nrOfPrintedTx;
    bool printed;
};

class Transcript{
 public:
    Transcript(string txID = "", Gene* gene=NULL){
	this->parent = gene;
	this->t_id = txID;
	this->tx_complete.first = false;
	this->tx_complete.second = false;
	this->tl_complete.first = false;
	this->tl_complete.second = false;
	this->isNotFrameCorrect = false;
	this->tss = -1;
	this->tts = -1;
    }

    list<string> supporter;		// at the moment these list only consist the pointer to transcripts, who are totaly equal to this transcript (doesnt matter whether their creation is related)
    list<Exon> exon_list;				// saves every feature (cds, utr, exon [cds+utr]
    list<Exon> intron_list;
    list<Exon> stop_list;
    Gene* parent;
    string source;
    string t_id;
    char strand;
    int tis;							// position of first base in translation sequence (translation initation site)
    int tes;							// position of last base in translation sequence (translation end site)
    int tss;
    int tts;
    int priority;
    pair<int,int> pred_range;			// prediction range borders: <low,high> sequence position
    pair<bool,bool> tx_complete;		// transcript end complete at <lower,higher> base number
    pair<bool,bool> tl_complete;		// translation end complete at <start,stop> codon
    pair<bool,bool> separated_codon;		// separated <start,stop> codon

    pair<string,string> joinpartner;	// pointer to transcript that was used to complete the <start,stop> codon side
    pair<string,string> utr_joinpartner;	// pointer to transcript that was used to fulfill the <downstream-UTR,upstream-UTR> codon side
    bool must_be_deleted;				// a transcript with this flag must be deleted
    pair<Exon*,Exon*> outerCds;				// CDS with the <lowest,highest> sequence position
    int qualityScore;
    float predictionScore;
    list<string> consistent;
    bool isNotFrameCorrect;
    string originalId;
    string inputFile;

    double penalty;

    compareType compareValue;

    boundaryStatusType startBoundaryType(){
	if (strand == '+'){
	    return exon_list.front().boundaryProblem;
	}else{
	    return exon_list.back().boundaryProblem;
	}
    }

    boundaryStatusType stopBoundaryType(){
	if (strand == '+'){
	    return exon_list.back().boundaryProblem;
	}else{
	    return exon_list.front().boundaryProblem;
	}
    }

    bool hasCommonExon(Transcript* tx){
	list<Exon>::const_iterator it1 = exon_list.begin();
	list<Exon>::const_iterator it2 = tx->exon_list.begin();

	while (it1 != exon_list.end() && it2 != tx->exon_list.end()){
	    if ((*it1).from == (*it2).from){
		if((*it1).to == (*it2).to){
		    if ((*it1).frame == (*it2).frame){
			return true;
		    }else{
			it1++; it2++;
		    }
		}else if ((*it1).to < (*it2).to){
		    it1++;
		}else{
		    it2++;
		}
	    }else if ((*it1).from > (*it2).from){
		it2++;
	    }else{
		it1++;
	    }
	}
	return false;
    }

    bool hasCommonTlStart(Transcript* tx){
	if (!this->tl_complete.first || !tx->tl_complete.first){return false;}
	if (tis == tx->tis && strand == tx->strand){return true;}else{return false;}
    }

    bool hasCommonTlStop(Transcript* tx){
	if (!this->tl_complete.second || !tx->tl_complete.second){return false;}
	if (tes == tx->tes && strand == tx->strand){return true;}{return false;}
    }

    string getChr(){
        return exon_list.front().chr;
    }

    int getTxStart(){
        return exon_list.front().from;
    }
    int getTxEnd(){
        return exon_list.back().to;
    }

    Exon* getCdsFront(){
	for (list<Exon>::iterator it = exon_list.begin(); it != exon_list.end(); it++){
	    if ((*it).feature == "CDS"){
		return &(*it);	// (*it).from
	    }
	}
	return NULL;
    }

    Exon* getCdsBack(){
	for (list<Exon>::reverse_iterator it = exon_list.rbegin(); it != exon_list.rend(); it++){
	    if ((*it).feature == "CDS"){
		return &(*it);	// (*it).to
	    }
	}
	return NULL;
    }

    void initiate(Properties &properties){
        exon_list.sort();				// essential for the other steps, so it should be the first one
	stop_list.sort();

	outerCds.first = getCdsFront();
	outerCds.second = getCdsBack();

        exontoutr();

	if (strand == '+'){
	    tis = (outerCds.first)->from;
	    tes = (outerCds.second)->to;
	}else{
	    tis = (outerCds.second)->to;
	    tes = (outerCds.first)->from;
	}
	if (properties.stopincoding){
            tes_to_cds();					// needs utr structure!
	}

	calcRangeToBoundary(properties);
    }

    void calcRangeToBoundary(Properties &properties){
	for (list<Exon>::iterator it = exon_list.begin(); it != exon_list.end(); it++){
	    (*it).rangeToBoundary = -1;
	}
//	exon_list.front().rangeToBoundary = -1;
//	exon_list.back().rangeToBoundary = -1;
	if (pred_range.second){
	    exon_list.back().rangeToBoundary = pred_range.second - getTxEnd();

	    if (exon_list.back().rangeToBoundary > properties.errordistance){
		//boundaryProblem.second = NO_PROB;
		exon_list.back().boundaryProblem = NO_PROB;
	    }else{
		//boundaryProblem.second = TOO_CLOSE;
		exon_list.back().boundaryProblem = TOO_CLOSE;
	    }
	}
	if (pred_range.first){
	    exon_list.front().rangeToBoundary = getTxStart() - pred_range.first;

	    if (exon_list.front().rangeToBoundary > properties.errordistance){
		//boundaryProblem.first = NO_PROB;
		exon_list.front().boundaryProblem = NO_PROB;
	    }else{
		//boundaryProblem.first = TOO_CLOSE;
		exon_list.front().boundaryProblem = TOO_CLOSE;
	    }
	}
    }

    void tes_to_cds(){
	/*int cdsfront;
	if (outerCds.first)
	    cdsfront = (outerCds.first)->from;
	else
	    cdsfront = -1;
	int cdsback;
	if (outerCds.second)
	    cdsback = (outerCds.second)->to;
	else
	    cdsback = -1;*/
	bool done = false;
//if (t_id == "g10104.t1"){cout << "HIER" << endl; sleep(1);}
	if (strand == '+'){
	    if (tl_complete.second){
		if (tes != stop_list.back().to){
		    for (list<Exon>::iterator it = exon_list.begin(); it != exon_list.end(); it++){
			if (((*it).feature == "UTR" || (*it).feature == "3'-UTR" || (*it).feature == "5'-UTR") && tes < (*it).from){
			    done = true;
			    if ((*it).from == tes + 1){
				int temp = (*it).to - (*it).from + 1;
				if (temp > 3){
				    (*it).from += 3;
				    it--;
				    (*it).to += 3;
				}else if (temp == 3){
				    it = exon_list.erase(it);
				    it--;
				    (*it).to += 3;
				}else{
				    it = exon_list.erase(it);
				    it--;
				    (*it).to += temp;
				    it++;
				    if ((*it).to - (*it).from + 1 > 3 - temp){
					Exon new_cds = (*it);
					new_cds.from = (*it).from;
					new_cds.to = (*it).from + (3 - (temp + 1));
					new_cds.feature = "CDS";
					new_cds.frame = (new_cds.to - new_cds.from + 1) % 3;
					(*it).from += 3 - temp;
					exon_list.insert(it, new_cds);
				    }else if ((*it).to - (*it).from + 1 == 3 - temp){
					(*it).feature = "CDS";
					(*it).frame = 0;
				    }else{
					cerr << "Fatal error: An exon with start/stop codon part is not even 3 bases long." << endl;
					exit( EXIT_FAILURE );
				    }
				}
			    }else{
				Exon new_cds = (*it);
				new_cds.from = (*it).from;
				new_cds.to = (*it).from + 2;
				new_cds.feature = "CDS";
				new_cds.frame = 0;
				(*it).from += 3;
				exon_list.insert(it, new_cds);
			    }
			    break;
			}
		    }
		    if (!done){
			if (tes == stop_list.back().to - 3){
			    (outerCds.second)->to += 3;
			}else{
			    Exon new_cds = *(outerCds.second);
			    new_cds.from = stop_list.back().to - 2;
			    new_cds.to = stop_list.back().to;
			    new_cds.feature = "CDS";
			    new_cds.frame = 0;
			    exon_list.push_back(new_cds);
			}
		    }
		    tes = stop_list.back().to;
		}
	    }/*else{
		tes = cdsback;
	    }
	    if (!tl_complete.first){
		tis = cdsfront;
	    }*/
	}else{	// strand == "-"
	    if (tl_complete.second){	// stop codon exists
//if (t_id == "g10104.t1"){cout << "HIER2" << endl; sleep(1);}
		if (tes != stop_list.front().from){		// stop codon is not in cds
//if (t_id == "g10104.t1"){cout << "HIER4" << endl; sleep(1);}
		    for (list<Exon>::iterator it = exon_list.begin(); it != exon_list.end(); it++){
			if ((*it).feature == "CDS" && tes == (*it).from){
			    if (it == exon_list.end()){
				cerr << "Fatal error: No exon for an existing codon." << endl;
				exit( EXIT_FAILURE );
			    }
			    //it--;
			    if (it != exon_list.begin() && ((*(it--)).feature == "UTR" || (*it).feature == "3'-UTR" || (*it).feature == "5'-UTR")){
				done = true;
				if ((*it).to == tes - 1){
				    int temp = (*it).to - (*it).from + 1;
				    if (temp > 3){
					(*it).to -= 3;
					it++;
					(*it).from -= 3;
				    }else if (temp == 3){
					it = exon_list.erase(it);
					(*it).from -= 3;
				    }else{
					it = exon_list.erase(it);
					(*it).from -= temp;
					it--;
					if ((*it).to - (*it).from + 1 > 3 - temp){
					    Exon new_cds = (*it);
					    new_cds.to = (*it).to;
					    new_cds.from = (*it).to - (3 - (temp + 1));
					    new_cds.feature = "CDS";
					    new_cds.frame = (new_cds.to - new_cds.from + 1) % 3;
					    (*it).to -= 3 - temp;
					    it++;
					    exon_list.insert(it, new_cds);
					}else if ((*it).to - (*it).from + 1 == 3 - temp){
					    (*it).feature = "CDS";
					    (*it).frame = 0;
					}else{
					    cerr << "Fatal error: An exon with start/stop codon part is not even 3 bases long." << endl;
					    exit( EXIT_FAILURE );
					}
				    }
				}else{
				    Exon new_cds = (*it);
				    new_cds.to = (*it).to;
				    new_cds.from = (*it).to - 2;
				    new_cds.feature = "CDS";
				    new_cds.frame = 0;
				    (*it).from -= 3;
				    it++;
				    exon_list.insert(it, new_cds);
				}
			    }
			    break;
			}
		    }
		    if (!done){
			if (tes == stop_list.front().from + 3){
			    (outerCds.first)->from -= 3;
			}else{
			    Exon new_cds = *(outerCds.first);
			    new_cds.from = stop_list.front().from;
			    new_cds.to = stop_list.front().from + 2;
			    new_cds.feature = "CDS";
			    new_cds.frame = 0;
			    exon_list.push_front(new_cds);
			}
		    }
//if (t_id == "g10104.t1"){cout << "HIER3" << endl; sleep(1);}
		    tes = stop_list.front().from;
		}
	    }/*else{
		tes = cdsfront;
	    }
	    if (!tl_complete.first){
		tis = cdsback;
	    }*/
	}
    }

    void exontoutr(){
	list<Exon>::iterator it_prev;
	list<Exon>::iterator it_next;
	for (list<Exon>::iterator it = exon_list.begin(); it != exon_list.end(); it++){
	    if ((*it).feature == "exon"){
		if (it != exon_list.begin()){
		    it_prev = it;
		    it_prev--;
		    if ((*it_prev).feature == "CDS"){
			if ((*it_prev).to >= (*it).from){
			    if ((*it_prev).from == (*it).from && (*it_prev).to == (*it).to){
				it = exon_list.erase(it);
				it--;
				continue;
			    }else if ((*it_prev).from == (*it).from){
				(*it).from = (*it_prev).to + 1;
				(*it).feature = "UTR";
				continue;
			    }else{cerr << "unexpected case in exontoutr() nr 1 with: " << originalId << " from " << inputFile << " because there is an exon from " << (*it).from << " to " << (*it).to << " and a CDS from " << (*it_prev).from << " to " << (*it_prev).to << endl;}
			}
		    }
		}
		it_next = it;
		it_next++;
		if (it_next != exon_list.end()){
		    if ((*it_next).feature == "CDS"){
			if ((*it).to >= (*it_next).from){
			    if ((*it_next).from == (*it).from && (*it_next).to == (*it).to){
				it = exon_list.erase(it);
				it--;
				continue;
			    }else if ((*it_next).to == (*it).to){
				(*it).to = (*it_next).from - 1;
				(*it).feature = "UTR";
				continue;
			    }else if ((*it).to > (*it_next).to){
				(*it).feature = "UTR";
				Exon new_utr = (*it);
				new_utr.from = (*it_next).to + 1;
				(*it).to = (*it_next).from - 1;
				it_next++;
				exon_list.insert(it_next, new_utr);
				continue;
			    }else{cerr << "unexpected case in exontoutr() nr 2" << originalId << " from " << inputFile << " because there is an exon from " << (*it).from << " to " << (*it).to << " and a CDS from " << (*it_next).from << " to " << (*it_next).to << endl;}
			}
		    }
		}
		(*it).feature = "UTR";
	    }
	}
    }
    bool operator<(Transcript const& rhs) const {
	if (min(tis,tes) != min(rhs.tis,rhs.tes))
	    return (min(tis,tes) < min(rhs.tis,rhs.tes));
	else
	    return (max(tis,tes) < max(rhs.tis,rhs.tes));
    }

    bool indirectBoundaryProblem(Transcript* tx){

	if (overlappingCdsWithAnything(this, tx) /*&& this->strand == tx->strand*/){
	    // does *it transcend the predictionBoundary of *itInside in tail direction
	    if (this->pred_range.second && tx->getTxEnd() > this->pred_range.second){
		return true;
	    }

	    // does *it transcend the predictionBoundary of *itInside in front direction
	    if (this->pred_range.first && tx->getTxStart() < this->pred_range.first){
		return true; // tx->exon_list.front().boundaryProblem = INDIRECT_PROBLEM;
	    }
	}
/*cout << "INSIDE: " << this->exon_list.front().indirectBoundProbEnemies.size() << " " << this->exon_list.front().indirectBoundProbEnemies.size() << endl;
	for(list<Transcript*>::iterator it = this->exon_list.front().indirectBoundProbEnemies.begin(); it != this->exon_list.front().indirectBoundProbEnemies.end(); it++){
	    if ( (*it) == tx ){
		return true;
	    }
	}
	for(list<Transcript*>::iterator it = this->exon_list.back().indirectBoundProbEnemies.begin(); it != this->exon_list.back().indirectBoundProbEnemies.end(); it++){
cout << (*it) << " " << tx << " " << (*it)->t_id << " " << tx->t_id << " " << (*it)->originalId << " " << tx->originalId << endl;
	    if ( (*it)->t_id == tx->t_id ){
		return true;
	    }
	}*/
	return false;
    }
};

class Point{
 public:
    double i_x;
    double i_y;
    int x;
    int y;

    bool operator<(Point const& rhs) const {
	if (y != rhs.y)
	    return (y > rhs.y);
	else
	    return (x < rhs.x);
    }
};

void divideInOverlapsAndConquer(list<Transcript*> &transcript_list, Properties &properties);
void workAtOverlap(list<Transcript*> &overlap, Properties &properties);
void selection(list<Transcript*> &overlap, Properties &properties);
void bactSelection(list<Transcript*> &overlap, Properties &properties);
void eukaSelection(list<Transcript*> &overlap, Properties &properties);
bool areOverlapping(Transcript* t1, Transcript* t2);
void joinCall(list<Transcript*> &overlap, Properties &properties);
void compareAndSplit(list<Transcript*> &overlap, Properties &properties);
double simpleProkScore(Transcript const* tx);
bool areSimilar(Transcript const* t1, Transcript const* t2);
void tooCloseToBoundary(list<Transcript*> &overlap, Properties &properties);
void search_n_destroy_doublings(list<Transcript*> &overlap, Properties &properties, bool ab_initio);
void search_n_destroy_parts(list<Transcript*> &overlap, Properties &properties);
bool compare_transcripts(Transcript const* t1, Transcript const* t2);
pair<bool,bool> is_part_of(Transcript const* t1, Transcript const* t2);
bool compare_priority(Transcript const* lhs, Transcript const* rhs);
void join(list<Transcript*> &overlap, char side, Properties &properties);
void joining(Transcript* t2, char strand, Transcript* txNew, int fittingCase, Properties &properties);
int is_combinable(Transcript const* t1, Transcript const* t2, char strand, bool frontSide, Properties &properties);
bool check_frame_annotation(Transcript const &transcript);
void eval_gtf(list<Transcript*> &overlap, Properties &properties);
bool strandeq(Exon ex1, Exon ex2, char strand);

void weight_info(list<Transcript*> &overlap);

bool overlapping(Transcript* t1, Transcript* t2);
// see top:  bool overlappingCdsWithAnything(Transcript* t1, Transcript* t2);
bool overlappingCdsWithCds(Transcript* t1, Transcript* t2);
bool overlappingCdsOnlyWithUtr(Transcript* t1, Transcript* t2);
bool overlappingUtrOnly(Transcript* t1, Transcript* t2);
int isCombinable(Transcript* t1, Transcript* t2, bool frontSide, Properties &properties);

bool alternativeVariants(Transcript* t1, Transcript* t2);
void eukaSelectionDevelopment(list<Transcript*> &overlap, Properties &properties);
void recursiv(list<list<Transcript*>> &groups, list<Transcript*> actGroup, list<string> openTx, Properties &properties, unordered_map<string,Transcript*> &txMap, list<string> selectedTx);
list<Transcript*> calcGeneStructur(list<Transcript*> overlap, Properties &properties);

void deleteTx(Transcript* tx, Properties &properties);
bool compare_quality(Transcript const* lhs, Transcript const* rhs);

void testInputAlternatives(Properties &properties);

void displayWarning(string const &warning, Properties &properties, string warningString);
void warningSummary(string const &warning, string const &warning2, Properties &properties, string warningString);

#endif

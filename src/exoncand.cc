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
#include <fstream>
#include <iostream>
#include <string>
#include <sys/time.h>

long int ExonCandidate::getStart() {
	return begin;
}

long int ExonCandidate::getEnd() {
	return end;
}

ExonType ExonCandidate::getExonType() {
	return type;
}

StateType ExonCandidate::getStateType(){

  if(type >= singleGene && type <= terminal_exon)
    return StateType (type + 1);

  if(type >= rsingleGene && type <= rterminal_2)
    return StateType (type + 28);

  return TYPE_UNKNOWN;
}

double ExonCandidate::getScore() {
	return score;
}

string ExonCandidate::createKey() {
	if(this == NULL) {
		return "no key";
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

char* randomDNAGenerator (double gc_content) {
	const char *gc = "cg"; // 0->c, 1->g
	const char *at = "at"; // 0->a, 1->t
	int dnaLength = 10;
	//int number = 5;
	double randomValue;
	char *sequence;

	//for (int i=0; i<number; i++) {
		sequence = new char[dnaLength+1];
		for (int j=0; j<dnaLength; j++) {
			randomValue=rand()/ (1.0 + RAND_MAX);
			if (randomValue<=gc_content) {
				sequence[j] = gc[rand() % 2]; // Do they have the same probability? G is part of all splice sites!
			} else {
				sequence[j] = at[rand() % 2];
			}
			//cout<<sequence[j];
		}
		//cout<<endl;
	//}
	return sequence;
}

// maybe I do not need it anymore
/*bool compType (ExonCandidate* a, ExonCandidate* b) {
	if (a->type != b->type) {
		return (a->type<b->type);
	} else{
		return (a->end>b->end);
	}
}*/

bool comp (ExonCandidate* a, ExonCandidate* b) {
	if (a->type != b->type) {
		return (a->type<b->type);
	} else  if (a->end != b->end) {
		return (a->end<b->end);
	} else {
		return (a->begin<b->begin);
	}
}

list<ExonCandidate* >  getExonCands(const char *dna, float assqthresh, float dssqthresh){
	int n = strlen(dna);
	int max_exon_length = 12000;
	int frame;
	/*int single_fw=0;
  	  int single_rv=0;
  	  int initial_fw=0;
  	  int initial_rv=0;
  	  int internal_fw=0;
  	  int internal_rv=0;
  	  int terminal_fw=0;
  	  int terminal_rv=0;*/
	list<int> exonStart;
	list<int> exonRCStop;
	list<int> exonASS;
	list<int> exonRDSS;
	ExonCandidate *ptr;
	list<ExonCandidate*> candidates;

	OpenReadingFrame orf(dna, max_exon_length, n);
  
	for (int i=0; i<=n-1; i++) {
	// positions of all startcodons "atg"
		if (onStart(dna+i)) {
	    	exonStart.push_back(i+1);
		}
	// positons of all ASSs "ag"
		if (onASS(dna+i)) {
			exonASS.push_back(i);
		}
    // positions of all reverse DSS "ac"
		if (onRDSS(dna+i)) {
			exonRDSS.push_back(i);
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

	for (int i=n-1; i>=2; i--) {
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
					ptr->begin=*ritStart;
					ptr->end=i+3;
					ptr->type=singleGene;
					candidates.push_back(ptr);
					//single_fw++;
				}
				ritStart++;
			};
		}

	// computing initials on the forward strand with at least startcodon plus base
		if (onDSS(dna+i)) {
			for (frame=0; frame<=2; frame++) {
				ritStart=ritStart_cur;
				while ((i<*ritStart)&&(ritStart!=exonStart.rend())){
					ritStart++;
				}
				ritStart_cur=ritStart;
				int lmb = orf.leftmostExonBegin(frame,i,true);
				while((lmb<=*ritStart)&&(i-*ritStart<=max_exon_length)&&(ritStart!=exonStart.rend())) {
					if ((i-*ritStart>=2)&&((i-*ritStart+1)%3==frame)) {
						ptr = new ExonCandidate;
						ptr->begin=*ritStart;
						ptr->end=i;
						if (frame==0) {
							ptr->type=initial_0;
						} else if (frame==1) {
							ptr->type=initial_1;
						} else {
							ptr->type=initial_2;
						}
						candidates.push_back(ptr);
						//initial_fw++;
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
					if (i-*ritASS>=5) {
						ptr = new ExonCandidate;
						ptr->begin=*ritASS+3;
						ptr->end=i;
						if (frame==0) {
							ptr->type=internal_0;
						} else if (frame==1) {
							ptr->type=internal_1;
						} else {
							ptr->type=internal_2;
						}
						//cout<<ptr->createKey()<<endl;
						candidates.push_back(ptr);
						//internal_fw++;
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
						ptr->begin=*ritASS+3;
						ptr->end=i+3;
						ptr->type=terminal_exon;
						candidates.push_back(ptr);
						//terminal_fw++;
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
					ptr->begin=*ritRCStop+1;
					ptr->end=i+3;
					ptr->type=rsingleGene;
					candidates.push_back(ptr);
					//single_rv++;
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
						ptr->begin=*ritRDSS+3;
						ptr->end=i+3;
						ptr->type=rinitial_exon;
						candidates.push_back(ptr);
						//initial_rv++;
					}
					ritRDSS++;
				};
			}
		}

	// computing internals on the reverse strand with at least a codon
		if (onRASS(dna+i)) {
			for (frame=0; frame<=2; frame++) {
				ritRDSS=ritRDSS_cur;
				while ((i<*ritRDSS)&&(ritRDSS!=exonRDSS.rend())){
					ritRDSS++;
				}
				ritRDSS_cur=ritRDSS;
				int lmb = orf.leftmostExonBegin(frame,i,false);
				while((lmb<=*ritRDSS)&&(i-*ritRDSS<=max_exon_length)&&(ritRDSS!=exonRDSS.rend())) {
					if (i-*ritRDSS>=5) {
						ptr = new ExonCandidate;
						ptr->begin=*ritRDSS+3;
						ptr->end=i;
						if (frame==0) {
							ptr->type=rinternal_0;
						} else if (frame==1) {
							ptr->type=rinternal_1;
						} else {
							ptr->type=rinternal_2;
						}
						candidates.push_back(ptr);
						//internal_rv++;
					}
					ritRDSS++;
				};
			}
		}

	// computing terminals on the reverse strand with at least one base plus stopcodon
		if (onRASS(dna+i)) {
			for (frame=0; frame<=2; frame++) {
		    	ritRCStop=ritRCStop_cur;
		    	while ((i<*ritRCStop)&&(ritRCStop!=exonRCStop.rend())){
		    		ritRCStop++;
		    	}
		    	ritRCStop_cur=ritRCStop;
		    	while ((i-*ritRCStop<=max_exon_length)&&(ritRCStop!=exonRCStop.rend())) {
		    		if ((i-*ritRCStop>=4)&&((i-*ritRCStop)%3==frame)) {
		    			ptr = new ExonCandidate;
		    			ptr->begin=*ritRCStop+1;
		    			ptr->end=i;
		    			if (frame==0) {
		    				ptr->type=rterminal_0;
		    			} else if (frame==1) {
		    				ptr->type=rterminal_1;
		    			} else {
		    				ptr->type=rterminal_2;
		    			}
		    			candidates.push_back(ptr);
		    	    	//terminal_rv++;
		    	    	break;
		    		} else {
		    			 ritRCStop++;
		    		}
		    	};
			}
		}
	}
	// to use computeIntervals correctly, you have to arrange the candidates according to statetype
	candidates.sort(comp);
	if (!candidates.empty()) {
		for (list<ExonCandidate*>::iterator lit=candidates.begin(); lit!=candidates.end(); lit++) {
			cout<<"start "<<(*lit)->begin;
			cout<<"     end "<<(*lit)->end;
			cout<<"     typ "<<(*lit)->type<<endl;
		}
	}

	// int anzahl=single_fw+single_rv+initial_fw+initial_rv+internal_fw+internal_rv+terminal_fw+terminal_rv;

	/*cout<<"Single Gene forward: "<<single_fw++<<endl;
	cout<<"Single Gene reverse: "<<single_rv++<<endl;
	cout<<"Initial forward: "<<initial_fw<<endl;
	cout<<"Initial reverse: "<<initial_rv<<endl;
	cout<<"Internal forward: "<<internal_fw<<endl;
	cout<<"Internal reverse: "<<internal_rv<<endl;
	cout<<"Terminal forward: "<<terminal_fw<<endl;
	cout<<"Terminal reverse: "<<terminal_rv<<endl;*/
	//cout<< "strand length: "<<n<<endl;

	// assqthresh, dssqthresh
	Double assminprob = IntronModel::assBinProbs.getMinProb(assqthresh);
	Double dssminprob = IntronModel::dssBinProbs.getMinProb(dssqthresh);
	cout << "thresholds: dssminprob=" << dssminprob << " assminprob=" << assminprob << endl;
	// demo of splice site scores
	Double p;
	for (int i=0; i<200; i++) { // limit for testing and demo of code
		if (i + Constant::ass_whole_size() + Constant::ass_upwindow_size < n){ // see class IntronModel for boundaries
			p = IntronModel::aSSProb(i, true); // second argument: plusstrand?
			if (p>0) {
			//	cout <<  ((p >= assminprob)? "strong" : "weak") << " intron may end at pos " << i+Constant::ass_upwindow_size+Constant::ass_start+ASS_MIDDLE-1
					// << " (0-based) on plusstrand, " << endl;
			}
		}
		/*if (i + Constant::ass_whole_size() + Constant::ass_upwindow_size < n){ // see class IntronModel for boundaries
			p = IntronModel::aSSProb(i, false); // second argument: minusstrand?
			if (p>0) {
				cout <<  ((p >= assminprob)? "strong" : "weak") << " intron may end at pos " << i+Constant::ass_upwindow_size+Constant::ass_start+ASS_MIDDLE-1
				<< " (0-based) on minusstrand, " << endl;
			}
		}*/
		if (i + Constant::dss_whole_size() < n){ // see class IntronModel for boundaries
			p = IntronModel::dSSProb(i, true);
		 	if (p>0) {
		 	//	cout << ((p >= dssminprob)? "strong" : "weak") << " intron may start at pos " << i+Constant::dss_start  << " (0-based) on plusstrand, "  << endl;
		 	}
		}
		/*if (i + Constant::dss_whole_size() + Constant::dss_upwindow_size < n){ // see class IntronModel for boundaries
			p = IntronModel::dSSProb(i, false); // second argument: minusstrand?
			if (p>0) {
				cout <<  ((p >= dssminprob)? "strong" : "weak") << " intron may start at pos " << i+Constant::dss_upwindow_size+Constant::dss_start+ASS_MIDDLE-1
				<< " (0-based) on minusstrand, " << endl;
			}
		}*/
	}
	//computeIntervals(candidates, n);
	return candidates;
}

// I probably do not need this anymore
 // computes intervals to the exon candidates vector so that access to the vector later is quicker
/*void computeIndices(list<ExonCandidate*> cand, int seqlen) {
	 firstEnd[i] holds an iterator to cand to the first element f in the vector such that
	 	f->end >= i*block_size, where i>=0 and i ist such that i*block_size <= seqlen
	//list<ExonCandidate*>::iterator **firstEnd;
	 lastBegin[i] holds an iterator to cand to the vector such that all following elements f have
	 	f->begin > i*block_size, where i>=0 and i ist such that i*block_size <= seqlen
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

#include "jg_transcript.h"
#include "jg_ios.h"

using namespace std;

int count1 = 0;				// only for semantic tests
int value = 0;				// only for semantic tests

void divideInOverlapsAndConquer(list<Transcript> &transcript_list, Properties &properties){
    // divide a transcript list in clusters of overlapping transcripts (transitive closure) and work on each cluster separately
    transcript_list.sort();

    list<Transcript> new_transcripts;
    list<Transcript*> overlap;

    /*    string filename2 = "/home/lars/lars/test_data/eval_test.txt";
	  fstream outfile2;
	  outfile2.open(filename2, ios::out);
	  outfile2.close();*/

    int max_base = max(transcript_list.front().tis,transcript_list.front().tes);

    for (list<Transcript>::iterator it = transcript_list.begin(); it != transcript_list.end(); it++){
	//cout << it->exon_list.begin()->chr << endl;
	//sleep(1);
        if (min((*it).tis,(*it).tes) < max_base){
            overlap.push_back(&*it);
            if (max_base < max((*it).tis,(*it).tes)){
                max_base = max((*it).tis,(*it).tes);
            }
        }
        else{
            //eval_gtf(overlap, errordistance);
	    if (!properties.onlyCompare){
		workAtOverlap(overlap, new_transcripts, properties);
		saveOverlap(overlap, properties.outFileName, properties);
	    }else{
		compareAndSplit(overlap, properties);
	    }
            overlap.clear();
            max_base = max((*it).tis,(*it).tes);
            overlap.push_front(&*it);
        }
    }
    workAtOverlap(overlap, new_transcripts, properties);
    saveOverlap(overlap, properties.outFileName, properties);
    overlap.clear();
}

void workAtOverlap(list<Transcript*> &overlap, list<Transcript> &new_transcripts, Properties &properties)
{
    // calls methods for joining transcripts with the target that most of transcripts are complete (have start and stop codon) and delete duplicates and other unwanted transcripts from the overlap

    // search doubling and delete the ones with the lowest priority/score
    search_n_destroy_doublings(overlap, properties.errordistance, true);

    if (properties.join){
	joinCall(overlap, new_transcripts, properties);
    }
    // searchs and destroys transcripts, which are fully in other transcripts
    search_n_destroy_parts(overlap, properties.errordistance);

    //........................

    if (properties.selecting){
	selection(overlap, properties);
    }


    //if (overlap.size() > 1)
    //	weight_info(overlap);
}

void selection(list<Transcript*> &overlap, Properties &properties){
    // seek for the higherest priority with a complete transcript in the current overlap
    overlap.sort(compare_priority);
    int highest_complete_priority = 0;
    list<Transcript*> completeHighPriorityTranscripts;
    for (list<Transcript*>::iterator it = overlap.begin(); it != overlap.end(); it++){
	if ((*it)->tl_complete.first && (*it)->tl_complete.second){
	    if (highest_complete_priority < (*it)->priority){
		completeHighPriorityTranscripts.clear();
		highest_complete_priority = (*it)->priority;
		completeHighPriorityTranscripts.push_back(*it);
	    }else if (highest_complete_priority == (*it)->priority){
		completeHighPriorityTranscripts.push_back(*it);
	    }
	}
    }

    // if there is such a highest priority delete transcripts with lower priority and some uncomplete ones and that ones with low score and so on...
    if (highest_complete_priority){
	if (properties.genemodel == "bacterium")
	    for (list<Transcript*>::iterator it = overlap.begin(); it != overlap.end(); it++){
		if (highest_complete_priority > (*it)->priority){
		    if ((*it)->tl_complete.second == false){
			it = overlap.erase(it);
			it--;
		    }else{
			for (list<Transcript*>::iterator itInside = completeHighPriorityTranscripts.begin(); itInside != completeHighPriorityTranscripts.end(); itInside++){
			    if ((*itInside)->tes == (*it)->tes && (*itInside)->strand == (*it)->strand){
				it = overlap.erase(it);
				it--;
				break;
			    }
			    int lowInside;
			    int highInside;
			    int high;
			    int low;
			    if ((*itInside)->strand == '+'){
				lowInside = (*itInside)->tis;
				highInside = (*itInside)->tes;
			    }else{
				lowInside = (*itInside)->tes;
				highInside = (*itInside)->tis;
			    }
			    if ((*it)->strand == '+'){
				low = (*it)->tis;
				high = (*it)->tes;
			    }else{
				low = (*it)->tes;
				high = (*it)->tis;
			    }
			    if ((low > lowInside && low < highInside) && (high < highInside || ((*itInside)->boha.first != 1 && (*itInside)->boha.second != 1))){
				it = overlap.erase(it);
				it--;
				break;
			    }
			    if ((high < highInside && high > lowInside) && (low > lowInside || ((*itInside)->boha.first != 1 && (*itInside)->boha.second != 1))){
				it = overlap.erase(it);
				it--;
				break;
			    }
			}
		    }
		}else if (highest_complete_priority == (*it)->priority){
		    if (!(*it)->tl_complete.first || !(*it)->tl_complete.second){
			it = overlap.erase(it);
			it--;		
		    }else{
			for (list<Transcript*>::iterator itInside = completeHighPriorityTranscripts.begin(); itInside != completeHighPriorityTranscripts.end(); itInside++){
			    if ((*itInside)->tes == (*it)->tes && (*itInside)->strand == (*it)->strand){
				if (simpleProkScore(*it) < simpleProkScore(*itInside)){
				    it = overlap.erase(it);
				    it--;
				    break;
				}
			    }
			}
		    }
		}
	    }
    }else{
	if (highest_complete_priority){
	    for (list<Transcript*>::iterator it = overlap.begin(); it != overlap.end(); it++){
		if (highest_complete_priority > (*it)->priority){
		    it = overlap.erase(it);
		    it--;
		}else if ((!(*it)->tl_complete.first || !(*it)->tl_complete.second) && (highest_complete_priority == (*it)->priority)){
		    for (list<Transcript*>::iterator itInside = completeHighPriorityTranscripts.begin(); itInside != completeHighPriorityTranscripts.end(); itInside++){
			if (areSimilar(*it,*itInside)){
			    it = overlap.erase(it);
			    it--;
			    break;
			}
		    }
		}
	    }
	}
    }
}

void joinCall(list<Transcript*> &overlap, list<Transcript> &new_transcripts, Properties &properties){
    // test if some exon is to close to prediction range border
    tooCloseToBorder(overlap, new_transcripts, '3', properties.errordistance);
    // join stop codons
    join(overlap, new_transcripts, '3');
    search_n_destroy_doublings(overlap, properties.errordistance, false);

    // delete ?special? transcripts
    for (list<Transcript*>::iterator it = overlap.begin(); it != overlap.end(); it++){
	if ((*it)->correction_ancestor.second && (*it)->tl_complete.second){
	    if (!((*it)->joinpartner.second->boha.second)){
		(*it)->must_be_deleted = true;
		for (list<Transcript*>::iterator iti = (*it)->descendant.second.begin(); iti != (*it)->descendant.second.end(); iti++){
		    (*iti)->must_be_deleted = true;
		}
	    }
	}
    }
    for (list<Transcript*>::iterator it = overlap.begin(); it != overlap.end(); it++){
	if ((*it)->must_be_deleted == true){
	    it = overlap.erase(it);
	    it--;
	}
    }

    // same for the start codon
    tooCloseToBorder(overlap, new_transcripts, '5', properties.errordistance);
    join(overlap, new_transcripts, '5');
    search_n_destroy_doublings(overlap, properties.errordistance, false);

    for (list<Transcript*>::iterator it = overlap.begin(); it != overlap.end(); it++){
	if ((*it)->correction_ancestor.first && (*it)->tl_complete.first){
	    if (!((*it)->joinpartner.first->boha.first)){
		(*it)->must_be_deleted = true;
		for (list<Transcript*>::iterator iti = (*it)->descendant.first.begin(); iti != (*it)->descendant.first.end(); iti++){
		    (*iti)->must_be_deleted = true;
		}
	    }
	}
    }
    for (list<Transcript*>::iterator it = overlap.begin(); it != overlap.end(); it++){
	if ((*it)->must_be_deleted == true){
	    it = overlap.erase(it);
	    it--;
	}
    }
}

void compareAndSplit(list<Transcript*> &overlap, Properties &properties){
    list<Transcript*> priorityTx1;
    list<Transcript*> priorityTx2;
    for (list<Transcript*>::iterator it = overlap.begin(); it != overlap.end(); it++){
	if ((*it)->priority == *properties.priorities.begin())
	    priorityTx1.push_back(*it);
	else
	    priorityTx2.push_back(*it);
    }

    list<Transcript*> equal1;
    list<Transcript*> equal2;
    list<Transcript*> same_stop1;
    list<Transcript*> same_stop2;
    list<Transcript*> unequal1;
    list<Transcript*> unequal2;
    for (list<Transcript*>::iterator it = priorityTx1.begin(); it != priorityTx1.end(); it++){
	list<Transcript*>::iterator it_temp = it;
	it_temp++;
	for (list<Transcript*>::iterator it_inside = priorityTx2.begin(); it_inside != priorityTx2.end(); it_inside++){
	    pair<bool,bool> who_is_part = is_part_of(*it, *it_inside);
	    if (who_is_part.first && who_is_part.second){
		(*it)->boha.first = 2;
		(*it_inside)->boha.first = 2;
	    }else if ((*it)->tes == (*it_inside)->tes){
		if ((*it)->boha.first != 2)
		    (*it)->boha.first = 1;
		if ((*it_inside)->boha.first != 2)
		    (*it_inside)->boha.first = 1;

	    }else{
		if ((*it)->boha.first != 2 && (*it)->boha.first != 1)
		    (*it)->boha.first = 0;
		if ((*it_inside)->boha.first != 2 && (*it_inside)->boha.first != 1)
		    (*it_inside)->boha.first = 0;
	    }
	}
	if ((*it)->boha.first == 2)
	    equal1.push_back(*it);
	if ((*it)->boha.first == 1)
	    same_stop1.push_back(*it);
	if ((*it)->boha.first == 0)
	    unequal1.push_back(*it);
    }
    for (list<Transcript*>::iterator it = priorityTx2.begin(); it != priorityTx2.end(); it++){
	if ((*it)->boha.first == 2)
	    equal2.push_back(*it);
	if ((*it)->boha.first == 1)
	    same_stop2.push_back(*it);
	if ((*it)->boha.first == 0)
	    unequal2.push_back(*it);
    }
    string filenameEqual1 = "equal1.gtf", filenameEqual2 = "equal2.gtf", filenameStop1 = "same_stop1.gtf", filenameStop2 = "same_stop2.gtf", filenameUne1 = "unequal1.gtf", filenameUne2 = "unequal2.gtf";
    saveOverlap(equal1, filenameEqual1, properties);
    saveOverlap(equal1, filenameEqual2, properties);
    saveOverlap(same_stop1, filenameStop1, properties);
    saveOverlap(same_stop2, filenameStop2, properties);
    saveOverlap(unequal1, filenameUne1, properties);
    saveOverlap(unequal2, filenameUne2, properties);
}

double simpleProkScore(Transcript const* tx){
    // calculates a score for comparision; completeness shouldnt be used
    if (tx->exon_list.size() != 1)
	return 0;
    else
	return tx->exon_list.front().score;
}

bool areSimilar(Transcript const* t1, Transcript const* t2){
    // is true, if transcripts are similar
    if (t1->strand != t2->strand){
        return false;
    }
    pair<unsigned int,unsigned int> overHangNumber =  make_pair(0, 0);
    unsigned int equalExonNumber = 0;
    pair<int,int> whoHangOver = make_pair(0, 0);
    pair<float,float> overHangScore = make_pair(0, 0);
    list<Exon>::const_iterator it1 = t1->exon_list.begin();
    list<Exon>::const_iterator it2 = t2->exon_list.begin();

    //int lowestEqualExonBase;
    //int highestEqualExonBase;
    //bool lowestEqualExonFound = false;
    while (it1 != t1->exon_list.end() || it2 != t2->exon_list.end()){
        if (it1 != t1->exon_list.end() && it2 != t2->exon_list.end()){
            if ((*it1).from == (*it2).from){
                if ((*it1).to == (*it2).to){
                    //if (!lowestEqualExonFound){lowestEqualExonFound = true; lowestEqualExonBase = (*it1).from;}
                    //highestEqualExonBase = (*it1).to;
                    equalExonNumber++; it1++; it2++;
                }else if ((*it1).to < (*it2).to){it1++;}else{it2++;}
            }else if ((*it1).from < (*it2).from){
                if(it2 == t2->exon_list.begin()){whoHangOver.first = 1;}
                if ((*it1).to < (*it2).from){
                    if(it2 == t2->exon_list.begin()){overHangNumber.first++; overHangScore.first += (*it1).score;}
                    it1++;
                }else{
                    if ((*it1).to == (*it2).to){it1++; it2++;}else if ((*it1).to < (*it2).to){it1++;}else{it2++;}
                }
            }else{
                if(it1 == t1->exon_list.begin()){whoHangOver.first = 2;}
                if ((*it1).from > (*it2).to){
                    if(it1 == t1->exon_list.begin()){overHangNumber.first++; overHangScore.first += (*it2).score;}
                    it2++;
                }else{
                    if ((*it1).to == (*it2).to){it1++; it2++;}else if ((*it1).to < (*it2).to){it1++;}else{it2++;}
                }
            }
        }else if (it1 == t1->exon_list.end()){
            whoHangOver.second = 2;
            list<Exon>::const_iterator itemp = it1; itemp--;
            if ((*itemp).to < (*it2).from){overHangNumber.second++; overHangScore.second += (*it2).score;} 
            it2++;
        }else{
            whoHangOver.second = 1;
            list<Exon>::const_iterator itemp = it2; itemp--;
            if ((*itemp).to < (*it1).from){overHangNumber.second++; overHangScore.second += (*it1).score;} 
            it1++;
        }
    }
    
    /*if (lowestEqualExonFound){
      float scoreSum = 0;
      int countLoopRuns = 0;
      for (list<Exon>::const_iterator it = t1->exon_list.begin(); (it != t1->exon_list.end() && (*it).from != lowestEqualExonBase); it++){
      scoreSum += (*it1).score;
      countLoopRuns++;
      }
      if (countLoopRuns > 0){
      cout << scoreSum/countLoopRuns << endl;
      }

      scoreSum = 0;
      countLoopRuns = 0;
      for (list<Exon>::const_iterator it = t2->exon_list.begin(); (it != t2->exon_list.end() && (*it).from != lowestEqualExonBase); it++){
      scoreSum += (*it1).score;
      countLoopRuns++;
      }
      if (countLoopRuns > 0){
      cout << scoreSum/countLoopRuns << endl;
      }

      scoreSum = 0;
      countLoopRuns = 0;
      bool behindHighestExon = false;
      for (list<Exon>::const_iterator it = t1->exon_list.begin(); it != t1->exon_list.end(); it++){
      if (behindHighestExon){
      scoreSum += (*it).score;
      countLoopRuns++;
      }
      if ((*it).to == highestEqualExonBase){behindHighestExon = true;}
      }
      if (countLoopRuns > 0){
      cout << scoreSum/countLoopRuns << endl;
      }

      scoreSum = 0;
      countLoopRuns = 0;
      behindHighestExon = false;
      for (list<Exon>::const_iterator it = t2->exon_list.begin(); it != t2->exon_list.end(); it++){
      if (behindHighestExon){
      scoreSum += (*it).score;
      countLoopRuns++;
      }
      if ((*it).to == highestEqualExonBase){behindHighestExon = true;}
      }
      if (countLoopRuns > 0){
      cout << scoreSum/countLoopRuns << endl;
      }
      }*/

    float equalBorder = 0.5;            // maybe as option
    float overHangScoreBorder = 0.5;    // maybe as option
                                        // maybe also as option: [border = score1 - score2]; maybe output this score to use it later
    overHangScore.first /= overHangNumber.first;
    overHangScore.second /= overHangNumber.second;
    if ((whoHangOver.first == 0 && overHangNumber.first != 0) || (whoHangOver.second == 0 && overHangNumber.second != 0)){cerr << "Something went wrong in areSimilar at Pos. 1" << endl;}
    float overlapExonNumberT1 = t1->exon_list.size();
    // float overlapExonNumberT2 = t2->exon_list.size();
    if (whoHangOver.first == 1){overlapExonNumberT1 -= overHangNumber.first;}//else if (whoHangOver.first == 2){overlapExonNumberT2 -= overHangNumber.first;}
    if (whoHangOver.second == 1){overlapExonNumberT1 -= overHangNumber.second;}//else if (whoHangOver.second == 2){overlapExonNumberT2 -= overHangNumber.second;}

    float overlapEqualRatioT1 = (float) (equalExonNumber / overlapExonNumberT1);
    // float overlapEqualRatioT2 = (float) (equalExonNumber / overlapExonNumberT2);

    float centerScoreT1 = 0;
    float centerScoreT2 = 0;
    unsigned int centerNumberT1 = t1->exon_list.size();
    unsigned int centerNumberT2 = t2->exon_list.size();
    for (list<Exon>::const_iterator it1 = t1->exon_list.begin(); it1 != t1->exon_list.end(); it1++){centerScoreT1 += (*it1).score;}
    for (list<Exon>::const_iterator it2 = t2->exon_list.begin(); it2 != t2->exon_list.end(); it2++){centerScoreT2 += (*it2).score;}
    if (whoHangOver.first == 1){
        centerScoreT1 -= overHangScore.first;
        centerNumberT1 -= overHangNumber.first;
    }else if (whoHangOver.first == 2){
        centerScoreT2 -= overHangScore.first;
        centerNumberT2 -= overHangNumber.first;
    }
    if (whoHangOver.second == 1){
        centerScoreT1 -= overHangScore.second;
        centerNumberT1 -= overHangNumber.second;
    }else if (whoHangOver.second == 2){
        centerScoreT2 -= overHangScore.second;
        centerNumberT2 -= overHangNumber.second;
    }


    if (overlapEqualRatioT1 < equalBorder){
        if (centerScoreT2/centerNumberT2 < 0.8 || centerScoreT2/centerNumberT2 < 2*(centerScoreT1/centerNumberT1)){
            return false;
        }
    }

    if ((overHangScore.first > overHangScoreBorder) && (overHangNumber.first > 2 || overHangNumber.first > centerNumberT1)){
        if (whoHangOver.first == 1){
            return false;
        }
    }
    if ((overHangScore.second > overHangScoreBorder) && (overHangNumber.second > 2 || overHangNumber.second > centerNumberT1)){
        if (whoHangOver.second == 1){
            return false;
        }
    }
    return true;
}

void tooCloseToBorder(list<Transcript*> &overlap, list<Transcript> &new_transcripts, char side, int errordistance){
    if (errordistance != -1){
	list<Transcript*> new_overlap_part;
	for (list<Transcript*>::iterator it = overlap.begin(); it != overlap.end(); it++){
	    if (((*it)->strand == '+' && side == '3') || ((*it)->strand == '-' && side == '5')){
		if ((*it)->pred_range.second){
		    int distance = (*it)->pred_range.second - (*it)->exon_list.back().to;
		    if (distance <= errordistance){
			if ((*it)->exon_list.size() == 1){
			    if (overlap.size() <= 1){(*it)->boha.first = 1; (*it)->boha.second = 1; return;}
			    it = overlap.erase(it);
			    it--;
			    continue;
			}
			Transcript tx_new = (*(*it));
			if (tx_new.strand == '+'){
			    (*it)->boha.second = 1;
			    tx_new.boha.second = 2;
			    tx_new.correction_ancestor.second = (*it);
			    if (tx_new.exon_list.back().feature == "CDS"){
				if (tx_new.tl_complete.second == true){
				    tx_new.tl_complete.second = false;
				}
			    }
			}else{
			    (*it)->boha.first = 1;
			    tx_new.boha.first = 2;
			    tx_new.correction_ancestor.first = (*it);
			    if (tx_new.exon_list.back().feature == "CDS"){
				if (tx_new.tl_complete.first == true){
				    tx_new.tl_complete.first = false;
				}
			    }
			}
			tx_new.exon_list.pop_back();
			new_transcripts.push_back(tx_new);
			new_overlap_part.push_back(&new_transcripts.back());
		    }
		}
	    }else if (((*it)->strand == '+' && side == '5') || ((*it)->strand == '-' && side == '3')){
		if ((*it)->pred_range.first){
		    int distance = (*it)->exon_list.front().from - (*it)->pred_range.first;
		    if (distance <= errordistance){
			if ((*it)->exon_list.size() == 1){
			    if (overlap.size() <= 1){(*it)->boha.first = 1; (*it)->boha.second = 1; return;}
			    it = overlap.erase(it);
			    it--;
			    continue;
			}
			Transcript tx_new = (*(*it));
			if (tx_new.strand == '+'){
			    (*it)->boha.first = 1;
			    tx_new.boha.first = 2;
			    tx_new.correction_ancestor.first = (*it);
			    if (tx_new.exon_list.front().feature == "CDS"){
				if (tx_new.tl_complete.first == true){
				    tx_new.tl_complete.first = false;
				}
			    }
			}else{
			    (*it)->boha.second = 1;
			    tx_new.boha.second = 2;
			    tx_new.correction_ancestor.second = (*it);
			    if (tx_new.exon_list.front().feature == "CDS"){
				if (tx_new.tl_complete.second == true){
				    tx_new.tl_complete.second = false;
				}
			    }
			}
			tx_new.exon_list.pop_front();
			new_transcripts.push_back(tx_new);
			new_overlap_part.push_back(&new_transcripts.back());
		    }
		}
	    }
	}
	overlap.merge(new_overlap_part);
    }
}

void search_n_destroy_doublings(list<Transcript*> &overlap, int errordistance, bool ab_initio){
    // delete all transcripts that are completly part of another transcript (in particular all exons are also in the other transcript and vice versa); the one with the lesser priority will be deleted
    if (overlap.size() > 1){
	for (list<Transcript*>::iterator it = overlap.begin(); it != overlap.end(); it++){
	    list<Transcript*>::iterator it_temp = it;
	    it_temp++;
	    for (list<Transcript*>::iterator it_inside = it_temp; it_inside != overlap.end(); it_inside++){
		if (overlap.size() <= 1){return;}
		pair<bool,bool> who_is_part = is_part_of(*it, *it_inside);
		if (who_is_part.first == true){
		    if (who_is_part.second == true){
			if (ab_initio){
			    (*it)->supporter.push_front(*it_inside);
			    (*it_inside)->supporter.push_front(*it);
			}
			int mindist1 = -1;
			int mindist2 = -1;
			if ((*it)->pred_range.first && (*it)->pred_range.second){
			    mindist1 = min((*it)->exon_list.front().from - (*it)->pred_range.first, (*it)->pred_range.second - (*it)->exon_list.back().to);
			}
			if ((*it_inside)->pred_range.first && (*it_inside)->pred_range.second){
			    mindist2 = min((*it_inside)->exon_list.front().from - (*it_inside)->pred_range.first, (*it_inside)->pred_range.second - (*it_inside)->exon_list.back().to);
			}
			if ((mindist1 == -1 || mindist1 > errordistance) && (mindist2 == -1 || mindist2 > errordistance)){
			    if ((*it)->priority < (*it_inside)->priority){
				it = overlap.erase(it);
				it_inside = it;
			    }else{
				it_inside = overlap.erase(it_inside);
				it_inside--;
			    }
			}else if ((mindist1 < mindist2) || mindist2 == -1){
			    it = overlap.erase(it);
			    it_inside = it;
			}else{
			    it_inside = overlap.erase(it_inside);
			    it_inside--;
			}
		    }else{
			if (ab_initio){
			    (*it_inside)->supporter.push_front(*it);
			}
			//it = overlap.erase(it);
			//it_inside = it;
		    }
		}else{
		    if (who_is_part.second == true){
			if (ab_initio){
			    (*it)->supporter.push_front(*it_inside);
			}
			//it_inside = overlap.erase(it_inside);
			//it_inside--;
		    }
		}
	    }
	}
    }
}

void search_n_destroy_parts(list<Transcript*> &overlap, int errordistance){
    // delete all transcripts that are part of another transcript
    if (overlap.size() > 1){
	for (list<Transcript*>::iterator it = overlap.begin(); it != overlap.end(); it++){
	    list<Transcript*>::iterator it_temp = it;
	    it_temp++;
	    for (list<Transcript*>::iterator it_inside = it_temp; it_inside != overlap.end(); it_inside++){
		if (overlap.size() <= 1){return;}
		pair<bool,bool> who_is_part = is_part_of(*it, *it_inside);
		if (who_is_part.first == true){
		    if (who_is_part.second == true){
		    }else{
			it = overlap.erase(it);
			it_inside = it;
		    }
		}else{
		    if (who_is_part.second == true){
			it_inside = overlap.erase(it_inside);
			it_inside--;
		    }
	    	}
	    }
	}
    }
}

void join(list<Transcript*> &overlap, list<Transcript> &new_transcripts, char side){
    // devides an overlap in a list of start/stop codon donors and in an other list of start/stop codon acceptors and joins every pair (one of each list) if they are combinable
    list<Transcript*> new_overlap_part;
    list<Transcript*> donor;
    list<Transcript*> acceptor;
    for (list<Transcript*>::iterator it = overlap.begin(); it != overlap.end(); it++){
	if (side == '3'){
	    if ((*it)->tl_complete.second){
		donor.push_front(*it);
	    }else{
		acceptor.push_front(*it);
	    }
	}else if (side == '5'){
	    if ((*it)->tl_complete.first){
		donor.push_front(*it);
	    }else{
		acceptor.push_front(*it);
	    }
	}
    }
    for (list<Transcript*>::iterator it = acceptor.begin(); it != acceptor.end(); it++){
	for (list<Transcript*>::iterator it_donor = donor.begin(); it_donor != donor.end(); it_donor++){
	    if ((*it)->strand == (*it_donor)->strand){
		int fitting_case = is_combinable(*it, *it_donor, (*it)->strand, side);
		if (fitting_case){
		    Transcript tx_new;
		    joining(*it, *it_donor, (*it)->strand, side, tx_new, fitting_case);
		    new_transcripts.push_back(tx_new);
		    new_overlap_part.push_back(&new_transcripts.back());
		    if (side == '3'){
			(*it)->descendant.second.push_back(&new_transcripts.back());
		    }else{
			(*it)->descendant.first.push_back(&new_transcripts.back());
		    }
		}
	    }
	}
    }
    overlap.merge(new_overlap_part);
}

void joining(Transcript* t1, Transcript* t2, char strand, char side, Transcript &tx_new, int fitting_case){
    // joins transcripts in one direction so that every suitable exon will be transferred and returns a new "joined" transcript without deleting the old ones
    int minimum_intron_length = 20;								// also defined in is_combinable
    tx_new = *t1;
    //	cout << "----------- " << fitting_case << endl;
    bool are_at_add_part;
    list<Exon> temp_exon_list;
    temp_exon_list.clear();
    switch (fitting_case){
    case 0:
	cerr << "WARNING: shouldnt happen (in joining())!" << endl;
	break;
    case 1:
	if (strand == '+'){
	    tx_new.tes = t2->tes;
	    tx_new.tl_complete.second = t2->tl_complete.second;
	    tx_new.joinpartner.second = t2;
	}else{
	    tx_new.tis = t2->tis;
	    tx_new.tl_complete.first = t2->tl_complete.first;
	    tx_new.joinpartner.first = t2;
	}
	are_at_add_part = false;
	for (list<Exon>::iterator it = t2->exon_list.begin(); it != t2->exon_list.end(); it++){
	    if (are_at_add_part){
		tx_new.exon_list.push_back(*it);
	    }
	    if ((*t1).exon_list.back().to >= (*it).from && (*t1).exon_list.back().to <= (*it).to){
		tx_new.exon_list.back().to = (*it).to;
		if (strand == '-'){tx_new.exon_list.back().frame = (*it).frame;}
		are_at_add_part = true;
	    }
	}
	break;
    case 2:
	if (strand == '+'){
	    tx_new.tes = t2->tes;
	    tx_new.tl_complete.second = t2->tl_complete.second;
	    tx_new.joinpartner.second = t2;
	}else{
	    tx_new.tis = t2->tis;
	    tx_new.tl_complete.first = t2->tl_complete.first;
	    tx_new.joinpartner.first = t2;
	}
	are_at_add_part = false;
	for (list<Exon>::iterator it = t2->exon_list.begin(); it != t2->exon_list.end(); it++){
	    if ((*t1).exon_list.back().to <= ((*it).from - minimum_intron_length)){
		are_at_add_part = true;
	    }
	    if (are_at_add_part){
		tx_new.exon_list.push_back(*it);
	    }
	}
	break;
    case 3:
	if (strand == '+'){
	    tx_new.tis = t2->tis;
	    tx_new.tl_complete.first = t2->tl_complete.first;
	    tx_new.joinpartner.first = t2;
	}else{
	    tx_new.tes = t2->tes;
	    tx_new.tl_complete.second = t2->tl_complete.second;
	    tx_new.joinpartner.second = t2;
	}
	//are_at_add_part = true;
	for (list<Exon>::iterator it = t2->exon_list.begin(); it != t2->exon_list.end(); it++){
	    if ((*t1).exon_list.front().from >= (*it).from && (*t1).exon_list.front().from <= (*it).to){
		tx_new.exon_list.front().from = (*it).from;
		if (strand == '+'){tx_new.exon_list.front().frame = (*it).frame;}
		tx_new.exon_list.merge(temp_exon_list);
		//are_at_add_part = false;
		break;
	    }
	    //if (are_at_add_part){
	    temp_exon_list.push_back(*it);
	    //}
	}
	break;
    case 4:
	if (strand == '+'){
	    tx_new.tis = t2->tis;
	    tx_new.tl_complete.first = t2->tl_complete.first;
	    tx_new.joinpartner.first = t2;
	}else{
	    tx_new.tes = t2->tes;
	    tx_new.tl_complete.second = t2->tl_complete.second;
	    tx_new.joinpartner.second = t2;
	}
	//are_at_add_part = true;
	for (list<Exon>::iterator it = t2->exon_list.begin(); it != t2->exon_list.end(); it++){
	    if ((*t1).exon_list.front().from < ((*it).to + minimum_intron_length)){
		tx_new.exon_list.merge(temp_exon_list);
		//are_at_add_part = false;
		break;
	    }
	    //if (are_at_add_part){
	    temp_exon_list.push_back(*it);
	    //}
	}
	break;
    default:
	cerr << "WARNING: unexpected case (in joining())!" << endl;
    }
    //cout << tx_new.exon_list.size() << endl;
}

int is_combinable(Transcript const* t1, Transcript const* t2, char strand, char side){
    // is true,	if the first exon which is minimal the minimum_intron_length away from the last exon of the other transcript in the appropriate direction
    // 			&& the exons at these positions are frame-compatible
    // 			&& the transcripts are overlapping					// maybe we can improve something here, that combinable non-overlaping transcripts gets true (but be carefull)
    int minimum_intron_length = 20;								// also defined in joining
    if ((max( (*t1).tis,(*t1).tes) < min((*t2).tis,(*t2).tes)) || (min((*t1).tis,(*t1).tes) > max((*t2).tis,(*t2).tes))){
	return 0;
    }
    if ((strand == '+' && side == '3') || (strand == '-' && side == '5')){
	for (list<Exon>::const_iterator it = t2->exon_list.begin(); it != t2->exon_list.end(); it++){
	    if (((*t1).exon_list.back().to >= (*it).from) && ((*t1).exon_list.back().to <= (*it).to)){
		if ((*t1).exon_list.back().frame == -1 && (*it).frame == -1){return true;}
		else if (((*t1).exon_list.back().frame != -1 && (*it).frame == -1) || ((*t1).exon_list.back().frame == -1 && (*it).frame != -1)){return false;}
		if ((strand == '+') && ((3 - (((*t1).exon_list.back().to - (*t1).exon_list.back().from) - (*t1).exon_list.back().frame) % 3) % 3) == ((3 - (((*t1).exon_list.back().to - (*it).from) - (*it).frame) % 3) % 3) ){
		    return 1;
		}else if ((strand == '-') && (*t1).exon_list.back().frame == ((3 - (((*it).to - (*t1).exon_list.back().to) - (*it).frame) % 3) % 3) ){
		    return 1;
		}else{
		    return 0;
		}
	    }else{
		if ((*t1).exon_list.back().to <= ((*it).from - minimum_intron_length)){
		    if ((strand == '+') && ((*it).frame == (3 - ( ((*t1).exon_list.back().to - (*t1).exon_list.back().from + 1) - (*t1).exon_list.back().frame) % 3) % 3)){
			return 2;
		    }else if ((strand == '-') && ((*t1).exon_list.back().frame == (3 - ( ((*it).to - (*it).from + 1) - (*it).frame) % 3) % 3)){
			return 2;
		    }else{
			return 0;
		    }
		}
	    }
	}
    }else if ((strand == '+' && side == '5') || (strand == '-' && side == '3')){
	for (list<Exon>::const_iterator it = t2->exon_list.begin(); it != t2->exon_list.end(); it++){
	    if (((*t1).exon_list.front().from >= (*it).from) && ((*t1).exon_list.front().from <= (*it).to)){
		if ((*t1).exon_list.front().frame == -1 && (*it).frame == -1){return true;}
		else if (((*t1).exon_list.back().frame != -1 && (*it).frame == -1) || ((*t1).exon_list.back().frame == -1 && (*it).frame != -1)){return false;}
		if ( (strand == '-') && ((3 - (((*t1).exon_list.front().to - (*t1).exon_list.front().from) - (*t1).exon_list.front().frame) % 3) % 3) == ((3 - (((*it).to - (*t1).exon_list.front().from) - (*it).frame) % 3) % 3) ){
		    return 3;
		}else if ( (strand == '+') && (*t1).exon_list.front().frame == ((3 - (((*t1).exon_list.front().from - (*it).from) - (*it).frame) % 3) % 3) ){
		    return 3;
		}else{
		    return 0;
		}
	    }else{
		if ((*t1).exon_list.front().from < ((*it).to + minimum_intron_length)){
		    if (it != t2->exon_list.begin()){
			it--;
		    }else{
			return 0;
		    }
		    if ((strand == '-') && ((*it).frame == (3 - ( ((*t1).exon_list.front().to - (*t1).exon_list.front().from + 1) - (*t1).exon_list.front().frame) % 3) % 3)){
			return 4;
		    }else if ((strand == '+') && ((*t1).exon_list.front().frame == (3 - ( ((*it).to - (*it).from + 1) - (*it).frame) % 3) % 3)){
			return 4;
		    }else{
			return 0;
		    }
		}
	    }
	}
    }else{									// no '+' or '-' strand
	return 0;						// here we have an unexpected case and so we dont wanna combine
    }
    return 0;
}

pair<bool,bool> is_part_of(Transcript const* t1, Transcript const* t2)
{
    // is true,false or false,true if one transcript contains the other completely
    // is true,true if the transcripts are equal in exons			// this case could completly replace compare_transcripts
    bool t1_is_part = true;
    bool t2_is_part = true;
    if (t1->strand != t2->strand){
        t1_is_part = false;
        t2_is_part = false;
    }
    if (t1->exon_list.size() > t2->exon_list.size()){
	t1_is_part = false;
    }
    if (t1->exon_list.size() < t2->exon_list.size()){
	t2_is_part = false;
    }
    if (t1->tl_complete.first && t2->tl_complete.first && t1->tis != t2->tis){
	t1_is_part = false;
	t2_is_part = false;
    }
    if (t1->tl_complete.second && t2->tl_complete.second && t1->tes != t2->tes){
	t1_is_part = false;
	t2_is_part = false;
    }
    list<Exon>::const_iterator it1 = t1->exon_list.begin();
    list<Exon>::const_iterator it2 = t2->exon_list.begin();
    while (t1_is_part == true || t2_is_part == true){
	if ((*it1).from == (*it2).from){
	    if((*it1).to == (*it2).to && (*it1).frame == (*it2).frame){
		it1++;
		it2++;
	    }
	    else{
		t1_is_part = false;
		t2_is_part = false;
		break;
	    }
	}else if ((*it1).from > (*it2).from){
	    t2_is_part = false;
	    it2++;
	}else{
	    t1_is_part = false;
	    it1++;
	}
	if (it1 == t1->exon_list.end() && it2 == t2->exon_list.end()){
	    break;
	}else{
	    if (it1 == t1->exon_list.end() && !(it2 == t2->exon_list.end())){
		t2_is_part = false;
		break;
	    }
	    if (it2 == t2->exon_list.end() && !(it1 == t1->exon_list.end())){
		t1_is_part = false;
		break;
	    }
	}
    }
    return make_pair(t1_is_part, t2_is_part);
}

bool compare_priority(Transcript const* lhs, Transcript const* rhs){
    // returns true, if the priority of the first element is higher; used to sort transcripts by priority
    return ( lhs->priority > rhs->priority );
}

bool check_frame_annotation(Transcript const &transcript){
    // returns true, if the frame annotation of the transcript is correct
    list<Exon>::const_iterator it2;
    if (transcript.strand == '+'){
	for (list<Exon>::const_iterator it1 = transcript.exon_list.begin(); it1 != transcript.exon_list.end(); it1++){
	    it2 = it1;
	    it2++;
	    if (it2 == transcript.exon_list.end())
		break;
	    if ((*it1).frame != -1 && (*it2).frame != -1){
		if ((*it2).frame != (3 - ( ((*it1).to - (*it1).from + 1) - (*it1).frame) % 3) % 3){
		    return false;
		}
	    }
	}
    }else if (transcript.strand == '-'){
	for (list<Exon>::const_iterator it1 = transcript.exon_list.end(); it1 != transcript.exon_list.begin(); it1--){
	    if (it1 == transcript.exon_list.end()){
		it1--;
		if (it1 == transcript.exon_list.begin())
		    break;
	    }
	    it2 = it1;
	    it2--;
	    if (it2 == transcript.exon_list.begin())
		break;
	    if ((*it1).frame != -1 && (*it2).frame != -1){
		if ((*it2).frame != (3 - ( ((*it1).to - (*it1).from + 1) - (*it1).frame) % 3) % 3){
		    return false;
		}else
		    if (it2 == transcript.exon_list.begin())
			break;
	    }
	}
    }
    return true;
}

void eval_gtf(list<Transcript*> &overlap, int errordistance){
    string filename = "/home/lars/lars/test_data/eval_test.txt";
    fstream outfile;
    outfile.open(filename, ios::out | ios::app);

    list<Transcript*> annotation;
    list<Transcript*> prediction;
    list<Transcript*> single;
    for (list<Transcript*>::iterator it = overlap.begin(); it != overlap.end(); it++){
	if ((*it)->priority == 3){
	    annotation.push_back(*it);
	}
	search_n_destroy_doublings(annotation, errordistance, false);
	if ((*it)->pred_range.first && (*it)->pred_range.second){
	    prediction.push_back(*it);
	}
	search_n_destroy_doublings(prediction, errordistance, false);
	if ((*it)->priority == 1){
	    single.push_back(*it);
	}
    }
    if (annotation.empty()){
	for (list<Transcript*>::iterator itp = prediction.begin(); itp != prediction.end(); itp++){
	    for (list<Exon>::iterator itpe = (*itp)->exon_list.begin(); itpe != (*itp)->exon_list.end(); itpe++){
		(*itpe).tooFew = 0;
		(*itpe).tooMany = (*itpe).to - (*itpe).from + 1;

		(*itpe).penalty = (*itpe).tooFew + (*itpe).tooMany;
		(*itpe).distance = min((*itpe).from - (*itp)->pred_range.first, (*itp)->pred_range.second - (*itpe).to);
		list<Exon>::iterator itpe_temp = itpe;
		itpe_temp++;
		if (!(itpe == (*itp)->exon_list.begin()) && !(itpe_temp == (*itp)->exon_list.end())){
		    //	outfile << (*itpe).distance << "\t" << (*itpe).penalty << endl;
		}else{
		    outfile << (*itpe).distance << "\t" << (*itpe).penalty << endl;
		}
	    }
	}
	for (list<Transcript*>::iterator its = single.begin(); its != single.end(); its++){
	    for (list<Exon>::iterator itpe = (*its)->exon_list.begin(); itpe != (*its)->exon_list.end(); itpe++){
		(*itpe).tooFew = 0;
		(*itpe).tooMany = (*itpe).to - (*itpe).from + 1;

		(*itpe).penalty = (*itpe).tooFew + (*itpe).tooMany;
		count1++;
		value += (*itpe).penalty;
	    }
	}
    }else{
	for (list<Transcript*>::iterator itp = prediction.begin(); itp != prediction.end(); itp++){
	    for (list<Exon>::iterator itpe = (*itp)->exon_list.begin(); itpe != (*itp)->exon_list.end(); itpe++){
		int min_pen = -1;
		for (list<Transcript*>::iterator ita = annotation.begin(); ita != annotation.end(); ita++){
		    for (list<Exon>::iterator itae = (*ita)->exon_list.begin(); itae != (*ita)->exon_list.end(); itae++){
			if ((*itpe).to >= (*itae).from && (*itpe).from <= (*itae).to && (*itp)->strand == (*ita)->strand && strandeq(*itpe,*itae,(*itp)->strand)){
			    int more = 0;
			    int lesser = 0;
			    if ((*itpe).to > (*itae).to){
				more += (*itpe).to - (*itae).to;
			    }else if ((*itpe).to < (*itae).to){
				lesser += (*itae).to - (*itpe).to;
			    }
			    if ((*itpe).from > (*itae).from){
				lesser += (*itpe).from - (*itae).from;
			    }else if ((*itpe).from < (*itae).from){
				more += (*itae).from - (*itpe).from;
			    }
			    if (min_pen == -1 || min_pen > more + lesser){
				(*itpe).tooFew = lesser;
				(*itpe).tooMany = more;
				min_pen = (*itpe).tooFew + (*itpe).tooMany;
				//min_pen = (*itpe).tooFew;
			    }
			}else{
			    if (min_pen == -1 || min_pen > (*itpe).to - (*itpe).from + 1){
				(*itpe).tooFew = 0;
				(*itpe).tooMany = (*itpe).to - (*itpe).from + 1;
				min_pen = (*itpe).tooFew + (*itpe).tooMany;
				//min_pen = (*itpe).tooFew;
			    }
			}
		    }
		}
		(*itpe).penalty = min_pen;
		(*itpe).distance = min((*itpe).from - (*itp)->pred_range.first, (*itp)->pred_range.second - (*itpe).to);

		list<Exon>::iterator itpe_temp = itpe;
		itpe_temp++;
		if (!(itpe == (*itp)->exon_list.begin()) && !(itpe_temp == (*itp)->exon_list.end())){
		    //	outfile << (*itpe).distance << "\t" << (*itpe).penalty << endl;
		}else{
		    outfile << (*itpe).distance << "\t" << (*itpe).penalty << endl;
		}
	    }
	}
	for (list<Transcript*>::iterator its = single.begin(); its != single.end(); its++){
	    for (list<Exon>::iterator itpe = (*its)->exon_list.begin(); itpe != (*its)->exon_list.end(); itpe++){
		int min_pen = -1;
		for (list<Transcript*>::iterator ita = annotation.begin(); ita != annotation.end(); ita++){
		    for (list<Exon>::iterator itae = (*ita)->exon_list.begin(); itae != (*ita)->exon_list.end(); itae++){
			if ((*itpe).to >= (*itae).from && (*itpe).from <= (*itae).to && (*its)->strand == (*ita)->strand && strandeq(*itpe,*itae,(*its)->strand)){
			    int more = 0;
			    int lesser = 0;
			    if ((*itpe).to > (*itae).to){
				more += (*itpe).to - (*itae).to;
			    }else if ((*itpe).to < (*itae).to){
				lesser += (*itae).to - (*itpe).to;
			    }
			    if ((*itpe).from > (*itae).from){
				lesser += (*itpe).from - (*itae).from;
			    }else if ((*itpe).from < (*itae).from){
				more += (*itae).from - (*itpe).from;
			    }
			    if (min_pen == -1 || min_pen > more + lesser){
				(*itpe).tooFew = lesser;
				(*itpe).tooMany = more;
				min_pen = (*itpe).tooFew + (*itpe).tooMany;
				//min_pen = (*itpe).tooFew;
			    }
			}else{
			    if (min_pen == -1 || min_pen > (*itpe).to - (*itpe).from + 1){
				(*itpe).tooFew = 0;
				(*itpe).tooMany = (*itpe).to - (*itpe).from + 1;
				min_pen = (*itpe).tooFew + (*itpe).tooMany;
				//min_pen = (*itpe).tooFew;
			    }
			}
		    }
		}
		(*itpe).penalty = min_pen;
		count1++;
		value += (*itpe).penalty;
	    }
	}
    }
    outfile.close();
}

bool strandeq(Exon ex1, Exon ex2, char strand){
    if (strand == '+'){
	if (ex1.frame == -1 && ex2.frame == -1){return true;}
	if (ex1.from >= ex2.from)
	    if (((ex1.from - ex2.from) + ex1.frame) % 3 == ex2.frame){return true;}else{return false;}
	else
	    if (((ex2.from - ex1.from) + ex2.frame) % 3 == ex1.frame){return true;}else{return false;}
    }else{
	if (ex1.frame == -1 && ex2.frame == -1){return true;}
	if (ex1.to >= ex2.to)
	    if (((ex1.to - ex2.to) + ex2.frame) % 3 == ex1.frame){return true;}else{return false;}
	else
	    if (((ex2.to - ex1.to) + ex1.frame) % 3 == ex2.frame){return true;}else{return false;}
    }
    return false;
}

void output_exon_list(Transcript const* tx){
    // just for semantic tests
    cout << ">> Exon_list of " << tx->t_id << ": " << endl;
    cout << "Priority: " << tx->priority << endl;
    cout << tx->exon_list.size() << " elements, complete start: " << tx->tl_complete.first << ", complete stop: " << tx->tl_complete.second << ", strand: " << tx->strand << endl;
    cout << "start codon of: ";
    if (tx->joinpartner.first){cout << tx->joinpartner.first->t_id;}else{cout << "noone";}
    cout << ", stop codon of: ";
    if (tx->joinpartner.second){cout << tx->joinpartner.second->t_id;}else{cout << "noone";}
    cout << endl;
    if (tx->pred_range.first && tx->pred_range.second){
	cout << "Pred_range: " << tx->pred_range.first << "\t" << tx->pred_range.second << endl;
	cout << "Minimum distance from pred_range borders: " << min(tx->exon_list.front().from - tx->pred_range.first, tx->pred_range.second - tx->exon_list.back().to) << endl;
    }
    cout << "from\t\tto\t\tframe\tfeature\tstart\t\tstop\tchr" << endl;
    for (list<Exon>::const_iterator it = tx->exon_list.begin(); it != tx->exon_list.end(); it++){
	cout << (*it).from << "\t" << (*it).to << "\t" << (*it).frame << "\t" << (*it).feature << "\t" << tx->tis << "\t" << tx->tes << "\t" << (*it).chr << endl;
    }
    cout << "Exon_list end----------------------------------------------------------------" << endl;
}

void weight_info(list<Transcript*> &overlap){
    cout << "ID (priority)\tStaBoha\tStoBoha\tStartjoin\tStopjoin\t#Supp.\ttis\t\ttes" << endl;
    for (list<Transcript*>::iterator it = overlap.begin(); it != overlap.end(); it++){
	cout << (*it)->t_id << " (" << (*it)->priority << ")\t" << (*it)->boha.first << "\t" << (*it)->boha.second << "\t";
	if ((*it)->joinpartner.first)
	    cout << (*((*it)->joinpartner.first)).t_id << " (" << (*((*it)->joinpartner.first)).priority << ")(" << (*((*it)->joinpartner.first)).boha.first << ")\t";
	else{
	    cout << "noone";
	    if ((*it)->tl_complete.first)
		cout << " + \t";
	    else
		cout << " - \t";
	}
	if ((*it)->joinpartner.second)
	    cout << (*((*it)->joinpartner.second)).t_id << " (" << (*((*it)->joinpartner.second)).priority << ")(" << (*((*it)->joinpartner.second)).boha.second << ")\t";
	else{
	    cout << "noone";
	    if ((*it)->tl_complete.second)
		cout << " + \t";
	    else
		cout << " - \t";
	}
	/*for (list<Transcript*>::iterator its = (*it)->supporter.begin(); its != (*it)->supporter.end(); its++){
	  cout << (*its)->t_id << " ";
	  }*/
	cout << (*it)->supporter.size() << "\t" << (*it)->tis << "\t" << (*it)->tes << endl;
    }
    cout << "---------------------------------------------------------------------------------------------------" << endl;
}

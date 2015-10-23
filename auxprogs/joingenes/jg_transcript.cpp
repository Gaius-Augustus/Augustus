#include "jg_transcript.h"
#include "jg_ios.h"

using namespace std;

int count1 = 0;				// only for semantic tests
int value = 0;				// only for semantic tests

bool compareGenomePositions(Transcript const* lhs, Transcript const* rhs){
    if (min(lhs->tis,lhs->tes) != min(rhs->tis,rhs->tes))
	return (min(lhs->tis,lhs->tes) < min(rhs->tis,rhs->tes));
    else
	return (max(lhs->tis,lhs->tes) < max(rhs->tis,rhs->tes));
}

void divideInOverlapsAndConquer(list<Transcript*> &transcript_list, Properties &properties){
    // divide a transcript list in clusters of overlapping transcripts (transitive closure) and work on each cluster separately
    transcript_list.sort(compareGenomePositions);
    list<Transcript*> overlap;
    /*    string filename2 = "/home/lars/lars/test_data/eval_test.txt";
	  fstream outfile2;
	  outfile2.open(filename2, ios::out);
	  outfile2.close();*/

    int max_base = max(transcript_list.front()->tis,transcript_list.front()->tes);
    for (list<Transcript*>::iterator it = transcript_list.begin(); it != transcript_list.end(); it++){
        if (min((*it)->tis,(*it)->tes) < max_base){
            overlap.push_back(*it);
            if (max_base < max((*it)->tis,(*it)->tes)){
                max_base = max((*it)->tis,(*it)->tes);
            }
        }else{
            //eval_gtf(overlap, errordistance);
	    if (!properties.onlyCompare){
		workAtOverlap(overlap, properties);
		saveOverlap(overlap, properties.outFileName, properties);
	    }else{
		compareAndSplit(overlap, properties);
	    }
            overlap.clear();
            max_base = max((*it)->tis,(*it)->tes);
            overlap.push_front(*it);
        }
    }
    if (!properties.onlyCompare){
	workAtOverlap(overlap, properties);
	saveOverlap(overlap, properties.outFileName, properties);
    }else{
	compareAndSplit(overlap, properties);
    }
    overlap.clear();
}

void workAtOverlap(list<Transcript*> &overlap, Properties &properties)
{
    if (overlap.size() > 300){
	cerr << "WARNING: the size of the actual overlap is very big (" << overlap.size() << "). This may lead to a long computation time." << endl; 
    }

    // calls methods for joining transcripts with the target that most of transcripts are complete (have start and stop codon) and delete duplicates and other unwanted transcripts from the overlap
    // search doubling and delete the ones with the lowest priority/score

    search_n_destroy_doublings(overlap, properties, true);
    if (properties.join && properties.genemodel != "bacterium"){
	joinCall(overlap, properties);
    }

    // searchs and destroys transcripts, which are fully in other transcripts
    if (properties.genemodel != "bacterium")
        search_n_destroy_parts(overlap, properties);

    //........................
    if (properties.selecting){
	selection(overlap, properties);
    }
}

void selection(list<Transcript*> &overlap, Properties &properties){
    // seek for the higherest priority with a complete transcript in the current overlap
    if (properties.genemodel == "bacterium"){
	bactSelection(overlap, properties);
    }else{
	eukaSelectionDevelopment(overlap, properties);
    }
}

void bactSelection(list<Transcript*> &overlap, Properties &properties){
    for (list<Transcript*>::iterator it = overlap.begin(); it != overlap.end(); it++){
	if ((*it)->exon_list.front().rangeToBoundary <= properties.errordistance && (*it)->exon_list.front().rangeToBoundary != -1){
	    it = overlap.erase(it);
	    it--;
	    continue;
	}
	if ((*it)->exon_list.back().rangeToBoundary <= properties.errordistance && (*it)->exon_list.back().rangeToBoundary != -1){
	    it = overlap.erase(it);
	    it--;
	    continue;
	}
	if (!(*it)->tl_complete.second || !(*it)->tl_complete.first){
	    it = overlap.erase(it);
	    it--;
	    continue;
	}
    }
    for (list<Transcript*>::iterator it = overlap.begin(); it != overlap.end(); it++){
	list<Transcript*>::iterator itInside = it;
	itInside++;
	while (itInside != overlap.end()){
	    if ((*it)->strand == (*itInside)->strand && (*it)->tes == (*itInside)->tes){
		if ((*it)->priority < (*itInside)->priority)
		{
		    it = overlap.erase(it);
		    it--;
		    break;
		}else if ((*it)->priority > (*itInside)->priority){
		    itInside = overlap.erase(itInside);
		    continue;
		}else{	// priorities equal
		    if (simpleProkScore(*it) < simpleProkScore(*itInside)){
			it = overlap.erase(it);
			it--;
			break;
		    }else{
			itInside = overlap.erase(itInside);
			continue;
		    }
		}
	    }else if((*it)->strand == (*itInside)->strand && (*it)->getCdsFront()->from > (*itInside)->getCdsFront()->from && (*it)->getCdsBack()->to < (*itInside)->getCdsBack()->to){
		if ((*it)->getCdsFront()->score < 0.1){
		    it = overlap.erase(it);
		    it--;
		    break;
		}
	    }else if((*it)->strand == (*itInside)->strand && (*it)->getCdsFront()->from < (*itInside)->getCdsFront()->from && (*it)->getCdsBack()->to > (*itInside)->getCdsBack()->to){
		if ((*itInside)->getCdsFront()->score < 0.1){
		    itInside = overlap.erase(itInside);
		    continue;
		}
	    }
	    itInside++;
	}
    }
}

void eukaSelection(list<Transcript*> &overlap, Properties &properties){
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
	for (list<Transcript*>::iterator it = overlap.begin(); it != overlap.end(); it++){
	    if (highest_complete_priority > (*it)->priority){
		it = overlap.erase(it);
		it--;
	    }else if ((!(*it)->tl_complete.first || !(*it)->tl_complete.second) && (highest_complete_priority == (*it)->priority)){			// TODO: not only for uncomplete transcripts... maybe handle these ones otherwise anyway!!!
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

void joinCall(list<Transcript*> &overlap, Properties &properties){
    // if there are transcripts with and without start/stop codon, they can be joined here, if they fit together; the old transcript remain for the moment
    // test if some exon is to close to prediction range boundary
    tooCloseToBoundary(overlap, properties);
    // join stop codons
    join(overlap, '3', properties);

    search_n_destroy_doublings(overlap, properties, false);

    join(overlap, '5', properties);

    search_n_destroy_doublings(overlap, properties, false);
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

    list<Transcript*> equal;
    list<Transcript*> alternatives1;
    list<Transcript*> alternatives2;
    list<Transcript*> unequal1;
    list<Transcript*> unequal2;

    for (list<Transcript*>::iterator it = priorityTx1.begin(); it != priorityTx1.end(); it++){
	list<Transcript*>::iterator it_temp = it;
	it_temp++;
	for (list<Transcript*>::iterator it_inside = priorityTx2.begin(); it_inside != priorityTx2.end(); it_inside++){
	    pair<bool,bool> who_is_part = is_part_of(*it, *it_inside);
	    if (who_is_part.first && who_is_part.second){
		(*it)->compareValue = EQUAL;
		(*it_inside)->compareValue = EQUAL;
	    }else if ((properties.genemodel == "bacterium" && (*it)->tes && (*it_inside)->tes && (*it)->tes == (*it_inside)->tes) || (properties.genemodel != "bacterium" && alternativeVariants((*it), (*it_inside)))){
		if ((*it)->compareValue != 2)
		    (*it)->compareValue = ALTERNATIVE;
		if ((*it_inside)->compareValue != 2)
		    (*it_inside)->compareValue = ALTERNATIVE;
	    }else{
		if ((*it)->compareValue != 2 && (*it)->compareValue != 1)
		    (*it)->compareValue = UNEQUAL;
		if ((*it_inside)->compareValue != 2 && (*it_inside)->compareValue != 1)
		    (*it_inside)->compareValue = UNEQUAL;
	    }
	}
	if ((*it)->compareValue == 2)
	    equal.push_back(*it);
	if ((*it)->compareValue == 1)
	    alternatives1.push_back(*it);
	if ((*it)->compareValue == 0)
	    unequal1.push_back(*it);
    }
    for (list<Transcript*>::iterator it = priorityTx2.begin(); it != priorityTx2.end(); it++){
	if ((*it)->compareValue == 1)
	    alternatives2.push_back(*it);
	if ((*it)->compareValue == 0)
	    unequal2.push_back(*it);
    }
    string filenameEqual = "cEqual.gtf", filenameStop1 = "cAlternative1.gtf", filenameStop2 = "cAlternative2.gtf", filenameUne1 = "cUnequal1.gtf", filenameUne2 = "cUnequal2.gtf";
    saveOverlap(equal, filenameEqual, properties);
    saveOverlap(alternatives1, filenameStop1, properties);
    saveOverlap(alternatives2, filenameStop2, properties);
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

Transcript* createCopyOf(Transcript* tx, Properties &properties, list<Transcript*> &overlap){
    Transcript* txNew = new Transcript;
    (*txNew) = *tx;

// either
    //txNew->parent->children.push_back(txNew);
// or
    Gene* geneNew = new Gene;
    geneNew->children.push_back(txNew);
    unordered_map<string,Gene*>* empty = NULL;
    geneNew->g_id = nextFreeGeneID(properties, empty);
    (*properties.geneMap)[geneNew->g_id] = geneNew;
    txNew->parent = geneNew;
// either end

    txNew->t_id = nextFreeTxID(txNew->parent, properties, NULL);
    (*properties.transcriptMap)[txNew->t_id] = txNew;
    overlap.push_back(txNew);
    return txNew;
}

void deleteGene(Gene* gene, Properties &properties){
    if (gene->children.size() > 0){cerr << "SEMANTIC ERROR: It is not allowed to delete a gene without deleting transcripts before." << endl; exit( EXIT_FAILURE );}
    (*properties.geneMap).erase(gene->g_id);
    delete gene;
}

void deleteTx(Transcript* tx, Properties &properties){
    (*properties.transcriptMap).erase(tx->t_id);
    for (list<Transcript*>::iterator it = tx->parent->children.begin(); it != tx->parent->children.end(); it++){
	if ((*it) == tx){
	    it = tx->parent->children.erase(it);
	    it--;
	    if (tx->parent->children.empty()){
		deleteGene(tx->parent, properties);
	    }
	}
    }
    delete tx;	// problem with pointer to these transcript in "isCopyOf" and "joinpartner" and others in that kind
}

void tooCloseToBoundary(list<Transcript*> &overlap, Properties &properties){

    // if an errordistance is set, solve problem cases
    if (properties.errordistance < 0){
	return;
    }
    // find out, if there are indirect boundary problems
    for (list<Transcript*>::iterator it = overlap.begin(); it != overlap.end(); it++){

	list<Transcript*>::iterator itInside = it;
	itInside++;
	while (itInside != overlap.end()){
	    if (overlappingCdsWithAnything((*it), (*itInside)) /*&& (*it)->strand == (*itInside)->strand*/){
		// does *it transcend the predictionBoundary of *itInside in tail direction
		if ((*itInside)->pred_range.second && (*it)->getTxEnd() > (*itInside)->pred_range.second){
//		    (*itInside)->boundaryProblem.second = INDIRECT_PROBLEM;
		    (*itInside)->exon_list.back().boundaryProblem = INDIRECT_PROBLEM;
//		    (*itInside)->exon_list.back().indirectBoundProbEnemies.push_back(*it);
		}
		// does *it transcend the predictionBoundary of *itInside in front direction
		if ((*itInside)->pred_range.first && (*it)->getTxStart() < (*itInside)->pred_range.first){
//		    (*itInside)->boundaryProblem.first = INDIRECT_PROBLEM;
		    (*itInside)->exon_list.front().boundaryProblem = INDIRECT_PROBLEM;
//		    (*itInside)->exon_list.front().indirectBoundProbEnemies.push_back(*it);
		}
		// does *itInside transcend the predictionBoundary of *it in tail direction
		if ((*it)->pred_range.second && (*itInside)->getTxEnd() > (*it)->pred_range.second){
//		    (*it)->boundaryProblem.second = INDIRECT_PROBLEM;
		    (*it)->exon_list.back().boundaryProblem = INDIRECT_PROBLEM;
//		    (*it)->exon_list.back().indirectBoundProbEnemies.push_back(*itInside);
		}
		// does *itInside transcend the predictionBoundary of *it in front direction
		if ((*it)->pred_range.first && (*itInside)->getTxStart() < (*it)->pred_range.first){
//		    (*it)->boundaryProblem.first = INDIRECT_PROBLEM;
		    (*it)->exon_list.front().boundaryProblem = INDIRECT_PROBLEM;
//		    (*it)->exon_list.front().indirectBoundProbEnemies.push_back(*itInside);
		}
	    }
	    itInside++;
	}
    }

    // create new transcripts without the problem exons at front and/or back, if there is a boundary problem (direct or indirect) 
    for (list<Transcript*>::iterator it = overlap.begin(); it != overlap.end(); it++){
	unsigned int countCDS = 0;
        for (list<Exon>::iterator ite = (*it)->exon_list.begin(); ite != (*it)->exon_list.end(); ite++){
	    if ((*ite).feature == "CDS"){
		countCDS++;
	    }
	    if (countCDS >= 2){break;}
	}


	if ((*it)->exon_list.front().boundaryProblem >= 2 && (*it)->exon_list.back().boundaryProblem < 2){
	    // erase the transcript, if it has problem boundary exons only
	    if ((*it)->exon_list.size() == 1 || ((*it)->exon_list.front().feature == "CDS" && countCDS <= 1)){
		deleteTx(*it, properties);
		it = overlap.erase(it);
		it--;
		continue;
	    }
	    // create txNew
	    Transcript* txNew = createCopyOf(*it, properties, overlap);

	    string tempFeatureFront = txNew->exon_list.front().feature;

	    // delete problem exon
	    txNew->exon_list.pop_front();

	    // change some values of new transcripts for adaption
	    if (tempFeatureFront == "CDS"){
		if (txNew->strand == '+'){
		    txNew->tl_complete.first = false;
		    txNew->tis = txNew->exon_list.front().from;
		}else{
		    txNew->tl_complete.second = false;
		    txNew->stop_list.clear();
		    txNew->tes = txNew->exon_list.front().from;
		}
	    }else{
		if (txNew->strand == '+'){
		    txNew->tx_complete.first = false;
		    txNew->tss = -1;
		}else{
		    txNew->tx_complete.second = false;
		    txNew->tts = -1;
		}
	    }
	    txNew->exon_list.front().boundaryProblem = FREED;

	}else if ((*it)->exon_list.front().boundaryProblem < 2 && (*it)->exon_list.back().boundaryProblem >= 2){
	    // erase the transcript, if it has problem boundary exons only
	    if ((*it)->exon_list.size() == 1 || ((*it)->exon_list.back().feature == "CDS" && countCDS <= 1)){
		deleteTx(*it, properties);
		it = overlap.erase(it);
		it--;
		continue;
	    }
	    // create txNew
	    Transcript* txNew = createCopyOf(*it, properties, overlap);

	    string tempFeatureBack = txNew->exon_list.back().feature;

	    // delete problem exon
	    txNew->exon_list.pop_back();

	    // change some values of new transcripts for adaption
	    if (tempFeatureBack == "CDS"){
		if (txNew->strand == '-'){
		    txNew->tl_complete.first = false;
		    txNew->tis = txNew->exon_list.back().to;
		}else{
		    txNew->tl_complete.second = false;
		    txNew->stop_list.clear();
		    txNew->tes = txNew->exon_list.back().to;
		}
	    }else{
		if (txNew->strand == '-'){
		    txNew->tx_complete.first = false;
		    txNew->tss = -1;
		}else{
		    txNew->tx_complete.second = false;
		    txNew->tts = -1;
		}
	    }
	    txNew->exon_list.back().boundaryProblem = FREED;

	}else if ((*it)->exon_list.front().boundaryProblem >= 2 && (*it)->exon_list.back().boundaryProblem >= 2){
	    // erase the transcript, if it has problem boundary exons only
	    if ((*it)->exon_list.size() <= 2  || ((*it)->exon_list.front().feature == "CDS" && countCDS <= 1) || ((*it)->exon_list.back().feature == "CDS" && countCDS <= 1) || ((*it)->exon_list.front().feature == "CDS" && (*it)->exon_list.back().feature == "CDS" && countCDS <= 2)){
		deleteTx(*it, properties);
		it = overlap.erase(it);
		it--;
		continue;
	    }
	    Transcript* txNew = createCopyOf(*it, properties, overlap);

	    string tempFeatureFront = txNew->exon_list.front().feature;

	    // delete problem exon
	    txNew->exon_list.pop_front();

	    // change some values of new transcripts for adaption
	    if (tempFeatureFront == "CDS"){
		if (txNew->strand == '+'){
		    txNew->tl_complete.first = false;
		    txNew->tis = txNew->exon_list.front().from;
		}else{
		    txNew->tl_complete.second = false;
		    txNew->stop_list.clear();
		    txNew->tes = txNew->exon_list.front().from;
		}
	    }else{
		if (txNew->strand == '+'){
		    txNew->tx_complete.first = false;
		    txNew->tss = -1;
		}else{
		    txNew->tx_complete.second = false;
		    txNew->tts = -1;
		}
	    }
	    txNew->exon_list.front().boundaryProblem = FREED;

	    string tempFeatureBack = txNew->exon_list.back().feature;

	    // delete problem exon
	    txNew->exon_list.pop_back();

	    // change some values of new transcripts for adaption
	    if (tempFeatureBack == "CDS"){
		if (txNew->strand == '-'){
		    txNew->tl_complete.first = false;
		    txNew->tis = txNew->exon_list.back().to;
		}else{
		    txNew->tl_complete.second = false;
		    txNew->stop_list.clear();
		    txNew->tes = txNew->exon_list.back().to;
		}
	    }else{
		if (txNew->strand == '-'){
		    txNew->tx_complete.first = false;
		    txNew->tss = -1;
		}else{
		    txNew->tx_complete.second = false;
		    txNew->tts = -1;
		}
	    }
	    txNew->exon_list.back().boundaryProblem = FREED;
	}
    }
}

void search_n_destroy_doublings(list<Transcript*> &overlap, Properties &properties, bool abInitio){
    // delete all transcripts that are completly part of another transcript (in particular all exons are also in the other transcript and vice versa); the one with the lesser priority will be deleted
    if (overlap.size() > 1){
	for (list<Transcript*>::iterator it = overlap.begin(); it != overlap.end(); it++){
	    list<Transcript*>::iterator it_temp = it;
	    it_temp++;
	    for (list<Transcript*>::iterator itInside = it_temp; itInside != overlap.end(); itInside++){
		if (overlap.size() <= 1){return;}
		pair<bool,bool> who_is_part = is_part_of(*it, *itInside);		// <a is in b; b is in a>; if both are true, they are identical
		if (who_is_part.first == true){
		    if (who_is_part.second == true){
			if (abInitio){							// fill the supporter list only while the first run of this method
			    (*it)->supporter.push_front((*itInside)->originalId);
			    (*itInside)->supporter.push_front((*it)->originalId);
			}

			if (compare_quality((*it), *itInside)){
			    deleteTx((*itInside), properties);
			    itInside = overlap.erase(itInside);
			    itInside--;
			    
			}else{
			    deleteTx((*it), properties);
			    it = overlap.erase(it);
			    itInside = it;
			}

			// calculate minimal distance to the closest pred range boundary
			/*int mindist = -1;
			int mindistInside = -1;
			if ((*it)->pred_range.first && (*it)->pred_range.second){
			    mindist = min((*it)->exon_list.front().rangeToBoundary, (*it)->exon_list.back().rangeToBoundary);
			}
			if ((*itInside)->pred_range.first && (*itInside)->pred_range.second){
			    mindistInside = min((*itInside)->exon_list.front().rangeToBoundary, (*itInside)->exon_list.back().rangeToBoundary);
			}

			// delete a transcript if: less trustful (errordistance), lower priority, smaller score
			if ((mindist == -1 || mindist > properties.errordistance) && !(mindistInside == -1 || mindistInside > properties.errordistance)){
			    itInside = overlap.erase(itInside);
			    itInside--;
			}else if ((mindistInside == -1 || mindistInside > properties.errordistance) && !(mindist == -1 || mindist > properties.errordistance)){
			    it = overlap.erase(it);
			    itInside = it;
			}else{		// both in our both not in the error distance
			    if ((*it)->priority < (*itInside)->priority){
				it = overlap.erase(it);
				itInside = it;
			    }else if ((*it)->priority > (*itInside)->priority){
			        itInside = overlap.erase(itInside);
			        itInside--;
			    }else{	// priorities equal
				if (mindist < mindistInside){
				    it = overlap.erase(it);
				    itInside = it;
				}else{
				    itInside = overlap.erase(itInside);
				    itInside--;
				}
			    }
			}*/
		    }else{
			if (abInitio){
			    (*itInside)->supporter.push_front((*it)->originalId);
			}
			//it = overlap.erase(it);		// deleting throws information away
			//itInside = it;
		    }
		}else{
		    if (who_is_part.second == true){
			if (abInitio){
			    (*it)->supporter.push_front((*itInside)->originalId);
			}
			//itInside = overlap.erase(itInside);	// deleting throws information away
			//itInside--;
		    }
		}
	    }
	}
    }
}

void search_n_destroy_parts(list<Transcript*> &overlap, Properties &properties){
    // delete all transcripts that are part of another transcript
    if (overlap.size() > 1){
	for (list<Transcript*>::iterator it = overlap.begin(); it != overlap.end(); it++){
	    list<Transcript*>::iterator it_temp = it;
	    it_temp++;
	    for (list<Transcript*>::iterator itInside = it_temp; itInside != overlap.end(); itInside++){
		if (overlap.size() <= 1){return;}
		pair<bool,bool> who_is_part = is_part_of(*it, *itInside);

		if (who_is_part.first == true){
		    if (who_is_part.second == true){
		    }else{
			deleteTx((*it), properties);
			it = overlap.erase(it);
			itInside = it;
		    }
		}else{
		    if (who_is_part.second == true){
			deleteTx((*itInside), properties);
			itInside = overlap.erase(itInside);
			itInside--;
		    }
	    	}
	    }
	}
    }
}

void join(list<Transcript*> &overlap, char side, Properties &properties){
    // devides an overlap in a list of start/stop codon donors and in an other list of start/stop codon acceptors and joins every pair (one of each list) if they are combinable
    // list<Transcript*> new_overlap_part;
    list<Transcript*> donor;
    list<Transcript*> acceptor;

    for (list<Transcript*>::iterator it = overlap.begin(); it != overlap.end(); it++){
	if ((*it)->isNotFrameCorrect){continue;}
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
		bool frontSide = ((*it)->strand == '+' && side == '5') || ((*it)->strand == '-' && side == '3');	// closer to the first base on the dna strand (NOT the start codon position!)
		int fittingCase = isCombinable(*it, *it_donor, frontSide, properties);
		if (fittingCase){
		    Transcript* txNew = createCopyOf((*it), properties, overlap);
		    joining(*it_donor, (*it)->strand, txNew, fittingCase, properties);
		}
	    }
	}
    }

    // UTR joining:
    list<Transcript*> donorUTR;
    list<Transcript*> acceptorUTR;
    for (list<Transcript*>::iterator it = overlap.begin(); it != overlap.end(); it++){
	if (side == '3'){
	    if ((*it)->tx_complete.second){
		donorUTR.push_front(*it);
	    }else if ((*it)->tl_complete.second){
		acceptorUTR.push_front(*it);
	    }
	}else if (side == '5'){
	    if ((*it)->tx_complete.first){
		donorUTR.push_front(*it);
	    }else if ((*it)->tl_complete.first){
		acceptorUTR.push_front(*it);
	    }
	}

    }

    for (list<Transcript*>::iterator it = acceptorUTR.begin(); it != acceptorUTR.end(); it++){
	for (list<Transcript*>::iterator it_donor = donorUTR.begin(); it_donor != donorUTR.end(); it_donor++){
	    if ((*it)->strand == (*it_donor)->strand){
		bool frontSide = ((*it)->strand == '+' && side == '5') || ((*it)->strand == '-' && side == '3');
		int fittingCase = 0;
		if (frontSide){
		    if (((*it)->strand == '+' && (*it)->tis == (*it_donor)->tis) || ((*it)->strand == '-' && (*it)->tes == (*it_donor)->tes)){

			Transcript temp = *(*it);
			Transcript tempD = *(*it_donor);
			bool eraseExon = false;
			for (list<Exon>::iterator iti = temp.exon_list.begin(); iti != temp.exon_list.end(); iti++){
			    if ((*iti).feature == "CDS"){
				eraseExon = true;
			    }
			    if (eraseExon){
				temp.exon_list.erase(iti);
				iti--;
			    }
			}
			eraseExon = false;
			for (list<Exon>::iterator iti = tempD.exon_list.begin(); iti != tempD.exon_list.end(); iti++){
			    if ((*iti).feature == "CDS"){
				eraseExon = true;
			    }
			    if (eraseExon){
				tempD.exon_list.erase(iti);
				iti--;
			    }
			}
			if ((*it)->strand == '+'){
			    temp.tes = (*it_donor)->tes;
			}else{
			    temp.tis = (*it_donor)->tis;
			}
			pair<bool,bool> partBool;
			if (temp.exon_list.size() > 0){
			    partBool = is_part_of(&temp, &tempD);
			}else{
			    partBool.first = true;
			    partBool.second = true;
			}
			if (partBool.first){
			    fittingCase = 4;
			}
		    }
		}else{
		    if (((*it)->strand == '-' && (*it)->tis == (*it_donor)->tis) || ((*it)->strand == '+' && (*it)->tes == (*it_donor)->tes)){

			Transcript temp = *(*it);
			Transcript tempD = *(*it_donor);
			string lastFeature = "";
			for (list<Exon>::iterator iti = temp.exon_list.begin(); iti != temp.exon_list.end(); iti++){
			    if (lastFeature == "CDS" && (*iti).feature != "CDS"){
				break;
			    }
			    lastFeature = (*iti).feature;
			    temp.exon_list.erase(iti);
			    iti--;
			}
			lastFeature = "";
			for (list<Exon>::iterator iti = tempD.exon_list.begin(); iti != tempD.exon_list.end(); iti++){
			    if (lastFeature == "CDS" && (*iti).feature != "CDS"){
				break;
			    }
			    lastFeature = (*iti).feature;
			    tempD.exon_list.erase(iti);
			    iti--;
			}
			if ((*it)->strand == '+'){
			    temp.tis = (*it_donor)->tis;
			}else{
			    temp.tes = (*it_donor)->tes;
			}
			pair<bool,bool> partBool;
			if (temp.exon_list.size() > 0){
			    partBool = is_part_of(&temp, &tempD);
			}else{
			    partBool.first = true;
			    partBool.second = true;
			}

			if (partBool.first){
			    fittingCase = 2;
			}
		    }
		}
		if (fittingCase){
		    Transcript* txNew = createCopyOf((*it), properties, overlap);
		    joining(*it_donor, (*it)->strand, txNew, fittingCase, properties);
		}
	    }
	}
    }
}

void joining(Transcript* t2, char strand, Transcript* txNew, int fittingCase, Properties &properties){

    // joins transcripts in one direction so that every suitable exon will be transferred and returns a new "joined" transcript without deleting the old ones

    int lastPositionInOriginal = txNew->exon_list.back().to;
    int firstPositionInOriginal = txNew->exon_list.front().from;
    string lastPosfeature = txNew->exon_list.back().feature;
    string firstPosfeature = txNew->exon_list.front().feature;

    int nrOfJoinedCDS = 0;
    int nrOfCoreCDS = 0;
    bool inside = false;
    for (list<Exon>::iterator it = txNew->exon_list.begin(); it != txNew->exon_list.end(); it++){
	if (!inside && txNew->outerCds.first && txNew->outerCds.first->from == (*it).from && txNew->outerCds.first->to == (*it).to && txNew->outerCds.first->frame == (*it).frame){
	    inside = true;
	}
	if (inside){
	    nrOfCoreCDS++;
	    if (txNew->outerCds.second->from == (*it).from && txNew->outerCds.second->to == (*it).to && txNew->outerCds.second->frame == (*it).frame){break;}
	}
    }

    bool are_at_add_part;
    list<Exon> temp_exon_list;
    temp_exon_list.clear();
    switch (fittingCase){
    case 0:
	cerr << "WARNING: shouldnt happen (in joining())!" << endl;
	break;
    case 1:
	if (strand == '+'){
	    if (txNew->stop_list.empty()){txNew->stop_list = t2->stop_list;}else{cerr << "WARNING: existing stop list should be replaced (in joining())!1" << endl;}
	    txNew->tts = t2->tts;
	    txNew->tx_complete.second = t2->tx_complete.second;
	    txNew->tes = t2->tes;
	    txNew->tl_complete.second = t2->tl_complete.second;
	    txNew->joinpartner.second = t2->originalId;
	}else{
	    txNew->tss = t2->tss;
	    txNew->tx_complete.first = t2->tx_complete.first;
	    txNew->tis = t2->tis;
	    txNew->tl_complete.first = t2->tl_complete.first;
	    txNew->joinpartner.first = t2->originalId;
	}
	are_at_add_part = false;
	for (list<Exon>::iterator it = t2->exon_list.begin(); it != t2->exon_list.end(); it++){
	    if (are_at_add_part){
		txNew->exon_list.push_back(*it);
		if ((*it).feature == "CDS"){nrOfJoinedCDS++;}
	    }
	    if (lastPositionInOriginal >= (*it).from && lastPositionInOriginal <= (*it).to){
		txNew->exon_list.back().to = (*it).to;
    		txNew->exon_list.back().score = (*it).score;
		txNew->exon_list.back().rangeToBoundary = (*it).rangeToBoundary;
		if (strand == '-'){txNew->exon_list.back().frame = (*it).frame;}
		are_at_add_part = true;
	    }
	}
	break;
    case 2:
	if (strand == '+'){
	    if (txNew->stop_list.empty()){txNew->stop_list = t2->stop_list;}else{/*cerr << "WARNING: existing stop list should be replaced (in joining())!2" << endl;*/}
	    txNew->tts = t2->tts;
	    txNew->tx_complete.second = t2->tx_complete.second;
	    txNew->tes = t2->tes;
	    txNew->tl_complete.second = t2->tl_complete.second;
	    txNew->joinpartner.second = t2->originalId;
	}else{
	    txNew->tss = t2->tss;
	    txNew->tx_complete.first = t2->tx_complete.first;
	    txNew->tis = t2->tis;
	    txNew->tl_complete.first = t2->tl_complete.first;
	    txNew->joinpartner.first = t2->originalId;
	}
	are_at_add_part = false;
	for (list<Exon>::iterator it = t2->exon_list.begin(); it != t2->exon_list.end(); it++){
	    if (lastPositionInOriginal <= ((*it).from - (int) properties.minimumIntronLength) || ( (*it).feature != "CDS" && lastPosfeature == "CDS" && lastPositionInOriginal <= ((*it).from - 1) )){
		are_at_add_part = true;
	    }
	    if (are_at_add_part){
		txNew->exon_list.push_back(*it);
		if ((*it).feature == "CDS"){nrOfJoinedCDS++;}
	    }
	}
	break;
    case 3:
	if (strand == '+'){
	    txNew->tss = t2->tss;
	    txNew->tx_complete.first = t2->tx_complete.first;
	    txNew->tis = t2->tis;
	    txNew->tl_complete.first = t2->tl_complete.first;
	    txNew->joinpartner.first = t2->originalId;
	}else{
	    if (txNew->stop_list.empty()){txNew->stop_list = t2->stop_list;}else{cerr << "WARNING: existing stop list should be replaced (in joining())!3" << endl;}
	    txNew->tts = t2->tts;
	    txNew->tx_complete.second = t2->tx_complete.second;
	    txNew->tes = t2->tes;
	    txNew->tl_complete.second = t2->tl_complete.second;
	    txNew->joinpartner.second = t2->originalId;
	}
	for (list<Exon>::iterator it = t2->exon_list.begin(); it != t2->exon_list.end(); it++){
	    if (firstPositionInOriginal >= (*it).from && firstPositionInOriginal <= (*it).to){
		txNew->exon_list.front().from = (*it).from;
    		txNew->exon_list.front().score = (*it).score;
		txNew->exon_list.front().rangeToBoundary = (*it).rangeToBoundary;
		if (strand == '+'){txNew->exon_list.front().frame = (*it).frame;}
		txNew->exon_list.merge(temp_exon_list);
		break;
	    }
	    temp_exon_list.push_back(*it);
	    if ((*it).feature == "CDS"){nrOfJoinedCDS++;}
	}
	break;
    case 4:
	if (strand == '+'){
	    txNew->tss = t2->tss;
	    txNew->tx_complete.first = t2->tx_complete.first;
	    txNew->tis = t2->tis;
	    txNew->tl_complete.first = t2->tl_complete.first;
	    txNew->joinpartner.first = t2->originalId;
	}else{
	    if (txNew->stop_list.empty()){txNew->stop_list = t2->stop_list;}else{/*cerr << "WARNING: existing stop list should be replaced (in joining())!4" << endl;*/}
	    txNew->tts = t2->tts;
	    txNew->tx_complete.second = t2->tx_complete.second;
	    txNew->tes = t2->tes;
	    txNew->tl_complete.second = t2->tl_complete.second;
	    txNew->joinpartner.second = t2->originalId;
	}
	for (list<Exon>::iterator it = t2->exon_list.begin(); it != t2->exon_list.end(); it++){
	    if ((!((*it).feature != "CDS" && firstPosfeature == "CDS") && firstPositionInOriginal < ((*it).to + (int) properties.minimumIntronLength)) || ( (*it).feature != "CDS" && firstPosfeature == "CDS" && firstPositionInOriginal < ((*it).to + 1) )){
		txNew->exon_list.merge(temp_exon_list);
		break;
	    }
	    temp_exon_list.push_back(*it);
	    if ((*it).feature == "CDS"){nrOfJoinedCDS++;}
	}
	break;
    default:
	cerr << "WARNING: unexpected case (in joining())!" << endl;
    }

    if (nrOfJoinedCDS > nrOfCoreCDS){
	txNew->priority = t2->priority;
    }
}

int is_combinable(Transcript const* t1, Transcript const* t2, char strand, char side, Properties &properties){
    // is true,	if the first exon which is minimal the minimumIntronLength away from the last exon of the other transcript in the appropriate direction
    // 			&& the exons at these positions are frame-compatible
    // 			&& the transcripts are overlapping					// maybe we can improve something here, that combinable non-overlaping transcripts gets true (but be carefull)

    if ((max( (*t1).tis,(*t1).tes) < min((*t2).tis,(*t2).tes)) || (min((*t1).tis,(*t1).tes) > max((*t2).tis,(*t2).tes))){
	return 0;
    }
    if ((strand == '+' && side == '3') || (strand == '-' && side == '5')){
	for (list<Exon>::const_iterator it = t2->exon_list.begin(); it != t2->exon_list.end(); it++){
	    if (((*t1).exon_list.back().to >= (*it).from) && ((*t1).exon_list.back().to <= (*it).to)){
		if ((*t1).exon_list.back().frame == -1 && (*it).frame == -1){return 1;}
		else if (((*t1).exon_list.back().frame != -1 && (*it).frame == -1) || ((*t1).exon_list.back().frame == -1 && (*it).frame != -1)){return 0;}
		if ((strand == '+') && ((3 - (((*t1).exon_list.back().to - (*t1).exon_list.back().from) - (*t1).exon_list.back().frame) % 3) % 3) == ((3 - (((*t1).exon_list.back().to - (*it).from) - (*it).frame) % 3) % 3) ){
		    return 1;
		}else if ((strand == '-') && (*t1).exon_list.back().frame == ((3 - (((*it).to - (*t1).exon_list.back().to) - (*it).frame) % 3) % 3) ){
		    return 1;
		}else{
		    return 0;
		}
	    }else{
		if ((*t1).exon_list.back().to <= ((*it).from - (int) properties.minimumIntronLength)){
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
		if ((*t1).exon_list.front().frame == -1 && (*it).frame == -1){return 3;}
		else if (((*t1).exon_list.back().frame != -1 && (*it).frame == -1) || ((*t1).exon_list.back().frame == -1 && (*it).frame != -1)){return 0;}
		if ( (strand == '-') && ((3 - (((*t1).exon_list.front().to - (*t1).exon_list.front().from) - (*t1).exon_list.front().frame) % 3) % 3) == ((3 - (((*it).to - (*t1).exon_list.front().from) - (*it).frame) % 3) % 3) ){
		    return 3;
		}else if ( (strand == '+') && (*t1).exon_list.front().frame == ((3 - (((*t1).exon_list.front().from - (*it).from) - (*it).frame) % 3) % 3) ){
		    return 3;
		}else{
		    return 0;
		}
	    }else{
		if ((*t1).exon_list.front().from < ((*it).to + (int) properties.minimumIntronLength)){
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

pair<bool,bool> is_part_of(Transcript const* t1, Transcript const* t2){
    // is true,false or false,true if one transcript contains the other completely
    // is true,true if the transcripts are equal in exons			// this case could completly replace compare_transcripts
    int UTRtoleranceWindow = 10;
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
	    if((*it1).to == (*it2).to && (*it1).frame == (*it2).frame){		// same start and same stop of exon
		it1++;
		it2++;
	    }
	    else{								// same start but diffrent stop of exon
		if ((*it1).feature != "CDS" && (*it2).feature != "CDS" && (*it1).to <= (*it2).to+UTRtoleranceWindow && (*it1).to >= (*it2).to-UTRtoleranceWindow){
		    it1++;
		    it2++;
		}else{
		    t1_is_part = false;
		    t2_is_part = false;
		    break;
		}
	    }
	}else if ((*it1).to == (*it2).to){					// different start but same stop of exon
	    if ((*it1).feature != "CDS" && (*it2).feature != "CDS" && (*it1).from <= (*it2).from+UTRtoleranceWindow && (*it1).from >= (*it2).from-UTRtoleranceWindow){
		it1++;
		it2++;
	    }else{
		t1_is_part = false;
		t2_is_part = false;
		break;
	    }
	}else if ((*it1).from > (*it2).from){					// different start and different stop of exon AND exon1 is behind exon2
	    if ((*it1).feature != "CDS" && (*it2).feature != "CDS" && (*it1).from <= (*it2).from+UTRtoleranceWindow && (*it1).from >= (*it2).from-UTRtoleranceWindow && (*it1).to <= (*it2).to+UTRtoleranceWindow && (*it1).to >= (*it2).to-UTRtoleranceWindow){
		it1++;
	    }else{
		t2_is_part = false;
	    }
	    it2++;
	}else{									// different start and different stop of exon AND exon1 is in front of exon2
	    if ((*it1).feature != "CDS" && (*it2).feature != "CDS" && (*it1).from <= (*it2).from+UTRtoleranceWindow && (*it1).from >= (*it2).from-UTRtoleranceWindow && (*it1).to <= (*it2).to+UTRtoleranceWindow && (*it1).to >= (*it2).to-UTRtoleranceWindow){
		it2++;
	    }else{
		t1_is_part = false;
	    }
	    it1++;
	}
	if (!t1_is_part && !t2_is_part){break;}					// both transcripts are not part of the other one
	if (it1 == t1->exon_list.end() && it2 == t2->exon_list.end()){		// no more exon in both exon lists
	    break;
	}else{
	    if (it1 == t1->exon_list.end() && !(it2 == t2->exon_list.end())){	// no more exon in exon1 list
		t2_is_part = false;
		break;
	    }
	    if (it2 == t2->exon_list.end() && !(it1 == t1->exon_list.end())){	// no more exon in exon2 list
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

void eval_gtf(list<Transcript*> &overlap, Properties &properties){
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
	search_n_destroy_doublings(annotation, properties, false);
	if ((*it)->pred_range.first && (*it)->pred_range.second){
	    prediction.push_back(*it);
	}
	search_n_destroy_doublings(prediction, properties, false);
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
    cout << "Original ID: " << tx->originalId << " from " << tx->inputFile << endl;
    cout << "Priority: " << tx->priority << endl;
    cout << tx->exon_list.size() << " elements, complete start: " << tx->tl_complete.first << ", complete stop: " << tx->tl_complete.second << endl;
    cout << "start codon of: ";
    if (!tx->joinpartner.first.empty()){cout << tx->joinpartner.first;}else{cout << "noone";}
    cout << ", stop codon of: ";
    if (!tx->joinpartner.second.empty()){cout << tx->joinpartner.second;}else{cout << "noone";}
    cout << endl;
    if (tx->pred_range.first && tx->pred_range.second){
	cout << "Pred_range: " << tx->pred_range.first << "\t" << tx->pred_range.second << endl;
	cout << "Pred_range EXON: " << tx->exon_list.front().predRange.first << "\t" << tx->exon_list.back().predRange.second << endl;
	cout << "Minimum distance from pred_range borders: " << min(tx->exon_list.front().from - tx->pred_range.first, tx->pred_range.second - tx->exon_list.back().to) << endl;
	cout << "Boundary problem: " << tx->exon_list.front().boundaryProblem << " " << tx->exon_list.back().boundaryProblem << endl;
	cout << "Distances from pred_range borders: " << tx->exon_list.front().rangeToBoundary << " " << tx->exon_list.back().rangeToBoundary << endl;
    }
    double x = 0.0;
    for (list<Exon>::const_iterator it = tx->exon_list.begin(); it != tx->exon_list.end(); it++){
	cout << (*it).chr << "\tTESTOUTPUT\t" << (*it).feature << "\t" << (*it).from << "\t" << (*it).to << "\t" << (*it).score << "\t" << tx->strand << "\t" << (*it).frame << "\ttranscript_id \"" << tx->t_id << "\"; gene_id \"" << tx->parent->g_id << "\";" << endl;
	x += (*it).score;
    }
    x /= tx->exon_list.size();
    cout << "Mean score: " << x << endl;
    cout << "Exon_list end----------------------------------------------------------------" << endl;
}

void calculatePredictionScore(Transcript* tx){
    float cumScoreCDS = 0;
    float countCDS = 0;
    float cumScoreUTR = 0;
    float countUTR = 0;
    float cumScoreIntron = 0;
    float countIntron = 0;

    int maxIntronSize = 0;

    for (list<Exon>::const_iterator iex = tx->exon_list.begin(); iex != tx->exon_list.end(); iex++){


	list<Exon>::const_iterator iexNext = iex;
	iexNext++;
	if (iexNext != tx->exon_list.end()){
	    int intronSize = (*iexNext).from - (*iex).to - 1;
	    if (maxIntronSize < intronSize){maxIntronSize = intronSize;}
	    if (intronSize > 0){
	        countIntron++;
	    }
	}


	if (iex->feature == "CDS"){
	    cumScoreCDS += iex->score;
	    countCDS++;
	    if (iex->score < 0 && iex->score > 1){
		cerr << "WARNING: A value of column 6 is above 1 (100%) or below 0 (0%)." << endl;
	    }
	}else if (iex->feature == "UTR"){
	    cumScoreUTR += iex->score;
	    countUTR++;
	    if (iex->score < 0 && iex->score > 1){
		cerr << "WARNING: A value of column 6 is above 1 (100%) or below 0 (0%)." << endl;
	    }
	}
    }

    for (list<Exon>::const_iterator iin = tx->intron_list.begin(); iin != tx->intron_list.end(); iin++){
	cumScoreIntron += iin->score;
	//countIntron++;
	if (iin->score < 0 && iin->score > 1){
	    cerr << "WARNING: A value of column 6 is above 1 (100%) or below 0 (0%)." << endl;
	}
    }

    if (countCDS == 0){
	cumScoreCDS = 0;
    }else{
	cumScoreCDS /= countCDS;
    }
    if (countUTR == 0){
	cumScoreUTR = 0;
    }else{
	cumScoreUTR /= countUTR;
    }
    if (countIntron == 0){
	cumScoreIntron = 0;
    }else{
	cumScoreIntron /= countIntron;
    }

    tx->predictionScore = (4*cumScoreCDS + 4*cumScoreIntron + cumScoreUTR)/9;
}

bool compare_quality(Transcript const* lhs, Transcript const* rhs){
    // completeness
    if ((lhs->tl_complete.first && lhs->tl_complete.second) && !(rhs->tl_complete.first && rhs->tl_complete.second)){
	return true;
    }
    if (!(lhs->tl_complete.first && lhs->tl_complete.second) && (rhs->tl_complete.first && rhs->tl_complete.second)){
	return false;
    }

    if ((lhs->tl_complete.first || lhs->tl_complete.second) && !(rhs->tl_complete.first || rhs->tl_complete.second)){
	return true;
    }
    if (!(lhs->tl_complete.first || lhs->tl_complete.second) && (rhs->tl_complete.first || rhs->tl_complete.second)){
	return false;
    }

    // closeness to prediction range border
    if ((lhs->exon_list.front().boundaryProblem <= 1 && lhs->exon_list.back().boundaryProblem <= 1) && !(rhs->exon_list.front().boundaryProblem <= 1 && rhs->exon_list.back().boundaryProblem <= 1)){
	return true;
    }
    if (!(lhs->exon_list.front().boundaryProblem <= 1 && lhs->exon_list.back().boundaryProblem <= 1) && (rhs->exon_list.front().boundaryProblem <= 1 && rhs->exon_list.back().boundaryProblem <= 1)){
	return false;
    }

    // priority
    if (lhs->priority > rhs->priority){
	return true;
    }
    if (lhs->priority < rhs->priority){
	return false;
    }

    // prediction score
    if (lhs->predictionScore >= rhs->predictionScore){
	return true;
    }else{
	return false;
    }
}

/*bool isCopyOfSameOriginal(Transcript* t1,Transcript* t2){
    Transcript* tmp1 = t1;
    while (tmp1->isCopyOf != NULL){
	tmp1 = tmp1->isCopyOf;
    }
    Transcript* tmp2 = t2;
    while (tmp2->isCopyOf != NULL){
	tmp2 = tmp2->isCopyOf;
    }
    if (tmp1->t_id == tmp2->t_id){
	return true;
    }else{
	return false;
    }
}*/

bool shareAlternativeVariant(Gene* g1, Gene* g2){
    for (list<Transcript*>::iterator it = g1->children.begin(); it != g1->children.end(); it++){
	for (list<Transcript*>::iterator itInside = g2->children.begin(); itInside != g2->children.end(); itInside++){
	    if (alternativeVariants((*it), (*itInside))){
		return true;
	    }
	}
    }
    return false;
}

bool alternativeVariants(Transcript* t1, Transcript* t2){
    if (overlappingCdsWithAnything(t1, t2) && t1->strand == t2->strand && (t1->hasCommonExon(t2) || (t1->hasCommonTlStart(t2) || t1->hasCommonTlStop(t2)))){
	return true;
    }
    return false;
}

void eukaSelectionDevelopment(list<Transcript*> &overlap, Properties &properties){

    search_n_destroy_doublings(overlap, properties, false);

    // delete transcripts that belong to input files with priorities that sho‎uld be suppressed
    for (list<Transcript*>::iterator it = overlap.begin(); it != overlap.end(); it++){
	if (find(properties.supprList.begin(),properties.supprList.end(),(*it)->priority) != properties.supprList.end()){
	    deleteTx(*it, properties);
	    it = overlap.erase(it);
	    it--;
	    continue;
	}
    	// calculate score to compare the quality of transcripts
	calculatePredictionScore(*it);
    }

    overlap.sort(compare_quality);

    // figure out, which genes are present in this overlap (#genes <= #transcripts)
    list<Gene*> genes;
    for (list<Transcript*>::iterator it = overlap.begin(); it != overlap.end(); it++){
	if (find(genes.begin(),genes.end(),(*it)->parent) == genes.end()){
	    genes.push_back((*it)->parent);
	}
    }

    // delete uncomplete transcripts from gene, if there is a complete one
    for (list<Gene*>::iterator it = genes.begin(); it != genes.end(); it++){
	bool completeExists = false;
	for (list<Transcript*>::iterator tit = (*it)->children.begin(); tit != (*it)->children.end(); tit++){
	    if ((*tit)->tl_complete.first && (*tit)->tl_complete.second){
		completeExists = true;
		break;
	    }
	}
	if (completeExists){
	    for (list<Transcript*>::iterator tit = (*it)->children.begin(); tit != (*it)->children.end(); tit++){
		if ( !( (*tit)->tl_complete.first && (*tit)->tl_complete.second ) ){
		    for (list<Transcript*>::iterator oit = overlap.begin(); oit != overlap.end(); oit++){
			if (*tit == *oit){
			    oit = overlap.erase(oit);
			    oit--;
			}
		    }
		    deleteTx(*tit, properties);
		    tit--;
		}
	    }
	}
    }

    // if flag alternatives is on, try to combine transcripts to one gene, if they look like alternative spliced variants
    if (properties.alternatives && genes.size() >= 2){


/*	for (list<Gene*>::iterator it = genes.begin(); it != genes.end(); it++){

	    list<Gene*>::iterator itInside = it;
	    itInside++;
	    while (itInside != genes.end()){
		if (shareAlternativeVariant((*it), (*itInside))){

		    // Transfer all transcripts of itInside to it and delete itInside
		    for (list<Transcript*>::iterator tit = (*itInside)->children.begin(); tit != (*itInside)->children.end(); tit++){
			(*tit)->parent = (*it);
			(*it)->children.push_back(*tit);
			(*it)->nrOfTx += (*itInside)->nrOfTx;
			(*it)->nrOfPrintedTx += (*itInside)->nrOfPrintedTx;
			tit = (*itInside)->children.erase(tit);
			tit--;
		    }
		    deleteGene((*itInside), properties);
		    itInside = genes.erase(itInside);
		    continue;
		}
		itInside++;
	    }

	    // sort children 
	    (*it)->children.sort(compare_quality);
	}*/


	for (list<Gene*>::iterator it = genes.begin(); it != genes.end(); it++){

	    list<Gene*>::iterator itInside = it;
	    itInside++;

	    while (itInside != genes.end()){

		if (shareAlternativeVariant((*it), (*itInside))){

		    // Transfer all transcripts of itInside to it and delete itInside
		    for (list<Transcript*>::iterator tit = (*it)->children.begin(); tit != (*it)->children.end(); tit++){
			(*tit)->parent = (*itInside);
			(*itInside)->children.push_back(*tit);
			(*itInside)->nrOfTx += (*it)->nrOfTx;
			(*itInside)->nrOfPrintedTx += (*it)->nrOfPrintedTx;
			tit = (*it)->children.erase(tit);
			tit--;
		    }
		    deleteGene((*it), properties);
		    it = genes.erase(it);
		    it--;
		    break;
		}
		itInside++;
	    }
	}
	// sort children
	for (list<Gene*>::iterator it = genes.begin(); it != genes.end(); it++){ 
	    (*it)->children.sort(compare_quality);
	}
    }

    // look for indirect boundary problems
    for (list<Transcript*>::iterator it = overlap.begin(); it != overlap.end(); it++){
	list<Transcript*>::iterator itInside = it;
	itInside++;
	while (itInside != overlap.end()){
	    if ((*it)->parent == (*itInside)->parent){

//cout << "HIER: " << (*it)->exon_list.back().boundaryProblem << " " << (*itInside)->exon_list.back().boundaryProblem << " " << (*it)->exon_list.back().indirectBoundProbEnemies.front() << " " << (*itInside)->exon_list.back().indirectBoundProbEnemies.front()->originalId << endl;

		bool b1 = (*it)->indirectBoundaryProblem(*itInside);
		bool b2 = (*itInside)->indirectBoundaryProblem(*it);

//cout << "RESULT: " << b1 << " " << b2 << " " << (*it)->originalId << " " << (*itInside)->originalId << " " << (*it)->exon_list.back().feature << " " << (*itInside)->exon_list.back().feature << endl;

		if ( b1 && !b2 ){
		    // DELETE *it
		    deleteTx(*it, properties);
		    it = overlap.erase(it);
		    it--;
		}
		if ( b2 && !b1 ){
		    deleteTx(*itInside, properties);
		    itInside = overlap.erase(itInside);
		    itInside--;
		}
	    }
	    itInside++;
	}
    }

    // find out, which transcripts can coexist in a gene structure ("consistent")
    for (list<Transcript*>::iterator it = overlap.begin(); it != overlap.end(); it++){
	list<Transcript*>::iterator itInside = it;
	itInside++;
	while (itInside != overlap.end()){
	    if (!overlappingCdsWithAnything((*it), (*itInside)) || (*it)->strand != (*itInside)->strand || ((*it)->parent == (*itInside)->parent && (*it)->originalId != (*itInside)->originalId)){
	        (*it)->consistent.push_back((*itInside)->t_id);
		(*itInside)->consistent.push_back((*it)->t_id);
	    }
	    itInside++;
	}
    }

/*    list<string> openTx;
    unordered_map<string,Transcript*> txMap;
    for (list<Transcript*>::iterator it = overlap.begin(); it != overlap.end(); it++){
	openTx.push_back((*it)->t_id);
	txMap[(*it)->t_id] = (*it);
    }*/
//    list<list<Transcript*>> groups;
//    list<Transcript*> actGroup;
//    list<string> selectedTx;

//    recursiv(groups, actGroup, openTx, properties, txMap, selectedTx);

    //overlap = groups.front();
    overlap.sort(compare_quality);
    overlap = calcGeneStructur(overlap, properties);

}

list<Transcript*> calcGeneStructur(list<Transcript*> overlap, Properties &properties){
    list<Transcript*> newOvlp;
    list<string> openTx;

    for (list<Transcript*>::iterator it = overlap.begin(); it != overlap.end(); it++){
	openTx.push_back((*it)->t_id);
    }

    for (list<Transcript*>::iterator it = overlap.begin(); it != overlap.end(); it++){
	if ((find(openTx.begin(), openTx.end(), (*it)->t_id) != openTx.end())){
	    newOvlp.push_back(*it);
	    for (list<string>::iterator iti = openTx.begin(); iti != openTx.end(); iti++){
		if ((find((*it)->consistent.begin(), (*it)->consistent.end(), (*iti)) == (*it)->consistent.end())){
		    openTx.erase(iti);
		    iti--;
		}
	    }
	}
    }
    return newOvlp;
}

void recursiv(list<list<Transcript*>> &groups, list<Transcript*> actGroup, list<string> openTx, Properties &properties, unordered_map<string,Transcript*> &txMap, list<string> selectedTx){
    // start at every possible
    for (list<string>::iterator it = openTx.begin(); it != openTx.end(); it++){

	// if a transcript should be add to the group, which was already the start point (selectedTx), then we skip that case in order to get no redundant solutions
	if (find(selectedTx.begin(), selectedTx.end(), (*it)) != selectedTx.end()){
	    continue;
	}
	selectedTx.push_back(txMap[*it]->t_id);

	list<Transcript*> tmpGroup = actGroup;
	tmpGroup.push_back(txMap[*it]);
	list<string> openTemp;
	openTemp.clear();
	for (list<string>::iterator itInside = openTx.begin(); itInside != openTx.end(); itInside++){
	    if (find((txMap[*it])->consistent.begin(), (txMap[*it])->consistent.end(), (*itInside)) != (txMap[*it])->consistent.end()){
		openTemp.push_back(*itInside);
	    }
	}

	// tests, whether there are no transcripts in this overlap left that are consistent to all already taken transcripts
	if (openTemp.empty()){
	    groups.push_back(tmpGroup);
	}
	else{
	    recursiv(groups, tmpGroup, openTemp, properties, txMap, selectedTx);
	}
    }
}

bool overlapping(Transcript* t1, Transcript* t2){
    int t1TxFirst = t1->getTxStart();
    int t1TxLast = t1->getTxEnd();

    int t2TxFirst = t2->getTxStart();
    int t2TxLast = t2->getTxEnd();

    return !((t1TxLast < t2TxFirst) || (t1TxFirst > t2TxLast));
}

bool overlappingCdsWithAnything(Transcript* t1, Transcript * t2){

    int t1TlFirst = min(t1->tis,t1->tes);
    int t1TlLast = max(t1->tis,t1->tes);
    int t1TxFirst = t1->getTxStart();
    int t1TxLast = t1->getTxEnd();

    int t2TlFirst = min(t2->tis,t2->tes);
    int t2TlLast = max(t2->tis,t2->tes);
    int t2TxFirst = t2->getTxStart();
    int t2TxLast = t2->getTxEnd();

    return !((t1TlLast < t2TxFirst) || (t1TlFirst > t2TxLast)) || !((t1TxLast < t2TlFirst) || (t1TxFirst > t2TlLast));
}

bool overlappingCdsWithCds(Transcript* t1, Transcript* t2){
    return (((t2->tis - t1->tis) * (t2->tis - t1->tes)) < 0) || (((t2->tes - t1->tis) * (t2->tes - t1->tes)) < 0);
}

bool overlappingCdsOnlyWithUtr(Transcript* t1, Transcript* t2){
    return overlappingCdsWithAnything(t1,t2) && !overlappingCdsWithCds(t1, t2);
}

bool overlappingUtrOnly(Transcript* t1, Transcript* t2){
    return overlapping(t1, t2) && !overlappingCdsWithAnything(t1, t2);
}

int isCombinable(Transcript* t1, Transcript* t2, bool frontSide, Properties &properties){
// return 0 if t1 and t2 are not combinable (joinable), otherwise return the joining case
    // if not overlapping, they are not combinable
    if (!overlapping(t1, t2)){ return 0; }

    // backSide ("+" && "3'" and "-" && "5'")
    if (!frontSide){
	if ( ( (*t1).tes < (*t2).tis && (*t1).strand == '+' ) || ( (*t1).tes > (*t2).tis && (*t1).strand == '-' ) ){ return 0; }
	// for every exon in t2
	for (list<Exon>::const_iterator it = t2->exon_list.begin(); it != t2->exon_list.end(); it++){
	    if ((*it).feature != "CDS"){continue;}
	    // return 1: if t1.back() ends in an exon of t2 such that they are combinable; return 2 if t1.back() does not end ...
	    if (((*t1).exon_list.back().to >= (*it).from) && ((*t1).exon_list.back().to <= (*it).to)){
		if ((*t1).exon_list.back().frame == -1 && (*it).frame == -1){return 1;}
		else if (((*t1).exon_list.back().frame != -1 && (*it).frame == -1) || ((*t1).exon_list.back().frame == -1 && (*it).frame != -1)){return 0;}
		if ((t1->strand == '+') && ((3 - (((*t1).exon_list.back().to - (*t1).exon_list.back().from) - (*t1).exon_list.back().frame) % 3) % 3) == ((3 - (((*t1).exon_list.back().to - (*it).from) - (*it).frame) % 3) % 3) ){
		    return 1;
		}else if ((t1->strand == '-') && (*t1).exon_list.back().frame == ((3 - (((*it).to - (*t1).exon_list.back().to) - (*it).frame) % 3) % 3) ){
		    return 1;
		}else{
		    return 0;
		}
	    }else{
		if ((*t1).exon_list.back().to <= ((*it).from - (int) properties.minimumIntronLength)){
		    if ((*t1).exon_list.back().frame == -1 && (*it).frame == -1){return 2;}
		    else if (((*t1).exon_list.back().frame != -1 && (*it).frame == -1) || ((*t1).exon_list.back().frame == -1 && (*it).frame != -1)){return 0;}
		    if ((t1->strand == '+') && ((*it).frame == (3 - ( ((*t1).exon_list.back().to - (*t1).exon_list.back().from + 1) - (*t1).exon_list.back().frame) % 3) % 3)){
			return 2;
		    }else if ((t1->strand == '-') && ((*t1).exon_list.back().frame == (3 - ( ((*it).to - (*it).from + 1) - (*it).frame) % 3) % 3)){
			return 2;
		    }else{
			return 0;
		    }
		}
	    }
	}
    }else{		// frontSide
	if ( ((*t1).tis > (*t2).tes && (*t1).strand == '+') || ((*t1).tis < (*t2).tes && (*t1).strand == '-') ){ return 0; }
	for (list<Exon>::reverse_iterator it = t2->exon_list.rbegin(); it != t2->exon_list.rend(); it++){
	    if ((*it).feature != "CDS"){continue;}
	    // return 3: if t1.front() ends in an exon of t2 such that they are combinable; return 4 if t1.front() does not end ...

	    if (((*t1).exon_list.front().from >= (*it).from) && ((*t1).exon_list.front().from <= (*it).to)){
		if ((*t1).exon_list.front().frame == -1 && (*it).frame == -1){return 3;}
		else if (((*t1).exon_list.back().frame != -1 && (*it).frame == -1) || ((*t1).exon_list.back().frame == -1 && (*it).frame != -1)){return 0;}
		if ( (t1->strand == '-') && ((3 - (((*t1).exon_list.front().to - (*t1).exon_list.front().from) - (*t1).exon_list.front().frame) % 3) % 3) == ((3 - (((*it).to - (*t1).exon_list.front().from) - (*it).frame) % 3) % 3) ){
		    return 3;
		}else if ( (t1->strand == '+') && (*t1).exon_list.front().frame == ((3 - (((*t1).exon_list.front().from - (*it).from) - (*it).frame) % 3) % 3) ){
		    return 3;
		}else{
		    return 0;
		}
	    }else{
		if ((*t1).exon_list.front().from >= ((*it).to + (int) properties.minimumIntronLength)){
		    if ((*t1).exon_list.front().frame == -1 && (*it).frame == -1){return 4;}
		    else if (((*t1).exon_list.back().frame != -1 && (*it).frame == -1) || ((*t1).exon_list.back().frame == -1 && (*it).frame != -1)){return 0;}
/*		    if (it != t2->exon_list.begin()){
			it--;
		    }else{
			return 0;
		    }*/
		    if ((t1->strand == '-') && ((*it).frame == (3 - ( ((*t1).exon_list.front().to - (*t1).exon_list.front().from + 1) - (*t1).exon_list.front().frame) % 3) % 3)){
			return 4;
		    }else if ((t1->strand == '+') && ((*t1).exon_list.front().frame == (3 - ( ((*it).to - (*it).from + 1) - (*it).frame) % 3) % 3)){
			return 4;
		    }else{
			return 0;
		    }
		}
	    }
	}
    }
    return 0;
}

void testInputAlternatives(Properties &properties){
    for (auto pointer = (*properties.geneMap).begin(); pointer != (*properties.geneMap).end(); pointer++){
	if (pointer->second->children.size() <= 1){continue;}
	for (list<Transcript*>::iterator it = pointer->second->children.begin(); it != pointer->second->children.end(); it++){
	    bool x = false;
            if (!(*it)->tl_complete.first || !(*it)->tl_complete.second){continue;}

	    for (list<Transcript*>::iterator itI = pointer->second->children.begin(); itI != pointer->second->children.end(); itI++){
		if (!(*itI)->tl_complete.first || !(*itI)->tl_complete.second){continue;}
		if ((*it)->t_id == (*itI)->t_id){continue;}
		if (alternativeVariants((*it), (*itI))){
		    x = true;
		}
	    }
    	    if (!x){
	       cerr << "The transkript " << (*it)->originalId << " from file " << (*it)->inputFile << " is not an alternative variant in its gene in the definition of this program!" << endl;
	    }
	}
    }
}

void displayWarning(string const &warning, Properties &properties, string warningString){
    int unsigned n = 2;
    properties.warningCount[warningString]++;
    if (properties.warningCount[warningString] <= n){
        cerr << "WARNING: " << warning << endl;
	if (properties.warningCount[warningString] == n){
	    cerr << "(This problem occured already " << n << " times and will not be printed further)..." << endl;
	}
    }
}

void warningSummary(string const &warning, string const &warning2, Properties &properties, string warningString){
    if (properties.warningCount[warningString] == 0){return;}
    if (warning.empty()){
	cerr << "The " << warningString << " problem occured " << properties.warningCount[warningString] << " times." << endl;
    }else{
	cerr << warning << properties.warningCount[warningString] << warning2 << endl;
    }
}

/**********************************************************************
 * file:    evaluation.cc
 * licence: Artistic Licence, see file LICENCE.TXT or 
 *          http://www.opensource.org/licenses/artistic-license.php
 * descr.:  
 * authors: Mario Stanke, mario@gobics.de
 *
 * date    |   author      |  changes
 * --------|---------------|------------------------------------------
 * 19.06.02| Mario Stanke  | creation of the file
 * 30.11.02| Mario Stanke  | incomplete genes
 * 11.04.07| Mario Stanke  | evaluation of the pred. of the tts
 **********************************************************************/

#include "evaluation.hh"

// project includes
#include "gene.hh"
#include "namgene.hh"

// standard C/C++ includes
#include <iomanip>  // for setw, setprecision


void Evaluation::addToEvaluation(Transcript *predictedGenes, Transcript *annotatedGenes, Strand strand, Double quotient){
    // find flanking regions on both sides (any strand)
    leftFlankEnd=-1, rightFlankBegin=-1;
    for (Transcript* a = annotatedGenes; a != NULL; a = a->next){
	if ((leftFlankEnd == -1 && a->geneBegin() > 0) || (leftFlankEnd >=0 && a->geneBegin() - 1 < leftFlankEnd))
	    leftFlankEnd = a->geneBegin() - 1;
	if (rightFlankBegin == -1 || a->geneEnd() + 1 > rightFlankBegin)
	    rightFlankBegin = a->geneEnd() + 1;
    }
    if (rightFlankBegin == -1)
	for (Transcript* p = predictedGenes; p != NULL; p = p->next)
	    if (p->geneEnd() + 1 > rightFlankBegin)
		rightFlankBegin = p->geneEnd() + 1;

    int oldgeneFN = geneFN;
    Transcript *predFW, *predBW, *annoFW, *annoBW;
    
    predFW = Transcript::getGenesOnStrand(predictedGenes, plusstrand);
    predBW = Transcript::getGenesOnStrand(predictedGenes, minusstrand);
    annoFW = Transcript::getGenesOnStrand(annotatedGenes, plusstrand);
    annoBW = Transcript::getGenesOnStrand(annotatedGenes, minusstrand);  
    if (strand == bothstrands || strand == plusstrand) {
	addToEvaluation(predFW, annoFW);
    }
    if (strand == bothstrands || strand == minusstrand) {
	addToEvaluation(predBW, annoBW);
    }

    // temporarily for analysis
    if (geneFN > oldgeneFN) {
	string seqname = annotatedGenes ? annotatedGenes->seqname : "";
	//if (seqname)
	// cerr << "# at least 1 gene not correctly predicted in sequence " << seqname << endl;
    }

    /*
     * add the quotient to the vector of quotients
     */
    if (!(quotient < 0.0))
	quotients.push_front(quotient);

    // delete the gene sequences of the strands
    Transcript::destroyGeneSequence(predFW);
    Transcript::destroyGeneSequence(predBW);
    Transcript::destroyGeneSequence(annoFW);
    Transcript::destroyGeneSequence(annoBW);
}


/*
 * add one dataset to the evaluation
 */
void Evaluation::addToEvaluation(Transcript* predictedGeneList, Transcript* annotatedGeneList){
#ifdef DEBUG
    // for analysis only
    if (!predictedGeneList && annotatedGeneList)
	cerr << "In strand containing " << annotatedGeneList->id << " nothing was found." << endl;
    // end for analysis only
    for (Transcript *pg = predictedGeneList; pg != NULL; pg = pg->next)
	for (State *intron = pg->introns; intron != NULL; intron = intron->next)
	    if (intron->length()>longestPredIntronLen)
		longestPredIntronLen = intron->length();
#endif

    /*
     * Initialize variables
     */
     nukTPBK = nukFPBK = nukFPBKinside = nukFNBK = 0;
     exonTPBK = exonFPBK = exonFNBK = 0;

    /*
     * First make two lists just of all exons
     * Then let them contain only the unique exons.
     */
    list<State> *predictedExons=new list<State>, *annotatedExons=new list<State>;
    State *st;
    const Transcript *g;

    g = predictedGeneList;
    while (g) {
	st = g->exons;
	while (st) {
	    predictedExons->push_back(*st);
	    st = st->next;
	}
	g = g->next;
    }
 
    g = annotatedGeneList;
    while (g) {
	st = g->exons;
	while (st) {
	    annotatedExons->push_back(*st);
	    st = st->next;
	}
 	g = g->next;
    }   
    numAnnoExons += annotatedExons->size();
    numPredExons += predictedExons->size();
    /*
     * make the list of unique
     */
    
    predictedExons->sort();
    annotatedExons->sort();
    predictedExons->unique();
    annotatedExons->unique();

    numUniqueAnnoExons += annotatedExons->size();
    numUniquePredExons += predictedExons->size();


    // evaluateQuickOnNukleotideLevel(predictedExons.first, -1, annotatedExons.first, -1);
    evaluateOnNucleotideLevel(predictedExons, annotatedExons); // coding bases only
    evaluateOnExonLevel(predictedExons, annotatedExons);
    evaluateOnGeneLevel(predictedGeneList, annotatedGeneList);
    evaluateOnUTRLevel(predictedGeneList, annotatedGeneList);

    /*
     * update the derived values
     */ 
    // nucleotide level
    nukSens = (double) nukTP/(nukTP + nukFN);
    nukSpec = (double) nukTP/(nukTP + nukFP);

    // exon level
    exonFP = exonFP_wrong + exonFP_partial + exonFP_overlapping;
    exonFN = exonFN_wrong + exonFN_partial + exonFN_overlapping;
    exonSens = (double) exonTP/(exonTP+exonFN);
    exonSpec = (double) exonTP/(exonTP+exonFP);

    // gene level
    geneFP = numPredGenes - geneTP;
    geneSens = (double) geneTP / (geneTP + geneFN);
    geneSpec = (double) geneTP / (geneTP + geneFP);

    // UTR level see print();
    UTRexonSens = (double) UTRexonTP/(UTRexonTP+UTRexonFN);
    UTRexonSpec = (double) UTRexonTP/(UTRexonTP+UTRexonFP);
    nucUSens = (double) nucUTP/(nucUTP + nucUFN);
    nucUSpec = (double) nucUTP/(nucUTP + nucUFP);

    delete predictedExons;
    delete annotatedExons;
    numDataSets++;
}

void Evaluation::evaluateQuickOnNucleotideLevel(State* const predictedExon, int curPredBegin, 
						State* const annotatedExon, int curAnnoBegin){
    if (curPredBegin < 0 && predictedExon) {
	curPredBegin = predictedExon->begin;
    }
    if (curAnnoBegin < 0 && annotatedExon) {
	curAnnoBegin = annotatedExon->begin;
    }
    if (predictedExon == NULL && annotatedExon == NULL) {
	return;
    } else if (annotatedExon == NULL) {
	// only predicted Exons
	nukFP   += predictedExon->end - curPredBegin + 1;
	nukFPBK += predictedExon->end - curPredBegin + 1;
	evaluateQuickOnNucleotideLevel(predictedExon->next, -1, annotatedExon, curAnnoBegin);
    } else if (predictedExon == NULL) {
	nukFN   += annotatedExon->end - curAnnoBegin + 1;
	nukFNBK += annotatedExon->end - curAnnoBegin + 1;
	evaluateQuickOnNucleotideLevel(predictedExon, curPredBegin, annotatedExon->next, -1);
    } else {
	// both a predicted and an annotated exon
	if (curPredBegin < curAnnoBegin) {
	    // predicted exon begins first
	    if (predictedExon->end < annotatedExon->begin) {
		// predicted exon completely before annotated exon
		nukFP   += predictedExon->end - curPredBegin + 1;
		nukFPBK += predictedExon->end - curPredBegin + 1;
		evaluateQuickOnNucleotideLevel(predictedExon->next, -1, annotatedExon, curAnnoBegin);
	    } else {
		// part of predicted exon before first annotated exon
		nukFP   += curAnnoBegin - curPredBegin;
		nukFPBK += curAnnoBegin - curPredBegin;
		evaluateQuickOnNucleotideLevel(predictedExon, curAnnoBegin, annotatedExon, curAnnoBegin);
	    }
	} else if (curAnnoBegin < curPredBegin){
	    // annotated exon begins first
	     if (annotatedExon->end < predictedExon->begin) {
		// annotated exon completely before predicted exon
		nukFN   += annotatedExon->end - curAnnoBegin + 1;
		nukFNBK += annotatedExon->end - curAnnoBegin + 1;
		evaluateQuickOnNucleotideLevel(predictedExon, curPredBegin, annotatedExon->next, -1);
	    } else {
		// part of annotated exon before first predicted exon
		nukFN   += curPredBegin - curAnnoBegin;
		nukFNBK += curPredBegin - curAnnoBegin;
		evaluateQuickOnNucleotideLevel(predictedExon, curPredBegin, annotatedExon, curPredBegin);
	    }
	} else {
	    // annotated and predicted exon begin at the same position
	    if (predictedExon->end < annotatedExon->end) {
		// predicted exon shorter
		nukTP   += predictedExon->end - curPredBegin + 1;
		nukTPBK += predictedExon->end - curPredBegin + 1;
		evaluateQuickOnNucleotideLevel(predictedExon->next, -1, annotatedExon, predictedExon->end + 1); 
	    } else if (annotatedExon->end < predictedExon->end) {
		// annotated exon shorter	
		nukTP   += annotatedExon->end - curAnnoBegin +1;
		nukTPBK += annotatedExon->end - curAnnoBegin +1;
		evaluateQuickOnNucleotideLevel(predictedExon, annotatedExon->end +1, annotatedExon->next, - 1);
	    } else {
		// predicted and annotated exon have the same end and beginning :-)
		nukTP   += predictedExon->end - curPredBegin + 1;
		nukTPBK += predictedExon->end - curPredBegin + 1;
		evaluateQuickOnNucleotideLevel(predictedExon->next, -1, annotatedExon->next, - 1);
	    }
	}
    }
}

void Evaluation::evaluateOnNucleotideLevel(list<State> *predictedExon, list<State> *annotatedExon, bool UTR){
    
  list<State>::iterator pred, anno;
  int n=0;
  char *nuc;
  // find the last state in each list
  for (anno = annotatedExon->begin(); anno != annotatedExon->end(); anno++){
      if (anno->end > n)
	  n = anno->end;
  }
  for (pred = predictedExon->begin(); pred != predictedExon->end(); pred++)
    if (pred->end > n)
      n = pred->end;

  /*
   * classify each nucleotide into these four cases
   * first bit (1): annotated?
   * second bit (2): predicted?
   * third bit (4): predicted and in flanking region?
   */
   try {
    nuc = new char[n+1];
  } catch (...) {
    throw ProjectError("Not enough memory for evaluation on base level.");
    }

   for (int i=0; i<=n; i++)
     nuc[i] = 0;
  
    /*
     * Loop over the predicted exons and find the predicted coding nucs 
     */
  for (pred = predictedExon->begin(); pred != predictedExon->end(); pred++) {
      for (int i = (pred->begin >=0)? pred->begin : 0; i <= pred->end; i++){
	  nuc[i] |= 2;
	  if (i > leftFlankEnd && i < rightFlankBegin)
	      nuc[i] |= 4;
      }
  }
  /*
   * Loop over the annotated exons and find the annotated coding nucs 
   */
  for (anno = annotatedExon->begin(); anno != annotatedExon->end(); anno++) {
    for (int i = (anno->begin >= 0)? anno->begin : 0; i <= anno->end; i++)
      nuc[i] |= 1;
  }
  
  /*
   * Now count the TP, FP and FN
   */
  for (int i=0; i<=n; i++){
      int f = nuc[i];
      if (f == 1){
	if (!UTR) {nukFNBK++;} else {nucUFN++;}
      }
      if((f & 1) == 0 && (f & 2) ){
	if (!UTR) {nukFPBK++;} else {nucUFP++;}
      }
      if((f & 1) == 0 && (f & 4) ){
	if (!UTR) {nukFPBKinside++;} else {nucUFPinside++;}
      }
      if ((f & 1) && (f & 2) ){
	if (!UTR) {nukTPBK++;} else {nucUTP++;}
      }
    }
  if (!UTR){
    nukFN += nukFNBK;
    nukFP += nukFPBK;
    nukFPinside += nukFPBKinside;
    nukTP += nukTPBK;
  }
  delete []nuc;
}

void Evaluation::evaluateQuickOnExonLevel(State* predictedExon, State* annotatedExon){
    
    State *examined, *firstOverlap, *lastOverlap;

    /*
     * First go through the predicted exons and check for each one whether it is
     * a true positive, or one of the three kinds of false positives.
     */

    firstOverlap = annotatedExon;
    for (examined = predictedExon; examined != NULL; examined = examined->next){
	while (firstOverlap && firstOverlap->end < examined->begin) {
	    firstOverlap = firstOverlap->next;
	}
	if (firstOverlap == NULL || firstOverlap->begin > examined->end) {
	    exonFP_wrong++;
	    exonFPBK++;
	} else if ((firstOverlap->begin == examined->begin) &&
		   (firstOverlap->end == examined->end)) {
	    exonTP++;
	    exonTPBK++;
	} else {
	    lastOverlap = firstOverlap;
	    while (lastOverlap && lastOverlap->end < examined->end) {
		lastOverlap = lastOverlap->next;
	    }
	    if ((firstOverlap->begin == examined->begin) || 
		(lastOverlap && (lastOverlap->end == examined->end))) {
		exonFP_partial++;
		exonFPBK++;
	    } else {
		exonFP_overlapping++;
		exonFPBK++;
	    }
	}
    }  
    /*
     * Now go through the annotated exons and check for each one whether it is
     * one of the three kinds of false negatives.
     */
    firstOverlap = predictedExon;
    for (examined = annotatedExon; examined != NULL; examined = examined->next){
	while (firstOverlap && firstOverlap->end < examined->begin) {
	    firstOverlap = firstOverlap->next;
	}
	if (firstOverlap == NULL || firstOverlap->begin > examined->end) {
	    exonFN_wrong++;
	    exonFNBK++;
	} else if ((firstOverlap->begin == examined->begin) &&
		   (firstOverlap->end == examined->end)) {
	    // do nothing: true positives have been counted above
	} else {
	    lastOverlap = firstOverlap;
	    while (lastOverlap && lastOverlap->end < examined->end) {
		lastOverlap = lastOverlap->next;
	    }
	    if ((firstOverlap->begin == examined->begin) || 
		(lastOverlap && (lastOverlap->end == examined->end))) {
		exonFN_partial++;
		exonFNBK++;
	    } else {
		exonFN_overlapping++;
		exonFNBK++;
	    }
	}
    }
}

void Evaluation::evaluateOnExonLevel(list<State> *predictedExon, list<State> *annotatedExon, bool UTR){
    
    list<State>::iterator examined, anno, pred;
    int klasse;
    /*
     * First go through the predicted exons and check for each one whether it is
     * a true positive, or one of the four kinds of false positives.
     * 0 : wrong     FP
     * 1 : overlap   FP
     * 2 : partial   FP, one end exact
     * 3 : almost    FP, one end exact, other end at most UTRoffThresh bp off
     * 4 : correct   TP
     */

    for (examined = predictedExon->begin(); examined != predictedExon->end(); examined++){
	klasse = 0;
	for (anno = annotatedExon->begin(); anno != annotatedExon->end(); anno++){
	    // overlaps with the annotated exon?
	    if (!(examined->begin > anno->end || examined->end < anno->begin)) 
		if (klasse<1) 
		    klasse = 1;
	    // one splice site correct?
	    if (examined->begin == anno->begin || examined->end == anno->end)
		if (klasse < 2)
		    klasse = 2;
	    // almost correct?
	    if (abs(examined->begin-anno->begin) <= UTRoffThresh && abs(examined->end-anno->end) <= UTRoffThresh)
	      if (klasse < 3)
		klasse = 3;
	    // correct?
	    if (examined->begin == anno->begin && examined->end == anno->end)
		if (klasse < 4)
		    klasse = 4;
	}
	switch (klasse) {
	  case 0: if (!UTR) {exonFP_wrong++;       exonFPBK++; } else {UTRexonFP++;} break;
	  case 1: if (!UTR) {exonFP_overlapping++; exonFPBK++; } else {UTRexonFP++;} break;
	  case 2: if (!UTR) {exonFP_partial++;     exonFPBK++; } else {UTRexonFP++;} break;
          case 3: if (!UTR) {exonFP_partial++;     exonFPBK++; } else {UTRexonTP++;} break; // UTR exons less stringent
	  case 4: if (!UTR) {exonTP++;             exonTPBK++; } else {UTRexonTP++;} break;
	    default: break;
	}
    }

    /*
     * Now go through the annotated exons and check for each one whether it is
     * one of the three kinds of false negatives.
     * 0 : wrong     FN
     * 1 : overlap   FN
     * 2 : partial   FN
     * 3 : almost    FN (UTR:TP)
     * 3 : correct   TP (not counted here, because exons could be annotated multiply)
     */

    for (examined = annotatedExon->begin(); examined != annotatedExon->end(); examined++){
	klasse = 0;
	for (pred = predictedExon->begin(); pred !=  predictedExon->end(); pred++) {
	    // overlaps with the predicted exon?
	    if (!(examined->begin > pred->end || examined->end < pred->begin)) 
		if (klasse<1) 
		    klasse = 1;
	    // one splice site correct?
	    if (examined->begin == pred->begin || examined->end == pred->end)
		if (klasse < 2)
		    klasse = 2;
	    if (abs(examined->begin - pred->begin) <= UTRoffThresh && abs(examined->end - pred->end) <= UTRoffThresh)
              if (klasse < 3)
                klasse = 3;
	    // correct?
	    if (examined->begin == pred->begin && examined->end == pred->end)
		if (klasse < 4)
		    klasse = 4;
	}
	switch (klasse) {
	    case 0: if (!UTR) {exonFN_wrong++; }       else {UTRexonFN++;} break;
	    case 1: if (!UTR) {exonFN_overlapping++; } else {UTRexonFN++;} break;
	    case 2: if (!UTR) {exonFN_partial++; }     else {UTRexonFN++;} break;
	    case 3: if (!UTR) {exonFN_partial++; }     else {            } break; // UTR exons less stringent
	    default: break;
	}
    }
}

void Evaluation::evaluateQuickOnGeneLevel(Transcript* const predictedGeneList, Transcript* const annotatedGeneList){
   Transcript *examined, *firstOverlap;
   State *aexon, *pexon;

   /*
    * Check for each annotated Gene, whether it is completely correct predicted
    */
   firstOverlap = predictedGeneList;
   for(examined = annotatedGeneList; examined != NULL; examined = (Gene*) examined->next) {
       numAnnoGenes++;
       while (firstOverlap && firstOverlap->exons->begin < examined->exons->begin)
	   firstOverlap = firstOverlap->next;

       if (firstOverlap) {
	   aexon = examined->exons;
	   pexon = firstOverlap->exons;
	   while (aexon && pexon && aexon->begin == pexon->begin && aexon->end == pexon->end) {
	       aexon = aexon->next;
	       pexon = pexon->next;
	   }
	   if (aexon == NULL && pexon == NULL && firstOverlap->complete) {
	       geneTP++;
	   } else {
	       geneFN++;
	   }
       } else {
	   geneFN++;
       }
   }

   /*
    * Count the number of ((complete)) predicted Genes
    */
   for (examined = predictedGeneList; examined != NULL; examined = examined->next)
       numPredGenes++;
}

void Evaluation::evaluateOnGeneLevel(Transcript* const predictedGeneList, Transcript* const annotatedGeneList){
   Transcript *examined, *pred, *previous;
   State *aexon, *pexon;
   bool correctPredicted;
   /*
    * Check for each annotated Gene, whether it is completely correct predicted
    */

   for(examined = annotatedGeneList; examined != NULL; examined = (Gene*) examined->next) {
       numAnnoGenes++;
       correctPredicted = false;
       for (pred = predictedGeneList; pred != NULL && !correctPredicted; pred = pred->next) {
	   // check whether examined is equal to pred
	   aexon = examined->exons;
	   pexon = pred->exons;
	   while (aexon && pexon && aexon->begin == pexon->begin && aexon->end == pexon->end) {
	       aexon = aexon->next;
	       pexon = pexon->next;
	   }
	   if (aexon == NULL && pexon == NULL && pred->complete) {
	       correctPredicted = true;
	   }
       }
       if (correctPredicted)
	   geneTP++;
       else
	   geneFN++;
   }

   /*
    * Count the number of ((complete)) predicted CDS
    */
   for (examined = predictedGeneList; examined != NULL; examined = examined->next) {
       bool seen = false;
       for (previous = predictedGeneList; previous != examined && previous != NULL && !seen; previous = previous->next) {
	   if (*previous == *examined)
	       seen = true;
       }
       if (!seen)
	   numPredGenes++;
   }
}

void Evaluation::evaluateOnUTRLevel(Transcript* const predictedGeneList, Transcript* const annotatedGeneList){
  Transcript *examined, *pred;
  Gene * examinedG, *predG;
  int predTSS, annoTSS, predTIS, annoTIS, diff;
  int predTTS, annoTTS, predSTP, annoSTP;
  /*
   * Check for each predicted gene, whether it has a TSS
   * and whether an annotated gene with the same translation start site exists
   */
  for (pred = predictedGeneList; pred != NULL; pred = pred->next) {
      predG = dynamic_cast<Gene*>(pred);
      if (!predG)
	  continue; // ignore non-coding genes;
      
      predTSS = (predG->strand == plusstrand)? predG->transstart : predG->transend;
      predTIS = (predG->strand == plusstrand)? predG->codingstart : predG->codingend;
      if (predTSS >= 0 && predG->complete5utr){
	  numTotalPredTSS++;
	  for (examined = annotatedGeneList; examined != NULL; examined = examined->next) {
	      examinedG = dynamic_cast<Gene*>(examined);
	      if (!examinedG)
		  continue; // ignore non-coding genes;
   
	      annoTSS = (examinedG->strand == plusstrand)? examinedG->transstart : examinedG->transend;
	      annoTIS = (examinedG->strand == plusstrand)? examinedG->codingstart : examinedG->codingend;
	      if (annoTIS == predTIS && annoTSS >= 0){
		  diff = predTSS - annoTSS;
		  if (diff < 0)
		      diff = -diff;
		  numTSS++;
		  if (diff <= MAXUTRDIST)
		      tssDist[diff]++;
	      }
	  }
      }
  }
  /*
   * Check for each predicted gene, whether it has a TTS
   * and whether an annotated gene with the same stop codon exists
   */
  for (pred = predictedGeneList; pred != NULL; pred = pred->next) {
      predG = dynamic_cast<Gene*>(pred);
      if (!predG)
	  continue; // ignore non-coding genes;

      predTTS = (predG->strand == plusstrand)? predG->transend : predG->transstart;
      predSTP = (predG->strand == plusstrand)? predG->codingend : predG->codingstart;
      if (predTTS >= 0 && predG->complete3utr){
	  numTotalPredTTS++;
	  for (examined = annotatedGeneList; examined != NULL; examined = examined->next) {
	      examinedG = dynamic_cast<Gene*>(examined);
	      if (!examinedG)
		  continue; // ignore non-coding genes;
	      annoTTS = (examinedG->strand == plusstrand)? examinedG->transend : examinedG->transstart;
	      annoSTP = (examinedG->strand == plusstrand)? examinedG->codingend : examinedG->codingstart;
	      if (annoSTP == predSTP && annoTTS >= 0){
		  diff = predTTS - annoTTS;
		  if (diff<0)
		      diff = -diff;
		  numTTS++;
		  if (diff <= MAXUTRDIST)
		      ttsDist[diff]++;
	      }
	  }
      }
  }

  /* evaluate on UTR exons
   * First make two lists just of all UTR exons
   * Then let them contain only the unique exons.
   */
  list<State> *predictedUTRExons=new list<State>, *annotatedUTRExons=new list<State>;
  State *st;
  const Transcript *t;
  const Gene *g;

  t = predictedGeneList;
  while (t) {
      g = dynamic_cast<const Gene*>(t);
      if (g) {
	  st = g->utr5exons;
	  while (st) {
	      predictedUTRExons->push_back(*st);
	      st = st->next;
	  }
	  st = g->utr3exons;
	  while (st) {
	      predictedUTRExons->push_back(*st);
	      st = st->next;
	  }
      }
      t = t->next;
  }

  t = annotatedGeneList;
  while (t) {
      g = dynamic_cast<const Gene*>(t);
      if (g) {
	  st = g->utr5exons;
	  while (st) {
	      annotatedUTRExons->push_back(*st);
	      st = st->next;
	  }
	  st = g->utr3exons;
	  while (st) {
	      annotatedUTRExons->push_back(*st);
	      st = st->next;
	  }
      }
      t = t->next;
  }
  numAnnoUTRExons += annotatedUTRExons->size();
  numPredUTRExons += predictedUTRExons->size();
  /*
   * make the lists of UTR exons unique
   */

  predictedUTRExons->sort();
  annotatedUTRExons->sort();
  predictedUTRExons->unique();
  annotatedUTRExons->unique();

  numUniqueAnnoUTRExons += annotatedUTRExons->size();
  numUniquePredUTRExons += predictedUTRExons->size();

  evaluateOnExonLevel(predictedUTRExons, annotatedUTRExons, true);
  evaluateOnNucleotideLevel(predictedUTRExons, annotatedUTRExons, true);
}

void Evaluation::finishEvaluation(){
    // UTR level evaluation
    // transcription start site
    if (numTSS>0){
	meanTssDist = 0.0;
	int numInRange=0;
	medianTssDist = -1;
	int cumNum=0;
	cout << "TSS distances " << endl;
	for (int i=0; i <= MAXUTRDIST; i++){
	    if (tssDist[i]>0) {
		cumNum += tssDist[i];
		if (2*cumNum >= numTSS && medianTssDist < 0)
		    medianTssDist = i;
		meanTssDist += i*tssDist[i];
		numInRange += tssDist[i];
		cout << i << "\ttimes:" << tssDist[i] << endl;
	    }
	}
	if (numTSS-numInRange > 0)
	    cout << "Warning: " << numTSS-numInRange << " TSS are off by more than " << MAXUTRDIST << endl;
	meanTssDist += MAXUTRDIST * (numTSS-numInRange); // assume that every distance > MAXUTRDIST is equal to MAXUTRDIST
	meanTssDist /= numTSS;
    } else {
	medianTssDist = -1;
	meanTssDist = -1;
    }
    // transcription termination site
    if (numTTS>0){
	meanTtsDist = 0.0;
	int numInRange = 0;
	medianTtsDist = -1;
	int cumNum=0;
	cout << "TTS distances " << endl;
	for (int i=0; i <= MAXUTRDIST; i++){
	    if (ttsDist[i]>0) {
		cumNum += ttsDist[i];
		if (2*cumNum >= numTTS && medianTtsDist < 0)
		    medianTtsDist = i;
		meanTtsDist += i*ttsDist[i];
		numInRange += ttsDist[i];
		cout << i << "\ttimes:" << ttsDist[i] << endl;
	    }
	}
	if (numTTS-numInRange > 0)
	    cout << "Warning: " << numTTS-numInRange << " TTS are off by more than " << MAXUTRDIST << endl;
	meanTtsDist += MAXUTRDIST * (numTTS-numInRange); // assume that every distance > MAXUTRDIST is equal to MAXUTRDIST
	meanTtsDist /= numTTS;
    } else {
	medianTtsDist = -1;
	meanTtsDist = -1;
    }
}

void Evaluation::print(){
    cout << "\n*******      Evaluation of gene prediction     *******" << endl << endl;

    // nucleotide level
    cout << "---------------------------------------------\\" << endl;
    cout << setw(16) << " " << " | " << setw(11) << "sensitivity" << " | " 
	 << setw(11) << "specificity" << " |" << endl;
    cout << "---------------------------------------------|" << endl;
    cout << setw(16) << "nucleotide level" << " | " << setw(11)
	 << setprecision(3) << nukSens 
	 << " | " <<  setw(11) << setprecision(3) << nukSpec << " |" << endl;
    cout << "---------------------------------------------/" << endl << endl;

    // exon level
    cout << "-----------------------------------------------------------------------------" 
	 << "-----------------------------\\" << endl;
    cout << setw(10) << " " << " | " << setw(6) << "#pred" << " | " << setw(6) << "#anno"
	 << " | " << setw(4) << " " << " | " << setw(18) << "FP = false pos." 
	 << " | " << setw(18) << "FN = false neg." << " | " 
	 << setw(11) << " " << " | " 
	 << setw(11) << " " << " |" << endl;
    cout << setw(10) << " " << " | " << setw(6) << "total/" << " | " << setw(6) << "total/"
	 << " | "  << setw(4) << "TP" << " |"<< setw(19) << "--------------------" 
	 << "|" << setw(19) << "--------------------" 
	 << "| " << setw(11) << "sensitivity" << " | " 
	 << setw(11) << "specificity" << " |" << endl;
    cout << setw(10) << " " << " | " << setw(6) << "unique" << " | " << setw(6) << "unique"
	 << " | " << setw(4) << " " << " | " << setw(4) << "part" << " | " << setw(4) 
	 << "ovlp" << " | "  << setw(4) << "wrng" << " | " <<  setw(4) << "part" 
	 << " | " << setw(4) << "ovlp" << " | " 
	 << setw(4) << "wrng" << " | " << setw(11) << " " << " | " 
	 << setw(11) << " " << " |" << endl;
    cout << "-----------------------------------------------------------------------------" 
	 << "-----------------------------|" << endl;
    cout << setw(10) << " " << " | " << setw(6) << " " << " | " << setw(6) << " "
	 << " | " << setw(4) << " " << " | " << setw(18) << exonFP << " | " << setw(18) << exonFN 
	 << " | " << setw(11) << " " << " | " 
	 << setw(11) << " " << " |" << endl;
    cout << setw(10) << "exon level" << " | " << setw(6) << numPredExons << " | " 
	 << setw(6) << numAnnoExons << " | " << setw(4) << exonTP << " | "
	 << setw(18) << "------------------" << " | " << setw(18) << "------------------" 
	 << " | " << setw(11) << setprecision(3) << exonSens << " | " 
	 << setw(11) << setprecision(3) << exonSpec << " |" << endl;
    cout << setw(10) << " " << " | " << setw(6) << numUniquePredExons << " | " << setw(6) << numUniqueAnnoExons
	 << " | " << setw(4) << " " << " | " << setw(4) << exonFP_partial << " | " << setw(4) 
	 << exonFP_overlapping << " | "  << setw(4) << exonFP_wrong << " | " <<  setw(4) 
	 << exonFN_partial << " | " << setw(4) << exonFN_overlapping << " | " 
	 << setw(4) << exonFN_wrong << " | " << setw(11) << " " << " | " 
	 << setw(11) << " " << " |" << endl;
    cout << "-----------------------------------------------------------------------------" 
	 << "-----------------------------/" << endl << endl;
    // transcript level
    cout << "-----------------------------------------------------" 
	 << "-----------------------\\" << endl;
    cout << setw(10) << "transcript" << " | " << setw(5) << "#pred" << " | " 
	 << setw(5) << "#anno" << " | " << setw(4) << "TP" << " | " 
	 << setw(4) << "FP" << " | " << setw(4) << "FN" << " | " 
	 << setw(11) << "sensitivity" << " | " 
	 << setw(9) << "specificity" << " |" << endl;
    cout << "------------------------------------------------------" 
	 << "----------------------|" << endl;
    cout << setw(10) << "gene level" << " | " << setw(5) << numPredGenes << " | " 
	 << setw(5) << numAnnoGenes << " | " << setw(4) << geneTP << " | " 
	 << setw(4) << geneFP << " | " << setw(4) << geneFN << " | " 
	 << setw(11) << setprecision(3) << geneSens << " | " 
	 << setw(11) << setprecision(3) << geneSpec << " |" << endl;
    cout << "-------------------------------------------------------" 
	 << "---------------------/" << endl;
    if (numTotalPredTSS > 0 || numTotalPredTTS > 0){
      // utr level
      cout << endl;
      cout << "------------------------------------------------------------------------\\" << endl;
      cout << setw(15) << "UTR" << " | "  << setw(10) << "total pred" << " | " << setw(14) << "CDS bnd. corr." << " | "
	   << setw(10) << "meanDiff" << " | "<< setw(10) << "medianDiff" << " |" << endl;
      cout << "------------------------------------------------------------------------|" << endl;
      cout << setw(15) << "TSS" << " | "  << setw(10) << numTotalPredTSS << " | " << setw(14) << numTSS << " | " 
	   << setw(10) << setprecision(3) << meanTssDist << " | " << setw(10) << medianTssDist << " |" << endl;
      cout << setw(15) << "TTS" << " | "  << setw(10) << numTotalPredTTS << " | " << setw(14) << numTTS << " | " 
	   << setw(10) << setprecision(3) << meanTtsDist << " | " << setw(10) << medianTtsDist << " |" << endl;
      cout << "------------------------------------------------------------------------|" << endl;
      cout << setw(15) << "UTR" << " | "  << setw(10) << "uniq. pred" << " | " << setw(14) << "unique anno" << " | "
           << setw(10) << "   sens." << " | "<< setw(10) << "     spec." << " |" << endl;
      cout << "------------------------------------------------------------------------|" << endl;
      cout << setw(15) << " " << " | "  << setw(45) << "true positive = 1 bound. exact, 1 bound. <= " << UTRoffThresh << "bp off" " |" << endl;
      cout << setw(15) << "UTR exon level" << " | "  << setw(10) << numUniquePredUTRExons << " | " << setw(14) << numUniqueAnnoUTRExons << " | "
           << setw(10) << setprecision(3) << UTRexonSens << " | " << setw(10) << UTRexonSpec << " |" << endl;
      cout << "------------------------------------------------------------------------|" << endl;
      cout << setw(15) << "UTR base level" << " | "  << setw(10) << (nucUTP+nucUFP) << " | " << setw(14) << (nucUTP+nucUFN) << " | "
           << setw(10) << setprecision(3) << nucUSens << " | " << setw(10) << nucUSpec << " |" << endl;
      cout << "------------------------------------------------------------------------/" << endl;
      cout << "nucUTP= " << nucUTP << " nucUFP=" << nucUFP << " nucUFPinside= " << nucUFPinside << " nucUFN=" << nucUFN << endl;
  
    }


#ifdef DEBUG
    cout << "longest predicted intron had length " << longestPredIntronLen << endl;
#endif
}

void Evaluation::printQuotients(){
    cout << "a-posteriori probability of viterbi path" << endl;
    cout << "----------------------------------------" << endl;
    cout << "a-posteriori probability of correct path" << endl << endl;

    quotients.sort();
    int one=0, ten=0;
    list<Double>::iterator it = quotients.begin();
    for (; it != quotients.end() && *it < 1.000001; it++)
	one++;
    cout << one << " times were the paths equally likely (identical)." << endl;
    cout << "sorted quotients of the rest:" << endl;
    for (; it != quotients.end(); it++){
	cout << *it << endl;
	if (*it < 10.0)
	    ten++;
    }
    cout << endl << ten << " quotients were between 1 and 10" << endl;
    cout << endl;
}


/*
 * predictAndEvaluate
 */

Evaluation* predictAndEvaluate(vector<AnnoSequence*> trainGeneList, FeatureCollection &extrinsicFeatures){
    Evaluation *eval = new Evaluation();
    NAMGene namgene;
    StatePath *viterbiPath, *condensedViterbiPath;
    bool crf = inCRFTraining;
    inCRFTraining = 0; // this way Viterbi will be much more efficient because of use of static variables
    for (int i=0; i < trainGeneList.size(); i++){
	AnnoSequence* evalGene = trainGeneList[i];
	SequenceFeatureCollection& sfc=extrinsicFeatures.getSequenceFeatureCollection(evalGene->seqname);
	sfc.prepare(evalGene, false);
	viterbiPath = namgene.getTrainViterbiPath(evalGene->sequence, &sfc);
	condensedViterbiPath = StatePath::condenseStatePath(viterbiPath);
	Transcript* viterbiGenes = condensedViterbiPath->projectOntoGeneSequence("v");
	// condensedViterbiPath->print();
	eval->addToEvaluation(viterbiGenes, evalGene->anno->genes, bothstrands);
	delete viterbiPath;
	delete condensedViterbiPath;
	Transcript::destroyGeneSequence(viterbiGenes);
    }
    inCRFTraining = crf;
    return eval;
}

/***************************************************************************\
 * Filename : evaluation.hh
 * Author   : Mario Stanke
 * Email    : stanke@math.uni-goettingen.de
 *
 * Copyright: (C) 2002 by Mario Stanke
 *
 * Description: 
 *
 *
 * Date       |   Author      |  Changes
 *------------+---------------+---------------------------------
 * 19.06.2002 | Mario Stanke  | Creation of the file
 * 11.04.2007 | Mario Stanke  | Evaluation of the pred. of the transcription termination site (tts)
 \**************************************************************************/

#ifndef _EVALUATION_HH
#define _EVALUATION_HH

// project includes
#include "types.hh"
#include "gene.hh"
#include "extrinsicinfo.hh"

// standard C/C++ includes
#include <list>
#include <vector>


#define MAXUTRDIST 5000
class Evaluation {
public:
  Evaluation(){
    nukTP = nukFP = nukFN = nukFPinside = 0;
    nucUTP = nucUFP = nucUFN = nucUFPinside = 0;
    exonTP = exonFP_partial = exonFP_overlapping = exonFP_wrong = 0;
    exonFN_partial = exonFN_overlapping = exonFN_wrong = 0;
    UTRexonTP = UTRexonFP = UTRexonFN = 0;
    UTRoffThresh = 20;
    geneTP = geneFN = 0;
    numPredExons = numAnnoExons = 0;
    numPredUTRExons = numAnnoUTRExons = 0;
    numUniquePredExons = numUniqueAnnoExons = 0;
    numUniquePredUTRExons = numUniqueAnnoUTRExons = 0;
    numPredGenes = numAnnoGenes = 0;
    numDataSets = 0;
    longestPredIntronLen = 0;
    tssDist = new int[MAXUTRDIST+1];
    for (int i=0; i<= MAXUTRDIST; i++)
      tssDist[i] = 0;
    numTotalPredTSS = numTSS = 0;
    ttsDist = new int[MAXUTRDIST+1];
    for (int i=0; i<= MAXUTRDIST; i++)
      ttsDist[i] = 0;
    numTotalPredTTS = numTTS = 0;
    leftFlankEnd = rightFlankBegin = -1;
  }
  ~Evaluation(){
      if (tssDist)
	  delete [] tssDist;
      if (ttsDist)
	  delete [] ttsDist;
  };

  void addToEvaluation(Transcript* prediction, Transcript *database, Strand strand, Double quotient = -1.0);
  void addToEvaluation(Transcript* predictedGeneList, Transcript* annotatedGeneList);
  void finishEvaluation();
  void print();
  void printQuotients();
private:
  /*
   * Quick evaluation is fast but requires that both gene lists 
   *
   */
  void evaluateQuickOnNucleotideLevel(State* const predictedExon, int curPredBegin, 
				      State* const annotatedExon, int curAnnoBegin);
  void evaluateQuickOnExonLevel(State* predictedExon, State* annotatedExon);
  void evaluateQuickOnGeneLevel(Transcript* const predictedGeneList, Transcript* const annotatedGeneList);

  void evaluateOnNucleotideLevel(list<State> *predictedExon, list<State> *annotatedExon, bool UTR=false);    
  void evaluateOnExonLevel(list<State> *predictedExon, list<State> *annotatedExon, bool UTR=false);
  void evaluateOnGeneLevel(Transcript* const predictedGeneList, Transcript* const annotatedGeneList);
  void evaluateOnUTRLevel(Transcript* const predictedGeneList, Transcript* const annotatedGeneList);
public:
  // nucleotide level
  int nukTP, nukFP, nukFN,
      nukFPinside; // false positive coding base inside gene area (as opposed to in flanking regions)
  int nucUTP, nucUFP, nucUFN, // UTR bases
    nucUFPinside; // false positive noncoding base inside gene area (as opposed to in flanking regions)
  double nukSens, nukSpec;         // coding base sensitivity and specifity
  double nucUSens, nucUSpec;       // non-coding base sensitivity and specifity
  double exonSens, exonSpec;       // coding exon sensitivity and specificity
  double UTRexonSens, UTRexonSpec; // exon sensitivity and specificity
  double geneSens, geneSpec;       // gene sensitivity and specifity
private:
  //TP = true positive, FP = false positive, FN = false negative
    int  leftFlankEnd, rightFlankBegin;
    list<Double> quotients;
    int longestPredIntronLen;

  // exon level
  int numPredExons, numAnnoExons;
  int numUniquePredExons, numUniqueAnnoExons;
  
  int exonTP, 
    exonFP, 
    exonFP_partial,        // predicted exon unequal to but included in an annotated exon
    exonFP_overlapping,
    exonFP_wrong,
    exonFN,
    exonFN_partial,        // annotated exon  unequal to but included in a predicted exon
    exonFN_overlapping,
    exonFN_wrong;

  // gene level
  int geneTP, geneFP, geneFN;
  int numPredGenes, numAnnoGenes;

  // UTR level
  int *tssDist; // array that holds for each distance the number of predicted TSS that is off by this distance
  int numTSS;   // number of gene pairs (anno, pred) with identical translation start and where both have an annotated TSS
  int numTotalPredTSS;
  double meanTssDist;
  int  medianTssDist;
  int *ttsDist; // array that holds for each distance the number of predicted TTS that is off by this distance
  int numTTS;   // number of gene pairs (anno, pred) with identical stop codon and where both have an annotated TTS
  int numTotalPredTTS;
  double meanTtsDist;
  int  medianTtsDist;
  int numPredUTRExons, numAnnoUTRExons;
  int numUniquePredUTRExons, numUniqueAnnoUTRExons;
  int UTRexonTP, UTRexonFP, UTRexonFN;
  int UTRoffThresh; // count UTR exon as correct, if one end is exact and the other end at most this many bp off
  /*
   * data members for the "Burge-Karlin"-Method computing first the 
   * specifity and sensitivity for each sequence and then taking their
   * means afterwards
   */

  int numDataSets;
  // nukleotide level
  int nukTPBK, nukFPBK, nukFPBKinside, nukFNBK;
           
  // exon level
  int exonTPBK, exonFPBK, exonFNBK;
};

/*
 * predictAndEvaluate
 *
 * Predict genes on the given set of annotated sequences given the current parameters.
 * Then evaluate the accuracy against the given annotation.
 */
Evaluation* predictAndEvaluate(vector<AnnoSequence*> trainGeneList, FeatureCollection &extrinsicFeatures);

#endif // _EVALUATION_HH
